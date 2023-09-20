"""BoltzTraP2 is a python software interpolating band structures and
computing materials properties from dft band structure using Boltzmann
semi-classical transport theory.
This module provides a pymatgen interface to BoltzTraP2.
Some of the code is written following the examples provided in BoltzTraP2.

BoltzTraP2 has been developed by Georg Madsen, Jesús Carrete, Matthieu J. Verstraete.

https://gitlab.com/sousaw/BoltzTraP2
https://www.sciencedirect.com/science/article/pii/S0010465518301632

References are:

    Georg K.H.Madsen, Jesús Carrete, Matthieu J.Verstraete
    BoltzTraP2, a program for interpolating band structures and
    calculating semi-classical transport coefficients
    Computer Physics Communications 231, 140-145, 2018

    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71

Todo:
- DONE: spin polarized bands
- read first derivative of the eigenvalues from vasprun.xml (mommat)
- handle magnetic moments (magmom)
"""

from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import numpy as np
from monty.serialization import dumpfn, loadfn
from tqdm import tqdm

from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine, Spin
from pymatgen.electronic_structure.boltztrap import BoltztrapError
from pymatgen.electronic_structure.dos import CompleteDos, Dos, Orbital
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Vasprun
from pymatgen.symmetry.bandstructure import HighSymmKpath

try:
    from BoltzTraP2 import bandlib as BL
    from BoltzTraP2 import fite, sphere, units
except ImportError:
    raise BoltztrapError("BoltzTraP2 has to be installed and working")

__author__ = "Francesco Ricci"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Fracesco Ricci"
__email__ = "frankyricci@gmail.com"
__status__ = "Development"
__date__ = "August 2018"


class VasprunBSLoader:
    """Loader for Bandstructure and Vasprun pmg objects."""

    def __init__(self, obj, structure=None, nelect=None):
        """
        Args:
            obj: Either a pmg Vasprun or a BandStructure object.
            structure: Structure object in case is not included in the BandStructure object.
            nelect: number of electrons in case a BandStructure obj is provided.

        Example:
            vrun = Vasprun('vasprun.xml')
            data = VasprunBSLoader(vrun)
        """
        if isinstance(obj, Vasprun):
            structure = obj.final_structure
            nelect = obj.parameters["NELECT"]
            bs_obj = obj.get_band_structure()
        elif isinstance(obj, BandStructure):
            bs_obj = obj
        else:
            raise BoltztrapError("The object provided is neither a Bandstructure nor a Vasprun.")

        self.kpoints = np.array([kp.frac_coords for kp in bs_obj.kpoints])

        if bs_obj.structure:
            self.structure = bs_obj.structure
        elif structure:
            self.structure = structure
        else:
            raise BoltztrapError("A structure must be given.")

        self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
        self.proj_all = None
        if bs_obj.projections:
            self.proj_all = {sp: p.transpose((1, 0, 3, 2)) for sp, p in bs_obj.projections.items()}

        e = np.array(list(bs_obj.bands.values()))
        e = e.reshape(-1, e.shape[-1])
        self.ebands_all = e * units.eV

        self.is_spin_polarized = bs_obj.is_spin_polarized

        if bs_obj.is_spin_polarized:
            self.dosweight = 1.0
        else:
            self.dosweight = 2.0

        self.lattvec = self.atoms.get_cell().T * units.Angstrom
        self.mommat_all = None  # not implemented yet
        self.mommat = None  # not implemented yet
        self.magmom = None  # not implemented yet
        self.fermi = bs_obj.efermi * units.eV
        self.UCvol = self.structure.volume * units.Angstrom**3

        if not bs_obj.is_metal():
            self.vbm_idx = max(bs_obj.get_vbm()["band_index"][Spin.up] + bs_obj.get_vbm()["band_index"][Spin.down])
            self.cbm_idx = min(bs_obj.get_cbm()["band_index"][Spin.up] + bs_obj.get_cbm()["band_index"][Spin.down])
            self.vbm = bs_obj.get_vbm()["energy"]
            self.cbm = bs_obj.get_cbm()["energy"]
        else:
            self.vbm_idx = self.cbm_idx = None
            self.vbm = self.fermi / units.eV
            self.cbm = self.fermi / units.eV

        if nelect:
            self.nelect_all = nelect
        elif self.vbm_idx:
            self.nelect_all = self.vbm_idx + self.cbm_idx + 1
        else:
            raise BoltztrapError("nelect must be given.")

    @classmethod
    def from_file(cls, vasprun_file):
        """Get a vasprun.xml file and return a VasprunBSLoader."""
        vrun_obj = Vasprun(vasprun_file, parse_projected_eigen=True)
        return cls(vrun_obj)

    def get_lattvec(self):
        """The lattice vectors."""
        try:
            return self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec

    def get_volume(self):
        """Volume."""
        try:
            return self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol

    def bandana(self, emin=-np.inf, emax=np.inf):
        """Cut out bands outside the range (emin,emax)."""
        bandmin = np.min(self.ebands_all, axis=1)
        bandmax = np.max(self.ebands_all, axis=1)
        ntoolow = np.count_nonzero(bandmax <= emin)
        accepted = np.logical_and(bandmin < emax, bandmax > emin)
        # self.data_bkp = np.copy(self.data.ebands)
        self.ebands = self.ebands_all[accepted]

        self.proj = {}
        if self.proj_all:
            if len(self.proj_all) == 2:
                h = int(len(accepted) / 2)
                self.proj[Spin.up] = self.proj_all[Spin.up][:, accepted[:h], :, :]
                self.proj[Spin.down] = self.proj_all[Spin.down][:, accepted[h:], :, :]
            elif len(self.proj_all) == 1:
                self.proj[Spin.up] = self.proj_all[Spin.up][:, accepted, :, :]

        if self.mommat_all:
            self.mommat = self.mommat[:, accepted, :]
        # Removing bands may change the number of valence electrons
        if self.nelect_all:
            self.nelect = self.nelect_all - self.dosweight * ntoolow

        return accepted


class BandstructureLoader:
    """Loader for Bandstructure object."""

    def __init__(self, bs_obj, structure=None, nelect=None, mommat=None, magmom=None):
        """
        Args:
            bs_obj: BandStructure object.
            structure: Structure object. It is needed if it is not contained in the BandStructure obj.
            nelect: Number of electrons in the calculation.
            mommat: Matrix of derivatives of energy eigenvalues. TODO Not implemented yet.
            magmom: Matrix of magnetic moments in non collinear calculations. Not implemented yet.

        Example:
            vrun = Vasprun('vasprun.xml')
            bs = vrun.get_band_structure()
            st = vrun.final_structure
            ne = vrun.parameters['NELECT']
            data = BandstructureLoader(bs,st,ne)
        """
        warnings.warn("Deprecated Loader. Use VasprunBSLoader instead.")

        self.kpoints = np.array([kp.frac_coords for kp in bs_obj.kpoints])

        if structure is None:
            self.structure = bs_obj.structure
        else:
            self.structure = structure

        self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
        self.proj_all = None
        if bs_obj.projections:
            self.proj_all = {sp: p.transpose((1, 0, 3, 2)) for sp, p in bs_obj.projections.items()}

        e = np.array(list(bs_obj.bands.values()))
        e = e.reshape(-1, e.shape[-1])
        self.ebands_all = e * units.eV

        self.is_spin_polarized = bs_obj.is_spin_polarized

        if bs_obj.is_spin_polarized:
            self.dosweight = 1.0
        else:
            self.dosweight = 2.0

        self.lattvec = self.atoms.get_cell().T * units.Angstrom
        self.mommat_all = mommat  # not implemented yet
        self.mommat = mommat  # not implemented yet
        self.magmom = magmom  # not implemented yet
        self.fermi = bs_obj.efermi * units.eV
        self.UCvol = self.structure.volume * units.Angstrom**3

        if not bs_obj.is_metal():
            self.vbm_idx = max(bs_obj.get_vbm()["band_index"][Spin.up] + bs_obj.get_vbm()["band_index"][Spin.down])
            self.cbm_idx = min(bs_obj.get_cbm()["band_index"][Spin.up] + bs_obj.get_cbm()["band_index"][Spin.down])
            self.vbm = bs_obj.get_vbm()["energy"]
            self.cbm = bs_obj.get_cbm()["energy"]
            self.nelect_all = self.vbm_idx * self.dosweight
        else:
            self.vbm_idx = self.cbm_idx = None
            self.vbm = self.fermi
            self.cbm = self.fermi
            self.nelect_all = nelect

    def get_lattvec(self):
        """The lattice vectors."""
        try:
            return self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec

    def bandana(self, emin=-np.inf, emax=np.inf):
        """Cut out bands outside the range (emin,emax)."""
        bandmin = np.min(self.ebands_all, axis=1)
        bandmax = np.max(self.ebands_all, axis=1)
        ntoolow = np.count_nonzero(bandmax <= emin)
        accepted = np.logical_and(bandmin < emax, bandmax > emin)
        # self.data_bkp = np.copy(self.data.ebands)
        self.ebands = self.ebands_all[accepted]

        self.proj = {}
        if self.proj_all:
            if len(self.proj_all) == 2:
                h = int(len(accepted) / 2)
                self.proj[Spin.up] = self.proj_all[Spin.up][:, accepted[:h], :, :]
                self.proj[Spin.down] = self.proj_all[Spin.down][:, accepted[h:], :, :]
            elif len(self.proj) == 1:
                self.proj[Spin.up] = self.proj_all[Spin.up][:, accepted, :, :]

        if self.mommat_all:
            self.mommat = self.mommat[:, accepted, :]
        # Removing bands may change the number of valence electrons
        if self.nelect_all:
            self.nelect = self.nelect_all - self.dosweight * ntoolow

        return accepted

    def set_upper_lower_bands(self, e_lower, e_upper):
        """Set fake upper/lower bands, useful to set the same energy
        range in the spin up/down bands when calculating the DOS.
        """
        warnings.warn(
            "This method does not work anymore in case of spin \
        polarized case due to the concatenation of bands !"
        )

        lower_band = e_lower * np.ones((1, self.ebands.shape[1]))
        upper_band = e_upper * np.ones((1, self.ebands.shape[1]))

        self.ebands = np.vstack((lower_band, self.ebands, upper_band))
        if self.proj:
            for sp, proj in self.proj.items():
                proj_lower = proj[:, 0:1, :, :]
                proj_upper = proj[:, -1:, :, :]
                self.proj[sp] = np.concatenate((proj_lower, proj, proj_upper), axis=1)

    def get_volume(self):
        """Volume."""
        try:
            return self.UCvol
        except AttributeError:
            lattvec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(lattvec))
        return self.UCvol


class VasprunLoader:
    """Loader for Vasprun object."""

    def __init__(self, vrun_obj=None):
        """vrun_obj: Vasprun object."""
        warnings.warn("Deprecated Loader. Use VasprunBSLoader instead.")

        if vrun_obj:
            self.kpoints = np.array(vrun_obj.actual_kpoints)
            self.structure = vrun_obj.final_structure
            self.atoms = AseAtomsAdaptor.get_atoms(self.structure)
            self.proj = None
            if len(vrun_obj.eigenvalues) == 1:
                e = next(iter(vrun_obj.eigenvalues.values()))
                self.ebands = e[:, :, 0].transpose() * units.eV
                self.dosweight = 2.0
                if vrun_obj.projected_eigenvalues:
                    self.proj = next(iter(vrun_obj.projected_eigenvalues.values()))

            elif len(vrun_obj.eigenvalues) == 2:
                raise BoltztrapError("spin bs case not implemented")

            self.lattvec = self.atoms.get_cell().T * units.Angstrom

            # TODO: read mommat from vasprun
            self.mommat = self.magmom = self.spin = None
            self.fermi = vrun_obj.efermi * units.eV
            self.nelect = vrun_obj.parameters["NELECT"]
            self.UCvol = self.structure.volume * units.Angstrom**3

            bs_obj = vrun_obj.get_band_structure()
            if not bs_obj.is_metal():
                self.vbm_idx = max(bs_obj.get_vbm()["band_index"][Spin.up] + bs_obj.get_vbm()["band_index"][Spin.down])
                self.cbm_idx = min(bs_obj.get_cbm()["band_index"][Spin.up] + bs_obj.get_cbm()["band_index"][Spin.down])
                self.vbm = bs_obj.get_vbm()["energy"]
                self.cbm = bs_obj.get_cbm()["energy"]
            else:
                self.vbm_idx = self.cbm_idx = None
                self.vbm = self.fermi
                self.cbm = self.fermi

    @classmethod
    def from_file(cls, vasprun_file):
        """Get a vasprun.xml file and return a VasprunLoader."""
        vrun_obj = Vasprun(vasprun_file, parse_projected_eigen=True)
        return VasprunLoader(vrun_obj)

    def get_lattvec(self):
        """Lattice vectors."""
        try:
            return self.lattvec
        except AttributeError:
            self.lattvec = self.atoms.get_cell().T * units.Angstrom
        return self.lattvec

    def bandana(self, emin=-np.inf, emax=np.inf):
        """Cut out bands outside the range (emin,emax)."""
        band_min = np.min(self.ebands, axis=1)
        band_max = np.max(self.ebands, axis=1)
        ii = np.nonzero(band_min < emax)
        n_emax = ii[0][-1]
        ii = np.nonzero(band_max > emin)
        n_emin = ii[0][0]
        # BoltzTraP2.misc.info("BANDANA output")
        # for iband in range(len(self.ebands)):
        # BoltzTraP2.misc.info(iband, bandmin[iband], bandmax[iband], (
        # (bandmin[iband] < emax) & (bandmax[iband] > emin)))
        self.ebands = self.ebands[n_emin : n_emax + 1]

        if isinstance(self.proj, np.ndarray):
            self.proj = self.proj[:, n_emin : n_emax + 1, :, :]

        if self.mommat is not None:
            self.mommat = self.mommat[:, n_emin : n_emax + 1, :]
        # Removing bands may change the number of valence electrons
        if self.nelect is not None:
            self.nelect -= self.dosweight * n_emin
        return n_emin, n_emax

    def get_volume(self):
        """Volume of cell."""
        try:
            return self.UCvol
        except AttributeError:
            latt_vec = self.get_lattvec()
            self.UCvol = np.abs(np.linalg.det(latt_vec))
        return self.UCvol


class BztInterpolator:
    """Interpolate the dft band structures."""

    def __init__(
        self,
        data,
        lpfac=10,
        energy_range=1.5,
        curvature=True,
        save_bztInterp=False,
        load_bztInterp=False,
        save_bands=False,
        fname="bztInterp.json.gz",
    ):
        """
        Args:
            data: A loader
            lpfac: the number of interpolation points in the real space. By
                default 10 gives 10 time more points in the real space than
                the number of kpoints given in reciprocal space.
            energy_range: usually the interpolation is not needed on the entire energy
                range but on a specific range around the Fermi level.
                This energy in eV fix the range around the Fermi level
                (E_fermi-energy_range,E_fermi+energy_range) of
                bands that will be interpolated
                and taken into account to calculate the transport properties.
            curvature: boolean value to enable/disable the calculation of second
                derivative related transport properties (Hall coefficient).
            save_bztInterp: Default False. If True coefficients and equivalences are
                saved in fname file.
            load_bztInterp: Default False. If True the coefficients and equivalences
                are loaded from fname file, not calculated. It can be faster than
                re-calculate them in some cases.
            save_bands: Default False. If True interpolated bands are also stored.
                It can be slower than interpolate them. Not recommended.
            fname: File path where to store/load from the coefficients and equivalences.

        Example:
            data = VasprunLoader().from_file('vasprun.xml')
            bztInterp = BztInterpolator(data)
        """
        bands_loaded = False
        self.data = data
        num_kpts = self.data.kpoints.shape[0]
        self.efermi = self.data.fermi
        middle_gap_en = (self.data.cbm + self.data.vbm) / 2
        self.accepted = self.data.bandana(
            emin=(middle_gap_en - energy_range) * units.eV,
            emax=(middle_gap_en + energy_range) * units.eV,
        )

        if load_bztInterp:
            bands_loaded = self.load(fname)
        else:
            self.equivalences = sphere.get_equivalences(self.data.atoms, self.data.magmom, num_kpts * lpfac)
            self.coeffs = fite.fitde3D(self.data, self.equivalences)

        if not bands_loaded:
            self.eband, self.vvband, self.cband = fite.getBTPbands(
                self.equivalences, self.coeffs, self.data.lattvec, curvature=curvature
            )

        if save_bztInterp:
            self.save(fname, save_bands)

    def load(self, fname="bztInterp.json.gz"):
        """Load the coefficient, equivalences, bands from fname."""
        d = loadfn(fname)
        if len(d) > 2:
            self.equivalences, coeffs, self.eband, self.vvband, self.cband = d
            bands_loaded = True
        elif len(d) == 2:
            self.equivalences, coeffs = loadfn(fname)
            bands_loaded = False
        else:
            raise BoltztrapError("Something wrong reading the data file!")
        self.coeffs = coeffs[0] + coeffs[1] * 1j
        return bands_loaded

    def save(self, fname="bztInterp.json.gz", bands=False):
        """Save the coefficient, equivalences to fname.
        If bands is True, also interpolated bands are stored.
        """
        if bands:
            dumpfn(
                [
                    self.equivalences,
                    [self.coeffs.real, self.coeffs.imag],
                    self.eband,
                    self.vvband,
                    self.cband,
                ],
                fname,
            )
        else:
            dumpfn([self.equivalences, [self.coeffs.real, self.coeffs.imag]], fname)

    def get_band_structure(self, kpaths=None, kpoints_lbls_dict=None, density=20):
        """Return a BandStructureSymmLine object interpolating bands along a
        High symmetry path calculated from the structure using HighSymmKpath
        function. If kpaths and kpoints_lbls_dict are provided, a custom
        path is interpolated.
        kpaths: List of lists of following kpoints labels defining
                the segments of the path. E.g. [['L','M'],['L','X']]
        kpoints_lbls_dict: Dict where keys are the kpoint labels used in kpaths
                and values are their fractional coordinates.
                E.g. {'L':np.array(0.5,0.5,0.5)},
                      'M':np.array(0.5,0.,0.5),
                      'X':np.array(0.5,0.5,0.)}
        density: Number of points in each segment.
        """
        if isinstance(kpaths, list) and isinstance(kpoints_lbls_dict, dict):
            kpoints = []
            for kpath in kpaths:
                for idx, k_pt in enumerate(kpath[:-1]):
                    sta = kpoints_lbls_dict[k_pt]
                    end = kpoints_lbls_dict[kpath[idx + 1]]
                    kpoints.append(np.linspace(sta, end, density))
            kpoints = np.concatenate(kpoints)
        else:
            kpath = HighSymmKpath(self.data.structure)
            kpoints = np.vstack(kpath.get_kpoints(density, coords_are_cartesian=False)[0])
            kpoints_lbls_dict = kpath.kpath["kpoints"]

        lattvec = self.data.get_lattvec()
        egrid, vgrid = fite.getBands(kpoints, self.equivalences, lattvec, self.coeffs)
        if self.data.is_spin_polarized:
            h = sum(np.array_split(self.accepted, 2)[0])
            egrid = np.array_split(egrid, [h], axis=0)
            bands_dict = {
                Spin.up: (egrid[0] / units.eV),
                Spin.down: (egrid[1] / units.eV),
            }
        else:
            bands_dict = {Spin.up: (egrid / units.eV)}

        return BandStructureSymmLine(
            kpoints,
            bands_dict,
            self.data.structure.lattice.reciprocal_lattice,
            self.efermi / units.eV,
            labels_dict=kpoints_lbls_dict,
        )

    def get_dos(self, partial_dos=False, npts_mu=10000, T=None, progress=False):
        """Return a Dos object interpolating bands.

        Args:
            partial_dos: if True, projections will be interpolated as well
                and partial doses will be return. Projections must be available
                in the loader.
            npts_mu: number of energy points of the Dos
            T: parameter used to smooth the Dos
            progress: Default False, If True a progress bar is shown when
                partial dos are computed.
        """
        dos_dict = {}
        enr = (self.eband.min(), self.eband.max())
        if self.data.is_spin_polarized:
            h = sum(np.array_split(self.accepted, 2)[0])
            eband_ud = np.array_split(self.eband, [h], axis=0)
            vvband_ud = np.array_split(self.vvband, [h], axis=0)
            spins = [Spin.up, Spin.down]
        else:
            eband_ud = [self.eband]
            vvband_ud = [self.vvband]
            spins = [Spin.up]

        for spin, eb, vvb in zip(spins, eband_ud, vvband_ud):
            energies, densities, vvdos, cdos = BL.BTPDOS(eb, vvb, npts=npts_mu, erange=enr)

            if T:
                densities = BL.smoothen_DOS(energies, densities, T)

            dos_dict.setdefault(spin, densities)

        tdos = Dos(self.efermi / units.eV, energies / units.eV, dos_dict)

        if partial_dos:
            tdos = self.get_partial_doses(tdos, eband_ud, spins, enr, npts_mu, T, progress)

        return tdos

    def get_partial_doses(self, tdos, eband_ud, spins, enr, npts_mu, T, progress):
        """Return a CompleteDos object interpolating the projections.

        tdos: total dos previously calculated
        npts_mu: number of energy points of the Dos
        T: parameter used to smooth the Dos
        progress: Default False, If True a progress bar is shown.
        """
        if not self.data.proj:
            raise BoltztrapError("No projections loaded.")

        bkp_data_ebands = np.copy(self.data.ebands)

        pdoss = {}
        if progress:
            n_iter = np.prod(np.sum([np.array(i.shape)[2:] for i in self.data.proj.values()]))
            t = tqdm(total=n_iter * 2)
        for spin, eb in zip(spins, eband_ud):
            for idx, site in enumerate(self.data.structure):
                if site not in pdoss:
                    pdoss[site] = {}
                for iorb, orb in enumerate(Orbital):
                    if progress:
                        t.update()
                    if iorb == self.data.proj[spin].shape[-1]:
                        break

                    if orb not in pdoss[site]:
                        pdoss[site][orb] = {}

                    self.data.ebands = self.data.proj[spin][:, :, idx, iorb].T
                    coeffs = fite.fitde3D(self.data, self.equivalences)
                    proj, vvproj, cproj = fite.getBTPbands(self.equivalences, coeffs, self.data.lattvec)

                    edos, pdos = BL.DOS(eb, npts=npts_mu, weights=np.abs(proj.real), erange=enr)

                    if T:
                        pdos = BL.smoothen_DOS(edos, pdos, T)

                    pdoss[site][orb][spin] = pdos

        self.data.ebands = bkp_data_ebands

        return CompleteDos(self.data.structure, total_dos=tdos, pdoss=pdoss)


class BztTransportProperties:
    """Compute Seebeck, Conductivity, Electrical part of thermal conductivity
    and Hall coefficient, conductivity effective mass, Power Factor tensors
    w.r.t. the chemical potential and temperatures, from dft band structure via
    interpolation.
    """

    def __init__(
        self,
        BztInterpolator,
        temp_r=None,
        doping=None,
        npts_mu=4000,
        CRTA=1e-14,
        margin=None,
        save_bztTranspProps=False,
        load_bztTranspProps=False,
        fname="bztTranspProps.json.gz",
    ):
        """
        Args:
            BztInterpolator: a BztInterpolator previously generated
            temp_r: numpy array of temperatures at which to calculate transport properties
            doping: doping levels at which to calculate transport properties. If provided,
                transport properties w.r.t. these doping levels are also computed. See
                compute_properties_doping() method for details.
            npts_mu: number of energy points at which to calculate transport properties
            CRTA: constant value of the relaxation time
            margin: The energy range of the interpolation is extended by this value on both sides.
                Defaults to 9 * units.BOLTZMANN * temp_r.max().
            save_bztTranspProps: Default False. If True all computed transport properties
                will be stored in fname file.
            load_bztTranspProps: Default False. If True all computed transport properties
                will be loaded from fname file.
            fname: File path where to save/load transport properties.

        Upon creation, it contains properties tensors w.r.t. the chemical potential
        of size (len(temp_r),npts_mu,3,3):
            Conductivity_mu (S/m), Seebeck_mu (microV/K), Kappa_mu (W/(m*K)),
            Power_Factor_mu (milliW/K m);
            cond_Effective_mass_mu (m_e) calculated as Ref.
        Also:
            Carrier_conc_mu: carrier concentration of size (len(temp_r),npts_mu)
            Hall_carrier_conc_trace_mu: trace of Hall carrier concentration of size
                (len(temp_r),npts_mu)
            mu_r_eV: array of energies in eV and with E_fermi at 0.0
                where all the properties are calculated.

        Example:
            bztTransp = BztTransportProperties(bztInterp,temp_r = np.arange(100,1400,100))
        """
        if temp_r is None:
            temp_r = np.arange(100, 1400, 100)

        self.dosweight = BztInterpolator.data.dosweight
        self.volume = BztInterpolator.data.get_volume()
        self.nelect = BztInterpolator.data.nelect
        self.efermi = BztInterpolator.data.fermi / units.eV

        if margin is None:
            margin = 9 * units.BOLTZMANN * temp_r.max()

        if load_bztTranspProps:
            self.load(fname)
        else:
            self.CRTA = CRTA
            self.temp_r = temp_r
            self.doping = doping

            self.epsilon, self.dos, self.vvdos, self.cdos = BL.BTPDOS(
                BztInterpolator.eband,
                BztInterpolator.vvband,
                npts=npts_mu,
                cband=BztInterpolator.cband,
            )

            mur_indices = np.logical_and(
                self.epsilon > self.epsilon.min() + margin,
                self.epsilon < self.epsilon.max() - margin,
            )

            self.mu_r = self.epsilon[mur_indices]  # mu range
            self.mu_r_eV = self.mu_r / units.eV - self.efermi

            N, L0, L1, L2, Lm11 = BL.fermiintegrals(
                self.epsilon,
                self.dos,
                self.vvdos,
                mur=self.mu_r,
                Tr=temp_r,
                dosweight=self.dosweight,
                cdos=self.cdos,
            )

            # Compute the Onsager coefficients from those Fermi integrals
            self.Conductivity_mu, self.Seebeck_mu, self.Kappa_mu, Hall_mu = BL.calc_Onsager_coefficients(
                L0, L1, L2, self.mu_r, temp_r, self.volume, Lm11=Lm11
            )

            # Common properties rescaling
            self.Conductivity_mu *= CRTA  # S / m
            self.Seebeck_mu *= 1e6  # microvolt / K
            self.Kappa_mu *= CRTA  # W / (m K)
            self.Hall_carrier_conc_trace_mu = (
                units.Coulomb
                * 1e-6
                / (np.abs(Hall_mu[:, :, 0, 1, 2] + Hall_mu[:, :, 2, 0, 1] + Hall_mu[:, :, 1, 2, 0]) / 3)
            )
            self.Carrier_conc_mu = (N + self.nelect) / (self.volume / (units.Meter / 100.0) ** 3)

            # Derived properties
            cond_eff_mass = np.zeros((len(self.temp_r), len(self.mu_r), 3, 3))
            for t in range(len(self.temp_r)):
                for i in range(len(self.mu_r)):
                    try:
                        cond_eff_mass[t, i] = (
                            np.linalg.inv(self.Conductivity_mu[t, i])
                            * self.Carrier_conc_mu[t, i]
                            * units.qe_SI**2
                            / units.me_SI
                            * 1e6
                        )
                    except np.linalg.LinAlgError:
                        pass

            self.Effective_mass_mu = cond_eff_mass * CRTA

            self.Power_Factor_mu = (self.Seebeck_mu @ self.Seebeck_mu) @ self.Conductivity_mu
            self.Power_Factor_mu *= 1e-9  # milliWatt / m / K**2

            self.contain_props_doping = False

            if isinstance(doping, np.ndarray):
                self.compute_properties_doping(doping, temp_r)

            if save_bztTranspProps:
                self.save(fname)

    def compute_properties_doping(self, doping, temp_r=None):
        """Calculate all the properties w.r.t. the doping levels in input.

        Args:
            doping: numpy array specifying the doping levels
            temp_r: numpy array specifying the temperatures

        When executed, it add the following variable at the BztTransportProperties
        object:
            Conductivity_doping, Seebeck_doping, Kappa_doping, Power_Factor_doping,
            cond_Effective_mass_doping are dictionaries with 'n' and 'p' keys and
            arrays of dim (len(temp_r),len(doping),3,3) as values.
            Carriers_conc_doping: carriers concentration for each doping level and T.
            mu_doping_eV: the chemical potential corrispondent to each doping level.
        """
        if temp_r is None:
            temp_r = self.temp_r

        self.Conductivity_doping, self.Seebeck_doping, self.Kappa_doping, self.Carriers_conc_doping = {}, {}, {}, {}

        self.Power_Factor_doping, self.Effective_mass_doping = {}, {}

        mu_doping = {}
        doping_carriers = [dop * (self.volume / (units.Meter / 100.0) ** 3) for dop in doping]

        for dop_type in ["n", "p"]:
            sbk = np.zeros((len(temp_r), len(doping), 3, 3))
            cond = np.zeros((len(temp_r), len(doping), 3, 3))
            kappa = np.zeros((len(temp_r), len(doping), 3, 3))
            hall = np.zeros((len(temp_r), len(doping), 3, 3, 3))
            dc = np.zeros((len(temp_r), len(doping)))

            if dop_type == "p":
                doping_carriers = [-dop for dop in doping_carriers]

            mu_doping[dop_type] = np.zeros((len(temp_r), len(doping)))
            for idx_t, temp in enumerate(temp_r):
                for idx_d, dop_car in enumerate(doping_carriers):
                    mu_doping[dop_type][idx_t, idx_d] = BL.solve_for_mu(
                        self.epsilon, self.dos, self.nelect + dop_car, temp, self.dosweight, True, False  # noqa: FBT003
                    )

                N, L0, L1, L2, Lm11 = BL.fermiintegrals(
                    self.epsilon,
                    self.dos,
                    self.vvdos,
                    mur=mu_doping[dop_type][idx_t],
                    Tr=np.array([temp]),
                    dosweight=self.dosweight,
                )

                cond[idx_t], sbk[idx_t], kappa[idx_t], hall[idx_t] = BL.calc_Onsager_coefficients(
                    L0, L1, L2, mu_doping[dop_type][idx_t], np.array([temp]), self.volume, Lm11
                )

                dc[idx_t] = self.nelect + N

            self.Conductivity_doping[dop_type] = cond * self.CRTA  # S / m
            self.Seebeck_doping[dop_type] = sbk * 1e6  # microVolt / K
            self.Kappa_doping[dop_type] = kappa * self.CRTA  # W / (m K)
            # self.Hall_doping[dop_type] = hall
            self.Carriers_conc_doping[dop_type] = dc / (self.volume / (units.Meter / 100.0) ** 3)

            self.Power_Factor_doping[dop_type] = (sbk @ sbk) @ cond * self.CRTA * 1e3

            cond_eff_mass = np.zeros((len(temp_r), len(doping), 3, 3))
            for idx_t in range(len(temp_r)):
                for idx_d, dop in enumerate(doping):
                    try:
                        cond_eff_mass[idx_t, idx_d] = (
                            np.linalg.inv(cond[idx_t, idx_d]) * dop * units.qe_SI**2 / units.me_SI * 1e6
                        )
                    except np.linalg.LinAlgError:
                        pass

            self.Effective_mass_doping[dop_type] = cond_eff_mass

        self.doping = doping
        self.mu_doping = mu_doping
        self.mu_doping_eV = {k: v / units.eV - self.efermi for k, v in mu_doping.items()}
        self.contain_props_doping = True

    # def find_mu_doping(self, epsilon, dos, N0, T, dosweight=2.):
    #     """
    #     Find the mu.

    #     :param epsilon:
    #     :param dos:
    #     :param N0:
    #     :param T:
    #     :param dosweight:
    #     """
    #     delta = np.empty_like(epsilon)
    #     for i, e in enumerate(epsilon):
    #         delta[i] = BL.calc_N(epsilon, dos, e, T, dosweight) + N0
    #     delta = np.abs(delta)
    #     # Find the position optimizing this distance
    #     pos = np.abs(delta).argmin()
    #     return epsilon[pos]

    def save(self, fname="bztTranspProps.json.gz"):
        """Save the transport properties to fname file."""
        lst_props = [
            self.temp_r,
            self.CRTA,
            self.epsilon,
            self.dos,
            self.vvdos,
            self.cdos,
            self.mu_r,
            self.mu_r_eV,
            self.Conductivity_mu,
            self.Seebeck_mu,
            self.Kappa_mu,
            self.Carrier_conc_mu,
            self.Hall_carrier_conc_trace_mu,
            self.Power_Factor_mu,
            self.Effective_mass_mu,
        ]

        if self.contain_props_doping:
            props = [
                self.Conductivity_doping,
                self.Seebeck_doping,
                self.Kappa_doping,
                self.Power_Factor_doping,
                self.Effective_mass_doping,
                self.Carriers_conc_doping,
                self.doping,
                self.mu_doping,
                self.mu_doping_eV,
            ]
            lst_props.extend(props)
        dumpfn(lst_props, fname)

    def load(self, fname="bztTranspProps.json.gz"):
        """Load the transport properties from fname file."""
        d = loadfn(fname)
        (
            self.temp_r,
            self.CRTA,
            self.epsilon,
            self.dos,
            self.vvdos,
            self.cdos,
            self.mu_r,
            self.mu_r_eV,
            self.Conductivity_mu,
            self.Seebeck_mu,
            self.Kappa_mu,
            self.Carrier_conc_mu,
            self.Hall_carrier_conc_trace_mu,
            self.Power_Factor_mu,
            self.Effective_mass_mu,
        ) = d[:15]
        if len(d) > 15:
            (
                self.Conductivity_doping,
                self.Seebeck_doping,
                self.Kappa_doping,
                self.Power_Factor_doping,
                self.Effective_mass_doping,
                self.Carriers_conc_doping,
                self.doping,
                self.mu_doping,
                self.mu_doping_eV,
            ) = d[15:]
            self.contains_doping_props = True

        return True


class BztPlotter:
    """Plotter to plot transport properties, interpolated bands along some high
    symmetry k-path, and DOS.

    Example:
        bztPlotter = BztPlotter(bztTransp,bztInterp)
        fig = self.bztPlotter.plot_props('S', 'mu', 'temp', temps=[300, 500])
        fig.show()
    """

    def __init__(self, bzt_transP=None, bzt_interp=None):
        """:param bzt_transP:
        :param bzt_interp:
        """
        self.bzt_transP = bzt_transP
        self.bzt_interp = bzt_interp

    def plot_props(
        self,
        prop_y,
        prop_x,
        prop_z="temp",
        output="avg_eigs",
        dop_type="n",
        doping=None,
        temps=None,
        xlim=(-2, 2),
        ax: plt.Axes = None,
    ):
        """Function to plot the transport properties.

        Args:
            prop_y: property to plot among ("Conductivity","Seebeck","Kappa","Carrier_conc",
                "Hall_carrier_conc_trace"). Abbreviations are possible, like "S" for "Seebeck"
            prop_x: independent variable in the x-axis among ('mu','doping','temp')
            prop_z: third variable to plot multiple curves ('doping','temp')
            output: 'avg_eigs' to plot the average of the eigenvalues of the properties
                tensors; 'eigs' to plot the three eigenvalues of the properties
                tensors.
            dop_type: 'n' or 'p' to specify the doping type in plots that use doping
                levels as prop_x or prop_z
            doping: list of doping level to plot, useful to reduce the number of curves
                when prop_z='doping'
            temps: list of temperatures to plot, useful to reduce the number of curves
                when prop_z='temp'
            xlim: chemical potential range in eV, useful when prop_x='mu'
            ax: figure.axes where to plot. If None, a new figure is produced.

        Example:
            bztPlotter.plot_props('S','mu','temp',temps=[600,900,1200]).show()
            more example are provided in the notebook
            "How to use Boltztra2 interface.ipynb".
        """
        props = (
            "Conductivity",
            "Seebeck",
            "Kappa",
            "Effective_mass",
            "Power_Factor",
            "Carrier_conc",
            "Hall_carrier_conc_trace",
        )
        props_lbl = (
            "Conductivity",
            "Seebeck",
            "$K_{el}$",
            "Effective mass",
            "Power Factor",
            "Carrier concentration",
            "Hall carrier conc.",
        )
        props_unit = (
            r"$(\mathrm{S\,m^{-1}})$",
            r"($\mu$V/K)",
            r"$(W / (m \cdot K))$",
            r"$(m_e)$",
            r"$( mW / (m\cdot K^2)$",
            r"$(cm^{-3})$",
            r"$(cm^{-3})$",
        )

        props_short = [p[: len(prop_y)] for p in props]

        if prop_y not in props_short:
            raise BoltztrapError("prop_y not valid")

        if prop_x not in ("mu", "doping", "temp"):
            raise BoltztrapError("prop_x not valid")

        if prop_z not in ("doping", "temp"):
            raise BoltztrapError("prop_z not valid")

        idx_prop = props_short.index(prop_y)

        leg_title = ""

        mu = self.bzt_transP.mu_r_eV

        if prop_z == "doping" and prop_x == "temp":
            p_array = eval("self.bzt_transP." + props[idx_prop] + "_" + prop_z)
        else:
            p_array = eval("self.bzt_transP." + props[idx_prop] + "_" + prop_x)

        if ax is None:
            plt.figure(figsize=(10, 8))

        temps_all = self.bzt_transP.temp_r.tolist()
        if temps is None:
            temps = self.bzt_transP.temp_r.tolist()

        if isinstance(self.bzt_transP.doping, np.ndarray):
            doping_all = self.bzt_transP.doping.tolist()
            if doping is None:
                doping = doping_all

        # special case of carrier and hall carrier concentration 2d arrays (temp,mu)
        if idx_prop in [5, 6]:
            if prop_z == "temp" and prop_x == "mu":
                for temp in temps:
                    ti = temps_all.index(temp)
                    prop_out = p_array[ti] if idx_prop == 6 else np.abs(p_array[ti])
                    plt.semilogy(mu, prop_out, label=f"{temp} K")

                plt.xlabel(r"$\mu$ (eV)", fontsize=30)
                plt.xlim(xlim)
            else:
                raise BoltztrapError(
                    "only prop_x=mu and prop_z=temp are \
                    available for c.c. and Hall c.c.!"
                )

        elif prop_z == "temp" and prop_x == "mu":
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[ti])[0]
                if output == "avg_eigs":
                    plt.plot(mu, prop_out.mean(axis=1), label=f"{temp} K")
                elif output == "eigs":
                    for i in range(3):
                        plt.plot(
                            mu,
                            prop_out[:, i],
                            label=f"eig {i} {temp} K",
                        )

            plt.xlabel(r"$\mu$ (eV)", fontsize=30)
            plt.xlim(xlim)

        elif prop_z == "temp" and prop_x == "doping":
            for temp in temps:
                ti = temps_all.index(temp)
                prop_out = np.linalg.eigh(p_array[dop_type][ti])[0]
                if output == "avg_eigs":
                    plt.semilogx(doping_all, prop_out.mean(axis=1), "s-", label=f"{temp} K")
                elif output == "eigs":
                    for i in range(3):
                        plt.plot(
                            doping_all,
                            prop_out[:, i],
                            "s-",
                            label=f"eig {i} {temp} K",
                        )
            plt.xlabel(r"Carrier conc. $cm^{-3}$", fontsize=30)
            leg_title = dop_type + "-type"

        elif prop_z == "doping" and prop_x == "temp":
            for dop in doping:
                di = doping_all.index(dop)
                prop_out = np.linalg.eigh(p_array[dop_type][:, di])[0]
                if output == "avg_eigs":
                    plt.plot(
                        temps_all,
                        prop_out.mean(axis=1),
                        "s-",
                        label=str(dop) + " $cm^{-3}$",
                    )
                elif output == "eigs":
                    for i in range(3):
                        plt.plot(
                            temps_all,
                            prop_out[:, i],
                            "s-",
                            label=f"eig {i} {dop} $cm^{{-3}}$",
                        )

            plt.xlabel(r"Temperature (K)", fontsize=30)
            leg_title = f"{dop_type}-type"

        plt.ylabel(f"{props_lbl[idx_prop]} {props_unit[idx_prop]}", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend(title=leg_title if leg_title != "" else "", fontsize=15)
        plt.tight_layout()
        plt.grid()
        return plt

    def plot_bands(self):
        """Plot a band structure on symmetry line using BSPlotter()."""
        if self.bzt_interp is None:
            raise BoltztrapError("BztInterpolator not present")

        sbs = self.bzt_interp.get_band_structure()

        return BSPlotter(sbs).get_plot()

    def plot_dos(self, T=None, npoints=10000):
        """Plot the total Dos using DosPlotter()."""
        if self.bzt_interp is None:
            raise BoltztrapError("BztInterpolator not present")

        tdos = self.bzt_interp.get_dos(T=T, npts_mu=npoints)
        dosPlotter = DosPlotter()
        dosPlotter.add_dos("Total", tdos)

        return dosPlotter


def merge_up_down_doses(dos_up, dos_dn):
    """Merge the up and down DOSs.

    Args:
    dos_up: Up DOS.
    dos_dn: Down DOS
    Return:
    CompleteDos object
    """
    warnings.warn(
        "This function is not useful anymore. VasprunBSLoader deals \
                   with spin case."
    )
    cdos = Dos(
        dos_up.efermi,
        dos_up.energies,
        {Spin.up: dos_up.densities[Spin.up], Spin.down: dos_dn.densities[Spin.down]},
    )

    if hasattr(dos_up, "pdos") and hasattr(dos_dn, "pdos"):
        pdoss = {}
        for site in dos_up.pdos:
            pdoss.setdefault(site, {})
            for orb in dos_up.pdos[site]:
                pdoss[site].setdefault(orb, {})
                pdoss[site][orb][Spin.up] = dos_up.pdos[site][orb][Spin.up]
                pdoss[site][orb][Spin.down] = dos_dn.pdos[site][orb][Spin.down]

        cdos = CompleteDos(dos_up.structure, total_dos=cdos, pdoss=pdoss)

    return cdos
