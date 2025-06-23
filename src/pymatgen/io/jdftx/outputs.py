"""JDFTx outputs parsing module.

Module for parsing outputs of JDFTx.

Note: JDFTXOutfile will be moved back to its own module once a more broad outputs
class is written.
"""

from __future__ import annotations

import pprint
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
from monty.dev import deprecated

from pymatgen.core import Lattice
from pymatgen.core.units import Ha_to_eV
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.jdftx._output_utils import (
    _get_nbands_from_bandfile_filepath,
    _get_orb_label_list,
    _get_u_to_oa_map,
    _parse_kptsfrom_bandprojections_file,
    get_proj_tju_from_file,
    orb_ref_list,
    read_outfile_slices,
)
from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from pymatgen.core.structure import Structure
    from pymatgen.core.trajectory import Trajectory
    from pymatgen.io.jdftx.inputs import JDFTXInfile
    from pymatgen.io.jdftx.jelstep import JElSteps
    from pymatgen.io.jdftx.jminsettings import (
        JMinSettingsElectronic,
        JMinSettingsFluid,
        JMinSettingsIonic,
        JMinSettingsLattice,
    )
    from pymatgen.io.jdftx.joutstructures import JOutStructures

__author__ = "Ben Rich, Jacob Clary"


# This will be a long list so keep alphabetical to make adding more files easier.
# Dump files that are redundant to data contained in the outfile attribute are included for completeness
# but commented out.
dump_file_names = (
    "bandProjections",
    "dosUp",
    "dosDn",
    "d_tot",
    "eigenvals",
    # "eigStats"
    # "Ecomponents",
    "fillings",
    # "force",
    "Gvectors",
    # "ionpos",
    # "lattice",
    "kPtsn",
    "n_up",
    "n_dn",
    "tau",
    "tau_up",
    "tau_dn",
    "wfns",
)

implemented_store_vars = (
    "bandProjections",
    "eigenvals",
    "kpts",
    "bandstructure",
)


@dataclass
class JDFTXOutputs:
    """JDFTx outputs parsing class.

    A class to read and process JDFTx outputs.

    Methods:
        from_calc_dir(calc_dir: str | Path, store_vars: list[str] = None) -> JDFTXOutputs:
            Return JDFTXOutputs object from the path to a directory containing JDFTx out files.

    Attributes:
        calc_dir (str | Path): The path to the directory containing the JDFTx output files.
        store_vars (list[str]): A list of the names of dump files to read and store.
        paths (dict[str, Path]): A dictionary of the paths to the dump files.
        outfile (JDFTXOutfile): The JDFTXOutfile object for the out file.
        bandProjections (np.ndarray): The band projections. Stored in shape (nstates, nbands, nproj) where nstates
            is nspin*nkpts (nkpts may not equal prod(kfolding) if symmetry reduction occurred), nbands is the number of
            bands, and nproj is the number of projections. This shape is chosen instead of the pymatgen convention of
            (nspin, nkpt, nbands, nion, nionproj) to save on memory as nonionproj is different depending on the ion
            type. This array may also be complex if specified in 'band-projections-params' in the JDFTx input, allowing
            for pCOHP analysis.
        eigenvals (np.ndarray): The eigenvalues. Stored in shape (nstates, nbands) where nstates is nspin*nkpts (nkpts
            may not equal prod(kfolding) if symmetry reduction occurred) and nbands is the number of bands.
        orb_label_list (tuple[str, ...]): A tuple of the orbital labels for the bandProjections file, where the i'th
            element describes the i'th orbital. Orbital labels are formatted as "<ion>#<ion-number>(<orbital>)",
            where <ion> is the element symbol of the ion, <ion-number> is the 1-based index of the ion-type in the
            structure (ie C#2 would be the second carbon atom, but not necessarily the second ion in the structure),
            and <orbital> is a string describing "l" and "ml" quantum numbers (ie "p_x" or "d_yz"). Note that while "z"
            corresponds to the "z" axis, "x" and "y" are arbitrary and may not correspond to the actual x and y axes of
            the structure. In the case where multiple shells of a given "l" are available within the projections, a
            0-based index will appear mimicking a principle quantum number (ie "0px" for first shell and "1px" for
            second shell). The actual principal quantum number is not stored in the JDFTx output files and must be
            inferred by the user.
        kpts (list[np.ndarray]): A list of the k-points used in the calculation. Each k-point is a 3D numpy array.
        wk_list (list[np.ndarray]): A list of the weights for the k-points used in the calculation.
    """

    calc_dir: str | Path = field(init=True)
    outfile_name: str | Path | None = field(init=True)
    store_vars: list[str] = field(default_factory=list, init=True)
    paths: dict[str, Path] = field(init=False)
    outfile: JDFTXOutfile = field(init=False)
    bandProjections: NDArray | None = field(init=False)
    eigenvals: NDArray | None = field(init=False)
    kpts: list[NDArray] | None = field(init=False)
    wk_list: list[NDArray] | None = field(init=False)
    # Misc metadata for interacting with the data
    orb_label_list: tuple[str, ...] | None = field(init=False)
    bandstructure: BandStructure | None = field(init=False)

    @classmethod
    def from_calc_dir(
        cls, calc_dir: str | Path, store_vars: list[str] | None = None, outfile_name: str | Path | None = None
    ) -> JDFTXOutputs:
        """
        Create a JDFTXOutputs object from a directory containing JDFTx out files.

        Args:
            calc_dir (str | Path): The path to the directory containing the JDFTx out files.
            is_bgw (bool): Mark True if data must be usable for BGW calculations. This will change the behavior of the
                parser to be stricter with certain criteria.
            none_slice_on_error (bool): If True, will return None if an error occurs while parsing a slice instead of
                halting the parsing process. This can be useful for parsing files with multiple slices where some slices
                may be incomplete or corrupted.
            outfile_name (str | Path): The name of the outfile to use. If None, will search for the outfile in the
                calc_dir. If provided, will concatenate with calc_dir as the outfile path. Use this if the calc_dir
                contains multiple files that may be mistaken for the outfile (ie multiple files with the '.out' suffix).
        Returns:
            JDFTXOutputs: The JDFTXOutputs object.
        """
        store_vars = cls._check_store_vars(store_vars)
        return cls(calc_dir=Path(calc_dir), store_vars=store_vars, outfile_name=outfile_name)

    def __post_init__(self):
        self._init_paths()
        self._store_vars()
        self._init_bandstructure()

    def _check_store_vars(store_vars: list[str] | None) -> list[str]:
        if store_vars is None:
            return []
        if "bandstructure" in store_vars:
            store_vars += ["kpts", "eigenvals"]
            store_vars.pop(store_vars.index("bandstructure"))
        return list(set(store_vars))

    def _init_paths(self):
        self.paths = {}
        if self.calc_dir is None:
            raise ValueError("calc_dir must be set as not None before initializing.")
        if self.outfile_name is None:
            outfile_path = _find_jdftx_out_file(self.calc_dir)
        else:
            outfile_path = self.calc_dir / self.outfile_name
            if not outfile_path.exists():
                raise FileNotFoundError(f"Provided outfile path {outfile_path} does not exist.")
        self.outfile = JDFTXOutfile.from_file(outfile_path)
        prefix = self.outfile.prefix
        for fname in dump_file_names:
            try:
                paths = _find_jdftx_dump_file(self.calc_dir, fname)
                path = _disambiguate_paths(paths, fname, prefix)
                self.paths[fname] = path
            except FileNotFoundError:
                pass

    def _store_vars(self):
        for var in implemented_store_vars:
            if var not in self.store_vars:
                setattr(self, var, None)
        for var in self.store_vars:
            self._store_var(var)

    def _store_var(self, var: str):
        if not hasattr(self, f"_store_{var}"):
            raise NotImplementedError(f"Storing {var} is not currently implemented.")
        if hasattr(self, f"_check_{var}"):
            check_method = getattr(self, f"_check_{var}")
            check_method()
        init_method = getattr(self, f"_store_{var}")
        init_method()

    def _check_bandProjections(self):
        """Check for misaligned data within bandProjections file."""
        if "bandProjections" in self.paths:
            if not self.paths["bandProjections"].exists():
                raise RuntimeError("Allocated path for bandProjections does not exist.")
            nbands = _get_nbands_from_bandfile_filepath(self.paths["bandProjections"])
            if not nbands == self.outfile.nbands:
                raise ValueError("Number of bands in bandProjections file does not match number of bands in outfile.")

    def _store_bandProjections(self):
        if "bandProjections" in self.paths:
            self.bandProjections = get_proj_tju_from_file(self.paths["bandProjections"])
            self.orb_label_list = _get_orb_label_list(self.paths["bandProjections"])

    def _check_eigenvals(self):
        """Check for misaligned data within eigenvals file."""
        if "eigenvals" in self.paths:
            if not self.paths["eigenvals"].exists():
                raise RuntimeError("Allocated path for eigenvals does not exist.")
            tj = len(np.fromfile(self.paths["eigenvals"]))
            nstates_float = tj / self.outfile.nbands
            if not np.isclose(nstates_float, int(nstates_float)):
                raise ValueError("Number of eigenvalues is not an integer multiple of number of bands.")

    def _store_eigenvals(self):
        if "eigenvals" in self.paths:
            self.eigenvals = np.fromfile(self.paths["eigenvals"])
            nstates = int(len(self.eigenvals) / self.outfile.nbands)
            self.eigenvals = self.eigenvals.reshape(nstates, self.outfile.nbands) * Ha_to_eV

    def _check_kpts(self):
        if "bandProjections" in self.paths and self.paths["bandProjections"].exists():
            # TODO: Write kpt file inconsistency checking
            return
        if "kPts" in self.paths and self.paths["kPts"].exists():
            raise NotImplementedError("kPts file parsing not yet implemented.")
        raise RuntimeError("No k-point data found in JDFTx output files.")

    def _store_kpts(self):
        if "bandProjections" in self.paths and self.paths["bandProjections"].exists():
            wk_list, kpts_list = _parse_kptsfrom_bandprojections_file(self.paths["bandProjections"])
            self.kpts = kpts_list
            self.wk_list = wk_list

    def _init_bandstructure(self):
        if True in [v is None for v in [self.kpts, self.eigenvals]]:
            return
        kpoints = np.array(self.kpts)
        eigenvals = self._get_pmg_eigenvals()
        projections = None
        if self.bandProjections is not None:
            projections = self._get_pmg_projections()
        self.bandstructure = BandStructure(
            kpoints,
            eigenvals,
            Lattice(self.outfile.structure.lattice.reciprocal_lattice.matrix),
            self.outfile.mu,
            projections=projections,
            structure=self.outfile.structure,
        )

    def _get_lmax(self) -> tuple[str | None, int | None]:
        """Get the maximum l quantum number and projection array orbital length.

        Returns:
            tuple[str | None, int | None]: The maximum l quantum number and projection array orbital length.
                Both are None if no bandProjections file is available.
        """
        if self.orb_label_list is None:
            return None, None
        orbs = [label.split("(")[1].split(")")[0] for label in self.orb_label_list]
        if orb_ref_list[-1][0] in orbs:
            return "f", 16
        if orb_ref_list[-2][0] in orbs:
            return "d", 9
        if orb_ref_list[-3][0] in orbs:
            return "p", 4
        if orb_ref_list[-4][0] in orbs:
            return "s", 1
        raise ValueError("Unrecognized orbital labels in orb_label_list.")

    def _get_pmg_eigenvals(self) -> dict | None:
        if self.eigenvals is None:
            return None
        _e_skj = self.eigenvals.copy().reshape(self.outfile.nspin, -1, self.outfile.nbands)
        _e_sjk = np.swapaxes(_e_skj, 1, 2)
        spins = [Spin.up, Spin.down]
        eigenvals = {}
        for i in range(self.outfile.nspin):
            eigenvals[spins[i]] = _e_sjk[i]
        return eigenvals

    def _get_pmg_projections(self) -> dict | None:
        """Return pymatgen-compatible projections dictionary.

        Converts the bandProjections array to a pymatgen-compatible dictionary

        Returns:
            dict | None:
        """
        lmax, norbmax = self._get_lmax()
        if norbmax is None:
            return None
        if self.orb_label_list is None:
            return None
        _proj_tju = self.bandProjections.copy()
        if _proj_tju.dtype is np.complex64:
            # Convert <orb|band(state)> to |<orb|band(state)>|^2
            _proj_tju = np.abs(_proj_tju) ** 2
        # Convert to standard datatype - using np.real to suppress warnings
        proj_tju = np.array(np.real(_proj_tju), dtype=float)
        proj_skju = proj_tju.reshape([self.outfile.nspin, -1, self.outfile.nbands, len(self.orb_label_list)])
        proj_sjku = np.swapaxes(proj_skju, 1, 2)
        nspin, nbands, nkpt, nproj = proj_sjku.shape
        u_to_oa_map = _get_u_to_oa_map(self.paths["bandProjections"])
        projections = {}
        spins = [Spin.up, Spin.down]
        for i in range(nspin):
            projections[spins[i]] = np.zeros([nbands, nkpt, norbmax, len(self.outfile.structure)])
            # TODO: Consider jitting this loop
            for u in range(nproj):
                projections[spins[i]][:, :, *u_to_oa_map[u]] += proj_sjku[i, :, :, u]
        return projections


_jof_atr_from_last_slice = (
    "prefix",
    "jstrucs",
    "jsettings_fluid",
    "jsettings_electronic",
    "jsettings_lattice",
    "jsettings_ionic",
    "xc_func",
    "lattice_initial",
    "lattice_final",
    "lattice",
    "a",
    "b",
    "c",
    "fftgrid",
    "geom_opt",
    "geom_opt_type",
    "efermi",
    "egap",
    "emin",
    "emax",
    "homo",
    "lumo",
    "homo_filling",
    "lumo_filling",
    "is_metal",
    "converged",
    "etype",
    "broadening_type",
    "broadening",
    "kgrid",
    "truncation_type",
    "truncation_radius",
    "pwcut",
    "rhocut",
    "pp_type",
    "total_electrons",
    "semicore_electrons",
    "valence_electrons",
    "total_electrons_uncharged",
    "semicore_electrons_uncharged",
    "valence_electrons_uncharged",
    "nbands",
    "atom_elements",
    "atom_elements_int",
    "atom_types",
    "spintype",
    "nspin",
    "nat",
    "atom_coords_initial",
    "atom_coords_final",
    "atom_coords",
    "structure",
    "has_solvation",
    "fluid",
    "is_gc",
    "eopt_type",
    "elecmindata",
    "stress",
    "strain",
    "forces",
    "nstep",
    "e",
    "grad_k",
    "alpha",
    "linmin",
    "abs_magneticmoment",
    "tot_magneticmoment",
    "mu",
    "elec_nstep",
    "elec_e",
    "elec_grad_k",
    "elec_alpha",
    "elec_linmin",
    "electronic_output",
    "t_s",
    "ecomponents",
    "infile",
)

# TODO: Remove references to the deprecated 'jsettings_*' attributes in `JDFTXOutfile` and `JDFTXOutfileSlice`
# TODO: JDFTXOutfile and JDFTXOutfileSlice are VERY bloated. The following attributes are redundant
# to fields of the `Structure` object and their removal likely won't cause any confusion.
# - lattice
# - lattice_final
# - atom_coords
# - atom_coords_final
# - a/b/c
# The following attributes are redundant to the internal `JDFTXInfile` object, but since the JDFTXInfile
# object is not a standard pymatgen object, it may be better to decorate them with a deprecated decorator
# until end of year. Additionally, it should be made certain which of these are dumped in the outfile's
# input section regardless of what was set in the input file - for those that are not, they need to be kept.
# - lattice_initial
# - atom_coords_initial
# - broadening
# - broadening_type
# - kgrid
# - truncation_type
# - truncation_radius
# - fluid
# - pwcut
# - rhocut
# The following are attributes that come from the electronic/fluid/ionic/lattice optimization step
# logs. I am not sure if it is actually helpful to have access to these from the `JDFTXOutfile` object,
# as they are only informative for analyzing optimization convergence. Additionally, the fields are
# less universal than I thought, and can change depending on the optimization settings, so it might
# be smarter to store them in a dictionary of arbitrary keys instead within the contained substructures.
# - grad_k
# - alpha
# - linmin
# - elec_grad_k
# - elec_alpha
# - elec_linmin


@dataclass
class JDFTXOutfile:
    """
    JDFTx out file parsing class.

    A class to read and process a JDFTx out file.

    Methods:
        from_file(file_path: str | Path) -> JDFTXOutfile:
            Return JDFTXOutfile object from the path to a JDFTx out file.

    Attributes:
        slices (list[JDFTXOutfileSlice]): A list of JDFTXOutfileSlice objects. Each slice corresponds to an individual
            call of the JDFTx executable. Subsequent JDFTx calls within the same directory and prefix will append
            outputs to the same out file. More than one slice may correspond to restarted calculations, geom + single
            point calculations, or optimizations done with 3rd-party wrappers like ASE.
        prefix (str): The prefix of the most recent JDFTx call.
        jstrucs (JOutStructures): The JOutStructures object from the most recent JDFTx call. This object contains a
            series of JOutStructure objects in its 'slices' attribute, each corresponding to a single structure
            (multiple iff performing a geometric optimization) as well as convergence data for the structures as a
            series.
        jsettings_fluid (JMinSettingsFluid): The JMinSettingsFluid object from the most recent JDFTx call. This object
            contains only a 'params' attribute, which is a dictionary of the input parameters for the fluid
            optimization.
        jsettings_electronic (JMinSettingsElectronic): The JMinSettingsElectronic object from the most recent JDFTx
            call. This object contains only a 'params' attribute, which is a dictionary of the input parameters for the
            electronic optimization.
        jsettings_lattice (JMinSettingsLattice): The JMinSettingsLattice object from the most recent JDFTx call. This
            object contains only a 'params' attribute, which is a dictionary of the input parameters for the lattice
            optimization.
        jsettings_ionic (JMinSettingsIonic): The JMinSettingsIonic object from the most recent JDFTx call. This object
            contains only a 'params' attribute, which is a dictionary of the input parameters for the ionic
            optimization.
        xc_func (str): The exchange-correlation functional used in the most recent JDFTx call. See documentation for
            JDFTx online for a list of available exchange-correlation functionals.
        lattice_initial (np.ndarray): The initial lattice vectors of the most recent JDFTx call as a 3x3 numpy array.
            In units of Angstroms.
        lattice_final (np.ndarray): The final lattice vectors of the most recent JDFTx call as a 3x3 numpy array. In
            units of Angstroms.
        lattice (np.ndarray): The lattice vectors of the most recent JDFTx call as a 3x3 numpy array (redundant to
            lattice_final).
        a (float): Length of the first lattice vector. In units of Angstroms.
        b (float): Length of the second lattice vector. In units of Angstroms.
        c (float): Length of the third lattice vector. In units of Angstroms.
        fftgrid (list[int]): The FFT grid shape used in the most recent JDFTx call. Can be used to properly shape
            densities dumped as binary files.
        geom_opt (bool): True if the most recent JDFTx call was a geometry optimization (lattice or ionic).
        geom_opt_type (str): The type of geometry optimization performed in the most recent JDFTx call. Options are
            'lattice' or 'ionic' if geom_opt, else "single point". ('lattice' optimizations perform ionic optimizations
            as well unless ion positions are given in direct coordinates).
        ecomponents (dict): The components of the total energy in eV of the most recent JDFTx call.
        efermi (float): The Fermi energy in eV of the most recent JDFTx call. Equivalent to "mu".
        egap (float): The band gap in eV of the most recent JDFTx call. (Only available if eigstats was dumped).
        emin (float): The minimum energy in eV (smallest Kohn-Sham eigenvalue) of the most recent JDFTx call. (Only
            available if eigstats was dumped).
        emax (float): The maximum energy in eV (largest Kohn-Sham eigenvalue) of the most recent JDFTx call. (Only
            available if eigstats was dumped).
        homo (float): The energy in eV of the band-gap lower bound (Highest Occupied Molecular Orbital) (Only available
            if eigstats was dumped).
        lumo (float): The energy in eV of the band-gap upper bound (Lowest Unoccupied Molecular Orbital) (Only
            available if eigstats was dumped).
        homo_filling (float): The electron filling at the homo band-state. (Only available if eigstats was dumped).
        lumo_filling (float): The electron filling at the lumo band-state. (Only available if eigstats was dumped).
        is_metal (bool): True if fillings of homo and lumo band-states are off-set by 1 and 0 by at least an arbitrary
            tolerance of 0.01 (ie 1 - 0.015 and 0.012 for homo/lumo fillings would be metallic, while 1-0.001 and 0
            would not be). (Only available if eigstats was dumped).
        converged (bool): True if most recent SCF cycle converged (and geom forces converged is calc is geom_opt)
        etype (str): String representation of total energy-type of system. Commonly "G" (grand-canonical potential) for
            GC calculations, and "F" for canonical (fixed electron count) calculations.
        broadening_type (str): Type of broadening for electronic filling about Fermi-level requested. Either "Fermi",
            "Cold", "MP1", or "Gauss".
        broadening (float): Magnitude of broadening for electronic filling.
        kgrid (list[int]): Shape of k-point grid used in calculation. (equivalent to k-point folding)
        truncation_type (str): Type of coulomb truncation used to prevent interaction between periodic images along
            certain directions. "periodic" means no coulomb truncation was used.
        truncation_radius (float | None): If spherical truncation_type, this is the radius of the coulomb truncation
            sphere.
        pwcut (float): The plane-wave cutoff energy in Hartrees used in the most recent JDFTx call.
        rhocut (float): The density cutoff energy in Hartrees used in the most recent JDFTx call.
        pp_type (str): The pseudopotential library used in the most recent JDFTx call. Currently only "GBRV" and "SG15"
            are supported by this output parser.
        total_electrons (float): The total number of electrons in the most recent JDFTx call (redundant to nelectrons).
        semicore_electrons (int): The number of semicore electrons in the most recent JDFTx call.
        valence_electrons (float): The number of valence electrons in the most recent JDFTx call.
        total_electrons_uncharged (int): The total number of electrons in the most recent JDFTx call, uncorrected for
            charge. (ie total_electrons + charge)
        semicore_electrons_uncharged (int): The number of semicore electrons in the most recent JDFTx call, uncorrected
            for charge. (ie semicore_electrons + charge)
        valence_electrons_uncharged (int): The number of valence electrons in the most recent JDFTx call, uncorrected
            for charge. (ie valence_electrons + charge)
        nbands (int): The number of bands used in the most recent JDFTx call.
        atom_elements (list[str]): The list of each ion's element symbol in the most recent JDFTx call.
        atom_elements_int (list[int]): The list of ion's atomic numbers in the most recent JDFTx call.
        atom_types (list[str]): Non-repeating list of each ion's element symbol in the most recent JDFTx call.
        spintype (str): The spin type used in the most recent JDFTx call. Options are "none", "collinear",
        nspin (int): The number of spins used in the most recent JDFTx call.
        nat (int): The number of atoms in the most recent JDFTx call.
        atom_coords_initial (list[list[float]]): The initial atomic coordinates of the most recent JDFTx call.
        atom_coords_final (list[list[float]]): The final atomic coordinates of the most recent JDFTx call.
        atom_coords (list[list[float]]): The atomic coordinates of the most recent JDFTx call.
        structure (Structure): The updated pymatgen Structure object of the most recent JDFTx call.
        trajectory (Trajectory): The Trajectory object of the most recent JDFTx call.
        has_solvation (bool): True if the most recent JDFTx call included a solvation calculation.
        fluid (str): The fluid used in the most recent JDFTx call.
        is_gc (bool): True if the most recent slice is a grand canonical calculation.
        eopt_type (str): The type of energy iteration used in the most recent JDFTx call.
        elecmindata (JElSteps): The JElSteps object from the most recent JDFTx call. This object contains a series of
            JElStep objects in its 'steps' attribute, each corresponding to a single energy iteration.
        stress (np.ndarray): The stress tensor of the most recent JDFTx call as a 3x3 numpy array. In units of
            eV/Angstrom^3.
        strain (np.ndarray): The strain tensor of the most recent JDFTx call as a 3x3 numpy array.
        nstep (int): The number of geometric optimization steps in the most recent JDFTx call.
        e (float): The final energy in eV of the most recent JDFTx call (equivalent to the call's etype).
        grad_k (float): The final norm of the preconditioned gradient for geometric optimization of the most recent
            JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
            (written as "|grad|_K" in JDFTx output).
        alpha (float): The step size of the final geometric step in the most recent JDFTx call.
        linmin (float): The final normalized projection of the geometric step direction onto the gradient for the most
            recent JDFTx call.
        abs_magneticmoment (float | None): The absolute magnetic moment of the most recent JDFTx call.
        tot_magneticmoment (float | None): The total magnetic moment of the most recent JDFTx call.
        mu (float): The Fermi energy in eV of the most recent JDFTx call.
        elec_e (float): The final energy in eV of the most recent electronic optimization step.
        elec_nstep (int): The number of electronic optimization steps in the most recent JDFTx call.
        elec_grad_k (float): The final norm of the preconditioned gradient for electronic optimization of the most
            recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
            (written as "|grad|_K" in JDFTx output).
        elec_alpha (float): The step size of the final electronic step in the most recent JDFTx call.
        elec_linmin (float): The final normalized projection of the electronic step direction onto the gradient for the
            most recent JDFTx call.
        infile (JDFTXInfile): The JDFTXInfile object representing the input parameters of the JDFTXOutfile.

    Magic Methods:
        __getitem__(key: str | int) -> Any: Decides behavior of how JDFTXOutfile objects are indexed. If the key is a
            string, it will return the value of the property with the same name. If the key is an integer, it will
            return the slice of the JDFTXOutfile object at that index.
        __len__() -> int: Returns the number of slices in the JDFTXOutfile object.
        __getattr__(name: str) -> Any: Returns the value of the property with the same name as the input string.
        __str__() -> str: Returns a string representation of the JDFTXOutfile object.
    """

    slices: list[JDFTXOutfileSlice | None] = field(default_factory=list)
    prefix: str = field(init=False)
    jstrucs: JOutStructures = field(init=False)
    jsettings_fluid: JMinSettingsFluid = field(init=False)
    jsettings_electronic: JMinSettingsElectronic = field(init=False)
    jsettings_lattice: JMinSettingsLattice = field(init=False)
    jsettings_ionic: JMinSettingsIonic = field(init=False)
    xc_func: str = field(init=False)
    lattice_initial: np.ndarray = field(init=False)
    lattice_final: np.ndarray = field(init=False)
    lattice: np.ndarray = field(init=False)
    a: float = field(init=False)
    b: float = field(init=False)
    c: float = field(init=False)
    fftgrid: list[int] = field(init=False)
    geom_opt: bool = field(init=False)
    geom_opt_type: str = field(init=False)
    efermi: float = field(init=False)
    egap: float = field(init=False)
    emin: float = field(init=False)
    emax: float = field(init=False)
    homo: float = field(init=False)
    lumo: float = field(init=False)
    homo_filling: float = field(init=False)
    lumo_filling: float = field(init=False)
    is_metal: bool = field(init=False)
    converged: bool = field(init=False)
    etype: str = field(init=False)
    broadening_type: str = field(init=False)
    broadening: float = field(init=False)
    kgrid: list[int] = field(init=False)
    truncation_type: str = field(init=False)
    truncation_radius: float = field(init=False)
    pwcut: float = field(init=False)
    rhocut: float = field(init=False)
    pp_type: str = field(init=False)
    total_electrons: float = field(init=False)
    semicore_electrons: int = field(init=False)
    valence_electrons: float = field(init=False)
    total_electrons_uncharged: int = field(init=False)
    semicore_electrons_uncharged: int = field(init=False)
    valence_electrons_uncharged: int = field(init=False)
    nbands: int = field(init=False)
    atom_elements: list[str] = field(init=False)
    atom_elements_int: list[int] = field(init=False)
    atom_types: list[str] = field(init=False)
    spintype: str = field(init=False)
    nspin: int = field(init=False)
    nat: int = field(init=False)
    atom_coords_initial: list[list[float]] = field(init=False)
    atom_coords_final: list[list[float]] = field(init=False)
    atom_coords: list[list[float]] = field(init=False)
    structure: Structure = field(init=False)
    trajectory: Trajectory = field(init=False)
    has_solvation: bool = field(init=False)
    fluid: str = field(init=False)
    is_gc: bool = field(init=False)
    eopt_type: str = field(init=False)
    elecmindata: JElSteps = field(init=False)
    stress: np.ndarray = field(init=False)
    strain: np.ndarray = field(init=False)
    forces: np.ndarray = field(init=False)
    nstep: int = field(init=False)
    e: float = field(init=False)
    grad_k: float = field(init=False)
    alpha: float = field(init=False)
    linmin: float = field(init=False)
    abs_magneticmoment: float = field(init=False)
    tot_magneticmoment: float = field(init=False)
    mu: float = field(init=False)
    elec_nstep: int = field(init=False)
    elec_e: float = field(init=False)
    elec_grad_k: float = field(init=False)
    elec_alpha: float = field(init=False)
    elec_linmin: float = field(init=False)
    electronic_output: float = field(init=False)
    infile: JDFTXInfile = field(init=False)

    @classmethod
    def from_calc_dir(
        cls, calc_dir: str | Path, is_bgw: bool = False, none_slice_on_error: bool | None = None
    ) -> JDFTXOutfile:
        """
        Create a JDFTXOutfile object from a directory containing JDFTx out files.

        Args:
            calc_dir (str | Path): The path to the directory containing the JDFTx out files.
            is_bgw (bool): Mark True if data must be usable for BGW calculations. This will change the behavior of the
                parser to be stricter with certain criteria.
            none_slice_on_error (bool): If True, will return None if an error occurs while parsing a slice instead of
                halting the parsing process. This can be useful for parsing files with multiple slices where some slices
                may be incomplete or corrupted.

        Returns:
            JDFTXOutfile: The JDFTXOutfile object.
        """
        file_path = _find_jdftx_out_file(Path(calc_dir))
        return cls.from_file(file_path=file_path, is_bgw=is_bgw, none_slice_on_error=none_slice_on_error)

    @classmethod
    def from_file(
        cls, file_path: str | Path, is_bgw: bool = False, none_slice_on_error: bool | None = None
    ) -> JDFTXOutfile:
        """
        Create a JDFTXOutfile object from a JDFTx out file.

        Args:
            file_path (str | Path): The path to the JDFTx out file.
            is_bgw (bool): Mark True if data must be usable for BGW calculations. This will change the behavior of the
                parser to be stricter with certain criteria.
            none_slice_on_error (bool | None): If True, will return None if an error occurs while parsing a slice
                instead of halting the parsing process. This can be useful for parsing files with multiple slices where
                some slices may be incomplete or corrupted. If False, all slices may raise errors. If None, only the
                final slice can raise an error upon parsing (default behavior)

        Returns:
            JDFTXOutfile: The JDFTXOutfile object.
        """
        texts = read_outfile_slices(str(file_path))
        if none_slice_on_error is None:
            none_slice_bools = [i != len(texts) - 1 for i in range(len(texts))]
        else:
            none_slice_bools = [none_slice_on_error for i in range(len(texts))]
        slices = [
            JDFTXOutfileSlice._from_out_slice(text, is_bgw=is_bgw, none_on_error=none_slice_bools[i])
            for i, text in enumerate(texts)
        ]
        return cls(slices=slices)

    # TODO: Write testing for this function
    def to_jdftxinfile(self) -> JDFTXInfile:
        """
        Convert the JDFTXOutfile object to a JDFTXInfile object with the most recent structure.
        If the input structure is desired, simply fetch JDFTXOutfile.infile

        Returns:
            JDFTXInfile: A JDFTXInfile object representing the input parameters of the JDFTXOutfile.
        """
        return self.slices[-1].to_jdftxinfile()

    def __post_init__(self):
        last_slice = None
        for slc in self.slices[::-1]:
            if slc is not None:
                last_slice = slc
                break
        if last_slice is not None:
            for var in _jof_atr_from_last_slice:
                setattr(self, var, getattr(last_slice, var))
            self.trajectory = self._get_trajectory()

    def _get_trajectory(self) -> Trajectory | None:
        """Set the trajectory attribute of the JDFTXOutfile object."""
        traj = None
        for outfile_slice in self.slices:
            if outfile_slice is not None:
                if traj is None:
                    traj = outfile_slice.trajectory
                elif outfile_slice.trajectory is not None:
                    traj.extend(outfile_slice.trajectory)
        return traj

    def as_dict(self) -> dict:
        """
        Convert the JDFTXOutfile object to a dictionary.

        Returns:
            dict: A dictionary representation of the JDFTXOutfile object.
        """
        dct = {}
        for fld in self.__dataclass_fields__:
            if fld == "slices":
                dct[fld] = [slc.as_dict() for slc in self.slices]
                continue
            value = getattr(self, fld)
            dct[fld] = value
        return dct

    @deprecated(as_dict, deadline=(2025, 10, 4))
    def to_dict(self):
        return self.as_dict()

    ###########################################################################
    # Magic methods
    ###########################################################################

    def __getitem__(self, key: int | str) -> JDFTXOutfileSlice | Any:
        """Return item.

        Args:
            key (int | str): The key of the item.

        Returns:
            JDFTXOutfileSlice | Any: The value of the item.

        Raises:
            TypeError: If the key type is invalid.
        """
        val = None
        if type(key) is int:
            val = self.slices[key]
        elif type(key) is str:
            val = getattr(self, key)
        else:
            raise TypeError(f"Invalid key type: {type(key)}")
        return val

    def __len__(self) -> int:
        """Return length of JDFTXOutfile object.

        Returns:
            int: The number of geometric optimization steps in the JDFTXOutfile object.
        """
        return len(self.slices)

    def __str__(self) -> str:
        """Return string representation of JDFTXOutfile object.

        Returns:
            str: The string representation of the JDFTXOutfile object.
        """
        return pprint.pformat(self)


def _find_jdftx_out_file(calc_dir: Path) -> Path:
    """
    Find the JDFTx out file in a directory.

    Args:
        calc_dir (Path): The directory containing the JDFTx out file.

    Returns:
        Path: The path to the JDFTx out file.
    """
    out_files = _find_jdftx_dump_file(calc_dir, "out")
    if len(out_files) > 1:
        raise FileNotFoundError("Multiple JDFTx out files found in directory.")
    return out_files[0]


def _find_jdftx_dump_file(calc_dir: Path, dump_fname: str) -> list[Path]:
    """
    Find the JDFTx out file in a directory.

    Args:
        calc_dir (Path): The directory containing the JDFTx out file.

    Returns:
        Path: The path to the JDFTx out file.
    """
    dump_files = list(calc_dir.glob(f"*.{dump_fname}")) + list(calc_dir.glob(f"{dump_fname}"))
    dump_files = [f for f in dump_files if f.is_file()]
    if len(dump_files) == 0:
        raise FileNotFoundError(f"No JDFTx {dump_fname} file found in directory.")
    return dump_files


def _disambiguate_paths(paths: list[Path], dump_fname: str, prefix: str | None) -> Path:
    """
    Determine the correct dump file path from a list of paths.

    Args:
        paths (list[Path]): The paths to disambiguate.

    Returns:
        Path: The most relevant path.
    """
    if len(paths) == 1:
        return paths[0]
    _fprefix = ""
    if prefix is not None:
        _fprefix = f".{prefix}"
    for path in paths:
        if path.name == f"{dump_fname}{_fprefix}":
            return path
    raise FileNotFoundError(
        f"Multiple JDFTx {dump_fname} files found in directory, but none match the prefix: {prefix}."
    )
