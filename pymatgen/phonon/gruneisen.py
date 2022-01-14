# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to define a Grueneisen band structure.
"""

import numpy as np
import scipy.constants as const
from monty.dev import requires
from monty.json import MSONable

from pymatgen.core import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.units import amu_to_kg
from pymatgen.phonon.bandstructure import PhononBandStructure, PhononBandStructureSymmLine
from pymatgen.phonon.dos import PhononDos

try:
    import phonopy
    from phonopy.phonon.dos import TotalDos
except ImportError as ex:
    print(ex)
    phonopy = None

__author__ = "A. Bonkowski, J. George, G. Petretto"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__maintainer__ = "A. Bonkowski, J. George"
__email__ = "alexander.bonkowski@rwth-aachen.de, janine.george@bam.de"
__status__ = "Production"
__date__ = "Apr 11, 2021"


class GruneisenParameter(MSONable):
    """
    Class for Grueneisen parameters on a regular grid.
    """

    def __init__(
        self,
        qpoints,
        gruneisen,
        frequencies,
        multiplicities=None,
        structure=None,
        lattice=None,
    ):
        """

        Args:
            qpoints: list of qpoints as numpy arrays, in frac_coords of the given lattice by default
            gruneisen: list of gruneisen parameters as numpy arrays, shape: (3*len(structure), len(qpoints))
            frequencies: list of phonon frequencies in THz as a numpy array with shape (3*len(structure), len(qpoints))
            multiplicities: list of multiplicities
            structure: The crystal structure (as a pymatgen Structure object) associated with the gruneisen parameters.
            lattice: The reciprocal lattice as a pymatgen Lattice object. Pymatgen uses the physics convention of
                     reciprocal lattice vectors WITH a 2*pi coefficient

        """

        self.qpoints = qpoints
        self.gruneisen = gruneisen
        self.frequencies = frequencies
        self.multiplicities = multiplicities
        self.lattice = lattice
        self.structure = structure

    def average_gruneisen(self, t=None, squared=True, limit_frequencies=None):
        """
        Calculates the average of the Gruneisen based on the values on the regular grid.
        If squared is True the average will use the squared value of the Gruneisen and a squared root
        is performed on the final result.
        Values associated to negative frequencies will be ignored.
        See Scripta Materialia 129, 88 for definitions.
        Adapted from classes in abipy that have been written by Guido Petretto (UCLouvain).

        Args:
            t: the temperature at which the average Gruneisen will be evaluated. If None the acoustic Debye
                temperature is used (see acoustic_debye_temp).
            squared: if True the average is performed on the squared values of the Grueneisen.
            limit_frequencies: if None (default) no limit on the frequencies will be applied.
                Possible values are "debye" (only modes with frequencies lower than the acoustic Debye
                temperature) and "acoustic" (only the acoustic modes, i.e. the first three modes).

        Returns:
            The average Gruneisen parameter

        """
        if t is None:
            t = self.acoustic_debye_temp

        w = self.frequencies
        wdkt = w * const.tera / (const.value("Boltzmann constant in Hz/K") * t)
        exp_wdkt = np.exp(wdkt)
        cv = np.choose(
            w > 0, (0, const.value("Boltzmann constant in eV/K") * wdkt ** 2 * exp_wdkt / (exp_wdkt - 1) ** 2)
        )  # in eV

        gamma = self.gruneisen

        if squared:
            gamma = gamma ** 2

        if limit_frequencies == "debye":
            adt = self.acoustic_debye_temp
            ind = np.where((w >= 0) & (w <= adt * const.value("Boltzmann constant in Hz/K")))
        elif limit_frequencies == "acoustic":
            w_acoustic = w[:, :3]
            ind = np.where(w_acoustic >= 0)
        elif limit_frequencies is None:
            ind = np.where(w >= 0)
        else:
            raise ValueError(f"{limit_frequencies} is not an accepted value for limit_frequencies")

        weights = self.multiplicities
        g = np.dot(weights[ind[0]], np.multiply(cv, gamma)[ind]).sum() / np.dot(weights[ind[0]], cv[ind]).sum()

        if squared:
            g = np.sqrt(g)

        return g

    def thermal_conductivity_slack(self, squared=True, limit_frequencies=None, theta_d=None, t=None):
        """
        Calculates the thermal conductivity at the acoustic Debye temperature with the Slack formula,
        using the average Gruneisen.
        Adapted from abipy.

        Args:
            squared (bool): if True the average is performed on the squared values of the Gruenisen
            limit_frequencies: if None (default) no limit on the frequencies will be applied.
                Possible values are "debye" (only modes with frequencies lower than the acoustic Debye
                temperature) and "acoustic" (only the acoustic modes, i.e. the first three modes).
            theta_d: the temperature used to estimate the average of the Gruneisen used in the
                Slack formula. If None the acoustic Debye temperature is used (see
                acoustic_debye_temp). Will also be considered as the Debye temperature in the
                Slack formula.
            t: temperature at which the thermal conductivity is estimated. If None the value at
                the calculated acoustic Debye temperature is given. The value is obtained as a
                simple rescaling of the value at the Debye temperature.

        Returns:
            The value of the thermal conductivity in W/(m*K)

        """
        average_mass = np.mean([s.specie.atomic_mass for s in self.structure]) * amu_to_kg
        if theta_d is None:
            theta_d = self.acoustic_debye_temp
        mean_g = self.average_gruneisen(t=theta_d, squared=squared, limit_frequencies=limit_frequencies)

        f1 = 0.849 * 3 * (4 ** (1.0 / 3.0)) / (20 * np.pi ** 3 * (1 - 0.514 * mean_g ** -1 + 0.228 * mean_g ** -2))
        f2 = (const.k * theta_d / const.hbar) ** 2
        f3 = const.k * average_mass * self.structure.volume ** (1.0 / 3.0) * 1e-10 / (const.hbar * mean_g ** 2)
        k = f1 * f2 * f3

        if t is not None:
            k *= theta_d / t

        return k

    @property  # type: ignore
    @requires(phonopy, "This method requires phonopy to be installed")
    def tdos(self):
        """
        The total DOS (re)constructed from the gruneisen.yaml file
        """

        # Here, we will reuse phonopy classes
        class TempMesh:
            """
            Temporary Class
            """

        a = TempMesh()
        a.frequencies = np.transpose(self.frequencies)
        a.weights = self.multiplicities

        b = TotalDos(a)
        b.run()

        return b

    @property
    def phdos(self):
        """

        Returns: PhononDos object

        """
        return PhononDos(self.tdos.frequency_points, self.tdos.dos)

    @property
    def debye_temp_limit(self):
        """
        Debye temperature in K. Adapted from apipy.
        """
        from scipy.interpolate import UnivariateSpline

        f_mesh = self.tdos.frequency_points * const.tera
        dos = self.tdos.dos

        i_a = UnivariateSpline(f_mesh, dos * f_mesh ** 2, s=0).integral(f_mesh[0], f_mesh[-1])
        i_b = UnivariateSpline(f_mesh, dos, s=0).integral(f_mesh[0], f_mesh[-1])

        integrals = i_a / i_b
        t_d = np.sqrt(5 / 3 * integrals) / const.value("Boltzmann constant in Hz/K")

        return t_d

    def debye_temp_phonopy(self, freq_max_fit=None):
        """
        Get Debye temperature in K as implemented in phonopy.
        Args:
            freq_max_fit: Maximum frequency to include for fitting.
                          Defaults to include first quartile of frequencies.

        Returns:
            Debye temperature in K.

        """
        # Use of phonopy classes to compute Debye frequency
        t = self.tdos
        t.set_Debye_frequency(num_atoms=self.structure.num_sites, freq_max_fit=freq_max_fit)
        f_d = t.get_Debye_frequency()  # in THz
        # f_d in THz is converted in a temperature (K)
        t_d = const.value("Planck constant") * f_d * const.tera / const.value("Boltzmann constant")

        return t_d

    @property
    def acoustic_debye_temp(self):
        """
        Acoustic Debye temperature in K, i.e. the Debye temperature divided by nsites**(1/3).
        Adapted from abipy.
        """
        return self.debye_temp_limit / self.structure.num_sites ** (1 / 3)


class GruneisenPhononBandStructure(PhononBandStructure):
    """
    This is the most generic phonon band structure data possible
    it's defined by a list of qpoints + frequencies for each of them.
    Additional information may be given for frequencies at Gamma, where
    non-analytical contribution may be taken into account.
    """

    def __init__(
        self,
        qpoints,
        frequencies,
        gruneisenparameters,
        lattice,
        eigendisplacements=None,
        labels_dict=None,
        coords_are_cartesian=False,
        structure=None,
    ):
        """
        Args:
            qpoints: list of qpoint as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in THz as a numpy array with shape
                (3*len(structure), len(qpoints)). The First index of the array
                refers to the band and the second to the index of the qpoint.
            gruneisenparameters: list of Grueneisen parameters with the same structure
                frequencies.
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient.
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                he First index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """

        PhononBandStructure.__init__(
            self,
            qpoints,
            frequencies,
            lattice,
            nac_frequencies=None,
            eigendisplacements=eigendisplacements,
            nac_eigendisplacements=None,
            labels_dict=labels_dict,
            coords_are_cartesian=coords_are_cartesian,
            structure=structure,
        )
        self.gruneisen = gruneisenparameters

    def as_dict(self):
        """

        Returns:
            MSONable (dict)

        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "lattice_rec": self.lattice_rec.as_dict(),
            "qpoints": [],
        }
        # qpoints are not Kpoint objects dicts but are frac coords. This makes
        # the dict smaller and avoids the repetition of the lattice
        for q in self.qpoints:
            d["qpoints"].append(q.as_dict()["fcoords"])
        d["bands"] = self.bands.tolist()
        d["labels_dict"] = {}
        for kpoint_letter, kpoint_object in self.labels_dict.items():
            d["labels_dict"][kpoint_letter] = kpoint_object.as_dict()["fcoords"]
        # split the eigendisplacements to real and imaginary part for serialization
        d["eigendisplacements"] = dict(
            real=np.real(self.eigendisplacements).tolist(), imag=np.imag(self.eigendisplacements).tolist()
        )
        d["gruneisen"] = self.gruneisen.tolist()
        if self.structure:
            d["structure"] = self.structure.as_dict()

        return d

    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d (dict): Dict representation

        Returns:
            GruneisenPhononBandStructure: Phonon band structure with Grueneisen parameters.

        """

        lattice_rec = Lattice(d["lattice_rec"]["matrix"])
        eigendisplacements = np.array(d["eigendisplacements"]["real"]) + np.array(d["eigendisplacements"]["imag"]) * 1j
        structure = Structure.from_dict(d["structure"]) if "structure" in d else None
        return cls(
            qpoints=d["qpoints"],
            frequencies=np.array(d["bands"]),
            gruneisenparameters=np.array(d["gruneisen"]),
            lattice=lattice_rec,
            eigendisplacements=eigendisplacements,
            labels_dict=d["labels_dict"],
            structure=structure,
        )


class GruneisenPhononBandStructureSymmLine(GruneisenPhononBandStructure, PhononBandStructureSymmLine):
    """
    This object stores a GruneisenPhononBandStructureSymmLine together with Grueneisen parameters
    for every frequency.
    """

    def __init__(
        self,
        qpoints,
        frequencies,
        gruneisenparameters,
        lattice,
        eigendisplacements=None,
        labels_dict=None,
        coords_are_cartesian=False,
        structure=None,
    ):
        """

        Args:
            qpoints: list of qpoints as numpy arrays, in frac_coords of the
                given lattice by default
            frequencies: list of phonon frequencies in eV as a numpy array with shape
                (3*len(structure), len(qpoints))
            gruneisenparameters: list of Grueneisen parameters as a numpy array with the
                shape (3*len(structure), len(qpoints))
            lattice: The reciprocal lattice as a pymatgen Lattice object.
                Pymatgen uses the physics convention of reciprocal lattice vectors
                WITH a 2*pi coefficient
            eigendisplacements: the phonon eigendisplacements associated to the
                frequencies in cartesian coordinates. A numpy array of complex
                numbers with shape (3*len(structure), len(qpoints), len(structure), 3).
                he First index of the array refers to the band, the second to the index
                of the qpoint, the third to the atom in the structure and the fourth
                to the cartesian coordinates.
            labels_dict: (dict) of {} this links a qpoint (in frac coords or
                cartesian coordinates depending on the coords) to a label.
            coords_are_cartesian: Whether the qpoint coordinates are cartesian.
            structure: The crystal structure (as a pymatgen Structure object)
                associated with the band structure. This is needed if we
                provide projections to the band structure
        """
        GruneisenPhononBandStructure.__init__(
            self,
            qpoints=qpoints,
            frequencies=frequencies,
            gruneisenparameters=gruneisenparameters,
            lattice=lattice,
            eigendisplacements=eigendisplacements,
            labels_dict=labels_dict,
            coords_are_cartesian=coords_are_cartesian,
            structure=structure,
        )

        PhononBandStructureSymmLine._reuse_init(
            self, eigendisplacements=eigendisplacements, frequencies=frequencies, has_nac=False, qpoints=qpoints
        )

    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d: Dict representation

        Returns: GruneisenPhononBandStructureSummLine

        """
        lattice_rec = Lattice(d["lattice_rec"]["matrix"])
        eigendisplacements = np.array(d["eigendisplacements"]["real"]) + np.array(d["eigendisplacements"]["imag"]) * 1j
        structure = Structure.from_dict(d["structure"]) if "structure" in d else None
        return cls(
            qpoints=d["qpoints"],
            frequencies=np.array(d["bands"]),
            gruneisenparameters=np.array(d["gruneisen"]),
            lattice=lattice_rec,
            eigendisplacements=eigendisplacements,
            labels_dict=d["labels_dict"],
            structure=structure,
        )
