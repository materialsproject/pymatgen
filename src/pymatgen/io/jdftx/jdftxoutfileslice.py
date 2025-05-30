"""JDFTx Outfile Slice Class.

This module defines the JDFTxOutfileSlice class, which is used to read and
process a JDFTx out file.
"""

from __future__ import annotations

import inspect
import math
import pprint
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from monty.dev import deprecated

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.io.jdftx.jelstep import JElSteps
    from pymatgen.io.jdftx.jminsettings import JMinSettings
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.core.trajectory import Trajectory
from pymatgen.core.units import Ha_to_eV, ang_to_bohr, bohr_to_ang
from pymatgen.io.jdftx._output_utils import (
    _init_dict_from_colon_dump_lines,
    find_all_key,
    find_first_range_key,
    find_key,
    find_key_first,
    get_colon_val,
    key_exists,
)
from pymatgen.io.jdftx.inputs import JDFTXInfile
from pymatgen.io.jdftx.joutstructures import JOutStructures, _get_joutstructures_start_idx

__author__ = "Ben Rich"

_jofs_atr_from_jstrucs = (
    "structure",
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
    "elec_nstep",
    "elec_e",
    "elec_grad_k",
    "elec_alpha",
    "elec_linmin",
    "t_s",
)


@dataclass
class JDFTXOutfileSlice:
    """A class to read and process a slice of a JDFTx out file.

    A class to read and process a slice of a JDFTx out file, where a "slice" is a segment of an out file corresponding
    to a single call of JDFTx.

    Methods:
        from_out_slice(text: list[str]): Read slice of out file into a JDFTXOutfileSlice instance.

    Attributes:
        prefix (str | None): Prefix of dump files for JDFTx calculation.
        jstrucs (JOutStructures | None): JOutStructures instance containing intermediate structures. Holds a "slices"
            attribute, which is a list of JOutStructure instances. (A JOutStructure instance functions as a Structure
            object, along with a JElSteps instance (stored as elecmindata) and other JDFTx-calculation-specific data.)
        jsettings_fluid (JMinSettings | None): JMinSettings instance containing fluid optimization settings.
        jsettings_electronic (JMinSettings | None): JMinSettings instance containing electronic optimization settings.
        jsettings_lattice (JMinSettings | None): JMinSettings instance containing lattice optimization settings.
        jsettings_ionic (JMinSettings | None): JMinSettings instance containing ionic optimization settings.
        xc_func (str | None): Exchange-correlation functional used in the calculation.
        lattice_initial (np.ndarray | None): Initial lattice matrix in Angstroms.
        lattice_final (np.ndarray | None): Final lattice matrix in Angstroms.
        lattice (np.ndarray | None): Current lattice matrix in Angstroms.
        a (float | None): Lattice parameter a in Angstroms.
        b (float | None): Lattice parameter b in Angstroms.
        c (float | None): Lattice parameter c in Angstroms.
        fftgrid (list[int] | None): Shape of FFT grid used in calculation (3 integers).
        geom_opt (bool | None): True if geometric (lattice or ionic) optimization was performed.
        geom_opt_type (str | None): Type of geometric optimization performed (lattice or ionic, where lattice implies
            ionic as well unless geometry was given in direct coordinates).
        efermi (float | None): Fermi energy in eV (may be None if eigstats are not dumped).
        egap (float | None): Band gap in eV (None if eigstats are not dumped).
        emin (float | None): Minimum energy in eV (None if eigstats are not dumped).
        emax (float | None): Maximum energy in eV (None if eigstats are not dumped).
        homo (float | None): Energy of last band-state before Fermi level (acronym for Highest Occupied Molecular
            Orbital, even though these are not molecular orbitals and this state may not be entirely occupied)
            (None if eigstats are not dumped).
        lumo (float | None): Energy of first band-state after Fermi level (acronym for Lowest Unoccupied Molecular
            Orbital, even though these are not molecular orbitals, and this state may not be entirely unoccupied)
            (None if eigstats are not dumped).
        homo_filling (float | None): Filling of "homo" band-state as calculated within this class object from the homo
            energy, Fermi level, electronic broadening type and electronic broadening parameter. (None if eigstats are
            not dumped).
        lumo_filling (float | None): Filling of "lumo" band-state as calculated within this class object from the homo
            energy, Fermi level, electronic broadening type and electronic broadening parameter. (None if eigstats are
            not dumped).
        is_metal (bool | None): True if fillings of homo and lumo band-states are off-set by 1 and 0 by at least an
            arbitrary tolerance of 0.01 (ie 1 - 0.015 and 0.012 for homo/lumo fillings would be metallic, while 1-0.001
            and 0 would not be). (Only available if eigstats was dumped).
        etype (str | None): String representation of total energy-type of system. Commonly "G" (grand-canonical
            potential) for GC calculations, and "F" for canonical (fixed electron count) calculations.
        broadening_type (str): Type of broadening for electronic filling about Fermi-level requested. Either "Fermi",
            "Cold", "MP1", or "Gauss".
        broadening (float): Magnitude of broadening for electronic filling.
        kgrid (list[int]): Shape of k-point grid used in calculation. (equivalent to k-point folding).
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
            charge. (ie total_electrons + charge).
        semicore_electrons_uncharged (int): The number of semicore electrons in the most recent JDFTx call, uncorrected
            for charge. (ie semicore_electrons + charge).
        valence_electrons_uncharged (int): The number of valence electrons in the most recent JDFTx call, uncorrected
            for charge. (ie valence_electrons + charge).
        nbands (int): The number of bands used in the most recent JDFTx call.
        atom_elements (list[str]): The list of each ion's element symbol in the most recent JDFTx call.
        atom_elements_int (list[int]): The list of ion's atomic numbers in the most recent JDFTx call.
        atom_types (list[str]): Non-repeating list of each ion's element symbol in the most recent JDFTx call.
        spintype (str): The spin type used in the most recent JDFTx call. Options are "none", "collinear".
        nspin (int): The number of spins used in the most recent JDFTx call.
        nat (int): The number of atoms in the most recent JDFTx call.
        atom_coords_initial (list[list[float]]): The initial atomic coordinates of the most recent JDFTx call.
        atom_coords_final (list[list[float]]): The final atomic coordinates of the most recent JDFTx call.
        atom_coords (list[list[float]]): The atomic coordinates of the most recent JDFTx call.
        has_solvation (bool): True if the most recent JDFTx call included a solvation calculation.
        fluid (str): The fluid used in the most recent JDFTx call.
        is_gc (bool): True if the most recent slice is a grand canonical calculation.
        is_bgw (bool): True if data must be usable for a BerkeleyGW calculation (user-set).
        has_eigstats (bool): True if eigstats were dumped in the most recent JDFTx call.
        has_parsable_pseudo (bool): True if the most recent JDFTx call used a pseudopotential that can be parsed by this
            output parser. Options are currently "GBRV" and "SG15".

    Properties:
        t_s (float | None): The total time in seconds for the calculation.
        converged (bool | None): True if calculation converged.
        trajectory (Trajectory): pymatgen Trajectory object containing intermediate Structure's of outfile slice
            calculation.
        electronic_output (dict): Dictionary with all relevant electronic information dumped from an eigstats log.
        structure (Structure): Calculation result as pymatgen Structure.
        initial_structure (Structure): Initial structure of the calculation, as read from inputs portion of out file.
        eopt_type (str | None): eopt_type from most recent JOutStructure.
        elecmindata (JElSteps): elecmindata from most recent JOutStructure.
        stress (np.ndarray | None): Stress tensor from most recent JOutStructure in units eV/Ang^3.
        strain (np.ndarray | None): Strain tensor from most recent JOutStructure (unitless).
        nstep (int | None): (geometric) nstep from most recent JOutStructure.
        e (float | None): Energy of system "etype" from most recent JOutStructure.
        grad_k (float): The final norm of the preconditioned gradient for geometric optimization of the most recent
            JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
            (written as "|grad|_K" in JDFTx output).
        alpha (float): The step size of the final geometric step in the most recent JDFTx call.
        linmin (float): The final normalized projection of the geometric step direction onto the gradient for the most
            recent JDFTx call.
        abs_magneticmoment (float | None): The absolute magnetic moment of the most recent JDFTx call.
        tot_magneticmoment (float | None): The total magnetic moment of the most recent JDFTx call.
        mu (float): The Fermi energy of the most recent JDFTx call.
        elec_e (float): The final energy of the most recent electronic optimization step.
        elec_nstep (int): The number of electronic optimization steps in the most recent JDFTx call.
        elec_grad_k (float): The final norm of the preconditioned gradient for electronic optimization of the most
            recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
            (written as "|grad|_K" in JDFTx output).
        elec_alpha (float): The step size of the final electronic step in the most recent JDFTx call.
        elec_linmin (float): The final normalized projection of the electronic step direction onto the gradient for the
            most recent JDFTx call.

    Magic Methods:
        __getattr__(name: str) -> Any: Overwrite of default __getattr__ method to allow for reference of un-defined
            attributes to the "jstrucs" class field. This referring behavior is ideally never used as (currently) all
            referrable attributes are defined properties, but is included to prevent errors in the case of future
            changes.
        __str__() -> str: Return a string representation of the class instance using pprint module.
        __repr__() -> str: Create string representation of the class instance. Overwritten from default behavior for
            dataclass so that properties are included in the string, and verbose attributes with redundant information
            are trimmed.
    """

    prefix: str | None = None

    jstrucs: JOutStructures | None = None
    jsettings_fluid: JMinSettings | None = None
    jsettings_electronic: JMinSettings | None = None
    jsettings_lattice: JMinSettings | None = None
    jsettings_ionic: JMinSettings | None = None
    constant_lattice: bool | None = None

    xc_func: str | None = None

    lattice_initial: np.ndarray | None = None
    lattice_final: np.ndarray | None = None
    lattice: np.ndarray | None = None
    a: float | None = None
    b: float | None = None
    c: float | None = None

    fftgrid: list[int] | None = None
    geom_opt: bool | None = None
    geom_opt_type: str | None = None

    # grouping fields related to electronic parameters.
    # Used by the get_electronic_output() method
    _electronic_output: ClassVar[list[str]] = [
        "efermi",
        "egap",
        "optical_egap",
        "emin",
        "emax",
        "homo",
        "lumo",
        "homo_filling",
        "lumo_filling",
        "is_metal",
    ]
    efermi: float | None = None
    egap: float | None = None
    optical_egap: float | None = None
    emin: float | None = None
    emax: float | None = None
    homo: float | None = None
    lumo: float | None = None
    homo_filling: float | None = None
    lumo_filling: float | None = None
    is_metal: bool | None = None
    etype: str | None = None

    broadening_type: str | None = None
    broadening: float | None = None
    kgrid: list | None = None
    truncation_type: str | None = None
    truncation_radius: float | None = None
    pwcut: float | None = None
    rhocut: float | None = None

    pp_type: str | None = None
    semicore_electrons: int | None = None
    valence_electrons: float | None = None
    total_electrons_uncharged: int | None = None
    semicore_electrons_uncharged: int | None = None
    valence_electrons_uncharged: int | None = None
    nbands: int | None = None

    atom_elements: list | None = None
    atom_elements_int: list | None = None
    atom_types: list | None = None
    spintype: str | None = None
    nspin: int | None = None
    nat: int | None = None
    atom_coords_initial: np.ndarray | None = None
    atom_coords_final: np.ndarray | None = None
    atom_coords: np.ndarray | None = None

    has_solvation: bool = False
    fluid: str | None = None
    is_gc: bool | None = None
    is_bgw: bool = False
    has_eigstats: bool = False
    parsable_pseudos: ClassVar[list[str]] = ["GBRV", "SG15"]
    has_parsable_pseudo: bool = False

    _total_electrons_backup: int | None = None
    total_electrons: float | None = None
    _mu_backup: float | None = None

    t_s: float | None = None
    converged: bool | None = None
    structure: Structure | None = None
    initial_structure: Structure | None = None
    trajectory: Trajectory | None = None
    electronic_output: dict | None = None
    eopt_type: str | None = None
    elecmindata: JElSteps | None = None
    stress: np.ndarray | None = None
    strain: np.ndarray | None = None
    forces: np.ndarray | None = None
    nstep: int | None = None
    e: float | None = None
    grad_k: float | None = None
    alpha: float | None = None
    linmin: float | None = None
    abs_magneticmoment: float | None = None
    tot_magneticmoment: float | None = None
    mu: float | None = None
    elec_nstep: int | None = None
    elec_e: float | None = None
    elec_grad_k: float | None = None
    elec_alpha: float | None = None
    elec_linmin: float | None = None

    def _get_mu(self) -> None | float:
        """Sets mu from most recent JOutStructure. (Equivalent to efermi)"""
        _mu = None
        if self.jstrucs is not None:
            _mu = self.jstrucs.mu
        if _mu is None:
            _mu = self._mu_backup
        return _mu

    ###########################################################################
    # Creation methods
    ###########################################################################

    # TODO: There are a littany of points in a JDFTx out file slice where an unexpected termination
    # (due to a partial calculation) could lead to a fatal error in the parser. As a fix for this, this
    # method contains a try-except block enabled by `none_on_error`. In the long term though
    # all the subobjects should be able to handle partial initialization while returning as much
    # information as possible.
    @classmethod
    def _from_out_slice(
        cls, text: list[str], is_bgw: bool = False, none_on_error: bool = True
    ) -> JDFTXOutfileSlice | None:
        """
        Read slice of out file into a JDFTXOutfileSlice instance.

        Args:
            text (list[str]): File to read.
            is_bgw (bool): True if data must be usable for a BerkeleyGW calculation.
            none_on_error (bool): If True, return None if an error occurs. If False, raise the error.

        Returns:
            JDFTXOutfileSlice | None: An instance of JDFTXOutfileSlice or None if an error occurs and
            none_on_error is True.
        """
        instance = cls()
        instance.is_bgw = is_bgw
        try:
            instance._from_out_slice_init_all(text)
        except (ValueError, IndexError, TypeError, KeyError, AttributeError):
            if none_on_error:
                return None
            raise
        return instance

    def _from_out_slice_init_all(self, text: list[str]) -> None:
        self._set_internal_infile(text)
        # self._set_min_settings(text)
        self._set_geomopt_vars(text)
        self._set_jstrucs(text)
        self._set_backup_vars(text)
        self.prefix = self._get_prefix(text)
        spintype, nspin = self._get_spinvars(text)
        self.xc_func = self._get_xc_func(text)
        self.spintype = spintype
        self.nspin = nspin
        broadening_type, broadening = self._get_broadeningvars(text)
        self.broadening_type = broadening_type
        self.broadening = broadening
        self.kgrid = self._get_kgrid(text)
        truncation_type, truncation_radius = self._get_truncationvars(text)
        self.truncation_type = truncation_type
        self.truncation_radius = truncation_radius
        self.pwcut = self._get_pw_cutoff(text)
        self.rhocut = self._get_rho_cutoff(text)
        self.fftgrid = self._get_fftgrid(text)
        self._set_eigvars(text)
        self._set_orb_fillings()
        self.is_metal = self.determine_is_metal()
        self._set_fluid(text)
        self._set_nbands(text)
        self._set_atom_vars(text)
        self._set_total_electrons()
        self._set_pseudo_vars(text)
        self._set_lattice_vars(text)
        self.has_solvation = self._check_solvation()

        # @ Cooper added @#
        self.is_gc = key_exists("target-mu", text)
        self._set_ecomponents(text)

        # Previously were properties, but are now set as attributes
        self._from_out_slice_init_all_post_init()

    def _set_internal_infile(self, text: list[str]) -> None:
        """Set the internal infile for the JDFTXOutfileSlice.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        start_line_idx = find_key("Input parsed successfully", text)
        if start_line_idx is None:
            raise ValueError("JDFTx input parsing failed on most recent call.")
        start_line_idx += 2
        end_line_idx = None
        for i in range(start_line_idx, len(text)):
            if not len(text[i].strip()):
                end_line_idx = i
                break
        if end_line_idx is None:
            raise ValueError("Calculation did not begin for this out file slice.")
        self.infile = JDFTXInfile.from_str(
            "\n".join(text[start_line_idx:end_line_idx]), validate_value_boundaries=False
        )
        self.constant_lattice = True
        if "lattice-minimize" in self.infile:
            latsteps = self.infile["lattice-minimize"]["nIterations"]
            self.constant_lattice = not (int(latsteps) > 0)

    def _set_t_s(self) -> None:
        """Return the total time in seconds for the calculation.

        Returns:
            float: The total time in seconds for the calculation.
        """
        _t_s = None
        if self.jstrucs:
            _t_s = self.jstrucs.t_s
        self.t_s = _t_s

    def _set_converged(self) -> None:
        """Return True if calculation converged.

        Returns:
            bool: True if the electronic and geometric optimization have converged (or only the former if a single-point
            calculation).
        """
        if self.jstrucs is None:
            return
        converged = self.jstrucs.elec_converged
        if self.geom_opt:
            converged = converged and self.jstrucs.geom_converged
        self.converged = converged

    def _set_trajectory(self) -> None:
        """Return pymatgen trajectory object.

        Returns:
            Trajectory: pymatgen Trajectory object containing intermediate Structure's of outfile slice calculation.
        """
        if self.jstrucs is not None:
            structures = [slc.structure for slc in self.jstrucs]
            constant_lattice = self.constant_lattice if self.constant_lattice is not None else False
            frame_properties = [slc.properties for slc in self.jstrucs]
            self.trajectory = Trajectory.from_structures(
                structures=structures,
                constant_lattice=constant_lattice,
                frame_properties=frame_properties,
            )

    def _set_electronic_output(self) -> None:
        """Return a dictionary with all relevant electronic information.

        Returns:
            dict: Dictionary with values corresponding to these keys in _electronic_output field.
        """
        dct = {}
        for field in self._electronic_output:
            if field in self.__dataclass_fields__:
                value = getattr(self, field)
                dct[field] = value
        self.electronic_output = dct

    def _from_out_slice_init_all_post_init(self) -> None:
        """Post init for running at end of "_from_out_slice_init_all" method.

        Sets class variables previously defined as properties.
        """
        self._set_t_s()
        self._set_converged()
        self._set_electronic_output()

    def _get_xc_func(self, text: list[str]) -> str | None:
        """Get the exchange-correlation functional used in the calculation.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            str: Exchange-correlation functional used.
        """
        line = find_key("elec-ex-corr", text)
        if line is None:
            return None
        return text[line].strip().split()[-1].strip()

    def _get_prefix(self, text: list[str]) -> str | None:
        """
        Get output prefix from the out file.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            str: Prefix of dump files for JDFTx calculation.
        """
        line = find_key("dump-name", text)
        if line is None:
            return None
        dumpname = text[line].split()[1]
        return dumpname.split(".")[0] if "." in dumpname else dumpname

    def _get_spinvars(self, text: list[str]) -> tuple[str, int]:
        """
        Set spintype and nspin from out file text for instance.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            tuple:
            spintype (str): Type of spin in calculation.
            nspin (int): Number of spin types in calculation.
        """
        line = find_key("spintype ", text)
        if line is None:
            return "no-spin", 1
        spintype = text[line].split()[1]
        if spintype == "no-spin":
            nspin = 1
        elif spintype == "z-spin":
            nspin = 2
        else:
            raise NotImplementedError("have not considered this spin yet")
        return spintype, nspin

    def _get_broadeningvars(self, text: list[str]) -> tuple[str | None, float]:
        """Get broadening type and value from out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            tuple[str, float]: Broadening type and parameter for electronic smearing.
        """
        line = find_key("elec-smearing ", text)
        if line is not None:
            broadening_type = text[line].split()[1]
            broadening = float(text[line].split()[2])
        else:
            broadening_type = None
            broadening = 0
        return broadening_type, broadening

    def _get_truncationvars(self, text: list[str]) -> tuple[str | None, Any | None]:
        """Get truncation type and value from out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            tuple[str, float] | tuple[None, None]: Truncation type and radius of truncation
            (if truncation_type is spherical).
        """
        maptypes = {
            "Periodic": "periodic",
            "Slab": "slab",
            "Cylindrical": "wire",
            "Wire": "wire",
            "Spherical": "spherical",
            "Isolated": "box",
        }
        line = find_key("coulomb-interaction", text)
        if line is None:
            return None, None
        truncation_type = None
        truncation_radius = None
        if line is not None:
            truncation_type = text[line].split()[1]
            if truncation_type not in maptypes:
                raise ValueError("Problem with this truncation!")
            truncation_type = maptypes[truncation_type]
            direc = None
            if len(text[line].split()) == 3:
                direc = text[line].split()[2]
            if self.is_bgw:
                if truncation_type == "slab" and direc != "001":
                    raise ValueError("BGW slab Coulomb truncation must be along z!")
                if truncation_type == "wire" and direc != "001":
                    raise ValueError("BGW wire Coulomb truncation must be periodic in z!")
            if truncation_type == "spherical":
                line = find_key("Initialized spherical truncation of radius", text)
                if line is None:
                    raise ValueError("No spherical truncation found in out file.")
                truncation_radius = float(text[line].split()[5]) / ang_to_bohr
        else:
            raise ValueError("No truncation type found in out file.")
        return truncation_type, truncation_radius

    def _get_pw_cutoff(self, text: list[str]) -> float | None:
        """Get the electron cutoff from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            float | None: Plane wave cutoff used in calculation.
        """
        line = find_key("elec-cutoff ", text)
        if line is None:
            return None
        return float(text[line].split()[1]) * Ha_to_eV

    def _get_rho_cutoff(self, text: list[str]) -> float | None:
        """Get the electron cutoff from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            float: Electron density cutoff used in calculation.
        """
        line = find_key("elec-cutoff ", text)
        if line is None:
            return None
        lsplit = text[line].split()
        if len(lsplit) == 3:
            rhocut = float(lsplit[2]) * Ha_to_eV
        else:
            if self.pwcut is None:
                self.pwcut = self._get_pw_cutoff(text)
            rhocut = None if self.pwcut is None else float(self.pwcut * 4)
        return rhocut

    def _get_fftgrid(self, text: list[str]) -> list[int] | None:
        """Get the FFT grid from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            list[int]: FFT grid used in calculation.
        """
        line = find_key_first("Chosen fftbox size", text)
        if line is None:
            return None
        return [int(x) for x in text[line].split()[6:9]]

    def _get_kgrid(self, text: list[str]) -> list[int] | None:
        """Get the kpoint grid from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            list[int]: Kpoint grid used in calculation.
        """
        line = find_key("kpoint-folding ", text)
        if line is None:
            return None
        return [int(x) for x in text[line].split()[1:4]]

    def _get_eigstats_varsdict(self, text: list[str], prefix: str | None) -> dict[str, float | None]:
        """Get the eigenvalue statistics from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.
            prefix (str): Prefix for the eigStats section in the out file.

        Returns:
            dict[str, float | None]: Dictionary of eigenvalue statistics.
        """
        varsdict: dict[str, float | None] = {}
        lines1 = find_all_key("Dumping ", text)
        lines2 = find_all_key("eigStats' ...", text)
        lines3 = [lines1[i] for i in range(len(lines1)) if lines1[i] in lines2]
        if not lines3:
            for key in list(eigstats_keymap.keys()):
                varsdict[eigstats_keymap[key]] = None
            self.has_eigstats = False
        else:
            line_start = lines3[-1]
            line_start_rel_idx = lines1.index(line_start)
            line_end = lines1[line_start_rel_idx + 1] if len(lines1) >= line_start_rel_idx + 2 else len(lines1) - 1
            _varsdict = _init_dict_from_colon_dump_lines([text[idx] for idx in range(line_start, line_end)])
            for key in _varsdict:
                varsdict[eigstats_keymap[key]] = float(_varsdict[key]) * Ha_to_eV
            self.has_eigstats = all(eigstats_keymap[key] in varsdict for key in eigstats_keymap) and all(
                eigstats_keymap[key] is not None for key in eigstats_keymap
            )
        return varsdict

    def _set_eigvars(self, text: list[str]) -> None:
        """Set the eigenvalue statistics variables.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        eigstats = self._get_eigstats_varsdict(text, self.prefix)
        for key, val in eigstats.items():
            setattr(self, key, val)
        if self.efermi is None:
            if self.mu is None:
                self.mu = self._get_mu()
            self.efermi = self.mu

    def _get_pp_type(self, text: list[str]) -> str | None:
        """Get the pseudopotential type used in calculation.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            str: Pseudopotential library used. Returns None if not GBRV or SG15
            (pseudopotentials parsable by this parser).
        """
        skey = "Reading pseudopotential file"
        line = find_key(skey, text)
        if line is None:
            raise ValueError("Unable to find pseudopotential dump lines")
        ppfile_example = text[line].split(skey)[1].split(":")[0].strip("'").strip()
        pptype = None
        readable = self.parsable_pseudos
        for _pptype in readable:
            if _pptype in ppfile_example:
                if pptype is not None:
                    if ppfile_example.index(pptype) < ppfile_example.index(_pptype):
                        pptype = _pptype
                    else:
                        pass
                else:
                    pptype = _pptype
        if pptype is not None:
            self.has_parsable_pseudo = True
        return pptype

    def _set_pseudo_vars(self, text: list[str]) -> None:
        """Set the pseudopotential variables.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        self.pp_type = self._get_pp_type(text)
        if self.has_parsable_pseudo and self.pp_type in ["GBRV", "SG15"]:
            self._set_pseudo_vars_t1(text)
        # Otherwise variables requiring parsing pseudopotential output will be kept as None

    def _set_pseudo_vars_t1(self, text: list[str]) -> None:
        """Set the pseudopotential variables for SG15 and GBRV pseudopotentials.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        all_val_lines = find_all_key("valence electrons", text)
        atom_total_elec = []
        bounds_list = get_pseudo_read_section_bounds(text)
        for bounds in bounds_list:
            startline = bounds[0]
            endline = bounds[1]
            val_lines = [x for x in all_val_lines if x < endline and x > startline]
            val_line = val_lines[0]
            val_elec = int(text[val_line].split("valence electrons")[0].strip().split()[-1])
            atom_total_elec.append(val_elec)
        total_elec_dict = {}
        if self.atom_types is not None:
            for i, atom in enumerate(self.atom_types):
                total_elec_dict[atom] = atom_total_elec[i]
        else:
            raise ValueError("Pseudopotential data cannot be allocated without atom types.")
        if self.atom_elements is None:
            raise ValueError("Atom elements not set yet.")
        # Explicit zipping due to pre-commit in three lines below
        element_total_electrons = np.array([total_elec_dict[x] for x in self.atom_elements])
        pmg_elements = [Element(x) for x in self.atom_elements]
        element_valence_electrons = np.array(
            [
                np.sum(np.array([int(v[2:]) for v in el.electronic_structure.split(".") if "]" not in v]))
                for el in pmg_elements
            ]
        )
        element_semicore_electrons = element_total_electrons - element_valence_electrons
        self.total_electrons_uncharged = np.sum(element_total_electrons)
        self.valence_electrons_uncharged = np.sum(element_valence_electrons)
        self.semicore_electrons_uncharged = np.sum(element_semicore_electrons)
        self.semicore_electrons = self.semicore_electrons_uncharged
        if (self.total_electrons is not None) and (self.semicore_electrons is not None):
            self.valence_electrons = self.total_electrons - self.semicore_electrons  # accounts for if system is charged
        else:
            raise ValueError("Total electrons and semicore electrons must be set.")

    def _collect_settings_lines(self, text: list[str], start_flag: str) -> list[int]:
        """Collect the lines of settings from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.
            start_flag (str): Key to start collecting settings lines.

        Returns:
            list[int]: List of line numbers where settings occur.
        """
        started = False
        line_texts = []
        for i, line_text in enumerate(text):
            if started:
                if line_text.strip().split()[-1].strip() == "\\":
                    line_texts.append(i)
                else:
                    started = False
            elif start_flag in line_text:
                started = True
            elif line_texts:
                break
        return line_texts

    def _create_settings_dict(self, text: list[str], start_flag: str) -> dict:
        """Get a dictionary of settings from the out file text.

        Create a dictionary of settings from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.
            start_flag (str): Key to start collecting settings lines.

        Returns:
            dict: Dictionary of settings.
        """
        line_texts = self._collect_settings_lines(text, start_flag)
        settings_dict = {}
        for line_text in line_texts:
            line_text_list = text[line_text].strip().split()
            key = line_text_list[0].lower()
            value = line_text_list[1]
            settings_dict[key] = value
        return settings_dict

    def _set_geomopt_vars(self, text: list[str]) -> None:
        """Set the geom_opt and geom_opt_type class variables.

        Set vars geom_opt and geom_opt_type for initializing self.jstrucs.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        found = False
        if "ionic-dynamics" in self.infile and int(self.infile["ionic-dynamics"]["nSteps"]) > 0:
            self.geom_opt = True
            self.geom_opt_type = "ionic"
            self.geom_opt_label = "IonicDynamics"
            found = True
        if (not found and "lattice-minimize" in self.infile) and int(
            self.infile["lattice-minimize"]["nIterations"] > 0
        ):
            self.geom_opt = True
            self.geom_opt_label = "LatticeMinimize"
            self.geom_opt_type = "lattice"
            found = True
        if (not found and "ionic-minimize" in self.infile) and int(self.infile["ionic-minimize"]["nIterations"] > 0):
            self.geom_opt = True
            self.geom_opt_label = "IonicMinimize"
            self.geom_opt_type = "ionic"
            found = True
        if not found:
            self.geom_opt = False
            self.geom_opt_type = "single point"
            self.geom_opt_label = "IonicMinimize"

    def _get_initial_structure(self, text: list[str]) -> Structure | None:
        """Get the initial structure from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.

        Returns:
            Structure | None: Initial structure of the calculation.
        """
        output_start_idx = _get_joutstructures_start_idx(text)
        if output_start_idx is None:
            init_struc = _get_init_structure(text)
        else:
            init_struc = _get_init_structure(text[:output_start_idx])
        if init_struc is None:
            raise ValueError("Provided out file slice's inputs preamble does not contain input structure data.")
        return init_struc

    def _set_jstrucs(self, text: list[str]) -> None:
        """Set the jstrucs class variable.

        Set the JStructures object to jstrucs from the out file text and all class attributes initialized from jstrucs.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        self.initial_structure = self._get_initial_structure(text)
        # In the case where no ion optimization updates are printed, JOutStructures
        if self.geom_opt_type is None:
            raise ValueError("geom_opt_type not set yet.")
        self.jstrucs = JOutStructures._from_out_slice(
            text, opt_type=self.geom_opt_type, init_struc=self.initial_structure
        )
        if self.etype is None:
            self.etype = self.jstrucs[-1].etype
        if self.jstrucs is not None:
            self._set_trajectory()
            self.mu = self._get_mu()
            for var in _jofs_atr_from_jstrucs:
                setattr(self, var, getattr(self.jstrucs, var))

    def _set_backup_vars(self, text: list[str]) -> None:
        """Set backups for important variables.

        Set backup versions of critical variables if missing from constructed jstrucs (that can be easily fetched
        through type-casting strings).

        Args:
            text (list[str]): Output of read_file for out file.
        """
        if self.total_electrons is None:
            lines = find_all_key("nElectrons", text)
            _val = None
            for line in lines[::-1]:
                _val = get_colon_val(text[line], "nElectrons:")
                if _val is not None:
                    break
            if _val is not None and np.isnan(_val):
                nval = None
            elif _val is not None:
                nval = int(_val)
            else:
                nval = None
            self._total_electrons_backup = nval

        if self.mu is None:
            lines = find_all_key("mu", text)
            _val = None
            for line in lines[::-1]:
                _val = get_colon_val(text[line], "mu:")
                if _val is not None:
                    break
            if _val is not None and np.isnan(_val):
                mval = None
            elif _val is not None:
                mval = float(_val)
            else:
                mval = None
            self._mu_backup = mval

    def _set_orb_fillings_nobroad(self, nspin: float) -> None:
        """Set the orbital fillings without broadening.

        Args:
            nspin (float): Number of spins in calculation.
        """
        self.homo_filling = 2 / nspin
        self.lumo_filling = 0

    def _set_orb_fillings_broad(
        self, nspin: float, ehomo: float, elumo: float, efermi: float, broadening_type: str, broadening: float
    ) -> None:
        """Set the orbital fillings with broadening.

        Args:
            nspin (float): Number of spins in calculation.
            ehomo (float): Energy of the highest occupied molecular orbital.
            elumo (float): Energy of the lowest unoccupied molecular orbital.
            efermi (float): Fermi energy.
            broadening_type (str): Type of broadening.
            broadening (float): Broadening parameter.
        """
        self.homo_filling = (2 / nspin) * self._calculate_filling(broadening_type, broadening, ehomo, efermi)
        self.lumo_filling = (2 / nspin) * self._calculate_filling(broadening_type, broadening, elumo, efermi)

    def _set_orb_fillings(self) -> None:
        """Set the orbital fillings.

        Calculate and set homo and lumo fillings.
        """
        if self.has_eigstats:
            if self.nspin is not None:
                if self.broadening_type is not None:
                    if self.broadening is not None:
                        if self.efermi is not None:
                            if self.homo is not None:
                                if self.lumo is not None:
                                    self._set_orb_fillings_broad(
                                        self.nspin,
                                        self.homo,
                                        self.lumo,
                                        self.efermi,
                                        self.broadening_type,
                                        self.broadening,
                                    )
                                else:
                                    raise ValueError(
                                        "Cannot set orbital fillings with broadening with self.lumo as None"
                                    )
                            else:
                                raise ValueError("Cannot set orbital fillings with broadening with self.homo as None")
                        else:
                            raise ValueError("Cannot set orbital fillings with broadening with self.efermi as None")
                    else:
                        raise ValueError("Cannot set orbital fillings with broadening with self.broadening as None")
                else:
                    self._set_orb_fillings_nobroad(self.nspin)
            else:
                raise ValueError("Cannot set homo/lumo filling with self.nspin as None")

    def _set_fluid(self, text: list[str]) -> None:
        """Set the fluid class variable.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        line = find_first_range_key("fluid ", text)
        val = text[line[0]].split()[1]
        self.fluid = None if val.lower() == "none" else val

    def _set_total_electrons(self) -> None:
        """Sets total_electrons from most recent JOutStructure."""
        tot_elec = None
        if self.jstrucs is not None:
            _tot_elec = self.jstrucs.nelectrons
            if _tot_elec is not None:
                tot_elec = _tot_elec
        if (tot_elec is None) and (self._total_electrons_backup is not None):
            tot_elec = self._total_electrons_backup
        self.total_electrons = tot_elec

    def _set_nbands(self, text: list[str]) -> None:
        """Set the Nbands class variable.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        lines = find_all_key("elec-n-bands", text)
        if len(lines):
            line = lines[0]
            nbands = int(text[line].strip().split()[-1].strip())
        else:
            lines = find_all_key("nBands:", text)
            line = lines[0]
            nbands = int(text[line].split("nBands:")[1].strip().split()[0].strip())
        self.nbands = nbands

    def _set_atom_vars(self, text: list[str]) -> None:
        """Set the atom variables.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        startline = find_key("Input parsed successfully", text)
        if startline is None:
            raise ValueError("Unable to find start line for atom lines")
        endline = find_key("---------- Initializing the Grid ----------", text)
        if endline is None:
            raise ValueError("Unable to find end line for atom lines")
        lines = find_first_range_key("ion ", text, startline=startline, endline=endline)
        atom_elements = [text[x].split()[1] for x in lines]
        self.nat = len(atom_elements)
        atom_coords = [text[x].split()[2:5] for x in lines]
        self.atom_coords_initial = np.array(atom_coords, dtype=float)
        atom_types = []
        for x in atom_elements:
            if x not in atom_types:
                atom_types.append(x)
        self.atom_elements = atom_elements
        self.atom_elements_int = [Element(x).Z for x in self.atom_elements]
        self.atom_types = atom_types
        if isinstance(self.structure, Structure):
            self.atom_coords = self.structure.cart_coords
            self.atom_coords_final = self.structure.cart_coords

    def _set_lattice_vars(self, text: list[str]) -> None:
        """Set the lattice variables.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        if self.jstrucs is not None:
            self.lattice_initial = self.jstrucs[0].lattice.matrix
            self.lattice_final = self.jstrucs[-1].lattice.matrix
            self.lattice = self.jstrucs[-1].lattice.matrix.copy()
            self.a, self.b, self.c = np.sum(self.jstrucs[-1].lattice.matrix ** 2, axis=1) ** 0.5
        else:
            raise ValueError("No structures found in out file.")

    def _set_ecomponents(self, text: list[str]) -> None:
        """
        Set the energy components dictionary from the out file text.

        Args:
            text (list[str]): Output of read_file for out file.
        """
        if self.jstrucs is not None:
            ecomp = self.jstrucs[-1].ecomponents
            if not isinstance(ecomp, dict):
                ecomp = {}
            if self.etype not in ecomp:
                if self.etype is None:
                    raise ValueError("etype not set yet.")
                ecomp[self.etype] = self.jstrucs[-1].e
            self.ecomponents = ecomp
        else:
            raise ValueError("No structures found in out file.")

    def _calculate_filling(self, broadening_type: str, broadening: float, eig: float, efermi: float) -> float:
        """
        Calculate the filling for a given eigenvalue.

        Use the broadening type, broadening value, eigenvalue, and fermi energy to calculate the filling at the
        eigenvalue.

        Args:
            broadening_type (str): Type of broadening to use.
            broadening (float): Broadening parameter.
            eig (float): Eigenvalue.
            efermi (float): Fermi energy.

        Returns:
            float: Filling at the eigenvalue.
        """
        # most broadening implementations do not have the denominator factor
        # of 2, but JDFTx does currently.
        x = (eig - efermi) / (2.0 * broadening)
        if broadening_type == "Fermi":
            filling = 0.5 * (1 - np.tanh(x))
        elif broadening_type == "Gauss":
            filling = 0.5 * (1 - math.erf(x))
        elif broadening_type == "MP1":
            filling = 0.5 * (1 - math.erf(x)) - x * np.exp(-1 * x**2) / (2 * np.pi**0.5)
        elif broadening_type == "Cold":
            filling = 0.5 * (1 - math.erf(x + 0.5**0.5)) + np.exp(-1 * (x + 0.5**0.5) ** 2) / (2 * np.pi) ** 0.5
        else:
            raise NotImplementedError("Have not added other broadening types")
        return filling

    def _determine_is_metal(self, tol_partial: float, nspin: int, homo_filling: float, lumo_filling: float) -> bool:
        """Return boolean for whether system is metallic.

        Return boolean for whether system is metallic. True if difference in filling in homo and lumo states
        exceed 0 and 1 by a given tol_partial parameter.

        Returns:
            bool: True if system is metallic
        """
        return not (homo_filling / (2 / nspin) > (1 - tol_partial) and lumo_filling / (2 / nspin) < tol_partial)

    def determine_is_metal(self) -> bool | None:
        """Determine if the system is a metal based on the fillings of homo and lumo.

        Returns:
            bool: True if system is metallic.
        """
        tol_partial = 0.01
        if self.has_eigstats:
            if self.nspin is not None:
                if self.homo_filling is not None:
                    if self.lumo_filling is not None:
                        return self._determine_is_metal(tol_partial, self.nspin, self.homo_filling, self.lumo_filling)
                    raise ValueError("Cannot determine if system is metal - self.lumo_filling undefined")
                raise ValueError("Cannot determine if system is metal - self.homo_filling undefined")
            raise ValueError("Cannot determine if system is metal - self.nspin undefined")
        return None

    def to_jdftxinfile(self) -> JDFTXInfile:
        """
        Convert the JDFTXOutfile object to a JDFTXInfile object with the most recent structure.
        If the input structure is desired, simply fetch JDFTXOutfile.infile

        Returns:
            JDFTXInfile: A JDFTXInfile object representing the input parameters of the JDFTXOutfile.
        """
        # Use internal infile as a reference for calculation parameters
        base_infile = self.infile.copy()
        # Strip references to the input
        base_infile.strip_structure_tags()
        if self.structure is None:
            return base_infile
        infile = JDFTXInfile.from_structure(self.structure)
        infile += base_infile
        return infile

    def _check_solvation(self) -> bool:
        """Check for implicit solvation.

        Returns:
            bool: True if calculation used implicit solvation.
        """
        return self.fluid is not None

    def write(self) -> None:
        """Return an error.

        Raises:
            NotImplementedError: There is no need to write a JDFTx out file.
        """
        raise NotImplementedError("There is no need to write a JDFTx out file")

    def as_dict(self) -> dict:
        """Convert dataclass to dictionary representation.

        Returns:
            dict: JDFTXOutfileSlice in dictionary format.
        """
        dct = {}
        for fld in self.__dataclass_fields__:
            value = getattr(self, fld)
            if hasattr(value, "as_dict"):
                dct[fld] = value.as_dict()
            else:
                dct[fld] = value
        return dct

    @deprecated(as_dict, deadline=(2025, 10, 4))
    def to_dict(self):
        return self.as_dict()

    # TODO: Re-do this now that there are no properties
    def __repr__(self) -> str:
        """Return string representation.

        Returns:
            str: String representation of the JDFTXOutfileSlice.
        """
        out_str = f"{self.__class__.__name__}("
        for cls in inspect.getmro(self.__class__):
            for key, value in cls.__dict__.items():
                if not key.startswith("_") and (not callable(value) or isinstance(value, property)):
                    pref = ""
                    suff = ", \n"
                    val = repr(getattr(self, key))
                    if "jsettings" in key:
                        pref = "\n    "
                    if key == "jstrucs":
                        val = "(... JOutStructures object ...)"
                    out_str += f"{pref}{key}={val}{suff}"
        out_str += ")"
        return out_str

    def __str__(self) -> str:
        """Return string representation.

        Returns:
            str: String representation of the JDFTXOutfileSlice.
        """
        return pprint.pformat(self)


eigstats_keymap = {
    "eMin": "emin",
    "HOMO": "homo",
    "mu": "efermi",
    "LUMO": "lumo",
    "eMax": "emax",
    "HOMO-LUMO gap": "egap",
    "Optical gap": "optical_egap",
}


def get_pseudo_read_section_bounds(text: list[str]) -> list[list[int]]:
    """Get the boundary line numbers for the pseudopotential read section.

    Args:
        text (list[str]): Output of read_file for out file.

    Returns:
        list[list[int]]: List of line numbers for the pseudopotential read sections.
    """
    start_lines = find_all_key("Reading pseudopotential file", text)
    section_bounds = []
    for start_line in start_lines:
        bounds = [start_line]
        for i in range(start_line, len(text)):
            if not len(text[i].strip()):
                bounds.append(i)
                break
        section_bounds.append(bounds)
    return section_bounds


def _get_init_structure(pre_out_slice: list[str]) -> Structure | None:
    """
    Return initial structure.

    Return the initial structure from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure structural data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
        contains the initial structure information.

    Returns:
        Structure | None: The initial structure if available, otherwise None.
    """
    try:
        lat_mat = _get_initial_lattice(pre_out_slice)
        coords = _get_initial_coords(pre_out_slice)
        species = _get_initial_species(pre_out_slice)
        return Structure(lattice=lat_mat, species=species, coords=coords, coords_are_cartesian=True)
    except AttributeError:
        return None


def _get_initial_lattice(pre_out_slice: list[str]) -> np.ndarray:
    """Return initial lattice.

    Return the initial lattice from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure lattice data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial lattice information.

    Returns:
        np.ndarray: The initial lattice matrix.
    """
    lat_lines = find_first_range_key("lattice  ", pre_out_slice)
    if len(lat_lines):
        lat_line = lat_lines[0]
        lat_mat = np.zeros([3, 3])
        for i in range(3):
            line_text = pre_out_slice[lat_line + i + 1].strip().split()
            for j in range(3):
                lat_mat[i, j] = float(line_text[j])
        return lat_mat.T * bohr_to_ang
    raise AttributeError("Lattice not found in pre_out_slice")


def _get_initial_coords(pre_out_slice: list[str]) -> np.ndarray:
    """Return initial coordinates.

    Return the initial coordinates from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure coordinate data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial coordinates information.

    Returns:
        np.ndarray: The initial coordinates.
    """
    lines = _get_ion_lines(pre_out_slice)
    coords = np.zeros([len(lines), 3])
    for i, line in enumerate(lines):
        line_text = pre_out_slice[line].strip().split()[2:]
        for j in range(3):
            coords[i, j] = float(line_text[j])
    coords_type_lines = find_first_range_key("coords-type", pre_out_slice)
    if len(coords_type_lines):
        coords_type = pre_out_slice[coords_type_lines[0]].strip().split()[1]
        if coords_type.lower() != "cartesian":
            coords = np.dot(coords, _get_initial_lattice(pre_out_slice))
        else:
            coords *= bohr_to_ang
    return coords


def _get_initial_species(pre_out_slice: list[str]) -> list[str]:
    """Return initial species.

    Return the initial species from the pre_out_slice, corresponding to all data cut from JOutStructure list
    initialization. This is needed to ensure species data that is not being updated (and therefore not being
    logged in the out file) is still available.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the initial species information.

    Returns:
        list[str]: The initial species.
    """
    lines = _get_ion_lines(pre_out_slice)
    species_strs = []
    for line in lines:
        species_strs.append(pre_out_slice[line].strip().split()[1])
    return species_strs


def _get_ion_lines(pre_out_slice: list[str]) -> list[int]:
    """Return ion lines.

    Return the ion lines from the pre_out_slice, ensuring that all the ion lines are consecutive.

    Args:
        pre_out_slice (list[str]): A slice of a JDFTx out file (individual call of JDFTx) that
            contains the ion lines information.

    Returns:
        list[int]: The ion lines.
    """
    _lines = find_first_range_key("ion ", pre_out_slice)
    if not len(_lines):
        raise AttributeError("Ion lines not found in pre_out_slice")
    gaps = [_lines[i + 1] - _lines[i] for i in range(len(_lines) - 1)]
    if not all(g == 1 for g in gaps):
        # TODO: Write the fix for this case
        raise AttributeError("Ion lines not consecutive in pre_out_slice")
    return _lines
