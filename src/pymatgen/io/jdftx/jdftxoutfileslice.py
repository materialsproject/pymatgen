"""JDFTx Outfile Slice Class.

This module defines the JDFTxOutfileSlice class, which is used to read and
process a JDFTx out file.

"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

if TYPE_CHECKING:
    from pymatgen.core import Structure
from pymatgen.core.trajectory import Trajectory
from pymatgen.core.units import Ha_to_eV, ang_to_bohr

from atomate2.jdftx.io.data import atom_valence_electrons
from atomate2.jdftx.io.jdftxoutfileslice_helpers import (
    find_all_key,
    find_first_range_key,
    find_key,
    find_key_first,
    get_pseudo_read_section_bounds,
    key_exists,
)
from atomate2.jdftx.io.jminsettings import (
    JMinSettingsElectronic,
    JMinSettingsFluid,
    JMinSettingsIonic,
    JMinSettingsLattice,
)
from atomate2.jdftx.io.joutstructures import JOutStructures


class ClassPrintFormatter:
    """Generic class object print formatter.

    Generic class object print formatter.
    """

    def __str__(self) -> str:
        """Return class object as str for readable format in command line."""
        return (
            str(self.__class__)
            + "\n"
            + "\n".join(
                str(item) + " = " + str(self.__dict__[item])
                for item in sorted(self.__dict__)
            )
        )


@dataclass
class JDFTXOutfileSlice(ClassPrintFormatter):
    """A class to read and process a JDFTx out file.

    A class to read and process a JDFTx out file.

    Attributes
    ----------
        see JDFTx documentation for tag info and typing
    """

    prefix: str = None

    jstrucs: JOutStructures = None
    jsettings_fluid: (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ) = None
    jsettings_electronic: (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ) = None
    jsettings_lattice: (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ) = None
    jsettings_ionic: (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ) = None

    xc_func: str = None

    # lattice_initial: list[list[float]] = None
    # lattice_final: list[list[float]] = None
    # lattice: list[list[float]] = None
    lattice_initial: np.ndarray = None
    lattice_final: np.ndarray = None
    lattice: np.ndarray = None
    a: float = None
    b: float = None
    c: float = None

    fftgrid: list[int] = None
    geom_opt: bool = None
    geom_opt_type: str = None

    # grouping fields related to electronic parameters.
    # Used by the get_electronic_output() method
    _electronic_output: ClassVar[list[str]] = [
        "efermi",
        "egap",
        "emin",
        "emax",
        "homo",
        "lumo",
        "homo_filling",
        "lumo_filling",
        "is_metal",
    ]
    efermi: float = None
    egap: float = None
    emin: float = None
    emax: float = None
    homo: float = None
    lumo: float = None
    homo_filling: float = None
    lumo_filling: float = None
    is_metal: bool = None
    etype: str = None

    broadening_type: str = None
    broadening: float = None
    kgrid: list = None
    truncation_type: str = None
    truncation_radius: float = None
    pwcut: float = None
    rhocut: float = None

    pp_type: str = None
    total_electrons: float = None
    semicore_electrons: int = None
    valence_electrons: float = None
    total_electrons_uncharged: int = None
    semicore_electrons_uncharged: int = None
    valence_electrons_uncharged: int = None
    nbands: int = None

    atom_elements: list = None
    atom_elements_int: list = None
    atom_types: list = None
    spintype: str = None
    nspin: int = None
    nat: int = None
    atom_coords_initial: list[list[float]] = None
    atom_coords_final: list[list[float]] = None
    atom_coords: list[list[float]] = None

    has_solvation: bool = False
    fluid: str = None
    is_gc: bool = None

    @property
    def t_s(self) -> float:
        """Return the total time in seconds for the calculation.

        Return the total time in seconds for the calculation.

        Returns
        -------
        t_s: float
            The total time in seconds for the calculation
        """
        t_s = None
        if self.jstrucs:
            t_s = self.jstrucs.t_s
        return t_s

    @property
    def is_converged(self) -> bool:
        """Return True if calculation converged.

        Return True if the electronic and geometric optimization have converged
        (or only the former if a single-point calculation)

        Returns
        -------
        converged: bool
            True if calculation converged
        """
        converged = self.jstrucs.elec_converged
        if self.geom_opt:
            converged = converged and self.jstrucs.geom_converged
        return converged

    @property
    def trajectory(self) -> Trajectory:
        """Return pymatgen trajectory object.

        Return pymatgen trajectory object containing intermediate Structure's
        of outfile slice calculation.

        Returns
        -------
        traj: Trajectory
            pymatgen Trajectory object
        """
        constant_lattice = self.jsettings_lattice.niterations == 0
        return Trajectory.from_structures(
            structures=self.jstrucs, constant_lattice=constant_lattice
        )

    @property
    def electronic_output(self) -> dict:
        """Return a dictionary with all relevant electronic information.

        Return dict with values corresponding to these keys in _electronic_output
        field.
        """
        dct = {}
        for field in self.__dataclass_fields__:
            if field in self._electronic_output:
                value = getattr(self, field)
                dct[field] = value
        return dct

    @property
    def structure(self) -> Structure:
        """Return calculation result as pymatgen Structure.

        Return calculation result as pymatgen Structure.

        Returns
        -------
        structure: Structure
            pymatgen Structure object
        """
        return self.jstrucs[-1]

    @classmethod
    def from_out_slice(cls, text: list[str]) -> JDFTXOutfileSlice:
        """Read slice of out file into a JDFTXOutfileSlice instance.

        Read slice of out file into a JDFTXOutfileSlice instance.

        Parameters
        ----------
        text: list[str]
            file to read
        """
        instance = cls()

        instance.set_min_settings(text)
        instance.set_geomopt_vars(text)
        instance.set_jstrucs(text)
        instance.prefix = instance.get_prefix(text)
        spintype, nspin = instance.get_spinvars(text)
        instance.xc_func = instance.get_xc_func(text)
        instance.spintype = spintype
        instance.nspin = nspin
        broadening_type, broadening = instance.get_broadeningvars(text)
        instance.broadening_type = broadening_type
        instance.broadening = broadening
        instance.kgrid = instance.get_kgrid(text)
        truncation_type, truncation_radius = instance.get_truncationvars(text)
        instance.truncation_type = truncation_type
        instance.truncation_radius = truncation_radius
        instance.pwcut = instance.get_pw_cutoff(text)
        instance.rhocut = instance.get_rho_cutoff(text)
        instance.fftgrid = instance.get_fftgrid(text)
        instance.set_eigvars(text)
        instance.set_orb_fillings()
        instance.is_metal = instance.determine_is_metal()
        instance.set_fluid(text)
        instance.set_total_electrons(text)
        instance.set_nbands(text)
        instance.set_atom_vars(text)
        instance.set_pseudo_vars(text)
        instance.set_lattice_vars(text)
        instance.has_solvation = instance.check_solvation()

        # @ Cooper added @#
        instance.is_gc = key_exists("target-mu", text)
        instance.set_ecomponents(text)
        # instance._build_trajectory(templines)

        return instance

    def get_xc_func(self, text: list[str]) -> str:
        """Get the exchange-correlation functional used in the calculation.

        Get the exchange-correlation functional used in the calculation.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        xc_func: str
            exchange-correlation functional used
        """
        line = find_key("elec-ex-corr", text)
        return text[line].strip().split()[-1].strip()

    def get_prefix(self, text: list[str]) -> str:
        """Get output prefix from the out file.

        Get output prefix from the out file.

        Parameters
        ----------
            text: list[str]
                output of read_file for out file

        Returns
        -------
            prefix: str
                prefix of dump files for JDFTx calculation
        """
        prefix = None
        line = find_key("dump-name", text)
        dumpname = text[line].split()[1]
        if "." in dumpname:
            prefix = dumpname.split(".")[0]
        return prefix

    def get_spinvars(self, text: list[str]) -> tuple[str, int]:
        """Set spintype and nspin from out file text for instance.

        Set spintype and nspin from out file text for instance.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        spintype: str
            type of spin in calculation
        nspin: int
            number of spin types in calculation
        """
        line = find_key("spintype ", text)
        spintype = text[line].split()[1]
        if spintype == "no-spin":
            spintype = None
            nspin = 1
        elif spintype == "z-spin":
            nspin = 2
        else:
            raise NotImplementedError("have not considered this spin yet")
        return spintype, nspin

    def get_broadeningvars(self, text: list[str]) -> tuple[str, float]:
        """Get broadening type and value from out file text.

        Get broadening type and value from out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        broadening_type: str
            type of electronic smearing
        broadening: float
            parameter for electronic smearing
        """
        line = find_key("elec-smearing ", text)
        if line is not None:
            broadening_type = text[line].split()[1]
            broadening = float(text[line].split()[2])
        else:
            broadening_type = None
            broadening = 0
        return broadening_type, broadening

    def get_truncationvars(self, text: list[str]) -> tuple[str, float]:
        """Get truncation type and value from out file text.

        Get truncation type and value from out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        truncation_type: str
            type of coulomb truncation
        truncation_radius: float | None
            radius of truncation (if truncation_type is spherical)
        """
        maptypes = {
            "Periodic": None,
            "Slab": "slab",
            "Cylindrical": "wire",
            "Wire": "wire",
            "Spherical": "spherical",
            "Isolated": "box",
        }
        line = find_key("coulomb-interaction", text)
        truncation_type = None
        truncation_radius = None
        if line is not None:
            truncation_type = text[line].split()[1]
            truncation_type = maptypes[truncation_type]
            direc = None
            if len(text[line].split()) == 3:
                direc = text[line].split()[2]
            if truncation_type == "slab" and direc != "001":
                raise ValueError("BGW slab Coulomb truncation must be along z!")
            if truncation_type == "wire" and direc != "001":
                raise ValueError(
                    "BGW wire Coulomb truncation must be periodic \
                                 in z!"
                )
            if truncation_type == "error":
                raise ValueError("Problem with this truncation!")
            if truncation_type == "spherical":
                line = find_key("Initialized spherical truncation of radius", text)
                truncation_radius = float(text[line].split()[5]) / ang_to_bohr
        return truncation_type, truncation_radius

    def get_pw_cutoff(self, text: list[str]) -> float:
        """Get the electron cutoff from the out file text.

        Get the electron cutoff from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        pwcut: float
            plane wave cutoff used in calculation
        """
        line = find_key("elec-cutoff ", text)
        return float(text[line].split()[1]) * Ha_to_eV

    def get_rho_cutoff(self, text: list[str]) -> float:
        """Get the electron cutoff from the out file text.

        Get the electron cutoff from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        rhocut: float
            electron density cutoff used in calculation
        """
        line = find_key("elec-cutoff ", text)
        lsplit = text[line].split()
        if len(lsplit) == 3:
            rhocut = float(lsplit[2]) * Ha_to_eV
        else:
            pwcut = self.pwcut
            if self.pwcut is None:
                pwcut = self.get_pw_cutoff(text)
            rhocut = float(pwcut * 4)
        return rhocut

    def get_fftgrid(self, text: list[str]) -> list[int]:
        """Get the FFT grid from the out file text.

        Get the FFT grid from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        fftgrid: list[int]
            FFT grid used in calculation
        """
        line = find_key_first("Chosen fftbox size", text)
        return [int(x) for x in text[line].split()[6:9]]

    def get_kgrid(self, text: list[str]) -> list[int]:
        """Get the kpoint grid from the out file text.

        Get the kpoint grid from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        kgrid: list[int]
            kpoint grid used in calculation
        """
        line = find_key("kpoint-folding ", text)
        return [int(x) for x in text[line].split()[1:4]]

    def get_eigstats_varsdict(
        self, text: list[str], prefix: str | None
    ) -> dict[str, float]:
        """Get the eigenvalue statistics from the out file text.

        Get the eigenvalue statistics from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        prefix: str
            prefix for the eigStats section in the out file

        Returns
        -------
        varsdict: dict[str, float]
            dictionary of eigenvalue statistics
        """
        varsdict = {}
        _prefix = ""
        if prefix is not None:
            _prefix = f"{prefix}."
        line = find_key(f"Dumping '{_prefix}eigStats' ...", text)
        if line is None:
            raise ValueError(
                'Must run DFT job with "dump End EigStats" to get summary gap\
                      information!'
            )
        varsdict["emin"] = float(text[line + 1].split()[1]) * Ha_to_eV
        varsdict["homo"] = float(text[line + 2].split()[1]) * Ha_to_eV
        varsdict["efermi"] = float(text[line + 3].split()[2]) * Ha_to_eV
        varsdict["lumo"] = float(text[line + 4].split()[1]) * Ha_to_eV
        varsdict["emax"] = float(text[line + 5].split()[1]) * Ha_to_eV
        varsdict["egap"] = float(text[line + 6].split()[2]) * Ha_to_eV
        return varsdict

    def set_eigvars(self, text: list[str]) -> None:
        """Set the eigenvalue statistics variables.

        Set the eigenvalue statistics variables.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        eigstats = self.get_eigstats_varsdict(text, self.prefix)
        self.emin = eigstats["emin"]
        self.homo = eigstats["homo"]
        self.efermi = eigstats["efermi"]
        self.lumo = eigstats["lumo"]
        self.emax = eigstats["emax"]
        self.egap = eigstats["egap"]

    def get_pp_type(self, text: list[str]) -> str:
        """Get the pseudopotential type used in calculation.

        Get the pseudopotential type used in calculation.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file

        Returns
        -------
        pptype: str
            Pseudopotential library used
        """
        skey = "Reading pseudopotential file"
        line = find_key(skey, text)
        ppfile_example = text[line].split(skey)[1].split(":")[0].strip("'")
        pptype = None
        readable = ["GBRV", "SG15"]
        for _pptype in readable:
            if _pptype in ppfile_example:
                if pptype is not None:
                    if ppfile_example.index(pptype) < ppfile_example.index(_pptype):
                        pptype = _pptype
                    else:
                        pass
                else:
                    pptype = _pptype
        if pptype is None:
            raise ValueError(
                f"Could not determine pseudopotential type from file name\
                      {ppfile_example}"
            )
        return pptype

    def set_pseudo_vars(self, text: list[str]) -> None:
        """Set the pseudopotential variables.

        Set the pseudopotential variables.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        self.pp_type = self.get_pp_type(text)
        if self.pp_type in ["SG15", "GBRV"]:
            self.set_pseudo_vars_t1(text)
        else:
            raise NotImplementedError(
                "Outfile parsing requires SG15 or\
                                       GBRV pseudos"
            )

    def set_pseudo_vars_t1(self, text: list[str]) -> None:
        """Set the pseudopotential variables for SG15 pseudopotentials.

        Set the pseudopotential variables for SG15 pseudopotentials.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        all_val_lines = find_all_key("valence electrons", text)
        atom_total_elec = []
        bounds_list = get_pseudo_read_section_bounds(text)
        for bounds in bounds_list:
            startline = bounds[0]
            endline = bounds[1]
            val_lines = [x for x in all_val_lines if x < endline and x > startline]
            val_line = val_lines[0]
            val_elec = int(
                text[val_line].split("valence electrons")[0].strip().split()[-1]
            )
            atom_total_elec.append(val_elec)
        total_elec_dict = dict(zip(self.atom_types, atom_total_elec))
        element_total_electrons = np.array(
            [total_elec_dict[x] for x in self.atom_elements]
        )
        element_valence_electrons = np.array(
            [atom_valence_electrons[x] for x in self.atom_elements]
        )
        element_semicore_electrons = element_total_electrons - element_valence_electrons
        self.total_electrons_uncharged = np.sum(element_total_electrons)
        self.valence_electrons_uncharged = np.sum(element_valence_electrons)
        self.semicore_electrons_uncharged = np.sum(element_semicore_electrons)
        self.semicore_electrons = self.semicore_electrons_uncharged
        self.valence_electrons = (
            self.total_electrons - self.semicore_electrons
        )  # accounts for if system is charged

    def _collect_settings_lines(self, text: list[str], start_flag: str) -> list[int]:
        """Collect the lines of settings from the out file text.

        Collect the lines of settings from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        start_flag: str
            key to start collecting settings lines

        Returns
        -------
        lines: list[int]
            list of line numbers where settings occur
        """
        started = False
        lines = []
        for i, line in enumerate(text):
            if started:
                if line.strip().split()[-1].strip() == "\\":
                    lines.append(i)
                else:
                    started = False
            elif start_flag in line:
                started = True
                # lines.append(i) # we DONT want to do this
            elif len(lines):
                break
        return lines

    def _create_settings_dict(self, text: list[str], start_flag: str) -> dict:
        """Get a dictionary of settings from the out file text.

        Create a dictionary of settings from the out file text

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        start_flag: str
            key to start collecting settings lines

        Returns
        -------
        settings_dict: dict
            dictionary of settings
        """
        lines = self._collect_settings_lines(text, start_flag)
        settings_dict = {}
        for line in lines:
            line_text_list = text[line].strip().split()
            key = line_text_list[0].lower()
            value = line_text_list[1]
            settings_dict[key] = value
        return settings_dict

    def get_settings_object(
        self,
        text: list[str],
        settings_class: type[
            JMinSettingsElectronic
            | JMinSettingsFluid
            | JMinSettingsIonic
            | JMinSettingsLattice
        ],
    ) -> (
        JMinSettingsElectronic
        | JMinSettingsFluid
        | JMinSettingsIonic
        | JMinSettingsLattice
    ):
        """Get appropriate JMinSettings mutant.

        Get the settings object from the out file text

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        settings_class: Type[JMinSettings]
            settings class to create object from

        Returns
        -------
        settings_obj: JMinSettings
            settings object
        """
        settings_dict = self._create_settings_dict(text, settings_class.start_flag)
        return settings_class(**settings_dict) if len(settings_dict) else None

    def set_min_settings(self, text: list[str]) -> None:
        """Set the settings objects from the out file text.

        Set the settings objects from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        self.jsettings_fluid = self.get_settings_object(text, JMinSettingsFluid)
        self.jsettings_electronic = self.get_settings_object(
            text, JMinSettingsElectronic
        )
        self.jsettings_lattice = self.get_settings_object(text, JMinSettingsLattice)
        self.jsettings_ionic = self.get_settings_object(text, JMinSettingsIonic)

    def set_geomopt_vars(self, text: list[str]) -> None:
        """Set the geom_opt and geom_opt_type class variables.

        Set vars geom_opt and geom_opt_type for initializing self.jstrucs

        Parameters
        ----------
            text: list[str]
                output of read_file for out file
        """
        if self.jsettings_ionic is None or self.jsettings_lattice is None:
            self.set_min_settings(text)
        if self.jsettings_ionic is None or self.jsettings_lattice is None:
            raise ValueError("Unknown issue in setting settings objects")
        if self.jsettings_lattice.niterations > 0:
            self.geom_opt = True
            self.geom_opt_type = "lattice"
        elif self.jsettings_ionic.niterations > 0:
            self.geom_opt = True
            self.geom_opt_type = "ionic"
        else:
            self.geom_opt = False
            self.geom_opt_type = "single point"

    def set_jstrucs(self, text: list[str]) -> None:
        """Set the jstrucs class variable.

        Set the JStructures object to jstrucs from the out file text.

        Parameters
        ----------
            text: list[str]
                output of read_file for out file
        """
        self.jstrucs = JOutStructures.from_out_slice(text, iter_type=self.geom_opt_type)
        if self.etype is None:
            self.etype = self.jstrucs[-1].etype

    def set_orb_fillings(self) -> None:
        """Set the orbital fillings.

        Calculate and set homo and lumo fillings
        """
        if self.broadening_type is not None:
            self.homo_filling = (2 / self.nspin) * self.calculate_filling(
                self.broadening_type, self.broadening, self.homo, self.efermi
            )
            self.lumo_filling = (2 / self.nspin) * self.calculate_filling(
                self.broadening_type, self.broadening, self.lumo, self.efermi
            )
        else:
            self.homo_filling = 2 / self.nspin
            self.lumo_filling = 0

    def set_fluid(
        self, text: list[str]
    ) -> None:  # Is this redundant to the fluid settings?
        """Set the fluid class variable.

        Set the fluid class variable.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        line = find_first_range_key("fluid ", text)
        self.fluid = text[line[0]].split()[1]
        if self.fluid == "None":
            self.fluid = None

    def set_total_electrons(self, text: list[str]) -> None:
        """Set the total_Electrons class variable.

        Set the total_electrons class variable.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        total_electrons = self.jstrucs[-1].elecmindata[-1].nelectrons
        self.total_electrons = total_electrons

    def set_nbands(self, text: list[str]) -> None:
        """Set the Nbands class variable.

        Set the nbands class variable.
        Prioritizes finding nBands from the reiteration of the input parameters.
        If this line is not found, then it pulls it from "Setting up k-points,
        bands, fillings" section.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
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

    def set_atom_vars(self, text: list[str]) -> None:
        """Set the atom variables.

        Set the atom variables from the out file text.

        Parameters
        ----------
        text: list[str]
        output of read_file for out file
        """
        startline = find_key("Input parsed successfully", text)
        endline = find_key("---------- Initializing the Grid ----------", text)
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
        mapping_dict = dict(zip(atom_types, range(1, len(atom_types) + 1)))
        self.atom_elements_int = [mapping_dict[x] for x in self.atom_elements]
        self.atom_types = atom_types
        line = find_key("# Ionic positions in", text) + 1
        coords = np.array(
            [text[i].split()[2:5] for i in range(line, line + self.nat)], dtype=float
        )
        self.atom_coords_final = coords
        self.atom_coords = self.atom_coords_final.copy()

    def set_lattice_vars(self, text: list[str]) -> None:
        """Set the lattice variables.

        Set the lattice variables from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        self.lattice_initial = self.jstrucs[0].lattice.matrix
        self.lattice_final = self.jstrucs[-1].lattice.matrix
        self.lattice = self.lattice_final.copy()
        self.a, self.b, self.c = np.sum(self.lattice**2, axis=1) ** 0.5

    def set_ecomponents(self, text: list[str]) -> None:
        """Set the energy components dictionary.

        Set the energy components dictionary from the out file text.

        Parameters
        ----------
        text: list[str]
            output of read_file for out file
        """
        ecomp = self.jstrucs[-1].ecomponents
        if self.etype not in ecomp:
            ecomp[self.etype] = self.jstrucs[-1].E
        # line = find_key("# Energy components:", text)
        self.ecomponents = ecomp

    def calculate_filling(
        self, broadening_type: str, broadening: float, eig: float, efermi: float
    ) -> float:
        """Calculate the filling for a given eigenvalue.

        Use the broadening type, broadening value, eigenvalue, and fermi energy
        to calculate the filling at the eigenvalue.

        Parameters
        ----------
        broadening_type: str
            type of broadening to use
        broadening: float
            broadening parameter
        eig: float
            eigenvalue
        efermi: float
            fermi energy

        Returns
        -------
        filling: float
            filling at the eigenvalue
        """
        # most broadening implementations do not have the denominator factor
        # of 2, but JDFTx does currently.
        # Remove if use this for other code outfile reading
        x = (eig - efermi) / (2.0 * broadening)
        if broadening_type == "Fermi":
            filling = 0.5 * (1 - np.tanh(x))
        elif broadening_type == "Gauss":
            filling = 0.5 * (1 - math.erf(x))
        elif broadening_type == "MP1":
            filling = 0.5 * (1 - math.erf(x)) - x * np.exp(-1 * x**2) / (2 * np.pi**0.5)
        elif broadening_type == "Cold":
            filling = (
                0.5 * (1 - math.erf(x + 0.5**0.5))
                + np.exp(-1 * (x + 0.5**0.5) ** 2) / (2 * np.pi) ** 0.5
            )
        else:
            raise NotImplementedError("Have not added other broadening types")

        return filling

    def determine_is_metal(self) -> bool:
        """Determine if the system is a metal based.

        Determine if the system is a metal based on the fillings of
        homo and lumo.

        Returns
        -------
        is_metal: bool
            True if system is metallic
        """
        tol_partial = 0.01
        is_metal = True
        if (
            self.homo_filling / (2 / self.nspin) > (1 - tol_partial)
            and self.lumo_filling / (2 / self.nspin) < tol_partial
        ):
            is_metal = False
        return is_metal

    def check_solvation(self) -> bool:
        """Check for implicit solvation.

        Check if calculation used implicit solvation

        Returns
        -------
        has_solvation: bool
            True if calculation used implicit solvation
        """
        return self.fluid is not None

    def write(self) -> NotImplementedError:
        """Return an error.

        Return an error. (pre-commit needs a docustring here)
        """
        # don't need a write method since will never do that
        return NotImplementedError("There is no need to write a JDFTx out file")

    def to_dict(self) -> dict:
        """Convert dataclass to dictionary representation.

        Convert dataclass to dictionary representation.

        Returns
        -------
        dict
            JDFTXOutfileSlice in dictionary format
        """
        # convert dataclass to dictionary representation
        dct = {}
        for field in self.__dataclass_fields__:
            value = getattr(self, field)
            dct[field] = value
        return dct
