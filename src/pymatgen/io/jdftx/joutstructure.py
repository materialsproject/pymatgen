"""Class object for storing a single JDFTx optimization step.

A mutant of the pymatgen Structure class for flexibility in holding JDFTx
"""

from __future__ import annotations

from typing import ClassVar, TypeVar

import numpy as np
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.units import Ha_to_eV, bohr_to_ang

from atomate2.jdftx.io.jeiters import JEiters
from atomate2.jdftx.io.joutstructure_helpers import (
    _get_colon_var_t1,
    correct_iter_type,
    is_charges_line,
    is_ecomp_start_line,
    is_forces_start_line,
    is_lattice_start_line,
    is_lowdin_start_line,
    is_moments_line,
    is_posns_start_line,
    is_strain_start_line,
    is_stress_start_line,
)

T = TypeVar("T", bound="JOutStructure")


class JOutStructure(Structure):
    """Class object for storing a single JDFTx optimization step.

    A mutant of the pymatgen Structure class for flexibility in holding JDFTx
    optimization data.
    """

    iter_type: str = None
    etype: str = None
    eiter_type: str = None
    emin_flag: str = None
    ecomponents: dict = None
    elecmindata: JEiters = None
    stress: np.ndarray = None
    strain: np.ndarray = None
    iter: int = None
    E: float = None
    grad_k: float = None
    alpha: float = None
    linmin: float = None
    t_s: float = None
    geom_converged: bool = False
    geom_converged_reason: str = None
    line_types: ClassVar[list[str]] = [
        "emin",
        "lattice",
        "strain",
        "stress",
        "posns",
        "forces",
        "ecomp",
        "lowdin",
        "opt",
    ]

    def __init__(
        self,
        lattice: np.ndarray,
        species: list[str],
        coords: list[np.ndarray],
        site_properties: dict[str, list],
    ) -> None:
        super().__init__(
            lattice=lattice,
            species=species,
            coords=coords,
            site_properties=site_properties,
        )

    @classmethod
    def from_text_slice(
        cls,
        text_slice: list[str],
        eiter_type: str = "ElecMinimize",
        iter_type: str = "IonicMinimize",
        emin_flag: str = "---- Electronic minimization -------",
    ) -> JOutStructure:
        """Return JOutStructure object.

        Create a JAtoms object from a slice of an out file's text corresponding
        to a single step of a native JDFTx optimization.

        Parameters
        ----------
        text_slice: list[str]
            A slice of text from a JDFTx out file corresponding to a single
            optimization step / SCF cycle
        eiter_type: str
            The type of electronic minimization step
        iter_type: str
            The type of optimization step
        emin_flag: str
            The flag that indicates the start of a log message for a JDFTx
            optimization step
        """
        instance = cls(lattice=np.eye(3), species=[], coords=[], site_properties={})
        if iter_type not in ["IonicMinimize", "LatticeMinimize"]:
            iter_type = correct_iter_type(iter_type)
        instance.eiter_type = eiter_type
        instance.iter_type = iter_type
        instance.emin_flag = emin_flag
        line_collections = instance.init_line_collections()
        for line in text_slice:
            read_line = False
            for line_type in line_collections:
                sdict = line_collections[line_type]
                if sdict["collecting"]:
                    lines, getting, got = instance.collect_generic_line(
                        line, sdict["lines"]
                    )
                    sdict["lines"] = lines
                    sdict["collecting"] = getting
                    sdict["collected"] = got
                    read_line = True
                    break
            if not read_line:
                for line_type in line_collections:
                    if (
                        not line_collections[line_type]["collected"]
                    ) and instance.is_generic_start_line(line, line_type):
                        line_collections[line_type]["collecting"] = True
                        line_collections[line_type]["lines"].append(line)
                        break

        # ecomponents needs to be parsed before emin to set etype
        instance.parse_ecomp_lines(line_collections["ecomp"]["lines"])
        instance.parse_emin_lines(line_collections["emin"]["lines"])
        # Lattice must be parsed before posns/forces in case of direct
        # coordinates
        instance.parse_lattice_lines(line_collections["lattice"]["lines"])
        instance.parse_posns_lines(line_collections["posns"]["lines"])
        instance.parse_forces_lines(line_collections["forces"]["lines"])
        # Strain and stress can be parsed in any order
        instance.parse_strain_lines(line_collections["strain"]["lines"])
        instance.parse_stress_lines(line_collections["stress"]["lines"])
        # Lowdin must be parsed after posns
        instance.parse_lowdin_lines(line_collections["lowdin"]["lines"])
        # Opt line must be parsed after ecomp
        instance.parse_opt_lines(line_collections["opt"]["lines"])

        return instance

    def init_line_collections(self) -> dict:
        """Initialize line collection dict.

        Initialize a dictionary of line collections for each type of line in a
        JDFTx out file.

        Returns
        -------
        line_collections: dict
            A dictionary of line collections for each type of line in a JDFTx
            out file
        """
        line_collections = {}
        for line_type in self.line_types:
            line_collections[line_type] = {
                "lines": [],
                "collecting": False,
                "collected": False,
            }
        return line_collections

    def is_emin_start_line(self, line_text: str) -> bool:
        """Return True if emin start line.

        Return True if the line_text is the start of a log message for a JDFTx
        optimization step.

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
        """
        return self.emin_flag in line_text

    def is_opt_start_line(self, line_text: str) -> bool:
        """Return True if opt start line.

        Return True if the line_text is the start of a log message for a JDFTx
        optimization step

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a log message for a JDFTx
            optimization step
        """
        is_line = f"{self.iter_type}:" in line_text
        return is_line and "Iter:" in line_text

    def get_etype_from_emin_lines(self, emin_lines: list[str]) -> str:
        """Return energy type string.

        Return the type of energy from the electronic minimization data of a
        JDFTx out file.

        Parameters
        ----------
        emin_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            electronic minimization data

        Returns
        -------
        etype: str
            The type of energy from the electronic minimization data of a JDFTx
            out file
        """
        etype = None
        for line in emin_lines:
            if "F:" in line:
                etype = "F"
                break
            if "G:" in line:
                etype = "G"
                break
        return etype

    def set_etype_from_emin_lines(self, emin_lines: list[str]) -> None:
        """Set etype class variable.

        Set the type of energy from the electronic minimization data of a
        JDFTx out file.

        Parameters
        ----------
        emin_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            electronic minimization data
        """
        self.etype = self.get_etype_from_emin_lines(emin_lines)
        if self.etype is None:
            raise ValueError(
                "Could not determine energy type from electronic minimization \
                    data"
            )

    def parse_emin_lines(self, emin_lines: list[str]) -> None:
        """Parse electronic minimization lines.

        Parse the lines of text corresponding to the electronic minimization
        data of a JDFTx out file.

        Parameters
        ----------
        emin_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            electronic minimization data
        """
        if len(emin_lines):
            if self.etype is None:
                self.set_etype_from_emin_lines(emin_lines)
            self.elecmindata = JEiters.from_text_slice(
                emin_lines, iter_type=self.eiter_type, etype=self.etype
            )

    def parse_lattice_lines(self, lattice_lines: list[str]) -> None:
        """Parse lattice lines.

        Parse the lines of text corresponding to the lattice vectors of a
        JDFTx out file.

        Parameters
        ----------
        lattice_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            lattice vectors
        """
        r = None
        if len(lattice_lines):
            r = self._brkt_list_of_3x3_to_nparray(lattice_lines, i_start=2)
            r = r.T * bohr_to_ang
            self.lattice = Lattice(r)

    def parse_strain_lines(self, strain_lines: list[str]) -> None:
        """Parse strain lines.

        Parse the lines of text corresponding to the strain tensor of a
        JDFTx out file.

        Parameters
        ----------
        strain_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            strain tensor
        """
        st = None
        if len(strain_lines):
            st = self._brkt_list_of_3x3_to_nparray(strain_lines, i_start=1)
            st = st.T * 1  # Conversion factor?
        self.strain = st

    def parse_stress_lines(self, stress_lines: list[str]) -> None:
        """Parse stress lines.

        Parse the lines of text corresponding to the stress tensor of a
        JDFTx out file.

        Parameters
        ----------
        stress_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            stress tensor
        """
        st = None
        if len(stress_lines):
            st = self._brkt_list_of_3x3_to_nparray(stress_lines, i_start=1)
            st = st.T * 1  # Conversion factor?
        self.stress = st

    def parse_posns_lines(self, posns_lines: list[str]) -> None:
        """Parse positions lines.

        Parse the lines of text corresponding to the positions of a
        JDFTx out file

        Parameters
        ----------
        posns_lines: list[str]
            A list of lines of text from a JDFTx out file
        """
        natoms = len(posns_lines) - 1
        coords_type = posns_lines[0].split("positions in")[1]
        coords_type = coords_type.strip().split()[0].strip()
        posns = []
        names = []
        for i in range(natoms):
            line = posns_lines[i + 1]
            name = line.split()[1].strip()
            posn = np.array([float(x.strip()) for x in line.split()[2:5]])
            names.append(name)
            posns.append(posn)
        posns = np.array(posns)
        if coords_type.lower() != "cartesian":
            posns = np.dot(posns, self.lattice.matrix)
        else:
            posns *= bohr_to_ang
        for i in range(natoms):
            self.append(species=names[i], coords=posns[i], coords_are_cartesian=True)

    def parse_forces_lines(self, forces_lines: list[str]) -> None:
        """Parse forces lines.

        Parse the lines of text corresponding to the forces of a
        JDFTx out file.

        Parameters
        ----------
        forces_lines: list[str]
            A list of lines of text from a JDFTx out file containing the forces
        """
        natoms = len(forces_lines) - 1
        coords_type = forces_lines[0].split("Forces in")[1]
        coords_type = coords_type.strip().split()[0].strip()
        forces = []
        for i in range(natoms):
            line = forces_lines[i + 1]
            force = np.array([float(x.strip()) for x in line.split()[2:5]])
            forces.append(force)
        forces = np.array(forces)
        if coords_type.lower() != "cartesian":
            forces = np.dot(
                forces, self.lattice.matrix
            )  # TODO: Double check this conversion
            # (since self.cell is in Ang, and I need the forces in eV/ang, how
            # would you convert forces from direct coordinates into cartesian?)
        else:
            forces *= 1 / bohr_to_ang
        forces *= Ha_to_eV
        self.forces = forces

    def parse_ecomp_lines(self, ecomp_lines: list[str]) -> None:
        """Parse energy component lines.

        Parse the lines of text corresponding to the energy components of a
        JDFTx out file

        Parameters
        ----------
        ecomp_lines: list[str]
            A list of lines of text from a JDFTx out file
        """
        self.ecomponents = {}
        for line in ecomp_lines:
            if " = " in line:
                lsplit = line.split(" = ")
                key = lsplit[0].strip()
                val = float(lsplit[1].strip())
                self.ecomponents[key] = val * Ha_to_eV
        if (self.etype is None) and (key in ["F", "G"]):
            self.etype = key

    def parse_lowdin_lines(self, lowdin_lines: list[str]) -> None:
        """Parse Lowdin lines.

        Parse the lines of text corresponding to a Lowdin population analysis
        in a JDFTx out file

        Parameters
        ----------
        lowdin_lines: list[str]
            A list of lines of text from a JDFTx out file
        """
        charges_dict: dict[str, list[float]] = {}
        moments_dict: dict[str, list[float]] = {}
        for line in lowdin_lines:
            if is_charges_line(line):
                charges_dict = self.parse_lowdin_line(line, charges_dict)
            elif is_moments_line(line):
                moments_dict = self.parse_lowdin_line(line, moments_dict)
        names = [s.name for s in self.species]
        charges = None
        moments = None
        if len(charges_dict):
            charges = np.zeros(len(names))
            for el in charges_dict:
                idcs = [int(i) for i in range(len(names)) if names[i] == el]
                for i, idx in enumerate(idcs):
                    charges[idx] += charges_dict[el][i]
        if len(moments_dict):
            moments = np.zeros(len(names))
            for el in moments_dict:
                idcs = [i for i in range(len(names)) if names[i] == el]
                for i, idx in enumerate(idcs):
                    moments[idx] += moments_dict[el][i]
        self.charges = charges
        self.magnetic_moments = moments

    def parse_lowdin_line(
        self, lowdin_line: str, lowdin_dict: dict[str, list[float]]
    ) -> dict[str, list[float]]:
        """Parse Lowdin line.

        Parse a line of text from a JDFTx out file corresponding to a
        Lowdin population analysis

        Parameters
        ----------
        lowdin_line: str
            A line of text from a JDFTx out file
        lowdin_dict: dict[str, list[float]]
            A dictionary of Lowdin population analysis data

        Returns
        -------
        lowdin_dict: dict[str, float]
            A dictionary of Lowdin population analysis data
        """
        tokens = [v.strip() for v in lowdin_line.strip().split()]
        name = tokens[2]
        vals = [float(x) for x in tokens[3:]]
        lowdin_dict[name] = vals
        return lowdin_dict

    def is_opt_conv_line(self, line_text: str) -> bool:
        """Return True if line_text is geom opt convergence line.

        Return True if the line_text is the end of a JDFTx optimization step

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file

        Returns
        -------
        is_line: bool
            True if the line_text is the end of a JDFTx optimization step
        """
        return f"{self.iter_type}: Converged" in line_text

    def parse_opt_lines(self, opt_lines: list[str]) -> None:
        """Parse optimization lines.

        Parse the lines of text corresponding to the optimization step of a
        JDFTx out file

        Parameters
        ----------
        opt_lines: list[str]
            A list of lines of text from a JDFTx out file
        """
        if len(opt_lines):
            for line in opt_lines:
                if self.is_opt_start_line(line):
                    n_iter = int(_get_colon_var_t1(line, "Iter:"))
                    self.iter = n_iter
                    en = _get_colon_var_t1(line, f"{self.etype}:")
                    self.E = en * Ha_to_eV
                    grad_k = _get_colon_var_t1(line, "|grad|_K: ")
                    self.grad_k = grad_k
                    alpha = _get_colon_var_t1(line, "alpha: ")
                    self.alpha = alpha
                    linmin = _get_colon_var_t1(line, "linmin: ")
                    self.linmin = linmin
                    t_s = _get_colon_var_t1(line, "t[s]: ")
                    self.t_s = t_s
                elif self.is_opt_conv_line(line):
                    self.geom_converged = True
                    self.geom_converged_reason = (
                        line.split("(")[1].split(")")[0].strip()
                    )

    def is_generic_start_line(self, line_text: str, line_type: str) -> bool:
        # I am choosing to map line_type to a function this way because
        # I've had horrible experiences with storing functions in dictionaries
        # in the past
        """Return True if the line_text is start of line_type log message.

        Return True if the line_text is the start of a section of the
        JDFTx out file corresponding to the line_type.

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file
        line_type: str
            The type of line to check for

        Returns
        -------
        is_line: bool
            True if the line_text is the start of a section of the
            JDFTx out file
        """
        if line_type == "lowdin":
            return is_lowdin_start_line(line_text)
        if line_type == "opt":
            return self.is_opt_start_line(line_text)
        if line_type == "ecomp":
            return is_ecomp_start_line(line_text)
        if line_type == "forces":
            return is_forces_start_line(line_text)
        if line_type == "posns":
            return is_posns_start_line(line_text)
        if line_type == "stress":
            return is_stress_start_line(line_text)
        if line_type == "strain":
            return is_strain_start_line(line_text)
        if line_type == "lattice":
            return is_lattice_start_line(line_text)
        if line_type == "emin":
            return self.is_emin_start_line(line_text)
        raise ValueError(f"Unrecognized line type {line_type}")

    def collect_generic_line(
        self, line_text: str, generic_lines: list[str]
    ) -> tuple[list[str], bool, bool]:
        """Collect generic log line.

        Collect a line of text into a list of lines if the line is not empty,
        and otherwise updates the collecting and collected flags.

        Parameters
        ----------
        line_text: str
            A line of text from a JDFTx out file
        generic_lines: list[str]
            A list of lines of text of the same type

        Returns
        -------
        generic_lines: list[str]
            A list of lines of text of the same type
        collecting: bool
            True if the line_text is not empty
        collected: bool
            True if the line_text is empty (end of section)
        """
        collecting = True
        collected = False
        if not len(line_text.strip()):
            collecting = False
            collected = True
        else:
            generic_lines.append(line_text)
        return generic_lines, collecting, collected

    def _brkt_list_of_3_to_nparray(self, line: str) -> np.ndarray:
        """Return 3x1 numpy array.

        Convert a string of the form "[ x y z ]" to a 3x1 numpy array

        Parameters
        ----------
        line: str
            A string of the form "[ x y z ]"
        """
        return np.array([float(x) for x in line.split()[1:-1]])

    def _brkt_list_of_3x3_to_nparray(
        self, lines: list[str], i_start: int = 0
    ) -> np.ndarray:
        """Return 3x3 numpy array.

        Convert a list of strings of the form "[ x y z ]" to a 3x3 numpy array

        Parameters
        ----------
        lines: list[str]
            A list of strings of the form "[ x y z ]"
        i_start: int
            The index of the first line in lines

        Returns
        -------
        out: np.ndarray
            A 3x3 numpy array
        """
        out = np.zeros([3, 3])
        for i in range(3):
            out[i, :] += self._brkt_list_of_3_to_nparray(lines[i + i_start])
        return out
