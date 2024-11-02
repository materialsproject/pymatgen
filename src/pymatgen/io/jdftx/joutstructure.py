"""Class object for storing a single JDFTx optimization step.

A mutant of the pymatgen Structure class for flexibility in holding JDFTx
"""

from __future__ import annotations

import inspect
from typing import Any, ClassVar

import numpy as np

from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pymatgen.io.jdftx.jelstep import JElSteps
from pymatgen.io.jdftx.utils import (
    _brkt_list_of_3x3_to_nparray,
    correct_geom_opt_type,
    get_colon_var_t1,
    is_charges_line,
    is_ecomp_start_line,
    is_forces_start_line,
    is_lattice_start_line,
    is_lowdin_start_line,
    is_magnetic_moments_line,
    is_posns_start_line,
    is_strain_start_line,
    is_stress_start_line,
)

__author__ = "Ben Rich"


class JOutStructure(Structure):
    """Class object for storing a single JDFTx optimization step.

    A mutant of the pymatgen Structure class for flexibility in holding JDFTx
    optimization data.
    """

    opt_type: str | None = None
    etype: str | None = None
    eopt_type: str | None = None
    emin_flag: str | None = None
    ecomponents: dict | None = None
    elecmindata: JElSteps | None = None
    stress: np.ndarray | None = None
    strain: np.ndarray | None = None
    nstep: int | None = None
    e: float | None = None
    grad_k: float | None = None
    alpha: float | None = None
    linmin: float | None = None
    t_s: float | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
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
    selective_dynamics: list[int] | None = None

    @property
    def mu(self) -> float | None:
        """Return the chemical potential.

        Return the chemical potential (Fermi level) in eV.

        Returns
        -------
        mu: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.mu
        return None

    @property
    def nelectrons(self) -> float | None:
        """Return the number of electrons.

        Return the total number of electrons in the electron density.

        Returns
        -------
        nelectrons: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.nelectrons
        return None

    @property
    def abs_magneticmoment(self) -> float | None:
        """Return the absolute magnetic moment.

        Return the absolute magnetic moment of the electron density.

        Returns
        -------
        abs_magneticmoment: float"""
        if self.elecmindata is not None:
            return self.elecmindata.abs_magneticmoment
        return None

    @property
    def tot_magneticmoment(self) -> float | None:
        """Return the total magnetic moment.

        Return the total magnetic moment of the electron density.

        Returns
        -------
        tot_magneticmoment: float"""
        if self.elecmindata is not None:
            return self.elecmindata.tot_magneticmoment
        return None

    @property
    def elec_nstep(self) -> int | None:
        """Return the most recent electronic step number.

        Return the nstep property of the electronic minimization data, where
        nstep corresponds to the SCF step number.

        Returns
        -------
        elec_nstep: int
        """
        if self.elecmindata is not None:
            return self.elecmindata.nstep
        return None

    @property
    def elec_e(self) -> int | None:
        """Return the most recent electronic energy.

        Return the e property of the electronic minimization, where e corresponds
        to the energy of the system's "etype".

        Returns
        -------
        elec_e: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.e
        return None

    @property
    def elec_grad_k(self) -> float | None:
        """Return the most recent electronic grad_k.

        Return the most recent electronic grad_k, where grad_k here corresponds
        to the gradient of the electronic density along the most recent line minimization.

        Returns
        -------
        grad_k: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.grad_k
        return None

    @property
    def elec_alpha(self) -> float | None:
        """Return the most recent electronic alpha.

        Return the most recent electronic alpha, where alpha here corresponds to the
        step size of the most recent SCF step along the line minimization.

        Returns
        -------
        alpha: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.alpha
        return None

    @property
    def elec_linmin(self) -> float | None:
        """Return the most recent electronic linmin.

        Return the most recent electronic linmin, where linmin here corresponds to
        the normalized alignment projection of the electronic energy gradient to
        the line minimization direction. (-1 perfectly anti-aligned, 0 orthogonal, 1 perfectly aligned)

        Returns
        -------
        linmin: float
        """
        if self.elecmindata is not None:
            return self.elecmindata.linmin
        return None

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
        eopt_type: str = "ElecMinimize",
        opt_type: str = "IonicMinimize",
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
        eopt_type: str
            The type of electronic minimization step
        opt_type: str
            The type of optimization step
        emin_flag: str
            The flag that indicates the start of a log message for a JDFTx
            optimization step
        """
        instance = cls(lattice=np.eye(3), species=[], coords=[], site_properties={})
        if opt_type not in ["IonicMinimize", "LatticeMinimize"]:
            opt_type = correct_geom_opt_type(opt_type)
        instance.eopt_type = eopt_type
        instance.opt_type = opt_type
        instance.emin_flag = emin_flag
        line_collections = instance.init_line_collections()
        line_collections = instance.gather_line_collections(line_collections, text_slice)

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

        # In case of single-point calculation
        if instance.e is None:  # This doesn't defer to elecmindata.e due to the existence of a class variable e
            if instance.etype is not None:
                if instance.ecomponents is not None:
                    if instance.etype in instance.ecomponents:
                        instance.e = instance.ecomponents[instance.etype]
                    elif instance.elecmindata is not None:
                        instance.e = instance.elecmindata.e
                    else:
                        raise ValueError("Could not determine total energy due to lack of elecmindata")
                else:
                    raise ValueError("Could not determine total energy due to lack of ecomponents")
            else:
                raise ValueError("Could not determine total energy due to lack of etype")

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

    def gather_line_collections(self, line_collections: dict, text_slice: list[str]) -> dict:
        """Gather line collections.

        Gather lines of text from a JDFTx out file into a dictionary of line collections.

        Parameters
        ----------
        line_collections: dict
            A dictionary of line collections for each type of line in a JDFTx
            out file. Assumed pre-initialized (there exists sdict["lines"]: list[str],
            sdict["collecting"]: bool, sdict["collected"]: bool for every sdict = line_collections[line_type])

        text_slice: list[str]
            A slice of text from a JDFTx out file corresponding to a single
            optimization step / SCF cycle
        """
        for line in text_slice:
            read_line = False
            for line_type in line_collections:
                sdict = line_collections[line_type]
                if sdict["collecting"]:
                    lines, getting, got = self.collect_generic_line(line, sdict["lines"])
                    sdict["lines"] = lines
                    sdict["collecting"] = getting
                    sdict["collected"] = got
                    read_line = True
                    break
            if not read_line:
                for line_type in line_collections:
                    if (not line_collections[line_type]["collected"]) and self.is_generic_start_line(line, line_type):
                        line_collections[line_type]["collecting"] = True
                        line_collections[line_type]["lines"].append(line)
                        break
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
        if self.emin_flag is None:
            raise ValueError("emin_flag is not set")
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
        is_line = f"{self.opt_type}:" in line_text
        return is_line and "Iter:" in line_text

    def get_etype_from_emin_lines(self, emin_lines: list[str]) -> str | None:
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
        # If not F or G, most likely given as Etot
        if etype is None:
            for line in emin_lines:
                if "Etot:" in line:
                    etype = "Etot"
                    break
        # Used as last-case scenario as the ordering of <etype> after <Iter> needs
        # to be checked by reading through the source code (TODO)
        if etype is None:
            for line in emin_lines:
                if "Iter:" in line:
                    # Assume line will have etype in "... Iter: n <etype>: num ..."
                    etype = line.split("Iter:")[1].split(":")[0].strip().split()[-1].strip()
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
            self.elecmindata = JElSteps.from_text_slice(emin_lines, opt_type=self.eopt_type, etype=self.etype)

    def parse_lattice_lines(self, lattice_lines: list[str]) -> None:
        """Parse lattice lines.

        Parse the lines of text corresponding to the lattice vectors of a
        JDFTx out file. Collects the lattice matrix "r" as a 3x3 numpy array first
        in column-major order (vec i = r[:,i]), then transposes it to row-major
        order (vec i = r[i,:]) and converts from Bohr to Angstroms.

        Parameters
        ----------
        lattice_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            lattice vectors
        """
        r = None
        if len(lattice_lines):
            r = _brkt_list_of_3x3_to_nparray(lattice_lines, i_start=2)
            r = r.T * bohr_to_ang
            self.lattice = Lattice(r)

    def parse_strain_lines(self, strain_lines: list[str]) -> None:
        """Parse strain lines.

        Parse the lines of text corresponding to the strain tensor of a
        JDFTx out file. Converts from column-major to row-major order.

        Parameters
        ----------
        strain_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            strain tensor
        """
        # TODO: Strain is a unitless quantity, so column to row-major conversion
        # should cover all unit conversion. Double check if this is true.
        st = None
        if len(strain_lines):
            st = _brkt_list_of_3x3_to_nparray(strain_lines, i_start=1)
            st = st.T
        self.strain = st

    def parse_stress_lines(self, stress_lines: list[str]) -> None:
        """Parse stress lines.

        Parse the lines of text corresponding to the stress tensor of a
        JDFTx out file and converts from Ha/Bohr^3 to eV/Ang^3.

        Parameters
        ----------
        stress_lines: list[str]
            A list of lines of text from a JDFTx out file containing the
            stress tensor
        """
        # TODO: Lattice optimizations dump stress in cartesian coordinates in units
        # "[Eh/a0^3]" (Hartree per bohr cubed). Check if this changes for direct
        # coordinates.
        st = None
        if len(stress_lines):
            st = _brkt_list_of_3x3_to_nparray(stress_lines, i_start=1)
            st = st.T
            st *= Ha_to_eV / (bohr_to_ang**3)
        self.stress = st

    def parse_posns_lines(self, posns_lines: list[str]) -> None:
        """Parse positions lines.

        Parse the lines of text corresponding to the positions of a
        JDFTx out file

        Parameters
        ----------
        posns_lines: list[str]
            A list of lines of text from a JDFTx out file. Collected lines will
            start with the string "# Ionic positions in ..." (specifying either
            cartesian or direct coordinates), followed by a line for each ion (atom)
            in the format "ion_name x y z sd", where ion_name is the name of the element,
            and sd is a flag indicating whether the ion is excluded from optimization (1)
             or not (0).
        """
        if len(posns_lines):
            natoms = len(posns_lines) - 1
            coords_type = posns_lines[0].split("positions in")[1]
            coords_type = coords_type.strip().split()[0].strip()
            posns: list[np.ndarray] = []
            names: list[str] = []
            selective_dynamics: list[int] = []
            for i in range(natoms):
                line = posns_lines[i + 1]
                name = line.split()[1].strip()
                posn = np.array([float(x.strip()) for x in line.split()[2:5]])
                sd = int(line.split()[5])
                names.append(name)
                posns.append(posn)
                selective_dynamics.append(sd)
            posns = np.array(posns)
            if coords_type.lower() != "cartesian":
                posns = np.dot(posns, self.lattice.matrix)
            else:
                posns *= bohr_to_ang
            for i in range(natoms):
                self.append(species=names[i], coords=posns[i], coords_are_cartesian=True)
            self.selective_dynamics = selective_dynamics

    def parse_forces_lines(self, forces_lines: list[str]) -> None:
        """Parse forces lines.

        Parse the lines of text corresponding to the forces of a
        JDFTx out file.

        Parameters
        ----------
        forces_lines: list[str]
            A list of lines of text from a JDFTx out file containing the forces
        """
        if len(forces_lines):
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
                # TODO: Double check conversion of forces from direct to cartesian
                forces = np.dot(forces, self.lattice.matrix)
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
            A list of lines of text from a JDFTx out file. All lines will either be
            the header line, a break line of only "...---...", or a line of the form
            "component = value" where component is the name of the energy component
            and value is the value of the energy component in Hartrees.
        """
        self.ecomponents = {}
        key = None
        for line in ecomp_lines:
            if " = " in line:
                lsplit = line.split(" = ")
                key = lsplit[0].strip()
                val = float(lsplit[1].strip())
                self.ecomponents[key] = val * Ha_to_eV
        if key is not None and (self.etype is None) and (key in ["F", "G"]):
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
            elif is_magnetic_moments_line(line):
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

    def parse_lowdin_line(self, lowdin_line: str, lowdin_dict: dict[str, list[float]]) -> dict[str, list[float]]:
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
        return f"{self.opt_type}: Converged" in line_text

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
                    nstep = int(get_colon_var_t1(line, "Iter:"))
                    self.nstep = nstep
                    en = get_colon_var_t1(line, f"{self.etype}:")
                    self.E = en * Ha_to_eV
                    grad_k = get_colon_var_t1(line, "|grad|_K: ")
                    self.grad_k = grad_k
                    alpha = get_colon_var_t1(line, "alpha: ")
                    self.alpha = alpha
                    linmin = get_colon_var_t1(line, "linmin: ")
                    self.linmin = linmin
                    t_s = get_colon_var_t1(line, "t[s]: ")
                    self.t_s = t_s
                elif self.is_opt_conv_line(line):
                    self.geom_converged = True
                    self.geom_converged_reason = line.split("(")[1].split(")")[0].strip()

    def is_generic_start_line(self, line_text: str, line_type: str) -> bool:
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

    def collect_generic_line(self, line_text: str, generic_lines: list[str]) -> tuple[list[str], bool, bool]:
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

    # This method is likely never going to be called as all (currently existing)
    # attributes of the most recent slice are explicitly defined as a class
    # property. However, it is included to reduce the likelihood of errors
    # upon future changes to downstream code.
    def __getattr__(self, name: str) -> Any:
        """Return attribute value.

        Return the value of an attribute.

        Parameters
        ----------
        name: str
            The name of the attribute

        Returns
        -------
        value
            The value of the attribute
        """
        # Only works for actual attributes of the class
        if name in self.__dict__:
            return self.__dict__[name]

        # Extended for properties
        for cls in inspect.getmro(self.__class__):
            if name in cls.__dict__ and isinstance(cls.__dict__[name], property):
                return cls.__dict__[name].__get__(self)

        # Check if the attribute is in self.jstrucs
        if hasattr(self.elecmindata, name):
            return getattr(self.elecmindata, name)

        # If the attribute is not found in either, raise an AttributeError
        raise AttributeError(f"{self.__class__.__name__} not found: {name}")
