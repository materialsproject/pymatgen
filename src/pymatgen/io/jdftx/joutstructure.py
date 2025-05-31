"""Class object for storing a single JDFTx geometric optimization step.

A mutant of the pymatgen Structure class for flexibility in holding JDFTx.
"""

from __future__ import annotations

import pprint
import warnings
from typing import TYPE_CHECKING, ClassVar

import numpy as np
from monty.dev import deprecated

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import ArrayLike, NDArray

    from pymatgen.util.typing import CompositionLike

from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.units import Ha_to_eV, bohr_to_ang
from pymatgen.io.jdftx._output_utils import (
    _brkt_list_of_3x3_to_nparray,
    correct_geom_opt_type,
    get_colon_val,
    is_lowdin_start_line,
)
from pymatgen.io.jdftx.jelstep import JElSteps

__author__ = "Ben Rich"

_jos_atrs_from_elecmindata = ("mu", "nelectrons", "abs_magneticmoment", "tot_magneticmoment")
_jos_atrs_elec_from_elecmindata = ("nstep", "e", "grad_k", "alpha", "linmin")


def _is_homogenous(val: list):
    try:
        np.array(val)
        return True
    except ValueError:
        return False


# TODO: Phase out redundant class attributes to data stored in `properties` and `site_properties` dicts.
class JOutStructure(Structure):
    """Class object for storing a single JDFTx optimization step.

    A mutant of the pymatgen Structure class for flexibility in holding JDFTx
    optimization data.

    Properties:
        charges (np.ndarray | None): The Lowdin charges of the atoms in the system.
        magnetic_moments (np.ndarray | None): The magnetic moments of the atoms in the system.
    Attributes:
        opt_type (str | None): The type of optimization step.
        etype (str | None): The type of energy from the electronic minimization data.
        eopt_type (str | None): The type of electronic minimization step.
        emin_flag (str | None): The flag that indicates the start of a log message for a JDFTx optimization step.
        ecomponents (dict | None): The energy components of the system.
        elecmindata (JElSteps | None): The electronic minimization data.
        stress (np.ndarray | None): The stress tensor.
        strain (np.ndarray | None): The strain tensor.
        nstep (int | None): The most recent step number.
        e (float | None): The total energy of the system.
        grad_k (float | None): The gradient of the electronic density along the most recent line minimization.
        alpha (float | None): The step size of the most recent SCF step along the line minimization.
        linmin (float | None): The normalized alignment projection of the electronic energy gradient to the line
                                minimization direction.
        t_s (float | None): The time in seconds for the optimization step.
        geom_converged (bool): Whether the geometry optimization has converged.
        geom_converged_reason (str | None): The reason for geometry optimization convergence.
        line_types (ClassVar[list[str]]): The types of lines in a JDFTx out file.
        selective_dynamics (list[int] | None): The selective dynamics flags for the atoms in the system.
        mu (float | None): The chemical potential (Fermi level) in eV.
        nelectrons (float | None): The total number of electrons in the electron density.
        abs_magneticmoment (float | None): The absolute magnetic moment of the electron density.
        tot_magneticmoment (float | None): The total magnetic moment of the electron density.
        elec_nstep (int | None): The most recent electronic step number.
        elec_e (float | None): The most recent electronic energy.
        elec_grad_k (float | None): The most recent electronic grad_k.
        elec_alpha (float | None): The most recent electronic alpha.
        elec_linmin (float | None): The most recent electronic linmin.
        structure (Structure | None): The Structure object of the system. (helpful for uses where the JOutStructure
                                      metadata causes issues)
        is_md (bool): Whether the optimization step is a molecular dynamics step.
    """

    opt_type: str | None = None
    etype: str | None = None
    eopt_type: str | None = None
    emin_flag: str | None = None
    ecomponents: dict | None = None
    elecmindata: JElSteps | None = None
    stress: NDArray[np.float64] | None = None
    kinetic_stress: NDArray[np.float64] | None = None
    strain: NDArray[np.float64] | None = None
    forces: NDArray[np.float64] | None = None
    nstep: int | None = None
    e: float | np.float64 | None = None
    grad_k: float | np.float64 | None = None
    alpha: float | np.float64 | None = None
    linmin: float | np.float64 | None = None
    t_s: float | np.float64 | None = None
    geom_converged: bool = False
    geom_converged_reason: str | None = None
    line_types: ClassVar[list[str]] = [
        "emin",
        "lattice",
        "strain",
        "kinetic_stress",
        "stress",
        "posns",
        "forces",
        "ecomp",
        "lowdin",
        "opt",
    ]
    mu: float | np.float64 | None = None
    nelectrons: float | np.float64 | None = None
    abs_magneticmoment: float | np.float64 | None = None
    tot_magneticmoment: float | np.float64 | None = None
    elec_nstep: int | None = None
    elec_e: float | np.float64 | None = None
    elec_grad_k: float | np.float64 | None = None
    elec_alpha: float | np.float64 | None = None
    elec_linmin: float | np.float64 | None = None
    structure: Structure | None = None
    is_md: bool = False
    thermostat_velocity: np.ndarray | None = None
    _velocities: list[NDArray[np.float64] | None] | None = None
    _constraint_vectors: list[NDArray[np.float64] | list[NDArray[np.float64]] | None] | None = None
    _constraint_types: list[str | None] | None = None
    _group_names: list[list[str] | None] | None = None

    def _elecmindata_postinit(self) -> None:
        """Post-initialization method for attributes taken from elecmindata."""
        if self.elecmindata is not None:
            for var in _jos_atrs_from_elecmindata:
                if hasattr(self.elecmindata, var):
                    setattr(self, var, getattr(self.elecmindata, var))
            for var in _jos_atrs_elec_from_elecmindata:
                if hasattr(self.elecmindata, var):
                    setattr(self, f"elec_{var}", getattr(self.elecmindata, var))
        else:
            for var in _jos_atrs_from_elecmindata:
                setattr(self, var, None)
            for var in _jos_atrs_elec_from_elecmindata:
                setattr(self, f"elec_{var}", None)

    @property
    def charges(self) -> NDArray[np.float64] | None:
        """Return the Lowdin charges.

        Returns:
            np.ndarray: The Lowdin charges of the atoms in the system.
        """
        if "charges" not in self.site_properties:
            return None
        return np.array(self.site_properties["charges"])

    @charges.setter
    def charges(self, charges: NDArray[np.float64] | None) -> None:
        """Set the Lowdin charges.

        Args:
            charges (np.ndarray): The Lowdin charges of the atoms in the system.
        """
        if charges is not None:
            self.add_site_property("charges", list(charges))
        elif "charges" in self.site_properties:
            self.remove_site_property("charges")

    @property
    def magnetic_moments(self) -> NDArray[np.float64] | None:
        """Return the magnetic moments.

        Returns:
            np.ndarray: The magnetic moments of the atoms in the system.
        """
        if "magmom" not in self.site_properties:
            return None
        return np.array(self.site_properties["magmom"])

    @magnetic_moments.setter
    def magnetic_moments(self, magnetic_moments: NDArray[np.float64]) -> None:
        """Set the magnetic moments.

        Args:
            magnetic_moments (np.ndarray): The magnetic moments of the atoms in the system.
        """
        if magnetic_moments is not None:
            self.add_site_property("magmom", list(magnetic_moments))
        elif "magmom" in self.site_properties:
            self.remove_site_property("magmom")

    @property
    def selective_dynamics(self) -> Sequence | None:
        """Return the selective dynamics.

        Returns:
            list[list[int]]: The selective dynamics of the atoms in the system.
        """
        if "selective_dynamics" not in self.site_properties:
            return None
        return self.site_properties["selective_dynamics"]

    @selective_dynamics.setter
    def selective_dynamics(self, selective_dynamics: np.ndarray | list[list[int]]) -> None:
        """Set the selective dynamics.

        Args:
            selective_dynamics (np.ndarray | list[list[int]]): The selective dynamics of the atoms in the system.
        """
        if selective_dynamics is not None:
            self.add_site_property("selective_dynamics", list(selective_dynamics))
        elif "selective_dynamics" in self.site_properties:
            self.remove_site_property("selective_dynamics")

    @property
    def velocities(self) -> Sequence | list[np.ndarray | None] | None:
        """Return the velocities.

        Returns:
            list[np.ndarray | None] | None: The velocities of the atoms in the system.
        """
        if "velocities" not in self.site_properties:
            return self._velocities
        return self.site_properties["velocities"]

    @velocities.setter
    def velocities(self, velocities: list[np.ndarray | None] | None) -> None:
        """Set the velocities.

        Args:
            velocities (list[np.ndarray | None] | None): The velocities of the atoms in the system.
        """
        self._velocities = velocities
        if (velocities is not None) and _is_homogenous(velocities):
            self.add_site_property("velocities", list(velocities))
        elif "velocities" in self.site_properties:
            self.remove_site_property("velocities")

    @property
    def constraint_types(self) -> Sequence | list[str | None] | None:
        """Return the constraint_types.

        Returns:
            list[str | None] | None: The constraint_types of the atoms in the system.
        """
        if "constraint_types" not in self.site_properties:
            return self._constraint_types
        return self.site_properties["constraint_types"]

    @constraint_types.setter
    def constraint_types(self, constraint_types: list[str | None] | None) -> None:
        """Set the constraint_types.

        Args:
            constraint_types (list[str | None] | None): The constraint_types of the atoms in the system.
        """
        self._constraint_types = constraint_types
        if (constraint_types is not None) and _is_homogenous(constraint_types):
            self.add_site_property("constraint_types", list(constraint_types))
        elif "constraint_types" in self.site_properties:
            self.remove_site_property("constraint_types")

    @property
    def constraint_vectors(self) -> Sequence | list[np.ndarray | list[np.ndarray] | None] | None:
        """Return the velocities.

        Returns:
            list[np.ndarray | list[np.ndarray] | None] | None: The constraint_vectors of the atoms in the system.
        """
        if "constraint_vectors" not in self.site_properties:
            return self._constraint_vectors
        return self.site_properties["constraint_vectors"]

    @constraint_vectors.setter
    def constraint_vectors(self, constraint_vectors: list[np.ndarray | list[np.ndarray] | None] | None) -> None:
        """Set the constraint_vectors.

        Args:
            constraint_vectors (list[np.ndarray | list[np.ndarray] | None] | None): The constraint_vectors of the
            atoms in the system.
        """
        if (constraint_vectors is not None) and _is_homogenous(constraint_vectors):
            self.add_site_property("constraint_vectors", list(constraint_vectors))
        elif "constraint_vectors" in self.site_properties:
            self.remove_site_property("constraint_vectors")

    @property
    def group_names(self) -> Sequence | list[list[str] | None] | None:
        """Return the group_names.

        Returns:
            list[list[str] | None] | None: The group_names of the atoms in the system.
        """
        if "group_names" not in self.site_properties:
            return self._group_names
        return self.site_properties["group_names"]

    @group_names.setter
    def group_names(self, group_names: list[list[str] | None] | None) -> None:
        """Set the group_names.

        Args:
            group_names (list[list[str] | None] | None): The group_names of the atoms in the system.
        """
        if (group_names is not None) and _is_homogenous(group_names):
            self.add_site_property("group_names", list(group_names))
        elif "group_names" in self.site_properties:
            self.remove_site_property("group_names")

    def __init__(
        self,
        lattice: ArrayLike | Lattice,
        species: Sequence[CompositionLike],
        coords: ArrayLike | Sequence[ArrayLike],
        site_properties: dict | None = None,
        **kwargs,
    ) -> None:
        super().__init__(
            lattice=lattice,
            species=species,
            coords=coords,
            site_properties=site_properties,
            **kwargs,
        )

    @classmethod
    def _from_text_slice(
        cls,
        text_slice: list[str],
        eopt_type: str | None = "ElecMinimize",
        opt_type: str | None = "IonicMinimize",
        emin_flag: str = "---- Electronic minimization -------",
        init_structure: Structure | None = None,
        is_md: bool = False,
    ) -> JOutStructure:
        """
        Return JOutStructure object.

        Create a JAtoms object from a slice of an out file's text corresponding
        to a single step of a native JDFTx optimization.

        Args:
            text_slice (list[str]): A slice of text from a JDFTx out file corresponding to a single
            optimization step / SCF cycle.
            eopt_type (str): The type of electronic minimization step.
            opt_type (str): The type of optimization step.
            emin_flag (str): The flag that indicates the start of a log message for a JDFTx
            optimization step.
            init_structure (Structure | None): The initial structure of the system.
            is_md (bool): Whether the optimization step is a molecular dynamics step.

        Returns:
            JOutStructure: The created JOutStructure object.
        """
        if init_structure is None:
            instance = cls(lattice=np.eye(3), species=[], coords=[], site_properties={})
        else:
            instance = cls(
                lattice=init_structure.lattice.matrix,
                species=init_structure.species,
                coords_are_cartesian=True,
                coords=init_structure.cart_coords,
                site_properties=init_structure.site_properties,
            )
        if opt_type not in ["IonicMinimize", "LatticeMinimize", "IonicDynamics"]:
            opt_type = correct_geom_opt_type(opt_type)
        instance.eopt_type = eopt_type
        instance.opt_type = opt_type
        instance.emin_flag = emin_flag
        instance.is_md = is_md
        line_collections = instance._init_line_collections()
        line_collections = instance._gather_line_collections(line_collections, text_slice)

        # ecomponents needs to be parsed before emin and opt to set etype
        instance._parse_ecomp_lines(line_collections["ecomp"]["lines"])
        if instance.is_md:
            instance._parse_md_opt_lines(line_collections["opt"]["lines"])
        else:
            instance._parse_opt_lines(line_collections["opt"]["lines"])
        instance._parse_emin_lines(line_collections["emin"]["lines"])
        # Lattice must be parsed before posns/forces in case of direct coordinates
        instance._parse_lattice_lines(line_collections["lattice"]["lines"])
        # Posns must be parsed before forces and lowdin analysis so that they can be stored in site_properties
        instance._parse_posns_lines(line_collections["posns"]["lines"])
        instance._parse_forces_lines(line_collections["forces"]["lines"])
        instance._parse_lowdin_lines(line_collections["lowdin"]["lines"])
        # Can be parsed at any point
        instance._parse_strain_lines(line_collections["strain"]["lines"])
        instance._parse_stress_lines(line_collections["stress"]["lines"])
        instance._parse_kinetic_stress_lines(line_collections["kinetic_stress"]["lines"])
        instance._parse_thermostat_line(line_collections["posns"]["lines"])  # Thermostat-velocity dumped with positions
        # In case of single-point calculation
        instance._init_e_sp_backup()
        # Setting attributes from elecmindata (set during _parse_emin_lines)
        instance._elecmindata_postinit()
        # Set relevant properties in self.properties
        instance._fill_properties()
        # Done last in case of any changes to site-properties
        instance._init_structure()
        return instance

    def _init_e_sp_backup(self) -> None:
        """Initialize self.e with coverage for single-point calculations."""
        err_str = None
        if self.e is None:  # This doesn't defer to elecmindata.e due to the existence of a class variable e
            if self.etype is not None:
                if self.ecomponents is not None:
                    if self.etype in self.ecomponents:
                        self.e = self.ecomponents[self.etype]
                    elif self.elecmindata is not None:
                        self.e = self.elecmindata.e
                    else:
                        err_str = "Could not determine total energy due to lack of elecmindata"
                else:
                    err_str = "Could not determine total energy due to lack of ecomponents"
            else:
                err_str = "Could not determine total energy due to lack of etype"
        if err_str is not None:
            warnings.warn(err_str, stacklevel=2)

    def _init_line_collections(self) -> dict:
        """Initialize line collection dict.

        Initialize a dictionary of line collections for each type of line in a
        JDFTx out file.

        Returns:
            dict: A dictionary of line collections for each type of line in a JDFTx
            out file.
        """
        line_collections = {}
        for line_type in self.line_types:
            line_collections[line_type] = {
                "lines": [],
                "collecting": False,
                "collected": False,
            }
        return line_collections

    def _gather_line_collections(self, line_collections: dict, text_slice: list[str]) -> dict:
        """Gather line collections.

        Gather lines of text from a JDFTx out file into a dictionary of line collections.

        Args:
            line_collections (dict): A dictionary of line collections for each type of line in a JDFTx
            out file. Assumed pre-initialized (there exists sdict["lines"]: list[str],
            sdict["collecting"]: bool, sdict["collected"]: bool for every sdict = line_collections[line_type]).
            text_slice (list[str]): A slice of text from a JDFTx out file corresponding to a single
            optimization step / SCF cycle.
        """
        for line in text_slice:
            read_line = False
            for sdict in line_collections.values():
                if sdict["collecting"]:
                    lines, getting, got = self._collect_generic_line(line, sdict["lines"])
                    sdict["lines"] = lines
                    sdict["collecting"] = getting
                    sdict["collected"] = got
                    read_line = True
                    break
            if not read_line:
                for line_type, sdict in line_collections.items():
                    if (not sdict["collected"]) and self._is_generic_start_line(line, line_type):
                        sdict["collecting"] = True
                        sdict["lines"].append(line)
                        break
        return line_collections

    def _is_emin_start_line(self, line_text: str) -> bool:
        """Return True if emin start line.

        Return True if the line_text is the start of a log message for a JDFTx
        optimization step.

        Args:
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the start of a log message for a JDFTx
            optimization step.
        """
        if self.emin_flag is None:
            raise ValueError("emin_flag is not set")
        return self.emin_flag in line_text

    def _is_opt_start_line(self, line_text: str) -> bool:
        """Return True if opt start line.

        Args:
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the start of a log message for a JDFTx optimization step.
        """
        is_line = f"{self.opt_type}:" in line_text
        if self.is_md:
            return is_line and "Step:" in line_text
        return is_line and "Iter:" in line_text

    def _get_etype_from_emin_lines(self, emin_lines: list[str]) -> str | None:
        """Return energy type string.

        Return the type of energy from the electronic minimization data of a
        JDFTx out file.

        Args:
            emin_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            electronic minimization data.

        Returns:
            str: The type of energy from the electronic minimization data of a JDFTx
            out file.
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

    def _set_etype_from_emin_lines(self, emin_lines: list[str]) -> None:
        """Set etype class variable.

        Set the type of energy from the electronic minimization data of a
        JDFTx out file.

        Args:
            emin_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            electronic minimization data.
        """
        self.etype = self._get_etype_from_emin_lines(emin_lines)

    def _parse_emin_lines(self, emin_lines: list[str]) -> None:
        """Parse electronic minimization lines.

        Args:
            emin_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            electronic minimization data.
        """
        if len(emin_lines):
            if self.etype is None:
                self._set_etype_from_emin_lines(emin_lines)
            if self.eopt_type is None:
                raise ValueError("eopt_type is not set")
            if self.etype is None:
                raise ValueError("etype is not set")
            self.elecmindata = JElSteps._from_text_slice(emin_lines, opt_type=self.eopt_type, etype=self.etype)
        else:
            if self.eopt_type is None:
                raise ValueError("eopt_type is not set")
            if self.etype is None:
                raise ValueError("etype is not set")
            self.elecmindata = JElSteps._from_nothing(opt_type=self.eopt_type, etype=self.etype)

    def _parse_lattice_lines(self, lattice_lines: list[str]) -> None:
        """Parse lattice lines.

        Args:
            lattice_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            lattice vectors. Collects the lattice matrix "r" as a 3x3 numpy array first
            in column-major order (vec i = r[:,i]), then transposes it to row-major
            order (vec i = r[i,:]) and converts from Bohr to Angstroms.
        """
        r = None
        if len(lattice_lines) >= 5:
            r = _brkt_list_of_3x3_to_nparray(lattice_lines, i_start=2)
            r = r.T * bohr_to_ang
            self.lattice = Lattice(np.array(r))

    def _parse_strain_lines(self, strain_lines: list[str]) -> None:
        """Parse strain lines.

        Args:
            strain_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            strain tensor. Converts from column-major to row-major order.
        """
        st = None
        if len(strain_lines) == 4:
            st = _brkt_list_of_3x3_to_nparray(strain_lines, i_start=1)
            st = st.T
        self.strain = st

    def _parse_stress_lines(self, stress_lines: list[str]) -> None:
        """Parse stress lines.

        Parse the lines of text corresponding to the stress tensor of a
        JDFTx out file and converts from Ha/Bohr^3 to eV/Ang^3.

        Args:
            stress_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            stress tensor.
        """
        # TODO: Lattice optimizations dump stress in cartesian coordinates in units
        # "[Eh/a0^3]" (Hartree per bohr cubed). Check if this changes for direct
        # coordinates.
        st = None
        if len(stress_lines) == 4:
            st = _brkt_list_of_3x3_to_nparray(stress_lines, i_start=1)
            st = st.T
            st *= Ha_to_eV / (bohr_to_ang**3)
        self.stress = st

    def _parse_kinetic_stress_lines(self, stress_lines: list[str]) -> None:
        """Parse stress (including kinetic component) lines.

        Parse the lines of text corresponding to the kinetic stress tensor of a
        JDFTx out file and converts from Ha/Bohr^3 to eV/Ang^3.

        Args:
            stress_lines (list[str]): A list of lines of text from a JDFTx out file containing the
            stress tensor.
        """
        # TODO: Lattice optimizations dump stress in cartesian coordinates in units
        # "[Eh/a0^3]" (Hartree per bohr cubed). Check if this changes for direct
        # coordinates.
        st = None
        if len(stress_lines) == 4:
            st = _brkt_list_of_3x3_to_nparray(stress_lines, i_start=1)
            st = st.T
            st *= Ha_to_eV / (bohr_to_ang**3)
        self.kinetic_stress = st

    def _parse_thermostat_line(self, posns_lines: list[str]) -> None:
        """Parse thermostat velocity line.

        Parse the lines of text corresponding to the thermostat velocity of a
        JDFTx out file.
        Args:
            posns_lines (list[str]): A list of lines of text from a JDFTx out file.
        """
        if len(posns_lines):
            for line in posns_lines[::-1]:
                if "thermostat-velocity" in line:
                    self.thermostat_velocity = np.array([float(x.strip()) for x in line.split()[1:]])
                elif "ion" in line:
                    break
        else:
            self.thermostat_velocity = None

    def _check_for_structure_consistency(self, names: list[str]) -> bool:
        # If JOutStructure was constructed with a reference init_structure
        if len(self.species):
            if len(names) != len(self.species):
                return False
            _names = list(set(names))
            _self_names = [s.symbol for s in self.species]
            for _name in _names:
                if names.count(_name) != _self_names.count(_name):
                    return False
        return True

    def _parse_posns_lines(self, posns_lines: list[str]) -> None:
        """Parse positions lines.

        Parse the lines of text corresponding to the positions of a
        JDFTx out file.

        Args:
            posns_lines (list[str]): A list of lines of text from a JDFTx out file.
            Collected lines will start with the string "# Ionic positions in ..."
            (specifying either cartesian or direct coordinates), followed by a line
            for each ion (atom) in the format "ion_name x y z sd", where ion_name is
            the name of the element, and sd is a flag indicating whether the ion is
            excluded from optimization (1) or not (0).
        """
        self.copy()
        if len(posns_lines):
            coords_type = posns_lines[0].split("positions in")[1]
            coords_type = coords_type.strip().split()[0].strip()
            _posns: list[NDArray[np.float64]] = []
            names: list[str] = []
            selective_dynamics: list[list[int]] = []
            velocities: list[NDArray[np.float64] | None] = []
            constraint_types: list[str | None] = []
            constraint_vectors: list[list[NDArray[np.float64]] | NDArray[np.float64] | None] = []
            group_names_list: list[list[str] | None] = []
            natoms = len(posns_lines) - 1
            if "thermostat-velocity" in posns_lines[-1]:
                natoms -= 1
            for i in range(natoms):
                line = posns_lines[i + 1].rstrip("\n")
                name, posn, sd, velocity, constraint_type, constraint_vector, group_names = _parse_posn_line(line)
                names.append(name)
                _posns.append(posn)
                selective_dynamics.append(sd)
                velocities.append(velocity)
                constraint_types.append(constraint_type)
                constraint_vectors.append(constraint_vector)
                group_names_list.append(group_names)
            is_good = self._check_for_structure_consistency(names)
            if not is_good and len(self.species):
                # Abort structure updating if we have a pre-existing structure
                return
            self.remove_sites(list(range(len(self.species))))
            posns = np.array(_posns)
            if coords_type.lower() != "cartesian":
                posns = np.dot(posns, self.lattice.matrix)
            else:
                posns *= bohr_to_ang
            for i in range(natoms):
                self.append(species=names[i], coords=posns[i], coords_are_cartesian=True)
            self.selective_dynamics = selective_dynamics
            self.velocities = velocities
            self.constraint_types = constraint_types
            self.constraint_vectors = constraint_vectors
            self.group_names = group_names_list

    def _parse_forces_lines(self, forces_lines: list[str]) -> None:
        """Parse forces lines.

        Args:
            forces_lines (list[str]): A list of lines of text from a JDFTx out file containing the forces.
        """
        if len(forces_lines):
            natoms = len(forces_lines) - 1
            coords_type = forces_lines[0].split("Forces in")[1]
            coords_type = coords_type.strip().split()[0].strip()
            _forces: list[np.ndarray] = []
            for i in range(natoms):
                line = forces_lines[i + 1]
                force = np.array([float(x.strip()) for x in line.split()[2:5]])
                _forces.append(force)
            forces = np.array(_forces)
            if coords_type.lower() != "cartesian":
                forces = np.dot(forces, np.linalg.inv(self.lattice.matrix))
            else:
                forces *= 1 / bohr_to_ang
            forces *= Ha_to_eV
            self.forces = forces
            self.add_site_property("forces", list(forces))

    def _parse_ecomp_lines(self, ecomp_lines: list[str]) -> None:
        """Parse energy component lines.

        Parse the lines of text corresponding to the energy components of a
        JDFTx out file.

        Args:
            ecomp_lines (list[str]): A list of lines of text from a JDFTx out file. All lines will either be
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
        if key is not None and (self.etype is None) and (key in ["F", "G", "Etot"]):
            self.etype = key

    def _parse_lowdin_lines(self, lowdin_lines: list[str]) -> None:
        """Parse Lowdin lines.

        Parse the lines of text corresponding to a Lowdin population analysis
        in a JDFTx out file.

        Args:
            lowdin_lines (list[str]): A list of lines of text from a JDFTx out file.
        """
        charges_dict: dict[str, list[float]] = {}
        moments_dict: dict[str, list[float]] = {}
        for line in lowdin_lines:
            if _is_charges_line(line):
                charges_dict = self._parse_lowdin_line(line, charges_dict)
            elif _is_magnetic_moments_line(line):
                moments_dict = self._parse_lowdin_line(line, moments_dict)
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

    def _parse_lowdin_line(self, lowdin_line: str, lowdin_dict: dict[str, list[float]]) -> dict[str, list[float]]:
        """Parse Lowdin line.

        Parse a line of text from a JDFTx out file corresponding to a
        Lowdin population analysis.

        Args:
            lowdin_line (str): A line of text from a JDFTx out file.
            lowdin_dict (dict[str, list[float]]): A dictionary of Lowdin population analysis data.

        Returns:
            dict[str, list[float]]: A dictionary of Lowdin population analysis data.
        """
        tokens = [v.strip() for v in lowdin_line.strip().split()]
        name = tokens[2]
        vals = [float(x) for x in tokens[3:]]
        lowdin_dict[name] = vals
        return lowdin_dict

    def _is_opt_conv_line(self, line_text: str) -> bool:
        """Return True if line_text is geom opt convergence line.

        Return True if the line_text is the end of a JDFTx optimization step.

        Args:
            line_text (str): A line of text from a JDFTx out file.

        Returns:
            bool: True if the line_text is the end of a JDFTx optimization step.
        """
        return f"{self.opt_type}: Converged" in line_text

    def _parse_opt_lines(self, opt_lines: list[str]) -> None:
        """Parse optimization lines.

        Parse the lines of text corresponding to the optimization step of a
        JDFTx out file.

        Args:
            opt_lines (list[str]): A list of lines of text from a JDFTx out file.
        """
        if len(opt_lines):
            for line in opt_lines:
                if self._is_opt_start_line(line):
                    _nstep = get_colon_val(line, "Iter:")
                    if _nstep is None:
                        raise ValueError("Could not find Iter in line")
                    if np.isnan(_nstep):
                        raise ValueError("Could not convert Iter to int")
                    nstep = int(_nstep)
                    self.nstep = nstep
                    en = get_colon_val(line, f"{self.etype}:")
                    self.e = en * Ha_to_eV
                    grad_k = get_colon_val(line, "|grad|_K: ")
                    self.grad_k = grad_k
                    alpha = get_colon_val(line, "alpha: ")
                    self.alpha = alpha
                    linmin = get_colon_val(line, "linmin: ")
                    self.linmin = linmin
                    t_s = get_colon_val(line, "t[s]: ")
                    self.t_s = t_s
                elif self._is_opt_conv_line(line):
                    self.geom_converged = True
                    self.geom_converged_reason = line.split("(")[1].split(")")[0].strip()

    def _parse_md_opt_lines(self, opt_lines: list[str]) -> None:
        """Parse MD optimization lines.

        Parse the lines of text corresponding to the optimization step of a MD
        JDFTx out file.

        Args:
            opt_lines (list[str]): A list of lines of text from a JDFTx out file.
        """
        if len(opt_lines):
            for line in opt_lines:
                if self._is_opt_start_line(line):
                    _nstep = get_colon_val(line, "Step:")
                    if _nstep is None:
                        raise ValueError("Could not find Step in line")
                    if np.isnan(_nstep):
                        raise ValueError("Could not convert Iter to int")
                    nstep = int(_nstep)
                    self.nstep = nstep
                    self.t_s = get_colon_val(line, "t[s]: ")
                    self.pe = get_colon_val(line, "PE:") * Ha_to_eV
                    self.ke = get_colon_val(line, "KE:") * Ha_to_eV
                    self.t_k = get_colon_val(line, "T[K]:")
                    self.p_bar = get_colon_val(line, "P[Bar]:")
                    self.tmd_fs = get_colon_val(line, "tMD[fs]:")

    def _is_generic_start_line(self, line_text: str, line_type: str) -> bool:
        """Return True if the line_text is start of line_type log message.

        Return True if the line_text is the start of a section of the
        JDFTx out file corresponding to the line_type.

        Args:
            line_text (str): A line of text from a JDFTx out file.
            line_type (str): The type of line to check for.

        Returns:
            bool: True if the line_text is the start of a section of the
            JDFTx out file.
        """
        if line_type == "lowdin":
            return is_lowdin_start_line(line_text)
        if line_type == "opt":
            return self._is_opt_start_line(line_text)
        if line_type == "ecomp":
            return _is_ecomp_start_line(line_text)
        if line_type == "forces":
            return _is_forces_start_line(line_text)
        if line_type == "posns":
            return _is_posns_start_line(line_text)
        if line_type == "stress":
            return _is_stress_start_line(line_text)
        if line_type == "kinetic_stress":
            return _is_kinetic_stress_start_line(line_text)
        if line_type == "strain":
            return _is_strain_start_line(line_text)
        if line_type == "lattice":
            return _is_lattice_start_line(line_text)
        if line_type == "emin":
            return self._is_emin_start_line(line_text)
        raise ValueError(f"Unrecognized line type {line_type}")

    def _collect_generic_line(self, line_text: str, generic_lines: list[str]) -> tuple[list[str], bool, bool]:
        """Collect generic log line.

        Collect a line of text into a list of lines if the line is not empty,
        and otherwise updates the collecting and collected flags.

        Args:
            line_text (str): A line of text from a JDFTx out file.
            generic_lines (list[str]): A list of lines of text of the same type.

        Returns:
            tuple: A tuple containing:
            - generic_lines (list[str]): A list of lines of text of the same type.
            - collecting (bool): True if the line_text is not empty.
            - collected (bool): True if the line_text is empty (end of section).
        """
        collecting = True
        collected = False
        if not len(line_text.strip()):
            collecting = False
            collected = True
        else:
            generic_lines.append(line_text)
        return generic_lines, collecting, collected

    def _fill_properties(self) -> None:
        """Fill properties attribute."""
        self.properties = {
            "eopt_type": self.eopt_type,
            "opt_type": self.opt_type,
            "etype": self.etype,
            "energy": self.e,
            "t_s": self.t_s,
            "geom_converged": self.geom_converged,
            "geom_converged_reason": self.geom_converged_reason,
            "stress": self.stress,
            "kinetic_stress": self.kinetic_stress,
            "strain": self.strain,
            "thermostat_velocity": self.thermostat_velocity,
        }

    def _init_structure(self) -> None:
        """Initialize structure attribute."""
        self.structure = Structure(
            lattice=self.lattice,
            species=self.species,
            coords=self.cart_coords,
            site_properties=self.site_properties,
            coords_are_cartesian=True,
            properties=self.properties,
        )

    def as_dict(self) -> dict:
        """
        Convert the JOutStructure object to a dictionary.

        Returns:
            dict: A dictionary representation of the JOutStructure object.
        """
        dct = {}
        for fld in self.__dict__:
            value = getattr(self, fld)
            if hasattr(value, "as_dict"):
                dct[fld] = value.as_dict()
            else:
                dct[fld] = value
        return dct

    @deprecated(as_dict, deadline=(2025, 10, 4))
    def to_dict(self):
        return self.as_dict()

    # TODO: Add string representation for JOutStructure-specific meta-data
    # This method currently only returns the Structure Summary as inherited from
    # the pymatgen Structure class.
    def __str__(self) -> str:
        """Return string representation.

        Returns:
            str: A string representation of the JOutStructure object.
        """
        return pprint.pformat(self)


# If constraint type is set to "None", it is not printed in out file
constraint_types = ["HyperPlane", "Linear", "Planar"]


def _parse_posn_line(
    posn_line: str,
    is_md: bool = False,
) -> tuple[
    str,
    NDArray[np.float64],
    list[int],
    NDArray[np.float64] | None,
    str | None,
    list[NDArray[np.float64]] | NDArray[np.float64] | None,
    list[str] | None,
]:
    """Parse a single position line.

    Parse a line of text corresponding to the position of an ion (atom) in a
    JDFTx out file. Returns ion name, ion position, selective dynamics, velocity (optional), constraint
    type (optional), constraint vector(s) (optional, may be a list if multiple HyperPlane constraints), and group names
    (if multiple HyperPlane constraints).

    tuple: A tuple containing the following elements:
        - str: Ion name (element symbol)
        - np.ndarray: Ion position coordinates as a 3D vector
        - list[int]: Selective dynamics flags (list of 0/1 values)
        - np.ndarray | None: Velocity as a 3D vector if present, otherwise None
        - str | None: Constraint type (e.g., "HyperPlane") if constrained, otherwise None
        - list[np.ndarray] | np.ndarray | None: Constraint vector(s), returns:
            - list[np.ndarray]: Multiple constraint vectors if multiple HyperPlane constraints
            - np.ndarray: Single constraint vector for a single constraint
            - None: If no constraints are present
        - list[str] | None: Group names for constraints if multiple HyperPlane constraints, otherwise None
    """
    psplit: list[str] = list(posn_line.split())
    name: str = psplit[1]
    posn: NDArray[np.float64] = np.array([float(x) for x in psplit[2:5]])
    velocity: NDArray[np.float64] | None = None
    sd: list[int] = []
    constraint_type: str | None = None
    constraint_vector: NDArray[np.float64] | list[NDArray[np.float64]] | None = None
    group_names: list[str] | None = None
    # Expecting a line to looks like "ion Ir   0.087   0.23   0.285 0"
    # If that last (6th) element is "v", we are dealing with an MD calculation and line might look like
    # "ion Ir   0.087481227766366   0.236770029479291   0.285789315528821 v   0.002  0.10   0.03 1"
    offset: int = 0
    if psplit[5] == "v":
        offset = 4
        velocity = np.array([float(x) for x in psplit[6:9]])
    # Convert the movescale tag to a redundant selective_dynamics tag to match the expected shape
    # (int(bool(v)) used since technically something like "0.1" can be passed to JDFTx to indicate non-freezing)
    sd = [int(bool(posn_line.split()[offset + 5])) for _ in range(3)]
    # Check for constraints (protected by try/except since its genuinely more likely that an accidental edit
    # triggers this since this is a rarely used feature)
    try:
        if len(psplit) > offset + 6:
            constraint_type = psplit[offset + 6]
            if constraint_type in constraint_types:
                # Check for multiple HyperPlane constraints
                if constraint_type == "HyperPlane":
                    if psplit.count("HyperPlane") > 1:
                        idcs = [i for i, v in enumerate(psplit) if v == "HyperPlane"]
                        constraint_vector = [np.array([float(x) for x in psplit[i + 1 : i + 4]]) for i in idcs]
                        group_names = [psplit[i + 4] for i in idcs]
                    else:
                        constraint_vector = np.array([float(x) for x in psplit[offset + 7 : offset + 10]])
                else:
                    constraint_vector = np.array([float(x) for x in psplit[offset + 7 : offset + 10]])
            else:
                # Ignore unrecognized constraint types
                constraint_type = None
    except (IndexError, ValueError, TypeError):
        pass
    return name, posn, sd, velocity, constraint_type, constraint_vector, group_names


def _is_stress_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of stress log message.

    Return True if the line_text is the start of stress log message.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of stress log message.
    """
    return "# Stress tensor in" in line_text


def _is_kinetic_stress_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of kinetic stress log message.

    Return True if the line_text is the start of kinetic stress log message.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of kinetic stress log message.
    """
    return "# Stress tensor including kinetic" in line_text


def _is_strain_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of strain log message.

    Return True if the line_text is the start of strain log message.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of strain log message.
    """
    return "# Strain tensor in" in line_text


def _is_posns_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of posns log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Args:
        line_text (str): A line of text from a JDFTx out file containing the positions of atoms.

    Returns:
        bool: True if the line_text is the start of a log message for a JDFTx optimization step.
    """
    return "# Ionic positions" in line_text


def _is_ecomp_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of ecomp log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of a log message for a JDFTx
        optimization step.
    """
    return "# Energy components" in line_text


def _is_charges_line(line_text: str) -> bool:
    """Return True if the line_text is start of charges log message.

    Return True if the line_text is a line of text from a JDFTx out file
    corresponding to a Lowdin population analysis.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is a line of text from a JDFTx out file
        corresponding to a Lowdin population.
    """
    return "oxidation-state" in line_text


def _is_magnetic_moments_line(line_text: str) -> bool:
    """Return True if the line_text is start of moments log message.

    Return True if the line_text is a line of text from a JDFTx out file
    corresponding to a Lowdin population analysis.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is a line of text from a JDFTx out file
        corresponding to a Lowdin population.
    """
    return "magnetic-moments" in line_text


def _is_forces_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of forces log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of a log message for a JDFTx
        optimization step.
    """
    return "# Forces in" in line_text


def _is_lattice_start_line(line_text: str) -> bool:
    """Return True if the line_text is the start of lattice log message.

    Return True if the line_text is the start of a log message for a JDFTx
    optimization step.

    Args:
        line_text (str): A line of text from a JDFTx out file.

    Returns:
        bool: True if the line_text is the start of a log message for a JDFTx
        optimization step.
    """
    return "# Lattice vectors:" in line_text
