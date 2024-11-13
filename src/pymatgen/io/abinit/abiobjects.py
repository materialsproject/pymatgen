"""Low-level classes and functions to work with Abinit input files."""

from __future__ import annotations

import abc
from collections.abc import Iterable
from enum import Enum, unique
from pprint import pformat
from typing import TYPE_CHECKING, NamedTuple, cast

import numpy as np
from monty.collections import AttrDict
from monty.design_patterns import singleton
from monty.json import MontyDecoder, MontyEncoder, MSONable

from pymatgen.core import ArrayWithUnit, Lattice, Species, Structure, units

if TYPE_CHECKING:
    from typing import Any, ClassVar

    from typing_extensions import Self


def lattice_from_abivars(cls=None, *args, **kwargs):
    """Get a `Lattice` object from a dictionary with the Abinit variables `acell`
    and either `rprim` in Bohr or `angdeg`. If acell is not given, the Abinit default
    of [1, 1, 1] Bohr is used.

    Args:
        cls: Lattice class to be instantiated. Defaults to pymatgen.core.Lattice.

    Example:
        lattice_from_abivars(acell=3*[10], rprim=np.eye(3))
    """
    cls = cls or Lattice
    kwargs.update(dict(*args))

    r_prim = kwargs.get("rprim")
    ang_deg = kwargs.get("angdeg")
    a_cell = kwargs["acell"]

    if r_prim is not None:
        if ang_deg is not None:
            raise ValueError("angdeg and rprimd are mutually exclusive")
        r_prim = np.reshape(r_prim, (3, 3))
        rprimd = [float(a_cell[i]) * r_prim[i] for i in range(3)]
        # Call pymatgen constructors (note that pymatgen uses Angstrom instead of Bohr).
        return cls(ArrayWithUnit(rprimd, "bohr").to("ang"))

    if ang_deg is not None:
        ang_deg = np.reshape(ang_deg, 3)

        if np.any(ang_deg <= 0.0):
            raise ValueError(f"Angles must be > 0 but got {ang_deg}")
        if ang_deg.sum() >= 360.0:
            raise ValueError(f"The sum of angdeg must be lower than 360, {ang_deg=}")

        # This code follows the implementation in ingeo.F90
        # See also https://docs.abinit.org/variables/basic/#angdeg
        tol12 = 1e-12
        pi, sin, cos, sqrt = np.pi, np.sin, np.cos, np.sqrt
        r_prim = np.zeros((3, 3))
        if (
            abs(ang_deg[0] - ang_deg[1]) < tol12
            and abs(ang_deg[1] - ang_deg[2]) < tol12
            and abs(ang_deg[0] - 90.0) + abs(ang_deg[1] - 90.0) + abs(ang_deg[2] - 90) > tol12
        ):
            # Treat the case of equal angles (except all right angles):
            # generates trigonal symmetry w.r.t. third axis
            cos_ang = cos(pi * ang_deg[0] / 180.0)
            a2 = 2.0 / 3.0 * (1.0 - cos_ang)
            aa = sqrt(a2)
            cc = sqrt(1.0 - a2)
            r_prim[0, 0] = aa
            r_prim[0, 1] = 0.0
            r_prim[0, 2] = cc
            r_prim[1, 0] = -0.5 * aa
            r_prim[1, 1] = sqrt(3.0) * 0.5 * aa
            r_prim[1, 2] = cc
            r_prim[2, 0] = -0.5 * aa
            r_prim[2, 1] = -sqrt(3.0) * 0.5 * aa
            r_prim[2, 2] = cc
        else:
            # Treat all the other cases
            r_prim[0, 0] = 1.0
            r_prim[1, 0] = cos(pi * ang_deg[2] / 180.0)
            r_prim[1, 1] = sin(pi * ang_deg[2] / 180.0)
            r_prim[2, 0] = cos(pi * ang_deg[1] / 180.0)
            r_prim[2, 1] = (cos(pi * ang_deg[0] / 180.0) - r_prim[1, 0] * r_prim[2, 0]) / r_prim[1, 1]
            r_prim[2, 2] = sqrt(1.0 - r_prim[2, 0] ** 2 - r_prim[2, 1] ** 2)

        # Call pymatgen constructors (note that pymatgen uses Angstrom instead of Bohr).
        rprimd = [float(a_cell[i]) * r_prim[i] for i in range(3)]
        return cls(ArrayWithUnit(rprimd, "bohr").to("ang"))

    raise ValueError(f"Don't know how to construct a Lattice from dict:\n{pformat(kwargs)}")


def structure_from_abivars(cls=None, *args, **kwargs) -> Structure:
    """
    Build a Structure object from a dictionary with ABINIT variables.

    Args:
        cls: Structure class to be instantiated. Defaults to Structure.

    Example:
        al_structure = structure_from_abivars(
            acell=3*[7.5],
            rprim=[0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0],
            typat=1,
            xred=[0.0, 0.0, 0.0],
            ntypat=1,
            znucl=13,
        )

    xred can be replaced with xcart or xangst.
    """
    kwargs.update(dict(*args))
    cls = cls or Structure

    # lattice = Lattice.from_dict(d, fmt="abivars")
    lattice = lattice_from_abivars(**kwargs)
    coords, coords_are_cartesian = kwargs.get("xred"), False

    if coords is None:
        coords = kwargs.get("xcart")
        if coords is not None:
            if "xangst" in kwargs:
                raise ValueError("xangst and xcart are mutually exclusive")
            coords = ArrayWithUnit(coords, "bohr").to("ang")
        else:
            coords = kwargs.get("xangst")
        coords_are_cartesian = True

    if coords is None:
        raise ValueError(f"Cannot extract coordinates from:\n {kwargs}")

    coords = np.reshape(coords, (-1, 3))

    znucl_type, typat = kwargs["znucl"], kwargs["typat"]

    if not isinstance(znucl_type, Iterable):
        znucl_type = [znucl_type]

    if not isinstance(typat, Iterable):
        typat = [typat]

    if len(typat) != len(coords):
        raise ValueError(f"{len(typat)=} must equal {len(coords)=}")

    # Note conversion to int and Fortran --> C indexing
    typat = np.array(typat, dtype=np.int64)
    species = [znucl_type[typ - 1] for typ in typat]

    return cls(
        lattice,
        species,
        coords,
        validate_proximity=False,
        to_unit_cell=False,
        coords_are_cartesian=coords_are_cartesian,
    )


def species_by_znucl(structure: Structure) -> list[Species]:
    """Get list of unique specie found in structure **ordered according to sites**.

    Example:
        Site0: 0.5 0 0 O
        Site1: 0   0 0 Si

    produces [Specie_O, Specie_Si] and not set([Specie_O, Specie_Si]) as in `types_of_specie`
    """
    # Please, do not change this algorithm and DO NOT USE set.
    # This logic produces a deterministic order of the species based on their first occurrence in sites.
    # This is something one can easily implement in Fortran to facilitate interoperability between pymatgen
    # and Abinit. Most importantly, we can reuse all the DFPT results produced so far in which the
    # old version of structure.types_of_specie (equivalent to this one) was used!
    types = []
    for site in structure:
        for sp, v in site.species.items():
            if sp not in types and v != 0:
                types.append(sp)
    return types


def structure_to_abivars(
    structure: Structure,
    enforce_znucl: list | None = None,
    enforce_typat: list | None = None,
    **kwargs,
):
    """
    Receives a structure and returns a dictionary with ABINIT variables.

    Args:
        enforce_znucl (list): ntypat entries with the value of Z for each type of atom.
            Used to change the default ordering. Defaults to None.
        enforce_typat (list): natom entries with the type index.
            Fortran conventions: start to count from 1. Used to change the default ordering.
    """
    if not structure.is_ordered:
        raise ValueError(
            "Received disordered structure with partial occupancies that cannot be converted into an "
            "Abinit input. Please use OrderDisorderedStructureTransformation or EnumerateStructureTransformation "
            "to build an appropriate supercell from partial occupancies or, alternatively, use the Rigid Band Model "
            "or the Virtual Crystal Approximation."
        )

    n_atoms = len(structure)
    enforce_order = False

    if enforce_znucl is not None or enforce_typat is not None:
        enforce_order = True
        # consistency check
        if enforce_znucl is None or enforce_typat is None:
            raise ValueError("Both enforce_znucl and enforce_typat are required!")

        if len(enforce_typat) != len(structure):
            raise ValueError(
                f"enforce_typat contains {len(enforce_typat)} entries while it should be {len(structure)=}"
            )

        if len(enforce_znucl) != structure.n_elems:
            raise ValueError(
                f"enforce_znucl contains {len(enforce_znucl)} entries while it should be {structure.n_elems=}"
            )

    if enforce_order:
        znucl_type = enforce_znucl
        typat = enforce_typat or []  # or [] added for mypy
    else:
        types_of_specie = species_by_znucl(structure)

        znucl_type = [specie.number for specie in types_of_specie]
        typat = np.zeros(n_atoms, int)
        for atm_idx, site in enumerate(structure):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

    r_prim = ArrayWithUnit(structure.lattice.matrix, "ang").to("bohr")
    ang_deg = structure.lattice.angles
    x_red = np.reshape([site.frac_coords for site in structure], (-1, 3))

    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    r_prim = np.where(np.abs(r_prim) > 1e-8, r_prim, 0.0)
    x_red = np.where(np.abs(x_red) > 1e-8, x_red, 0.0)

    # Info on atoms.
    dct = {
        "natom": n_atoms,
        "ntypat": structure.n_elems,
        "typat": typat,
        "znucl": znucl_type,
        "xred": x_red,
    }

    # Add info on the lattice.
    # Should we use (rprim, acell) or (ang_deg, acell) to specify the lattice?
    geo_mode = kwargs.pop("geomode", "rprim")
    if geo_mode == "automatic":
        geo_mode = "rprim"
        if structure.lattice.is_hexagonal():  # or structure.lattice.is_rhombohedral
            geo_mode = "angdeg"
            ang_deg = structure.lattice.angles
            # Here one could polish a bit the numerical values if they are not exact.
            # Note that in pmg the angles are 12, 20, 01 while in Abinit 12, 02, 01
            # One should make sure that the orientation is preserved (see Curtarolo's settings)

    if geo_mode == "rprim":
        dct.update(acell=3 * [1.0], rprim=r_prim)

    elif geo_mode == "angdeg":
        dct.update(
            acell=ArrayWithUnit(structure.lattice.abc, "ang").to("bohr"),
            angdeg=ang_deg,
        )
    else:
        raise ValueError(f"Wrong value for {geo_mode=}")

    return dct


def contract(string):
    """
    Examples:
        assert contract("1 1 1 2 2 3") == "3*1 2*2 1*3"
        assert contract("1 1 3 2 3") == "2*1 1*3 1*2 1*3".
    """
    if not string:
        return string

    tokens = string.split()
    old = tokens[0]
    count = [[1, old]]

    for tok in tokens[1:]:
        if tok == old:
            count[-1][0] += 1
        else:
            old = tok
            count.append([1, tok])

    return " ".join(f"{c}*{t}" for c, t in count)


class AbivarAble(abc.ABC):
    """An AbivarAble object provides a method to_abivars that returns a dictionary with the abinit variables."""

    @abc.abstractmethod
    def to_abivars(self):
        """Get a dictionary with the abinit variables."""

    # @abc.abstractmethod
    # def from_abivars(cls, vars):
    #    """Build the object from a dictionary with Abinit variables."""

    def __str__(self):
        return pformat(self.to_abivars(), indent=1, width=80, depth=None)

    def __contains__(self, key):
        return key in self.to_abivars()


@singleton
class MandatoryVariable:
    """
    Singleton used to tag mandatory variables, just because I can use
    the cool syntax: variable is MANDATORY!
    """

    def as_dict(self):
        return {}


@singleton
class DefaultVariable:
    """Singleton used to tag variables that will have the default value."""


MANDATORY = MandatoryVariable()
DEFAULT = DefaultVariable()


class SpinModeTuple(NamedTuple):
    mode: str
    nsppol: int
    nspinor: int
    nspden: int


class SpinMode(SpinModeTuple, AbivarAble, MSONable):
    """
    Different configurations of the electron density as implemented in abinit:
    One can use as_spinmode to construct the object via SpinMode.as_spinmode
    (string) where string can assume the values:

        - polarized
        - unpolarized
        - afm (anti-ferromagnetic)
        - spinor (non-collinear magnetism)
        - spinor_nomag (non-collinear, no magnetism)
    """

    __slots__ = ()

    @classmethod
    def as_spinmode(cls, obj):
        """Convert obj into a `SpinMode` instance."""
        if isinstance(obj, cls):
            return obj

        # Assume a string with mode
        try:
            return _mode_to_spin_vars[obj]
        except KeyError:
            raise KeyError(f"Wrong value for spin_mode: {obj}")

    def to_abivars(self):
        """Dictionary with Abinit input variables."""
        return {"nsppol": self.nsppol, "nspinor": self.nspinor, "nspden": self.nspden}

    def as_dict(self):
        """JSON-friendly dict representation of SpinMode."""
        out = {k: getattr(self, k) for k in self._fields}
        out |= {"@module": type(self).__module__, "@class": type(self).__name__}
        return out

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        return cls(**{key: dct[key] for key in dct if key in cls._fields})


_mode_to_spin_vars = {
    "unpolarized": SpinMode("unpolarized", 1, 1, 1),
    "polarized": SpinMode("polarized", 2, 1, 2),
    "afm": SpinMode("afm", 1, 1, 2),
    "spinor": SpinMode("spinor", 1, 2, 4),
    "spinor_nomag": SpinMode("spinor_nomag", 1, 2, 1),
}


class Smearing(AbivarAble, MSONable):
    """
    Variables defining the smearing technique. The preferred way to instantiate
    a `Smearing` object is via the class method Smearing.as_smearing(string).
    """

    # Map string_mode to occopt
    _mode2occopt: ClassVar[dict[str, int]] = {
        "nosmearing": 1,
        "fermi_dirac": 3,
        "marzari4": 4,
        "marzari5": 5,
        "methfessel": 6,
        "gaussian": 7,
    }

    def __init__(self, occopt, tsmear):
        """
        Args:
            occopt: Integer specifying the smearing technique.
            tsmear: Smearing parameter in Hartree units.
        """
        self.occopt = occopt
        self.tsmear = tsmear

    def __str__(self):
        string = f"occopt {self.occopt} # {self.mode} Smearing\n"
        if self.tsmear:
            string += f"tsmear {self.tsmear}"
        return string

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("occopt", "tsmear")
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        other = cast(Smearing, other)

        return self.occopt == other.occopt and np.allclose(self.tsmear, other.tsmear)

    def __bool__(self):
        return self.mode != "nosmearing"

    @classmethod
    def as_smearing(cls, obj):
        """
        Constructs an instance of `Smearing` from obj. Accepts obj in the form:

            * Smearing instance
            * "name:tsmear"  e.g. "gaussian:0.004"  (Hartree units)
            * "name:tsmear units" e.g. "gaussian:0.1 eV"
            * None --> no smearing
        """
        if obj is None:
            return Smearing.nosmearing()

        if isinstance(obj, cls):
            return obj

        # obj is a string
        if obj == "nosmearing":
            return cls.nosmearing()

        obj, tsmear = obj.split(":")
        obj.strip()

        occopt = cls._mode2occopt[obj]
        try:
            tsmear = float(tsmear)
        except ValueError:
            tsmear, unit = tsmear.split()
            tsmear = units.Energy(float(tsmear), unit).to("Ha")

        return cls(occopt, tsmear)

    @property
    def mode(self):
        """String with smearing technique."""
        for mode_str, occopt in self._mode2occopt.items():
            if occopt == self.occopt:
                return mode_str
        raise AttributeError(f"Unknown occopt {self.occopt}")

    @staticmethod
    def nosmearing():
        """For calculations without smearing."""
        return Smearing(1, 0.0)

    def to_abivars(self):
        """Return dictionary with Abinit variables."""
        if self.mode == "nosmearing":
            return {"occopt": 1, "tsmear": 0.0}
        return {"occopt": self.occopt, "tsmear": self.tsmear}

    def as_dict(self):
        """JSON-friendly dict representation of Smearing."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "occopt": self.occopt,
            "tsmear": self.tsmear,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        return cls(dct["occopt"], dct["tsmear"])


class ElectronsAlgorithm(dict, AbivarAble, MSONable):
    """Variables controlling the SCF/NSCF algorithm."""

    # None indicates that we use abinit defaults
    _DEFAULT: ClassVar[dict[str, int | None]] = {
        "iprcell": None,
        "iscf": None,
        "diemac": None,
        "diemix": None,
        "diemixmag": None,
        "dielam": None,
        "diegap": None,
        "dielng": None,
        "diecut": None,
        "nstep": 50,
    }

    def __init__(self, *args, **kwargs):
        """
        Args:
            iprcell: 1 if the cell is fixed, 2 if the cell is relaxed.
            iscf: SCF algorithm.
            diemac: Macroscopic dielectric constant.
            diemix: Mixing parameter for the electric field.
            diemixmag: Mixing parameter for the magnetic field.
            dielam: Damping factor for the electric field.
            diegap: Energy gap for the dielectric function.
            dielng: Length of the electric field.
            diecut: Cutoff for the dielectric function.
            nstep: Maximum number of SCF iterations.
        """
        super().__init__(*args, **kwargs)

        for key in self:
            if key not in self._DEFAULT:
                raise ValueError(f"{type(self).__name__}: No default value has been provided for {key=}")

    def to_abivars(self):
        """Dictionary with Abinit input variables."""
        return self.copy()

    def as_dict(self):
        """Get JSON-able dict representation."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            **self.copy(),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        dct = dct.copy()
        dct.pop("@module", None)
        dct.pop("@class", None)
        return cls(**dct)


class Electrons(AbivarAble, MSONable):
    """The electronic degrees of freedom."""

    def __init__(
        self,
        spin_mode="polarized",
        smearing="fermi_dirac:0.1 eV",
        algorithm=None,
        nband=None,
        fband=None,
        charge=0.0,
        comment=None,
    ):  # occupancies=None,
        """
        Constructor for Electrons object.

        Args:
            comment: String comment for Electrons
            charge: Total charge of the system. Default is 0.
        """
        super().__init__()

        self.comment = comment
        self.smearing = Smearing.as_smearing(smearing)
        self.spin_mode = SpinMode.as_spinmode(spin_mode)
        self.nband = nband
        self.fband = fband
        self.charge = charge
        self.algorithm = algorithm

    @property
    def nsppol(self) -> int:
        """Number of independent spin polarizations."""
        return self.spin_mode.nsppol

    @property
    def nspinor(self):
        """Number of independent spinor components."""
        return self.spin_mode.nspinor

    @property
    def nspden(self):
        """Number of independent density components."""
        return self.spin_mode.nspden

    def as_dict(self):
        """JSON friendly dict representation."""
        dct = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["spin_mode"] = self.spin_mode.as_dict()
        dct["smearing"] = self.smearing.as_dict()
        dct["algorithm"] = self.algorithm.as_dict() if self.algorithm else None
        dct["nband"] = self.nband
        dct["fband"] = self.fband
        dct["charge"] = self.charge
        dct["comment"] = self.comment
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        dct = dct.copy()
        dct.pop("@module", None)
        dct.pop("@class", None)
        dct["spin_mode"] = MontyDecoder().process_decoded(dct["spin_mode"])
        dct["smearing"] = MontyDecoder().process_decoded(dct["smearing"])
        dct["algorithm"] = MontyDecoder().process_decoded(dct["algorithm"]) if dct["algorithm"] else None
        return cls(**dct)

    def to_abivars(self):
        """Return dictionary with Abinit variables."""
        abivars = self.spin_mode.to_abivars()

        abivars |= {"nband": self.nband, "fband": self.fband, "charge": self.charge}

        if self.smearing:
            abivars.update(self.smearing.to_abivars())

        if self.algorithm:
            abivars.update(self.algorithm)

        # abivars["#comment"] = self.comment
        return abivars


@unique
class KSamplingModes(Enum):
    """Enum if the different samplings of the BZ."""

    monkhorst = 1
    path = 2
    automatic = 3


class KSampling(AbivarAble, MSONable):
    """Input variables defining the K-point sampling."""

    def __init__(
        self,
        mode=KSamplingModes.monkhorst,
        num_kpts=0,
        kpts=((1, 1, 1),),
        kpt_shifts=(0.5, 0.5, 0.5),
        kpts_weights=None,
        use_symmetries=True,
        use_time_reversal=True,
        chksymbreak=None,
        comment=None,
    ):
        """
        Highly flexible constructor for KSampling objects. The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility. For most usage cases, the object be constructed
        far more easily using the convenience static constructors:

            #. gamma_only
            #. gamma_centered
            #. monkhorst
            #. monkhorst_automatic
            #. path

        and it is recommended that you use those.

        Args:
            mode: Mode for generating k-poits. Use one of the KSamplingModes enum types.
            num_kpts: Number of kpoints if mode is "automatic"
                Number of division for the sampling of the smallest segment if mode is "path".
                Not used for the other modes
            kpts: Number of divisions. Even when only a single specification is
                required, e.g. in the automatic scheme, the kpts should still
                be specified as a 2D array. e.g. [[20]] or [[2,2,2]].
            kpt_shifts: Shifts for Kpoints.
            use_symmetries: False if spatial symmetries should not be used
                to reduce the number of independent k-points.
            use_time_reversal: False if time-reversal symmetry should not be used
                to reduce the number of independent k-points.
            kpts_weights: Optional weights for kpoints. For explicit kpoints.
            chksymbreak: Abinit input variable: check whether the BZ sampling preserves the symmetry of the crystal.
            comment: String comment for Kpoints

        Note:
            The default behavior of the constructor is monkhorst.
        """
        if isinstance(mode, str):
            mode = KSamplingModes[mode]

        super().__init__()

        self.mode = mode
        self.comment = comment

        self.num_kpts = num_kpts
        self.kpts = kpts
        self.kpt_shifts = kpt_shifts
        self.kpts_weights = kpts_weights
        self.use_symmetries = use_symmetries
        self.use_time_reversal = use_time_reversal
        self.chksymbreak = chksymbreak

        abivars = {}

        if mode == KSamplingModes.monkhorst:
            if num_kpts != 0:
                raise ValueError(f"expect num_kpts to be zero, got {num_kpts}")
            ngkpt = np.reshape(kpts, 3)
            shiftk = np.reshape(kpt_shifts, (-1, 3))

            if use_symmetries and use_time_reversal:
                kptopt = 1
            elif not use_symmetries and use_time_reversal:
                kptopt = 2
            elif not use_symmetries and not use_time_reversal:
                kptopt = 3
            else:  # use_symmetries and not use_time_reversal
                kptopt = 4

            abivars |= {
                "ngkpt": ngkpt,
                "shiftk": shiftk,
                "nshiftk": len(shiftk),
                "kptopt": kptopt,
                "chksymbreak": chksymbreak,
            }

        elif mode == KSamplingModes.path:
            if num_kpts <= 0:
                raise ValueError("For Path mode, num_kpts must be specified and >0")

            kptbounds = np.reshape(kpts, (-1, 3))

            abivars |= {
                "ndivsm": num_kpts,
                "kptbounds": kptbounds,
                "kptopt": -len(kptbounds) + 1,
            }

        elif mode == KSamplingModes.automatic:
            kpts = np.reshape(kpts, (-1, 3))
            if len(kpts) != num_kpts:
                raise ValueError("For Automatic mode, num_kpts must be specified.")

            abivars |= {
                "kptopt": 0,
                "kpt": kpts,
                "nkpt": num_kpts,
                "kptnrm": np.ones(num_kpts),
                "wtk": kpts_weights,  # for iscf/=-2, wtk.
                "chksymbreak": chksymbreak,
            }

        else:
            raise ValueError(f"Unknown {mode=}")

        self.abivars = abivars
        # self.abivars["#comment"] = comment

    @property
    def is_homogeneous(self) -> bool:
        """Homogeneous sampling."""
        return self.mode != "path"

    @classmethod
    def gamma_only(cls):
        """Gamma-only sampling."""
        return cls(kpt_shifts=(0.0, 0.0, 0.0), comment="Gamma-only sampling")

    @classmethod
    def gamma_centered(cls, kpts=(1, 1, 1), use_symmetries=True, use_time_reversal=True):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            use_symmetries: False if spatial symmetries should not be used
                to reduce the number of independent k-points.
            use_time_reversal: False if time-reversal symmetry should not be used
                to reduce the number of independent k-points.

        Returns:
            KSampling object.
        """
        return cls(
            kpts=[kpts],
            kpt_shifts=(0.0, 0.0, 0.0),
            use_symmetries=use_symmetries,
            use_time_reversal=use_time_reversal,
            comment="gamma-centered mode",
        )

    @classmethod
    def monkhorst(
        cls,
        ngkpt,
        shiftk=(0.5, 0.5, 0.5),
        chksymbreak=None,
        use_symmetries=True,
        use_time_reversal=True,
        comment=None,
    ):
        """
        Convenient static constructor for a Monkhorst-Pack mesh.

        Args:
            ngkpt: Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            shiftk: Shift to be applied to the kpoints.
            use_symmetries: Use spatial symmetries to reduce the number of k-points.
            use_time_reversal: Use time-reversal symmetry to reduce the number of k-points.

        Returns:
            KSampling object.
        """
        return cls(
            kpts=[ngkpt],
            kpt_shifts=shiftk,
            use_symmetries=use_symmetries,
            use_time_reversal=use_time_reversal,
            chksymbreak=chksymbreak,
            comment=comment or "Monkhorst-Pack scheme with user-specified shiftk",
        )

    @classmethod
    def monkhorst_automatic(
        cls,
        structure,
        ngkpt,
        use_symmetries=True,
        use_time_reversal=True,
        chksymbreak=None,
        comment=None,
    ):
        """
        Convenient static constructor for an automatic Monkhorst-Pack mesh.

        Args:
            structure: Structure object.
            ngkpt: Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            use_symmetries: Use spatial symmetries to reduce the number of k-points.
            use_time_reversal: Use time-reversal symmetry to reduce the number of k-points.

        Returns:
            KSampling object.
        """
        # TODO
        shiftk = 3 * (0.5,)

        # if lattice.ishexagonal:
        # elif lattice.isbcc
        # elif lattice.isfcc

        return cls.monkhorst(
            ngkpt,
            shiftk=shiftk,
            use_symmetries=use_symmetries,
            use_time_reversal=use_time_reversal,
            chksymbreak=chksymbreak,
            comment=comment or "Automatic Monkhorst-Pack scheme",
        )

    @classmethod
    def _path(cls, ndivsm, structure=None, kpath_bounds=None, comment=None):
        """Static constructor for path in k-space.

        Args:
            structure: Structure object.
            kpath_bounds: List with the reduced coordinates of the k-points defining the path.
            ndivsm: Number of division for the smallest segment.
            comment: Comment string.

        Returns:
            KSampling object.
        """
        if kpath_bounds is None:
            # Compute the boundaries from the input structure.
            from pymatgen.symmetry.bandstructure import HighSymmKpath

            sp = HighSymmKpath(structure)

            # Flat the array since "path" is a list of lists!
            kpath_labels = []
            for labels in sp.kpath["path"]:
                kpath_labels.extend(labels)

            kpath_bounds = []
            for label in kpath_labels:
                red_coord = sp.kpath["kpoints"][label]
                kpath_bounds.append(red_coord)

        return cls(
            mode=KSamplingModes.path,
            num_kpts=ndivsm,
            kpts=kpath_bounds,
            comment=comment or "K-Path scheme",
        )

    @classmethod
    def path_from_structure(cls, ndivsm, structure) -> Self:
        """See _path for the meaning of the variables."""
        return cls._path(
            ndivsm,
            structure=structure,
            comment="K-path generated automatically from structure",
        )

    @classmethod
    def explicit_path(cls, ndivsm, kpath_bounds):
        """See _path for the meaning of the variables."""
        return cls._path(ndivsm, kpath_bounds=kpath_bounds, comment="Explicit K-path")

    @classmethod
    def automatic_density(
        cls,
        structure,
        kppa,
        chksymbreak=None,
        use_symmetries=True,
        use_time_reversal=True,
        shifts=(0.5, 0.5, 0.5),
    ):
        """Get an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes for hexagonal cells and Monkhorst-Pack grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure: Input structure
            kppa: Grid density
        """
        lattice = structure.lattice
        lengths = lattice.abc
        shifts = np.reshape(shifts, (-1, 3))
        ngrid = kppa / len(structure) / len(shifts)

        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3.0)

        num_div = [int(round(1.0 / lengths[i] * mult)) for i in range(3)]
        # ensure that num_div[i] > 0
        num_div = [i if i > 0 else 1 for i in num_div]

        comment = f"pymatge.io.abinit generated KPOINTS with grid density = {kppa} / atom"

        return cls(
            mode="monkhorst",
            num_kpts=0,
            kpts=[num_div],
            kpt_shifts=shifts,
            use_symmetries=use_symmetries,
            use_time_reversal=use_time_reversal,
            chksymbreak=chksymbreak,
            comment=comment,
        )

    def to_abivars(self):
        """Dictionary with Abinit variables."""
        return self.abivars

    def as_dict(self):
        """Get JSON-able dict representation."""
        enc = MontyEncoder()
        return {
            "mode": self.mode.name,
            "comment": self.comment,
            "num_kpts": self.num_kpts,
            "kpts": enc.default(np.array(self.kpts)),
            "kpt_shifts": self.kpt_shifts,
            "kpts_weights": self.kpts_weights,
            "use_symmetries": self.use_symmetries,
            "use_time_reversal": self.use_time_reversal,
            "chksymbreak": self.chksymbreak,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        dct = dct.copy()
        dct.pop("@module", None)
        dct.pop("@class", None)
        dct["kpts"] = MontyDecoder().process_decoded(dct["kpts"])
        return cls(**dct)


class Constraints(AbivarAble):
    """Define the constraints for structural relaxation."""

    def to_abivars(self):
        """Dictionary with Abinit variables."""
        raise NotImplementedError("")


class RelaxationMethod(AbivarAble, MSONable):
    """
    This object stores the variables for the (constrained) structural optimization
    ionmov and optcell specify the type of relaxation.
    The other variables are optional and their use depend on ionmov and optcell.
    A None value indicates that we use abinit default. Default values can
    be modified by passing them to the constructor.
    The set of variables are constructed in to_abivars depending on ionmov and optcell.
    """

    _default_vars: ClassVar[dict[str, Any]] = {
        "ionmov": MANDATORY,
        "optcell": MANDATORY,
        "ntime": 80,
        "dilatmx": 1.05,
        "ecutsm": 0.5,
        "strfact": None,
        "tolmxf": None,
        "strtarget": None,
        "atoms_constraints": {},  # Constraints are stored in a dictionary. {} means if no constraint is enforced.
    }

    IONMOV_DEFAULT = 3
    OPTCELL_DEFAULT = 2

    def __init__(self, *args, **kwargs):
        """
        Args:
            ionmov: The type of relaxation for the ions.
            optcell: The type of relaxation for the unit cell.
            ntime: Maximum number of iterations.
            dilatmx: Maximum allowed cell volume change.
            ecutsm: Energy convergence criterion.
            strfact: Stress convergence criterion.
            tolmxf: Force convergence criterion.
            strtarget: Target stress.
            atoms_constraints: Constraints for the atoms.
        """
        # Initialize abivars with the default values.
        self.abivars = {**self._default_vars}

        # Overwrite the keys with the args and kwargs passed to constructor.
        self.abivars.update(*args, **kwargs)

        self.abivars = AttrDict(self.abivars)

        for key in self.abivars:
            if key not in self._default_vars:
                raise ValueError(f"{type(self).__name__}: No default value has been provided for {key=}")

        for key in self.abivars:
            if key is MANDATORY:
                raise ValueError(f"{type(self).__name__}: No default value has been provided for the mandatory {key=}")

    @classmethod
    def atoms_only(cls, atoms_constraints=None):
        """Relax atomic positions, keep unit cell fixed."""
        if atoms_constraints is None:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=0)
        return cls(ionmov=cls.IONMOV_DEFAULT, optcell=0, atoms_constraints=atoms_constraints)

    @classmethod
    def atoms_and_cell(cls, atoms_constraints=None):
        """Relax atomic positions as well as unit cell."""
        if atoms_constraints is None:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=cls.OPTCELL_DEFAULT)
        return cls(
            ionmov=cls.IONMOV_DEFAULT,
            optcell=cls.OPTCELL_DEFAULT,
            atoms_constraints=atoms_constraints,
        )

    @property
    def move_atoms(self):
        """True if atoms must be moved."""
        return self.abivars.ionmov != 0

    @property
    def move_cell(self):
        """True if lattice parameters must be optimized."""
        return self.abivars.optcell != 0

    def to_abivars(self):
        """Get a dictionary with the abinit variables."""
        # These variables are always present.
        out_vars = {
            "ionmov": self.abivars.ionmov,
            "optcell": self.abivars.optcell,
            "ntime": self.abivars.ntime,
        }

        # Atom relaxation.
        if self.move_atoms:
            out_vars["tolmxf"] = self.abivars.tolmxf

        if self.abivars.atoms_constraints:
            # Add input variables for constrained relaxation.
            raise NotImplementedError("")
            out_vars.update(self.abivars.atoms_constraints.to_abivars())

        # Cell relaxation.
        if self.move_cell:
            out_vars.update(
                dilatmx=self.abivars.dilatmx,
                ecutsm=self.abivars.ecutsm,
                strfact=self.abivars.strfact,
                strtarget=self.abivars.strtarget,
            )

        return out_vars

    def as_dict(self):
        """Convert to dictionary."""
        dct = dict(self._default_vars)
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dictionary."""
        dct = dct.copy()
        dct.pop("@module", None)
        dct.pop("@class", None)

        return cls(**dct)


@unique
class PPModelModes(Enum):
    """Different kind of plasmon-pole models."""

    noppmodel = 0
    godby = 1
    hybersten = 2
    linden = 3
    farid = 4


class PPModel(AbivarAble, MSONable):
    """
    Parameters defining the plasmon-pole technique.
    The common way to instantiate a PPModel object is via the class method PPModel.as_ppmodel(string).
    """

    @classmethod
    def as_ppmodel(cls, obj):
        """
        Constructs an instance of PPModel from obj.

        Accepts obj in the form:
            * PPmodel instance
            * string. e.g "godby:12.3 eV", "linden".
        """
        if isinstance(obj, cls):
            return obj

        # obj is a string
        if ":" not in obj:
            mode, plasmon_freq = obj, None
        else:
            # Extract mode and plasmon_freq
            mode, plasmon_freq = obj.split(":")
            try:
                plasmon_freq = float(plasmon_freq)
            except ValueError:
                plasmon_freq, unit = plasmon_freq.split()
                plasmon_freq = units.Energy(float(plasmon_freq), unit).to("Ha")

        return cls(mode=mode, plasmon_freq=plasmon_freq)

    def __init__(self, mode="godby", plasmon_freq=None):
        """
        Args:
            mode: ppmodel type
            plasmon_freq: Plasmon frequency in Ha.
        """
        if isinstance(mode, str):
            mode = PPModelModes[mode]
        self.mode = mode
        self.plasmon_freq = plasmon_freq

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("mode", "plasmon_freq")
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented
        other = cast(PPModel, other)

        if self.mode != other.mode:
            return False

        if self.plasmon_freq is None:
            return other.plasmon_freq is None

        return np.allclose(self.plasmon_freq, other.plasmon_freq)

    def __bool__(self):
        return self.mode != PPModelModes.noppmodel

    def __repr__(self):
        return f"<{type(self).__name__} at {id(self)}, mode = {self.mode}>"

    def to_abivars(self):
        """Return dictionary with Abinit variables."""
        if self:
            return {"ppmodel": self.mode.value, "ppmfrq": self.plasmon_freq}
        return {}

    @classmethod
    def get_noppmodel(cls):
        """Calculation without plasmon-pole model."""
        return cls(mode="noppmodel", plasmon_freq=None)

    def as_dict(self):
        """Get JSON-able dict representation."""
        return {
            "mode": self.mode.name,
            "plasmon_freq": self.plasmon_freq,
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build from dict."""
        return cls(mode=dct["mode"], plasmon_freq=dct["plasmon_freq"])


class HilbertTransform(AbivarAble):
    """
    Parameters for the Hilbert-transform method (Screening code)
    i.e. the parameters defining the frequency mesh used for the spectral function
    and the frequency mesh used for the polarizability.
    """

    def __init__(
        self,
        nomegasf,
        domegasf=None,
        spmeth=1,
        nfreqre=None,
        freqremax=None,
        nfreqim=None,
        freqremin=None,
    ):
        """
        Args:
            nomegasf: Number of points for sampling the spectral function along the real axis.
            domegasf: Step in Ha for the linear mesh used for the spectral function.
            spmeth: Algorithm for the representation of the delta function.
            nfreqre: Number of points along the real axis (linear mesh).
            freqremax: Maximum frequency for W along the real axis (in hartree).
            nfreqim: Number of point along the imaginary axis (Gauss-Legendre mesh).
            freqremin: Minimum frequency for W along the real axis (in hartree).
        """
        # Spectral function
        self.nomegasf = nomegasf
        self.domegasf = domegasf
        self.spmeth = spmeth

        # Mesh for the contour-deformation method used for the integration of the self-energy
        self.nfreqre = nfreqre
        self.freqremax = freqremax
        self.freqremin = freqremin
        self.nfreqim = nfreqim

    def to_abivars(self):
        """Get a dictionary with the abinit variables."""
        return {
            # Spectral function
            "nomegasf": self.nomegasf,
            "domegasf": self.domegasf,
            "spmeth": self.spmeth,
            # Frequency mesh for the polarizability
            "nfreqre": self.nfreqre,
            "freqremax": self.freqremax,
            "nfreqim": self.nfreqim,
            "freqremin": self.freqremin,
        }


class ModelDielectricFunction(AbivarAble):
    """Model dielectric function used for BSE calculation."""

    def __init__(self, mdf_epsinf):
        """
        Args:
            mdf_epsinf: Value of epsilon_infinity.
        """
        self.mdf_epsinf = mdf_epsinf

    def to_abivars(self):
        """Return dictionary with abinit variables."""
        return {"mdf_epsinf": self.mdf_epsinf}


##########################################################################################
# WORK IN PROGRESS ######################################
##########################################################################################


class Screening(AbivarAble):
    """
    This object defines the parameters used for the
    computation of the screening function.
    """

    # Approximations used for W
    _WTYPES: ClassVar[dict[str, int]] = {"RPA": 0}

    # Self-consistency modes
    _SC_MODES: ClassVar[dict[str, int]] = {
        "one_shot": 0,
        "energy_only": 1,
        "wavefunctions": 2,
    }

    def __init__(
        self,
        ecuteps,
        nband,
        w_type="RPA",
        sc_mode="one_shot",
        hilbert=None,
        ecutwfn=None,
        inclvkb=2,
    ):
        """
        Args:
            ecuteps: Cutoff energy for the screening (Ha units).
            nband Number of bands for the Green's function
            w_type: Screening type
            sc_mode: Self-consistency mode.
            hilbert: Instance of HilbertTransform defining the parameters for the Hilber transform method.
            ecutwfn: Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).
            inclvkb: Option for the treatment of the dipole matrix elements (NC pseudos).
        """
        if w_type not in self._WTYPES:
            raise ValueError(f"W_TYPE: {w_type} is not supported")

        if sc_mode not in self._SC_MODES:
            raise ValueError(f"Self-consistecy mode {sc_mode} is not supported")

        self.ecuteps = ecuteps
        self.nband = nband
        self.w_type = w_type
        self.sc_mode = sc_mode

        self.ecutwfn = ecutwfn
        self.inclvkb = inclvkb

        if hilbert is not None:
            raise NotImplementedError("Hilber transform not coded yet")
            self.hilbert = hilbert

        # Default values (equivalent to those used in Abinit8)
        self.gwpara = 2
        self.awtr = 1
        self.symchi = 1

        self.optdriver = 3

    @property
    def use_hilbert(self):
        """True if we are using the Hilbert transform method."""
        return hasattr(self, "hilbert")

    # @property
    # def gwcalctyp(self):
    #    "Return the value of the gwcalctyp input variable"
    #    dig0 = str(self._SIGMA_TYPES[self.type])
    #    dig1 = str(self._SC_MODES[self.sc_mode]
    #    return dig1.strip() + dig0.strip()

    def to_abivars(self):
        """Get a dictionary with the abinit variables."""
        abivars = {
            "ecuteps": self.ecuteps,
            "ecutwfn": self.ecutwfn,
            "inclvkb": self.inclvkb,
            "gwpara": self.gwpara,
            "awtr": self.awtr,
            "symchi": self.symchi,
            "nband": self.nband,
            # "gwcalctyp": self.gwcalctyp,
            # "fftgw"    : self.fftgw,
            "optdriver": self.optdriver,
        }

        # Variables for the Hilber transform.
        if self.use_hilbert:
            abivars.update(self.hilbert.to_abivars())

        return abivars


class SelfEnergy(AbivarAble):
    """Define the parameters used for the computation of the self-energy."""

    _SIGMA_TYPES: ClassVar[dict[str, int]] = {
        "gw": 0,
        "hartree_fock": 5,
        "sex": 6,
        "cohsex": 7,
        "model_gw_ppm": 8,
        "model_gw_cd": 9,
    }

    _SC_MODES: ClassVar[dict[str, int]] = {
        "one_shot": 0,
        "energy_only": 1,
        "wavefunctions": 2,
    }

    def __init__(
        self,
        se_type,
        sc_mode,
        nband,
        ecutsigx,
        screening,
        gw_qprange=1,
        ppmodel=None,
        ecuteps=None,
        ecutwfn=None,
        gwpara=2,
    ):
        """
        Args:
            se_type: Type of self-energy (str)
            sc_mode: Self-consistency mode.
            nband: Number of bands for the Green's function
            ecutsigx: Cutoff energy for the exchange part of the self-energy (Ha units).
            screening: Screening instance.
            gw_qprange: Option for the automatic selection of k-points and bands for GW corrections.
                See Abinit docs for more detail. The default value makes the code computie the
                QP energies for all the point in the IBZ and one band above and one band below the Fermi level.
            ppmodel: PPModel instance with the parameters used for the plasmon-pole technique.
            ecuteps: Cutoff energy for the screening (Ha units).
            ecutwfn: Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).
        """
        if se_type not in self._SIGMA_TYPES:
            raise ValueError(f"SIGMA_TYPE: {se_type} is not supported")

        if sc_mode not in self._SC_MODES:
            raise ValueError(f"Self-consistecy mode {sc_mode} is not supported")

        self.type = se_type
        self.sc_mode = sc_mode
        self.nband = nband
        self.ecutsigx = ecutsigx
        self.screening = screening
        self.gw_qprange = gw_qprange
        self.gwpara = gwpara

        if ppmodel is not None:
            if screening.use_hilbert:
                raise ValueError("cannot use hilbert for screening")
            self.ppmodel = PPModel.as_ppmodel(ppmodel)

        self.ecuteps = ecuteps if ecuteps is not None else screening.ecuteps
        self.ecutwfn = ecutwfn
        self.optdriver = 4

        # band_mode in ["gap", "full"]

        # if isinstance(kptgw, str) and kptgw == "all":
        #    self.kptgw = self.nkptgw = None
        # else:
        #    self.kptgw = np.reshape(kptgw, (-1,3))
        #    self.nkptgw =  len(self.kptgw)

        # if bdgw is None:
        #    raise ValueError("bdgw must be specified")

        # if isinstance(bdgw, str):
        #    # TODO add new variable in Abinit so that we can specify
        #    # an energy interval around the KS gap.
        #    homo = float(nele) / 2.0
        #    #self.bdgw =

        # else:
        #    self.bdgw = np.reshape(bdgw, (-1,2))

        # self.freq_int = freq_int

    @property
    def use_ppmodel(self):
        """True if we are using the plasmon-pole approximation."""
        return hasattr(self, "ppmodel")

    @property
    def gwcalctyp(self):
        """The value of the gwcalctyp input variable."""
        dig0 = str(self._SIGMA_TYPES[self.type])
        dig1 = str(self._SC_MODES[self.sc_mode])
        return dig1.strip() + dig0.strip()

    @property
    def symsigma(self):
        """1 if symmetries can be used to reduce the number of q-points."""
        return 1 if self.sc_mode == "one_shot" else 0

    def to_abivars(self):
        """Get a dictionary with the abinit variables."""
        abivars = {
            "gwcalctyp": self.gwcalctyp,
            "ecuteps": self.ecuteps,
            "ecutsigx": self.ecutsigx,
            "symsigma": self.symsigma,
            "gw_qprange": self.gw_qprange,
            "gwpara": self.gwpara,
            "optdriver": self.optdriver,
            "nband": self.nband,
            # "ecutwfn": self.ecutwfn,
            # "kptgw": self.kptgw,
            # "nkptgw": self.nkptgw,
            # "bdgw": self.bdgw,
        }

        # TODO: problem with the spin
        # if len(self.bdgw) != self.nkptgw:
        #    raise ValueError("lengths of bdgw and nkptgw mismatch")

        # ppmodel variables
        if self.use_ppmodel:
            abivars |= self.ppmodel.to_abivars()

        return abivars


class ExcHamiltonian(AbivarAble):
    """Contain parameters for the solution of the Bethe-Salpeter equation."""

    # Types of excitonic Hamiltonian.
    _EXC_TYPES: ClassVar[dict[str, int]] = {
        "TDA": 0,  # Tamm-Dancoff approximation.
        "coupling": 1,  # Calculation with coupling.
    }

    # Algorithms used to compute the macroscopic dielectric function
    # and/or the exciton wavefunctions.
    _ALGO2VAR: ClassVar[dict[str, int]] = {
        "direct_diago": 1,
        "haydock": 2,
        "cg": 3,
    }

    # Options specifying the treatment of the Coulomb term.
    _COULOMB_MODES = ("diago", "full", "model_df")

    def __init__(
        self,
        bs_loband,
        nband,
        mbpt_sciss,
        coulomb_mode,
        ecuteps,
        spin_mode="polarized",
        mdf_epsinf=None,
        exc_type="TDA",
        algo="haydock",
        with_lf=True,
        bs_freq_mesh=None,
        zcut=None,
        **kwargs,
    ):
        r"""
        Args:
            bs_loband: Lowest band index (Fortran convention) used in the e-h  basis set.
                Can be scalar or array of shape (nsppol,). Must be >= 1 and <= nband
            nband: Max band index used in the e-h  basis set.
            mbpt_sciss: Scissors energy in Hartree.
            coulomb_mode: Treatment of the Coulomb term.
            ecuteps: Cutoff energy for W in Hartree.
            mdf_epsinf: Macroscopic dielectric function :math:`\\epsilon_\\inf` used in
                the model dielectric function.
            exc_type: Approximation used for the BSE Hamiltonian
            with_lf: True if local field effects are included <==> exchange term is included
            bs_freq_mesh: Frequency mesh for the macroscopic dielectric function (start, stop, step) in Ha.
            zcut: Broadening parameter in Ha.
            **kwargs:
                Extra keywords.
        """
        spin_mode = SpinMode.as_spinmode(spin_mode)

        # We want an array bs_loband(nsppol).
        try:
            bs_loband = np.reshape(bs_loband, spin_mode.nsppol)
        except ValueError:
            bs_loband = np.array(spin_mode.nsppol * [int(bs_loband)])

        self.bs_loband = bs_loband
        self.nband = nband
        self.mbpt_sciss = mbpt_sciss
        self.coulomb_mode = coulomb_mode
        if coulomb_mode not in self._COULOMB_MODES:
            raise ValueError("coulomb_mode not in _COULOMB_MODES")
        self.ecuteps = ecuteps

        self.mdf_epsinf = mdf_epsinf
        self.exc_type = exc_type
        if exc_type not in self._EXC_TYPES:
            raise ValueError("exc_type not in _EXC_TYPES")
        self.algo = algo
        if algo not in self._ALGO2VAR:
            raise ValueError(f"{algo=} not in {self._ALGO2VAR=}")
        self.with_lf = with_lf

        # If bs_freq_mesh is not given, abinit will select its own mesh.
        self.bs_freq_mesh = np.array(bs_freq_mesh) if bs_freq_mesh is not None else bs_freq_mesh
        self.zcut = zcut
        self.optdriver = 99

        # Extra options.
        self.kwargs = kwargs
        # if "chksymbreak" not in self.kwargs:
        #    self.kwargs["chksymbreak"] = 0

        # Consistency check
        if any(bs_loband < 0):
            raise ValueError(f"bs_loband <= 0 while it is {bs_loband}")
        if any(bs_loband >= nband):
            raise ValueError(f"({bs_loband=}) >= ({nband=})")

    @property
    def inclvkb(self):
        """Treatment of the dipole matrix element (NC pseudos, default is 2)."""
        return self.kwargs.get("inclvkb", 2)

    @property
    def use_haydock(self):
        """True if we are using the Haydock iterative technique."""
        return self.algo == "haydock"

    @property
    def use_cg(self):
        """True if we are using the conjugate gradient method."""
        return self.algo == "cg"

    @property
    def use_direct_diago(self):
        """True if we are performing the direct diagonalization of the BSE Hamiltonian."""
        return self.algo == "direct_diago"

    def to_abivars(self):
        """Get a dictionary with the abinit variables."""
        abivars = {
            "bs_calctype": 1,
            "bs_loband": self.bs_loband,
            # nband=self.nband,
            "mbpt_sciss": self.mbpt_sciss,
            "ecuteps": self.ecuteps,
            "bs_algorithm": self._ALGO2VAR[self.algo],
            "bs_coulomb_term": 21,
            "mdf_epsinf": self.mdf_epsinf,
            "bs_exchange_term": 1 if self.with_lf else 0,
            "inclvkb": self.inclvkb,
            "zcut": self.zcut,
            "bs_freq_mesh": self.bs_freq_mesh,
            "bs_coupling": self._EXC_TYPES[self.exc_type],
            "optdriver": self.optdriver,
        }

        if self.use_haydock:
            abivars.update(
                bs_haydock_niter=100,  # No. of iterations for Haydock
                bs_hayd_term=0,  # No terminator
                bs_haydock_tol=[0.05, 0],  # Stopping criteria
            )

        elif self.use_direct_diago or self.use_cg:
            raise NotImplementedError

        else:
            raise ValueError(f"Unknown algorithm for EXC: {self.algo}")

        # Add extra kwargs
        abivars.update(self.kwargs)

        return abivars
