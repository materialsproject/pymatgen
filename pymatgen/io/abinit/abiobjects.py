# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Low-level objects providing an abstraction for the objects involved in the calculation.
"""
from __future__ import unicode_literals, division, print_function

import collections
import abc
import six
import numpy as np
import pymatgen.core.units as units

from pprint import pformat
from monty.design_patterns import singleton
from monty.collections import AttrDict
from enum import Enum
from monty.json import MSONable
from pymatgen.core.units import ArrayWithUnit
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.serializers.json_coders import pmg_serialize
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.json import MontyEncoder, MontyDecoder


def lattice_from_abivars(cls=None, *args, **kwargs):
    """
    Returns a `Lattice` object from a dictionary
    with the Abinit variables `acell` and `rprim` in Bohr.
    If acell is not given, the Abinit default is used i.e. [1,1,1] Bohr

    Args:
        cls: Lattice class to be instantiated. pymatgen.core.lattice.Lattice if cls is None

    Example:

    lattice_from_abivars(acell=3*[10], rprim=np.eye(3))
    """
    cls = Lattice if cls is None else cls
    kwargs.update(dict(*args))
    d = kwargs

    rprim = d.get("rprim", None)
    angdeg = d.get("angdeg", None)
    acell = d["acell"]

    # Call pymatgen constructors (note that pymatgen uses Angstrom instead of Bohr).
    if rprim is not None:
        assert angdeg is None
        rprim = np.reshape(rprim, (3,3))
        rprimd = [float(acell[i]) * rprim[i] for i in range(3)]
        return cls(ArrayWithUnit(rprimd, "bohr").to("ang"))

    elif angdeg is not None:
        # angdeg(0) is the angle between the 2nd and 3rd vectors,
        # angdeg(1) is the angle between the 1st and 3rd vectors,
        # angdeg(2) is the angle between the 1st and 2nd vectors,
        raise NotImplementedError("angdeg convention should be tested")
        angles = angdeg
        angles[1] = -angles[1]
        l = ArrayWithUnit(acell, "bohr").to("ang")
        return cls.from_lengths_and_angles(l, angdeg)
    else:
        raise ValueError("Don't know how to construct a Lattice from dict: %s" % str(d))


def structure_from_abivars(cls=None, *args, **kwargs):
    """
    Build a :class:`Structure` object from a dictionary with ABINIT variables.

    Args:
        cls: Structure class to be instantiated. pymatgen.core.structure.Structure if cls is None

    example:

        al_structure = structure_from_abivars(
            acell=3*[7.5],
            rprim=[0.0, 0.5, 0.5,
                   0.5, 0.0, 0.5,
                   0.5, 0.5, 0.0],
            typat=1,
            xred=[0.0, 0.0, 0.0],
            ntypat=1,
            znucl=13,
        )

    `xred` can be replaced with `xcart` or `xangst`.
    """
    kwargs.update(dict(*args))
    d = kwargs

    cls = Structure if cls is None else cls

    lattice = Lattice.from_dict(d, fmt="abivars")
    coords, coords_are_cartesian = d.get("xred", None), False

    if coords is None:
        coords = d.get("xcart", None)
        if coords is not None:
            if "xangst" in d:
                raise ValueError("xangst and xcart are mutually exclusive")
            coords = ArrayWithUnit(coords, "bohr").to("ang")
        else:
            coords = d.get("xangst", None)
        coords_are_cartesian = True

    if coords is None:
        raise ValueError("Cannot extract coordinates from:\n %s" % str(d))

    coords = np.reshape(coords, (-1,3))

    znucl_type, typat = d["znucl"], d["typat"]

    if not isinstance(znucl_type, collections.Iterable):
        znucl_type = [znucl_type]

    if not isinstance(typat, collections.Iterable):
        typat = [typat]

    assert len(typat) == len(coords)

    # Note Fortran --> C indexing
    #znucl_type = np.rint(znucl_type)
    species = [znucl_type[typ-1] for typ in typat]

    return cls(lattice, species, coords, validate_proximity=False,
               to_unit_cell=False, coords_are_cartesian=coords_are_cartesian)


def structure_to_abivars(structure, **kwargs):
    """Returns a dictionary with the ABINIT variables."""
    types_of_specie = structure.types_of_specie
    natom = structure.num_sites

    znucl_type = [specie.number for specie in types_of_specie]

    znucl_atoms = structure.atomic_numbers

    typat = np.zeros(natom, np.int)
    for atm_idx, site in enumerate(structure):
        typat[atm_idx] = types_of_specie.index(site.specie) + 1

    rprim = ArrayWithUnit(structure.lattice.matrix, "ang").to("bohr")
    xred = np.reshape([site.frac_coords for site in structure], (-1,3))

    # Set small values to zero. This usually happens when the CIF file
    # does not give structure parameters with enough digits.
    rprim = np.where(np.abs(rprim) > 1e-8, rprim, 0.0)
    xred = np.where(np.abs(xred) > 1e-8, xred, 0.0)

    # Info on atoms.
    d = dict(
        natom=natom,
        ntypat=len(types_of_specie),
        typat=typat,
        znucl=znucl_type,
        xred=xred,
    )

    # Add info on the lattice.
    # Should we use (rprim, acell) or (angdeg, acell) to specify the lattice?
    geomode = kwargs.pop("geomode", "rprim")
    #latt_dict = structure.lattice.to_abivars(geomode=geomode)

    if geomode == "rprim":
        d.update(dict(
            acell=3 * [1.0],
            rprim=rprim))

    elif geomode == "angdeg":
        d.update(dict(
            acell=3 * [1.0],
            angdeg=angdeg))
    else:
        raise ValueError("Wrong value for geomode: %s" % geomode)

    return d


def contract(s):
    """
    >>> assert contract("1 1 1 2 2 3") == "3*1 2*2 1*3"
    >>> assert contract("1 1 3 2 3") == "2*1 1*3 1*2 1*3"
    """
    if not s: return s

    tokens = s.split()
    old = tokens[0]
    count = [[1, old]]

    for t in tokens[1:]:
        if t == old:
            count[-1][0] += 1
        else:
            old = t
            count.append([1, t])

    return " ".join("%d*%s" % (c, t) for c, t in count)


class AbivarAble(six.with_metaclass(abc.ABCMeta, object)):
    """
    An `AbivarAble` object provides a method `to_abivars`
    that returns a dictionary with the abinit variables.
    """

    @abc.abstractmethod
    def to_abivars(self):
        """Returns a dictionary with the abinit variables."""

    #@abc.abstractmethod
    #def from_abivars(cls, vars):
    #    """Build the object from a dictionary with Abinit variables."""

    def __str__(self):
        return pformat(self.to_abivars(), indent=1, width=80, depth=None)

    def __contains__(self, key):
        return key in self.to_abivars()


@singleton
class MandatoryVariable(object):
    """
    Singleton used to tag mandatory variables, just because I can use
    the cool syntax: variable is MANDATORY!
    """


@singleton
class DefaultVariable(object):
    """Singleton used to tag variables that will have the default value"""

MANDATORY = MandatoryVariable()
DEFAULT = DefaultVariable()


class SpinMode(collections.namedtuple('SpinMode', "mode nsppol nspinor nspden"), AbivarAble, MSONable):
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
    @classmethod
    def as_spinmode(cls, obj):
        """Converts obj into a `SpinMode` instance"""
        if isinstance(obj, cls):
            return obj
        else:
            # Assume a string with mode
            try:
                return _mode2spinvars[obj]
            except KeyError:
                raise KeyError("Wrong value for spin_mode: %s" % str(obj))

    def to_abivars(self):
        return {
            "nsppol": self.nsppol,
            "nspinor": self.nspinor,
            "nspden": self.nspden,
        }

    @pmg_serialize
    def as_dict(self):
        return {k: getattr(self, k) for k in self._fields}

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: d[k] for k in d if k in cls._fields})


# An handy Multiton
_mode2spinvars = {
    "unpolarized": SpinMode("unpolarized", 1, 1, 1),
    "polarized": SpinMode("polarized", 2, 1, 2),
    "afm": SpinMode("afm", 1, 1, 2),
    "spinor": SpinMode("spinor", 1, 2, 4),
    "spinor_nomag": SpinMode("spinor_nomag", 1, 2, 1),
}


class Smearing(AbivarAble, MSONable):
    """
    Variables defining the smearing technique. The preferred way to instanciate
    a `Smearing` object is via the class method Smearing.as_smearing(string)
    """
    #: Mapping string_mode --> occopt
    _mode2occopt = {
        'nosmearing': 1,
        'fermi_dirac': 3,
        'marzari4': 4,
        'marzari5': 5,
        'methfessel': 6,
        'gaussian': 7}

    def __init__(self, occopt, tsmear):
        self.occopt = occopt
        self.tsmear = tsmear

    def __str__(self):
        s = "occopt %d # %s Smearing\n" % (self.occopt, self.mode)
        if self.tsmear:
            s += 'tsmear %s' % self.tsmear
        return s

    def __eq__(self, other):
        return (self.occopt == other.occopt and
                np.allclose(self.tsmear, other.tsmear))

    def __ne__(self, other):
        return not self == other

    def __bool__(self):
        return self.mode != "nosmearing"

    # py2 old version
    __nonzero__ = __bool__

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
        else:
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
        for (mode_str, occopt) in self._mode2occopt.items():
            if occopt == self.occopt:
                return mode_str
        raise AttributeError("Unknown occopt %s" % self.occopt)

    @staticmethod
    def nosmearing():
        return Smearing(1, 0.0)

    def to_abivars(self):
        if self.mode == "nosmearing":
            return {"occopt": 1, "tsmear": 0.0}
        else:
            return {"occopt": self.occopt, "tsmear": self.tsmear,}

    @pmg_serialize
    def as_dict(self):
        """json friendly dict representation of Smearing"""
        return {"occopt": self.occopt, "tsmear": self.tsmear}

    @staticmethod
    def from_dict(d):
        return Smearing(d["occopt"], d["tsmear"])


class ElectronsAlgorithm(dict, AbivarAble, MSONable):
    """Variables controlling the SCF/NSCF algorithm."""
    # None indicates that we use abinit defaults.
    _DEFAULT = dict(
        iprcell=None, iscf=None, diemac=None, diemix=None, diemixmag=None,
        dielam=None, diegap=None, dielng=None, diecut=None, nstep=50)

    def __init__(self, *args, **kwargs):
        super(ElectronsAlgorithm, self).__init__(*args, **kwargs)

        for k in self:
            if k not in self._DEFAULT:
                raise ValueError("%s: No default value has been provided for "
                                 "key %s" % (self.__class__.__name__, k))

    def to_abivars(self):
        return self.copy()

    @pmg_serialize
    def as_dict(self):
        return self.copy()

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        d.pop("@module", None)
        d.pop("@class", None)
        return cls(**d)


class Electrons(AbivarAble, MSONable):
    """The electronic degrees of freedom"""
    def __init__(self, spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                 algorithm=None, nband=None, fband=None, charge=0.0, comment=None):  # occupancies=None,
        """
        Constructor for Electrons object.

        Args:
            comment: String comment for Electrons
            charge: Total charge of the system. Default is 0.
        """
        super(Electrons, self).__init__()

        self.comment = comment
        self.smearing = Smearing.as_smearing(smearing)
        self.spin_mode = SpinMode.as_spinmode(spin_mode)

        self.nband = nband
        self.fband = fband
        self.charge = charge
        self.algorithm = algorithm

    @property
    def nsppol(self):
        return self.spin_mode.nsppol

    @property
    def nspinor(self):
        return self.spin_mode.nspinor

    @property
    def nspden(self):
        return self.spin_mode.nspden

    def as_dict(self):
        "json friendly dict representation"
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["spin_mode"] = self.spin_mode.as_dict()
        d["smearing"] = self.smearing.as_dict()
        d["algorithm"] = self.algorithm.as_dict() if self.algorithm else None
        d["nband"] = self.nband
        d["fband"] = self.fband
        d["charge"] = self.charge
        d["comment"] = self.comment
        return d

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        d.pop("@module", None)
        d.pop("@class", None)
        dec = MontyDecoder()
        d["spin_mode"] = dec.process_decoded(d["spin_mode"])
        d["smearing"] = dec.process_decoded(d["smearing"])
        d["algorithm"] = dec.process_decoded(d["algorithm"]) if d["algorithm"] else None
        return cls(**d)


    def to_abivars(self):
        abivars = self.spin_mode.to_abivars()

        abivars.update({
            "nband"  : self.nband,
            "fband"  : self.fband,
            "charge" : self.charge,
        })

        if self.smearing:
            abivars.update(self.smearing.to_abivars())

        if self.algorithm:
            abivars.update(self.algorithm)

        #abivars["#comment"] = self.comment
        return abivars

class KSamplingModes(Enum):
    monkhorst = 1
    path = 2
    automatic = 3


class KSampling(AbivarAble, MSONable):
    """
    Input variables defining the K-point sampling.
    """
    # Modes supported by the constructor.


    def __init__(self, mode=KSamplingModes.monkhorst, num_kpts= 0,
                 kpts=((1, 1, 1),),
                 kpt_shifts=(0.5, 0.5, 0.5),
                 kpts_weights=None, use_symmetries=True, use_time_reversal=True, chksymbreak=None,
                 comment=None):
        """
        Highly flexible constructor for KSampling objects.  The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility.  For most usage cases, the object be constructed
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
                  be specified as a 2D array. e.g., [[20]] or [[2,2,2]].
            kpt_shifts: Shifts for Kpoints.
            use_symmetries: False if spatial symmetries should not be used
                to reduce the number of independent k-points.
            use_time_reversal: False if time-reversal symmetry should not be used
                to reduce the number of independent k-points.
            kpts_weights: Optional weights for kpoints. For explicit kpoints.
            chksymbreak: Abinit input variable: check whether the BZ sampling preserves the symmetry of the crystal.
            comment: String comment for Kpoints

        .. note::
            The default behavior of the constructor is monkhorst.
        """
        if isinstance(mode, six.string_types):
            mode = KSamplingModes[mode]

        super(KSampling, self).__init__()

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
            assert num_kpts == 0
            ngkpt  = np.reshape(kpts, 3)
            shiftk = np.reshape(kpt_shifts, (-1,3))

            if use_symmetries and use_time_reversal: kptopt = 1
            if not use_symmetries and use_time_reversal: kptopt = 2
            if not use_symmetries and not use_time_reversal: kptopt = 3
            if use_symmetries and not use_time_reversal: kptopt = 4

            abivars.update({
                "ngkpt"      : ngkpt,
                "shiftk"     : shiftk,
                "nshiftk"    : len(shiftk),
                "kptopt"     : kptopt,
                "chksymbreak": chksymbreak,
            })

        elif mode == KSamplingModes.path:
            if num_kpts <= 0:
                raise ValueError("For Path mode, num_kpts must be specified and >0")

            kptbounds = np.reshape(kpts, (-1,3))
            #print("in path with kptbound: %s " % kptbounds)

            abivars.update({
                "ndivsm"   : num_kpts,
                "kptbounds": kptbounds,
                "kptopt"   : -len(kptbounds)+1,
            })

        elif mode == KSamplingModes.automatic:
            kpts = np.reshape(kpts, (-1,3))
            if len(kpts) != num_kpts:
                raise ValueError("For Automatic mode, num_kpts must be specified.")

            kptnrm = np.ones(num_kpts)

            abivars.update({
                "kptopt"     : 0,
                "kpt"        : kpts,
                "nkpt"       : num_kpts,
                "kptnrm"     : kptnrm,
                "wtk"        : kpts_weights,  # for iscf/=-2, wtk.
                "chksymbreak": chksymbreak,
            })

        else:
            raise ValueError("Unknown mode %s" % mode)

        self.abivars = abivars
        #self.abivars["#comment"] = comment

    @property
    def is_homogeneous(self):
        return self.mode not in ["path"]

    @classmethod
    def gamma_only(cls):
        """Gamma-only sampling"""
        return cls(kpt_shifts=(0.0,0.0,0.0), comment="Gamma-only sampling")

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
            :class:`KSampling` object.
        """
        return cls(kpts=[kpts], kpt_shifts=(0.0, 0.0, 0.0),
                   use_symmetries=use_symmetries, use_time_reversal=use_time_reversal,
                   comment="gamma-centered mode")

    @classmethod
    def monkhorst(cls, ngkpt, shiftk=(0.5, 0.5, 0.5), chksymbreak=None, use_symmetries=True,
                  use_time_reversal=True, comment=None):
        """
        Convenient static constructor for a Monkhorst-Pack mesh.

        Args:
            ngkpt: Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            shiftk: Shift to be applied to the kpoints.
            use_symmetries: Use spatial symmetries to reduce the number of k-points.
            use_time_reversal: Use time-reversal symmetry to reduce the number of k-points.

        Returns:
            :class:`KSampling` object.
        """
        return cls(
            kpts=[ngkpt], kpt_shifts=shiftk,
            use_symmetries=use_symmetries, use_time_reversal=use_time_reversal, chksymbreak=chksymbreak,
            comment=comment if comment else "Monkhorst-Pack scheme with user-specified shiftk")

    @classmethod
    def monkhorst_automatic(cls, structure, ngkpt,
                            use_symmetries=True, use_time_reversal=True, chksymbreak=None, comment=None):
        """
        Convenient static constructor for an automatic Monkhorst-Pack mesh.

        Args:
            structure: :class:`Structure` object.
            ngkpt: Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            use_symmetries: Use spatial symmetries to reduce the number of k-points.
            use_time_reversal: Use time-reversal symmetry to reduce the number of k-points.

        Returns:
            :class:`KSampling` object.
        """
        sg = SpacegroupAnalyzer(structure)
        #sg.get_crystal_system()
        #sg.get_point_group_symbol()
        # TODO
        nshiftk = 1
        #shiftk = 3*(0.5,) # this is the default
        shiftk = 3*(0.5,)

        #if lattice.ishexagonal:
        #elif lattice.isbcc
        #elif lattice.isfcc

        return cls.monkhorst(
            ngkpt, shiftk=shiftk, use_symmetries=use_symmetries, use_time_reversal=use_time_reversal,
            chksymbreak=chksymbreak, comment=comment if comment else "Automatic Monkhorst-Pack scheme")

    @classmethod
    def _path(cls, ndivsm, structure=None, kpath_bounds=None, comment=None):
        """
        Static constructor for path in k-space.

        Args:
            structure: :class:`Structure` object.
            kpath_bounds: List with the reduced coordinates of the k-points defining the path.
            ndivsm: Number of division for the smallest segment.
            comment: Comment string.

        Returns:
            :class:`KSampling` object.
        """
        if kpath_bounds is None:
            # Compute the boundaries from the input structure.
            from pymatgen.symmetry.bandstructure import HighSymmKpath
            sp = HighSymmKpath(structure)

            # Flat the array since "path" is a a list of lists!
            kpath_labels = []
            for labels in sp.kpath["path"]:
                kpath_labels.extend(labels)

            kpath_bounds = []
            for label in kpath_labels:
                red_coord = sp.kpath["kpoints"][label]
                #print("label %s, red_coord %s" % (label, red_coord))
                kpath_bounds.append(red_coord)

        return cls(mode=KSamplingModes.path, num_kpts=ndivsm, kpts=kpath_bounds,
                   comment=comment if comment else "K-Path scheme")

    @classmethod
    def path_from_structure(cls, ndivsm, structure):
        """See _path for the meaning of the variables"""
        return cls._path(ndivsm,  structure=structure, comment="K-path generated automatically from structure")

    @classmethod
    def explicit_path(cls, ndivsm, kpath_bounds):
        """See _path for the meaning of the variables"""
        return cls._path(ndivsm, kpath_bounds=kpath_bounds, comment="Explicit K-path")

    @classmethod
    def automatic_density(cls, structure, kppa, chksymbreak=None, use_symmetries=True, use_time_reversal=True,
                          shifts=(0.5, 0.5, 0.5)):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
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
        ngrid = kppa / structure.num_sites

        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3.)

        num_div = [int(round(1.0 / lengths[i] * mult)) for i in range(3)]
        # ensure that num_div[i] > 0
        num_div = [i if i > 0 else 1 for i in num_div]

        angles = lattice.angles
        hex_angle_tol = 5      # in degrees
        hex_length_tol = 0.01  # in angstroms

        right_angles = [i for i in range(3) if abs(angles[i] - 90) < hex_angle_tol]

        hex_angles = [i for i in range(3)
                      if abs(angles[i] - 60) < hex_angle_tol or
                      abs(angles[i] - 120) < hex_angle_tol]

        is_hexagonal = (len(right_angles) == 2 and len(hex_angles) == 1
                        and abs(lengths[right_angles[0]] -
                                lengths[right_angles[1]]) < hex_length_tol)

        #style = KSamplingModes.gamma
        #if not is_hexagonal:
        #    num_div = [i + i % 2 for i in num_div]
        #    style = KSamplingModes.monkhorst

        comment = "pymatge.io.abinit generated KPOINTS with grid density = " + "{} / atom".format(kppa)

        shifts = np.reshape(shifts, (-1, 3))

        return cls(
            mode="monkhorst", num_kpts=0, kpts=[num_div], kpt_shifts=shifts,
            use_symmetries=use_symmetries, use_time_reversal=use_time_reversal, chksymbreak=chksymbreak,
            comment=comment)

    def to_abivars(self):
        return self.abivars

    def as_dict(self):
        enc = MontyEncoder()
        return {'mode': self.mode.name, 'comment': self.comment,
                'num_kpts': self.num_kpts,
                'kpts': enc.default(np.array(self.kpts)), 'kpt_shifts': self.kpt_shifts,
                'kpts_weights': self.kpts_weights, 'use_symmetries': self.use_symmetries,
                'use_time_reversal': self.use_time_reversal, 'chksymbreak': self.chksymbreak,
                '@module': self.__class__.__module__, '@class': self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        d.pop('@module', None)
        d.pop('@class', None)
        dec = MontyDecoder()
        d['kpts'] = dec.process_decoded(d['kpts'])
        return cls(**d)


class Constraints(AbivarAble):
    """This object defines the constraints for structural relaxation"""
    def to_abivars(self):
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
    _default_vars =  {
        "ionmov"           : MANDATORY,
        "optcell"          : MANDATORY,
        "ntime"            : 80,
        "dilatmx"          : 1.05,
        "ecutsm"           : 0.5,
        "strfact"          : None,
        "tolmxf"           : None,
        "strtarget"        : None,
        "atoms_constraints": {}, # Constraints are stored in a dictionary. {} means if no constraint is enforced.
    }

    IONMOV_DEFAULT = 3
    OPTCELL_DEFAULT = 2

    def __init__(self, *args, **kwargs):

        # Initialize abivars with the default values.
        self.abivars = self._default_vars

        # Overwrite the keys with the args and kwargs passed to constructor.
        self.abivars.update(*args, **kwargs)

        self.abivars = AttrDict(self.abivars)

        for k in self.abivars:
            if k not in self._default_vars:
                raise ValueError("%s: No default value has been provided for key %s" % (self.__class__.__name__, k))

        for k in self.abivars:
            if k is MANDATORY:
                raise ValueError("%s: No default value has been provided for the mandatory key %s" %
                                 (self.__class__.__name__, k))

    @classmethod
    def atoms_only(cls, atoms_constraints=None):
        if atoms_constraints is None:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=0)
        else:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=0, atoms_constraints=atoms_constraints)

    @classmethod
    def atoms_and_cell(cls, atoms_constraints=None):
        if atoms_constraints is None:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=cls.OPTCELL_DEFAULT)
        else:
            return cls(ionmov=cls.IONMOV_DEFAULT, optcell=cls.OPTCELL_DEFAULT, atoms_constraints=atoms_constraints)

    @property
    def move_atoms(self):
        """True if atoms must be moved."""
        return self.abivars.ionmov != 0

    @property
    def move_cell(self):
        """True if lattice parameters must be optimized."""
        return self.abivars.optcell != 0

    def to_abivars(self):
        """Returns a dictionary with the abinit variables"""
        # These variables are always present.
        out_vars = {
            "ionmov" : self.abivars.ionmov,
            "optcell": self.abivars.optcell,
            "ntime"  : self.abivars.ntime,
        }

        # Atom relaxation.
        if self.move_atoms:
            out_vars.update({
                "tolmxf": self.abivars.tolmxf,
            })

        if self.abivars.atoms_constraints:
            # Add input variables for constrained relaxation.
            raise NotImplementedError("")
            out_vars.update(self.abivars.atoms_constraints.to_abivars())

        # Cell relaxation.
        if self.move_cell:
            out_vars.update({
                "dilatmx"  : self.abivars.dilatmx,
                "ecutsm"   : self.abivars.ecutsm,
                "strfact"  : self.abivars.strfact,
                "strtarget": self.abivars.strtarget,
            })

        return out_vars

    def as_dict(self):
        d = dict(self._default_vars)
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        d.pop('@module', None)
        d.pop('@class', None)

        return cls(**d)


class PPModelModes(Enum):
    noppmodel = 0
    godby = 1
    hybersten = 2
    linden = 3
    farid = 4


class PPModel(AbivarAble, MSONable):
    """
    Parameters defining the plasmon-pole technique.
    The common way to instanciate a PPModel object is via the class method PPModel.as_ppmodel(string)
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
        if isinstance(mode, six.string_types):
            mode = PPModelModes[mode]
        self.mode = mode
        self.plasmon_freq = plasmon_freq

    def __eq__(self, other):
        if other is None:
            return False
        else:
            if self.mode != other.mode:
                return False

            if self.plasmon_freq is None:
                return other.plasmon_freq is None
            else:
                return np.allclose(self.plasmon_freq, other.plasmon_freq)

    def __ne__(self, other):
        return not self == other

    def __bool__(self):
        return self.mode != PPModelModes.noppmodel

    # py2 old version
    __nonzero__ = __bool__

    def __repr__(self):
        return "<%s at %s, mode = %s>" % (self.__class__.__name__, id(self),
                                          str(self.mode))

    def to_abivars(self):
        if self:
            return {"ppmodel": self.mode.value,
                    "ppmfrq": self.plasmon_freq}
        else:
            return {}

    @classmethod
    def get_noppmodel(cls):
        return cls(mode="noppmodel", plasmon_freq=None)

    def as_dict(self):
        return {"mode": self.mode.name, "plasmon_freq": self.plasmon_freq,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @staticmethod
    def from_dict(d):
        return PPModel(mode=d["mode"], plasmon_freq=d["plasmon_freq"])


class HilbertTransform(AbivarAble):
    """
    Parameters for the Hilbert-transform method (Screening code)
    i.e. the parameters defining the frequency mesh used for the spectral function
    and the frequency mesh used for the polarizability
    """
    def __init__(self, nomegasf, domegasf=None, spmeth=1, nfreqre=None, freqremax=None, nfreqim=None, freqremin=None):
        """
        Args:
            nomegasf: Number of points for sampling the spectral function along the real axis.
            domegasf: Step in Ha for the linear mesh used for the spectral function.
            spmeth: Algorith for the representation of the delta function.
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
        """Returns a dictionary with the abinit variables"""
        return {
                # Spectral function
                "nomegasf": self.nomegasf,
                "domegasf": self.domegasf,
                "spmeth"  : self.spmeth,
                 # Frequency mesh for the polarizability
                "nfreqre"  : self.nfreqre,
                "freqremax": self.freqremax,
                "nfreqim"  : self.nfreqim,
                "freqremin": self.freqremin,
                }


class ModelDielectricFunction(AbivarAble):
    """Model dielectric function used for BSE calculation"""
    def __init__(self, mdf_epsinf):
        self.mdf_epsinf = mdf_epsinf

    def to_abivars(self):
        return {"mdf_epsinf": self.mdf_epsinf}

##########################################################################################
#################################  WORK IN PROGRESS ######################################
##########################################################################################


class Screening(AbivarAble):
    """
    This object defines the parameters used for the
    computation of the screening function.
    """
    # Approximations used for W
    _WTYPES = {
        "RPA": 0,
    }

    # Self-consistecy modes
    _SC_MODES = {
        "one_shot"     : 0,
        "energy_only"  : 1,
        "wavefunctions": 2,
    }

    def __init__(self, ecuteps, nband, w_type="RPA", sc_mode="one_shot",
                 hilbert=None, ecutwfn=None, inclvkb=2):
        """
        Args:
            ecuteps: Cutoff energy for the screening (Ha units).
            nband Number of bands for the Green's function
            w_type: Screening type
            sc_mode: Self-consistency mode.
            hilbert: Instance of :class:`HilbertTransform` defining the parameters for the Hilber transform method.
            ecutwfn: Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).
            inclvkb: Option for the treatment of the dipole matrix elements (NC pseudos).
        """
        if w_type not in self._WTYPES:
            raise ValueError("W_TYPE: %s is not supported" % w_type)

        if sc_mode not in self._SC_MODES:
            raise ValueError("Self-consistecy mode %s is not supported" % sc_mode)

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
        self.gwpara=2
        self.awtr  =1
        self.symchi=1

        self.optdriver = 3

    @property
    def use_hilbert(self):
        return hasattr(self, "hilbert")

    #@property
    #def gwcalctyp(self):
    #    "Return the value of the gwcalctyp input variable"
    #    dig0 = str(self._SIGMA_TYPES[self.type])
    #    dig1 = str(self._SC_MODES[self.sc_mode]
    #    return dig1.strip() + dig0.strip()

    def to_abivars(self):
        """Returns a dictionary with the abinit variables"""
        abivars = {
            "ecuteps"   : self.ecuteps,
            "ecutwfn"   : self.ecutwfn,
            "inclvkb"   : self.inclvkb,
            "gwpara"    : self.gwpara,
            "awtr"      : self.awtr,
            "symchi"    : self.symchi,
            #"gwcalctyp": self.gwcalctyp,
            #"fftgw"    : self.fftgw,
            "optdriver" : self.optdriver,
        }

        # Variables for the Hilber transform.
        if self.use_hilbert:
            abivars.update(self.hilbert.to_abivars())

        return abivars


class SelfEnergy(AbivarAble):
    """
    This object defines the parameters used for the computation of the self-energy.
    """
    _SIGMA_TYPES = {
        "gw"          : 0,
        "hartree_fock": 5,
        "sex"         : 6,
        "cohsex"      : 7,
        "model_gw_ppm": 8,
        "model_gw_cd" : 9,
    }

    _SC_MODES = {
        "one_shot"     : 0,
        "energy_only"  : 1,
        "wavefunctions": 2,
    }

    def __init__(self, se_type, sc_mode, nband, ecutsigx, screening,
                 gw_qprange=1, ppmodel=None, ecuteps=None, ecutwfn=None, gwpara=2):
        """
        Args:
            se_type: Type of self-energy (str)
            sc_mode: Self-consistency mode.
            nband: Number of bands for the Green's function
            ecutsigx: Cutoff energy for the exchange part of the self-energy (Ha units).
            screening: :class:`Screening` instance.
            gw_qprange: Option for the automatic selection of k-points and bands for GW corrections.
                See Abinit docs for more detail. The default value makes the code computie the
                QP energies for all the point in the IBZ and one band above and one band below the Fermi level.
            ppmodel: :class:`PPModel` instance with the parameters used for the plasmon-pole technique.
            ecuteps: Cutoff energy for the screening (Ha units).
            ecutwfn: Cutoff energy for the wavefunctions (Default: ecutwfn == ecut).
        """
        if se_type not in self._SIGMA_TYPES:
            raise ValueError("SIGMA_TYPE: %s is not supported" % se_type)

        if sc_mode not in self._SC_MODES:
            raise ValueError("Self-consistecy mode %s is not supported" % sc_mode)

        self.type = se_type
        self.sc_mode = sc_mode
        self.nband = nband
        self.ecutsigx = ecutsigx
        self.screening = screening
        self.gw_qprange = gw_qprange
        self.gwpara = gwpara

        if ppmodel is not None:
            assert not screening.use_hilbert
            self.ppmodel = PPModel.as_ppmodel(ppmodel)

        self.ecuteps = ecuteps if ecuteps is not None else screening.ecuteps
        self.ecutwfn = ecutwfn
        self.optdriver = 4

        #band_mode in ["gap", "full"]

        #if isinstance(kptgw, str) and kptgw == "all":
        #    self.kptgw = None
        #    self.nkptgw = None
        #else:
        #    self.kptgw = np.reshape(kptgw, (-1,3))
        #    self.nkptgw =  len(self.kptgw)

        #if bdgw is None:
        #    raise ValueError("bdgw must be specified")

        #if isinstance(bdgw, str):
        #    # TODO add new variable in Abinit so that we can specify
        #    # an energy interval around the KS gap.
        #    homo = float(nele) / 2.0
        #    #self.bdgw =

        #else:
        #    self.bdgw = np.reshape(bdgw, (-1,2))

        #self.freq_int = freq_int

    @property
    def use_ppmodel(self):
        """True if we are using the plasmon-pole approximation."""
        return hasattr(self, "ppmodel")

    @property
    def gwcalctyp(self):
        """Returns the value of the gwcalctyp input variable."""
        dig0 = str(self._SIGMA_TYPES[self.type])
        dig1 = str(self._SC_MODES[self.sc_mode])
        return dig1.strip() + dig0.strip()

    @property
    def symsigma(self):
        """1 if symmetries can be used to reduce the number of q-points."""
        return 1 if self.sc_mode == "one_shot" else 0

    def to_abivars(self):
        """Returns a dictionary with the abinit variables."""
        abivars = dict(
            gwcalctyp=self.gwcalctyp,
            ecuteps=self.ecuteps,
            ecutsigx=self.ecutsigx,
            symsigma=self.symsigma,
            gw_qprange=self.gw_qprange,
            gwpara=self.gwpara,
            optdriver=self.optdriver,
            #"ecutwfn"  : self.ecutwfn,
            #"kptgw"    : self.kptgw,
            #"nkptgw"   : self.nkptgw,
            #"bdgw"     : self.bdgw,
        )

        # FIXME: problem with the spin
        #assert len(self.bdgw) == self.nkptgw

        # ppmodel variables
        if self.use_ppmodel:
            abivars.update(self.ppmodel.to_abivars())

        return abivars


class ExcHamiltonian(AbivarAble):
    """This object contains the parameters for the solution of the Bethe-Salpeter equation."""
    # Types of excitonic Hamiltonian.
    _EXC_TYPES = {
        "TDA": 0,          # Tamm-Dancoff approximation.
        "coupling": 1,     # Calculation with coupling.
    }

    # Algorithms used to compute the macroscopic dielectric function
    # and/or the exciton wavefunctions.
    _ALGO2VAR = {
        "direct_diago": 1,
        "haydock"     : 2,
        "cg"          : 3,
    }

    # Options specifying the treatment of the Coulomb term.
    _COULOMB_MODES = [
        "diago",
        "full",
        "model_df"
        ]

    def __init__(self, bs_loband, nband, mbpt_sciss, coulomb_mode, ecuteps, spin_mode="polarized", mdf_epsinf=None,
                 exc_type="TDA", algo="haydock", with_lf=True, bs_freq_mesh=None, zcut=None, **kwargs):
        """
        Args:
            bs_loband: Lowest band index (Fortran convention) used in the e-h  basis set.
                Can be scalar or array of shape (nsppol,). Must be >= 1 and <= nband
            nband: Max band index used in the e-h  basis set.
            mbpt_sciss: Scissors energy in Hartree.
            coulomb_mode: Treatment of the Coulomb term.
            ecuteps: Cutoff energy for W in Hartree.
            mdf_epsinf: Macroscopic dielectric function :math:`\epsilon_\inf` used in
                the model dielectric function.
            exc_type: Approximation used for the BSE Hamiltonian
            with_lf: True if local field effects are included <==> exchange term is included
            bs_freq_mesh: Frequency mesh for the macroscopic dielectric function (start, stop, step) in Ha.
            zcut: Broadening parameter in Ha.
            **kwargs:
                Extra keywords
        """
        spin_mode = SpinMode.as_spinmode(spin_mode)

        # We want an array bs_loband(nsppol).
        try:
            bs_loband = np.reshape(bs_loband, spin_mode.nsppol)
        except ValueError:
            bs_loband = np.array(spin_mode.nsppol * [int(bs_loband)])

        self.bs_loband = bs_loband
        self.nband  = nband
        self.mbpt_sciss = mbpt_sciss
        self.coulomb_mode = coulomb_mode
        assert coulomb_mode in self._COULOMB_MODES
        self.ecuteps = ecuteps

        self.mdf_epsinf = mdf_epsinf
        self.exc_type = exc_type
        assert exc_type in self._EXC_TYPES
        self.algo = algo
        assert algo in self._ALGO2VAR
        self.with_lf = with_lf

        # if bs_freq_mesh is not given, abinit will select its own mesh.
        self.bs_freq_mesh = np.array(bs_freq_mesh) if bs_freq_mesh is not None else bs_freq_mesh
        self.zcut = zcut
        self.optdriver = 99

        # Extra options.
        self.kwargs = kwargs
        #if "chksymbreak" not in self.kwargs:
        #    self.kwargs["chksymbreak"] = 0

        # Consistency check
        if any(bs_loband < 0):
            raise ValueError("bs_loband <= 0 while it is %s" % bs_loband)
        if any(bs_loband >= nband):
            raise ValueError("bs_loband (%s) >= nband (%s)" % (bs_loband, nband))

    @property
    def inclvkb(self):
        """Treatment of the dipole matrix element (NC pseudos, default is 2)"""
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
        """Returns a dictionary with the abinit variables."""
        abivars = dict(
            bs_calctype=1,
            bs_loband=self.bs_loband,
            #nband=self.nband,
            mbpt_sciss=self.mbpt_sciss,
            ecuteps=self.ecuteps,
            bs_algorithm=self._ALGO2VAR[self.algo],
            bs_coulomb_term=21,
            mdf_epsinf=self.mdf_epsinf,
            bs_exchange_term=1 if self.with_lf else 0,
            inclvkb=self.inclvkb,
            zcut=self.zcut,
            bs_freq_mesh=self.bs_freq_mesh,
            bs_coupling=self._EXC_TYPES[self.exc_type],
            optdriver=self.optdriver,
        )

        if self.use_haydock:
            # FIXME
            abivars.update(
                bs_haydock_niter=100,      # No. of iterations for Haydock
                bs_hayd_term=0,            # No terminator
                bs_haydock_tol=[0.05, 0],  # Stopping criteria
            )

        elif self.use_direct_diago:
            raise NotImplementedError("")

        elif self.use_cg:
            raise NotImplementedError("")

        else:
            raise ValueError("Unknown algorithm for EXC: %s" % self.algo)

        # Add extra kwargs
        abivars.update(self.kwargs)

        return abivars
