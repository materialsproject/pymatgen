"Low-level classes providing an abstraction for the objects involved in the calculation."
from __future__ import division, print_function

import collections
import numpy as np
import os.path

from pprint import pprint, pformat

from pymatgen.util.decorators import singleton
from pymatgen.core.design_patterns import Enum
from pymatgen.core.units import any2Ha
from pymatgen.core.physical_constants import Ang2Bohr 
from pymatgen.serializers.json_coders import MSONable 
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure import Structure

##########################################################################################

class AbivarAble(object):
    """
    An AbivarAble object provides a method to_abivars that returns a dictionary with the abinit variables.
    """
    def to_abivars(self):
        raise RuntimeError("%s: must implement the method to_abivars that returns a dictionary with the abinit variables" % (
            self.__class__.__name__))

    def __str__(self):
        return pformat(self.to_abivars(), indent=1, width=80, depth=None)

@singleton
class MandatoryVariable(object):
    """
    Singleton used to tag mandatory variables, just because I can use 
    the cool syntax: variable is MANDATORY!
    """

@singleton
class DefaultVariable(object):
    """
    Singleton used to tag variables that will have the default value
    """

MANDATORY = MandatoryVariable()
DEFAULT   = DefaultVariable()

##########################################################################################

class SpinMode(collections.namedtuple('SpinMode', "mode nsppol nspinor nspden"), AbivarAble):
    """
    Different configurations of the electron density as implemented in abinit:
    One can use asspimode to construct the object via SpinMode.asspinmode(string) 
    where string can assume the values:

        - polarized
        - unpolarized
        - afm (anti-ferromagnetic)
        - spinor (non-collinear magnetism)
        - spinor_nomag (non-collinear, no magnetism)
    """
    @classmethod
    def asspinmode(cls, obj):
        "Converts obj into a SpinMode instance"
        if isinstance(obj, cls):
            return obj
        else:
            # Assume a string with mode
            try:
                return _mode2spinvars[obj]
            except KeyError:
                raise KeyError("Wrong value for spin_mode: %s" % str(obj))

    def to_abivars(self):
        return {"nsppol" : self.nsppol,
                "nspinor": self.nspinor,
                "nspden" : self.nspden,
                }

# An handy Multiton
_mode2spinvars = {
    "unpolarized" : SpinMode("unpolarized" , 1, 1, 1),
    "polarized"   : SpinMode("polarized"   , 2, 1, 2),
    "afm"         : SpinMode("afm"         , 1, 1, 2),
    "spinor"      : SpinMode("spinor"      , 1, 2, 4),
    "spinor_nomag": SpinMode("spinor_nomag", 1, 2, 1),
}

##########################################################################################

class Smearing(AbivarAble, MSONable):
    """
    Variables defining the smearing technique. The preferred way to instanciate 
    a Smearing object is via the class method Smearing.assmearing(string)
    """
    #: Mapping string_mode --> occopt
    _mode2occopt = {  
        'nosmearing' : 1,
        'fermi_dirac': 3,
        'marzari4'   : 4,
        'marzari5'   : 5,
        'methfessel' : 6,
        'gaussian'   : 7,
    }

    modes = Enum(k for k in _mode2occopt)

    def __init__(self, occopt, tsmear): 
        self.occopt = occopt
        self.tsmear = tsmear

    def __str__(self):
        s  = "occopt %d # %s Smearing\n" % (self.occopt, self.mode)
        if self.tsmear:
            s += 'tsmear %s' % self.tsmear
        return s

    def __eq__(self, other):
        if other is None: 
            return False
        else:
            return (self.occopt == other.occopt and self.tsmear == other.tsmear)

    def __ne__(self, other):
        return not self == other

    def __bool__(self):
        return self.mode != "nosmearing"

    # py2 old version
    __nonzero__ = __bool__

    @classmethod
    def assmearing(cls, obj):
        """
        Constructs an instance of Smearing from obj . 
        Accepts obj in the form:

            * Smearing instance
            * "name:tsmear"  e.g. "gaussian:0.004"  (Hartree units)
            * "name:tsmear units" e.g. "gaussian:0.1 eV"
        """
        if isinstance(obj, cls):
            return obj

        # obj is a string
        obj, tsmear = obj.split(":")
        obj.strip()

        if obj == "nosmearing":
            return cls.nosmearing()
        else:
            occopt = cls._mode2occopt[obj]
            try:
                tsmear = float(tsmear)
            except ValueError:
                tsmear, units = tsmear.split()
                tsmear = any2Ha(units)(float(tsmear))

            return cls(occopt, tsmear)

    @property
    def mode(self):
        for (mode_str, occopt) in self._mode2occopt.items():
            if occopt == self.occopt: 
                return mode_str
        raise AttributeError("Unknown occopt %s" % self.occopt)

    @staticmethod
    def nosmearing():
        return Smearing(1, None)

    def to_abivars(self):
        if self.mode == "nosmearing":
            return {}
        else:
            return {"occopt": self.occopt,
                    "tsmear": self.tsmear,}

    @property
    def to_dict(self):
        """json friendly dict representation of Smearing"""
        d = {"occopt": self.occopt, 
             "tsmear": self.tsmear,}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
                                                                     
    @staticmethod
    def from_dict(d):
        return Smearing(d["occopt"], d["tsmear"])

##########################################################################################

class SCFSolver(dict, AbivarAble):
    "Variables controlling the SCF algorithm."

    # None indicates that we use abinit defaults.
    _default_vars = {
        "iprcell"  : None,  
        "iscf"     : None, 
        "diemac"   : None, 
        "diemix"   : None,  
        "diemixmag": None, 
        "dielam"   : None,  
        "diegap"   : None, 
        "dielng"   : None, 
        "diecut"   : None, 
        "nstep"    : 70,
        #"ixc"   ??
        #"toldfe"
        #"tolwfr"
        #"nline"
    }

    def __init__(self, *args, **kwargs):
        super(SCFSolver, self).__init__(*args, **kwargs)

        for k in self.vars:
            if k not in self._default_vars:
                raise ValueError("%s: No default value has been provided for key %s" % (self.__class__.__name__, k))

    def to_abivars(self):
        return self.copy()

##########################################################################################

class AbiStructure(Structure, AbivarAble):
    "Patches the pymatgen structure adding the method to_abivars"

    @classmethod
    def asabistructure(cls, obj):
        """
        Convert obj into an AbiStructure object.
        Accepts:
            - AbiStructure instance
            - Subinstances of pymatgen.
            - File paths
        """
        if isinstance(obj, cls):
            return obj

        if isinstance(obj, Structure):
            # Promote 
            return cls(obj)

        if isinstance(obj, str):
            # Handle file path.
            if os.path.isfile(obj):

                if obj.endswith(".nc"):
                    from .netcdf import structure_from_etsf_file
                    structure = structure_from_etsf_file(obj)
                    print(structure._sites)
                else:
                    from pymatgen.io.smartio import read_structure
                    structure = read_structure(obj)

                # Promote 
                return cls(structure)

        raise ValueError("Don't know how to convert object %s to a AbiStructure structure" % str(obj))

    def __new__(cls, structure):
        new = structure
        new.__class__ = cls
        return new

    def __init__(self, structure):
        pass

    def to_abivars(self):
        "Returns a dictionary with the abinit variables."
        types_of_specie = self.types_of_specie
        natom = self.num_sites

        znucl_type = [specie.number for specie in types_of_specie]
                                                                     
        znucl_atoms = self.atomic_numbers

        typat = np.zeros(natom, np.int)
        for (atm_idx, site) in enumerate(self):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

        significant_figures = 12
        format_str = "{{:.{0}f}}".format(significant_figures)
        fmt = format_str.format

        lines = []
        for vec in Ang2Bohr(self.lattice.matrix):
            lines.append(" ".join([fmt(c) for c in vec]))
        rprim = "\n" + "\n".join(lines)

        lines = []
        for (i, site) in enumerate(self):
            coords = site.frac_coords 
            lines.append( " ".join([fmt(c) for c in coords]) + " # " + site.species_string )
        xred = '\n' + "\n".join(lines)

        abivars = {
            "acell" : 3 * [1.0],
            "rprim" : rprim,
            "natom" : natom,
            "ntypat": len(types_of_specie),
            "typat" : typat,
            "xred"  : xred,
            "znucl" : znucl_type,
        }

        return abivars

##########################################################################################

class KSampling(AbivarAble):
    """
    Input variables defining the K-point sampling.
    """
    #: Modes supported by the constructor.
    modes = Enum(('monkhorst', 'path', 'automatic',))

    def __init__(self, 
                 mode              = "monkhorst",
                 num_kpts          = 0, 
                 kpts              = ((1, 1, 1),), 
                 kpt_shifts        = (0.5, 0.5, 0.5), 
                 use_symmetries    = True,
                 use_time_reversal = True,
                 kpts_weights      = None, 
                 labels            = None, 
                 chksymbreak       = None,
                 comment           = None,
                 ):
        """
        Highly flexible constructor for KSampling objects.  The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility.  For most usage cases, the object be constructed 
        far more easily using the convenience static constructors:
        
            # gamma_only
            # gamma_centered
            # monkhorst
            # monkhorst_automatic
            # path 
        
        and it is recommended that you use those.

        Args:
            mode:
                Mode for generating k-poits. Use one of the KSampling.modes enum types.
            num_kpts:
                Number of kpoints if mode is "automatic"
                Number of division for the sampling of the smallest segment if mode is "path"
                Not used for the other modes
            kpts:
                Number of divisions. Even when only a single specification is
                required, e.g. in the automatic scheme, the kpts should still
                be specified as a 2D array. e.g., [[20]] or [[2,2,2]].
            kpt_shifts:
                Shifts for Kpoints.
            use_symmetries:
                False if spatial symmetries should not be used to reduced the number of independent k-points.
            use_time_reversal:
                False if time-reversal symmetry should not be used to reduced the number of independent k-points.
            kpts_weights:
                Optional weights for kpoints. For explicit kpoints.
            labels:
                In path-mode, this should provide a list of labels for each kpt.
            chksymbreak:
                Abinit input variable: check whether the BZ sampling preserves the symmetry of the crystal.
            comment:
                String comment for Kpoints

        The default behavior of the constructor is monkhorst.
        """
        if mode not in KSampling.modes:
            raise ValueError("Unknown kpoint mode %s" % mode)

        super(KSampling, self).__init__()

        self.comment = comment

        # FIXME
        abivars = {"chksymbreak": 0}

        if mode in ("monkhorst",):
            assert num_kpts == 0
            ngkpt  = np.reshape(kpts, (3,)),
            shiftk = np.reshape(kpt_shifts, (-1,3))

            if use_symmetries and use_time_reversal: kptopt = 1
            if not use_symmetries and use_time_reversal: kptopt = 2
            if not use_symmetries and not use_time_reversal: kptopt = 3
            if use_symmetries and not use_time_reversal: kptopt = 4

            abivars.update({
                "ngkpt"       : ngkpt,
                "shiftk"      : shiftk,
                "nshiftk"     : len(shiftk),
                "kptopt"      : kptopt,
                "chksymbreak" : chksymbreak,
            })

        elif mode in ("path",):
            if num_kpts <= 0:
                raise ValueError("For Path mode, num_kpts must be specified and >0")

            kptbounds = np.reshape(kpts, (-1,3,))

            abivars.update({
                "kptbounds" : kptbounds,
                "ndivsm"    : num_kpts,
                "kptopt"    : -len(kptbounds)+1,
                #"labels"   : labels,  # TODO
            })
        
        elif mode in ("automatic",):
            kpts = np.reshape(kpts, (-1,3))
            if len(kpts) != num_kpts:
                raise ValueError("For Automatic mode, num_kpts must be specified.")

            kptnrm = np.ones(num_kpts) 

            abivars.update({
                "kptopt"      : 0,
                "kpt"         : kpts,
                "nkpt"        : num_kpts,
                "kptnrm"      : kptnrm,
                "wtk"         : kpts_weights,  # for iscf/=-2, wtk.
                "chksymbreak" : chksymbreak,
            })

        else:
            raise ValueError("Unknown mode %s" % mode)

        self.abivars = abivars
        self.abivars["_comment"] = comment

    @classmethod
    def gamma_only(cls):
        "Gamma-only sampling"
        return cls(kpt_shifts=(0.0,0.0,0.0), comment="Gamma-only sampling")

    @classmethod
    def gamma_centered(cls, kpts=(1, 1, 1), use_symmetries=True, use_time_reversal=True):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint grid.

        Args:
            kpts:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
        use_symmetries:
            False if spatial symmetries should not be used to reduced the number of independent k-points.
        use_time_reversal:
            False if time-reversal symmetry should not be used to reduced the number of independent k-points.

        Returns:
            KSampling object
        """
        return cls(kpts              = [kpts],
                   use_symmetries    = use_symmetries,
                   use_time_reversal = use_time_reversal,
                   comment           = "gamma-centered mode", 
                   )

    #@classmethod
    #def symetry_breaking(cls, kpts, shift, use_symmetries=True, use_time_reversal=True):
    #    return cls(
    #                   kpts=[kpts],
    #                   kpt_shifts=shift, 
    #                   use_symmetries = use_symmetries,
    #                   use_time_reversal = use_time_reversal,
    #                   chksymbreak = 0,
    #                   comment="shifted mode", 
    #                   )

    @classmethod
    def monkhorst(cls, ngkpt, shiftk=(0.5, 0.5, 0.5), chksymbreak=None, use_symmetries=True, use_time_reversal=True, comment=None):
        """
        Convenient static constructor for a Monkhorst-Pack mesh.

        Args:
            ngkpt:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            shiftk:
                Shift to be applied to the kpoints. Defaults to (0,0,0).
            use_symmetries:
                Use spatial symmetries to reduce the number of k-points.
            use_time_reversal:
                Use time-reversal symmetry to reduce the number of k-points.

        Returns:
            KSampling object
        """
        return cls(kpts              = [ngkpt],
                   kpt_shifts        = shiftk,
                   use_symmetries    = use_symmetries,
                   use_time_reversal = use_time_reversal,
                   chksymbreak       = chksymbreak,
                   comment           = comment if comment else "Monkhorst-Pack scheme with user-specified shiftk", 
                  )

    @classmethod
    def monkhorst_automatic(cls, structure, ngkpt, chksymbreak=None, use_symmetries=True, use_time_reversal=True, comment=None):
        """
        Convenient static constructor for an automatic Monkhorst-Pack mesh.
                                                                                   
        Args:
            structure:
                paymatgen structure object.
            ngkpt:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            use_symmetries:
                Use spatial symmetries to reduce the number of k-points.
            use_time_reversal:
                Use time-reversal symmetry to reduce the number of k-points.
                                                                                   
        Returns:
            KSampling object
        """
        sg = SymmetryFinder(structure)
        #sg.get_crystal_system()
        #sg.get_point_group()
        # TODO
        nshiftk = 1
        #shiftk = 3*(0.5,) # this is the default
        shiftk = 3*(0.5,)

        #if lattice.ishexagonal:
        #elif lattice.isbcc
        #elif lattice.isfcc

        return cls.monkhorst(ngkpt, 
                             shiftk            = shiftk, 
                             chksymbreak       = chksymbreak,
                             use_symmetries    = use_symmetries, 
                             use_time_reversal = use_time_reversal,
                             comment           = comment if comment else "Automatic Monkhorst-Pack scheme",
                            )

    @classmethod
    def _path(cls, ndivsm, structure=None, kpath_bounds=None, comment=None):
        """
        Static constructor for path in k-space. 
                                                                                   
        Args:
            structure:
                pymatgen structure.
            kpath_bounds:
                List with the reduced coordinates of the k-points defining the path.
            ndivsm:
                Number of division for the smallest segment.
            comment: 
                Comment string.
                                                                                   
        Returns:
            KSampling object
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
                print("label %s, red_coord %s" % (label, red_coord))
                kpath_bounds.append(red_coord)

        return cls(mode     = KSampling.modes.path, 
                   num_kpts = ndivsm,
                   kpts     = kpath_bounds,
                   comment  = comment if comment else "K-Path scheme", 
                  )

    @classmethod
    def path_from_structure(cls, ndivsm, structure):
        "See _path for the meaning of the variables"
        return cls._path(ndivsm,  structure=structure, comment="K-path generated automatically from pymatgen structure")
                                                                                                   
                                                                                                   
    @classmethod
    def explicit_path(cls, ndivsm, kpath_bounds):
        "See _path for the meaning of the variables"
        return cls._path(ndivsm, kpath_bounds=kpath_bounds, comment="Explicit K-path")

    @classmethod
    def automatic_density(cls, structure, kppa, chksymbreak=None, use_symmetries=True, use_time_reversal=True):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes for hexagonal cells and Monkhorst-Pack grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure:
                Input structure
            kppa:
                Grid density
        """
        #raise NotImplementedError()
        #rec_lattice = structure.lattice.reciprocal_lattice

        #min_idx, min_abc = minloc(rec_lattice.abc)
        # See np.argmax
        #ratios = rec_lattice.abc / min_abc

        #kpt_shifts = [0.5, 0.5, 0.5]
        #kpt_shifts = np.atleast_2d(kpt_shifts)

        #num_shifts = len(kpt_shifts)

        #ndiv, num_points = 0, 0

        #while num_points < min_npoints:
        #    ndiv += 1
        #    trial_divs = [int(round(n)) for n in ratios * ndiv]
        #    # ensure that trial_divs  > 0
        #    trial_divs = [i if i > 0 else 1 for i in trial_divs]
        #    num_points = num_shifts * np.product(trial_divs)

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

        right_angles = [i for i in xrange(3) if abs(angles[i] - 90) < hex_angle_tol]

        hex_angles = [i for i in xrange(3)
                      if abs(angles[i] - 60) < hex_angle_tol or
                      abs(angles[i] - 120) < hex_angle_tol]

        is_hexagonal = (len(right_angles) == 2 and len(hex_angles) == 1
                        and abs(lengths[right_angles[0]] -
                                lengths[right_angles[1]]) < hex_length_tol)

        #style = Kpoints.modes.gamma
        #if not is_hexagonal:
        #    num_div = [i + i % 2 for i in num_div]
        #    style = Kpoints.modes.monkhorst

        comment = "pymatgen generated KPOINTS with grid density = " + "{} / atom".format(kppa)

        return cls(mode              = "monkhorst",
                   num_kpts          = 0, 
                   kpts              = [num_div], 
                   kpt_shifts        = [0.5, 0.5, 0.5],
                   chksymbreak       = chksymbreak,
                   use_symmetries    = use_symmetries   , 
                   use_time_reversal = use_time_reversal,
                   comment           = comment, 
                  )

    @property
    def to_dict(self):
        """json friendly dict representation of KSampling"""
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        raise NotImplementedError("")
        return d

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")

    def to_abivars(self):
        return self.abivars.copy()

##########################################################################################

class Electrons(AbivarAble):
    """
    The electronic degrees of freedom
    """
    def __init__(self, 
                 spin_mode = "polarized", 
                 smearing  = "fermi_dirac:0.1 eV",
                 scf_solver = None,
                 nband     = None,
                 fband     = None,
                 charge    = 0.0,
                 comment   = None,
                 #occupancies = None,
                ):
        """
        Constructor for Electrons object.
                                                                                
        Args:
            comment:
                String comment for Electrons

            charge: float
                    Total charge of the system. Default is 0.
        """
        super(Electrons, self).__init__()

        self.comment = comment

        self.smearing = Smearing.assmearing(smearing)

        self.spin_mode = SpinMode.asspinmode(spin_mode)

        self.nband = nband
        self.fband = fband
        self.charge = charge

        self.scf_solver = self.scf_solver

        # FIXME
        if nband is None:
            self.fband = 4

    @property 
    def nsppol(self): 
        return self.spin_mode.nsppol

    @property 
    def nspinor(self): 
        return self.spin_mode.nspinor

    @property 
    def nspden(self): 
        return self.spin_mode.nspden

    @property
    def to_dict(self):
        "json friendly dict representation"
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        raise NotImplementedError("")
        return d

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")

    def to_abivars(self):
        abivars = self.spin_mode.to_abivars()
                                                 
        abivars.update({
            "nband"  : self.nband,
            "fband"  : self.fband,
            "charge" : self.charge,
        })

        if self.smearing:
            abivars.update(self.smearing.to_abivars())

        if self.scf_solver:
            abivars.update(self.scf_solver.to_abivars())

        abivars["_comment"] = comment
        return abivars

##########################################################################################

class Constraints(AbivarAble):
    "Object defining the constraints for structural relaxation"

    def to_abivars(self):
        raise NotImplementedError("")

class RelaxationMethod(AbivarAble):
    """
    Container to store the variables for (constrained) structural optimization
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
        "dilatmx"          : 1.1, 
        "ecutsm"           : 0.5, 
        "strfact"          : None,
        "tolmxf"           : None,  
        "strtarget"        : None,
        "atoms_constraints": {}, # Constraints are stored in a dictionary. Empty if no constraint is enforced"
    }

    IONMOV_DEFAULT = 3
    OPTCELL_DEFAULT = 2

    def __init__(self, **kwargs):

        # Initialize abiars with the default values.
        self.abivars = self._default_vars
                                                                                                                     
        # Overwrite keys with the args and kwargs passed to constructor.
        self.abivars.update(*args, **kwargs)

        self.abivars = AttrDict(self.abivars)
                                                                                                                     
        for k in self.abivars:
            if k not in self._default_vars:
                raise ValueError("%s: No default value has been provided for key %s" % (self.__class__.__name__, k))

        for k in self.abivars:
            if k is MANDATORY:
                raise ValueError("%s: No default value has been provided for the mandatory key %s" % (self.__class__.__name__, k))

    @classmethod
    def atoms_only(cls, atoms_constraints=None):
        if atoms_constraints is None:
            new = cls(ionmov=cls.IONMOV_DEFAULT, optcell=0)
        else:
            new = cls(ionmov=cls.IONMOV_DEFAULT, optcell=0, atoms_constraints=atoms_constraints)
        return new

    @classmethod
    def atoms_and_cell(cls, atoms_constraints=None):
        if atoms_constraints is None:
            new = cls(ionmov=cls.IONMOV_DEFAULT, optcell=cls.OPTCELL_DEFAULT)
        else:
            new = cls(ionmov=cls.IOMOV_DEFAULT, optcell=cls.OPTCELL_DEFAULT, atoms_constraints=atoms_constraints)
        return new

    @property
    def move_atoms(self):
        "True if atoms must be moved"
        return self.vars.ionmov != 0
                                                       
    @property
    def move_cell(self):
        "True if lattice parameters must be optimized"
        return self.vars.optcell != 0

    def to_abivars(self):
        "Returns a dictionary with the abinit variables"
        vars = self.vars

        # These variables are always present.
        abivars = { 
            "ionmov" : vars.ionmov,
            "optcell": vars.optcell,
            "ntime"  : vars.ntime,
        }

        # Atom relaxation.
        if self.move_atoms:
            abivars.update({
                "tolmxf": vars.tolmxf,
            })

        # Cell relaxation.
        if self.move_cell:
            abivars.update({
                "dilatmx"  : vars.dilatmx,
                "ecutsm"   : vars.ecutsm,
                "strfact"  : vars.strfact,
                "strtarget": vars.strtarget,
            })

        if vars.atoms_constraints:
            # Add input variables for constrained relaxation.
            raise NotImplementedError("")
            abivars.update(vars.atoms_constraints.to_abivars())

        return abivars

##########################################################################################

class PPModel(AbivarAble, MSONable):
    """
    Parameters defining the plasmon-pole technique.
    The common way to instanciate a PPModel object is via the class method PPModel.asppmodel(string)
    """

    _mode2ppmodel = {  
        "noppmodel": 0,
        "godby"    : 1,
        "hybersten": 2,
        "linden"   : 3,
        "farid"    : 4,
    }

    modes = Enum(k for k in _mode2ppmodel)

    @classmethod
    def asppmodel(cls, obj):
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
                plasmon_freq, units = plasmon_freq.split()
                plasmon_freq = any2Ha(units)(float(plasmon_freq))

        return cls(mode=mode, plasmon_freq=plasmon_freq)

    def __init__(self, mode="godby", plasmon_freq=None):
        assert mode in PPModel.modes
        self.mode = mode
        self.plasmon_freq = plasmon_freq

    def __bool__(self):
        return self.mode != "noppmodel"

    # py2 old version
    __nonzero__ = __bool__

    def __repr__(self):
        return "<%s at %s, mode = %s>" % (self.__class__.__name__, id(self), str(self.mode))

    def to_abivars(self):
        if self:
            return {"ppmodel": self._mode2ppmodel[self.mode],
                    "ppmfrq" : self.plasmon_freq,}
        else:
            return {}

    @classmethod
    def noppmodel(cls):
        return cls(mode=modes.noppmodel, plasmon_freq=None)

    @property
    def to_dict(self):
        d = {"mode": self.mode, "plasmon_freq": plasmon_freq}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return PPmodel(mode=d["mode"], plasmon_freq=d["plasmon_freq"])

##########################################################################################
#################################  WORK IN PROGRESS ######################################
##########################################################################################

class Screening(AbivarAble):
    # FIXME
    _types = {
        "gw"           : 0,
        "hartree_fock" : 5,
        "sex"          : 6,
        "cohsex"       : 7,
        "model_gw_ppm" : 8,
        "model_gw_cd"  : 9,
    }
                                    
    types = Enum(_types)
                                     
    _sc_modes = {
        "one_shot"     : 0,
        "energy_only"  : 1,
        "wavefunctions": 2,
    }
                                    
    scmodes = Enum(_sc_modes.keys())

    # Example
    #ecuteps = 2
    #scr_strategy = ScreeningStrategy(fftgw=11)
    #w = Screening(w_type, sc_mode, ecuteps, ecutwfn=None, ppmodel=None, freq_mesh=None)
    #scr_strategy.set_w(w)

    def __init__(self, w_type, sc_mode, ecuteps, ecutwfn=None, ppmodel=None, freq_mesh=None):
        self.type = w_type
        self.sc_mode = sc_mode
        self.ecuteps = ecuteps
        self.ecutwfn = ecutwfn

        if [obj is not None for obj in [ppmodel, freq_mesh]].count(True) != 1:
            raise ValueError("Either ppmodel or freq_mesh can be specified") 

        if ppmodel is not None:
            self.ppmodel = PPmodel.asppmodel(ppmodel)
        else:
            raise NotImplementedError("")
            self.wmesh = WMesh.aswmesh(wmesh)

    @property
    def ppmfreq(self):
        pass

    def to_abivars(self):
        "Returns a dictionary with the abinit variables"
        abivars = { 
            "ecuteps"  : vars.ecuteps,
            "ecutwfn"  : vars.ecutwfn,
            #"gwpara"   : vars.gwpara,
            #"awtr"     : vars.awtr,
            #"symchi"   : vars.symchi,
            #"gwmem"    : vars.gwmem,
            #"gwcalctyp": vars.gwcalctyp,
            #"fftgw"    : vars.fftgw,
            #"inclvkb"  : vars.inclvkb,
        }

        if self.has_ppmodel:
            abivars.update(vars.ppmodel.to_abivars())

        # Variables for the Hilber transform.
        if self.hilbert_transform:
            hilbert = {
                "spmeth"  : vars.spmeth,
                "nomegasf": vars.nomegasf,
            }
            abivars.update(hilbert)

        return abivars                                                                                       

##########################################################################################

class SelfEnergy(AbivarAble):
    _types = {
        "gw"          : 0,
        "hartree_fock": 5,
        "sex"         : 6,
        "cohsex"      : 7,
        "model_gw_ppm": 8,
        "model_gw_cd" : 9,
    }
                                    
    types = Enum(_types)
                                     
    _sc_modes = {
        "one_shot"     : 0,
        "energy_only"  : 1,
        "wavefunctions": 2,
    }
                                    
    scmodes = Enum(_sc_modes.keys())

    # Example
    #ecuteps, nband = 2, 20
    #self_energy = SelfEnergy("gw", "one_shot", ppmodel="godby", ecuteps, nband)
    #se_strategy = SelfEnergyStrategy(fftgw=11)
    #se_strategy.learn_data(data)
    #se_strategy.set_self_energy(self_energy)

    def __init__(self, se_type, sc_mode, ecuteps, nband, 
                 ecutwfn=None, ecutsigx=None, kptgw="all", bdgw=None, ppmodel="noppmodel", freq_int=None):
        "freq_int: cd for contour deformation"

        self.type     = se_type
        self.sc_mode  = sc_mode
        self.ecutwfn  = ecutwfn
        self.ecutsigx = ecutsigx
        self.ecuteps  = ecuteps
        self.nband    = nband
        #band_mode in ["gap", "full"]

        if isinstance(kptgw, str) and kptgw == "all":
            self.kptgw = None
            self.nkptgw = None
        else:
            self.kptgw = np.reshape(kptgw, (-1,3))
            self.nkptgw =  len(self.kptgw)

        if bdgw is None:
            raise ValueError("bdgw must be specified")

        if isinstance(bdgw, str):
            # TODO add new variable in Abinit so that we can specify 
            # an energy interval around the KS gap.
            homo = float(nele) / 2.0
            #self.bdgw = 

        else:
            self.bdgw = np.reshape(bdgw, (-1,2))

        self.freq_int = freq_int

        self.ppmodel = PPModel.asppmodel(ppmodel)

    @property
    def gwcalctyp(self):
        "Return the value of the gwcalctyp input variable"
        dig0 = str(SelfEnergy._types[self.type])
        dig1 = str(SelfEnergy._sc_modes[self.sc_mode])
        return dig1.strip() + dig0.strip()

    @property
    def symsigma(self):
        return 1 if self.sc_mode == "one_shot" else 0

    def nband_green(self):
        return self.nband

    def to_abivars(self):
        "Returns a dictionary with the abinit variables"

        abivars = { 
            "gwcalctyp": self.gwcalctyp,
            "ecuteps"  : self.ecuteps,
            "ecutwfn"  : self.ecutwfn,
            "ecutsigx" : self.ecutsigx,
            "symsigma" : self.symsigma,
            "kptgw"    : self.kptgw, 
            "nkptgw"   : self.nkptgw,
            "bdgw"     : self.bdgw,
        }

        # FIXME: problem with the spin
        #assert len(self.bdgw) == self.nkptgw

        # ppmodel variables
        if self.ppmodel:
            abivars.update(self.ppmodel.to_abivars())

        return abivars                                                                                       

##########################################################################################

class CDFrequencyMesh(AbivarAble):
    "Mesh for the contour-deformation method used for the integration of the self-energy"

    #    #def from_slice
    #    #def from_linspace
    #    #FrequencyMesh(real_slice=slice(0,15,1), units="Ha")
    #
    #    def __init__(self, r_slice, nfreqim=10, freqim_alpha=5, units="Ha"):
             #np.arange(start, s.stop, step))
    #        self.r_slice = r_slice
    #        self.units = units

    #
    #    @property
    #    def num_rpoints(self)
    #        "Number of points along the real axis"
    #
    #    @property
    #    def num_ipoint(self)
    #        "Number of points along the imaginary axis"
    #
    #    @property
    #    def r_step(self)
    #        "Step used to sample the real axis"
    #

    def to_abivars(self):
        "Returns a dictionary with the abinit variables"
        abivars = {
            "freqremax"   : self.freqremax,
            "freqremin"   : self.freqremin,
            "nfreqre"     : self.nfreqre,
            "nfreqim"     : self.nfreqim,
            "freqim_alpha": self.freqim_alpha
        }
        return abivars

##########################################################################################

class HilbertTransform(AbivarAble):
    "Parameters for the Hilbert-transform method (RPA, Screening code)"

    def __init__(self, nomegasf, spmeth=1):
        self.spmeth = spmeth
        self.nomegasf = nomegasf

    def to_abivars(self):
        return {"spmeth"  : self.spmeth,
                "nomegasf": self.nomegasf,}

##########################################################################################

class ModelDielectricFunction(AbivarAble):
    "Model dielectric function used for BSE calculation"

    def __init__(self, mdf_epsinf):
        self.mdf_epsinf = mdf_epsinf

    def to_abivars(self):
        return {"mdf_epsinf": self.mdf_epsinf}

##########################################################################################
