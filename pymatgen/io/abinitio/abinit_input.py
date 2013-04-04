"""
Classes for manipulating/writing ABINIT input files. 
"""
from __future__ import division, print_function

import sys
import os
import os.path
import warnings
import collections
import numpy as np

from pprint import pprint

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.design_patterns import Enum
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr
from pymatgen.core.units import any2Ha
from pymatgen.util.string_utils import str_aligned, str_delimited
from pymatgen.serializers.json_coders import MSONable #, PMGJSONDecoder
from pymatgen.symmetry.finder import SymmetryFinder

from .pseudos import PseudoTable

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

##########################################################################################

class AbinitCard(dict, MSONable):
    """
    Base class representing a set of Abinit input parameters associated to a particular topic. 
    Essentially consists of a dictionary with some helper functions.

    The input variables are stored in the dictionary as {"varname" : value}, whereas instance
    attributes are stored in self[_attributes].
    The __str__ method returns a string with the corresponding section of the abinit input file.
    """
    VALID_VARS = tuple()

    def __init__(self, attributes=None, **kwargs):
        """
        Creates an AbinitCard object.

        Args:
            attributes:
                dictionary with extra attributes to be stored in self["_attributes"]
            kwargs:
                The input variables as dictionary {varname : varvalue}
        """
        super(AbinitCard, self).__init__(**kwargs)

        self["_attributes"] = {}
        if attributes is not None:
            self["_attributes"].update(attributes)

    def __str__(self):
        "String representation (return the section of the abinit input file)"
        return self.get_string()
                                                                                          
    def __getattribute__(self, name):
        try:
            # Default behaviour
            return super(AbinitCard, self).__getattribute__(name)
        except AttributeError:
            try:
                # Try in the dictionary.
                return self[name]
            except KeyError:
                # Try in _attributes.
                try:
                    return self["_attributes"][name]
                except KeyError as exc:
                    raise AttributeError(str(exc))
                                                                                          
    def __setattr__(self, name, value):
        # Don't break private attributes e.g self.__class__
        if name.startswith("__"): 
            super(AbinitCard, self).__setattr__(name, value)
        else:
            self["_attributes"].update({name : value})

    def __add__(self, other):
        return self.merge(other)

    def copy(self):
        "Return a shallow copy of self."
        kwargs =  super(AbinitCard, self).copy()
        attributes = kwargs.pop("_attributes", None)
        new = AbinitCard(attributes, **kwargs)
        new.__class__ = self.__class__
        return new

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def comment(self):
        "String comment"
        try:
            return str(self._comment)
        except:
            return "No comment available"

    #def __setitem__(self, key, val):
    #    """
    #    Add parameter-value pair to self.  Warns if parameter is not in list of
    #    valid tags. Also cleans the parameter and val by stripping
    #    leading and trailing white spaces.
    #    """
    #    if self.VALID_VARS and key.strip() not in self.VALID_VARS:
    #        warnings.warn(key + " not in VALID_VARS")
    #    super(AbinitCard, self).__setitem__(key.strip(),
    #                                        self.proc_val(key.strip(), val.strip())
    #                                        if isinstance(val, str) else val)

    @property
    def to_dict(self):
        raise NotImplementedError("Subclasses should provide the to_dict method")

    @classmethod
    def from_dict(cls, d):
        raise NotImplementedError("Subclasses should provide the from_dict method")

    def get_string(self, sort_keys=False, pretty=False):
        """
        Returns a string representation of self. The reason why this
        method is different from the __str__ method is to provide options for pretty printing.

        Args:
            sort_keys:
                Set to True to sort the AbinitCard parameters alphabetically. Defaults to False.
            pretty:
                Set to True for pretty aligned output. Defaults to False.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)

        lines = []
        if self.comment: 
            lines.append(["# comment: ", self.comment])
            lines.append(["#", "\n"])

        for k in keys:
            if k == "_attributes": continue # Don't print the attributes.

            if self.VALID_VARS and k.strip() not in self.VALID_VARS:
                warnings.warn("key:" + k + " not in VALID_VARS")

            value = self[k]
            if value is None: continue

            if isinstance(value, collections.Iterable) and not isinstance(value, basestring):
                arr = np.array(value)

                if len(arr.shape) in [0,1]: 
                    # scalar or vector.
                    lines.append([k, " ".join([str(i) for i in arr])])

                else: 
                    # array --> matrix 
                    matrix = np.reshape(arr, (-1, arr.shape[-1])) 
                    for (idx, row) in enumerate(matrix):
                        kname = k +"\n" if idx == 0 else ""
                        lines.append([kname, " ".join([str(i) for i in row])])

            else:
                lines.append([k, value])

        if pretty:
            return str_aligned(lines, header=None)
        else:
            return str_delimited(lines, header=None, delimiter=5*" ")

    def diff(self, other):
        """
        Diff function for AbinitCard.  Compares two objects and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other:
                The other AbinitCard object to compare to.

        Returns:
            Dict of the following format:
            {"Same" : parameters_that_are_the_same,
            "Different": parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {"natom":3}
        """
        similar_param, different_param = {}, {}
        for (k1, v1) in self.items():
            if k1 not in other:
                different_param[k1] = {self.name : v1, other.name : "Default"}
            elif v1 != other[k1]:
                different_param[k1] = {self.name : v1, other.name : other[k1]}
            else:
                similar_param[k1] = v1

        for (k2, v2) in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {self.name : "Default", other.name : v2}

        return {"Same": similar_param, "Different": different_param}

    def merge(self, other):
        """
        Add all the values of another AbinitCard object to this object.
        """
        if self.name != other.name:
            raise ValueError("Cannot merge cards of different type: %s, %s!" % (self.name, other.name) )

        s_items = {k : v for (k, v) in self.items()}
        s_attributes = s_items.pop("_attributes")

        o_items = {k : v for (k, v) in other.items()}
        o_attributes = o_items.pop("_attributes")

        new_dict = s_items.copy()
        for (k, v) in o_items:
            if k in new_dict and new_dict[k] != v:
                raise ValueError("AbinitCards have conflicting variables!")
            else:
                new_dict[k] = v

        new_attributes = s_attributes.copy()
        for (k, v) in o_attributes:
            if k in new_attributes and new_attributes[k] != v:
                raise ValueError("AbinitCards have conflicting attributes!")
            else:
                new_attributes[k] = v

        new = AbinitCard(new_attributes, **new_dict)
        # Return same type as self.
        new.__class__ = self.__class__
        return new

##########################################################################################

class System(AbinitCard):
    """
    Object for representing the geometry of the system and the list of pseudopotentials

    .. attribute:: structure
        Associated Structure.

    .. attribute:: comment
        Optional comment string.
    """
    VALID_VARS = [
        "acell", 
        "rprim",
        "natom",
        "ntypat",
        "typat",
        "znucl",
        "xred",
        "nsym",
    ]

    def __init__(self, structure, pseudos, nsym=None, comment=None): 
        """
        Args:
            structure:
                Input pymatgen structure.
            pseudos:
                List of pseudopotentials.
            nsym:
                Number of space group symmetries used by abinit. Set it to 1 to disable the use of symmetries.
            comment:
                Comment string.
        """
        if not structure.is_ordered:
            raise ValueError("Structure with partial occupancies cannot be converted into System!")

        super(System, self).__init__()

        self.structure = structure
        self.nsym = nsym

        self._comment = structure.formula if comment is None else comment

        if not isinstance(pseudos, PseudoTable):
            pseudos = PseudoTable(pseudos)

        # Extract pseudos for this calculation (needed when we pass an entier periodic table)
        table = pseudos
        pseudos = []

        for typ in structure.types_of_specie:
            # Get list of pseudopotentials in table from atom symbol.
            pseudos_for_type = table.pseudos_with_symbol(typ)
                                                                                 
            if pseudos_for_type is None or len(pseudos_for_type) != 1:
                raise ValueError("Cannot find unique pseudo for type %s" % typ)
                                                                                 
            pseudos.append(pseudos_for_type[0])
                                                                                     
        self.pseudos = PseudoTable(pseudos)

        # TODO
        #self.num_valence_electrons 
                                                                              
        types_of_specie = structure.types_of_specie
                                                                     
        natom = structure.num_sites, 
        ntypat = len(types_of_specie)
                                                                     
        znucl_type = [specie.number for specie in types_of_specie]
                                                                     
        znucl_atoms = structure.atomic_numbers

        typat = np.zeros(natom, np.int)
        for (atm_idx, site) in enumerate(structure):
            typat[atm_idx] = types_of_specie.index(site.specie) + 1

        significant_figures = 12
        format_str = "{{:.{0}f}}".format(significant_figures)
        fmt = format_str.format

        lines = []
        for vec in Ang2Bohr(structure.lattice.matrix):
            lines.append(" ".join([fmt(c) for c in vec]))
        rprim = "\n" + "\n".join(lines)

        lines = []
        for (i, site) in enumerate(structure):
            coords = site.frac_coords 
            lines.append( " ".join([fmt(c) for c in coords]) + " # " + site.species_string )
        xred = '\n' + "\n".join(lines)

        self.update({
            "acell"  : 3 * [1.0],
            "rprim"  : rprim,
            "natom"  : natom, 
            "ntypat" : ntypat,
            "typat"  : typat,
            "znucl"  : znucl_type,
            "xred"   : xred,
            "nsym"   : nsym,
            #"pseudo_list" : ", ".join(p.path for p in self.pseudos),
        })


    @property
    def isnc(self):
        "True if norm-conserving calculation"
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        "True if PAW calculation"
        return all(p.ispaw for p in self.pseudos)

    # TODO
    #@property
    #def num_valence_electrons(self):
    #    for specie in self.structure:
    #        print specie
    #        specie.pseudo
    #    return sum(p.num_valence_electrons for p in self.pseudos)

    @staticmethod
    def boxed_molecule(pseudos, cart_coords, acell=3*(10,)):
        """
        Creates a molecule in a periodic box of lengths acell [Bohr]

        Args:
            pseudos:
                List of pseudopotentials
            cart_coords: 
                Cartesian coordinates
            acell:
                Lengths of the box in *Bohr*
        """
        cart_coords = np.atleast_2d(cart_coords)
                                                                                 
        molecule = Molecule([p.symbol for p in pseudos], cart_coords)
                                                                                 
        l = Bohr2Ang(acell)
                                                                                 
        structure = molecule.get_boxed_structure(l[0], l[1], l[2])
                                                                                 
        comment = structure.formula + " in a periodic box of size %s [Bohr]" % str(acell)
                                                                                 
        return System(structure, pseudos, comment=comment)

    @staticmethod
    def boxed_atom(pseudo, cart_coords=3*(0,), acell=3*(10,)):
        """
        Creates an atom in a periodic box of lengths acell [Bohr]
                                                                     
        Args:
            pseudo:
                Pseudopotential object.
            cart_coords: 
                Cartesian coordinates
            acell:
                Lengths of the box in *Bohr*
        """
        return System.boxed_molecule([pseudo], cart_coords, acell=acell)

    @property
    def to_dict(self):
        return {"@module"   : self.__class__.__module__,
                "@class"    : self.__class__.__name__,
                "structure" : self.structure,
                "pseudos"   : self.pseudos,
                "nsym"      : self.nsym,
                "comment"   : self._comment,
        }

    @staticmethod
    def from_dict(d):
        return System(
            d["structure"],
            d["pseudos"],
            d["nsym"],
            comment=d["comment"],
        )

##########################################################################################

class Kpoints(AbinitCard):
    """
    Input variables defining the K-point sampling.
    """
    #: Modes supported by the constructor.
    modes = Enum(('monkhorst', 'path', 'automatic',))

    VALID_VARS = [
        "kptopt",
        "ngkpt",
        "nshiftk",
        "shiftk",
        "kptbounds",
        "ndivsm",
        "nkpt",
        "kpt",
        "kptnrm",
        "wtk",
        "chksymbreak",
    ]

    def __init__(self, 
                 mode              = modes.monkhorst,
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
        Highly flexible constructor for Kpoints object.  The flexibility comes
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
                Mode for generating k-poits. Use one of the Kpoints.modes enum types.
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
        if mode not in Kpoints.modes:
            raise ValueError("Unknown kpoint mode %s" % mode)

        super(Kpoints, self).__init__()

        self._comment = comment

        # FIXME
        chksymbreak = 0

        if mode in ("monkhorst",):
            assert num_kpts == 0
            ngkpt  = np.reshape(kpts, (3,)),
            shiftk = np.reshape(kpt_shifts, (-1,3))

            if use_symmetries and use_time_reversal: kptopt = 1
            if not use_symmetries and use_time_reversal: kptopt = 2
            if not use_symmetries and not use_time_reversal: kptopt = 3
            if use_symmetries and not use_time_reversal: kptopt = 4

            self.update({
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

            self.update({
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

            self.update({
                "kptopt"      : 0,
                "kpt"         : kpts,
                "nkpt"        : num_kpts,
                "kptnrm"      : kptnrm,
                "wtk"         : kpts_weights,  # for iscf/=-2, wtk.
                "chksymbreak" : chksymbreak,
            })

        else:
            raise ValueError("Unknown mode %s" % mode)

    @staticmethod
    def gamma_only():
        "Gamma-only sampling"
        return Kpoints(kpt_shifts=(0.0,0.0,0.0), comment="Gamma-only sampling")

    @staticmethod
    def gamma_centered(kpts=(1, 1, 1), use_symmetries=True, use_time_reversal=True):
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
            Kpoints object
        """
        return Kpoints(
                       kpts              = [kpts],
                       use_symmetries    = use_symmetries,
                       use_time_reversal = use_time_reversal,
                       comment           = "gamma-centered mode", 
                       )

    #@staticmethod
    #def symetry_breaking(kpts, shift, use_symmetries=True, use_time_reversal=True):
    #    return Kpoints(
    #                   kpts=[kpts],
    #                   kpt_shifts=shift, 
    #                   use_symmetries = use_symmetries,
    #                   use_time_reversal = use_time_reversal,
    #                   chksymbreak = 0,
    #                   comment="shifted mode", 
    #                   )

    @staticmethod
    def monkhorst(ngkpt, shiftk=(0.5, 0.5, 0.5), chksymbreak=None, use_symmetries=True, use_time_reversal=True, comment=None):
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
            Kpoints object
        """
        return Kpoints(
                       kpts              = [ngkpt],
                       kpt_shifts        = shiftk,
                       use_symmetries    = use_symmetries,
                       use_time_reversal = use_time_reversal,
                       chksymbreak       = chksymbreak,
                       comment           = comment if comment else "Monkhorst-Pack scheme with user-specified shiftk", 
                       )

    @staticmethod
    def monkhorst_automatic(structure, ngkpt, chksymbreak=None, use_symmetries=True, use_time_reversal=True, comment=None):
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
            Kpoints object
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

        return Kpoints.monkhorst(ngkpt, 
                                 shiftk            = shiftk, 
                                 chksymbreak       = chksymbreak,
                                 use_symmetries    = use_symmetries, 
                                 use_time_reversal = use_time_reversal,
                                 comment           = comment if comment else "Automatic Monkhorst-Pack scheme",
                                 )

    @staticmethod
    def _path(ndivsm, structure=None, kpath_bounds=None, comment=None):
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
            Kpoints object
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

        return Kpoints(
                      mode      = Kpoints.modes.path, 
                      num_kpts  = ndivsm,
                      kpts      = kpath_bounds,
                      comment   = comment if comment else "K-Path scheme", 
                      )

    @staticmethod
    def path_from_structure(ndivsm, structure):
        "See _path for the meaning of the variables"
        return Kpoints._path(ndivsm, 
                             structure=structure, 
                             comment="K-path generated automatically from pymatgen structure")
                                                                                                   
                                                                                                   
    @staticmethod
    def explicit_path(ndivsm, kpath_bounds):
        "See _path for the meaning of the variables"
        return Kpoints._path(ndivsm, 
                             kpath_bounds=kpath_bounds,
                             comment="Explicit K-path")


    @staticmethod
    def automatic_density(structure, kppa, chksymbreak=None, use_symmetries=True, use_time_reversal=True):
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

        return Kpoints(mode              = "monkhorst",
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
        """json friendly dict representation of Kpoints"""
        raise NotImplementedError("")
        #d = {"comment": self.comment, 
        #     "nkpoints": self.num_kpts,
        #     "generation_style": self.style, 
        #     "kpoints": self.kpts,
        #     "usershift": self.kpt_shift
        #    }
        #d["@module"] = self.__class__.__module__
        #d["@class"] = self.__class__.__name__
        #return d

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")
        #return Kpoints(
        #         mode              = modes.monkhorst,
        #         num_kpts          = 0, 
        #         kpts              = ((1, 1, 1),), 
        #         kpt_shifts        = (0, 0, 0), 
        #         use_symmetries    = True,
        #         use_time_reversal = True,
        #         kpts_weights      = None, 
        #         labels            = None, 
        #         chksymbreak       = None,
        #         comment           = None,
        #         ):

##########################################################################################

class Smearing(MSONable):
    "Container with the abinit variables defining the smearing technique."

    #: Mapping string_mode --> occopt
    _mode2occopt = {  
        'None'        : 1,
        'fermi_dirac' : 3,
        'marzari4'    : 4,
        'marzari5'    : 5,
        'methfessel'  : 6,
        'gaussian'    : 7,
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
        if other is None: return False 
        return ( self.occopt == other.occopt and
                 self.tsmear == other.tsmear)

    def __ne__(self, other):
        return not self == other

    @classmethod
    def assmearing(cls, obj):
        """
        Constructs an instance of Smearing from obj . 
        Accepts obj in the form:

            * "name:tsmear"  e.g. "gaussian:0.004"  (Hartree units)
            * "name:tsmear units" e.g. "gaussian:0.1 eV"
            * None
            * Smearing instance
        """
        if isinstance(obj, cls):
            return obj

        if obj is None:
            return Smearing.NoSmearing()

        # obj is a string
        obj, tsmear = obj.split(":")
        obj.strip()

        if obj == "None" or obj is None:
            return cls.NoSmearing()
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
    def NoSmearing():
        return Smearing(1, None)

    @staticmethod
    def FermiDirac(tsmear):
        return Smearing(3, tsmear)

    @staticmethod
    def Gaussian(tsmear):
        return Smearing(7, tsmear)

    @property
    def to_dict(self):
        """json friendly dict representation of Kpoints"""
        raise NotImplementedError("")
        d = {"occopt": self.occopt, 
             "tsmear": self.tsmear,
            }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
                                                                     
    @staticmethod
    def from_dict(d):
        return Smearing(d["occopt"], d["tsmear"])

##########################################################################################

class Electrons(AbinitCard):
    """
    Electron section.
    """
    #: Different configurations of the electron density.
    SM = collections.namedtuple('SpinMode', "nsppol nspinor nspden")

    _mode2vars = {
        "unpolarized"       : SM(1, 1, 1),
        "polarized"         : SM(2, 1, 2),
        "afm"               : SM(1, 1, 2),
        "spinor"            : SM(1, 2, 4),
        "spinor_nomagnetic" : SM(1, 2, 1),
    }
    del SM

    modes = Enum(k for k in _mode2vars)

    VALID_VARS = []

    def __init__(self, 
        spin_mode = modes.polarized, 
        nband     = None,
        fband     = None,
        smearing  = None,
        iscf      = None,
        diemac    = None,
        charge    = None,
        #occupancies = None,
        comment   = None,
        ):
        """
        Highly flexible constructor for Electrons object.  The flexibility comes
        at the cost of usability and in general, it is recommended that you use
        the default constructor only if you know exactly what you are doing and
        requires the flexibility.  For most usage cases, the object can be constructed 
        far more easily using the convenience static constructors: 
        and it is recommended that you use those.
                                                                                
        Args:
            comment:
                String comment for Kpoints

            charge: float
                    Total charge of the system. Default is 0.
        """
        if spin_mode not in Electrons.modes:
            raise ValueError("Unknow spin_mode %s" % spin_mode)

        super(Electrons, self).__init__()

        self._comment = comment
        self.spin_mode = spin_mode

        spin_vars = Electrons._mode2vars[spin_mode]

        # FIXME
        if nband is None:
            fband = 4

        self.update({
            "nsppol" : spin_vars.nsppol,
            "nspinor": spin_vars.nspinor,
            "nspden" : spin_vars.nspden,
            "nband"  : nband,
            "fband"  : fband,
            "iscf"   : iscf,
            "diemac" : diemac,
            #diemac  : None if self..ismetal else ... TODO
        })

        self.smearing = smearing

        if smearing is not None:
            self.update({
                "occopt": smearing.occopt,
                "tsmear": smearing.tsmear,
            })

    @property 
    def nsppol(self): return self["nsppol"]

    @property 
    def nspinor(self): return self["nspinor"]

    @property 
    def nspden(self): return self["nspden"]

    @property 
    def nband(self): return self["nband"]

    @property
    def to_dict(self):
        """json friendly dict representation of Kpoints"""
        raise NotImplementedError("")

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")
        return Electrons(
            spin_mode = modes.polarized, 
            nband     = None,
            fband     = None,
            smearing  = None,
            iscf      = None,
            diemac    = None,
            #occupancies = None,
            comment   = None,
            )


##########################################################################################

class RelaxStrategy(object):
    "Container to store the variables for (constrained) optimization"

    def __init__(self, ionmov=3, optcell=2, **vars):
        self.ionmov = ionmov
        self.optcell = optcell

        self.extra_vars = vars
        #        ionmov      = 3, 
        #        optcell     = 2, 
        #        dilatmx     = 1.1, 
        #        ecutsm      = 0.5, 
        #        ntime       = 80, 
        #        strfact     = 100,
        #        tolmxf      = None,
        #        strtarget   = None,
        #        ):

    @classmethod
    def asstrategy(cls, obj):
        """
        Constructs an instance of RelaxStrategy from obj. 
        Accepts obj in the form:

            * RelaxStrategy object.
            * String, e.g. "ions", "cell", "ions-cell"
        """
        if isinstance(obj, cls):
            return obj

        # obj is a String
        if "-" in obj:
            items = obj.split(":")
            assert len(items) == 2
        else:
            items = [obj,]

        ionmov = 3 if "ions" in items else 0
        optcell = 2 if "cell" in items else 0
        return cls(ionmov, optcell)

    @property
    def move_ions(self):
        "True if ions must be moved"
        return self.ionmov != 0

    @property
    def move_cell(self):
        "True if lattice parameters must be optimized"
        return self.optcell != 0

    @property
    def with_constraints(self):
        "Dictionary with the constraints. Empty dict if no constraint is enforced"
        return {}

class GeoOptimization(AbinitCard):
    "Card containing the parameters for the relaxation of the atomic position and of the unit cell."

    def __init__(self, strategy,
                ionmov      = 3, 
                optcell     = 2, 
                dilatmx     = 1.1, 
                ecutsm      = 0.5, 
                ntime       = 80, 
                tolmxf      = None,
                strtarget   = None,
                strfact     = 100,
                ):

        super(GeoOptimization, self).__init__()

        self.update({ 
            "ionmov"    : ionmov,
            "optcell"   : optcell,
            "ntime"     : ntime,
            "tolmxf"    : tolmxf,
            "strfact"   : 100 if optcell != 0 else None,
            "ecutsm"    : ecutsm  if optcell != 0 else None,
            "dilatmx"   : dilatmx if optcell != 0 else None,
            "strtarget" : strtarget,
            #"tolrff"   : 0.02,
            #"toldff"   : 1.0D-7 (or tolvrs 1.0D-12 if all ions are on special positions)
        })

        pprint(self)

    @property
    def to_dict(self):
        """json friendly dict representation of Kpoints"""
        raise NotImplementedError("")

    @staticmethod
    def from_dict(d):
        return GeoOptimization(
                iomov      = 3, 
                optcell    = 2, 
                dilatmx    = 1.1, 
                ecutsm     = 0.5, 
                ntime      = 80, 
                tolmxf     = None,
                strtarget  = None,
                strfact    = 100,
                constraint = None,
                )

##########################################################################################
#class Phonons(AbinitCard):
#    "Card containing the parameters for the computation of phonons."
#
#    _optdriver = 1
#
#    def __init__(self,  qpoint):
#        raise NotImplementedError("")
#        super(Screening, self).__init__()

#class DDK(AbinitCard):

##########################################################################################

class PPModel(MSONable):
    "Pameters defining the plasmon-pole technique."

    _mode2ppmodel = {  
        'None'      : 0,
        'hybersten' : 1,
        'godby'     : 2,
        'linden'    : 3,
        'farid'     : 4,
    }

    modes = Enum(k for k in _mode2ppmodel)

    @classmethod
    def asppmodel(cls, obj):
        """
        Constructs an instance of PPModel from obj. 
        Accepts obj in the form:
            * PPmodel instance
            * None
            * string. e.g "godby:12.3 eV", "linden".
        """
        if isinstance(obj, cls):
            return obj
                                                           
        if obj is None:
            return Smearing.NoSmearing()

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

    def __init__(self, mode=modes.godby, plasmon_freq=None):
        assert mode in PPModel.modes
        self.mode = mode
        self.plasmon_freq = plasmon_freq

    def tovariables(self):
        return {"ppmodel" : self._mode2ppmodel[self.mode],
                "ppmfrq"  : self.plasmon_freq}

    def __repr__(self):
        return "<%s at %s, mode = %s>" % (self.__class__.__name__, id(self), str(self.mode))

    def __eq__(self, other):
        if other is None: return False 
        return (self.mode == other.mode and
                self.plasmon_freq == other.plasmon_freq)

    def __ne__(self, other):
        return not self == other

    @staticmethod
    def Hybertsen():
        return PPModel(mode="hybersten")

    @staticmethod
    def Godby(plasmon_freq=None):
        return PPModel(mode="godby", plasmon_freq=plasmon_freq)

    @staticmethod
    def Linden():
        return PPModel(mode="linden")

    @staticmethod
    def Engel():
        return PPModel(mode="engel")

    @property
    def to_dict(self):
        d = {"mode" : self.mode,
             "plasmon_freq" : plasmon_freq}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return PPmodel(mode=d["mode"], plasmon_freq=d["plasmon_freq"])

##########################################################################################

class ScreeningFrequencyMesh(MSONable):

    #"freqremax"
    #"freqremin"
    #"nfreqre"
    #"nfreqim",

    def __init__(self, nomega_real=None, maxomega_real=None, nomega_imag=None):
        self.nomega_real = nomega_real
        self.maxomega_real = maxomega_real
        self.nomega_imag = nomega_imag

    def tovariables(self):
        raise NotImplementedError()
        return {"nomega" : self._mode2ppmodel[self.mode],
                "ppfreq"  : self.plasmon_freq}

    @property
    def to_dict(self):
        d = {"nomega_real"   : self.nomega_real,
             "maxomega_real" : self.maxomega_real,
             "nomega_imag"   : self.nomega_imag,
             }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return ScreeningFrequencyMesh(
            nomega_real   = d["nomega_real"],
            maxomega_real = d["maxomega_real"],
            nomega_imag   = d["nomega_imag"], 
            )

##########################################################################################

class Screening(AbinitCard):
    "Card containing the parameters used for the computation of the screening function."

    _modes = [
        "automatic",
        "adler_wiser",
        "hilbert_transform",
    ]

    modes = Enum(k for k in _modes)

    _optdriver = 3
                                              
    def __init__(self, ecuteps, ppmodel_or_freqmesh, 
                 mode    = "automatic",
                 symchi  = 1,
                 inclvkb = 2,
                 ecutwfn = None,
                 gwmem   = None,
                 fftgw   = None,
                 comment = None,
                ):

        assert mode in self.modes

        super(Screening, self).__init__()

        self._comment = comment

        gwcalctyp, spmeth, nomegasf = 3 * (None,)

        self.ppmodel_or_freqmesh = ppmodel_or_freqmesh

        #self.update(ppmodel_or_freqmesh.tovariables())

        self.update({
            "ecuteps"    : ecuteps,
            "gwcalctyp"  : gwcalctyp,
            "ecutwfn"    : ecutwfn,
            "awtr"       : 1,
            "symchi"     : symchi,
            "spmeth"     : None,
            "nomegasf"   : None,
            "gwmem"      : gwmem,
            "fftgw"      : fftgw,
            "inclvkb"    : inclvkb,
        })

    @staticmethod
    def static():
        raise NotImplementedError("")

    @staticmethod
    def for_ppmodel(ppmodel):
        raise NotImplementedError("")

    @staticmethod
    def for_contour_deformation(frequency_mesh):
        raise NotImplementedError("")

    @property
    def to_dict(self):
        d = {"ecuteps"   : self.ecuteps,
             }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return Screening(ecuteps, ppmodel_or_freqmesh, 
                 mode    = "automatic",
                 symchi  = 1,
                 inclvkb = 2,
                 ecutwfn = None,
                 gwmem   = None,
                 fftgw   = None,
                 comment = None,
                )

##########################################################################################

class SelfEnergy(AbinitCard):
    "Parameters for the computation of the self-energy."

    _types = {
        "gw"           : 0,
        "hartree_fock" : 5,
        "sex"          : 6,
        "cohsex"       : 7,
        "model_gw_ppm" : 8,
        "model_gw_cd"  : 9,
    }

    types = Enum(_types)

    _scmodes = {
        "one_shot"      : 0,
        "energy_only"   : 1,
        "wavefunctions" : 2,
    }

    scmodes = Enum(_scmodes.keys())

    _optdriver = 4

    def __init__(self, nband_sigma, ecuteps, ecutsigx, kptgw, bdgw, 
                 type                = types.gw,
                 scmode              = scmodes.one_shot,
                 ppmodel_of_freqmesh = None,
                 gwmem               = None,
                 fftgw               = None,
                 comment             = None
                ):

        super(SelfEnergy, self).__init__()

        self._comment = comment
                               
        assert type in SelfEnergy.types
        assert scmode in SelfEnergy.scmodes 
        self.type = type
        self.scmode = scmode

        gwcalctyp = None
                                                                  
        #self.update(type_scmode_tovariables())

        # TODO: check nsppol = 2 case
        kptgw = np.reshape(kptgw, (-1,3))
        nkptgw = len(kptgw)

        bdgw  = np.reshape(bdgw, (-1,2))
        assert len(bdgw == nkptgw)

        symsigma = 0 if scmode == SelfEnergy.scmodes.wavefunctions else 1

        self.update({
            "ecuteps"  : ecuteps,
            "ecutsigx" : ecutsigx,
            "nkptgw"   : nkptgw,
            "kptgw"    : kptgw,
            "bdgw"     : bdgw,
            "gwmem"    : gwmem,
            "fftgw"    : fftgw,
            "symsigma" : symsigma
        })

    def type_scmode_tovars(self):
        d = { "gwcalctyp"  : gwcalctyp,
            }
        return d

##########################################################################################

class BetheSalpeter(AbinitCard):
    "Card containing the parameters for the solution of the Bethe-Salpeter equation."

    _optdriver = 99

    #types = Enum(_types)

    #_modes = {
    #    "Tamm_Dandcoff" : 0,
    #    "Coupling"      : 1,
    #}

    #scmodes = Enum(_scmodes.keys())

    _algo2var = {
        "direct_diago" : 1, 
        "haydock"      : 2,
        "cg"           : 3,
    }

    #algorithms = Enum(_algo2var)

    _coulomb_mode = {
        "diago",
        "full",
        "model_df"
    }

    def __init__(self, loband, nband,
                 with_lf             = True,
                 gwmem               = None,
                 fftgw               = None,
                 comment             = None
                ):

        raise NotImplementedError("")

        super(SelfEnergy, self).__init__()

        self._comment = comment
                               
        #self.type = type
        #self.scmode = scmode
        #gwcalctyp = None
                                                                  
        self.update({
            "bs_exchange_term" : 1 if with_lf else 0,
            "inclvkb"  : 2,
            "ecuteps"  : ecuteps,
            "soenergy" : soenergy,
            "gwmem"    : gwmem,
            "fftgw"    : fftgw,
            "bs_freq_mesh" : freq_mesh
            #bs_haydock_niter5 60      # No. of iterations for Haydock
            #bs_hayd_term5      0      # No terminator
            #bs_haydock_tol5 0.05 0
            #mdf_epsinf         12.0
        })

##########################################################################################

class Control(AbinitCard):

    #: Basic variables needed for the different runlevels
    _mode2optdriver = {
        "scf"       : 0 , 
        "nscf"      : 0 ,
        "dfpt"      : 1 ,
        "screening" : 3 ,
        "sigma"     : 4 ,
        "bse"       : 99,
    }

    modes = Enum(k for k in _mode2optdriver)

    #: Tolerance used by the runlevels.
    _mode2tolname = {
        "scf"       : 'tolvrs', 
        "nscf"      : 'tolwfr', 
        "dfpt"      : 'toldfe',   # ?
        "screening" : 'toldfe',   # dummy
        "sigma"     : 'toldfe',   # dummy
        "bse"       : 'toldfe',   # ?
    }

    # Tolerances for the different levels of accuracy.
    T = collections.namedtuple('Tolerance', "low normal high")
    _tolerances = {
        "toldfe" : T(1.e-9,  1.e-10, 1.e-11), 
        "tolvrs" : T(1.e-7,  1.e-8,  1.e-9),
        "tolwfr" : T(1.e-15, 1.e-17, 1.e-19),
        "tolrdf" : T(0.04,   0.02,   0.01),
        }
    del T

    accuracies = Enum(("low", "normal", "high",))

    #paral_modes = Enum("automatic", "explicit",)

    _mandatory_cards = ['System', 'Kpoints',  'Electrons']

    def __init__(self, *cards, **kwargs):

        super(Control, self).__init__()

        # Keep a reference to the input kwargs so that we can call the constructor via self.
        self._cards = cards
        self._input_kwargs = kwargs

        # Get optdriver from the input cards by inspecting the class attribute _optdriver.
        optdrivers, names_found = [], []
        for card in cards:
            names_found.append(card.name)
            try:
                optdrivers.append(card._optdriver)
            except AttributeError:
                pass

        # Sanity check: no card type duplication is allowed.
        counter = collections.Counter(names_found)
        duplicated = filter(lambda t : t[1] > 1, counter.items())

        if duplicated:
            raise ValueError("Found duplicated card names: %s" % str(duplicated))

        # Check for presence of mandatory cards.
        if not set(names_found).issuperset(Control._mandatory_cards):
            raise ValueError("Mandatory cards are missing, got %s" % str(names_found))

        # Extract useful cards.
        system    = [card for card in cards if card.name == "System"][0]
        electrons = [card for card in cards if card.name == "Electrons"][0]

        # Default is Ground-State calculation (SCF or NSCF run).
        self.optdriver = 0

        if optdrivers:
            self.optdriver = optdrivers[0]
            if len(optdrivers) > 1:
                err_msg = "Got more than one value for optdriver! %s" % optdrivers
                raise ValueError(err_msg)

        # Find mode from optdriver and iscf.
        match = [k for (k,v) in self._mode2optdriver.items() if v == self.optdriver] 
        if not match:
            raise ValueError("Cannot find mode string for optdriver %s" % self.optdriver)
                                                                                          
        if len(match) == 1:
            self.mode = match[0]
        else:
            # It's a GS run. Look at the value of iscf.
            assert self.optdriver == 0
            if electrons.iscf in [-2 ,-3]:
                self.mode = "nscf"
            else:
                self.mode = "scf"

        # Setup accuracy, cutoff and other basic parameters.
        self.accuracy = kwargs.pop("accuracy",    "normal")
        ecut          = kwargs.pop("ecut",        None)
        pawecutdg     = kwargs.pop("pawecutdg",   None)
        prtwf         = kwargs.pop("prtwf",       None)
        prtden        = kwargs.pop("prtden",       None)
        boxcutmin     = kwargs.pop("boxcutmin",   None)
        want_forces   = kwargs.pop("want_forces", False)
        want_stress   = kwargs.pop("want_stress", False)

        self._comment = "mode = %s, accuracy = %s" % (self.mode, self.accuracy)

        hints = [p.hint_for_accuracy(self.accuracy) for p in system.pseudos]

        if system.isnc:
            if ecut is None:
                # Find the cutoff energy from the hints.
                ecut = max(hint.ecut for hint in hints)

        elif system.ispaw:
            raise NotImplementedError("")
            #ratio = max(p.suggested_augratio(accuracy) for p in self.pseudos])
            #ratio = augration_high if high else augratio_norm 
            #pawecutdg = ecut * ratio

        else:
            raise RuntimeError("Neither NC nor PAW!!")

        # Default tolerances
        tol_varname = self._mode2tolname[self.mode]
        tol = self._tolerances[tol_varname]

        #min_nband = self.System.num_valence_electrons // 2

        nbdbuf, nband = None, electrons.nband

        if tol_varname == "tolwfr" and nband is not None:
            nbdbuf = max(4, nband * 0.02) 

        prtgkk = 1 if self.mode == "dfpt" else None

        aidx = ["low", "normal", "high"].index(self.accuracy)

        self.update( {
            "optdriver" : self.optdriver,
            "ecut"      : ecut,
            "nbdbuf"    : nbdbuf,
            "nstep"     : [50,  75, 100][aidx],
            tol_varname : getattr(tol, self.accuracy),
            "prtwf"     : prtwf,
            "prtden"    : prtden,
            "prtgkk"    : prtgkk,
            "boxcutmin" : boxcutmin,
            "optforces" : None if (self.need_forces or want_forces) else 2,
            "optstress" : None if (self.need_stress or want_stress) else 0,
        })

        if system.ispaw:
            self.update({
                "pawecutdg" : pawecutdg,
            })

        # Add what is left in kwargs.
        self.update(**kwargs)

    @property
    def need_forces(self):
        "True if forces are required at each SCF step (like the stresses)."
        return self.has_card_names("GeoOptimization")
                                                                            
    @property
    def need_stress(self):
        "True if the computation of the stress is required"
        # TODO: here it's easier to check if optcell != 0
        return self.has_card_names("GeoOptimization")

    def has_card_names(self, *card_names):
        self_card_names = [card.name for card in self._cards]
        for card_name in card_names: 
            if card_name not in self_card_names: return False
        return True

    @property
    def to_dict(self):
        raise NotimplementedError("")
        #d = {k: v.to_dict for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        raise NotimplementedError("")
        return Control(*cards, **kwargs)

##########################################################################################

class Input(dict, MSONable):
    """
    Class to contain the set of ABINIT input variable defining the calculation.
    Essentially an ordered dictionary of predefined AbinitCard
    """

    #: Cards that must be passed to the constructor.
    _mandatory_cards = ['System', 'Kpoints',  'Electrons', 'Control',]

    def __init__(self, *cards):
        """
        Args:
            cards:
                List of AbinitCards. Mandatory cards: system, kpoints, electrons
        """
        super(Input, self).__init__()

        for card in cards:
            if card.name in self:
                raise ValueError("Card %s is already in the dictionary" % card.name)
            self[card.name] = card 

        # Sanity check: no variable duplication is allowed.
        counter = collections.Counter(self.varnames)
        duplicated = filter(lambda t : t[1] > 1, counter.items())
                                                                                  
        if duplicated:
            raise ValueError("Found duplicated variables: %s" % str(duplicated))

        # Check for presence of mandatory cards.
        card_names = set([card.name for card in cards])
        if not set(card_names).issuperset(Input._mandatory_cards):
            raise ValueError("Mandatory cards are missing, got %s" % str(card_names))

    def __getattribute__(self, name):
        try:
            # Default behaviour
            return super(Input, self).__getattribute__(name)
        except AttributeError:
            # Try in self
            return self[name]

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), str([k for k in self.card_names]))

    def __str__(self):
        output = [
            "# This file has been generated by pymatgen",
            "# *** Do not edit. Any change will be lost if the input if regenerated ***"
            ]
        keys = self.card_names
        keys.sort()
        for k in keys:
            card = self[k]
            output.append( k.center(len(k)+2).center(90, "#") )
            output.append(str(card))
            output.append("")
        return "\n".join(output)

    def __eq__(self, other):
        if other is None: return False
        if set(self.card_names) != set(other.card_names):
            return False
        for (card_name, card) in self:
            if card != other[card_name]: return False
        return True

    def __ne__(self, other): 
        return not self == other

    @property
    def pseudos(self):
        "List of pseudopotentials"
        return self.System.pseudos

    @property
    def structure(self):
        "Pymatgen structure"
        return self.System.structure

    @property
    def description(self):
        "Description string"
        lines = []
        app = lines.append
        for (card_name, card) in self.items():
            card_comment = getattr(card, "comment", "No comment avaiable")
            app(card_name + ": " + str(card_comment))
        return "\n".join(lines)

    @property
    def to_dict(self):
        d = {"cards" : self.cards}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return Input(*d["cards"])

    @property
    def card_names(self):
        return self.keys()

    @property
    def cards(self):
        return self.values()

    def copy(self):
        return Input(*self.cards)

    def iter_variables(self):
        "Iterator over (varname, varvalue)"
        for card in self.cards:
            for item in card.items(): 
                if item[0] != "_attributes": yield item

    @property
    def varnames(self):
        return [v[0] for v in self.variables]

    @property
    def variables(self):
        return list(self.iter_variables())

    @property
    def varnames2card(self):
        d = {}
        for (card_name, card) in self.items():
            d[card_name] = card.keys()

        vnames2card = {}
        for (card_name, vname_list) in d.items():
            for var_name in vname_list:

                if var_name == "_attributes": 
                    continue

                if var_name in vnames2card:
                    raise ValueError("variable %s occurs more than once" % var_name)

                vnames2card[var_name] = self[card_name]

        return vnames2card

    def get_variable(self, variable_name):
        try:
            card = self.varnames2card[variable_name]
            return card[variable_name]
        except KeyError:
            raise KeyError("variable %s not found" % variable_name)

    def set_variable(self, variable_name, variable_value):
        try:
            card = self.varnames2card[variable_name]
            card[variable_name] = variable_value
        except KeyError:
            raise KeyError("variable %s not found" % variable_name)

    def delete_variable(self, variable_name):
        try:
            card = self.varnames2card[variable_name]
            del card[variable_name]
        except KeyError:
            raise KeyError("variable %s not found" % variable_name)

    def copy_cards(self, exclude=None, aslist=True):
        if exclude is None: 
            exclude = []
        else:
            # Check for typos!
            if not set(exclude).issubset(self.card_names):
                print("exclude_list: %s" % str(exclude))
                print("card_names: %s" % str(self.card_names))
                raise ValueError("Wrong value for exclude, likely a typo")
            
        # The entry "Control" is automatically generated when we instanciate 
        # Control hence we don't add it to the cards.
        new_dict = {}
        for card_name in self.card_names:
            if card_name not in exclude + ["Control"]: 
                new_dict[card_name] = self[card_name].copy()

        if aslist:
            return new_dict.values()
        else:
            return new_dict

    # These methods are obsolete!
    def add_card(self, card):
        if card.name in self.card_names:
            raise ValueError("Card %s is already in %s" % (card.name, repr(self)))
        new = Input(*self.cards)
        new[card.name] = card
        return new

    def replace_card(self, card_name, new_card, check_name=True):
        "Modify the value of self[card_name] with new_card"

        old_card = self[card_name]
        if check_name and old_card.name != new_card.name:
            raise ValueError("Cannot replace %s with Card of type %s" % (card_name, new_card.name))

        new = Input(*self.cards)
        new[card_name] = new_card
        return new

    def replace_control_with(self, new_card):   
        return self.replace_card("Control", new_card)

    def replace_kpoints_with(self, new_card):   
        return self.replace_card("Kpoints", new_card)

    def add_vars_to_card(self, vars, card_name):
        "Add new_vars to the AbinitCard self[card_name]. Return new instance of Input."

        if card_name not in self.card_names:
            raise ValueError("card_name %s not present in Input" % card_name)

        cards = {}
        for (k, v) in self.items():
            if k == card_name:
                cards[k] = self[k].copy()
                cards[k].update(vars)
            else:
                cards[k] = self[k]

        return Input(*cards.values())

    def add_vars_to_control(self, vars):
        return self.add_vars_to_card(vars, "Control")

    def add_vars_to_electrons(self, vars):
        return self.add_vars_to_card(vars, "Electrons")

    def add_vars_to_kpoints(self, vars):
        return self.add_vars_to_card(vars, "Kpoints")

    def add_vars_to_system(self, vars):
        return self.add_vars_to_card(vars, "System")

    @staticmethod
    def SCF_groundstate(structure, pseudos, 
                        ngkpt = None,
                        kppa  = None,
                        spin_mode = "polarized", 
                        smearing  = None,
                        **kwargs
                        ):
        """
        Constructor for self-consistent ground-state calculations.

        Args:
            structure:
                pymatgen structure
            pseudos:
                List of pseudopotentials.
            ngkpt:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            kppa:
                Grid density (kppa / natom). 
            spin_mode: 
                Flag defining the spin polarization (nsppol, nspden, nspinor). Defaults to "polarized"
            smearing: 
                Smearing instance. None if smearing technique is not used. 
            **kwargs:
                Extra variables added directly to the input file

        Returns:
            AbinitInput instance.
        """

        # Initialize system section from structure.
        system = System(structure, pseudos)

        # Variabled for electrons
        electrons = Electrons(spin_mode=spin_mode, smearing=smearing)

        # K-point sampling.
        if [obj is not None for obj in [ngkpt, kppa,]].count(True) != 1:
            raise ValueError("only one variable among ngkpt, kppa should be supplied")

        if ngkpt is not None:
            kmesh = Kpoints.monkhorst_automatic(structure, ngkpt)
        elif kppa is not None:
            kmesh = Kpoints.automatic_density(structure, kppa)

        scf_control = Control(system, electrons, kmesh, **kwargs)

        return Input(system, electrons, kmesh, scf_control)

    @staticmethod
    def NSCF_kpath_from_SCF(scf_input, nband, 
                            ndivsm       = 20, 
                            kpath_bounds = None, 
                            **kwargs
                           ):
        """
        Constructor for non-self-consistent ground-state calculations (band structure calculations)
                                                                                                       
        Args:
            scf_input:
                Input for self-consistent calculations.
            nband:
                Number of bands to compute
            ndivsm:
                Number of division for the smallest segment
            kpath_bounds:
                Extrema of the k-path.
            **kwargs:
                Extra variables added directly to the input file
                                                                                                       
        Returns:
            AbinitInput instance.
        """
        # Change Kpoints and Electrons.
        nscf_cards = scf_input.copy_cards(exclude=["Kpoints", "Electrons",])

        if kpath_bounds is not None:
            nscf_cards.append(Kpoints.explicit_path(ndivsm, kpath_bounds))
        else:
            nscf_cards.append(Kpoints.path_from_structure(ndivsm, scf_input.structure))

        spin_mode = scf_input.Electrons.spin_mode

        nscf_cards.append(Electrons(spin_mode=spin_mode, nband=nband, iscf=-2))

        nscf_cards.append(Control(*nscf_cards, **kwargs))

        return Input(*nscf_cards)

    @staticmethod
    def NSCF_kmesh_from_SCF(scf_input, nband, 
                            ngkpt = None, 
                            kppa  = None,
                            **kwargs
                           ):
        """
        Constructor for non-self-consistent ground-state calculations (DOS calculations).
                                                                                                       
        Args:
            scf_input:
                Input for self-consistent calculations.
            nband:
                Number of bands to compute
            ngkpt:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            kppa:
                Grid density (kppa / natom). 
            **kwargs:
                Extra variables added directly to the input file
                                                                                                       
        Returns:
            AbinitInput instance.
        """
        # Change Kpoints and Electrons.
        nscf_cards = scf_input.copy_cards(exclude=["Kpoints", "Electrons",])

        # K-point sampling.
        if ngkpt is not None:
            kmesh = Kpoints.monkhorst_automatic(scf_input.structure, ngkpt)
        elif kppa is not None:
            kmesh = Kpoints.automatic_density(scf_input.structure, kppa)

        nscf_cards.append(kmesh)

        spin_mode = scf_input.Electrons.spin_mode

        nscf_cards.append(Electrons(spin_mode=spin_mode, nband=nband, iscf=-3))

        nscf_cards.append(Control(*nscf_cards, **kwargs))

        return Input(*nscf_cards)

    @staticmethod
    def GeometryRelax(structure, pseudos, strategy,
              ngkpt     = None,
              kppa      = None,
              spin_mode = "polarized", 
              smearing  = None,
              **kwargs
             ):
        """
        Constructor for structure relaxation..
                                                                                                       
        Args:
            structure:
                pymatgen structure
            pseudos:
                List of pseudopotentials.
            strategy:
                RelaxStrategy instance
            ngkpt:
                Subdivisions N_1, N_2 and N_3 along reciprocal lattice vectors.
            kppa:
                Grid density (kppa / natom). 
            spin_mode: 
                Flag defining the spin polarization (nsppol, nspden, nspinor). Defaults to "polarized"
            smearing: 
                Smearing instance. None if smearing technique is not used. 
            **kwargs:
                Extra variables added directly to the input file

        Returns:
            AbinitInput instance.
        """
        # Initialize geometry section from structure.
        system = System(structure, pseudos)

        # Variables for electrons
        electrons = Electrons(spin_mode=spin_mode, smearing=smearing)

        # K-point sampling.
        if ngkpt is not None:
            kmesh = Kpoints.monkhorst_automatic(structure, ngkpt)
        elif kppa is not None:
            kmesh = Kpoints.automatic_density(structure, kppa)

        geo_optim = GeoOptimization(strategy)

        control = Control(system, electrons, kmesh, geo_optim, **kwargs)

        return Input(system, electrons, kmesh, geo_optim, control)

    @staticmethod
    def SCR_from_NSCF(nscf_input, ecuteps, ppmodel_or_freqmesh, nband_screening, 
                      smearing = None, 
                      inclvkb  = None,
                      **kwargs
                     ):
        """
        Constructor for screening calculations.
                                                                                                       
        Args:
            nscf_input:
                Input used for the non-self consistent calculation
            ecuteps:
                Cutoff energy [Ha] for the dielectric matrix
            ppmodel_or_freqmesh:
                object describing the frequency dependende of the screening function.
                Either PPmodel instance or ScreeningFrequencyMesh instance.
            nband_screening:
                Number of bands for the computation of the screening. 
            smearing: 
                Smearing instance. None if smearing technique is not used. 
            **kwargs:
                Extra variables added directly to the input file

        Returns:
            AbinitInput instance.
        """
        scr_cards = nscf_input.copy_cards(exclude=["Electrons",])

        spin_mode = nscf_input.Electrons.spin_mode

        scr_cards.append(Electrons(spin_mode=spin_mode, nband=nband_screening, smearing=smearing))

        scr_cards.append(Screening(ecuteps, ppmodel_or_freqmesh, inclvkb=inclvkb))

        scr_cards.append(Control(*scr_cards, **kwargs))

        return Input(*scr_cards)

    @staticmethod
    def SIGMA_from_SCR(scr_input, nband_sigma, ecuteps, ecutsigx, kptgw, bdgw, 
                       type     = "gw",
                       scmode   = "one_shot",
                       smearing = None,
                       **kwargs
                       ):
        """
        Constructor for sigma calculations.
                                                                                                       
        Args:
            scr_input:
                Input used for the screening calculation
            nband_sigma:
                Number of bands for the self-energy
            ecuteps:
                Cutoff energy [Ha] for the dielectric matrix
            ecutsigx:
                Cutoff energy [Ha] for the exchange part of the self-energy
            type:
            scmode:
            smearing: 
                Smearing instance. None if smearing technique is not used. 
            **kwargs:
                Extra variables added directly to the input file

        Returns:
            AbinitInput instance.
        """

        ppmodel_or_freqmesh = scr_input.Screening.ppmodel_or_freqmesh

        se_card = SelfEnergy(nband_sigma, ecuteps, ecutsigx, kptgw, bdgw, 
                             type                = type,
                             scmode              = scmode,
                             ppmodel_of_freqmesh = ppmodel_or_freqmesh,
                             )

        spin_mode = scr_input.Electrons.spin_mode
        assert smearing == scr_input.Electrons.smearing

        sigma_cards = scr_input.copy_cards(exclude=["Screening", "Electrons",])

        sigma_cards.append(Electrons(spin_mode=spin_mode, nband=nband_sigma, smearing=smearing))

        sigma_cards.append(se_card)

        sigma_cards.append(Control(*sigma_cards, **kwargs))

        return Input(*sigma_cards)
