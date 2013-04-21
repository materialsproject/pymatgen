"""
Classes for manipulating/writing ABINIT input files. 
"""
from __future__ import division, print_function

import sys
import os
import os.path
import warnings
import collections
import abc
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.design_patterns import Enum, AttrDict
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr
from pymatgen.core.units import any2Ha
from pymatgen.util.string_utils import str_aligned, str_delimited
from pymatgen.util.decorators import singleton
from pymatgen.serializers.json_coders import MSONable #, PMGJSONDecoder
from pymatgen.symmetry.finder import SymmetryFinder

from .pseudos import PseudoTable
from .abiobjects import AbivarAble, SpinMode, Smearing

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

##########################################################################################

class NewCard(AbivarAble):
    """
    Base class representing a set of Abinit input parameters associated to a particular topic. 
    Essentially consists of a dictionary with some helper functions.

    The objects are stored in the dictionary as {obj.__class__.__name__: obj}

    str(self) returns a string with the corresponding section of the abinit input file.
    """
    def __init__(self, *args, **kwargs):
        "Creates a Card object."
        self.obj_dict = dict()

        self._comment = kwargs.get("comment", "")

    def __str__(self):
        "String representation (the section of the abinit input file)"
        return self.get_string()
                                                                                          
    @property
    def name(self):
        return self.__class__.__name__

    @property
    def objects(self):
        "List of objects stored in self"
        return self.obj_dict.values()

    def get_obj(self, classname):
        "Return the object from its class name"
        return self.obj_dict[classname]

    def add_obj(self, obj):
        "Add the object to the card"
        cname = obj.__class__.__name__
        if cname in self.obj_dict:
            raise ValueError("%s is already in the Card" % cname)
        self.obj_dict[cname] = obj

    def pop_obj(self, obj, *d):
        "Remove the object obj and returns it"
        return self.obj_dict.pop(obj.__class__.__name__, d)

    def to_abivars():
        "Returns a dictionary with the abinit variables defined by the Card"
        abivars = {}
        for obj in self.objects:
            abivars.update(obj.to_abivars())
        return abivars

    @property
    def comment(self):
        "String comment (comment of self + comments of the objects, if any)"
        lines = [self._comment,]

        for obj in self.objects:
            if hasattr(obj, "comment"):
                lines.append("%s: %s" % (obj.__class__.__name__, obj.comment))

        return "\n".join(lines)

    def get_string(self, pretty=False):
        """
        Returns a string representation of self. The reason why this
        method is different from the __str__ method is to provide options for pretty printing.

        Args:
            pretty:
                Set to True for pretty aligned output. Defaults to False.
        """
        abivars = self.to_abivars()
        sorted_keys = sorted([k for k in abivars])

        lines = []
        if self.comment: 
            lines.append(["# comment: ", self.comment.replace("\n", "\n#")])
            lines.append(["#", "\n"])

        for k in sorted_keys:
            value = abivars[k]

            # Use default values if value is None.
            if value is None: 
                continue

            if isinstance(value, collections.Iterable) and not isinstance(value, str):
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

    #def merge(self, other):
    #    """
    #    Add all the values of another Card object to this object.
    #    """
    #    if self.name != other.name:
    #        raise ValueError("Cannot merge cards of different type: %s, %s!" % (self.name, other.name) )

    #    s_items = {k : v for (k, v) in self.items()}
    #    s_attributes = s_items.pop("_attributes")

    #    o_items = {k : v for (k, v) in other.items()}
    #    o_attributes = o_items.pop("_attributes")

    #    new_dict = s_items.copy()
    #    for (k, v) in o_items:
    #        if k in new_dict and new_dict[k] != v:
    #            raise ValueError("Cards have conflicting variables!")
    #        else:
    #            new_dict[k] = v

    #    new_attributes = s_attributes.copy()
    #    for (k, v) in o_attributes:
    #        if k in new_attributes and new_attributes[k] != v:
    #            raise ValueError("Cards have conflicting attributes!")
    #        else:
    #            new_attributes[k] = v

    #    new = Card(new_attributes, **new_dict)
    #    # Return same type as self.
    #    new.__class__ = self.__class__
    #    return new

    #def diff(self, other):
    #    """
    #    Diff function for Card.  Compares two objects and indicates which
    #    parameters are the same and which are not. Useful for checking whether
    #    two runs were done using the same parameters.
                                                                                    
    #    Args:
    #        other:
    #            The other Card object to compare to.
                                                                                    
    #    Returns:
    #        Dict of the following format:
    #        {"Same" : parameters_that_are_the_same,
    #        "Different": parameters_that_are_different}
    #        Note that the parameters are return as full dictionaries of values.
    #        E.g. {"natom":3}
    #    """
    #    similar_param, different_param = {}, {}
    #    for (k1, v1) in self.items():
    #        if k1 not in other:
    #            different_param[k1] = {self.name : v1, other.name : "Default"}
    #        elif v1 != other[k1]:
    #            different_param[k1] = {self.name : v1, other.name : other[k1]}
    #        else:
    #            similar_param[k1] = v1
                                                                                    
    #    for (k2, v2) in other.items():
    #        if k2 not in similar_param and k2 not in different_param:
    #            if k2 not in self:
    #                different_param[k2] = {self.name : "Default", other.name : v2}
                                                                                    
    #    return {"Same": similar_param, "Different": different_param}

class Card(dict, MSONable):
    """
    Base class representing a set of Abinit input parameters associated to a particular topic. 
    Essentially consists of a dictionary with some helper functions.

    The input variables are stored in the dictionary as {"varname": value}, whereas instance
    attributes are stored in self[_attributes].
    str(object) returns a string with the corresponding section of the abinit input file.
    """
    def __init__(self, attributes=None, **kwargs):
        """
        Creates a Card object.

        Args:
            attributes:
                dictionary with extra attributes to be stored in self["_attributes"]
            kwargs:
                The input variables as dictionary {varname : varvalue}
        """
        super(Card, self).__init__(**kwargs)

        self["_attributes"] = {}
        if attributes is not None:
            self["_attributes"].update(attributes)

    def __str__(self):
        "String representation (return the section of the abinit input file)"
        return self.get_string()
                                                                                          
    def __getattribute__(self, name):
        try:
            # Default behaviour
            return super(Card, self).__getattribute__(name)
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
            super(Card, self).__setattr__(name, value)
        else:
            self["_attributes"].update({name : value})

    def __add__(self, other):
        return self.merge(other)

    def copy(self):
        "Return a shallow copy of self."
        kwargs =  super(Card, self).copy()
        attributes = kwargs.pop("_attributes", None)
        new = Card(attributes, **kwargs)
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
                Set to True to sort the Card parameters alphabetically. Defaults to False.
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

            value = self[k]

            # Use default values if value is None.
            if value is None: 
                continue

            if isinstance(value, collections.Iterable) and not isinstance(value, str):
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
        Diff function for Card.  Compares two objects and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.

        Args:
            other:
                The other Card object to compare to.

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
        Add all the values of another Card object to this object.
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
                raise ValueError("Cards have conflicting variables!")
            else:
                new_dict[k] = v

        new_attributes = s_attributes.copy()
        for (k, v) in o_attributes:
            if k in new_attributes and new_attributes[k] != v:
                raise ValueError("Cards have conflicting attributes!")
            else:
                new_attributes[k] = v

        new = Card(new_attributes, **new_dict)
        # Return same type as self.
        new.__class__ = self.__class__
        return new

##########################################################################################

class SystemCard(Card):
    """
    Object for representing the geometry of the system and the list of pseudopotentials

    .. attribute:: structure
        Associated Structure.

    .. attribute:: comment
        Optional comment string.
    """
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
            raise ValueError("Structure with partial occupancies cannot be converted into SystemCard!")

        super(SystemCard, self).__init__()

        self.structure = structure
        self.nsym = nsym

        self._comment = structure.formula if comment is None else comment

        # Extract pseudos for this calculation (needed when we pass an entire periodic table)
        table = PseudoTable.astable(pseudos)

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
                                                                     
        natom = structure.num_sites
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
            "xred"   : xred,
            "znucl"  : znucl_type,
            "nsym"   : nsym,
            #"pseudo_list": ", ".join(p.path for p in self.pseudos),
        })

    @property
    def isnc(self):
        "True if norm-conserving calculation"
        return self.pseudos.allnc

    @property
    def ispaw(self):
        "True if PAW calculation"
        return self.pseudos.allpaw

    # TODO
    #@property
    #def num_valence_electrons(self):
    #    for specie in self.structure:
    #        print specie
    #        specie.pseudo
    #    return sum(p.num_valence_electrons for p in self.pseudos)

    @classmethod
    def boxed_molecule(cls, pseudos, cart_coords, acell=3*(10,)):
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
                                                                                 
        return cls(structure, pseudos, comment=comment)

    @classmethod
    def boxed_atom(cls, pseudo, cart_coords=3*(0,), acell=3*(10,)):
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
        return cls.boxed_molecule([pseudo], cart_coords, acell=acell)

    @property
    def to_dict(self):
        return {"@module"   : self.__class__.__module__,
                "@class"    : self.__class__.__name__,
                "structure" : self.structure.to_dict,
                "pseudos"   : self.pseudos.to_dict,
                "nsym"      : self.nsym,
                "comment"   : self._comment,
        }

    @staticmethod
    def from_dict(d):
        return SystemCard(
            d["structure"],
            d["pseudos"],
            d["nsym"],
            comment=d["comment"],
        )

##########################################################################################

class KpointsCard(Card):
    """
    Input variables defining the K-point sampling.
    """
    #: Modes supported by the constructor.
    modes = Enum(('monkhorst', 'path', 'automatic',))

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
        Highly flexible constructor for KpointsCard object.  The flexibility comes
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
                Mode for generating k-poits. Use one of the KpointsCard.modes enum types.
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
        if mode not in KpointsCard.modes:
            raise ValueError("Unknown kpoint mode %s" % mode)

        super(KpointsCard, self).__init__()

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
            KpointsCard object
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
            KpointsCard object
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
            KpointsCard object
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
            KpointsCard object
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

        return cls(mode     = KpointsCard.modes.path, 
                   num_kpts = ndivsm,
                   kpts     = kpath_bounds,
                   comment  = comment if comment else "K-Path scheme", 
                  )

    @classmethod
    def path_from_structure(cls, ndivsm, structure):
        "See _path for the meaning of the variables"
        return cls._path(ndivsm, 
                             structure=structure, 
                             comment="K-path generated automatically from pymatgen structure")
                                                                                                   
                                                                                                   
    @classmethod
    def explicit_path(cls, ndivsm, kpath_bounds):
        "See _path for the meaning of the variables"
        return cls._path(ndivsm, 
                         kpath_bounds=kpath_bounds,
                         comment="Explicit K-path")


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
        """json friendly dict representation of KpointsCard"""
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        raise NotImplementedError("")
        return d

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")

##########################################################################################

class ElectronsCard(Card):
    """
    Electrons section.
    """
    def __init__(self, 
        spin_mode = "polarized", 
        smearing  = "fermi_dirac:0.1 eV",
        nband     = None,
        fband     = None,
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
                String comment for ElectronsCard

            charge: float
                    Total charge of the system. Default is 0.
        """
        super(ElectronsCard, self).__init__()

        self._comment = comment

        self.spin_mode = SpinMode.asspinmode(spin_mode)

        self.smearing = Smearing.assmearing(smearing)

        # FIXME
        if nband is None:
            fband = 4

        self.update(self.spin_mode.to_abivars())

        self.update({
            "nband"  : nband,
            "fband"  : fband,
            "iscf"   : iscf,
            "diemac" : diemac,
            #diemac  : None if self..ismetal else ... TODO
        })

        self.update(self.smearing.to_abivars())

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
    def nband(self): 
        return self["nband"]

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

##########################################################################################

class RelaxCard(Card):
    "Card containing the parameters for the relaxation of the atomic position and of the unit cell."

    def __init__(self, strategy):

        super(RelaxCard, self).__init__()

        self.strategy = strategy

        self.update(strategy.to_abivars())

    @property
    def to_dict(self):
        "json friendly dict representation"
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        raise NotImplementedError("")

    @staticmethod
    def from_dict(d):
        raise NotImplementedError("")

##########################################################################################
#class PhononsCard(Card):
#    "Card containing the parameters for the computation of phonons."
#    _optdriver = 1
#
#    def __init__(self,  qpoint):
#        raise NotImplementedError("")
#        super(PhononsCard, self).__init__()

#class DDKCard(Card):

##########################################################################################

class ScreeningCard(Card):
    "Card containing the parameters used for the computation of the screening function."
    _optdriver = 3
                                              
    def __init__(self, strategy, comment=None):

        super(ScreeningCard, self).__init__()

        self._comment = comment

        self.strategy = strategy 

        self.update(strategy.to_abivars())

    @property
    def to_dict(self):
        d = {"strategy": self.strategy.to_dict,
             "comment" : self._comment,
            }
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        return cls(d["strategy"], comment=d["comment"])

##########################################################################################

class SelfEnergyCard(Card):
    "Parameters for the computation of the self-energy."
    #_types = {
    #    "gw"           : 0,
    #    "hartree_fock" : 5,
    #    "sex"          : 6,
    #    "cohsex"       : 7,
    #    "model_gw_ppm" : 8,
    #    "model_gw_cd"  : 9,
    #}

    #types = Enum(_types)

    #_scmodes = {
    #    "one_shot"      : 0,
    #    "energy_only"   : 1,
    #    "wavefunctions" : 2,
    #}

    #scmodes = Enum(_scmodes.keys())

    _optdriver = 4

    def __init__(self, strategy, comment=None):

        super(SelfEnergyCard, self).__init__()

        self._comment = comment

        self.strategy = strategy
                               
        self.update(strategy.to_abivars())

##########################################################################################

class BSECard(Card):
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

    def __init__(self, strategy, comment=None):

        raise NotImplementedError("")

        super(BSECard, self).__init__()

        self._comment = comment

        self.strategy = strategy
                               
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
            #loband  = None, 
            #nband   = None,
            #with_lf = True,
            #gwmem   = None,
            #fftgw   = None,
            #comment = None
        })

        self.update(strategy.to_abivars())

##########################################################################################

class ControlCard(Card):

    #: Basic variables needed for the different runlevels
    _mode2optdriver = {
        "scf"      : 0 , 
        "nscf"     : 0 ,
        "dfpt"     : 1 ,
        "screening": 3 ,
        "sigma"    : 4 ,
        "bse"      : 99,
    }

    modes = Enum(k for k in _mode2optdriver)

    #: Tolerance used by the runlevels.
    _mode2tolname = {
        "scf"      : 'tolvrs', 
        "nscf"     : 'tolwfr', 
        "dfpt"     : 'toldfe',   # ?
        "screening": 'toldfe',   # dummy
        "sigma"    : 'toldfe',   # dummy
        "bse"      : 'toldfe',   # ?
    }

    # Tolerances for the different levels of accuracy.
    T = collections.namedtuple('Tolerance', "low normal high")
    _tolerances = {
        "toldfe": T(1.e-9,  1.e-10, 1.e-11), 
        "tolvrs": T(1.e-7,  1.e-8,  1.e-9),
        "tolwfr": T(1.e-15, 1.e-17, 1.e-19),
        "tolrdf": T(0.04,   0.02,   0.01),
        }
    del T

    accuracies = Enum(("low", "normal", "high",))

    #paral_modes = Enum("automatic", "explicit",)

    _mandatory_cards = ['SystemCard', 'KpointsCard',  'ElectronsCard']

    def __init__(self, *cards, **kwargs):

        super(ControlCard, self).__init__()

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
        if not set(names_found).issuperset(ControlCard._mandatory_cards):
            raise ValueError("Mandatory cards are missing, got %s" % str(names_found))

        # Extract useful cards.
        system    = [card for card in cards if card.name == "SystemCard"][0]
        electrons = [card for card in cards if card.name == "ElectronsCard"][0]

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
            "optdriver": self.optdriver,
            "ecut"     : ecut,
            "nbdbuf"   : nbdbuf,
            "nstep"    : [50,  75, 100][aidx],
            tol_varname: getattr(tol, self.accuracy),
            "prtwf"    : prtwf,
            "prtden"   : prtden,
            "prtgkk"   : prtgkk,
            "boxcutmin": boxcutmin,
            "optforces": None if (self.need_forces or want_forces) else 2,
            "optstress": None if (self.need_stress or want_stress) else 0,
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
        return self.has_card_names("RelaxCard")
                                                                            
    @property
    def need_stress(self):
        "True if the computation of the stress is required"
        # TODO: here it's easier to check if optcell != 0
        return self.has_card_names("RelaxCard")

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
        return ControlCard(*cards, **kwargs)

##########################################################################################

class Input(object):
    """
    Abstraction for the abinit input file: an input contains a dictionary of cards. 
    Each card is associated to a particular topic/calculation and gathers several objects 
    indexed by the class name e.g. card[__class__.__name__] --> object
    """
    # Mandatory cards that must be passed to the constructor.
    _mandatory_cards = ['SystemCard', 'KpointsCard',  'ElectronsCard', 'ControlCard',]

    def __init__(self, *cards):
        """
        Args:
            cards:
                List of Cards. Mandatory cards: SystemCard, KpointsCard, ElectronsCard, ControlCard
        """
        self._cards_dict = {}

        for card in cards:
            key = card.name
            if key in self._cards_dict:
                raise ValueError("Card %s is already in the dictionary" % card.name)

            self._cards_dict[key] = card 

        # Check for presence of mandatory cards.
        card_names = set([card.name for card in cards])
        if not set(card_names).issuperset(Input._mandatory_cards):
            raise ValueError("Mandatory cards are missing, got %s" % str(card_names))

        # Sanity check: no variable duplication is allowed.
        # TODO
        #counter = collections.Counter(self.varnames)
        #duplicated = filter(lambda t: t[1] > 1, counter.items())
        #                                                                          
        #if duplicated:
        #    raise ValueError("Found duplicated variables: %s" % str(duplicated))

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), str(self.card_names))

    def __str__(self):
        lines = [
            "# This input file has been generated by pymatgen",
            "# *** Do not edit. Any change will be lost if the input if regenerated ***"
            ]

        for k in sorted(self.card_names):
            card = self._cards_dict[k]
            lines.append( k.center(len(k)+2).center(90, "#") )
            lines.append(str(card))
            lines.append("")
                                                                                         
        return "\n".join(lines)

    @property
    def card_names(self):
        return [k for k in self._cards_dict]
                                         
    @property
    def cards(self):
        return self._cards_dict.values()

    @property
    def system_card(self):
        return self.cards_dict["SystemCard"]

    @property
    def electrons_card(self):
        return self.cards_dict["ElectronsCard"]

    @property
    def kpoints_card(self):
        return self.cards_dict["KpointsCard"]

    @property
    def control_card(self):
        return self.cards_dict["ControlCard"]

    @property
    def structure(self):
        "Pymatgen structure"
        return self.system_card.structure

    @property
    def pseudos(self):
        "List of pseudopotentials"
        return self.system_card.pseudos

    @property
    def spin_mode(self):
        "SpinMode object"
        return self.electrons_card.spin_mode

    @property
    def smearing(self):
        "Object defining the smearing technique"
        return self.electrons_card.smearing

    #@property
    #def description(self):
    #    "Description string"
    #    lines = []
    #    app = lines.append
    #    for (card_name, card) in self.items():
    #        card_comment = getattr(card, "comment", "No comment avaiable")
    #        app(card_name + ": " + str(card_comment))
    #    return "\n".join(lines)

    def copy(self):
        "Shallow copy"
        return Input(*self.cards)

    #def iter_variables(self):
    #    "Iterator over (varname, varvalue)"
    #    for card in self.cards:
    #        for obj in card.objects: 
    #            for (k, v) in obj.to_abivars():
    #                yield k,v

    #@property
    #def varnames(self):
    #    return [v[0] for v in self.variables]

    #@property
    #def variables(self):
    #    return list(self.iter_variables())

    #@property
    #def varnames2card(self):
    #    d = {}
    #    for (card_name, card) in self.items():
    #        d[card_name] = card.keys()

    #    vnames2card = {}
    #    for (card_name, vname_list) in d.items():
    #        for var_name in vname_list:

    #            if var_name == "_attributes": 
    #                continue

    #            if var_name in vnames2card:
    #                raise ValueError("variable %s occurs more than once" % var_name)

    #            vnames2card[var_name] = self[card_name]

    #    return vnames2card

    #def copy_cards(self, exclude=None, aslist=True):
    #    if exclude is None: 
    #        exclude = []
    #    else:
    #        # Check for typos!
    #        if not set(exclude).issubset(self.card_names):
    #            print("exclude_list: %s" % str(exclude))
    #            print("card_names: %s" % str(self.card_names))
    #            raise ValueError("Wrong value for exclude, likely a typo")
    #        
    #    # The entry "ControlCard" is automatically generated when we instanciate 
    #    # hence we don't add it to the cards.
    #    new_dict = {}
    #    for card_name in self.card_names:
    #        if card_name not in exclude + ["ControlCard"]: 
    #            new_dict[card_name] = self[card_name].copy()

    #    if aslist:
    #        return new_dict.values()
    #    else:
    #        return new_dict

    def add_vars_to_card(self, vars, card_name):
        "Add new_vars to the Card self[card_name]. Return new instance of Input."
        if card_name not in self.card_names:
            raise ValueError("card_name %s not present in Input" % card_name)

        print(card_name, vars)
                                                                                  
        cards = {}
        for (cname, card) in self._cards_dict.items():
            if cname == card_name:
                cards[cname] = card.copy()
                # FIXME
                cards[cname].update(vars)
            else:
                cards[cname] = card
                                                                                  
        return Input(*cards.values())
                                                                                  
    def add_vars_to_control(self, vars):
        return self.add_vars_to_card(vars, "ControlCard")
