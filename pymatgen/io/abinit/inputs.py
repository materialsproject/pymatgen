"""
This module defines a simplified interface for generating ABINIT input files.
The syntax is similar to the one used in ABINIT with small differences.
"""
from __future__ import print_function, division, unicode_literals, absolute_import

import os
import copy
import six
import abc
import json
import numpy as np

from collections import OrderedDict
from collections.abc import MutableMapping, Mapping
from monty.collections import dict2namedtuple
from monty.string import list_strings, is_string
from monty.json import MontyDecoder, MSONable
from pymatgen.util.serialization import pmg_serialize
from pymatgen.io.abinit.pseudos import PseudoTable, Pseudo
from pymatgen.io.abinit import abiobjects as aobj
from pymatgen.core.structure import Structure
from pymatgen.io.abinit.variable import InputVariable
#from abipy.core.structure import Structure
#from abipy.flowtk.abiobjects import structure_from_abivars, structure_to_abivars

import logging
logger = logging.getLogger(__file__)


# List of Abinit variables used to specify the structure.
# This variables should not be passed to set_vars since
# they will be generated with structure.to_abivars()
GEOVARS = set([
    "acell",
    "rprim",
    "rprimd"
    "angdeg",
    "xred",
    "xcart",
    "xangst",
    "znucl",
    "typat",
    "ntypat",
    "natom",
])

# Variables defining tolerances (used in pop_tolerances)
_TOLVARS = set([
    'toldfe',
    'tolvrs',
    'tolwfr',
    'tolrff',
    "toldff",
    "tolimg", # ?
    "tolmxf",
    "tolrde",
])

# Variables defining tolerances for the SCF cycle that are mutally exclusive
_TOLVARS_SCF = set([
    'toldfe',
    'tolvrs',
    'tolwfr',
    'tolrff',
    "toldff",
])

# Variables determining if data files should be read in input
_IRDVARS = set([
    "irdbseig",
    "irdbsreso",
    "irdhaydock",
    "irdddk",
    "irdden",
    "ird1den",
    "irdqps",
    "irdkss",
    "irdscr",
    "irdsuscep",
    "irdvdw",
    "irdwfk",
    "irdwfkfine",
    "irdwfq",
    "ird1wf",
])


class AbstractInput(six.with_metaclass(abc.ABCMeta, MutableMapping, object)):
    """
    Abstract class defining the methods that must be implemented by Input objects.
    """

    # ABC protocol: __delitem__, __getitem__, __iter__, __len__, __setitem__
    def __delitem__(self, key):
        return self.vars.__delitem__(key)

    def __getitem__(self, key):
        return self.vars.__getitem__(key)

    def __iter__(self):
        return self.vars.__iter__()

    def __len__(self):
        return len(self.vars)

    def __setitem__(self, key, value):
        self._check_varname(key)
        return self.vars.__setitem__(key, value)

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        return self.to_string()

    def write(self, filepath="run.abi"):
        """
        Write the input file to file to ``filepath``.
        """
        dirname = os.path.dirname(os.path.abspath(filepath))
        if not os.path.exists(dirname): os.makedirs(dirname)

        # Write the input file.
        with open(filepath, "wt") as fh:
            fh.write(str(self))

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def set_vars(self, *args, **kwargs):
        """
        Set the value of the variables.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue
        return kwargs

    def set_vars_ifnotin(self, *args, **kwargs):
        """
        Set the value of the variables but only if the variable is not already present.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        added = {}
        for varname, varvalue in kwargs.items():
            if varname not in self:
                self[varname] = varvalue
                added[varname] = varvalue
        return added

    def pop_vars(self, keys):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.
        Unlike remove_vars, no exception is raised if the variables are not in the input.

        Args:
            keys: string or list of strings with variable names.

        Example:
            inp.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])
        """
        return self.remove_vars(keys, strict=False)

    def remove_vars(self, keys, strict=True):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.

        Args:
            keys: string or list of strings with variable names.
            strict: If True, KeyError is raised if at least one variable is not present.
        """
        removed = {}
        for key in list_strings(keys):
            if strict and key not in self:
                raise KeyError("key: %s not in self:\n %s" % (key, list(self.keys())))
            if key in self:
                removed[key] = self.pop(key)

        return removed

    @abc.abstractproperty
    def vars(self):
        """Dictionary with the input variables. Used to implement dict-like interface."""

    @abc.abstractmethod
    def _check_varname(self, key):
        """Check if key is a valid name. Raise self.Error if not valid."""

    @abc.abstractmethod
    def to_string(self):
        """Returns a string with the input."""


class AbinitInputError(Exception):
    """Base error class for exceptions raised by ``BasicAbinitInput``."""


class BasicAbinitInput(six.with_metaclass(abc.ABCMeta, AbstractInput, MSONable, object)):
    """
    This object stores the ABINIT variables for a single dataset.
    """
    Error = AbinitInputError

    def __init__(self, structure, pseudos, pseudo_dir=None, comment=None, abi_args=None,
                 abi_kwargs=None):
        """
        Args:
            structure: Parameters defining the crystalline structure. Accepts |Structure| object
            file with structure (CIF, netcdf file, ...) or dictionary with ABINIT geo variables.
            pseudos: Pseudopotentials to be used for the calculation. Accepts: string or list of strings
                with the name of the pseudopotential files, list of |Pseudo| objects
                or |PseudoTable| object.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            abi_args: list of tuples (key, value) with the initial set of variables. Default: Empty
            abi_kwargs: Dictionary with the initial set of variables. Default: Empty
        """
        # Internal dict with variables. we use an ordered dict so that
        # variables will be likely grouped by `topics` when we fill the input.
        abi_args = [] if abi_args is None else abi_args
        for key, value in abi_args:
            self._check_varname(key)

        abi_kwargs = {} if abi_kwargs is None else abi_kwargs
        for key in abi_kwargs:
            self._check_varname(key)

        args = list(abi_args)[:]
        args.extend(list(abi_kwargs.items()))

        self._vars = OrderedDict(args)
        self.set_structure(structure)

        if pseudo_dir is not None:
            pseudo_dir = os.path.abspath(pseudo_dir)
            if not os.path.exists(pseudo_dir):
                raise self.Error("Directory %s does not exist" % pseudo_dir)
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        try:
            self._pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(self.structure)
        except ValueError as exc:
            raise self.Error(str(exc))

        if comment is not None: self.set_comment(comment)

    @pmg_serialize
    def as_dict(self):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        # Use a list of (key, value) to serialize the OrderedDict
        abi_args = []
        for key, value in self.items():
            if isinstance(value, np.ndarray): value = value.tolist()
            abi_args.append((key, value))

        return dict(structure=self.structure.as_dict(),
                    pseudos=[p.as_dict() for p in self.pseudos],
                    comment=self.comment,
                    abi_args=abi_args)

    @property
    def vars(self):
        return self._vars

    @classmethod
    def from_dict(cls, d):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        pseudos = [Pseudo.from_file(p['filepath']) for p in d['pseudos']]
        return cls(d["structure"], pseudos,
                   comment=d["comment"], abi_args=d["abi_args"])

    def add_abiobjects(self, *abi_objects):
        """
        This function receive a list of ``AbiVarable`` objects and add
        the corresponding variables to the input.
        """
        d = {}
        for aobj in abi_objects:
            if not hasattr(aobj, "to_abivars"):
                raise TypeError("type %s: %s does not have `to_abivars` method" % (type(aobj), repr(aobj)))
            d.update(self.set_vars(aobj.to_abivars()))
        return d

    def __setitem__(self, key, value):
        if key in _TOLVARS_SCF and hasattr(self, '_vars') and any(t in self._vars and t != key for t in _TOLVARS_SCF):
            logger.info("Replacing previously set tolerance variable: {0}."
                        .format(self.remove_vars(_TOLVARS_SCF, strict=False)))

        return super().__setitem__(key, value)

    def _check_varname(self, key):
        if key in GEOVARS:
            raise self.Error("You cannot set the value of a variable associated to the structure.\n"
                             "Use Structure objects to prepare the input file.")

    def to_string(self, post=None, mode="text",
                  with_structure=True, with_pseudos=True, exclude=None, verbose=0):
        r"""
        String representation.

        Args:
            mode: Either `text` or `html` if HTML output with links is wanted.
            post: String that will be appended to the name of the variables
                Note that post is usually autodetected when we have multiple datatasets
                It is mainly used when we have an input file with a single dataset
                so that we can prevent the code from adding "1" to the name of the variables
                (In this case, indeed, Abinit complains if ndtset=1 is not specified
                and we don't want ndtset=1 simply because the code will start to add
                _DS1_ to all the input and output files.
            with_structure: False if section with structure variables should not be printed.
            with_pseudos: False if JSON section with pseudo data should not be added.
            exclude: List of variable names that should be ignored.
        """
        lines = []
        app = lines.append

        if self.comment: app("# " + self.comment.replace("\n", "\n#"))

        post = post if post is not None else ""
        exclude = set(exclude) if exclude is not None else set()

        # Default is no sorting else alphabetical order.
        keys = sorted([k for k, v in self.items() if k not in exclude and v is not None])

        # Extract the items from the dict and add the geo variables at the end
        items = [(k, self[k]) for k in keys]
        if with_structure:
            items.extend(list(aobj.structure_to_abivars(self.structure).items()))

        for name, value in items:
            # Build variable, convert to string and append it
            vname = name + post
            app(str(InputVariable(vname, value)))

        s = "\n".join(lines)
        if not with_pseudos:
            return s

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")

        s += "\n#".join(ppinfo)
        return s

    @property
    def comment(self):
        """Optional string with comment. None if comment is not set."""
        try:
            return self._comment
        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    @property
    def structure(self):
        """The |Structure| object associated to this input."""
        return self._structure

    def set_structure(self, structure):
        """Set structure."""
        #self._structure = Structure.as_structure(structure)
        self._structure = as_structure(structure)

        # Check volume
        m = self.structure.lattice.matrix
        if np.dot(np.cross(m[0], m[1]), m[2]) <= 0:
            raise self.Error("The triple product of the lattice vector is negative. Use structure.abi_sanitize.")

        return self._structure

    # Helper functions to facilitate the specification of several variables.
    def set_kmesh(self, ngkpt, shiftk, kptopt=1):
        """
        Set the variables for the sampling of the BZ.

        Args:
            ngkpt: Monkhorst-Pack divisions
            shiftk: List of shifts.
            kptopt: Option for the generation of the mesh.
        """
        shiftk = np.reshape(shiftk, (-1, 3))
        return self.set_vars(ngkpt=ngkpt, kptopt=kptopt, nshiftk=len(shiftk), shiftk=shiftk)

    def set_gamma_sampling(self):
        """Gamma-only sampling of the BZ."""
        return self.set_kmesh(ngkpt=(1, 1, 1), shiftk=(0, 0, 0))

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2):
        """
        Set the variables for the computation of the electronic band structure.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            kptbounds: k-points defining the path in k-space.
                If None, we use the default high-symmetry k-path defined in the pymatgen database.
        """
        if kptbounds is None: kptbounds = self.structure.calc_kptbounds()
        kptbounds = np.reshape(kptbounds, (-1,3))
        #self.pop_vars(["ngkpt", "shiftk"]) ??

        return self.set_vars(kptbounds=kptbounds, kptopt=-(len(kptbounds)-1), ndivsm=ndivsm, iscf=iscf)

    def set_kptgw(self, kptgw, bdgw):
        """
        Set the variables (k-points, bands) for the computation of GW corrections.

        Args:
            kptgw: List of k-points in reduced coordinates.
            bdgw: Specifies the range of bands for the GW corrections.
                Accepts iterable that be reshaped to (nkptgw, 2)
                or a tuple of two integers if the extrema are the same for each k-point.
        """
        kptgw = np.reshape(kptgw, (-1,3))
        nkptgw = len(kptgw)
        if len(bdgw) == 2: bdgw = len(kptgw) * bdgw

        return self.set_vars(kptgw=kptgw, nkptgw=nkptgw, bdgw=np.reshape(bdgw, (nkptgw, 2)))

    def set_spin_mode(self, spin_mode):
        """
        Set the variables used to the treat the spin degree of freedom.
        Return dictionary with the variables that have been removed.

        Args:
            spin_mode: :class:`SpinMode` object or string. Possible values for string are:

            - polarized
            - unpolarized
            - afm (anti-ferromagnetic)
            - spinor (non-collinear magnetism)
            - spinor_nomag (non-collinear, no magnetism)
        """
        # Remove all variables used to treat spin
        old_vars = self.pop_vars(["nsppol", "nspden", "nspinor"])
        self.add_abiobjects(aobj.SpinMode.as_spinmode(spin_mode))
        return old_vars

    @property
    def pseudos(self):
        """List of |Pseudo| objects."""
        return self._pseudos

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    def new_with_vars(self, *args, **kwargs):
        """
        Return a new input with the given variables.

        Example:
            new = input.new_with_vars(ecut=20)
        """
        # Avoid modifications in self.
        new = self.deepcopy()
        new.set_vars(*args, **kwargs)
        return new

    def pop_tolerances(self):
        """
        Remove all the tolerance variables present in self.
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_TOLVARS, strict=False)

    def pop_irdvars(self):
        """
        Remove all the `ird*` variables present in self.
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_IRDVARS, strict=False)



def as_structure(obj):
    """
    Convert obj into a Structure. Accepts:

        - Structure object.
        - Filename
        - Dictionaries (MSONable format or dictionaries with abinit variables).
    """
    if isinstance(obj, Structure):
        return obj

    if is_string(obj):
        return Structure.from_file(obj)

    if isinstance(obj, Mapping):
        if "@module" in obj:
            return Structure.from_dict(obj)
        else:
            return aobj.structure_from_abivars(cls=None, **obj)

    raise TypeError("Don't know how to convert %s into a structure" % type(obj))
