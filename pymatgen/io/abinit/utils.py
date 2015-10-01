# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""Tools and helper functions for abinit calculations"""
from __future__ import unicode_literals, division, print_function

import os
import six
import re
import collections
import shutil
import operator
import numpy as np

from fnmatch import fnmatch
from six.moves import filter
from monty.collections import dict2namedtuple
from monty.string import list_strings
from monty.fnmatch import WildCard

import logging
logger = logging.getLogger(__name__)


def as_bool(s):
    """
    Convert a string into a boolean.

    >>> assert as_bool(True) is True and as_bool("Yes") is True and as_bool("false") is False
    """
    if s in (False, True): return s
    # Assume string
    s = s.lower()
    if s in ("yes", "true"): 
        return True
    elif s in ("no", "false"): 
        return False
    else:
        raise ValueError("Don't know how to convert type %s: %s into a boolean" % (type(s), s))


class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.
    Provides wrappers for the most commonly used functions defined in os.path.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    def __str__(self):
        return "<%s, %s>" % (self.__class__.__name__, self.path)

    def __eq__(self, other):
        return False if other is None else self.path == other.path
                                       
    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def path(self):
        """Absolute path of the file."""
        return self._path

    @property
    def basename(self):
        """File basename."""
        return os.path.basename(self.path)

    @property
    def relpath(self):
        """Relative path."""
        try:
            return os.path.relpath(self.path)
        except OSError:
            # current working directory may not be defined!
            return self.path

    @property
    def dirname(self):
        """Absolute path of the directory where the file is located."""
        return os.path.dirname(self.path)

    @property
    def exists(self):
        """True if file exists."""
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        """True if self is a NetCDF file"""
        return self.basename.endswith(".nc")

    def chmod(self, mode):
        """Change the access permissions of a file."""
        os.chmod(self.path, mode)

    def read(self):
        """Read data from file."""
        with open(self.path, "r") as f:
            return f.read()

    def readlines(self):
        """Read lines from files."""
        with open(self.path, "r") as f:
            return f.readlines()

    def write(self, string):
        """Write string to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            if not string.endswith("\n"):
                return f.write(string + "\n")
            else:
                return f.write(string)

    def writelines(self, lines):
        """Write a list of strings to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            return f.writelines(lines)

    def make_dir(self):
        """Make the directory where the file is located."""
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)

    def remove(self):
        """Remove the file."""
        try:
            os.remove(self.path)
        except:
            pass

    def move(self, dst):
        """
        Recursively move a file or directory to another location. This is
        similar to the Unix "mv" command.
        """
        shutil.move(self.path, dst)

    def get_stat(self):
        """Results from os.stat"""
        return os.stat(self.path)


class Directory(object):
    """
    Very simple class that provides helper functions
    wrapping the most commonly used functions defined in os.path.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    def __str__(self):
        return self.path

    def __eq__(self, other):
        return False if other is None else self.path == other.path

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def path(self):
        """Absolute path of the directory."""
        return self._path

    @property
    def relpath(self):
        """Relative path."""
        return os.path.relpath(self.path)

    @property
    def basename(self):
        """Directory basename."""
        return os.path.basename(self.path)

    def path_join(self, *p):
        """
        Join two or more pathname components, inserting '/' as needed.
        If any component is an absolute path, all previous path components will be discarded.
        """
        return os.path.join(self.path, *p)

    @property
    def exists(self):
        """True if file exists."""
        return os.path.exists(self.path)

    def makedirs(self):
        """
        Super-mkdir; create a leaf directory and all intermediate ones.
        Works like mkdir, except that any intermediate path segment (not
        just the rightmost) will be created if it does not exist.
        """
        if not self.exists:
            os.makedirs(self.path)

    def rmtree(self):
        """Recursively delete the directory tree"""
        shutil.rmtree(self.path, ignore_errors=True)

    def clean(self):
        """Remove all files in the directory tree while preserving the directory"""
        for path in self.list_filepaths():
            try:
                os.remove(path)
            except:
                pass

    def path_in(self, file_basename):
        """Return the absolute path of filename in the directory."""
        return os.path.join(self.path, file_basename)

    def list_filepaths(self, wildcard=None):
        """
        Return the list of absolute filepaths in the directory.

        Args:
            wildcard: String of tokens separated by "|". Each token represents a pattern.
                If wildcard is not None, we return only those files that match the given shell pattern (uses fnmatch).
                Example:
                  wildcard="*.nc|*.pdf" selects only those files that end with .nc or .pdf
        """
        # Select the files in the directory.
        fnames = [f for f in os.listdir(self.path)]
        filepaths = filter(os.path.isfile, [os.path.join(self.path, f) for f in fnames])

        # Filter using the shell patterns.
        if wildcard is not None:
            filepaths = WildCard(wildcard).filter(filepaths)

        return filepaths

    def has_abiext(self, ext):
        """
        Returns the absolute path of the ABINIT file with extension ext.
        Support both Fortran files and netcdf files. In the later case,
        we check whether a file with extension ext + ".nc" is present 
        in the directory. Returns empty string is file is not present.

        Raises:
            `ValueError` if multiple files with the given ext are found.
            This implies that this method is not compatible with multiple datasets.
        """

        ext = ext if ext.startswith('_') else '_'+ext

        files = []
        for f in self.list_filepaths():
            if f.endswith(ext) or f.endswith(ext + ".nc"):
                files.append(f)

        # This should fix the problem with the 1WF files in which the file extension convention is broken
        if not files:
            files = [f for f in self.list_filepaths() if fnmatch(f, "*%s*" % ext)]

        if not files:
            return ""

        if len(files) > 1:
            # ABINIT users must learn that multiple datasets are bad!
            err_msg = "Found multiple files with the same extensions:\n %s\nPlease avoid the use of mutiple datasets!" % files
            raise ValueError(err_msg)

        return files[0]

    def symlink_abiext(self, inext, outext):
        """Create a simbolic link"""
        infile = self.has_abiext(inext)
        if not infile:
            raise RuntimeError('no file with extension %s in %s' % (inext, self))

        for i in range(len(infile) - 1, -1, -1):
            if infile[i] == '_':
                break
        else:
            raise RuntimeError('Extension %s could not be detected in file %s' % (inext, infile))

        outfile = infile[:i] + '_' + outext
        os.symlink(infile, outfile)
        return 0

    def rename_abiext(self, inext, outext):
        """Rename the Abinit file with extension inext with the new extension outext"""
        infile = self.has_abiext(inext)
        if not infile:
            raise RuntimeError('no file with extension %s in %s' % (inext, self))

        for i in range(len(infile) - 1, -1, -1):
            if infile[i] == '_':
                break
        else:
            raise RuntimeError('Extension %s could not be detected in file %s' % (inext, infile))

        outfile = infile[:i] + '_' + outext
        shutil.move(infile, outfile)
        return 0

    def copy_abiext(self, inext, outext):
        """Copy the Abinit file with extension inext to a new file withw extension outext"""
        infile = self.has_abiext(inext)
        if not infile:
            raise RuntimeError('no file with extension %s in %s' % (inext, self))

        for i in range(len(infile) - 1, -1, -1):
            if infile[i] == '_':
                break
        else:
            raise RuntimeError('Extension %s could not be detected in file %s' % (inext, infile))

        outfile = infile[:i] + '_' + outext
        shutil.copy(infile, outfile)
        return 0

    def remove_exts(self, exts):
        """
        Remove the files with the given extensions. Unlike rmtree, this function preserves the directory path.
        Return list with the absolute paths of the files that have been removed.
        """
        paths = []

        for ext in list_strings(exts):
            path = self.has_abiext(ext)
            if not path: continue
            try:
                os.remove(path)
                paths.append(path)
            except IOError:
                logger.warning("Exception while trying to remove file %s" % path)

        return paths

    def find_last_timden_file(self):
        """
        ABINIT produces lots of out_TIM1_DEN files for each step and we need to find the lat
        one in order to prepare the restart or to connect other tasks to the structural relaxation.

        This function finds all the TIM?_DEN files in self and return a namedtuple (path, step)
        where `path` is the path of the last TIM?_DEN file and step is the iteration number.
        Returns None if the directory does not contain TIM?_DEN files.
        """
        regex = re.compile("out_TIM(\d+)_DEN(.nc)?$")

        timden_paths = [f for f in self.list_filepaths() if regex.match(os.path.basename(f))]
        if not timden_paths: return None

        # Build list of (step, path) tuples.
        stepfile_list = []
        for path in timden_paths:
            name = os.path.basename(path)
            match = regex.match(name)
            step, ncext = match.groups()
            stepfile_list.append((int(step), path))
                                                                                        
        # DSU sort.
        last = sorted(stepfile_list, key=lambda t: t[0])[-1]
        return dict2namedtuple(step=last[0], path=last[1])

    def find_1wf_files(self):
        """
        Abinit adds the idir-ipert index at the end of the 1WF file and this breaks the extension 
        e.g. out_1WF4. This method scans the files in the directories and returns a list of namedtuple
        Each named tuple gives the `path` of the 1FK file and the `pertcase` index.
        """
        regex = re.compile("out_1WF(\d+)(.nc)?$")

        wf_paths = [f for f in self.list_filepaths() if regex.match(os.path.basename(f))]
        if not wf_paths: return None

        # Build list of (pertcase, path) tuples.
        pertfile_list = []
        for path in wf_paths:
            name = os.path.basename(path)
            match = regex.match(name)
            pertcase, ncext = match.groups()
            pertfile_list.append((int(pertcase), path))
                                                                                        
        # DSU sort.
        pertfile_list = sorted(pertfile_list, key=lambda t: t[0])
        return [dict2namedtuple(pertcase=item[0], path=item[1]) for item in pertfile_list]

    def find_1den_files(self):
        """
        Abinit adds the idir-ipert index at the end of the 1DEN file and this breaks the extension 
        e.g. out_DEN1. This method scans the files in the directories and returns a list of namedtuple
        Each named tuple gives the `path` of the 1DEN file and the `pertcase` index.
        """
        regex = re.compile("out_DEN(\d+)(.nc)?$")
        den_paths = [f for f in self.list_filepaths() if regex.match(os.path.basename(f))]
        if not den_paths: return None

        # Build list of (pertcase, path) tuples.
        pertfile_list = []
        for path in den_paths:
            name = os.path.basename(path)
            match = regex.match(name)
            pertcase, ncext = match.groups()
            pertfile_list.append((int(pertcase), path))
                                                                                        
        # DSU sort.
        pertfile_list = sorted(pertfile_list, key=lambda t: t[0])
        return [dict2namedtuple(pertcase=item[0], path=item[1]) for item in pertfile_list]

# This dictionary maps ABINIT file extensions to the variables that must be used to read the file in input.
#
# TODO: It would be nice to pass absolute paths to abinit with getden_path
# so that I can avoid creating symbolic links before running but
# the presence of the C-bindings complicates the implementation
# (gfortran SIGFAULTs if I add strings to dataset_type!
_EXT2VARS = {
    "DEN": {"irdden": 1},
    "WFK": {"irdwfk": 1},
    "WFQ": {"irdwfq": 1},
    "SCR": {"irdscr": 1},
    "QPS": {"irdqps": 1},
    "1WF": {"ird1wf": 1},
    "1DEN": {"ird1den": 1},
    "BSR": {"irdbsreso": 1},
    "BSC": {"irdbscoup": 1},
    "HAYDR_SAVE": {"irdhaydock": 1},
    "DDK": {"irdddk": 1},
    "DDB": {},
    "GKK": {},
    "DKK": {},
}


def irdvars_for_ext(ext):
    """
    Returns a dictionary with the ABINIT variables 
    that must be used to read the file with extension ext.
    """
    return _EXT2VARS[ext].copy()


def abi_extensions():
    """List with all the ABINIT extensions that are registered."""
    return list(_EXT2VARS.keys())[:]


def abi_splitext(filename):
    """
    Split the ABINIT extension from a filename.
    "Extension" are found by searching in an internal database.

    Returns "(root, ext)" where ext is the registered ABINIT extension 
    The final ".nc" is included (if any) 

    >>> assert abi_splitext("foo_WFK") == ('foo_', 'WFK')
    >>> assert abi_splitext("/home/guido/foo_bar_WFK.nc") == ('foo_bar_', 'WFK.nc')
    """
    filename = os.path.basename(filename)
    is_ncfile = False
    if filename.endswith(".nc"):
        is_ncfile = True
        filename = filename[:-3]

    known_extensions = abi_extensions()

    # This algorith fails if we have two files 
    # e.g. HAYDR_SAVE, ANOTHER_HAYDR_SAVE
    for i in range(len(filename)-1, -1, -1):
        ext = filename[i:]
        if ext in known_extensions:
            break

    else:
        raise ValueError("Cannot find a registered extension in %s" % filename)

    root = filename[:i]
    if is_ncfile: 
        ext += ".nc"

    return root, ext


class FilepathFixer(object):
    """
    This object modifies the names of particular output files
    produced by ABINIT so that the file extension is preserved.
    Having a one-to-one mapping between file extension and data format
    is indeed fundamental for the correct behaviour of abinitio since:

        - We locate the output file by just inspecting the extension

        - We select the variables that must be added to the input file
          on the basis of the extension specified by the user during 
          the initialization of the `AbinitFlow`.

    Unfortunately, ABINIT developers like to append extra stuff 
    to the initial extension and therefore we have to call 
    `FilepathFixer` to fix the output files produced by the run.

    Example:
    
    >>> fixer = FilepathFixer()
    >>> assert fixer.fix_paths('/foo/out_1WF17') == {'/foo/out_1WF17': '/foo/out_1WF'}
    >>> assert fixer.fix_paths('/foo/out_1WF5.nc') == {'/foo/out_1WF5.nc': '/foo/out_1WF.nc'}
    """
    def __init__(self):
        # dictionary mapping the *official* file extension to
        # the regular expression used to tokenize the basename of the file
        # To add a new fix it's sufficient to add a new regexp and 
        # a static method _fix_EXTNAME
        self.regs = regs = {}
        import re
        regs["1WF"] = re.compile("(\w+_)1WF(\d+)(.nc)?$")
        regs["1DEN"] = re.compile("(\w+_)1DEN(\d+)(.nc)?$")

    @staticmethod
    def _fix_1WF(match):
        root, pert, ncext = match.groups()
        if ncext is None: ncext = ""
        return root + "1WF" + ncext

    @staticmethod
    def _fix_1DEN(match):
        root, pert, ncext = match.groups()
        if ncext is None: ncext = ""
        return root + "1DEN" + ncext

    def _fix_path(self, path):
        for ext, regex in self.regs.items():
            head, tail = os.path.split(path)

            match = regex.match(tail)
            if match:
                newtail = getattr(self, "_fix_" + ext)(match)
                newpath = os.path.join(head, newtail)
                return newpath, ext

        return None, None

    def fix_paths(self, paths):
        """
        Fix the filenames in the iterable paths

        Returns:
            old2new: Mapping old_path --> new_path
        """
        old2new, fixed_exts = {}, []

        for path in list_strings(paths):
            newpath, ext = self._fix_path(path)

            if newpath is not None:
                #if ext not in fixed_exts:
                #    if ext == "1WF": continue
                #    raise ValueError("Unknown extension %s" % ext)
                #print(ext, path, fixed_exts)
                #if ext != '1WF':
                #    assert ext not in fixed_exts
                if ext not in fixed_exts:
                    if ext == "1WF": continue
                    raise ValueError("Unknown extension %s" % ext)
                fixed_exts.append(ext)
                old2new[path] = newpath

        return old2new


def _bop_not(obj):
    """Boolean not."""
    return not bool(obj)


def _bop_and(obj1, obj2):
    """Boolean and."""
    return bool(obj1) and bool(obj2)


def _bop_or(obj1, obj2):
    """Boolean or."""
    return bool(obj1) or bool(obj2)


def _bop_divisible(num1, num2):
    """Return True if num1 is divisible by num2."""
    return (num1 % num2) == 0.0

# Mapping string --> operator.
_UNARY_OPS = {
    "$not": _bop_not,
}

_BIN_OPS = {
    "$eq": operator.eq,
    "$ne": operator.ne,
    "$gt": operator.gt,
    "$ge": operator.ge,
    "$lt": operator.lt,
    "$le": operator.le,
    "$divisible": _bop_divisible,
    "$and": _bop_and,
    "$or":  _bop_or,
    }


_ALL_OPS = list(_UNARY_OPS.keys()) + list(_BIN_OPS.keys())


def map2rpn(map, obj):
    """
    Convert a Mongodb-like dictionary to a RPN list of operands and operators.

    Reverse Polish notation (RPN) is a mathematical notation in which every 
    operator follows all of its operands, e.g.

    3 - 4 + 5 -->   3 4 - 5 + 

    >>> d = {2.0: {'$eq': 1.0}}
    >>> assert map2rpn(d, None) == [2.0, 1.0, '$eq']
    """
    rpn = []

    for k, v in map.items():

        if k in _ALL_OPS:
            if isinstance(v, collections.Mapping):
                # e.g "$not": {"$gt": "one"}
                # print("in op_vmap",k, v)
                values = map2rpn(v, obj)
                rpn.extend(values)
                rpn.append(k)

            elif isinstance(v, (list, tuple)):
                # e.g "$and": [{"$not": {"one": 1.0}}, {"two": {"$lt": 3}}]}
                # print("in_op_list",k, v)
                for d in v:
                    rpn.extend(map2rpn(d, obj))

                rpn.append(k)

            else: 
                # Examples
                # 1) "$eq"": "attribute_name"
                # 2) "$eq"": 1.0
                try:
                    #print("in_otherv",k, v)
                    rpn.append(getattr(obj, v))
                    rpn.append(k)

                except TypeError:
                    #print("in_otherv, raised",k, v)
                    rpn.extend([v, k])
        else:
            try:
                k = getattr(obj, k)
            except TypeError:
                k = k

            if isinstance(v, collections.Mapping):
                # "one": {"$eq": 1.0}}
                values = map2rpn(v, obj)
                rpn.append(k)
                rpn.extend(values) 
            else:
                #"one": 1.0 
                rpn.extend([k, v, "$eq"])

    return rpn


def evaluate_rpn(rpn):
    """
    Evaluates the RPN form produced my map2rpn.

    Returns:
        bool
    """
    vals_stack = []
                                                                              
    for item in rpn:
                                                                              
        if item in _ALL_OPS:
            # Apply the operator and push to the task.
            v2 = vals_stack.pop()
                                                                              
            if item in _UNARY_OPS:
                res = _UNARY_OPS[item](v2)
                                                                              
            elif item in _BIN_OPS:
                v1 = vals_stack.pop()
                res = _BIN_OPS[item](v1, v2)
            else:
                raise ValueError("%s not in unary_ops or bin_ops" % str(item))
                                                                              
            vals_stack.append(res)
                                                                              
        else:
            # Push the operand
            vals_stack.append(item)
                                                                              
        #print(vals_stack)
                                                                              
    assert len(vals_stack) == 1
    assert isinstance(vals_stack[0], bool)
                                                                              
    return vals_stack[0]


class Condition(object):
    """
    This object receives a dictionary that defines a boolean condition whose syntax is similar
    to the one used in mongodb (albeit not all the operators available in mongodb are supported here).

    Example:

    $gt: {field: {$gt: value} }

    $gt selects those documents where the value of the field is greater than (i.e. >) the specified value.

    $and performs a logical AND operation on an array of two or more expressions (e.g. <expression1>, <expression2>, etc.) 
    and selects the documents that satisfy all the expressions in the array. 

    { $and: [ { <expression1> }, { <expression2> } , ... , { <expressionN> } ] }

    Consider the following example:

    db.inventory.find( { qty: { $gt: 20 } } )
    This query will select all documents in the inventory collection where the qty field value is greater than 20.
    Consider the following example:

    db.inventory.find( { qty: { $gt: 20 } } )
    db.inventory.find({ $and: [ { price: 1.99 }, { qty: { $lt: 20 } }, { sale: true } ] } )
    """
    @classmethod
    def as_condition(cls, obj):
        """Convert obj into :class:`Condition`"""
        if isinstance(obj, cls):
            return obj
        else:
            return cls(cmap=obj)

    def __init__(self, cmap=None):
        self.cmap = {} if cmap is None else cmap

    def __str__(self):
        return str(self.cmap)

    def __bool__(self):
        return bool(self.cmap)

    __nonzero__ = __bool__

    def __call__(self, obj):
        if not self: return True
        try:
            return evaluate_rpn(map2rpn(self.cmap, obj))
        except Exception as exc:
            logger.warning("Condition(%s) raised Exception:\n %s" % (type(obj), str(exc)))
            return False


class Editor(object):
    """
    Wrapper class that calls the editor specified by the user 
    or the one specified in the $EDITOR env variable.
    """
    def __init__(self, editor=None):
        """If editor is None, $EDITOR is used."""
        self.editor = os.getenv("EDITOR", "vi") if editor is None else str(editor)

    def edit_files(self, fnames, ask_for_exit=True):
        exit_status = 0
        for idx, fname in enumerate(fnames):
            exit_status = self.edit_file(fname)
            if ask_for_exit and idx != len(fnames)-1 and self.user_wants_to_exit():
                break
        return exit_status

    def edit_file(self, fname):
        from subprocess import call
        retcode = call([self.editor, fname])

        if retcode != 0:
            import warnings
            warnings.warn("Error while trying to edit file: %s" % fname)

        return retcode 

    @staticmethod
    def user_wants_to_exit():
        """Show an interactive prompt asking if exit is wanted."""
        # Fix python 2.x.
        if six.PY2:
            my_input = raw_input
        else:
            my_input = input

        try:
            answer = my_input("Do you want to continue [Y/n]")
        except EOFError:
            return True

        return answer.lower().strip() in ["n", "no"]


class SparseHistogram(object):

    def __init__(self, items, key=None, num=None, step=None):
        if num is None and step is None:
            raise ValueError("Either num or step must be specified")

        from collections import defaultdict, OrderedDict

        values = [key(item) for item in items] if key is not None else items
        start, stop = min(values), max(values)
        if num is None:
            num = int((stop - start) / step)
            if num == 0: num = 1
        mesh = np.linspace(start, stop, num, endpoint=False)

        from monty.bisect import find_le

        hist = defaultdict(list)
        for item, value in zip(items, values):
            # Find rightmost value less than or equal to x.
            # hence each bin contains all items whose value is >= value
            pos = find_le(mesh, value)
            hist[mesh[pos]].append(item)

        #new = OrderedDict([(pos, hist[pos]) for pos in sorted(hist.keys(), reverse=reverse)])
        self.binvals = sorted(hist.keys())
        self.values = [hist[pos] for pos in self.binvals]
        self.start, self.stop, self.num = start, stop, num

    from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
    @add_fig_kwargs
    def plot(self, ax=None, **kwargs):
        """
        Plot the histogram with matplotlib, returns `matplotlib figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        yy = [len(v) for v in self.values]
        ax.plot(self.binvals, yy, **kwargs)

        return fig
