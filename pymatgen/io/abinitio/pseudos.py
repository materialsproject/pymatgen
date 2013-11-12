"""
This module provides objects describing the basic parameters of the 
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..
"""
from __future__ import division, print_function

import sys
import os
import abc
import collections
import json
import warnings
import cPickle as pickle
import cStringIO as StringIO
import numpy as np

from pymatgen.core.design_patterns import FrozenDict, AttrDict
from pymatgen.core.periodic_table import PeriodicTable
from pymatgen.util.num_utils import iterator_from_slice
from pymatgen.util.string_utils import list_strings, is_string

__all__ = [
    "Pseudo",
    "PseudoTable",
]

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"

# Tools and helper functions.
def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def _read_nlines(filename, nlines):
    """
    Read at most nlines lines from file filename.
    If nlines is < 0, the entire file is read.
    """
    if nlines < 0:
        with open(filename, 'r') as fh:
            return fh.readlines()

    lines = []
    with open(filename, 'r') as fh:
        for (lineno, line) in enumerate(fh):
            if lineno == nlines: break
            lines.append(line)
        return lines

_l2str = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
    4: "g",
    5: "h",
    6: "i",
}

_str2l = {v: k for k, v in _l2str.items()}

def l2str(l):
    """Convert the angular momentum l (int) to string."""
    try:
        return _l2str[l]
    except KeyError:
        return "Unknown angular momentum, received l = %s" % l


def str2l(s):
    """Convert a string to the angular momentum l (int)"""
    return _str2l[s]


def read_dojo_report(filename):
    """Helper function to read the DOJO_REPORT from file."""
    with open(filename, "r") as fh:
        lines = fh.readlines()
        try:
            start = lines.index("<DOJO_REPORT>\n")
        except ValueError:
            return {}

        stop = lines.index("</DOJO_REPORT>\n")
        return json.loads("".join(lines[start+1:stop]))

#class DojoReport(dict):
#    _LATEST_VERSION = 1.0
#    _START_LINE = "<DOJO_REPORT>\n"
#    _END_LINE = "</DOJO_REPORT>\n"
#
#    @classmethod
#    def from_file(cls, path):
#        new = read_dojo_report(path)
#        new.__class__ = cls
#        return new
#
#    #def to_file(self, path):


_PTABLE = PeriodicTable()

class Pseudo(object):
    """
    Abstract base class defining the methods that must be 
    implemented by the concrete pseudopotential classes.
    """
    __metaclass__ = abc.ABCMeta

    #def __init__(self, filepath):
    #    self.filepath = os.path.abspath(filepath)
    #    self._dojo_report = {}

    @classmethod
    def aspseudo(cls, obj):
        """
        Convert obj into a pseudo. Accepts:

            * Pseudo object.
            * string defining a valid path.
        """
        if isinstance(obj, cls):
            return obj
        else:
            # Assumes path.
            return cls.from_file(obj)

    @staticmethod
    def from_file(filename):
        """
        Return a pseudopotential object from filename.
        Note: the parser knows the concrete class that should be instanciated
        """
        return PseudoParser().parse(filename)

    def __repr__(self):
        return "<%s at %s, name = %s>" % (
            self.__class__.__name__, id(self), self.name)

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("<%s: %s>" % (self.__class__.__name__, self.name))
        app("  summary: " + self.summary.strip())
        app("  number of valence electrons: %s" % self.Z_val)
        #FIXME: rewrite the treatment of xc, use XML specs as starting point
        #app("  XC correlation (ixc): %s" % self._pspxc)  #FIXME
        app("  maximum angular momentum: %s" % l2str(self.l_max))
        app("  angular momentum for local part: %s" % l2str(self.l_local))
        if self.isnc:
            app("  radius for non-linear core correction: %s" % self.nlcc_radius)
        app("")

        hint_normal = self.hint_for_accuracy()
        if hint_normal is not None:
            app("  hint for normal accuracy: %s" % str(hint_normal))

        return "\n".join(lines)

    @abc.abstractproperty
    def summary(self):
        """String summarizing the most important properties."""

    @property
    def filepath(self):
        return os.path.abspath(self.path)

    @property
    def name(self):
        """File basename."""
        return os.path.basename(self.filepath)

    @abc.abstractproperty
    def Z(self):
        """The atomic number of the atom."""

    @abc.abstractproperty
    def Z_val(self):
        """Valence charge"""

    @property
    def element(self):
        """Pymatgen `Element`."""
        try:
            return _PTABLE[self.Z]
        except (KeyError, IndexError):
            return _PTABLE[int(self.Z)]

    @property
    def symbol(self):
        """Element symbol."""
        return self.element.symbol

    @abc.abstractproperty
    def l_max(self):
        """Maximum angular momentum."""

    @abc.abstractproperty
    def l_local(self):
        """Angular momentum used for the local part."""

    @property
    def isnc(self):
        """True if norm-conserving pseudopotential."""
        return isinstance(self, NcPseudo)

    @property
    def ispaw(self):
        """True if PAW pseudopotential."""
        return isinstance(self, PawPseudo)

    #@abc.abstractproperty
    #def xc_type(self):
    #    """XC family e.g LDA, GGA, MGGA."""

    #@abc.abstractproperty
    #def xc_flavor(self):
    #    """XC flavor e.g PW, PW91, PBE."""

    #@property
    #def xc_functional(self):
    #    """XC identifier e.g LDA-PW91, GGA-PBE, GGA-revPBE."""
    #    return "-".join([self.xc_type, self.xc_flavor])

    #@abc.abstractproperty
    #def has_soc(self):
    #    """True if pseudo contains spin-orbit coupling."""

    #@abc.abstractmethod
    #def num_of_projectors(self, l='s'):
    #    """Number of projectors for the angular channel l"""

    #@abc.abstractmethod
    #def generation_mode
    #    """scalar scalar-relativistic, relativistic."""

    @property
    def has_dojo_report(self):
        """True if self contains the DOJO_REPORT section."""
        return bool(self.dojo_report)

    def delta_factor(self, accuracy="normal"):
        """
        Returns the deltafactor [meV/natom] computed with the given accuracy.
        None if self does not have info on the deltafactor.
        """
        if not self.has_dojo_report:
            return None
        try:
            return self.dojo_report["delta_factor"][accuracy]["dfact"]
        except KeyError:
            return None

    def read_dojo_report(self):
        """Read the DOJO_REPORT section, returns {} if section is not present.""" 
        return read_dojo_report(self.path)

    def write_dojo_report(self, report):
        """Write a new DOJO_REPORT section to the pseudopotential file."""
        path = self.path

        # Create JSON string from report.
        jstring = json.dumps(report, indent=4, sort_keys=True) + "\n"

        # Read lines from file and insert jstring between the tags.
        with open(path, "r") as fh:
            lines = fh.readlines()
            try:
                start = lines.index("<DOJO_REPORT>\n")
            except ValueError:
                start = -1

            if start == -1:
               # DOJO_REPORT was not present.
               lines += ["\n", "<DOJO_REPORT>\n", jstring , "</DOJO_REPORT>\n",]
            else:
               stop = lines.index("</DOJO_REPORT>\n")
               lines.insert(stop, jstring)
               del lines[start+1:stop]

        #  Write new file.
        with open(path, "w") as fh:
            fh.writelines(lines)

    def remove_dojo_report(self):
        """Remove the DOJO_REPORT section from the pseudopotential file."""
        # Read lines from file and insert jstring between the tags.
        path = self.path
        with open(path, "r") as fh:
            lines = fh.readlines()
            try:
                start = lines.index("<DOJO_REPORT>\n")
            except ValueError:
                start = -1

            if start == -1:
               return

            stop = lines.index("</DOJO_REPORT>\n")
            if stop == -1:
               return

            del lines[start+1:stop]

        # Write new file.
        with open(path, "w") as fh:
            fh.writelines(lines)

    def hint_for_accuracy(self, accuracy="normal"):
        """
        Returns an hint object with parameters such as ecut [Ha] and 
        aug_ratio for given accuracy. Returns None if no hint is available.

        Args:
            accuracy: ["low", "normal", "high"]
        """
        if self.has_dojo_report:
            return Hint.from_dict(self.dojo_report["hints"][accuracy])
        else:
            return None

    @property
    def has_hints(self):
        """True if self provides hints on the cutoff energy."""
        for acc in ["low", "normal", "high"]:
            if self.hint_for_accuracy(acc) is None:
                return False
        return True

    #@property
    #def md5(self):
    #    """
    #    Return the checksum of the pseudopotential file.
    #    """
    #    import hashlib
    #    hasher = hashlib.md5()
    #    with open(self.filepath, "r") as fh:
    #        hasher.update(fh.read())
    #        return hasher.hexdigest()


class NcPseudo(object):
    """
    Abstract class defining the methods that must be implemented
    by the concrete classes representing norm-conserving pseudopotentials.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def nlcc_radius(self):
        """
        Radius at which the core charge vanish (i.e. cut-off in a.u.).
        Returns 0.0 if nlcc is not used.
        """

    @property
    def has_nlcc(self):
        """True if the pseudo is generated with non-linear core correction."""
        return self.nlcc_radius > 0.0


class PawPseudo(object):
    """
    Abstract class that defines the methods that must be implemented
    by the concrete classes representing PAW pseudopotentials.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def paw_radius(self):
        """Radius of the PAW sphere in a.u."""


class AbinitPseudo(Pseudo):
    """
    An AbinitPseudo is a pseudopotential whose file contains an abinit header.
    """
    def __init__(self, path, header):
        """
        Args:
            path:
                Filename.
            header:
                `AbinitHeader` instance.
        """
        self.path = path
        self._summary = header.summary
        if hasattr(self, "dojo_report"):
            self.dojo_report = header.dojo_report
        else:
            self.dojo_report = {}

        #self.pspcod  = header.pspcod

        for (attr_name, desc) in header.items():
            value = header.get(attr_name, None)

            # Hide these attributes since one should always use the public interface.
            setattr(self, "_" + attr_name, value)

    @property
    def summary(self):
        """Summary line reported in the ABINIT header."""
        return self._summary.strip()

    @property
    def Z(self):
        return self._zatom

    @property
    def Z_val(self):
        return self._zion

    @property
    def l_max(self):
        return self._lmax

    @property
    def l_local(self):
        return self._lloc


class NcAbinitPseudo(NcPseudo, AbinitPseudo):
    """
    Norm-conserving pseudopotential in the Abinit format.
    """
    @property
    def summary(self):
        return self._summary.strip()

    @property
    def Z(self):
        return self._zatom

    @property
    def Z_val(self):
        """Number of valence electrons."""
        return self._zion

    @property
    def l_max(self):
        return self._lmax

    @property
    def l_local(self):
        return self._lloc

    @property
    def nlcc_radius(self):
        return self._rchrg


class PawAbinitPseudo(PawPseudo, AbinitPseudo):
    """Paw pseudopotential in the Abinit format."""

    @property
    def paw_radius(self):
        return self._r_cut

    #def orbitals(self):


class Hint(collections.namedtuple("Hint", "ecut aug_ratio")):
    """
    Suggested value for the cutoff energy [Hartree units] and the augmentation ratio (PAW pseudo)
    """
    @property
    def to_dict(self):
        return {f: getattr(self, f) for f in self._fields}

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k,v in d.items() if not k.startswith("@")})


def _dict_from_lines(lines, key_nums, sep=None):
    """
    Helper function to parse formatted text structured like:

    value1 value2 ... sep key1, key2 ...

    key_nums is a list giving the number of keys for each line. 0 if line should be skipped.
    sep is a string denoting the character that separates the keys from the value (None if
    no separator is present).

    Returns:
        dict{key1 : value1, key2 : value2, ...}

    Raises:
        ValueError if parsing fails.
    """
    if is_string(lines):
        lines = [lines]

    if not isinstance(key_nums, collections.Iterable):
        key_nums = list(key_nums)

    if len(lines) != len(key_nums):
        err_msg = "lines = %s\n key_num =  %s" % (str(lines), str(key_nums))
        raise ValueError(err_msg)

    kwargs = FrozenDict()

    for (i, nk) in enumerate(key_nums):
        if nk == 0: continue
        line = lines[i]

        tokens = [t.strip() for t in line.split()]
        values, keys = tokens[:nk], "".join(tokens[nk:])
        # Sanitize keys: In some case we might string in for  foo[,bar]
        keys.replace("[", "").replace("]", "")
        keys = keys.split(",")

        if sep is not None:
            check = keys[0][0]
            if check != sep:
                raise ValueError("Expecting separator %s, got %s" % (sep, check))
            keys[0] = keys[0][1:]

        if len(values) != len(keys):
            msg = "line: %s\n len(keys) != len(value)\nkeys: %s\n values:  %s" % (line, keys, values)
            warnings.warn(msg)
            #raise ValueError(msg)

        kwargs.update(zip(keys, values))

    return kwargs


class AbinitHeader(dict):
    """Dictionary whose keys can be also accessed as attributes."""
    def __getattr__(self, name):
        try:
            # Default behaviour
            return super(AbinitHeader, self).__getattribute__(name)
        except AttributeError:
            try:
                # Try in the dictionary.
                return self[name]
            except KeyError as exc:
                raise AttributeError(str(exc))


def _int_from_str(string):
    """
    Convert string into integer

    Raise:
        TypeError if string is not a valid integer
    """
    float_num = float(string)
    int_num = int(float_num)
    if float_num == int_num:
        return int_num
    else:
        raise TypeError("Cannot convert string %s to int" % string)

class NcAbinitHeader(AbinitHeader):
    """
    The abinit header found in the NC pseudopotential files.
    """
    _attr_desc = collections.namedtuple("att", "default astype")

    _VARS = {
        # Mandatory
        "zatom"        : _attr_desc(None, _int_from_str),
        "zion"         : _attr_desc(None, float),
        "pspdat"       : _attr_desc(None, float),
        "pspcod"       : _attr_desc(None, int),
        "pspxc"        : _attr_desc(None, int),
        "lmax"         : _attr_desc(None, int),
        "lloc"         : _attr_desc(None, int),
        "r2well"       : _attr_desc(None, float),
        "mmax"         : _attr_desc(None, float),
        # Optional variables for non linear-core correction. HGH does not have it.
        "rchrg"        : _attr_desc(0.0,  float), # radius at which the core charge vanish (i.e. cut-off in a.u.)
        "fchrg"        : _attr_desc(0.0,  float),
        "qchrg"        : _attr_desc(0.0,  float),
    }
    del _attr_desc

    def __init__(self, summary, **kwargs):
        super(NcAbinitHeader, self).__init__()

        # APE uses llocal instead of lloc.
        if "llocal" in kwargs:
            kwargs["lloc"] = kwargs.pop("llocal")

        self.summary = summary.strip()

        for (key, desc) in NcAbinitHeader._VARS.items():
            default, astype = desc.default, desc.astype

            value = kwargs.pop(key, None)

            if value is None:
                value = default
                if default is None:
                    raise RuntimeError("Attribute %s must be specified" % key)
            else:
                try:
                    value = astype(value)
                except:
                    raise RuntimeError("Conversion Error for key, value %s" % (key, value))

            self[key] = value

        # Add dojo_report
        self["dojo_report"] = kwargs.pop("dojo_report", {})

        #if kwargs:
        #    raise RuntimeError("kwargs should be empty but got %s" % str(kwargs))

    @staticmethod
    def fhi_header(filename, ppdesc):
        """Parse the FHI abinit header."""
        # Example:
        # Troullier-Martins psp for element  Sc        Thu Oct 27 17:33:22 EDT 1994
        #  21.00000   3.00000    940714                zatom, zion, pspdat
        #    1    1    2    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
        # 1.80626423934776     .22824404341771    1.17378968127746   rchrg,fchrg,qchrg
        lines = _read_nlines(filename, -1)

        try:
            header = _dict_from_lines(lines[:4], [0, 3, 6, 3])
        except ValueError:
            # The last record with rchrg ... seems to be optional.
            header = _dict_from_lines(lines[:3], [0, 3, 6])

        summary = lines[0]

        header["dojo_report"] = read_dojo_report(filename)

        #print(header)
        return NcAbinitHeader(summary, **header)

    @staticmethod
    def hgh_header(filename, ppdesc):
        """Parse the HGH abinit header."""
        # Example:
        #Hartwigsen-Goedecker-Hutter psp for Ne,  from PRB58, 3641 (1998)
        #   10   8  010605 zatom,zion,pspdat
        # 3 1   1 0 2001 0  pspcod,pspxc,lmax,lloc,mmax,r2well
        lines = _read_nlines(filename, -1)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        header["dojo_report"] = read_dojo_report(filename)

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def tm_header(filename, ppdesc):
        """Parse the TM abinit header."""
        # Example:
        #Troullier-Martins psp for element Fm         Thu Oct 27 17:28:39 EDT 1994
        #100.00000  14.00000    940714                zatom, zion, pspdat
        #   1    1    3    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
        #   0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
        #   .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        #   1   3.116   4.632    1   3.4291849        l,e99.0,e99.9,nproj,rcpsp
        #   .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        #   2   4.557   6.308    1   2.1865358        l,e99.0,e99.9,nproj,rcpsp
        #   .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        #   3  23.251  29.387    1   2.4776730        l,e99.0,e99.9,nproj,rcpsp
        #   .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        #   3.62474762267880     .07409391739104    3.07937699839200   rchrg,fchrg,qchrg
        lines = _read_nlines(filename, -1)
        header = []

        for (lineno, line) in enumerate(lines):
            header.append(line)
            if lineno == 2: 
                # Read lmax.
                tokens = line.split()
                pspcod, pspxc, lmax, lloc = map(int, tokens[:4])
                mmax, r2well = map(float, tokens[4:6])
                #if tokens[-1].strip() != "pspcod,pspxc,lmax,lloc,mmax,r2well":
                #    raise RuntimeError("%s: Invalid line\n %s"  % (filename, line))

                lines = lines[3:]
                break

        # TODO
        # Parse the section with the projectors.
        #0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
        #.00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        projectors = collections.OrderedDict()
        for idx in range(2*(lmax+1)):
            line = lines[idx]
            if idx % 2 == 0: proj_info = [line,]
            if idx % 2 == 1:
                proj_info.append(line)
                d = _dict_from_lines(proj_info, [5,4])
                projectors[int(d["l"])] = d

        # Add the last line with info on nlcc.
        header.append(lines[idx+1])
        summary = header[0]

        header = _dict_from_lines(header, [0,3,6,3])

        header["dojo_report"] = read_dojo_report(filename)

        return NcAbinitHeader(summary, **header)


class PawAbinitHeader(AbinitHeader):
    """
    The abinit header found in the PAW pseudopotential files.
    """
    _attr_desc = collections.namedtuple("att", "default astype")

    _VARS = {
        "zatom"           : _attr_desc(None, _int_from_str),
        "zion"            : _attr_desc(None, float),
        "pspdat"          : _attr_desc(None, float),
        "pspcod"          : _attr_desc(None, int),
        "pspxc"           : _attr_desc(None, int),
        "lmax"            : _attr_desc(None, int),
        "lloc"            : _attr_desc(None, int),
        "mmax"            : _attr_desc(None, int),
        "r2well"          : _attr_desc(None, float),
        "pspfmt"          : _attr_desc(None, str),
        "creatorID"       : _attr_desc(None, int),
        "basis_size"      : _attr_desc(None, int),
        "lmn_size"        : _attr_desc(None, int),
        "orbitals"        : _attr_desc(None, list),
        "number_of_meshes": _attr_desc(None, int),
        "r_cut"           : _attr_desc(None, float), # r_cut(PAW) in the header
        "shape_type"      : _attr_desc(None, int),
        "rshape"          : _attr_desc(None, float),
    }
    del _attr_desc

    def __init__(self, summary, **kwargs):
        super(PawAbinitHeader, self).__init__()

        self.summary = summary.strip()

        for (key, desc) in self._VARS.items():
            default, astype = desc.default, desc.astype

            value = kwargs.pop(key, None)

            if value is None:
                value = default
                if default is None:
                    raise RuntimeError("Attribute %s must be specified" % key)
            else:
                try:
                    value = astype(value)
                except:
                    raise RuntimeError("Conversion Error for key %s, with value %s" % (key, value))

            self[key] = value

        if kwargs:
            raise RuntimeError("kwargs should be empty but got %s" % str(kwargs))

    @staticmethod
    def paw_header(filename, ppdesc):
        """Parse the PAW abinit header."""
        #Paw atomic data for element Ni - Generated by AtomPAW (N. Holzwarth) + AtomPAW2Abinit v3.0.5
        #  28.000  18.000 20061204               : zatom,zion,pspdat
        #  7  7  2 0   350 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well
        # paw3 1305                              : pspfmt,creatorID
        #  5 13                                  : basis_size,lmn_size
        # 0 0 1 1 2                              : orbitals
        # 3                                      : number_of_meshes
        # 1 3  350 1.1803778368E-05 3.5000000000E-02 : mesh 1, type,size,rad_step[,log_step]
        # 2 1  921 2.500000000000E-03                : mesh 2, type,size,rad_step[,log_step]
        # 3 3  391 1.1803778368E-05 3.5000000000E-02 : mesh 3, type,size,rad_step[,log_step]
        #  2.3000000000                          : r_cut(SPH)
        # 2 0.                

        # Example
        #C  (US d-loc) - PAW data extracted from US-psp (D.Vanderbilt) - generated by USpp2Abinit v2.3.0
        #   6.000   4.000 20090106               : zatom,zion,pspdat
        #  7 11  1 0   560 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well
        # paw4 2230                              : pspfmt,creatorID
        #  4  8                                  : basis_size,lmn_size
        # 0 0 1 1                                : orbitals
        # 5                                      : number_of_meshes
        # 1 2  560 1.5198032759E-04 1.6666666667E-02 : mesh 1, type,size,rad_step[,log_step]
        # 2 2  556 1.5198032759E-04 1.6666666667E-02 : mesh 2, type,size,rad_step[,log_step]
        # 3 2  576 1.5198032759E-04 1.6666666667E-02 : mesh 3, type,size,rad_step[,log_step]
        # 4 2  666 1.5198032759E-04 1.6666666667E-02 : mesh 4, type,size,rad_step[,log_step]
        # 5 2  673 1.5198032759E-04 1.6666666667E-02 : mesh 5, type,size,rad_step[,log_step]
        #  1.5550009124                          : r_cut(PAW)
        # 3 0.                                   : shape_type,rshape

        #Paw atomic data for element Si - Generated by atompaw v3.0.1.3 & AtomPAW2Abinit v3.3.1
        #  14.000   4.000 20120814               : zatom,zion,pspdat
        #  7      11  1 0   663 0.               : pspcod,pspxc,lmax,lloc,mmax,r2well
        # paw5 1331                              : pspfmt,creatorID
        #  4  8                                  : basis_size,lmn_size
        # 0 0 1 1                                : orbitals
        # 5                                      : number_of_meshes
        # 1 2  663 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 1, type,size,rad_step[,log_step]
        # 2 2  658 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 2, type,size,rad_step[,log_step]
        # 3 2  740 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 3, type,size,rad_step[,log_step]
        # 4 2  819 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 4, type,size,rad_step[,log_step]
        # 5 2  870 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 5, type,size,rad_step[,log_step]
        #  1.5669671236                          : r_cut(PAW)
        # 2 0.                                   : shape_type,rshape
        supported_formats = ["paw3", "paw4", "paw5"]
        if ppdesc.format not in supported_formats:
            raise NotImplementedError("format %s not in %s" % (ppdesc.format, supported_formats))

        lines = _read_nlines(filename, -1)

        summary = lines[0]
        header = _dict_from_lines(lines[:5], [0, 3, 6, 2, 2], sep=":")

        lines = lines[5:]
        # TODO
        # Parse orbitals and number of meshes.
        header["orbitals"] = [int(t) for t in lines[0].split(":")[0].split()]
        header["number_of_meshes"] = num_meshes = int(lines[1].split(":")[0])

        #print filename, header

        # Skip meshes =
        lines = lines[2+num_meshes:]
        #for midx in range(num_meshes):
        #    l = midx + 1

        #print lines[0]
        header["r_cut"] = float(lines[0].split(":")[0])
        #print lines[1]
        header.update(_dict_from_lines(lines[1], [2], sep=":"))

        report = read_dojo_report(filename)
        if report:
            header["dojo_report"] = report

        #print("PAW header\n", header)
        return PawAbinitHeader(summary, **header)


class PseudoParserError(Exception):
    """Base Error class for the exceptions raised by `PseudoParser`"""


class PseudoParser(object):
    """
    Responsible for parsing pseudopotential files and returning pseudopotential objects.

    Usage::
        pseudo = PseudoParser().parse("filename")
    """
    Error = PseudoParserError

    # Supported values of pspcod
    ppdesc = collections.namedtuple("ppdesc", "pspcod name psp_type format")

    # TODO Recheck
    _PSPCODES = collections.OrderedDict( {
        1 : ppdesc(1, "TM",  "NC", None),
        3 : ppdesc(3, "HGH", "NC", None),
        #4 : ppdesc(4, "NC",     , None),
        #5 : ppdesc(5, "NC",     , None),
        6 : ppdesc(6, "FHI", "NC", None),
        7 : ppdesc(6, "PAW_abinit_text", "PAW", None),
        #8 : ppdesc(8, "NC", None),
       10 : ppdesc(10, "HGHK", "NC", None),
    })
    del ppdesc

    def __init__(self):
        # List of files that have been parsed succesfully.
        self._parsed_paths = []

        # List of files that could not been parsed.
        self._wrong_paths  = []

    def scan_directory(self, dirname, exclude_exts=(), exclude_fnames=()):
        """
        Analyze the files contained in directory dirname.

        Args:
            dirname:
                directory path
            exclude_exts:
                list of file extensions that should be skipped.
            exclude_fnames:
                list of file names that should be skipped.

        returns: 
            List of pseudopotential objects.
        """
        for (i, ext) in enumerate(exclude_exts):
            if not ext.strip().startswith("."):
                exclude_exts[i] =  "." + ext.strip()

        # Exclude files depending on the extension.
        paths = []
        for fname in os.listdir(dirname):
            root, ext = os.path.splitext(fname)
            if ext in exclude_exts or fname in exclude_fnames or fname.startswith("."):
                continue
            paths.append(os.path.join(dirname, fname))

        pseudos = []
        for path in paths:
            # parse the file and generate the pseudo.
            pseudo = self.parse(path)

            if pseudo is not None:
                pseudos.append(pseudo)
                self._parsed_paths.extend(path)
            else:
                self._wrong_paths.extend(path)

        return pseudos

    def read_ppdesc(self, filename):
        """
        Read the pseudopotential descriptor from file filename.

        Returns:
            Pseudopotential descriptor. None if filename is not a valid pseudopotential file.

        Raises:
            `PseudoParserError` if fileformat is not supported.
        """
        if filename.endswith(".xml"):
            raise self.Error("XML pseudo not supported yet")

        else:
            # Assume file with the abinit header.
            lines = _read_nlines(filename, -1)

            for (lineno, line) in enumerate(lines):

                if lineno == 2:
                    try:
                        tokens = line.split()
                        pspcod, pspxc = map(int, tokens[:2])
                    except:
                        msg = "%s: Cannot parse pspcod, pspxc in line\n %s" % (filename, line)
                        sys.stderr.write(msg)
                        return None

                    #if tokens[-1].strip().replace(" ","") not in ["pspcod,pspxc,lmax,lloc,mmax,r2well",
                    #                              "pspcod,pspxc,lmax,llocal,mmax,r2well"]:
                    #    raise self.Error("%s: Invalid line\n %s"  % (filename, line))
                    #    return None

                    if pspcod not in self._PSPCODES:
                        raise self.Error("%s: Don't know how to handle pspcod %s\n"  % (filename, pspcod))

                    ppdesc = self._PSPCODES[pspcod]

                    if pspcod == 7:
                        # PAW -> need to know the format pspfmt
                        tokens = lines[lineno+1].split()
                        pspfmt, creatorID = tokens[:2]
                        #if tokens[-1].strip() != "pspfmt,creatorID":
                        #    raise self.Error("%s: Invalid line\n %s" % (filename, line))
                        #    return None

                        ppdesc = ppdesc._replace(format = pspfmt)

                    return ppdesc

            return None

    def parse(self, filename):
        """
        Read and parse a pseudopotential file. Main entry point for client code.

        Returns: 
            pseudopotential object or None if filename is not a valid pseudopotential file.
        """
        path = os.path.abspath(filename)

        # Only PAW supports XML at present.
        if filename.endswith(".xml"):
            return PawXmlSetup(path)

        ppdesc = self.read_ppdesc(path)

        if ppdesc is None: 
            return None

        psp_type = ppdesc.psp_type

        parsers = {
            "FHI"            : NcAbinitHeader.fhi_header,
            "TM"             : NcAbinitHeader.tm_header,
            "HGH"            : NcAbinitHeader.hgh_header,
            "HGHK"           : NcAbinitHeader.hgh_header,
            "PAW_abinit_text": PawAbinitHeader.paw_header,
        }

        try:
            header = parsers[ppdesc.name](path, ppdesc)

        except Exception as exc:
            raise self.Error(path + ":\n" + straceback())

        root, ext = os.path.splitext(path)

        # Add the content of input file (if present).
        # The name of the input is name + ".ini"
        #input = None
        #input_path = root + ".ini"
        #if os.path.exists(input_path):
        #    with open(input_path, 'r') as fh:
        #        input = fh.read()

        if psp_type == "NC":
            pseudo = NcAbinitPseudo(path, header)
        elif psp_type == "PAW":
            pseudo = PawAbinitPseudo(path, header)
        else:
            raise NotImplementedError("psp_type not in [NC, PAW]")

        return pseudo


#TODO use RadialFunction from pseudo_dojo.
class RadialFunction(collections.namedtuple("RadialFunction", "mesh values")):
    pass


class PawXmlSetup(Pseudo, PawPseudo):
    def __init__(self, filepath):
        # FIXME
        self.dojo_report = {}
        self.path = os.path.abspath(filepath)

        # Get the XML root (this trick is used to that the object is pickleable).
        root = self.root

        # Get the version of the XML format
        self.paw_setup_version = root.get("version")

        # Info on the atom.
        atom_attrib = root.find("atom").attrib

        #self._symbol = atom_attrib["symbol"]
        self._zatom = int(float(atom_attrib["Z"]))
        self.core, self.valence = map(float, [atom_attrib["core"], atom_attrib["valence"]])

        #xc_info = root.find("atom").attrib
        #self.xc_type, self.xc_name  = xc_info["type"], xc_info["name"]
        #self.ae_energy = {k: float(v) for k,v in root.find("ae_energy").attrib.items()}

        self._paw_radius = float(root.find("PAW_radius").attrib["rpaw"])

        #<valence_states>
        #  <state n="2" l="0" f="2"  rc="1.10" e="-0.6766" id="N-2s"/>
        #  <state n="2" l="1" f="3"  rc="1.10" e="-0.2660" id="N-2p"/>
        #  <state       l="0"        rc="1.10" e=" 0.3234" id="N-s1"/>
        #  <state       l="1"        rc="1.10" e=" 0.7340" id="N-p1"/>
        #  <state       l="2"        rc="1.10" e=" 0.0000" id="N-d1"/>
        #</valence_states>
        #
        # The valence_states element contains several state elements.
        # For this setup, the first two lines describe bound eigenstates
        # with occupation numbers and principal quantum numbers.
        # Notice, that the three additional unbound states should have no f and n attributes.
        # In this way, we know that only the first two bound states (with f and n attributes)
        # should be used for constructing an initial guess for the wave functions.

        self.valence_states = {}
        for node in root.find("valence_states"):
            attrib = AttrDict(node.attrib)
            assert attrib.id not in self.valence_states
            self.valence_states[attrib.id] = attrib
        #print(self.valence_states)

        # Parse the radial grids
        self.rad_grids = {}
        for node in root.findall("radial_grid"):
            grid_params = node.attrib
            id = grid_params["id"]
            assert id not in self.rad_grids

            self.rad_grids[id] = self._eval_grid(grid_params)

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.
                                                                                      
        In this case we just remove the XML root element process since Element object cannot be pickled.
        """
        return {k:v for k,v in self.__dict__.items() if k not in ["_root",]}

    @property
    def root(self):
        try:
            return self._root
        except AttributeError:
            from xml.etree import cElementTree  as ET
            tree = ET.parse(self.filepath)
            self._root = tree.getroot()
            return self._root

    @property
    def Z(self):
        return self._zatom

    @property
    def Z_val(self):
        """Number of valence electrons."""
        return self.valence

    # FIXME
    @property
    def l_max(self):
        """Maximum angular momentum."""
        return None
                                                        
    @property
    def l_local(self):
        """Angular momentum used for the local part."""
        return None

    @property
    def summary(self):
        """String summarizing the most important properties."""
        return ""

    @property
    def paw_radius(self):
        return self._paw_radius

    @staticmethod
    def _eval_grid(grid_params):
        """
        This function receives a dictionary with the parameters defining the
        radial mesh and returns a `ndarray` with the mesh
        """
        eq = grid_params.get("eq").replace(" " ,"")
        istart, iend = int(grid_params.get("istart")), int(grid_params.get("iend"))
        indices = range(istart, iend+1)

        if eq == 'r=a*exp(d*i)':
            a, d = float(grid_params['a']), float(grid_params['d'])
            mesh = [a * np.exp(d * i) for i in indices]

        elif eq == 'r=a*i/(n-i)':
            a, n = float(grid_params['a']), float(grid_params['n'])
            mesh = [a * i / (n - i) for i in indices]

        elif eq == 'r=a*(exp(d*i)-1)':
            a, d = float(grid_params['a']), float(grid_params['d'])
            mesh = [a * (np.exp(d * i) - 1.0) for i in indices]

        elif eq == 'r=d*i':
            d = float(grid_params['d'])
            mesh = [d * i for i in indices]

        elif eq == 'r=(i/n+a)^5/a-a^4':
            a, n = float(grid_params['a']), float(grid_params['n'])
            mesh = [(i / n + a)**5 / a - a**4 for i in indices]

        else:
            raise ValueError('Unknown grid type: %s' % eq)

        return np.array(mesh)

    def _parse_radfunc(self, func_name):
        """Parse the first occurence of func_name in the XML file."""
        node = self.root.find(func_name)
        grid = node.attrib["grid"]
        values = np.array([float(s) for s in node.text.split()])

        return self.rad_grids[grid], values, node.attrib

    def _parse_all_radfuncs(self, func_name):
        """Parse all the nodes with tag func_name in the XML file."""
        for node in self.root.findall(func_name):
            grid = node.attrib["grid"]
            values = np.array([float(s) for s in node.text.split()])

            yield self.rad_grids[grid], values, node.attrib

    @property
    def ae_core_density(self):
        """The all-electron radial density."""
        try:
            return self._ae_core_density

        except AttributeError:
            mesh, values, attrib = self._parse_radfunc("ae_core_density")
            self._ae_core_density = RadialFunction(mesh, values)
            return self._ae_core_density

    @property
    def pseudo_core_density(self):
        """The pseudized radial density."""
        try:
            return self._pseudo_core_density

        except AttributeError:
            mesh, values, attrib = self._parse_radfunc("pseudo_core_density")
            self._pseudo_core_density = RadialFunction(mesh, values)
            return self._pseudo_core_density

    @property
    def ae_partial_waves(self):
        """Dictionary with the AE partial waves indexed by state."""
        try:
            return self._ae_partial_waves

        except AttributeError:
            self._ae_partial_waves = {}
            for (mesh, values, attrib) in self._parse_all_radfuncs("ae_partial_wave"):
                state = attrib["state"]
                val_state = self.valence_states[state]
                self._ae_partial_waves[state] = RadialFunction(mesh, values)
                #print("val_state", val_state)

            return self._ae_partial_waves

    @property
    def pseudo_partial_waves(self):
        """Dictionary with the pseudo partial waves indexed by state."""
        try:
            return self._pseudo_partial_waves

        except AttributeError:
            self._pseudo_partial_waves = {}
            for (mesh, values, attrib) in self._parse_all_radfuncs("pseudo_partial_wave"):
                state = attrib["state"]
                val_state = self.valence_states[state]
                self._pseudo_partial_waves[state] = RadialFunction(mesh, values)

            return self._pseudo_partial_waves

    @property
    def projector_functions(self):
        """Dictionary with the PAW projectors indexed by state."""
        try:
            return self._projector_functions

        except AttributeError:
            self._projector_functions = {}
            for (mesh, values, attrib) in self._parse_all_radfuncs("projector_function"):
                state = attrib["state"]
                val_state = self.valence_states[state]
                self._projector_functions[state] = RadialFunction(mesh, values)

            return self._projector_functions

    def plot_densities(self, **kwargs):
        """
        Plot the PAW densities.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        title             Title of the plot (Default: "Densities").
        show              True to show the figure (Default).
        savefig           'abc.png' or 'abc.eps' to save the figure to a file.
        ================  ==============================================================

        Returns:
            `matplotlib` figure
        """
        title = kwargs.pop("title", "Densities")
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)
        ax.grid(True)
        ax.set_xlabel('r [Bohr]')
        #ax.set_ylabel('density')

        for i, den_name in enumerate(["ae_core_density", "pseudo_core_density"]):
            rden = getattr(self, den_name)
            label = "$n_c$" if i == 1 else "$\\tilde{n}_c$"
            ax.plot(rden.mesh, rden.mesh * rden.values, label=label, lw=2)

        plt.legend(loc="best")

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig:
            fig.savefig(savefig)

        return fig

    def plot_waves(self, **kwargs):
        """
        Plot the AE and the pseudo partial waves.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        title             Title of the plot (Default: "Partial Waves").
        show              True to show the figure (Default).
        savefig           'abc.png' or 'abc.eps' to save the figure to a file.
        ================  ==============================================================

        Returns:
            `matplotlib` figure
        """
        title = kwargs.pop("title", "Partial Waves")
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)
        ax.grid(True)
        ax.set_xlabel("r [Bohr]")
        ax.set_ylabel("$r\phi,\\, r\\tilde\phi\, [Bohr]^{-\\frac{1}{2}}$")

        ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        #ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.pseudo_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label="PS-WAVE: " + state)

        for state, rfunc in self.ae_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label="AE-WAVE: " + state)

        plt.legend(loc="best")

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig:
            fig.savefig(savefig)

        return fig

    def plot_projectors(self, **kwargs):
        """
        Plot the PAW projectors.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        title             Title of the plot (Default: "Projectors").
        show              True to show the figure (Default).
        savefig           'abc.png' or 'abc.eps' to save the figure to a file.
        ================  ==============================================================

        Returns:
            `matplotlib` figure
        """
        title = kwargs.pop("title", "Projectors")
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax = fig.add_subplot(1,1,1)
        ax.grid(True)
        ax.set_xlabel('r [Bohr]')
        ax.set_ylabel("$r\\tilde p\, [Bohr]^{-\\frac{1}{2}}$")

        ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        #ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.projector_functions.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, label="TPROJ: " + state)

        plt.legend(loc="best")

        if title is not None:
            fig.suptitle(title)

        if show:
            plt.show()

        if savefig:
            fig.savefig(savefig)

        return fig

    #def plot_potentials(self, **kwargs):
    #    """
    #        ================  ==============================================================
    #        kwargs            Meaning
    #        ================  ==============================================================
    #        title             Title of the plot (Default: None).
    #        show              True to show the figure (Default).
    #        savefig           'abc.png' or 'abc.eps' to save the figure to a file.
    #        ================  ==============================================================

    #    Returns:
    #        `matplotlib` figure
    #    """
    #    title = kwargs.pop("title", "Potentials")
    #    show = kwargs.pop("show", True)
    #    savefig = kwargs.pop("savefig", None)

    #    import matplotlib.pyplot as plt

    #    fig = plt.figure()

    #    ax = fig.add_subplot(1,1,1)
    #    ax.grid(True)
    #    ax.set_xlabel('r [Bohr]')
    #    ax.set_ylabel('density')
    #    ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
    #    ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

    #    for state, rfunc in self.potentials.items():
    #        ax.plot(rfunc.mesh, rfunc.values, label="TPROJ: " + state)

    #    plt.legend(loc="best")

    #    if title is not None:
    #        fig.suptitle(title)

    #    if show:
    #        plt.show()

    #    if savefig:
    #        fig.savefig(savefig)

    #    return fig


class PseudoTable(collections.Sequence):
    """
    Define the pseudopotentials from the element table.
    Individidual elements are accessed by name, symbol or atomic number.

    For example, the following all retrieve iron:

    print elements[26]
    Fe
    print elements.Fe
    Fe
    print elements.symbol('Fe')
    Fe
    print elements.name('iron')
    Fe
    print elements.isotope('Fe')
    Fe
    """
    @classmethod
    def astable(cls, items):
        """
        Return an instance of `PseudoTable` from the iterable items.
        """ 
        if isinstance(items, cls): return items
        return cls(items)

    def __init__(self, pseudos):
        """
        Args:
            pseudos:
                List of pseudopotentials or filepaths
        """
        # Store pseudos in a default dictionary with z as key.
        # Note that we can have more than one pseudo for given z.
        # hence the values are lists of pseudos.
        if not isinstance(pseudos, collections.Iterable):
            pseudos = [pseudos]

        if is_string(pseudos[0]):
            pseudos = list_strings(pseudos)

        self._pseudos_with_z = collections.defaultdict(list)

        for pseudo in pseudos:
            p = pseudo
            if not isinstance(pseudo, Pseudo):
                p = Pseudo.from_file(pseudo)

            self._pseudos_with_z[p.Z].append(p)

        for z in self.zlist:
            pseudo_list = self._pseudos_with_z[z]
            symbols = [p.symbol for p in pseudo_list]
            symbol = symbols[0]
            if any(symb != symbol for symb in symbols):
                raise ValueError("All symbols must be equal while they are: %s" % str(symbols))

            setattr(self, symbol, pseudo_list)

    def __getitem__(self, Z):
        """
        Retrieve pseudos for the atomic number z.
        Accepts both int and slice objects.
        """
        if isinstance(Z, slice):
            assert Z.stop is not None
            pseudos = []
            for znum in iterator_from_slice(Z):
                pseudos.extend(self._pseudos_with_z[znum])
            return pseudos
        else:
            return self._pseudos_with_z[Z]

    def __len__(self):
        return len(list(self.__iter__()))

    def __iter__(self):
        """Process the elements in Z order."""
        for z in self.zlist:
            for pseudo in self._pseudos_with_z[z]:
                yield pseudo

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, id(self))

    def __str__(self):
        lines = []
        app = lines.append
        app("<%s, len=%d>" % (self.__class__.__name__, len(self)))

        for pseudo in self:
            app(str(pseudo))

        return "\n".join(lines)

    @property
    def allnc(self):
        """True if all pseudos are norm-conserving."""
        return all(p.isnc for p in self)

    @property
    def allpaw(self):
        """True if all pseudos are PAW."""
        return all(p.ispaw for p in self)

    @property
    def zlist(self):
        """Ordered list with the atomic numbers available in the table."""
        zlist = list(self._pseudos_with_z.keys())
        zlist.sort()
        return zlist

    def iscomplete(self, zmax=118):
        """
        True if table is complete i.e. all elements with Z < zmax
        have at least on pseudopotential
        """
        for z in range(1, zmax):
            if not self[z]: return False
        return True

    def pseudos_with_symbol(self, symbol):
        """
        Return the list of pseudopotentials in the table the with given symbol.
        Return an empty list if no pseudo is avaiable
        """
        try:
            return getattr(self, str(symbol))
        except AttributeError:
            #raise
            return []

    def pseudo_from_name(self, name):
        """Return the pseudo in the table with the given name"""
        for pseudo in self:
            if pseudo.name == name:
                return pseudo
        return None

    def list_properties(self, *props, **kw):
        """
        Print a list of elements with the given set of properties.

        Args:
            *prop1*, *prop2*, ... : string
                Name of the properties to print
            *format*: string
                Template for displaying the element properties, with one
                % for each property.

        For example, print a table of mass and density.

        from periodictable import elements
        elements.list_properties('symbol','mass','density', format="%-2s: %6.2f u %5.2f g/cm^3")
        H :   1.01 u   0.07 g/cm^3
        He:   4.00 u   0.12 g/cm^3
        Li:   6.94 u   0.53 g/cm^3
        ...
        Bk: 247.00 u  14.00 g/cm^3
        """
        format = kw.pop('format',None)
        assert len(kw) == 0

        for pseudo in self:
            try:
                values = tuple(getattr(pseudo, p) for p in props)
            except AttributeError:
                # Skip elements which don't define all the attributes
                continue

            # Skip elements with a value of None
            if any(v is None for v in values):
                continue

            if format is None:
                print(" ".join(str(p) for p in values))
            else:
                try:
                    print(format % values)
                except:
                    print("format",format,"args",values)
                    raise

    #def print_table(self, stream=sys.stdout, filter_function=None):
    #    """
    #    A pretty ASCII printer for the periodic table, based on some filter_function.

    #    Args:
    #        filter_function:
    #            A filtering function that take a Pseudo as input and returns a boolean.
    #            For example, setting filter_function = lambda el: el.Z_val > 2 will print
    #            a periodic table containing only pseudos with Z_val > 2.
    #    """
    #    for row in range(1, 10):
    #        rowstr = []
    #        for group in range(1, 19):
    #            el = Element.from_row_and_group(row, group)
    #            if el and ((not filter_function) or filter_function(el)):
    #                rowstr.append("{:3s}".format(el.symbol))
    #            else:
    #                rowstr.append("   ")
    #        print(" ".join(rowstr))

    def sorted(self, attrname, reverse=False):
        """Sort the table according to the value of attribute attrname."""
        attrs = []
        for i, pseudo in self:
            try:
                a = getattr(pseudo, attrname)
            except AttributeError:
                a = np.inf
            attrs.append((i, a))

        # Sort attrs, and build new table with sorted pseudos.
        attrs = sorted(attrs, key=lambda t:t[1], reverse=reverse)
        return PseudoTable([self[a[0]] for a in attrs])

    def select(self, condition):
        """
        Select only those pseudopotentials for which condition is True.

        Args:
            condition:
                Function that accepts a `Pseudo` object and returns True or False.
        """
        return PseudoTable([p for p in self if condition(p)])

    def with_dojo_report(self):
        """Select pseudos containing the DOJO_REPORT section."""
        return self.select(condition=lambda p : p.has_dojo_report)
