# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides objects describing the basic parameters of the 
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..
"""
from __future__ import unicode_literals, division, print_function

import sys
import os
import abc
import collections
import json
import six
#import pprint
import numpy as np

from warnings import warn
from collections import OrderedDict, defaultdict, namedtuple
from monty.string import list_strings, is_string
from monty.itertools import iterator_from_slice
from monty.io import FileLock
from monty.collections import AttrDict, Namespace
from monty.functools import lazy_property
from monty.os.path import find_exts
from monty.dev import deprecated
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt
from pymatgen.core.periodic_table import PeriodicTable, Element
from pymatgen.serializers.json_coders import PMGSONable, pmg_serialize
from .eos import EOS
from monty.json import MontyDecoder

import logging
logger = logging.getLogger(__name__)


__all__ = [
    "Pseudo",
    "PseudoTable",
]

__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


_PTABLE = PeriodicTable()

# Tools and helper functions.


def straceback():
    """Returns a string with the traceback."""
    import sys
    import traceback
    return "\n".join((traceback.format_exc(), str(sys.exc_info()[0])))


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


class Pseudo(six.with_metaclass(abc.ABCMeta, PMGSONable, object)):
    """
    Abstract base class defining the methods that must be 
    implemented by the concrete pseudopotential classes.
    """

    @classmethod
    def as_pseudo(cls, obj):
        """
        Convert obj into a pseudo. Accepts:

            * Pseudo object.
            * string defining a valid path.
        """
        return obj if isinstance(obj, cls) else cls.from_file(obj)

    @staticmethod
    def from_file(filename):
        """
        Return a :class:`Pseudo` object from filename.
        Note: the parser knows the concrete class that should be instanciated
        """
        return PseudoParser().parse(filename)

    def __eq__(self, other):
        if other is None: return False
        # TODO
        # For the time being we check the filepath
        # A more robust algorithm would use md5
        #return self.filepath == other.filepath
        return (self.md5 == other.md5 and
                self.__class__ == other.__class__ and 
                self.Z == other.Z and 
                self.Z_val == other.Z_val and
                self.l_max == other.l_max )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "<%s at %s>" % (self.__class__.__name__, self.filepath)

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        app("<%s: %s>" % (self.__class__.__name__, self.basename))
        app("  summary: " + self.summary.strip())
        app("  number of valence electrons: %s" % self.Z_val)
        #FIXME: rewrite the treatment of xc, use XML specs as starting point
        #app("  XC correlation (ixc): %s" % self._pspxc)  #FIXME
        app("  maximum angular momentum: %s" % l2str(self.l_max))
        app("  angular momentum for local part: %s" % l2str(self.l_local))
        if self.isnc:
            app("  radius for non-linear core correction: %s" % self.nlcc_radius)
        app("")

        if self.has_hints:
            hint_normal = self.hint_for_accuracy()
            if hint_normal is not None:
                app("  hint for normal accuracy: %s" % str(hint_normal))
        else:
                app("  hints on cutoff-energy are not available")

        return "\n".join(lines)

    @abc.abstractproperty
    def summary(self):
        """String summarizing the most important properties."""

    @property
    def filepath(self):
        return os.path.abspath(self.path)

    @property
    def basename(self):
        """File basename."""
        return os.path.basename(self.filepath)

    @abc.abstractproperty
    def Z(self):
        """The atomic number of the atom."""

    @abc.abstractproperty
    def Z_val(self):
        """Valence charge"""

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def element(self):
        """Pymatgen :class:`Element`."""
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

    @lazy_property
    def md5(self):
        """MD5 hash value."""
        if self.has_dojo_report:
            if "md5" in self.dojo_report:
                return self.dojo_report["md5"]
            else:
                warn("Dojo report without md5 entry")

        return self.compute_md5()

    def compute_md5(self):
        """Compute MD5 hash value."""
        import hashlib

        if self.path.endswith(".xml"):
            # TODO: XML + DOJO_REPORT
            #raise NotImplementedError("md5 for XML files!")
            with open(self.path, "rt") as fh:
                text = fh.read()

        else:
            # If we have a pseudo with a dojo_report at the end.
            # we compute the hash from the data before DOJO_REPORT.
            # else all the lines are taken.
            with open(self.path, "rt") as fh:
                lines = fh.readlines()
                try:
                    start = lines.index("<DOJO_REPORT>\n")
                except ValueError:
                    start = len(lines)
                text = "".join(lines[:start])

        m = hashlib.md5(text.encode("utf-8"))
        return m.hexdigest()

    def check_and_fix_dojo_md5(self):
        report = self.read_dojo_report()

        if "md5" in report:
            if report["md5"] != self.md5:
                raise ValueError("md5 found in dojo_report does not agree\n"
                    "with the computed value\nFound %s\nComputed %s"  % (report["md5"], hash))
        else:
            report["md5"] = self.compute_md5()
            self.write_dojo_report(report=report)

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

    @pmg_serialize
    def as_dict(self, **kwargs):
        return dict(
            basename=self.basename,
            type=self.type,
            symbol=self.symbol,
            Z=self.Z,
            Z_val=self.Z_val,
            l_max=self.l_max,
            md5=self.md5,
            filepath=self.filepath
        )

    @classmethod
    def from_dict(cls, d):
        new = cls.from_file(d['filepath'])

        # Consistency test based on md5
        if "md5" in d and d["md5"] != new.md5:
            raise ValueError("The md5 found in file does not agree with the one in dict\n"
            "Received %s\nComputed %s" % (d["md5"], new.md5))

        return new

    def as_tmpfile(self):
        """
        Copy the pseudopotential to a temporary a file and returns a new pseudopotential object.
        """
        import tempfile, shutil
        _, dst = tempfile.mkstemp(suffix=self.basename, text=True)
        shutil.copy(self.path, dst)
        return self.__class__.from_file(dst)

    @property
    def has_dojo_report(self):
        """True if self contains the `DOJO_REPORT` section."""
        return bool(self.dojo_report)

    def delta_factor(self, accuracy="normal"):
        """
        Returns the deltafactor [meV/natom] computed with the given accuracy.
        None if the `Pseudo` does not have info on the deltafactor.
        """
        if not self.has_dojo_report:
            return None
        try:
            return self.dojo_report["delta_factor"][accuracy]["dfact"]
        except KeyError:
            return None

    def read_dojo_report(self):
        """
        Read the `DOJO_REPORT` section and set the `dojo_report` attribute.
        returns {} if section is not present.
        """ 
        self.dojo_report = DojoReport.from_file(self.path)
        return self.dojo_report

    def write_dojo_report(self, report=None):
        """Write a new `DOJO_REPORT` section to the pseudopotential file."""
        if report is None:
            report = self.dojo_report

        report["symbol"] = self.symbol

        if "md5" not in report:
            report["md5"] = self.md5

        if report["md5"] != self.md5:
            raise ValueError("md5 found in dojo_report does not agree\n"
               "with the computed value\nreport: %s\npseudo %s" % (report["md5"], self.md5))

        # Create JSON string from report.
        jstring = json.dumps(report, indent=4, sort_keys=True) + "\n"

        # Read lines from file and insert jstring between the tags.
        with open(self.path, "r") as fh:
            lines = fh.readlines()
            try:
                start = lines.index("<DOJO_REPORT>\n")
            except ValueError:
                start = -1

            if start == -1:
                # DOJO_REPORT was not present.
                lines += ["<DOJO_REPORT>\n", jstring , "</DOJO_REPORT>\n",]
            else:
                stop = lines.index("</DOJO_REPORT>\n")
                lines.insert(stop, jstring)
                del lines[start+1:stop]

        #  Write new file.
        with FileLock(self.path):
            with open(self.path, "w") as fh:
                fh.writelines(lines)

    def remove_dojo_report(self):
        """Remove the `DOJO_REPORT` section from the pseudopotential file."""
        # Read lines from file and insert jstring between the tags.
        with open(self.path, "r") as fh:
            lines = fh.readlines()
            try:
                start = lines.index("<DOJO_REPORT>\n")
            except ValueError:
                start = -1

            if start == -1: return

            stop = lines.index("</DOJO_REPORT>\n")
            if stop == -1: return

            del lines[start+1:stop]

        # Write new file.
        with FileLock(self.path):
            with open(self.path, "w") as fh:
                fh.writelines(lines)

    def hint_for_accuracy(self, accuracy="normal"):
        """
        Returns an hint object with parameters such as ecut [Ha] and 
        aug_ratio for given accuracy. Returns None if no hint is available.

        Args:
            accuracy: ["low", "normal", "high"]
        """
        if self.has_dojo_report:
            try:
                return Hint.from_dict(self.dojo_report["hints"][accuracy])
            except KeyError:
                return None
        else:
            return None

    @property
    def has_hints(self):
        """True if self provides hints on the cutoff energy."""
        for acc in ["low", "normal", "high"]:
            try:
                if self.hint_for_accuracy(acc) is None:
                    return False
            except KeyError:
                return False
        return True

    def open_pspsfile(self, ecut=20, pawecutdg=None):
        """
        Calls Abinit to compute the internal tables for the application of the 
        pseudopotential part. Returns :class:`PspsFile` object providing methods
        to plot and analyze the data or None if file is not found or it's not readable.

        Args:
            ecut: Cutoff energy in Hartree.
            pawecutdg: Cutoff energy for the PAW double grid.
        """
        from pymatgen.io.abinit.tasks import AbinitTask
        from abipy.core.structure import Structure
        from abipy.abio.factories import gs_input
        from abipy.electrons.psps import PspsFile

        # Build fake structure.
        lattice = 10 * np.eye(3)
        structure = Structure(lattice, [self.element], coords=[[0, 0, 0]])

        if self.ispaw and pawecutdg is None: pawecudg = ecut * 4
        inp = gs_input(structure, pseudos=[self], ecut=ecut, pawecutdg=pawecutdg, 
                       spin_mode="unpolarized", kppa=1)
        # Add prtpsps = -1 to make Abinit print the PSPS.nc file and stop.
        inp["prtpsps"] = -1

        # Build temporary task and run it.
        task = AbinitTask.temp_shell_task(inp)
        retcode = task.start_and_wait()

        filepath = task.outdir.has_abiext("_PSPS.nc")
        if not filepath:
            logger.critical("Cannot find PSPS.nc file in %s" % task.outdir)
            return None

        # Open the PSPS.nc file.
        try:
            return PspsFile(filepath)
        except Exception as exc:
            logger.critical("Exception while reading PSPS file at %s:\n%s" % (filepath, str(exc)))
            return None


class NcPseudo(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract class defining the methods that must be implemented
    by the concrete classes representing norm-conserving pseudopotentials.
    """

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

    @property
    def rcore(self):
        """Radius of the pseudization sphere in a.u."""
        try:
            return self._core
        except AttributeError:
            return None


class PawPseudo(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract class that defines the methods that must be implemented
    by the concrete classes representing PAW pseudopotentials.
    """

    #def nlcc_radius(self):
    #    """
    #    Radius at which the core charge vanish (i.e. cut-off in a.u.).
    #    Returns 0.0 if nlcc is not used.
    #    """
    #    return 0.0
    #                                                                           

    #@property
    #def has_nlcc(self):
    #    """True if the pseudo is generated with non-linear core correction."""
    #    return True

    @abc.abstractproperty
    def paw_radius(self):
        """Radius of the PAW sphere in a.u."""

    @property
    def rcore(self):
        """Alias of paw_radius."""
        return self.paw_radius


class AbinitPseudo(Pseudo):
    """
    An AbinitPseudo is a pseudopotential whose file contains an abinit header.
    """
    def __init__(self, path, header):
        """
        Args:
            path: Filename.
            header: :class:`AbinitHeader` instance.
        """
        self.path = path
        self._summary = header.summary

        if hasattr(header, "dojo_report"):
            self.dojo_report = header.dojo_report
        else:
            self.dojo_report = {}

        #self.pspcod  = header.pspcod

        for attr_name, desc in header.items():
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
    """Norm-conserving pseudopotential in the Abinit format."""
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


#class Hint(namedtuple("Hint", "ecut aug_ratio")):
class Hint(object):
    """
    Suggested value for the cutoff energy [Hartree units] 
    and the cutoff energy for the dense grid (only for PAW pseudos)
    """
    def __init__(self, ecut, pawecutdg=None):
        self.ecut = ecut
        self.pawecutdg = ecut if pawecutdg is None else pawecutdg
        
    @pmg_serialize
    def as_dict(self):
        return dict(ecut=self.ecut, pawecutdg=self.pawecutdg)

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

    kwargs = Namespace()

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
            raise ValueError(msg)

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
        # Needed to handle pseudos with fractional charge
        int_num = np.rint(float_num)
        warn("Converting float %s to int %s" % (float_num, int_num))
        return int_num
        #raise TypeError("Cannot convert string %s to int" % string)


class NcAbinitHeader(AbinitHeader):
    """The abinit header found in the NC pseudopotential files."""
    _attr_desc = namedtuple("att", "default astype")

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
        "rchrg"        : _attr_desc(0.0,  float),  # radius at which the core charge vanish (i.e. cut-off in a.u.)
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
                    raise RuntimeError("Conversion Error for key %s, value %s" % (key, value))

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

        header["dojo_report"] = DojoReport.from_file(filename)

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

        header["dojo_report"] = DojoReport.from_file(filename)

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def gth_header(filename, ppdesc):
        """Parse the GTH abinit header."""
        # Example:
        #Goedecker-Teter-Hutter  Wed May  8 14:27:44 EDT 1996
        #1   1   960508                     zatom,zion,pspdat
        #2   1   0    0    2001    0.       pspcod,pspxc,lmax,lloc,mmax,r2well
        #0.2000000 -4.0663326  0.6778322 0 0     rloc, c1, c2, c3, c4
        #0 0 0                              rs, h1s, h2s
        #0 0                                rp, h1p
        #  1.36 .2   0.6                    rcutoff, rloc
        lines = _read_nlines(filename, -1)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        header["dojo_report"] = DojoReport.from_file(filename)

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def oncvpsp_header(filename, ppdesc):
        """Parse the ONCVPSP abinit header."""
        # Example
        #Li    ONCVPSP  r_core=  2.01  3.02
        #      3.0000      3.0000      140504    zatom,zion,pspd
        #     8     2     1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well
        #  5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
        #     2     2     0     0     0    nproj
        #     0                 extension_switch
        #   0                        -2.5000025868368D+00 -1.2006906995331D+00
        #     1  0.0000000000000D+00  0.0000000000000D+00  0.0000000000000D+00
        #     2  1.0000000000000D-02  4.4140499497377D-02  1.9909081701712D-02
        lines = _read_nlines(filename, -1)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        header.update({'pspdat': header['pspd']})
        header.pop('pspd')
        try:
            header["dojo_report"] = DojoReport.from_file(filename)
        except DojoReport.Error:
            logger.warning('failed to read the dojo report for %s' % filename)
            header["dojo_report"] = None

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
        projectors = OrderedDict()
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

        header["dojo_report"] = DojoReport.from_file(filename)

        return NcAbinitHeader(summary, **header)


class PawAbinitHeader(AbinitHeader):
    """The abinit header found in the PAW pseudopotential files."""
    _attr_desc = namedtuple("att", "default astype")

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

        report = DojoReport.from_file(filename)
        if report:
            header["dojo_report"] = report

        #print("PAW header\n", header)
        return PawAbinitHeader(summary, **header)


class PseudoParserError(Exception):
    """Base Error class for the exceptions raised by :class:`PseudoParser`"""


class PseudoParser(object):
    """
    Responsible for parsing pseudopotential files and returning pseudopotential objects.

    Usage::

        pseudo = PseudoParser().parse("filename")
    """
    Error = PseudoParserError

    # Supported values of pspcod
    ppdesc = namedtuple("ppdesc", "pspcod name psp_type format")

    # TODO Recheck
    _PSPCODES = OrderedDict( {
        1: ppdesc(1, "TM",  "NC", None),
        2: ppdesc(2, "GTH",  "NC", None),
        3: ppdesc(3, "HGH", "NC", None),
        #4: ppdesc(4, "NC",     , None),
        #5: ppdesc(5, "NC",     , None),
        6: ppdesc(6, "FHI", "NC", None),
        7: ppdesc(6, "PAW_abinit_text", "PAW", None),
        8: ppdesc(8, "ONCVPSP", "NC", None),
       10: ppdesc(10, "HGHK", "NC", None),
    })
    del ppdesc
    # renumber functionals from oncvpsp todo confrim that 3 is 2
    _FUNCTIONALS = {1: {'n': 4, 'name': 'Wigner'},
                    2: {'n': 5, 'name': 'HL'},
                    3: {'n': 2, 'name': 'PWCA'},
                    4: {'n': 11, 'name': 'PBE'}}

    def __init__(self):
        # List of files that have been parsed succesfully.
        self._parsed_paths = []

        # List of files that could not been parsed.
        self._wrong_paths  = []

    def scan_directory(self, dirname, exclude_exts=(), exclude_fnames=()):
        """
        Analyze the files contained in directory dirname.

        Args:
            dirname: directory path
            exclude_exts: list of file extensions that should be skipped.
            exclude_fnames: list of file names that should be skipped.

        Returns:
            List of pseudopotential objects.
        """
        for (i, ext) in enumerate(exclude_exts):
            if not ext.strip().startswith("."):
                exclude_exts[i] =  "." + ext.strip()

        # Exclude files depending on the extension.
        paths = []
        for fname in os.listdir(dirname):
            root, ext = os.path.splitext(fname)
            path = os.path.join(dirname, fname)
            if (ext in exclude_exts or fname in exclude_fnames or 
                fname.startswith(".") or not os.path.isfile(path)): continue
            paths.append(path)

        pseudos = []
        for path in paths:
            # Parse the file and generate the pseudo.
            try:
                pseudo = self.parse(path)
            except:
                pseudo = None

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
            lines = _read_nlines(filename, 80)

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
                        raise self.Error("%s: Don't know how to handle pspcod %s\n" % (filename, pspcod))

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
            "GTH"            : NcAbinitHeader.gth_header,
            "TM"             : NcAbinitHeader.tm_header,
            "HGH"            : NcAbinitHeader.hgh_header,
            "HGHK"           : NcAbinitHeader.hgh_header,
            "ONCVPSP"        : NcAbinitHeader.oncvpsp_header,
            "PAW_abinit_text": PawAbinitHeader.paw_header,
        }

        try:
            header = parsers[ppdesc.name](path, ppdesc)
        except Exception as exc:
            raise self.Error(path + ":\n" + straceback())

        root, ext = os.path.splitext(path)

        if psp_type == "NC":
            pseudo = NcAbinitPseudo(path, header)
        elif psp_type == "PAW":
            pseudo = PawAbinitPseudo(path, header)
        else:
            raise NotImplementedError("psp_type not in [NC, PAW]")

        return pseudo


#TODO use RadialFunction from pseudo_dojo.
class RadialFunction(namedtuple("RadialFunction", "mesh values")):
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

        # Old XML files do not define this field!
        # In this case we set the PAW radius to None.
        #self._paw_radius = float(root.find("PAW_radius").attrib["rpaw"])

        pawr_element = root.find("PAW_radius")
        self._paw_radius = None
        if pawr_element is not None:
            self._paw_radius = float(pawr_element.attrib["rpaw"])

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
            gid = grid_params["id"]
            assert gid not in self.rad_grids

            self.rad_grids[id] = self._eval_grid(grid_params)

    def __getstate__(self):
        """
        Return state is pickled as the contents for the instance.
                                                                                      
        In this case we just remove the XML root element process since Element object cannot be pickled.
        """
        return {k: v for k, v in self.__dict__.items() if k not in ["_root"]}

    @property
    def root(self):
        try:
            return self._root
        except AttributeError:
            from xml.etree import cElementTree as Et
            tree = Et.parse(self.filepath)
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
        eq = grid_params.get("eq").replace(" ", "")
        istart, iend = int(grid_params.get("istart")), int(grid_params.get("iend"))
        indices = list(range(istart, iend+1))

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

    @add_fig_kwargs
    def plot_densities(self, ax=None, **kwargs):
        """
        Plot the PAW densities.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        ax.set_xlabel('r [Bohr]')
        #ax.set_ylabel('density')

        for i, den_name in enumerate(["ae_core_density", "pseudo_core_density"]):
            rden = getattr(self, den_name)
            label = "$n_c$" if i == 1 else "$\\tilde{n}_c$"
            ax.plot(rden.mesh, rden.mesh * rden.values, label=label, lw=2)

        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_waves(self, ax=None, **kwargs):
        """
        Plot the AE and the pseudo partial waves.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        ax.grid(True)
        ax.set_xlabel("r [Bohr]")
        ax.set_ylabel("$r\phi,\\, r\\tilde\phi\, [Bohr]^{-\\frac{1}{2}}$")

        ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        #ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.pseudo_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label="PS-WAVE: " + state)

        for state, rfunc in self.ae_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label="AE-WAVE: " + state)

        ax.legend(loc="best")
        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax=None, **kwargs):
        """
        Plot the PAW projectors.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig, plt = get_ax_fig_plt(ax)
        title = kwargs.pop("title", "Projectors")
        ax.grid(True)
        ax.set_xlabel('r [Bohr]')
        ax.set_ylabel("$r\\tilde p\, [Bohr]^{-\\frac{1}{2}}$")

        ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        #ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.projector_functions.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, label="TPROJ: " + state)

        ax.legend(loc="best")

        return fig

    #@add_fig_kwargs
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

    #    ax.legend(loc="best")

    #    if title is not None: fig.suptitle(title)
    #    if show: plt.show()
    #    if savefig: fig.savefig(savefig)
    #    return fig


class PseudoTable(six.with_metaclass(abc.ABCMeta, collections.Sequence, PMGSONable, object)):
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
    def as_table(cls, items):
        """
        Return an instance of :class:`PseudoTable` from the iterable items.
        """ 
        if isinstance(items, cls): return items
        return cls(items)

    @classmethod
    def from_dir(cls, top, exts=None, exclude_dirs="_*"):
        """
        Find all pseudos in the directory tree starting from top.

        Args:
            top: Top of the directory tree
            exts: List of files extensions. if exts == "all_files"
                    we try to open all files in top
            exclude_dirs: Wildcard used to exclude directories.
        
        return: :class:`PseudoTable` sorted by atomic number Z.
        """
        pseudos = []

        if exts == "all_files":
            for f in [os.path.join(top, fn) for fn in os.listdir(top)]:
                if os.path.isfile(f):
                    try:
                        p = Pseudo.from_file(f)
                        if p:
                            pseudos.append(p)
                        else:
                            logger.info('Skipping file %s' % f)
                    except:
                        logger.info('Skipping file %s' % f)
            if not pseudos:
                logger.warning('No pseudopotentials parsed from folder %s' % top)
                return None
            logger.info('Creating PseudoTable with %i pseudopotentials' % len(pseudos))

        else:
            if exts is None: exts=("psp8",)

            for p in find_exts(top, exts, exclude_dirs=exclude_dirs):
                try:
                    pseudos.append(Pseudo.from_file(p))
                except Exception as exc:
                    logger.critical("Error in %s:\n%s" % (p, exc))

        return cls(pseudos).sort_by_z()

    def __init__(self, pseudos):
        """
        Args:
            pseudos: List of pseudopotentials or filepaths
        """
        # Store pseudos in a default dictionary with z as key.
        # Note that we can have more than one pseudo for given z.
        # hence the values are lists of pseudos.
        if not isinstance(pseudos, collections.Iterable):
            pseudos = [pseudos]

        if len(pseudos) and is_string(pseudos[0]):
            pseudos = list_strings(pseudos)

        self._pseudos_with_z = defaultdict(list)

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
        Retrieve pseudos for the atomic number z. Accepts both int and slice objects.
        """
        if isinstance(Z, slice):
            assert Z.stop is not None
            pseudos = []
            for znum in iterator_from_slice(Z):
                pseudos.extend(self._pseudos_with_z[znum])
            return self.__class__(pseudos)
        else:
            return self.__class__(self._pseudos_with_z[Z])

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
        return sorted(list(self._pseudos_with_z.keys()))

    def as_dict(self, **kwargs):
        d = {}
        for p in self:
            k, count = p.element, 1
            # Handle multiple-pseudos with the same name!
            while k in d:
                k += k.split("#")[0] + "#" + str(count)
                count += 1
            d.update({k: p.as_dict()})
        d['@module'] = self.__class__.__module__
        d['@class'] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        pseudos = []
        dec = MontyDecoder()
        for k, v in d.items():
            if not k.startswith('@'):
                pseudos.append(dec.process_decoded(v))
        return cls(pseudos)

    def is_complete(self, zmax=118):
        """
        True if table is complete i.e. all elements with Z < zmax have at least on pseudopotential
        """
        for z in range(1, zmax):
            if not self[z]: return False
        return True

    def all_combinations_for_elements(self, element_symbols):
        """
        Return a list with all the the possible combination of pseudos 
        for the given list of element_symbols.
        Each item is a list of pseudopotential objects.

        Example::

            table.all_combinations_for_elements(["Li", "F"])
        """
        d = OrderedDict()
        for symbol in element_symbols:
            d[symbol] = self.select_symbols(symbol, ret_list=True)

        from itertools import product
        all = product(*d.values())
        return list(all)

    def pseudo_with_symbol(self, symbol, allow_multi=False):
        """
        Return the pseudo with the given chemical symbol.

        Args:
            symbols: String with the chemical symbol of the element
            allow_multi: By default, the method raises ValueError
                if multiple occurrences are found. Use allow_multi to prevent this.

        Raises:
            ValueError if symbol is not found or multiple occurences are present and not allow_multi
        """
        pseudos = self.select_symbols(symbol, ret_list=True)
        if not pseudos or (len(pseudos) > 1 and not allow_multi):
            raise ValueError("Found %d occurrences of symbol %s" % (len(pseudos), symbol))

        if not allow_multi: 
            return pseudos[0]
        else:
            return pseudos

    def pseudos_with_symbols(self, symbols):
        """
        Return the pseudos with the given chemical symbols.

        Raises:
            ValueError if one of the symbols is not found or multiple occurences are present.
        """
        pseudos = self.select_symbols(symbols, ret_list=True)
        found_symbols = [p.symbol for p in pseudos]
        duplicated_elements = [s for s, o in collections.Counter(found_symbols).items() if o > 1]

        if duplicated_elements:
            raise ValueError("Found multiple occurrences of symbol(s) %s" % ', '.join(duplicated_elements))
        missing_symbols = [s for s in symbols if s not in found_symbols]

        if missing_symbols:
            raise ValueError("Missing data for symbol(s) %s" % ', '.join(missing_symbols))

        return pseudos

    def select_symbols(self, symbols, ret_list=False):
        """
        Return a :class:`PseudoTable` with the pseudopotentials with the given list of chemical symbols.

        Args:
            symbols: str or list of symbols
                Prepend the symbol string with "-", to exclude pseudos.
            ret_list: if True a list of pseudos is returned instead of a :class:`PseudoTable`
        """
        symbols = list_strings(symbols)
        exclude = symbols[0].startswith("-")

        if exclude: 
            if not all(s.startswith("-") for s in symbols):
                raise ValueError("When excluding symbols, all strings must start with `-`")
            symbols = [s[1:] for s in symbols]
            #print(symbols)

        symbols = set(symbols)
        pseudos = []
        for p in self:
            if exclude:
                if p.symbol in symbols: continue
            else:
                if p.symbol not in symbols: continue

            pseudos.append(p)
    
        if ret_list:
            return pseudos
        else:
            return self.__class__(pseudos)

    def get_pseudos_for_structure(self, structure):
        """
        Return the list of :class:`Pseudo` objects to be used for this :class:`Structure`.

        Args:
            structure: pymatgen :class:`Structure`.

        Raises:
            `ValueError` if one of the chemical symbols is not found or 
            multiple occurences are present in the table.
        """
        symbols = structure.symbol_set
        return self.pseudos_with_symbols(symbols)

    #def list_properties(self, *props, **kw):
    #    """
    #    Print a list of elements with the given set of properties.

    #    Args:
    #        *prop1*, *prop2*, ... : string
    #            Name of the properties to print
    #        *format*: string
    #            Template for displaying the element properties, with one
    #            % for each property.

    #    For example, print a table of mass and density.

    #    from periodictable import elements
    #    elements.list_properties('symbol','mass','density', format="%-2s: %6.2f u %5.2f g/cm^3")
    #    H :   1.01 u   0.07 g/cm^3
    #    He:   4.00 u   0.12 g/cm^3
    #    Li:   6.94 u   0.53 g/cm^3
    #    ...
    #    Bk: 247.00 u  14.00 g/cm^3
    #    """
    #    format = kw.pop('format', None)
    #    assert len(kw) == 0

    #    for pseudo in self:
    #        try:
    #            values = tuple(getattr(pseudo, p) for p in props)
    #        except AttributeError:
    #            # Skip elements which don't define all the attributes
    #            continue

    #        # Skip elements with a value of None
    #        if any(v is None for v in values):
    #            continue

    #        if format is None:
    #            print(" ".join(str(p) for p in values))
    #        else:
    #            try:
    #                print(format % values)
    #            except:
    #                print("format",format,"args",values)
    #                raise

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
        """
        Sort the table according to the value of attribute attrname.

        Return:
            New class:`PseudoTable` object
        """
        attrs = []
        for i, pseudo in self:
            try:
                a = getattr(pseudo, attrname)
            except AttributeError:
                a = np.inf
            attrs.append((i, a))

        # Sort attrs, and build new table with sorted pseudos.
        return self.__class__([self[a[0]] for a in sorted(attrs, key=lambda t: t[1], reverse=reverse)])

    def sort_by_z(self):
        """Return a new :class:`PseudoTable` with pseudos sorted by Z"""
        return self.__class__(sorted(self, key=lambda p: p.Z))

    def select(self, condition):
        """
        Select only those pseudopotentials for which condition is True.
        Return new class:`PseudoTable` object.

        Args:
            condition:
                Function that accepts a :class:`Pseudo` object and returns True or False.
        """
        return self.__class__([p for p in self if condition(p)])

    def with_dojo_report(self):
        """Select pseudos containing the DOJO_REPORT section. Return new class:`PseudoTable` object."""
        return self.select(condition=lambda p: p.has_dojo_report)

    def get_dojo_dataframe(self, **kwargs):
        """
        Buid a pandas :class:`DataFrame` with the most important parameters extracted from the 
        `DOJO_REPORT` section of each pseudo in the table.

        Returns:
            frame, errors

            where frame is the pandas :class:`DataFrame` and errors is a list of errors
            encountered while trying to read the `DOJO_REPORT` from the pseudopotential file.
        """
        accuracies = ["low", "normal", "high"]

        trial2keys = {
            "deltafactor": ["dfact_meV", "dfactprime_meV"] + ["v0", "b0_GPa", "b1"], 
            "gbrv_bcc": ["a0_rel_err"],
            "gbrv_fcc": ["a0_rel_err"],
            "phonon": "all",
            #"phwoa": "all"
        }

        rows, names, errors = [], [], []

        for p in self:
            report = p.dojo_report
            #assert "version"  in report
            if "version" not in report:
                print("ignoring old report in ", p.basename)
                continue

            d = {"symbol": p.symbol, "Z": p.Z}
            names.append(p.basename)

            #read hints
            for acc in accuracies:
                try:
                    d.update({acc + "_ecut_hint": report['hints'][acc]['ecut']})
                except KeyError:
                    d.update({acc + "_ecut_hint": -1.0 })

            # FIXME
            try:
                ecut_acc = dict(
                    low=report.ecuts[2],
                    normal=report.ecuts[int(len(report.ecuts)/2)],
                    high=report.ecuts[-2],
                )
            except IndexError:
                ecut_acc = dict(
                    low=report.ecuts[0],
                    normal=report.ecuts[-1],
                    high=report.ecuts[-1],
                )

            for acc in accuracies:
                d[acc + "_ecut"] = ecut_acc[acc]

            try:
                for trial, keys in trial2keys.items():
                    data = report.get(trial, None)
                    if data is None: continue
                    # if the current trial has an entry for this ecut notting changes, else we take the ecut closes 
                    ecut_acc = dict(
                        low=sorted(data.keys())[0],
                        normal=sorted(data.keys())[int(len(data.keys())/2)],
                        high=sorted(data.keys())[-1],
                    )
                    for acc in accuracies:
                        ecut = ecut_acc[acc]
                        if keys is 'all':
                            ecuts = data
                            d.update({acc + "_" + trial: data[ecut]})
                        else:
                            if trial.startswith("gbrv"):
                                d.update({acc + "_" + trial + "_" + k: float(data[ecut][k]) for k in keys}) 
                            else:
                                d.update({acc + "_" + k: float(data[ecut][k]) for k in keys}) 

            except Exception as exc:
                logger.warning("%s raised %s" % (p.basename, exc))
                errors.append((p.basename, str(exc)))

            #print(d)
            rows.append(d)

        # Build sub-class of pandas.DataFrame
        return DojoDataFrame(rows, index=names), errors

    def select_rows(self, rows):
        """
        Return new class:`PseudoTable` object with pseudos in the given rows of the periodic table.
        rows can be either a int or a list of integers.
        """
        if not isinstance(rows, (list, tuple)): rows = [rows]
        return self.__class__([p for p in self if p.element.row in rows])

    def select_family(self, family):
        # e.g element.is_alkaline
        return self.__class__([p for p in self if getattr(p.element, "is_" + family)])

    def dojo_compare(self, what="all", **kwargs):
        """Compare ecut convergence and Deltafactor, GBRV results"""
        import matplotlib.pyplot as plt
        show = kwargs.pop("show", True)
        what = list_strings(what)
        figs = []

        if all(p.dojo_report.has_trial("deltafactor") for p in self) and \
               any(k in what for k in ("all", "ecut")):

            fig_etotal, ax_list = plt.subplots(nrows=len(self), ncols=1, sharex=True, squeeze=True)
            #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(self), ncols=1, sharex=True, squeeze=True)
            figs.append(fig_etotal)

            for ax, pseudo in zip(ax_list, self):
                pseudo.dojo_report.plot_etotal_vs_ecut(ax=ax, show=False, label=pseudo.basename)
            if show: plt.show()

        if all(p.dojo_report.has_trial("deltafactor") for p in self) and \
               any(k in what for k in ("all", "df", "deltafactor")):

            fig_deltafactor, ax_grid = plt.subplots(nrows=5, ncols=len(self), sharex=True, sharey="row", squeeze=False)
            #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=5, ncols=len(self), sharex=True, sharey="row", squeeze=False))
            figs.append(fig_deltafactor)

            for ax_list, pseudo in zip(ax_grid.T, self):
                pseudo.dojo_report.plot_deltafactor_convergence(ax_list=ax_list, show=False)

            fig_deltafactor.suptitle(" vs ".join(p.basename for p in self))
            if show: plt.show()

        # Compare GBRV results
        if all(p.dojo_report.has_trial("gbrv_bcc") for p in self) and \
           any(k in what for k in ("all", "gbrv")):

            fig_gbrv, ax_grid = plt.subplots(nrows=2, ncols=len(self), sharex=True, sharey="row", squeeze=False)
            figs.append(fig_gbrv)
            #ax_list, fig, plt = get_axarray_fig_plt(ax_list, ncols=len(self), sharex=True, sharey="row", squeeze=False))

            for ax_list, pseudo in zip(ax_grid.T, self):
                pseudo.dojo_report.plot_gbrv_convergence(ax_list=ax_list, show=False)

            fig_gbrv.suptitle(" vs ".join(p.basename for p in self))
            if show: plt.show()

        return figs

    @classmethod
    @deprecated(replacement=from_dir)
    def from_directory(cls, path):
        pseudos = []
        for f in [os.path.join(path, fn) for fn in os.listdir(path)]:
            if os.path.isfile(f):
                try:
                    p = Pseudo.from_file(f)
                    if p:
                        pseudos.append(p)
                    else:
                        logger.info('Skipping file %s' % f)
                except:
                    logger.info('Skipping file %s' % f)
        if not pseudos:
            logger.warning('No pseudopotentials parsed from folder %s' % path)
            return None
        logger.info('Creating PseudoTable with %i pseudopotentials' % len(pseudos))
        return cls(pseudos)

try:
    from pandas import DataFrame
except ImportError:
    DataFrame = object


class DojoDataFrame(DataFrame):
    """Extends pandas DataFrame adding helper functions."""
    ALL_ACCURACIES = ("low", "normal", "high")

    ALL_TRIALS = (
        "ecut",
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",
        "phonon",
        #"phwoa"
    )

    _TRIALS2KEY = {
        "ecut": "ecut",
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "gbrv_bcc_a0_rel_err",
        "gbrv_fcc": "gbrv_fcc_a0_rel_err",
        "phonon": "all",
        #"phwoa": "all"
    }

    _TRIALS2YLABEL = {
        "ecut": "Ecut [Ha]",
        "deltafactor": "$\Delta$-factor [meV]",
        "gbrv_bcc": "BCC $\Delta a_0$ (%)",
        "gbrv_fcc": "FCC $\Delta a_0$ (%)",
        "phonon": "Phonons with ASR",
        #"phwoa": "Phonons without ASR"
    }

    ACC2PLTOPTS = dict(
        low=dict(color="red"),
        normal=dict(color="blue"),
        high=dict(color="green"),
    )

    for v in ACC2PLTOPTS.values():
        v.update(linewidth=2, linestyle='dashed', marker='o', markersize=8)

    def tabulate(self, columns=None, stream=sys.stdout):
        from tabulate import tabulate
        if columns is None:
            accuracies = self.ALL_ACCURACIES
            columns = [acc + "_dfact_meV" for acc in accuracies] 
            columns += [acc + "_ecut" for acc in accuracies] 
            columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies] 
            columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies] 

        #return self[columns].to_html()
        tablefmt = "grid"
        floatfmt=".2f"
        stream.write(tabulate(self[columns], headers="keys", tablefmt=tablefmt, floatfmt=floatfmt))

    def get_accuracy(self, accuracy):
        columns = [c for c in self if c.startswith(accuracy)]
        return self.__class__(data=self[columns])

    def get_trials(self, accuracies="all"):
        accuracies = self.ALL_ACCURACIES if accuracies == "all" else list_strings(accuracies)

        columns = [acc + "_dfact_meV" for acc in accuracies] 
        columns += [acc + "_ecut" for acc in accuracies] 
        columns += [acc + "_gbrv_fcc_a0_rel_err" for acc in accuracies] 
        columns += [acc + "_gbrv_bcc_a0_rel_err" for acc in accuracies] 
        return self.__class__(data=self[columns])

    def select_rows(self, rows):
        if not isinstance(rows, (list, tuple)): rows = [rows]
        
        data = []
        for index, entry in self.iterrows():
            element = _PTABLE[entry.Z]
            if element.row in rows:
                data.append(entry)

        return self.__class__(data=data)

    def select_family(self, family):
        data = []
        for index, entry in self.iterrows():
            element = _PTABLE[entry.Z]
            # e.g element.is_alkaline
            if getattr(element, "is_" + family):
                data.append(entry)
        return self.__class__(data=data)

    @add_fig_kwargs
    def plot_hist(self, what="dfact_meV", bins=400, **kwargs):
        import matplotlib.pyplot as plt
        fig, ax_list = plt.subplots(nrows=len(self.ALL_ACCURACIES), ncols=1, sharex=True, sharey=False, squeeze=True)
        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(self.ALL_ACCURACIES), ncols=1, sharex=True, sharey=False, squeeze=True)

        for acc, ax in zip(self.ALL_ACCURACIES, ax_list):
            col = acc + "_" + what
            #print(col)
            #self[col].hist(ax=ax, bins=bins, label=col)
            self[col].plot(ax=ax, kind="bar", label=col)

        return fig

    @add_fig_kwargs
    def plot_trials(self, trials="all", accuracies="all", **kwargs):
        import matplotlib.pyplot as plt
        trials = self.ALL_TRIALS if trials == "all" else list_strings(trials)
        accuracies = self.ALL_ACCURACIES if accuracies == "all" else list_strings(accuracies)

        fig, ax_list = plt.subplots(nrows=len(trials), ncols=1, sharex=True, sharey=False, squeeze=True)
        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(trials), ncols=1, sharex=True, sharey=False, squeeze=True)
                                                                                                                      
        # See also http://matplotlib.org/examples/pylab_examples/barchart_demo.html
        for i, (trial, ax) in enumerate(zip(trials, ax_list)):
            what = self._TRIALS2KEY[trial]
            ax.set_ylabel(self._TRIALS2YLABEL[trial])
            minval, maxval = np.inf, -np.inf
            for acc in accuracies:
                col = acc + "_" + what
                legend = i == 0 
                data = self[col]
                minval, maxval = min(minval, data.min()), max(maxval, data.max())
                data.plot(ax=ax, legend=legend, use_index=True, label=acc, **self.ACC2PLTOPTS[acc])
                #data.plot(ax=ax, kind="bar") 

                if i == 0:
                    ax.legend(loc='best', shadow=True, frameon=True) #fancybox=True)

            ax.set_xticks(range(len(data.index)))
            ax.set_xticklabels(data.index)
            #ax.set_xticklabels([root for root, ext in map(os.path.splitext, data.index)])

            # Set ylimits
            #stepsize = None
            #if "gbrv" in trial: 
            #    ax.hlines(0.0, 0, len(data.index))
            #    #start, end = -0.6, +0.6
            #    start, end = max(-0.6, minval), min(+0.6, maxval)
            #    if end - start < 0.05: end = start + 0.1
            #    ax.set_ylim(start, end)
            #    ax.yaxis.set_ticks(np.arange(start, end, 0.05))

            if trial == "deltafactor":
                #start, end = 0.0, 15
                start, end  = 0.0, min(15, maxval)
                ax.set_ylim(start, end)
                #ax.yaxis.set_ticks(np.arange(start, end, 0.1))

            #if stepsize is not None:
            #    start, end = ax.get_ylim()
            #    ax.yaxis.set_ticks(np.arange(start, end, stepsize))

            plt.setp(ax.xaxis.get_majorticklabels(), rotation=25)

        return fig


class DojoReportError(Exception):
    """Exception raised by DoJoReport."""


class DojoReport(dict):
    """Dict-like object with the dojo report."""

    _TRIALS2KEY = {
        "deltafactor": "dfact_meV",
        "gbrv_bcc": "a0_rel_err",
        "gbrv_fcc": "a0_rel_err",
        #"phwoa": "all",
        "phonon": "all"
    }

    ALL_ACCURACIES = ("low", "normal", "high")

    # List of dojo_trials
    # Remember to update the list if you add a new test to the DOJO_REPORT
    ALL_TRIALS = (
        "deltafactor",
        "gbrv_bcc",
        "gbrv_fcc",
        "phonon",
        #"phwoa"
    )

    # Tolerances on the deltafactor prime (in eV) used for the hints.
    #ATOLS = (1.0, 0.2, 0.04)
    ATOLS = (0.5, 0.1, 0.02)

    Error = DojoReportError

    @classmethod
    def from_file(cls, filepath):
        """Read the DojoReport from file."""
        with open(filepath, "rt") as fh:
            lines = fh.readlines()
            try:
                start = lines.index("<DOJO_REPORT>\n")
            except ValueError:
                return {}

            stop = lines.index("</DOJO_REPORT>\n")
            #print("start, stop" ,start, stop)
            #print("".join(lines[start+1:stop]))

            d = json.loads("".join(lines[start+1:stop]))
            return cls(**d)

    @classmethod
    def from_hints(cls, ppgen_ecut, symbol):
        """Initialize the DojoReport from an initial value of ecut in Hartree."""
        dense_right = np.arange(ppgen_ecut, ppgen_ecut + 6*2, step=2)
        dense_left = np.arange(max(ppgen_ecut-6, 2), ppgen_ecut, step=2)
        coarse_high = np.arange(ppgen_ecut + 15, ppgen_ecut + 35, step=5)

        ecut_list = list(dense_left) + list(dense_right) + list(coarse_high)
        return cls(ecut_list=ecut_list, symbol=symbol) #, **{k: {}: for k in self.ALL_TRIALS})

    def __init__(self, *args, **kwargs): 
        super(DojoReport, self).__init__(*args, **kwargs)

        try:
            for trial in self.ALL_TRIALS:
                # Convert ecut to float and build an OrderedDict (results are indexed by ecut in ascending order)
                try:
                    d = self[trial]
                except KeyError:
                    continue
                ecuts_keys = sorted([(float(k), k) for k in d], key=lambda t:t[0])
                ord = OrderedDict([(t[0], d[t[1]]) for t in ecuts_keys])
                self[trial] = ord

        except ValueError:
            raise self.Error('Error while initializing the dojo report')

    #def __str__(self):
    #    stream = six.moves.StringIO()
    #    pprint.pprint(self, stream=stream, indent=2, width=80)
    #    return stream.getvalue()

    @property
    def symbol(self):
        """Chemical symbol."""
        return self["symbol"]

    @property
    def element(self):
        """Element object."""
        return Element(self.symbol)

    @property
    def has_hints(self):
        """True if hints on cutoff energy are present."""
        return "hints" in self

    @property
    def ecuts(self):
        """Numpy array with the list of ecuts that should be present in the dojo_trial sub-dicts"""
        return self["ecuts"]

    @property
    def trials(self):
        """Set of strings with the trials present in the report."""
        return set(list(self.keys())).intersection(self.ALL_TRIALS)

    def has_trial(self, dojo_trial, ecut=None):
        """
        True if the dojo_report contains dojo_trial with the given ecut.
        If ecut is not, we test if dojo_trial is present.
        """
        if dojo_trial not in self.ALL_TRIALS:
            raise self.Error("dojo_trial `%s` is not a registered DOJO TRIAL" % dojo_trial)

        if ecut is None:
            return dojo_trial in self
        else:
            #key = self._ecut2key(ecut)
            key = ecut
            try:
                self[dojo_trial][key]
                return True
            except KeyError:
                return False

    def add_ecuts(self, new_ecuts):
        """Add a list of new ecut values."""
        # Be careful with the format here! it should be %.1f
        # Select the list of ecuts reported in the DOJO section.
        prev_ecuts = self["ecuts"]

        for i in range(len(prev_ecuts)-1):
            if prev_ecuts[i] >= prev_ecuts[i+1]:
                raise self.Error("Ecut list is not ordered:\n %s" % prev_ecuts)

        from monty.bisect import find_le
        for e in new_ecuts:
            # Find rightmost value less than or equal to x.
            if e < prev_ecuts[0]:
                i = 0
            elif e > prev_ecuts[-1]:
                i = len(prev_ecuts)
            else:
                i = find_le(prev_ecuts, e)
                assert prev_ecuts[i] != e
                i += 1

            prev_ecuts.insert(i, e)

    def add_hints(self, hints):
        hints_dict = {
           "low": {'ecut': hints[0]},
           "normal" : {'ecut': hints[1]},
           "high" : {'ecut': hints[2]}
                     }
        self["hints"] = hints_dict

    #def validate(self, hints):
    #    Add md5 hash value
    #    self["validated"] = True

    @staticmethod
    def _ecut2key(ecut):
        """Convert ecut to a valid key. ecut can be either a string or a float."""
        if is_string(ecut):
            # Validate string
            i = ecut.index(".")
            if len(ecut[i+1:]) != 1:
                raise ValueError("string %s must have one digit")
            return ecut

        else:
            # Assume float
            return "%.1f" % ecut

    def add_entry(self, dojo_trial, ecut, d, overwrite=False):
        """
        Add an entry computed with the given ecut to the sub-dictionary associated to dojo_trial.

        Args:
            dojo_trial:
            ecut:
            d:
            overwrite:
        """
        if dojo_trial not in self.ALL_TRIALS:
            raise ValueError("%s is not a registered trial")
        section = self.get(dojo_trial, {})

        key = self._ecut2key(ecut)
        if key in section and not overwrite:
            raise self.Error("Cannot overwrite key %s in dojo_trial %s" % (key, dojo_trial))

        section[key] = d

    def find_missing_entries(self):
        """
        check the DojoReport. 
        This function tests if each trial contains an ecut entry. 
        Return a dictionary {trial_name: [list_of_missing_ecuts]}
        mapping the name of the Dojo trials to the list of ecut values that are missing 
        """
        d = {}

        for trial in self.ALL_TRIALS:
            data = self.get(trial, None)
            if data is None:
                # Gbrv results do not contain noble gases so ignore the error
                if "gbrv" in trial and self.element.is_noble_gas: 
                    assert data is None
                    continue
                d[trial] = self.ecuts

            else:
                computed_ecuts = self[trial].keys()
                for e in self.ecuts:
                    if e not in computed_ecuts:
                        if trial not in d: d[trial] = []
                        d[trial].append(e)

        if not d:
            assert len(computed_ecuts) == len(self.ecuts)

        return d

    def print_table(self, stream=sys.stdout):
        from monty.pprint import pprint_table
        pprint_table(self.get_dataframe(), out=stream)

    @add_fig_kwargs
    def plot_etotal_vs_ecut(self, ax=None, inv_ecut=False, **kwargs):
        """
        plot the convergence of the total energy as function of the energy cutoff ecut

        Args:
            ax: matplotlib Axes, if ax is None a new figure is created.

        Returns:
            `matplotlib` figure.
        """
        # Extract the total energy of the AE relaxed structure (4).
        d = OrderedDict([(ecut, data["etotals"][4]) for ecut, data in self["deltafactor"].items()])

        # Ecut mesh in Ha
        ecuts = np.array(list(d.keys()))
        ecut_min, ecut_max = np.min(ecuts), np.max(ecuts)

        # Energies per atom in meV and difference wrt 'converged' value
        num_sites = [v["num_sites"] for v in self["deltafactor"].values()][0]
        etotals_mev = np.array([d[e] for e in ecuts]) * 1000  / num_sites
        ediffs = etotals_mev - etotals_mev[-1]

        ax, fig, plt = get_ax_fig_plt(ax)
        #ax.yaxis.set_view_interval(-5, 5)

        lines, legends = [], []

        xs = 1/ecuts if inv_ecut else ecuts
        ys = etotals_mev if inv_ecut else ediffs

        line, = ax.plot(xs, ys, "-o", color="blue") #, linewidth=3.0, markersize=15)
        lines.append(line)

        label = kwargs.pop("label", None)
        if label is not None: ax.legend(lines, [label], loc='best', shadow=True)

        high_hint = self["ppgen_hints"]["high"]["ecut"]
        #ax.vlines(high_hint, min(ediffs), max(ediffs))
        #ax.vlines(high_hint, 0.5, 1.5)
        #ax.scatter([high_hint], [1.0], s=20) #, c='b', marker='o', cmap=None, norm=None)
        #ax.arrow(high_hint, 1, 0, 0.2, head_width=0.05, head_length=0.1, fc='k', ec='k',head_starts_at_zero=False)

        #ax.hlines(5, ecut_min, ecut_max, label="5.0")
        #ax.hlines(1, ecut_min, ecut_max, label="1.0")
        #ax.hlines(0.5, ecut_min, ecut_max, label="0.2")

        # Set xticks and labels.
        ax.grid(True)
        ax.set_xlabel("Ecut [Ha]")
        ax.set_xticks(xs)
        ax.set_ylabel("Delta Etotal/natom [meV]")
        #ax.set_xlim(0, max(xs))

        # Use logscale if possible.
        if all(ediffs[:-1] > 0): 
            ax.set_yscale("log")
            ax.set_xlim(xs[0]-1, xs[-2]+1)

        return fig

    @add_fig_kwargs
    def plot_deltafactor_eos(self, ax=None, **kwargs):
        """
        plot the EOS computed with the deltafactor setup.

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.
            
        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure.
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        trial = "deltafactor"
        ecuts = self[trial].keys()
        num_ecuts = len(ecuts)

        cmap = kwargs.pop("cmap", None)
        if cmap is None: cmap = plt.get_cmap("jet")

        for i, ecut in enumerate(ecuts):
            d = self[trial][ecut]
            num_sites, volumes, etotals = d["num_sites"], np.array(d["volumes"]), np.array(d["etotals"])

            # Use same fit as the one employed for the deltafactor.
            eos_fit = EOS.DeltaFactor().fit(volumes/num_sites, etotals/num_sites)

            label = "ecut %.1f" % ecut if i % 2 == 0 else ""
            label = "ecut %.1f" % ecut 
            eos_fit.plot(ax=ax, text=False, label=label, color=cmap(i/num_ecuts, alpha=1), show=False)

        return fig

    def get_ecut_dfactprime(self):
        data = self["deltafactor"]
        ecuts, values= data.keys(), []
        values = np.array([data[e]["dfactprime_meV"] for e in ecuts])
        return np.array(ecuts), values

    def compute_hints(self):
        ecuts, dfacts = self.get_ecut_dfactprime()
        abs_diffs = np.abs((dfacts - dfacts[-1]))
        #print(list(zip(ecuts, dfacts)))
        #print(abs_diffs)

        hints = 3 * [None]
        for ecut, adiff in zip(ecuts, abs_diffs):
            for i in range(3):
                if adiff <= self.ATOLS[i] and hints[i] is None:
                    hints[i] = ecut
                if adiff > self.ATOLS[i]: 
                    hints[i] = None
        return hints

    def check(self):
        """
        Check the dojo report for inconsistencies.
        Return a string with the errors found in the DOJO_REPORT.
        """
        errors = []
        app = errors.append

        if "version" not in self:
            app("version is missing")

        if "ppgen_hints" not in self:
            app("version is missing")

        if "md5" not in self:
            app("md5 checksum is missing!")

        # Check if we have computed each trial for the full set of ecuts in global_ecuts
        global_ecuts = self.ecuts

        missing = defaultdict(list)
        for trial in self.ALL_TRIALS:
            for ecut in global_ecuts:
                if not self.has_trial(trial, ecut=ecut):
                    missing[trial].append(ecut)

        if missing:
            app("The following list of ecut energies is missing:")
            for trial, ecuts in missing.items():
                app("%s: %s" % (trial, ecuts))
            
        return "\n".join(errors)

    @add_fig_kwargs
    def plot_deltafactor_convergence(self, code="WIEN2k", what=None, ax_list=None, **kwargs):
        """
        plot the convergence of the deltafactor parameters wrt ecut.

        Args:
            code: Reference code
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        all = ["dfact_meV", "dfactprime_meV", "v0", "b0_GPa", "b1"]
        if what is None:
            keys = all
        else:
            what = list_strings(what)
            if what[0].startswith("-"):
                # Exclude keys
                #print([type(w) for w in what])
                what = [w[1:] for w in what]
                keys = [k for k in all if k not in what]
            else:
                keys = what
            
        # get reference entry
        from pseudo_dojo.refdata.deltafactor import df_database
        reference = df_database().get_entry(symbol=self.symbol, code=code)

        d = self["deltafactor"]
        ecuts = list(d.keys())

        import matplotlib.pyplot as plt
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(keys), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(keys), ncols=1, sharex=True, squeeze=False)

        if len(keys) != len(ax_list): 
            raise ValueError("len(keys)=%s != len(ax_list)=%s" %  (len(keys), len(ax_list)))

        for i, (ax, key) in enumerate(zip(ax_list, keys)):
            values = np.array([float(d[ecut][key]) for ecut in ecuts])
            #try:
            refval = getattr(reference, key)
            #except AttributeError:
            #    refval = 0.0

            # Plot difference pseudo - ref.
            ax.plot(ecuts, values - refval, "o-")

            ax.grid(True)
            ax.set_ylabel("$\Delta$" + key)
            if i == len(keys) - 1: ax.set_xlabel("Ecut [Ha]")

            if key == "dfactprime_meV":
                # Add horizontal lines (used to find hints for ecut).
                last = values[-1]
                xmin, xmax = min(ecuts), max(ecuts)
                for pad, color in zip(self.ATOLS, ("blue", "red", "violet")):
                    ax.hlines(y=last + pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1, linestyles='dashed')
                    ax.hlines(y=last - pad, xmin=xmin, xmax=xmax, colors=color, linewidth=1, linestyles='dashed')
              
                # Set proper limits so that we focus on the relevant region.
                ax.set_ylim(last - 1.1*self.ATOLS[0], last + 1.1*self.ATOLS[0])

        return fig

    @add_fig_kwargs
    def plot_gbrv_eos(self, struct_type, ax=None, **kwargs):
        """
        Uses Matplotlib to plot the EOS computed with the GBRV setup

        Args:
            ax: matplotlib :class:`Axes` or None if a new figure should be created.

        ================  ==============================================================
        kwargs            Meaning
        ================  ==============================================================
        cmap              Color map. default `jet`
        ================  ==============================================================

        Returns:
            `matplotlib` figure or None if the GBRV test is not present
        """
        ax, fig, plt = get_ax_fig_plt(ax)

        trial = "gbrv_" + struct_type
        # Handle missing entries: noble gases, Hg ...
        if trial not in self: return None
        ecuts = self[trial].keys()
        num_ecuts = len(ecuts)

        cmap = kwargs.pop("cmap", None)
        if cmap is None: cmap = plt.get_cmap("jet")

        for i, ecut in enumerate(ecuts):
            d = self[trial][ecut]
            volumes, etotals = np.array(d["volumes"]), np.array(d["etotals"])

            eos_fit = EOS.Quadratic().fit(volumes, etotals)
            label = "ecut %.1f" % ecut if i % 2 == 0 else ""
            label = "ecut %.1f" % ecut 
            eos_fit.plot(ax=ax, text=False, label=label, color=cmap(i/num_ecuts, alpha=1), show=False)

        return fig

    @add_fig_kwargs
    def plot_gbrv_convergence(self, ax_list=None, **kwargs):
        """
        Uses Matplotlib to plot the convergence of the GBRV parameters wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        import matplotlib.pyplot as plt
        stypes = ("fcc", "bcc")
        if ax_list is None:
            fig, ax_list = plt.subplots(nrows=len(stypes), ncols=1, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            fig = plt.gcf()

        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(stypes), ncols=1, sharex=True, squeeze=False)

        if len(stypes) != len(ax_list): 
            raise ValueError("len(stypes)=%s != len(ax_list)=%s" %  (len(stypes), len(ax_list)))

        for i, (ax, stype) in enumerate(zip(ax_list, stypes)):
            trial = "gbrv_" + stype
            d = self[trial]
            ecuts = list(d.keys())
            values = np.array([float(d[ecut]["a0_rel_err"]) for ecut in ecuts])

            ax.grid(True)
            ax.set_ylabel("$\Delta$" + trial + "a0_rel_err")

            # Plot difference pseudo - ref.
            ax.plot(ecuts, values, "bo-")
            #ax.hlines(y=0.0, xmin=min(ecuts), xmax=max(ecuts), color="red")
            if i == len(ax_list) - 1: ax.set_xlabel("Ecut [Ha]")

        return fig

    @add_fig_kwargs
    def plot_phonon_convergence(self, ax_list=None, **kwargs):
        """
        Plot the convergence of the phonon modes wrt ecut.

        Args:
            ax_list: List of matplotlib Axes, if ax_list is None a new figure is created

        Returns:
            `matplotlib` figure.
        """
        d = self["phonon"]
        ecuts = list(d.keys())

        l = [(ecut, float(ecut)) for ecut in ecuts]
        s = sorted(l, key=lambda t: t[1])
        max_ecut = s[-1][0]
        s_ecuts = [ecut[0] for ecut in s]

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(nrows=2, sharex=True)
        #ax_list, fig, plt = get_axarray_fig_plt(ax_list, nrows=len(keys), ncols=1, sharex=True, squeeze=False)

        fmin, fmax = np.inf, -np.inf
        for i, v in enumerate(d[ecuts[0]]):
            values1 = np.array([float(d[ecut][i]) for ecut in s_ecuts])
            fmin = min(fmin, values1.min())
            fmax = max(fmax, values1.max())

            ax[0].plot(s_ecuts, values1, "o-")
            ax[0].grid(True)
            ax[0].set_ylabel("phonon modes [meV] (asr==2)")
            ax[0].set_xlabel("Ecut [Ha]")

            values2 = np.array([float(d[ecut][i]) - float(d[max_ecut][i]) for ecut in s_ecuts])

            ax[1].plot(s_ecuts, values2, "o-")
            ax[1].grid(True)
            ax[1].set_ylabel("w - w(ecut_max) [meV]")
            ax[1].set_xlabel("Ecut [Ha]")

        # Adjust limits.
        fmin -= 10 
        fmax += 10 
        ax[0].set_ylim(fmin, fmax)

        return fig
