# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides objects describing the basic parameters of the
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..
"""
from __future__ import unicode_literals, division, print_function

import abc
import collections
import json
import logging
import os
import sys
import numpy as np
import six

from collections import OrderedDict, defaultdict, namedtuple
from monty.collections import AttrDict, Namespace
from tabulate import tabulate
#from monty.dev import deprecated
from monty.functools import lazy_property
from monty.itertools import iterator_from_slice
from monty.json import MSONable, MontyDecoder
from monty.os.path import find_exts
from monty.string import list_strings, is_string
from pymatgen.core.periodic_table import Element
from pymatgen.core.xcfunc import XcFunc
from pymatgen.util.serialization import pmg_serialize
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig_plt

logger = logging.getLogger(__name__)


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
        for lineno, line in enumerate(fh):
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


class Pseudo(six.with_metaclass(abc.ABCMeta, MSONable, object)):
    """
    Abstract base class defining the methods that must be
    implemented by the concrete pseudopotential sub-classes.
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
        Build an instance of a concrete Pseudo subclass from filename.
        Note: the parser knows the concrete class that should be instantiated
        Client code should rely on the abstract interface provided by Pseudo.
        """
        return PseudoParser().parse(filename)

    def __eq__(self, other):
        if other is None: return False
        return (self.md5 == other.md5 and
                self.__class__ == other.__class__ and
                self.Z == other.Z and
                self.Z_val == other.Z_val and
                self.l_max == other.l_max )

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        try:
            return "<%s at %s>" % (self.__class__.__name__, os.path.relpath(self.filepath))
        except:
            # relpath can fail if the code is executed in demon mode.
            return "<%s at %s>" % (self.__class__.__name__, self.filepath)

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        app("<%s: %s>" % (self.__class__.__name__, self.basename))
        app("  summary: " + self.summary.strip())
        app("  number of valence electrons: %s" % self.Z_val)
        app("  maximum angular momentum: %s" % l2str(self.l_max))
        app("  angular momentum for local part: %s" % l2str(self.l_local))
        app("  XC correlation: %s" % self.xc)
        app("  supports spin-orbit: %s" % self.supports_soc)

        if self.isnc:
            app("  radius for non-linear core correction: %s" % self.nlcc_radius)

        if self.has_hints:
            for accuracy in ("low", "normal", "high"):
                hint = self.hint_for_accuracy(accuracy=accuracy)
                app("  hint for %s accuracy: %s" % (accuracy, str(hint)))

        return "\n".join(lines)

    @property
    @abc.abstractmethod
    def summary(self):
        """String summarizing the most important properties."""

    @property
    def filepath(self):
        return os.path.abspath(self.path)

    @property
    def basename(self):
        """File basename."""
        return os.path.basename(self.filepath)

    @property
    @abc.abstractmethod
    def Z(self):
        """The atomic number of the atom."""

    @property
    @abc.abstractmethod
    def Z_val(self):
        """Valence charge."""

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def element(self):
        """Pymatgen :class:`Element`."""
        try:
            return Element.from_Z(self.Z)
        except (KeyError, IndexError):
            return Element.from_Z(int(self.Z))

    @property
    def symbol(self):
        """Element symbol."""
        return self.element.symbol

    @property
    @abc.abstractmethod
    def l_max(self):
        """Maximum angular momentum."""

    @property
    @abc.abstractmethod
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
        #if self.has_dojo_report and "md5" in self.dojo_report: return self.dojo_report["md5"]
        return self.compute_md5()

    def compute_md5(self):
        """Compute and erturn MD5 hash value."""
        import hashlib
        with open(self.path, "rt") as fh:
            text = fh.read()
            m = hashlib.md5(text.encode("utf-8"))
            return m.hexdigest()

    @property
    @abc.abstractmethod
    def supports_soc(self):
        """
        True if the pseudo can be used in a calculation with spin-orbit coupling.
        Base classes should provide a concrete implementation that computes this value.
        """

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
            filepath=self.filepath,
            #xc=self.xc.as_dict(),
        )

    @classmethod
    def from_dict(cls, d):
        new = cls.from_file(d['filepath'])

        # Consistency test based on md5
        if "md5" in d and d["md5"] != new.md5:
            raise ValueError("The md5 found in file does not agree with the one in dict\n"
            "Received %s\nComputed %s" % (d["md5"], new.md5))

        return new

    def as_tmpfile(self, tmpdir=None):
        """
        Copy the pseudopotential to a temporary a file and returns a new pseudopotential object.
        Useful for unit tests in which we have to change the content of the file.

        Args:
            tmpdir: If None, a new temporary directory is created and files are copied here
                else tmpdir is used.
        """
        import tempfile, shutil
        tmpdir = tempfile.mkdtemp() if tmpdir is None else tmpdir
        new_path = os.path.join(tmpdir, self.basename)
        shutil.copy(self.filepath, new_path)

        # Copy dojoreport file if present.
        root, ext = os.path.splitext(self.filepath)
        djrepo = root + ".djrepo"
        if os.path.exists(djrepo):
            shutil.copy(djrepo, os.path.join(tmpdir, os.path.basename(djrepo)))

        # Build new object and copy dojo_report if present.
        new = self.__class__.from_file(new_path)
        if self.has_dojo_report: new.dojo_report = self.dojo_report.deepcopy()

        return new

    @property
    def has_dojo_report(self):
        """True if the pseudo has an associated `DOJO_REPORT` section."""
        return hasattr(self, "dojo_report") and bool(self.dojo_report)

    @property
    def djrepo_path(self):
        """The path of the djrepo file. None if file does not exist."""
        root, ext = os.path.splitext(self.filepath)
        path = root + ".djrepo"
        return path
        #if os.path.exists(path): return path
        #return None

    def hint_for_accuracy(self, accuracy="normal"):
        """
        Returns a :class:`Hint` object with the suggested value of ecut [Ha] and
        pawecutdg [Ha] for the given accuracy.
        ecut and pawecutdg are set to zero if no hint is available.

        Args:
            accuracy: ["low", "normal", "high"]
        """
        if not self.has_dojo_report:
            return Hint(ecut=0., pawecutdg=0.)

        # Get hints from dojoreport. Try first in hints then in ppgen_hints.
        if "hints" in self.dojo_report:
            return Hint.from_dict(self.dojo_report["hints"][accuracy])
        elif "ppgen_hints" in self.dojo_report:
            return Hint.from_dict(self.dojo_report["ppgen_hints"][accuracy])
        return Hint(ecut=0., pawecutdg=0.)

    @property
    def has_hints(self):
        """
        True if self provides hints on the cutoff energy.
        """
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

        if self.ispaw and pawecutdg is None: pawecutdg = ecut * 4
        inp = gs_input(structure, pseudos=[self], ecut=ecut, pawecutdg=pawecutdg,
                       spin_mode="unpolarized", kppa=1)
        # Add prtpsps = -1 to make Abinit print the PSPS.nc file and stop.
        inp["prtpsps"] = -1

        # Build temporary task and run it (ignore retcode because we don't exit cleanly)
        task = AbinitTask.temp_shell_task(inp)
        task.start_and_wait()

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

    @property
    @abc.abstractmethod
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

    @property
    @abc.abstractmethod
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
        self.header = header
        self._summary = header.summary

        # Build xc from header.
        self.xc = XcFunc.from_abinit_ixc(header["pspxc"])

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

    @property
    def supports_soc(self):
        # Treate ONCVPSP pseudos
        if self._pspcod == 8:
            switch = self.header["extension_switch"]
            if switch in (0, 1): return False
            if switch in (2, 3): return True
            raise ValueError("Don't know how to handle extension_switch: %s" % switch)

        # TODO Treat HGH HGHK pseudos

        # As far as I know, other Abinit pseudos do not support SOC.
        return False


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

    @property
    def supports_soc(self):
        return True


class Hint(object):
    """
    Suggested value for the cutoff energy [Hartree units]
    and the cutoff energy for the dense grid (only for PAW pseudos).
    """
    def __init__(self, ecut, pawecutdg=None):
        self.ecut = ecut
        self.pawecutdg = ecut if pawecutdg is None else pawecutdg

    def __str__(self):
        if self.pawecutdg is not None:
            return "ecut: %s, pawecutdg: %s" % (self.ecut, self.pawecutdg)
        else:
            return "ecut: %s" % (self.ecut)

    @pmg_serialize
    def as_dict(self):
        return dict(ecut=self.ecut, pawecutdg=self.pawecutdg)

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k, v in d.items() if not k.startswith("@")})


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
        # Sanitize keys: In some case we might get strings in the form: foo[,bar]
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
        logger.warning("Converting float %s to int %s" % (float_num, int_num))
        return int_num


class NcAbinitHeader(AbinitHeader):
    """The abinit header found in the NC pseudopotential files."""
    _attr_desc = namedtuple("att", "default astype")

    _VARS = {
        # Mandatory
        "zatom": _attr_desc(None, _int_from_str),
        "zion": _attr_desc(None, float),
        "pspdat": _attr_desc(None, float),
        "pspcod": _attr_desc(None, int),
        "pspxc": _attr_desc(None, int),
        "lmax": _attr_desc(None, int),
        "lloc": _attr_desc(None, int),
        "r2well": _attr_desc(None, float),
        "mmax": _attr_desc(None, float),
        # Optional variables for non linear-core correction. HGH does not have it.
        "rchrg": _attr_desc(0.0,  float),  # radius at which the core charge vanish (i.e. cut-off in a.u.)
        "fchrg": _attr_desc(0.0,  float),
        "qchrg": _attr_desc(0.0,  float),
    }
    del _attr_desc

    def __init__(self, summary, **kwargs):
        super(NcAbinitHeader, self).__init__()

        # pseudos generated by APE use llocal instead of lloc.
        if "llocal" in kwargs:
            kwargs["lloc"] = kwargs.pop("llocal")

        self.summary = summary.strip()

        for key, desc in NcAbinitHeader._VARS.items():
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

        # Add remaining arguments, e.g. extension_switch
        if kwargs:
            self.update(kwargs)

    @staticmethod
    def fhi_header(filename, ppdesc):
        """
        Parse the FHI abinit header. Example:

        Troullier-Martins psp for element  Sc        Thu Oct 27 17:33:22 EDT 1994
         21.00000   3.00000    940714                zatom, zion, pspdat
           1    1    2    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
        1.80626423934776     .22824404341771    1.17378968127746   rchrg,fchrg,qchrg
        """
        lines = _read_nlines(filename, 4)

        try:
            header = _dict_from_lines(lines[:4], [0, 3, 6, 3])
        except ValueError:
            # The last record with rchrg ... seems to be optional.
            header = _dict_from_lines(lines[:3], [0, 3, 6])

        summary = lines[0]

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def hgh_header(filename, ppdesc):
        """
        Parse the HGH abinit header. Example:

        Hartwigsen-Goedecker-Hutter psp for Ne,  from PRB58, 3641 (1998)
           10   8  010605 zatom,zion,pspdat
         3 1   1 0 2001 0  pspcod,pspxc,lmax,lloc,mmax,r2well
        """
        lines = _read_nlines(filename, 3)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def gth_header(filename, ppdesc):
        """
        Parse the GTH abinit header. Example:

        Goedecker-Teter-Hutter  Wed May  8 14:27:44 EDT 1996
        1   1   960508                     zatom,zion,pspdat
        2   1   0    0    2001    0.       pspcod,pspxc,lmax,lloc,mmax,r2well
        0.2000000 -4.0663326  0.6778322 0 0     rloc, c1, c2, c3, c4
        0 0 0                              rs, h1s, h2s
        0 0                                rp, h1p
          1.36 .2   0.6                    rcutoff, rloc
        """
        lines = _read_nlines(filename, 7)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def oncvpsp_header(filename, ppdesc):
        """
        Parse the ONCVPSP abinit header. Example:

        Li    ONCVPSP  r_core=  2.01  3.02
              3.0000      3.0000      140504    zatom,zion,pspd
             8     2     1     4   600     0    pspcod,pspxc,lmax,lloc,mmax,r2well
          5.99000000  0.00000000  0.00000000    rchrg fchrg qchrg
             2     2     0     0     0    nproj
             0                 extension_switch
           0                        -2.5000025868368D+00 -1.2006906995331D+00
             1  0.0000000000000D+00  0.0000000000000D+00  0.0000000000000D+00
             2  1.0000000000000D-02  4.4140499497377D-02  1.9909081701712D-02
        """
        lines = _read_nlines(filename, 6)

        header = _dict_from_lines(lines[:3], [0, 3, 6])
        summary = lines[0]

        # Replace pspd with pspdata
        header.update({'pspdat': header['pspd']})
        header.pop('pspd')

        # Read extension switch
        header["extension_switch"] = int(lines[5].split()[0])

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def tm_header(filename, ppdesc):
        """
        Parse the TM abinit header. Example:

        Troullier-Martins psp for element Fm         Thu Oct 27 17:28:39 EDT 1994
        100.00000  14.00000    940714                zatom, zion, pspdat
           1    1    3    0      2001    .00000      pspcod,pspxc,lmax,lloc,mmax,r2well
           0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
           .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
           1   3.116   4.632    1   3.4291849        l,e99.0,e99.9,nproj,rcpsp
           .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
           2   4.557   6.308    1   2.1865358        l,e99.0,e99.9,nproj,rcpsp
           .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
           3  23.251  29.387    1   2.4776730        l,e99.0,e99.9,nproj,rcpsp
           .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
           3.62474762267880     .07409391739104    3.07937699839200   rchrg,fchrg,qchrg
        """
        lines = _read_nlines(filename, -1)
        header = []

        for lineno, line in enumerate(lines):
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

        return NcAbinitHeader(summary, **header)


class PawAbinitHeader(AbinitHeader):
    """The abinit header found in the PAW pseudopotential files."""
    _attr_desc = namedtuple("att", "default astype")

    _VARS = {
        "zatom": _attr_desc(None, _int_from_str),
        "zion": _attr_desc(None, float),
        "pspdat": _attr_desc(None, float),
        "pspcod": _attr_desc(None, int),
        "pspxc": _attr_desc(None, int),
        "lmax": _attr_desc(None, int),
        "lloc": _attr_desc(None, int),
        "mmax": _attr_desc(None, int),
        "r2well": _attr_desc(None, float),
        "pspfmt": _attr_desc(None, str),
        "creatorID": _attr_desc(None, int),
        "basis_size": _attr_desc(None, int),
        "lmn_size": _attr_desc(None, int),
        "orbitals": _attr_desc(None, list),
        "number_of_meshes": _attr_desc(None, int),
        "r_cut": _attr_desc(None, float), # r_cut(PAW) in the header
        "shape_type": _attr_desc(None, int),
        "rshape": _attr_desc(None, float),
    }
    del _attr_desc

    def __init__(self, summary, **kwargs):
        super(PawAbinitHeader, self).__init__()

        self.summary = summary.strip()

        for key, desc in self._VARS.items():
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
        """
        Parse the PAW abinit header. Examples:

        Paw atomic data for element Ni - Generated by AtomPAW (N. Holzwarth) + AtomPAW2Abinit v3.0.5
          28.000  18.000 20061204               : zatom,zion,pspdat
          7  7  2 0   350 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well
         paw3 1305                              : pspfmt,creatorID
          5 13                                  : basis_size,lmn_size
         0 0 1 1 2                              : orbitals
         3                                      : number_of_meshes
         1 3  350 1.1803778368E-05 3.5000000000E-02 : mesh 1, type,size,rad_step[,log_step]
         2 1  921 2.500000000000E-03                : mesh 2, type,size,rad_step[,log_step]
         3 3  391 1.1803778368E-05 3.5000000000E-02 : mesh 3, type,size,rad_step[,log_step]
          2.3000000000                          : r_cut(SPH)
         2 0.

        Another format:

        C  (US d-loc) - PAW data extracted from US-psp (D.Vanderbilt) - generated by USpp2Abinit v2.3.0
           6.000   4.000 20090106               : zatom,zion,pspdat
          7 11  1 0   560 0.                    : pspcod,pspxc,lmax,lloc,mmax,r2well
         paw4 2230                              : pspfmt,creatorID
          4  8                                  : basis_size,lmn_size
         0 0 1 1                                : orbitals
         5                                      : number_of_meshes
         1 2  560 1.5198032759E-04 1.6666666667E-02 : mesh 1, type,size,rad_step[,log_step]
         2 2  556 1.5198032759E-04 1.6666666667E-02 : mesh 2, type,size,rad_step[,log_step]
         3 2  576 1.5198032759E-04 1.6666666667E-02 : mesh 3, type,size,rad_step[,log_step]
         4 2  666 1.5198032759E-04 1.6666666667E-02 : mesh 4, type,size,rad_step[,log_step]
         5 2  673 1.5198032759E-04 1.6666666667E-02 : mesh 5, type,size,rad_step[,log_step]
          1.5550009124                          : r_cut(PAW)
         3 0.                                   : shape_type,rshape

        Yet nnother one:

        Paw atomic data for element Si - Generated by atompaw v3.0.1.3 & AtomPAW2Abinit v3.3.1
          14.000   4.000 20120814               : zatom,zion,pspdat
          7      11  1 0   663 0.               : pspcod,pspxc,lmax,lloc,mmax,r2well
         paw5 1331                              : pspfmt,creatorID
          4  8                                  : basis_size,lmn_size
         0 0 1 1                                : orbitals
         5                                      : number_of_meshes
         1 2  663 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 1, type,size,rad_step[,log_step]
         2 2  658 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 2, type,size,rad_step[,log_step]
         3 2  740 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 3, type,size,rad_step[,log_step]
         4 2  819 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 4, type,size,rad_step[,log_step]
         5 2  870 8.2129718540404674E-04 1.1498160595656655E-02 : mesh 5, type,size,rad_step[,log_step]
          1.5669671236                          : r_cut(PAW)
         2 0.                                   : shape_type,rshape
        """
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
        1: ppdesc(1, "TM", "NC", None),
        2: ppdesc(2, "GTH", "NC", None),
        3: ppdesc(3, "HGH", "NC", None),
        4: ppdesc(4, "Teter", "NC", None),
        #5: ppdesc(5, "NC",     , None),
        6: ppdesc(6, "FHI", "NC", None),
        7: ppdesc(6, "PAW_abinit_text", "PAW", None),
        8: ppdesc(8, "ONCVPSP", "NC", None),
       10: ppdesc(10, "HGHK", "NC", None),
    })
    del ppdesc
    # renumber functionals from oncvpsp todo confrim that 3 is 2
    #_FUNCTIONALS = {1: {'n': 4, 'name': 'Wigner'},
    #                2: {'n': 5, 'name': 'HL'},
    #                3: {'n': 2, 'name': 'PWCA'},
    #                4: {'n': 11, 'name': 'PBE'}}

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
        for i, ext in enumerate(exclude_exts):
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

            for lineno, line in enumerate(lines):

                if lineno == 2:
                    try:
                        tokens = line.split()
                        pspcod, pspxc = map(int, tokens[:2])
                    except:
                        msg = "%s: Cannot parse pspcod, pspxc in line\n %s" % (filename, line)
                        logger.critical(msg)
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
            logger.critical("Cannot find ppdesc in %s" % path)
            return None

        psp_type = ppdesc.psp_type

        parsers = {
            "FHI": NcAbinitHeader.fhi_header,
            "GTH": NcAbinitHeader.gth_header,
            "TM": NcAbinitHeader.tm_header,
            "Teter": NcAbinitHeader.tm_header,
            "HGH": NcAbinitHeader.hgh_header,
            "HGHK": NcAbinitHeader.hgh_header,
            "ONCVPSP": NcAbinitHeader.oncvpsp_header,
            "PAW_abinit_text": PawAbinitHeader.paw_header,
        }

        try:
            header = parsers[ppdesc.name](path, ppdesc)
        except Exception:
            raise self.Error(path + ":\n" + straceback())

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

        # Build xc from header.
        xc_info = root.find("xc_functional").attrib
        self.xc = XcFunc.from_type_name(xc_info["type"], xc_info["name"])

        # Old XML files do not define this field!
        # In this case we set the PAW radius to None.
        #self._paw_radius = float(root.find("PAW_radius").attrib["rpaw"])

        #self.ae_energy = {k: float(v) for k,v in root.find("ae_energy").attrib.items()}
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

            self.rad_grids[gid] = self._eval_grid(grid_params)

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

    @property
    def supports_soc(self):
        """
        Here I assume that the ab-initio code can treat the SOC within the on-site approximation
        """
        return True

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
            for mesh, values, attrib in self._parse_all_radfuncs("ae_partial_wave"):
                state = attrib["state"]
                #val_state = self.valence_states[state]
                self._ae_partial_waves[state] = RadialFunction(mesh, values)

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
                #val_state = self.valence_states[state]
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
                #val_state = self.valence_states[state]
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
        ax.set_ylabel("$r\\phi,\\, r\\tilde\\phi\\, [Bohr]^{-\\frac{1}{2}}$")

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
        ax.set_ylabel("$r\\tilde p\\, [Bohr]^{-\\frac{1}{2}}$")

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


class PseudoTable(six.with_metaclass(abc.ABCMeta, collections.Sequence, MSONable, object)):
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
            if not isinstance(pseudo, Pseudo):
                pseudo = Pseudo.from_file(pseudo)
            if pseudo is not None:
                self._pseudos_with_z[pseudo.Z].append(pseudo)

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
        return self.to_table()

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

    #def max_ecut_pawecutdg(self, accuracy):
    #"""Return the maximum value of ecut and pawecutdg based on the hints available in the pseudos."""
    #    ecut = max(p.hint_for_accuracy(accuracy=accuracy).ecut for p in self)
    #    pawecutdg = max(p.hint_for_accuracy(accuracy=accuracy).pawecutdg for p in self)
    #    return ecut, pawecutdg

    def as_dict(self, **kwargs):
        d = {}
        for p in self:
            k, count = p.element.name, 1
            # k, count = p.element, 1
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
        return list(product(*d.values()))

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
        return self.pseudos_with_symbols(structure.symbol_set)

    def print_table(self, stream=sys.stdout, filter_function=None):
        """
        A pretty ASCII printer for the periodic table, based on some filter_function.

        Args:
            stream: file-like object
            filter_function:
                A filtering function that take a Pseudo as input and returns a boolean.
                For example, setting filter_function = lambda p: p.Z_val > 2 will print
                a periodic table containing only pseudos with Z_val > 2.
        """
        print(self.to_table(filter_function=filter_function), file=stream)

    def to_table(self, filter_function=None):
        """Return string with data in tabular form."""
        table = []
        for p in self:
            if filter_function is not None and filter_function(p): continue
            table.append([p.basename, p.symbol, p.Z_val, p.l_max, p.l_local, p.xc, p.type])
        return tabulate(table, headers= ["basename", "symbol", "Z_val", "l_max", "l_local", "XC", "type"],
                        tablefmt="grid")

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
