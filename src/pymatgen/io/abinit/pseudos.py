"""
This module provides objects describing the basic parameters of the
pseudopotentials used in Abinit, and a parser to instantiate pseudopotential objects..
"""

from __future__ import annotations

import abc
import collections
import hashlib
import logging
import os
import shutil
import sys
import tempfile
import traceback
from collections import defaultdict
from typing import TYPE_CHECKING, NamedTuple
from xml.etree import ElementTree as ET

import numpy as np
from monty.collections import AttrDict
from monty.functools import lazy_property
from monty.itertools import iterator_from_slice
from monty.json import MontyDecoder, MSONable
from monty.os.path import find_exts
from tabulate import tabulate

from pymatgen.core import Element
from pymatgen.core.xcfunc import XcFunc
from pymatgen.io.core import ParseError
from pymatgen.util.plotting import add_fig_kwargs, get_ax_fig

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence
    from typing import Any, ClassVar

    import matplotlib.pyplot as plt
    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.core import Structure

logger = logging.getLogger(__name__)


__author__ = "Matteo Giantomassi"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


# Tools and helper functions.


def _read_nlines(filename: str, n_lines: int) -> list[str]:
    """Read at most nlines from filename. If nlines is < 0, the entire file is read."""
    if n_lines < 0:
        with open(filename, encoding="utf-8") as file:
            return file.readlines()

    lines = []
    with open(filename, encoding="utf-8") as file:
        for lineno, line in enumerate(file):
            if lineno == n_lines:
                break
            lines.append(line)
        return lines


_l2str = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i"}

_str2l = {v: k for k, v in _l2str.items()}


def l2str(l_ang_mom):
    """Convert the angular momentum l (int) to string."""
    try:
        return _l2str[l_ang_mom]
    except KeyError:
        return f"Unknown angular momentum, received {l_ang_mom = }"


def str2l(s):
    """Convert a string to the angular momentum l (int)."""
    return _str2l[s]


class Pseudo(MSONable, abc.ABC):
    """
    Abstract base class defining the methods that must be
    implemented by the concrete pseudo-potential sub-classes.
    """

    @classmethod
    def as_pseudo(cls, obj: Pseudo | str) -> Pseudo:
        """Convert obj into a Pseudo.

        Args:
            obj (str | Pseudo): Path to the pseudo file or a Pseudo object.
        """
        return obj if isinstance(obj, cls) else cls.from_file(obj)  # type:ignore[arg-type]

    @classmethod
    def from_file(cls, filename: str) -> Self:
        """
        Build an instance of a concrete Pseudo subclass from filename.
        Note: the parser knows the concrete class that should be instantiated
        Client code should rely on the abstract interface provided by Pseudo.
        """
        return PseudoParser().parse(filename)

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("md5", "Z", "Z_val", "l_max")
        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented
        return (
            all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)
            and self.__class__ == other.__class__
        )

    def __repr__(self) -> str:
        try:
            return f"<{type(self).__name__} at {os.path.relpath(self.filepath)}>"
        except Exception:
            # relpath can fail if the code is executed in demon mode.
            return f"<{type(self).__name__} at {self.filepath}>"

    def __str__(self) -> str:
        return self.to_str()

    def to_str(self, verbose=0) -> str:
        """String representation."""
        lines: list[str] = []
        lines += (
            f"<{type(self).__name__}: {self.basename}>",
            f"  summary: {self.summary.strip()}",
            f"  number of valence electrons: {self.Z_val}",
            f"  maximum angular momentum: {l2str(self.l_max)}",
            f"  angular momentum for local part: {l2str(self.l_local)}",
            f"  XC correlation: {self.xc}",
            f"  supports spin-orbit: {self.supports_soc}",
        )

        if self.isnc:
            lines.append(f"  radius for non-linear core correction: {self.nlcc_radius}")

        if self.has_hints:
            for accuracy in ("low", "normal", "high"):
                hint = self.hint_for_accuracy(accuracy=accuracy)
                lines.append(f"  hint for {accuracy} accuracy: {hint}")

        return "\n".join(lines)

    @property
    @abc.abstractmethod
    def summary(self) -> str:
        """String summarizing the most important properties."""

    @property
    def filepath(self) -> str:
        """Absolute path to pseudopotential file."""
        return os.path.abspath(self.path)

    @property
    def basename(self) -> str:
        """File basename."""
        return os.path.basename(self.filepath)

    @property
    @abc.abstractmethod
    def Z(self) -> int:
        """The atomic number of the atom."""

    @property
    @abc.abstractmethod
    def Z_val(self) -> int:
        """Valence charge."""

    @property
    def type(self) -> str:
        """Type of pseudo."""
        return type(self).__name__

    @property
    def element(self) -> Element:
        """Pymatgen Element."""
        try:
            return Element.from_Z(self.Z)
        except (KeyError, IndexError):
            return Element.from_Z(int(self.Z))

    @property
    def symbol(self) -> str:
        """Element symbol."""
        return self.element.symbol

    @property
    @abc.abstractmethod
    def l_max(self) -> int:
        """Maximum angular momentum."""

    @property
    @abc.abstractmethod
    def l_local(self) -> int:
        """Angular momentum used for the local part."""

    @property
    def isnc(self) -> bool:
        """True if norm-conserving pseudopotential."""
        return isinstance(self, NcPseudo)

    @property
    def ispaw(self) -> bool:
        """True if PAW pseudopotential."""
        return isinstance(self, PawPseudo)

    @lazy_property
    def md5(self):
        """MD5 hash value."""
        # if self.has_dojo_report and "md5" in self.dojo_report: return self.dojo_report["md5"]
        return self.compute_md5()

    def compute_md5(self):
        """Compute and return MD5 hash value."""
        with open(self.path, encoding="utf-8") as file:
            text = file.read()
            # usedforsecurity=False needed in FIPS mode (Federal Information Processing Standards)
            # https://github.com/materialsproject/pymatgen/issues/2804
            md5 = hashlib.md5(usedforsecurity=False)
            md5.update(text.encode("utf-8"))
            return md5.hexdigest()

    @property
    @abc.abstractmethod
    def supports_soc(self):
        """
        True if the pseudo can be used in a calculation with spin-orbit coupling.
        Base classes should provide a concrete implementation that computes this value.
        """

    def as_dict(self, **kwargs):
        """Return dictionary for MSONable protocol."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "basename": self.basename,
            "type": self.type,
            "symbol": self.symbol,
            "Z": self.Z,
            "Z_val": self.Z_val,
            "l_max": self.l_max,
            "md5": self.md5,
            "filepath": self.filepath,
            # "xc": self.xc.as_dict(),
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build instance from dictionary (MSONable protocol)."""
        new = cls.from_file(dct["filepath"])

        # Consistency test based on md5
        if dct.get("md5") != new.md5:
            raise ValueError(
                f"The md5 found in file does not agree with the one in dict\nReceived {dct['md5']}\nComputed {new.md5}"
            )

        return new

    def as_tmpfile(self, tmpdir=None):
        """
        Copy the pseudopotential to a temporary a file and returns a new pseudopotential object.
        Useful for unit tests in which we have to change the content of the file.

        Args:
            tmpdir: If None, a new temporary directory is created and files are copied here
                else tmpdir is used.
        """
        tmpdir = tempfile.mkdtemp() if tmpdir is None else tmpdir
        new_path = os.path.join(tmpdir, self.basename)
        shutil.copy(self.filepath, new_path)

        # Copy dojo report file if present.
        root, _ext = os.path.splitext(self.filepath)
        dj_report = f"{root}.djrepo"
        if os.path.isfile(dj_report):
            shutil.copy(dj_report, os.path.join(tmpdir, os.path.basename(dj_report)))

        # Build new object and copy dojo_report if present.
        new = type(self).from_file(new_path)
        if self.has_dojo_report:
            new.dojo_report = self.dojo_report.deepcopy()

        return new

    @property
    def has_dojo_report(self):
        """True if the pseudo has an associated `DOJO_REPORT` section."""
        return hasattr(self, "dojo_report") and self.dojo_report

    @property
    def djrepo_path(self):
        """The path of the djrepo file. None if file does not exist."""
        root, _ext = os.path.splitext(self.filepath)
        return f"{root}.djrepo"
        # if os.path.isfile(path): return path
        # return None

    def hint_for_accuracy(self, accuracy="normal"):
        """Get a Hint object with the suggested value of ecut [Ha] and
        pawecutdg [Ha] for the given accuracy.
        ecut and pawecutdg are set to zero if no hint is available.

        Args:
            accuracy: ["low", "normal", "high"]
        """
        if not self.has_dojo_report:
            return Hint(ecut=0.0, pawecutdg=0.0)

        # Get hints from dojoreport. Try first in hints then in ppgen_hints.
        if "hints" in self.dojo_report:
            return Hint.from_dict(self.dojo_report["hints"][accuracy])
        if "ppgen_hints" in self.dojo_report:
            return Hint.from_dict(self.dojo_report["ppgen_hints"][accuracy])
        return Hint(ecut=0.0, pawecutdg=0.0)

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
        pseudopotential part. Returns PspsFile object providing methods
        to plot and analyze the data or None if file is not found or it's not readable.

        Args:
            ecut: Cutoff energy in Hartree.
            pawecutdg: Cutoff energy for the PAW double grid.
        """
        from abipy.abio.factories import gs_input
        from abipy.core.structure import Structure
        from abipy.electrons.psps import PspsFile
        from abipy.flowtk import AbinitTask

        # Build fake structure.
        lattice = 10 * np.eye(3)
        structure = Structure(lattice, [self.element], coords=[[0, 0, 0]])

        if self.ispaw and pawecutdg is None:
            pawecutdg = ecut * 4
        inp = gs_input(
            structure,
            pseudos=[self],
            ecut=ecut,
            pawecutdg=pawecutdg,
            spin_mode="unpolarized",
            kppa=1,
        )
        # Add prtpsps = -1 to make Abinit print the PSPS.nc file and stop.
        inp["prtpsps"] = -1

        # Build temporary task and run it (ignore retcode because we don't exit cleanly)
        task = AbinitTask.temp_shell_task(inp)
        task.start_and_wait()

        filepath = task.outdir.has_abiext("_PSPS.nc")
        if not filepath:
            logger.critical(f"Cannot find PSPS.nc file in {task.outdir}")
            return None

        # Open the PSPS.nc file.
        try:
            return PspsFile(filepath)
        except Exception as exc:
            logger.critical(f"Exception while reading PSPS file at {filepath}:\n{exc}")
            return None


class NcPseudo(abc.ABC):
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


class PawPseudo(abc.ABC):
    """
    Abstract class that defines the methods that must be implemented
    by the concrete classes representing PAW pseudopotentials.
    """

    # def nlcc_radius(self):
    #    """
    #    Radius at which the core charge vanish (i.e. cut-off in a.u.).
    #    Returns 0.0 if nlcc is not used.
    #    """
    #    return 0.0
    #

    # @property
    # def has_nlcc(self):
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
    """An AbinitPseudo is a pseudopotential whose file contains an abinit header."""

    def __init__(self, path, header):
        """
        Args:
            path: Filename.
            header: AbinitHeader instance.
        """
        self.path = path
        self.header = header
        self._summary = header.summary

        # Build xc from header.
        self.xc = XcFunc.from_abinit_ixc(header["pspxc"])

        for attr_name in header:
            value = header.get(attr_name)

            # Hide these attributes since one should always use the public interface.
            setattr(self, f"_{attr_name}", value)

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
        # Treat ONCVPSP pseudos

        if self._pspcod == 8:
            switch = self.header["extension_switch"]
            if switch in (0, 1):
                return False
            if switch in (2, 3):
                return True
            raise ValueError(f"Don't know how to handle extension_switch: {switch}")

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

    # def orbitals(self):

    @property
    def supports_soc(self):
        return True


class Hint:
    """
    Suggested value for the cutoff energy [Hartree units]
    and the cutoff energy for the dense grid (only for PAW pseudos).
    """

    def __init__(self, ecut, pawecutdg=None):
        self.ecut = ecut
        self.pawecutdg = ecut if pawecutdg is None else pawecutdg

    def __str__(self):
        if self.pawecutdg is not None:
            return f"ecut: {self.ecut}, pawecutdg: {self.pawecutdg}"
        return f"ecut: {self.ecut}"

    def as_dict(self):
        """Return dictionary for MSONable protocol."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "ecut": self.ecut,
            "pawecutdg": self.pawecutdg,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build instance from dictionary (MSONable protocol)."""
        return cls(**{k: v for k, v in dct.items() if not k.startswith("@")})


def _dict_from_lines(lines, key_nums, sep=None) -> dict:
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
    if isinstance(lines, str):
        lines = [lines]

    if not isinstance(key_nums, collections.abc.Iterable):
        key_nums = list(key_nums)

    if len(lines) != len(key_nums):
        raise ValueError(f"{lines = }\n{key_nums = }")

    # TODO: PR 4223: kwargs was using `monty.collections.Namespace`,
    # revert to original implementation if needed
    kwargs: dict = {}

    for idx, nk in enumerate(key_nums):
        if nk == 0:
            continue
        line = lines[idx]

        tokens = [tok.strip() for tok in line.split()]
        values = tokens[:nk]
        # Sanitize keys: In some case we might get strings in the form: foo[,bar]
        keys = "".join(tokens[nk:]).replace("[", "").replace("]", "").split(",")

        if sep is not None:
            check = keys[0][0]
            if check != sep:
                raise ValueError(f"Expecting separator {sep}, got {check}")
            keys[0] = keys[0][1:]

        if len(values) != len(keys):
            raise ValueError(f"{line=}\n {len(keys)=} must equal {len(values)=}")

        kwargs.update(zip(keys, values, strict=True))

    return kwargs


class AbinitHeader(dict):
    """Dictionary whose keys can be also accessed as attributes."""

    def __getattr__(self, name):
        try:
            return super().__getattribute__(name)  # this is just default behavior
        except AttributeError:
            try:
                return self[name]  # if above failed, try the dictionary
            except KeyError as exc:
                raise AttributeError(str(exc))


def _int_from_str(string):
    """
    Convert string into integer.

    Raise:
        TypeError if string is not a valid integer
    """
    float_num = float(string)
    int_num = int(float_num)
    if float_num == int_num:
        return int_num
    # Needed to handle pseudos with fractional charge
    int_num = np.rint(float_num)
    logger.warning(f"Converting float {float_num} to int {int_num}")
    return int_num


class NcAbinitHeader(AbinitHeader):
    """The abinit header found in the NC pseudopotential files."""

    _VARS: ClassVar[dict[str, tuple]] = {
        "zatom": (None, _int_from_str),
        "zion": (None, float),
        "pspdat": (None, float),
        "pspcod": (None, int),
        "pspxc": (None, int),
        "lmax": (None, int),
        "lloc": (None, int),
        "r2well": (None, float),
        "mmax": (None, float),
        "rchrg": (0.0, float),
        "fchrg": (0.0, float),
        "qchrg": (0.0, float),
    }

    def __init__(self, summary, **kwargs):
        super().__init__()

        # pseudos generated by APE use llocal instead of lloc.
        if "llocal" in kwargs:
            kwargs["lloc"] = kwargs.pop("llocal")

        self.summary = summary.strip()

        for key, desc in NcAbinitHeader._VARS.items():
            default, astype = desc
            value = kwargs.pop(key, None)

            if value is None:
                value = default
                if default is None:
                    raise RuntimeError(f"Attribute {key} must be specified")
            else:
                try:
                    value = astype(value)
                except Exception:
                    raise RuntimeError(f"Conversion Error for {key=}, {value=}")

            self[key] = value

        # Add remaining arguments, e.g. extension_switch
        if kwargs:
            self.update(kwargs)

    @staticmethod
    def fhi_header(filename, ppdesc):
        """Parse the FHI abinit header. Example:

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
        """Parse the HGH abinit header. Example:

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
        """Parse the GTH abinit header. Example:

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
        """Parse the ONCVPSP abinit header. Example:

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
        header["pspdat"] = header["pspd"]
        header.pop("pspd")

        # Read extension switch
        header["extension_switch"] = int(lines[5].split()[0])

        return NcAbinitHeader(summary, **header)

    @staticmethod
    def tm_header(filename, ppdesc):
        """Parse the TM abinit header. Example:

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
        lmax = None

        for lineno, line in enumerate(lines):
            header.append(line)
            if lineno == 2:
                # Read lmax
                tokens = line.split()
                _pspcod, _pspxc, lmax, _lloc = map(int, tokens[:4])
                _mmax, _r2well = map(float, tokens[4:6])
                # if tokens[-1].strip() != "pspcod,pspxc,lmax,lloc,mmax,r2well":
                #    raise RuntimeError("%s: Invalid line\n %s"  % (filename, line))

                lines = lines[3:]
                break

        # TODO
        # Parse the section with the projectors.
        # 0   4.085   6.246    0   2.8786493        l,e99.0,e99.9,nproj,rcpsp
        # .00000000    .0000000000    .0000000000    .00000000   rms,ekb1,ekb2,epsatm
        projectors = {}
        proj_info = []
        idx = None
        for idx in range(2 * (lmax + 1)):
            line = lines[idx]
            if idx % 2 == 0:
                proj_info = [
                    line,
                ]
            else:
                proj_info.append(line)
                d = _dict_from_lines(proj_info, [5, 4])
                projectors[int(d["l"])] = d

        # Add the last line with info on nlcc.
        header.append(lines[idx + 1])
        summary = header[0]

        header = _dict_from_lines(header, [0, 3, 6, 3])

        return NcAbinitHeader(summary, **header)


class PawAbinitHeader(AbinitHeader):
    """The abinit header found in the PAW pseudopotential files."""

    _VARS: ClassVar[dict[str, tuple]] = {
        "zatom": (None, _int_from_str),
        "zion": (None, float),
        "pspdat": (None, float),
        "pspcod": (None, int),
        "pspxc": (None, int),
        "lmax": (None, int),
        "lloc": (None, int),
        "mmax": (None, int),
        "r2well": (None, float),
        "pspfmt": (None, str),
        "creatorID": (None, int),
        "basis_size": (None, int),
        "lmn_size": (None, int),
        "orbitals": (None, list),
        "number_of_meshes": (None, int),
        "r_cut": (None, float),  # r_cut(PAW) in the header
        "shape_type": (None, int),
        "rshape": (None, float),
    }

    def __init__(self, summary, **kwargs):
        super().__init__()

        self.summary = summary.strip()

        for key, desc in self._VARS.items():
            default, astype = desc

            value = kwargs.pop(key, None)

            if value is None:
                value = default
                if default is None:
                    raise RuntimeError(f"Attribute {key} must be specified")
            else:
                try:
                    value = astype(value)
                except Exception:
                    raise RuntimeError(f"Conversion Error for {key=}, with {value=}")

            self[key] = value

        if kwargs:
            raise RuntimeError(f"kwargs should be empty but got {kwargs}")

    @staticmethod
    def paw_header(filename, ppdesc):
        """Parse the PAW abinit header. Examples:

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
            raise NotImplementedError(f"format {ppdesc.format} not in {supported_formats}")

        lines = _read_nlines(filename, -1)

        summary = lines[0]
        header = _dict_from_lines(lines[:5], [0, 3, 6, 2, 2], sep=":")

        lines = lines[5:]
        # TODO
        # Parse orbitals and number of meshes.
        header["orbitals"] = [int(tok) for tok in lines[0].split(":")[0].split()]
        header["number_of_meshes"] = num_meshes = int(lines[1].split(":")[0])

        # Skip meshes =
        lines = lines[2 + num_meshes :]
        # for midx in range(num_meshes):
        #    l = midx + 1

        header["r_cut"] = float(lines[0].split(":")[0])
        header.update(_dict_from_lines(lines[1], [2], sep=":"))

        return PawAbinitHeader(summary, **header)


class PseudoParseError(ParseError):
    """Base Error class for the exceptions raised by PseudoParser."""


class PseudoParser:
    """
    Responsible for parsing pseudopotential files and returning pseudopotential objects.

    Usage:
        pseudo = PseudoParser().parse("filename")
    """

    Error = PseudoParseError

    class ppdesc(NamedTuple):
        """Supported values of pspcod."""

        pspcod: int
        name: str
        psp_type: str
        format: None

    # TODO Recheck
    _PSPCODES: ClassVar[dict[int, ppdesc]] = {
        1: ppdesc(1, "TM", "NC", None),
        2: ppdesc(2, "GTH", "NC", None),
        3: ppdesc(3, "HGH", "NC", None),
        4: ppdesc(4, "Teter", "NC", None),
        # 5: ppdesc(5, "NC",     , None),
        6: ppdesc(6, "FHI", "NC", None),
        7: ppdesc(6, "PAW_abinit_text", "PAW", None),
        8: ppdesc(8, "ONCVPSP", "NC", None),
        10: ppdesc(10, "HGHK", "NC", None),
    }

    # renumber functionals from oncvpsp todo confirm that 3 is 2
    # _FUNCTIONALS = {1: {'n': 4, 'name': 'Wigner'},
    #                2: {'n': 5, 'name': 'HL'},
    #                3: {'n': 2, 'name': 'PWCA'},
    #                4: {'n': 11, 'name': 'PBE'}}

    def __init__(self):
        # List of files that have been parsed successfully.
        self._parsed_paths: list = []

        # List of files that could not been parsed.
        self._wrong_paths: list = []

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
                exclude_exts[i] = f".{ext.strip()}"

        # Exclude files depending on the extension.
        paths = []
        for fname in os.listdir(dirname):
            _root, ext = os.path.splitext(fname)
            path = os.path.join(dirname, fname)
            if ext in exclude_exts or fname in exclude_fnames or fname.startswith(".") or not os.path.isfile(path):
                continue
            paths.append(path)

        pseudos = []
        for path in paths:
            # Parse the file and generate the pseudo.
            try:
                pseudo = self.parse(path)
            except Exception:
                pseudo = None

            if pseudo is not None:
                pseudos.append(pseudo)
                self._parsed_paths.extend(path)
            else:
                self._wrong_paths.extend(path)

        return pseudos

    def read_ppdesc(self, filename):
        """
        Read the pseudopotential descriptor from filename.

        Returns:
            Pseudopotential descriptor. None if filename is not a valid pseudopotential file.

        Raises:
            `PseudoParseError` if fileformat is not supported.
        """
        if filename.endswith(".xml"):
            raise self.Error("XML pseudo not supported yet")

        # Assume file with the abinit header.
        lines = _read_nlines(filename, 80)

        for lineno, line in enumerate(lines):
            if lineno == 2:
                try:
                    tokens = line.split()
                    pspcod, _pspxc = map(int, tokens[:2])
                except Exception:
                    msg = f"{filename}: Cannot parse pspcod, pspxc in line\n {line}"
                    logger.critical(msg)
                    return None

                if pspcod not in self._PSPCODES:
                    raise self.Error(f"{filename}: Don't know how to handle {pspcod=}\n")

                ppdesc = self._PSPCODES[pspcod]

                if pspcod == 7:
                    # PAW -> need to know the format pspfmt
                    tokens = lines[lineno + 1].split()
                    pspfmt, _creatorID = tokens[:2]

                    ppdesc = ppdesc._replace(format=pspfmt)

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
            logger.critical(f"Cannot find ppdesc in {path}")
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
            str_traceback = "\n".join((traceback.format_exc(), str(sys.exc_info()[0])))

            raise self.Error(f"{path}:\n{str_traceback}")

        if psp_type == "NC":
            pseudo = NcAbinitPseudo(path, header)
        elif psp_type == "PAW":
            pseudo = PawAbinitPseudo(path, header)
        else:
            raise NotImplementedError("psp_type not in [NC, PAW]")

        return pseudo


class RadialFunction(NamedTuple):
    """Radial Function class.

    TODO: use RadialFunction from pseudo_dojo.
    """

    mesh: Any
    values: NDArray


class PawXmlSetup(Pseudo, PawPseudo):
    """Setup class for PawXml."""

    def __init__(self, filepath):
        """
        Args:
            filepath (str): Path to the XML file.
        """
        self.path = os.path.abspath(filepath)

        # Get the XML root (this trick is used to that the object is pickleable).
        root = self.root

        # Get the version of the XML format
        self.paw_setup_version = root.get("version")

        # Info on the atom.
        atom_attrib = root.find("atom").attrib

        # self._symbol = atom_attrib["symbol"]
        self._zatom = int(float(atom_attrib["Z"]))
        self.core, self.valence = map(float, [atom_attrib["core"], atom_attrib["valence"]])

        # Build xc from header.
        xc_info = root.find("xc_functional").attrib
        self.xc = XcFunc.from_type_name(xc_info["type"], xc_info["name"])

        # Old XML files do not define this field!
        # In this case we set the PAW radius to None.
        # self._paw_radius = float(root.find("PAW_radius").attrib["rpaw"])

        # self.ae_energy = {k: float(v) for k,v in root.find("ae_energy").attrib.items()}
        pawr_element = root.find("PAW_radius")
        self._paw_radius = None
        if pawr_element is not None:
            self._paw_radius = float(pawr_element.attrib["rpaw"])

        # <valence_states>
        #  <state n="2" l="0" f="2"  rc="1.10" e="-0.6766" id="N-2s"/>
        #  <state n="2" l="1" f="3"  rc="1.10" e="-0.2660" id="N-2p"/>
        #  <state       l="0"        rc="1.10" e=" 0.3234" id="N-s1"/>
        #  <state       l="1"        rc="1.10" e=" 0.7340" id="N-p1"/>
        #  <state       l="2"        rc="1.10" e=" 0.0000" id="N-d1"/>
        # </valence_states>
        #
        # The valence_states element contains several state elements.
        # For this setup, the first two lines describe bound eigenstates
        # with occupation numbers and principal quantum numbers.
        # Notice, that the three additional unbound states should have no f and n attributes.
        # In this way, we know that only the first two bound states (with f and n attributes)
        # should be used for constructing an initial guess for the wave functions.

        self.valence_states: dict = {}
        for node in root.find("valence_states"):
            attrib = AttrDict(node.attrib)
            if attrib.id in self.valence_states:
                raise ValueError(f"{attrib.id=} should not be in {self.valence_states=}")
            self.valence_states[attrib.id] = attrib

        # Parse the radial grids
        self.rad_grids: dict = {}
        for node in root.findall("radial_grid"):
            grid_params = node.attrib
            gid = grid_params["id"]
            if gid in self.rad_grids:
                raise ValueError(f"{gid=} should not be in {self.rad_grids=}")

            self.rad_grids[gid] = self._eval_grid(grid_params)

    def __getstate__(self):
        """Get state is pickled as the contents for the instance.

        In this case we just remove the XML root element process since Element object cannot be pickled.
        """
        return {k: v for k, v in self.__dict__.items() if k != "_root"}

    @lazy_property
    def root(self):
        """Root tree of XML."""
        tree = ET.parse(self.filepath)
        return tree.getroot()

    @property
    def Z(self):
        return self._zatom

    @property
    def Z_val(self):
        """Number of valence electrons."""
        return self.valence

    @property
    def l_max(self):
        """Maximum angular momentum."""
        # TODO return an actual value
        return

    @property
    def l_local(self):
        """Angular momentum used for the local part."""
        return

    @property
    def summary(self):
        """String summarizing the most important properties."""
        return ""

    @property
    def paw_radius(self):
        return self._paw_radius

    @property
    def supports_soc(self):
        """Here I assume that the ab-initio code can treat the SOC within the on-site approximation."""
        return True

    @staticmethod
    def _eval_grid(grid_params):
        """For a dictionary with the parameters defining the
        radial mesh, get a `ndarray` with the mesh.
        """
        eq = grid_params.get("eq").replace(" ", "")
        istart, iend = int(grid_params.get("istart")), int(grid_params.get("iend"))
        indices = list(range(istart, iend + 1))

        if eq == "r=a*exp(d*i)":
            a, d = float(grid_params["a"]), float(grid_params["d"])
            mesh = [a * np.exp(d * i) for i in indices]

        elif eq == "r=a*i/(n-i)":
            a, n = float(grid_params["a"]), float(grid_params["n"])
            mesh = [a * i / (n - i) for i in indices]

        elif eq == "r=a*(exp(d*i)-1)":
            a, d = float(grid_params["a"]), float(grid_params["d"])
            mesh = [a * (np.exp(d * i) - 1.0) for i in indices]

        elif eq == "r=d*i":
            d = float(grid_params["d"])
            mesh = [d * i for i in indices]

        elif eq == "r=(i/n+a)^5/a-a^4":
            a, n = float(grid_params["a"]), float(grid_params["n"])
            mesh = [(i / n + a) ** 5 / a - a**4 for i in indices]

        else:
            raise ValueError(f"Unknown grid type: {eq}")

        return np.array(mesh)

    def _parse_radfunc(self, func_name):
        """Parse the first occurrence of func_name in the XML file."""
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

    @lazy_property
    def ae_core_density(self):
        """The all-electron radial density."""
        mesh, values, _attrib = self._parse_radfunc("ae_core_density")
        return RadialFunction(mesh, values)

    @lazy_property
    def pseudo_core_density(self):
        """The pseudized radial density."""
        mesh, values, _attrib = self._parse_radfunc("pseudo_core_density")
        return RadialFunction(mesh, values)

    @lazy_property
    def ae_partial_waves(self):
        """Dictionary with the AE partial waves indexed by state."""
        ae_partial_waves = {}
        for mesh, values, attrib in self._parse_all_radfuncs("ae_partial_wave"):
            state = attrib["state"]
            # val_state = self.valence_states[state]
            ae_partial_waves[state] = RadialFunction(mesh, values)

        return ae_partial_waves

    @property
    def pseudo_partial_waves(self):
        """Dictionary with the pseudo partial waves indexed by state."""
        pseudo_partial_waves = {}
        for mesh, values, attrib in self._parse_all_radfuncs("pseudo_partial_wave"):
            state = attrib["state"]
            # val_state = self.valence_states[state]
            pseudo_partial_waves[state] = RadialFunction(mesh, values)

        return pseudo_partial_waves

    @lazy_property
    def projector_functions(self):
        """Dictionary with the PAW projectors indexed by state."""
        projector_functions = {}
        for mesh, values, attrib in self._parse_all_radfuncs("projector_function"):
            state = attrib["state"]
            # val_state = self.valence_states[state]
            projector_functions[state] = RadialFunction(mesh, values)

        return projector_functions

    def yield_figs(self, **kwargs):  # pragma: no cover
        """This function *generates* a predefined list of matplotlib figures with minimal input from the user."""
        yield self.plot_densities(title="PAW densities", show=False)
        yield self.plot_waves(title="PAW waves", show=False)
        yield self.plot_projectors(title="PAW projectors", show=False)
        # yield self.plot_potentials(title="potentials", show=False)

    @add_fig_kwargs
    def plot_densities(self, ax: plt.Axes = None, **kwargs):
        """
        Plot the PAW densities.

        Args:
            ax: matplotlib Axes or None if a new figure should be created.

        Returns:
            `matplotlib` figure
        """
        ax, fig = get_ax_fig(ax)

        ax.grid(visible=True)
        ax.set_xlabel("r [Bohr]")
        # ax.set_ylabel('density')

        for idx, density_name in enumerate(["ae_core_density", "pseudo_core_density"]):
            rden = getattr(self, density_name)
            label = "$n_c$" if idx == 1 else r"$\tilde{n}_c$"
            ax.plot(rden.mesh, rden.mesh * rden.values, label=label, lw=2)

        ax.legend(loc="best")

        return fig

    @add_fig_kwargs
    def plot_waves(self, ax: plt.Axes = None, fontsize=12, **kwargs):
        """
        Plot the AE and the pseudo partial waves.

        Args:
            ax: matplotlib Axes or None if a new figure should be created.
            fontsize: fontsize for legends and titles

        Returns:
            plt.Figure: matplotlib figure
        """
        ax, fig = get_ax_fig(ax)

        ax.grid(visible=True)
        ax.set_xlabel("r [Bohr]")
        ax.set_ylabel(r"$r\phi,\, r\tilde\phi\, [Bohr]^{-\frac{1}{2}}$")

        # ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        # ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.pseudo_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label=f"PS-WAVE: {state}")

        for state, rfunc in self.ae_partial_waves.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, lw=2, label=f"AE-WAVE: {state}")

        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    @add_fig_kwargs
    def plot_projectors(self, ax: plt.Axes = None, fontsize=12, **kwargs):
        """
        Plot the PAW projectors.

        Args:
            ax: matplotlib Axes or None if a new figure should be created.

        Returns:
            plt.Figure: matplotlib figure
        """
        ax, fig = get_ax_fig(ax)
        ax.grid(visible=True)
        ax.set_xlabel("r [Bohr]")
        ax.set_ylabel(r"$r\tilde p\, [Bohr]^{-\frac{1}{2}}$")

        # ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
        # ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

        for state, rfunc in self.projector_functions.items():
            ax.plot(rfunc.mesh, rfunc.mesh * rfunc.values, label=f"TPROJ: {state}")

        ax.legend(loc="best", shadow=True, fontsize=fontsize)

        return fig

    # @add_fig_kwargs
    # def plot_potentials(self, **kwargs):
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

    #    fig = plt.figure()

    #    ax = fig.add_subplot(1,1,1)
    #    ax.grid(visible=True)
    #    ax.set_xlabel('r [Bohr]')
    #    ax.set_ylabel('density')
    #    ax.axvline(x=self.paw_radius, linewidth=2, color='k', linestyle="--")
    #    ax.annotate("$r_c$", xy=(self.paw_radius + 0.1, 0.1))

    #    for state, rfunc in self.potentials.items():
    #        ax.plot(rfunc.mesh, rfunc.values, label=f"TPROJ: {state}")

    #    ax.legend(loc="best")

    #    if title is not None: fig.suptitle(title)
    #    if show: plt.show()
    #    if savefig: fig.savefig(savefig)
    #    return fig


class PseudoTable(collections.abc.Sequence, MSONable):
    """
    Define the pseudopotentials from the element table.
    Individidual elements are accessed by name, symbol or atomic number.

    For example, the following all retrieve iron:

    print(elements[26])
    print(elements.Fe)
    print(elements.symbol('Fe'))
    print(elements.name('iron'))
    print(elements.isotope('Fe'))
    """

    @classmethod
    def as_table(cls, items):
        """Return an instance of PseudoTable from the iterable items."""
        return items if isinstance(items, cls) else cls(items)

    @classmethod
    def from_dir(cls, top, exts=None, exclude_dirs="_*") -> Self | None:
        """Find all pseudos in the directory tree starting from top.

        Args:
            top: Top of the directory tree
            exts: List of files extensions. if exts == "all_files"
                    we try to open all files in top
            exclude_dirs: Wildcard used to exclude directories.

        Returns:
            PseudoTable sorted by atomic number Z.
        """
        pseudos = []

        if exts == "all_files":
            for filepath in [os.path.join(top, fn) for fn in os.listdir(top)]:
                if os.path.isfile(filepath):
                    try:
                        if pseudo := Pseudo.from_file(filepath):
                            pseudos.append(pseudo)
                        else:
                            logger.info(f"Skipping file {filepath}")
                    except Exception:
                        logger.info(f"Skipping file {filepath}")
            if not pseudos:
                logger.warning(f"No pseudopotentials parsed from folder {top}")
                return None
            logger.info(f"Creating PseudoTable with {len(pseudos)} pseudopotentials")

        else:
            if exts is None:
                exts = ("psp8",)

            for pseudo in find_exts(top, exts, exclude_dirs=exclude_dirs):  # type:ignore[assignment]
                try:
                    pseudos.append(Pseudo.from_file(pseudo))  # type:ignore[arg-type]
                except Exception as exc:
                    logger.critical(f"Error in {pseudo}:\n{exc}")

        return cls(pseudos).sort_by_z()

    def __init__(self, pseudos: Sequence[Pseudo]) -> None:
        """
        Args:
            pseudos: List of pseudopotentials or filepaths.
        """
        # Store pseudos in a default dictionary with z as key.
        # Note that we can have more than one pseudo for given z.
        # hence the values are lists of pseudos.
        if not isinstance(pseudos, collections.abc.Iterable):
            pseudos = [pseudos]

        if isinstance(pseudos, str):
            pseudos = [pseudos]

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
                raise ValueError(f"All symbols must be equal while they are: {symbols}")

            setattr(self, symbol, pseudo_list)

    def __getitem__(self, Z):
        """Retrieve pseudos for the atomic number z. Accepts both int and slice objects."""
        if isinstance(Z, slice):
            if Z.stop is None:
                raise ValueError("Z.stop is None")
            pseudos = []
            for znum in iterator_from_slice(Z):
                pseudos.extend(self._pseudos_with_z[znum])
            return type(self)(pseudos)
        return type(self)(self._pseudos_with_z[Z])

    def __len__(self) -> int:
        return len(list(iter(self)))

    def __iter__(self) -> Iterator[Pseudo]:
        """Process the elements in Z order."""
        for z in self.zlist:
            yield from self._pseudos_with_z[z]

    def __repr__(self) -> str:
        return f"<{type(self).__name__} at {id(self)}>"

    def __str__(self) -> str:
        return self.to_table()

    @property
    def allnc(self) -> bool:
        """True if all pseudos are norm-conserving."""
        return all(p.isnc for p in self)

    @property
    def allpaw(self):
        """True if all pseudos are PAW."""
        return all(p.ispaw for p in self)

    @property
    def zlist(self):
        """Ordered list with the atomic numbers available in the table."""
        return sorted(self._pseudos_with_z)

    # def max_ecut_pawecutdg(self, accuracy):
    # """Return the maximum value of ecut and pawecutdg based on the hints available in the pseudos."""
    #    ecut = max(p.hint_for_accuracy(accuracy=accuracy).ecut for p in self)
    #    pawecutdg = max(p.hint_for_accuracy(accuracy=accuracy).pawecutdg for p in self)
    #    return ecut, pawecutdg

    def as_dict(self, **kwargs):
        """Return dictionary for MSONable protocol."""
        dct = {}
        for p in self:
            k, count = p.element.name, 1
            # k, count = p.element, 1
            # Handle multiple-pseudos with the same name!
            while k in dct:
                k += f"{k.split('#')[0]}#{count}"
                count += 1
            dct[k] = p.as_dict()
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Build instance from dictionary (MSONable protocol)."""
        pseudos = []
        for k, v in dct.items():
            if not k.startswith("@"):
                pseudos.append(MontyDecoder().process_decoded(v))
        return cls(pseudos)

    def is_complete(self, zmax=118) -> bool:
        """True if table is complete i.e. all elements with Z < zmax have at least on pseudopotential."""
        return all(self[z] for z in range(1, zmax))

    def all_combinations_for_elements(self, element_symbols):
        """Get a list with all the possible combination of pseudos
        for the given list of element_symbols.
        Each item is a list of pseudopotential objects.

        Example:
            table.all_combinations_for_elements(["Li", "F"])
        """
        dct = {}
        for symbol in element_symbols:
            dct[symbol] = self.select_symbols(symbol, ret_list=True)

        from itertools import product

        return list(product(*dct.values()))

    def pseudo_with_symbol(self, symbol, allow_multi=False):
        """Get the pseudo with the given chemical symbol.

        Args:
            symbols: String with the chemical symbol of the element
            allow_multi: By default, the method raises ValueError
                if multiple occurrences are found. Use allow_multi to prevent this.

        Raises:
            ValueError if symbol is not found or multiple occurrences are present and not allow_multi
        """
        pseudos = self.select_symbols(symbol, ret_list=True)
        if not pseudos or (len(pseudos) > 1 and not allow_multi):
            raise ValueError(f"Found {len(pseudos)} occurrences of {symbol=}")

        return pseudos if allow_multi else pseudos[0]

    def pseudos_with_symbols(self, symbols):
        """Get the pseudos with the given chemical symbols.

        Raises:
            ValueError if one of the symbols is not found or multiple occurrences are present.
        """
        pseudos = self.select_symbols(symbols, ret_list=True)
        found_symbols = [p.symbol for p in pseudos]

        if duplicated_elements := [s for s, o in collections.Counter(found_symbols).items() if o > 1]:
            raise ValueError(f"Found multiple occurrences of symbol(s) {', '.join(duplicated_elements)}")

        if missing_symbols := [s for s in symbols if s not in found_symbols]:
            raise ValueError(f"Missing data for symbol(s) {', '.join(missing_symbols)}")

        return pseudos

    def select_symbols(self, symbols, ret_list=False):
        """Get a PseudoTable with the pseudopotentials with the given list of chemical symbols.

        Args:
            symbols: str or list of symbols
                Prepend the symbol string with "-", to exclude pseudos.
            ret_list: if True a list of pseudos is returned instead of a PseudoTable
        """
        if isinstance(symbols, str):
            symbols = [symbols]

        exclude = symbols[0].startswith("-")

        if exclude:
            if not all(s.startswith("-") for s in symbols):
                raise ValueError("When excluding symbols, all strings must start with `-`")
            symbols = [s[1:] for s in symbols]

        symbols = set(symbols)
        pseudos = []
        for p in self:
            if exclude:
                if p.symbol in symbols:
                    continue
            elif p.symbol not in symbols:
                continue

            pseudos.append(p)

        if ret_list:
            return pseudos
        return type(self)(pseudos)

    def get_pseudos_for_structure(self, structure: Structure):
        """Get the list of Pseudo objects to be used for this Structure.

        Args:
            structure: pymatgen Structure.

        Raises:
            `ValueError` if one of the chemical symbols is not found or
            multiple occurrences are present in the table.
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
            if filter_function is not None and filter_function(p):
                continue
            table.append([p.basename, p.symbol, p.Z_val, p.l_max, p.l_local, p.xc, p.type])
        return tabulate(
            table,
            headers=["basename", "symbol", "Z_val", "l_max", "l_local", "XC", "type"],
            tablefmt="grid",
        )

    def sorted(self, attrname, reverse=False):
        """Sort the table according to the value of attribute attrname.

        Returns:
            New class: `PseudoTable` object
        """
        attrs = []
        for i, pseudo in self:
            try:
                a = getattr(pseudo, attrname)
            except AttributeError:
                a = np.inf
            attrs.append((i, a))

        # Sort attrs, and build new table with sorted pseudos.
        return type(self)([self[a[0]] for a in sorted(attrs, key=lambda t: t[1], reverse=reverse)])

    def sort_by_z(self):
        """Return a new PseudoTable with pseudos sorted by Z."""
        return type(self)(sorted(self, key=lambda p: p.Z))

    def select(self, condition) -> PseudoTable:
        """Select only those pseudopotentials for which condition is True.

        Args:
            condition:
                Function that accepts a Pseudo object and returns True or False.

        Returns:
            PseudoTable: New PseudoTable instance with pseudos for which condition is True.
        """
        return type(self)([p for p in self if condition(p)])

    def with_dojo_report(self):
        """Select pseudos containing the DOJO_REPORT section. Return new class:`PseudoTable` object."""
        return self.select(condition=lambda p: p.has_dojo_report)

    def select_rows(self, rows):
        """Get new class:`PseudoTable` object with pseudos in the given rows of the periodic table.
        rows can be either a int or a list of integers.
        """
        if not isinstance(rows, list | tuple):
            rows = [rows]
        return type(self)([p for p in self if p.element.row in rows])

    def select_family(self, family):
        """Return PseudoTable with element belonging to the specified family, e.g. family="alkaline"."""
        # e.g element.is_alkaline
        return type(self)([p for p in self if getattr(p.element, f"is_{family}")])
