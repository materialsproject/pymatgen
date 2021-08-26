# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides
"""

from collections import OrderedDict, namedtuple

from monty.functools import lazy_property
from monty.json import MSONable  # , MontyEncoder
from monty.string import is_string

from pymatgen.core.libxcfunc import LibxcFunc

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "3.0.0"  # The libxc version used to generate this file!
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo@gmail.com"
__status__ = "Production"
__date__ = "May 16, 2016"


class XcFunc(MSONable):
    """
    This object stores information about the XC correlation functional.
    Client code usually creates the object by calling the class methods:

        - from_name
        - from_type_name

    or code-specific methods such as:

        - from_abinit_ixc

    Ax XcFunc instance is hashable and can therefore be used as key in dictionaries.

    The implementation is based on the libxc conventions
    and is inspired to the XML specification for atomic PAW datasets documented at:

        https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html

    For convenience, part of the pawxml documentation is reported here.

    The xc_functional element defines the exchange-correlation functional used for
    generating the dataset. It has the two attributes type and name.

    The type attribute can be LDA, GGA, MGGA or HYB.
    The name attribute designates the exchange-correlation functional
    and can be specified in the following ways:

    [1] Taking the names from the LibXC library. The correlation and exchange names
        are stripped from their XC_ part and combined with a + sign.
        Here is an example for an LDA functional:

        <xc_functional type="LDA", name="LDA_X+LDA_C_PW"/>

        and this is what PBE will look like:

        <xc_functional type="GGA", name="GGA_X_PBE+GGA_C_PBE"/>

    [2] Using one of the following pre-defined aliases:

    type    name    LibXC equivalent             Reference
    LDA     PW      LDA_X+LDA_C_PW               LDA exchange; Perdew, Wang, PRB 45, 13244 (1992)
    GGA     PW91    GGA_X_PW91+GGA_C_PW91        Perdew et al PRB 46, 6671 (1992)
    GGA     PBE     GGA_X_PBE+GGA_C_PBE          Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996)
    GGA     RPBE    GGA_X_RPBE+GGA_C_PBE         Hammer, Hansen, Nørskov, PRB 59, 7413 (1999)
    GGA     revPBE  GGA_X_PBE_R+GGA_C_PBE        Zhang, Yang, PRL 80, 890 (1998)
    GGA     PBEsol  GGA_X_PBE_SOL+GGA_C_PBE_SOL  Perdew et al, PRL 100, 136406 (2008)
    GGA     AM05    GGA_X_AM05+GGA_C_AM05        Armiento, Mattsson, PRB 72, 085108 (2005)
    GGA     BLYP    GGA_X_B88+GGA_C_LYP          Becke, PRA 38, 3098 (1988); Lee, Yang, Parr, PRB 37, 785
    """

    type_name = namedtuple("type_name", "type, name")

    xcf = LibxcFunc
    defined_aliases = OrderedDict(
        [  # (x, c) --> type_name
            # LDAs
            ((xcf.LDA_X, xcf.LDA_C_PW), type_name("LDA", "PW")),  # ixc 7
            ((xcf.LDA_X, xcf.LDA_C_PW_MOD), type_name("LDA", "PW_MOD")),
            ((xcf.LDA_X, xcf.LDA_C_PZ), type_name("LDA", "PZ")),  # ixc 2
            ((xcf.LDA_X, xcf.LDA_C_WIGNER), type_name("LDA", "W")),  # ixc 4
            ((xcf.LDA_X, xcf.LDA_C_HL), type_name("LDA", "HL")),  # ixc 5
            ((xcf.LDA_X, xcf.LDA_C_GL), type_name("LDA", "GL")),
            ((xcf.LDA_X, xcf.LDA_C_VWN), type_name("LDA", "VWN")),
            # GGAs
            ((xcf.GGA_X_PW91, xcf.GGA_C_PW91), type_name("GGA", "PW91")),
            ((xcf.GGA_X_PBE, xcf.GGA_C_PBE), type_name("GGA", "PBE")),
            ((xcf.GGA_X_RPBE, xcf.GGA_C_PBE), type_name("GGA", "RPBE")),  # ixc 15
            ((xcf.GGA_X_PBE_R, xcf.GGA_C_PBE), type_name("GGA", "revPBE")),  # ixc 14
            ((xcf.GGA_X_PBE_SOL, xcf.GGA_C_PBE_SOL), type_name("GGA", "PBEsol")),
            ((xcf.GGA_X_AM05, xcf.GGA_C_AM05), type_name("GGA", "AM05")),
            ((xcf.GGA_X_B88, xcf.GGA_C_LYP), type_name("GGA", "BLYP")),
        ]
    )
    del type_name

    # Correspondence between Abinit ixc notation and libxc notation.
    # see: http://www.abinit.org/doc/helpfiles/for-v7.8/input_variables/varbas.html#ixc
    # and 42_libpaw/m_pawpsp.F90 for the implementation.
    # Fortunately, all the other cases are handled with libxc.
    abinitixc_to_libxc = {
        1: dict(xc=xcf.LDA_XC_TETER93),
        2: dict(x=xcf.LDA_X, c=xcf.LDA_C_PZ),  # PZ  001009
        4: dict(x=xcf.LDA_X, c=xcf.LDA_C_WIGNER),  # W
        5: dict(x=xcf.LDA_X, c=xcf.LDA_C_HL),  # HL
        7: dict(x=xcf.LDA_X, c=xcf.LDA_C_PW),  # PW 001012
        11: dict(x=xcf.GGA_X_PBE, c=xcf.GGA_C_PBE),  # PBE
        14: dict(x=xcf.GGA_X_PBE_R, c=xcf.GGA_C_PBE),  # revPBE
        15: dict(x=xcf.GGA_X_RPBE, c=xcf.GGA_C_PBE),  # RPBE
    }
    del xcf

    @classmethod
    def aliases(cls):
        """List of registered names."""
        return [nt.name for nt in cls.defined_aliases.values()]

    @classmethod
    def asxc(cls, obj):
        """Convert object into Xcfunc."""
        if isinstance(obj, cls):
            return obj
        if is_string(obj):
            return cls.from_name(obj)
        raise TypeError("Don't know how to convert <%s:%s> to Xcfunc" % (type(obj), str(obj)))

    @classmethod
    def from_abinit_ixc(cls, ixc):
        """Build the object from Abinit ixc (integer)"""
        ixc = int(ixc)
        if ixc == 0:
            return None
        if ixc > 0:
            return cls(**cls.abinitixc_to_libxc[ixc])

        # libxc notation employed in Abinit: a six-digit number in the form XXXCCC or CCCXXX
        # ixc = str(ixc)
        # assert len(ixc[1:]) == 6
        # first, last = ixc[1:4], ixc[4:]
        ixc = abs(ixc)
        first = ixc // 1000
        last = ixc - first * 1000
        x, c = LibxcFunc(int(first)), LibxcFunc(int(last))
        if not x.is_x_kind:
            x, c = c, x  # Swap
        assert x.is_x_kind and c.is_c_kind
        return cls(x=x, c=c)

    @classmethod
    def from_name(cls, name):
        """Build the object from one of the registered names"""
        return cls.from_type_name(None, name)

    @classmethod
    def from_type_name(cls, typ, name):
        """Build the object from (type, name)."""
        # Try aliases first.
        for k, nt in cls.defined_aliases.items():
            if typ is not None and typ != nt.type:
                continue
            # print(name, nt.name)
            if name == nt.name:
                if len(k) == 1:
                    return cls(xc=k)
                if len(k) == 2:
                    return cls(x=k[0], c=k[1])
                raise ValueError("Wrong key: %s" % k)

        # At this point, we should have something in the form
        # name="GGA_X_PBE+GGA_C_PBE" or  name=""LDA_XC_TETER93"
        if "+" in name:
            # if typ is not None: raise ValueError("typ: `%s` but name: `%s`" % (typ, name))
            x, c = (s.strip() for s in name.split("+"))
            x, c = LibxcFunc[x], LibxcFunc[c]
            return cls(x=x, c=c)
        # if typ is not None: raise ValueError("typ: `%s` but name: `%s`" % (typ, name))
        xc = LibxcFunc[name]
        return cls(xc=xc)

    @classmethod
    def from_dict(cls, d):
        """
        Makes XcFunc obey the general json interface used in pymatgen for easier serialization.
        """
        return cls(xc=d.get("xc"), x=d.get("x"), c=d.get("c"))

    def as_dict(self):
        """
        Makes XcFunc obey the general json interface used in pymatgen for easier serialization.
        """
        d = {"@module": self.__class__.__module__, "@class": self.__class__.__name__}
        # print("in as_dict", type(self.x), type(self.c), type(self.xc))
        if self.x is not None:
            d["x"] = self.x.as_dict()
        if self.c is not None:
            d["c"] = self.c.as_dict()
        if self.xc is not None:
            d["xc"] = self.xc.as_dict()
        return d

    # def to_json(self):
    #    """
    #    Returns a json string representation of the MSONable object.
    #    """
    #    return json.dumps(self.as_dict()) #, cls=MontyEncoder)

    def __init__(self, xc=None, x=None, c=None):
        """
        Args:
            xc: LibxcFunc for XC functional.
            x, c: LibxcFunc for exchange and correlation part. Mutually exclusive with xc.
        """
        # Consistency check
        if xc is None:
            if x is None or c is None:
                raise ValueError("x or c must be specified when xc is None")
        else:
            if x is not None or c is not None:
                raise ValueError("x and c should be None when xc is specified")

        self.xc, self.x, self.c = xc, x, c

    @lazy_property
    def type(self):
        """The type of the functional."""
        if self.xc in self.defined_aliases:
            return self.defined_aliases[self.xc].type
        xc = (self.x, self.c)
        if xc in self.defined_aliases:
            return self.defined_aliases[xc].type

        # If self is not in defined_aliases, use LibxcFunc family
        if self.xc is not None:
            return self.xc.family
        return "+".join([self.x.family, self.c.family])

    @lazy_property
    def name(self):
        """
        The name of the functional. If the functional is not found in the aliases,
        the string has the form X_NAME+C_NAME
        """
        if self.xc in self.defined_aliases:
            return self.defined_aliases[self.xc].name
        xc = (self.x, self.c)
        if xc in self.defined_aliases:
            return self.defined_aliases[xc].name
        if self.xc is not None:
            return self.xc.name
        return "+".join([self.x.name, self.c.name])

    def __repr__(self):
        return "%s" % self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, XcFunc):
            return self.name == other.name
        # assume other is a string
        return self.name == other

    def __ne__(self, other):
        return not self == other

    # @property
    # def refs(self):

    # def info_dict()
    #    if self.xc is not None:
    #        return {"xc", self.xc.info_dict}
    #    else:
    #        return {"x", self.x.info_dict, "c", self.c.info_dict}
