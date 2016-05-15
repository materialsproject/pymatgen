# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
"""
from __future__ import unicode_literals, division, print_function

class XcFunctional(object):
    """
    https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html

    The xc_functional element defines the exchange-correlation functional used for 
    generating the dataset. It has the two attributes type and name.

    The type attribute can be LDA, GGA, MGGA or HYB.
    The name attribute designates the exchange-correlation functional 
    and can be specified in the following ways:

    [1] Taking the names from the LibXC library. The correlation and exchange names 
        are stripped from their XC_ part and combined with a + sign. Here is an example for an LDA functional:

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

    For the Abinit conventions see:  http://www.abinit.org/doc/helpfiles/for-v7.8/input_variables/varbas.html#ixc
    """

    from collections import namedtuple
    type_name = namedtuple("type_name", "stype, name")

    from pymatgen.io.abinit.libxc import LibxcEnum as xc
    aliases = {
	(xc.LDA_X, xc.LDA_C_PW): type_name("LDA", "PW"),
	(xc.GGA_X_PW91, xc.GGA_C_PW91): type_name("GGA", "PW91"),
	(xc.GGA_X_PBE, xc.GGA_C_PBE): type_name("GGA", "PBE"),
	(xc.GGA_X_RPBE, xc.GGA_C_PBE): type_name("GGA", "RPBE"),
	(xc.GGA_X_PBE_R, xc.GGA_C_PBE): type_name("GGA", "revPBE"), 
	(xc.GGA_X_PBE_SOL, xc.GGA_C_PBE_SOL): type_name("GGA", "PBEsol"),
	(xc.GGA_X_AM05, xc.GGA_C_AM05): type_name("GGA", "AM05"),
	(xc.GGA_X_B88, xc.GGA_C_LYP): type_name("GGA", "BLYP"), 
    }

    @classmethod
    def from_abinit_ixc(cls, ixc_string):
        """Build XC from the value of the Abinit variable ixc (string)"""
        ixc = ixc_string.strip()
        if not ixc.startswith("-"):
            return cls(**ixc2libxc[ixc])
        else:
            # libxc notation employed in Abinit: a six-digit number in the form XXXCCC or CCCXXX
            first, last = ixc[1:4], ixc[4:]
            x, c = LibxcEnum(first), LibxcEnum(last)
	    if not x.is_xonly: # Swap
                x, c = c, x
            assert x.is_xonly and c.isconly
            return cls(x=x, c=c)

    #@classmethod
    #def from_type_name(cls, stype, name):

    def __init__(self, xc=None, x=None, c=None):
        # Consistency check
        if xc is None:
             if x is not None and c is not None:
                raise ValueError("x or c must be specified when xc is None")
        else:
             if x is not None or c is not None:
                raise ValueError("x and c should be None when xc is specified")

        self.xc, self.x, self.c = xc, x, c

    def __repr__(self):
        if self.xc is not None:
            return self.xc.name
        else:
            return "+".join([self.x.name, self.c.name])
    __str__ = __repr__

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self == other

    #@property
    #def xtype(self):
    #    """Exchange family"""

    #@property
    #def xflavor(self):
    #    """Exchange flavor"""

    #@property
    #def ctype(self):
    #    """Correlation family"""
    #                          
    #@property
    #def cflavor(self):
    #    """Correlation flavor"""

    #@property
    #def info(self):

