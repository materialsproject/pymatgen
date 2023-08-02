---
layout: default
title: pymatgen.core.xcfunc.md
nav_exclude: true
---

# pymatgen.core.xcfunc module

This module provides.


### _class_ pymatgen.core.xcfunc.XcFunc(xc=None, x=None, c=None)
Bases: `MSONable`

This object stores information about the XC correlation functional.

Client code usually creates the object by calling the class methods:

>
> * from_name


> * from_type_name

or code-specific methods such as:

>
> * from_abinit_ixc

Ax XcFunc instance is hashable and can therefore be used as key in dictionaries.

The implementation is based on the libxc conventions
and is inspired to the XML specification for atomic PAW datasets documented at:

> [https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html](https://wiki.fysik.dtu.dk/gpaw/setups/pawxml.html)

For convenience, part of the pawxml documentation is reported here.

The xc_functional element defines the exchange-correlation functional used for
generating the dataset. It has the two attributes type and name.

The type attribute can be LDA, GGA, MGGA or HYB.
The name attribute designates the exchange-correlation functional
and can be specified in the following ways:

[1] Taking the names from the LibXC library. The correlation and exchange names

    are stripped from their

    ```
    XC_
    ```

     part and combined with a + sign.
    Here is an example for an LDA functional:

    <xc_functional type=”LDA”, name=”LDA_X+LDA_C_PW”/>

    and this is what PBE will look like:

    <xc_functional type=”GGA”, name=”GGA_X_PBE+GGA_C_PBE”/>

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


* **Parameters**


    * **xc** – LibxcFunc for XC functional.


    * **x** – LibxcFunc for exchange part. Mutually exclusive with xc.


    * **c** – LibxcFunc for correlation part. Mutually exclusive with xc.



#### abinitixc_to_libxc(_ = {1: {'xc': LibxcFunc.LDA_XC_TETER93}, 2: {'c': LibxcFunc.LDA_C_PZ, 'x': LibxcFunc.LDA_X}, 4: {'c': LibxcFunc.LDA_C_WIGNER, 'x': LibxcFunc.LDA_X}, 5: {'c': LibxcFunc.LDA_C_HL, 'x': LibxcFunc.LDA_X}, 7: {'c': LibxcFunc.LDA_C_PW, 'x': LibxcFunc.LDA_X}, 11: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_PBE}, 14: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_PBE_R}, 15: {'c': LibxcFunc.GGA_C_PBE, 'x': LibxcFunc.GGA_X_RPBE}_ )

#### _classmethod_ aliases()
List of registered names.


#### as_dict()
Makes XcFunc obey the general json interface used in pymatgen for easier serialization.


#### _classmethod_ asxc(obj)
Convert object into Xcfunc.


#### defined_aliases(_ = {(<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_PW: 12>): type_name(type='LDA', name='PW'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_PW_MOD: 13>): type_name(type='LDA', name='PW_MOD'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_PZ: 9>): type_name(type='LDA', name='PZ'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_WIGNER: 2>): type_name(type='LDA', name='W'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_HL: 4>): type_name(type='LDA', name='HL'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_GL: 5>): type_name(type='LDA', name='GL'), (<LibxcFunc.LDA_X: 1>, <LibxcFunc.LDA_C_VWN: 7>): type_name(type='LDA', name='VWN'), (<LibxcFunc.GGA_X_PW91: 109>, <LibxcFunc.GGA_C_PW91: 134>): type_name(type='GGA', name='PW91'), (<LibxcFunc.GGA_X_PBE: 101>, <LibxcFunc.GGA_C_PBE: 130>): type_name(type='GGA', name='PBE'), (<LibxcFunc.GGA_X_RPBE: 117>, <LibxcFunc.GGA_C_PBE: 130>): type_name(type='GGA', name='RPBE'), (<LibxcFunc.GGA_X_PBE_R: 102>, <LibxcFunc.GGA_C_PBE: 130>): type_name(type='GGA', name='revPBE'), (<LibxcFunc.GGA_X_PBE_SOL: 116>, <LibxcFunc.GGA_C_PBE_SOL: 133>): type_name(type='GGA', name='PBEsol'), (<LibxcFunc.GGA_X_AM05: 120>, <LibxcFunc.GGA_C_AM05: 135>): type_name(type='GGA', name='AM05'), (<LibxcFunc.GGA_X_B88: 106>, <LibxcFunc.GGA_C_LYP: 131>): type_name(type='GGA', name='BLYP')_ )

#### _classmethod_ from_abinit_ixc(ixc)
Build the object from Abinit ixc (integer).


#### _classmethod_ from_dict(d)
Makes XcFunc obey the general json interface used in pymatgen for easier serialization.


#### _classmethod_ from_name(name)
Build the object from one of the registered names.


#### _classmethod_ from_type_name(typ, name)
Build the object from (type, name).


#### name()
The name of the functional. If the functional is not found in the aliases,
the string has the form X_NAME+C_NAME.


#### type()
The type of the functional.