---
layout: default
title: pymatgen.core.libxcfunc.md
nav_exclude: true
---

# pymatgen.core.libxcfunc module

Enumerator with the libxc identifiers.
This is a low level object, client code should not interact with LibxcFunc directly
but use the API provided by the Xcfunc object defined in core.xcfunc.py.
Part of this module is automatically generated so be careful when refactoring stuff.
Use the script ~pymatgen/dev_scripts/regen_libxcfunc.py to regenerate the enum values.


### _class_ pymatgen.core.libxcfunc.LibxcFunc(value)
Bases: `Enum`

Enumerator with the identifiers. This object is used by Xcfunc
declared in xcfunc.py to create an internal representation of the XC functional.
This is a low level object, client code should not interact with LibxcFunc directly
but use the API provided by Xcfunc.


* **Parameters**

    **num** – Number for the xc.



#### GGA_C_AM05(_ = 13_ )

#### GGA_C_APBE(_ = 18_ )

#### GGA_C_BGCP(_ = 3_ )

#### GGA_C_FT97(_ = 8_ )

#### GGA_C_GAM(_ = 3_ )

#### GGA_C_HCTH_A(_ = 9_ )

#### GGA_C_LM(_ = 13_ )

#### GGA_C_LYP(_ = 13_ )

#### GGA_C_N12(_ = 8_ )

#### GGA_C_N12_SX(_ = 7_ )

#### GGA_C_OPTC(_ = 20_ )

#### GGA_C_OP_B88(_ = 8_ )

#### GGA_C_OP_G96(_ = 8_ )

#### GGA_C_OP_PBE(_ = 8_ )

#### GGA_C_OP_PW91(_ = 26_ )

#### GGA_C_OP_XALPHA(_ = 8_ )

#### GGA_C_P86(_ = 13_ )

#### GGA_C_PBE(_ = 13_ )

#### GGA_C_PBEFE(_ = 25_ )

#### GGA_C_PBEINT(_ = 6_ )

#### GGA_C_PBELOC(_ = 24_ )

#### GGA_C_PBE_JRGX(_ = 13_ )

#### GGA_C_PBE_SOL(_ = 13_ )

#### GGA_C_PW91(_ = 13_ )

#### GGA_C_Q2D(_ = 4_ )

#### GGA_C_REGTPSS(_ = 8_ )

#### GGA_C_REVTCA(_ = 9_ )

#### GGA_C_RGE2(_ = 14_ )

#### GGA_C_SOGGA11(_ = 15_ )

#### GGA_C_SOGGA11_X(_ = 15_ )

#### GGA_C_SPBE(_ = 8_ )

#### GGA_C_TCA(_ = 10_ )

#### GGA_C_WI(_ = 14_ )

#### GGA_C_WI0(_ = 15_ )

#### GGA_C_WL(_ = 14_ )

#### GGA_C_XPBE(_ = 13_ )

#### GGA_C_ZPBEINT(_ = 6_ )

#### GGA_C_ZPBESOL(_ = 6_ )

#### GGA_K_ABSP1(_ = 50_ )

#### GGA_K_ABSP2(_ = 50_ )

#### GGA_K_APBE(_ = 18_ )

#### GGA_K_APBEINT(_ = 5_ )

#### GGA_K_BALTIN(_ = 50_ )

#### GGA_K_DK(_ = 51_ )

#### GGA_K_ERNZERHOF(_ = 52_ )

#### GGA_K_FR_B88(_ = 51_ )

#### GGA_K_FR_PW86(_ = 51_ )

#### GGA_K_GE2(_ = 50_ )

#### GGA_K_GOLDEN(_ = 50_ )

#### GGA_K_GP85(_ = 51_ )

#### GGA_K_GR(_ = 50_ )

#### GGA_K_LC94(_ = 52_ )

#### GGA_K_LIEB(_ = 50_ )

#### GGA_K_LLP(_ = 52_ )

#### GGA_K_LUDENA(_ = 50_ )

#### GGA_K_MEYER(_ = 5_ )

#### GGA_K_OL1(_ = 51_ )

#### GGA_K_OL2(_ = 51_ )

#### GGA_K_PEARSON(_ = 51_ )

#### GGA_K_PERDEW(_ = 51_ )

#### GGA_K_REVAPBE(_ = 5_ )

#### GGA_K_REVAPBEINT(_ = 5_ )

#### GGA_K_TFVW(_ = 5_ )

#### GGA_K_THAKKAR(_ = 52_ )

#### GGA_K_TW1(_ = 18_ )

#### GGA_K_TW2(_ = 18_ )

#### GGA_K_TW3(_ = 18_ )

#### GGA_K_TW4(_ = 19_ )

#### GGA_K_VJKS(_ = 51_ )

#### GGA_K_VSK(_ = 51_ )

#### GGA_K_VW(_ = 50_ )

#### GGA_K_YT65(_ = 50_ )

#### GGA_XC_B97_D(_ = 17_ )

#### GGA_XC_B97_GGA1(_ = 9_ )

#### GGA_XC_EDF1(_ = 16_ )

#### GGA_XC_HCTH_120(_ = 16_ )

#### GGA_XC_HCTH_147(_ = 16_ )

#### GGA_XC_HCTH_407(_ = 16_ )

#### GGA_XC_HCTH_407P(_ = 9_ )

#### GGA_XC_HCTH_93(_ = 16_ )

#### GGA_XC_HCTH_P14(_ = 9_ )

#### GGA_XC_HCTH_P76(_ = 9_ )

#### GGA_XC_KT2(_ = 14_ )

#### GGA_XC_MOHLYP(_ = 19_ )

#### GGA_XC_MOHLYP2(_ = 19_ )

#### GGA_XC_MPWLYP1W(_ = 17_ )

#### GGA_XC_OBLYP_D(_ = 6_ )

#### GGA_XC_OPBE_D(_ = 6_ )

#### GGA_XC_OPWLYP_D(_ = 6_ )

#### GGA_XC_PBE1W(_ = 17_ )

#### GGA_XC_PBELYP1W(_ = 17_ )

#### GGA_XC_TH1(_ = 15_ )

#### GGA_XC_TH2(_ = 15_ )

#### GGA_XC_TH3(_ = 15_ )

#### GGA_XC_TH4(_ = 15_ )

#### GGA_XC_TH_FC(_ = 19_ )

#### GGA_XC_TH_FCFO(_ = 19_ )

#### GGA_XC_TH_FCO(_ = 19_ )

#### GGA_XC_TH_FL(_ = 19_ )

#### GGA_XC_VV10(_ = 25_ )

#### GGA_XC_XLYP(_ = 16_ )

#### GGA_X_2D_B86(_ = 12_ )

#### GGA_X_2D_B86_MGC(_ = 12_ )

#### GGA_X_2D_B88(_ = 12_ )

#### GGA_X_2D_PBE(_ = 12_ )

#### GGA_X_AIRY(_ = 19_ )

#### GGA_X_AK13(_ = 5_ )

#### GGA_X_AM05(_ = 12_ )

#### GGA_X_APBE(_ = 18_ )

#### GGA_X_B86(_ = 10_ )

#### GGA_X_B86_MGC(_ = 10_ )

#### GGA_X_B86_R(_ = 4_ )

#### GGA_X_B88(_ = 10_ )

#### GGA_X_BAYESIAN(_ = 12_ )

#### GGA_X_BGCP(_ = 3_ )

#### GGA_X_BPCCAC(_ = 9_ )

#### GGA_X_C09X(_ = 15_ )

#### GGA_X_CAP(_ = 27_ )

#### GGA_X_DK87_R1(_ = 11_ )

#### GGA_X_DK87_R2(_ = 11_ )

#### GGA_X_EV93(_ = 3_ )

#### GGA_X_FT97_A(_ = 11_ )

#### GGA_X_FT97_B(_ = 11_ )

#### GGA_X_G96(_ = 10_ )

#### GGA_X_GAM(_ = 3_ )

#### GGA_X_HCTH_A(_ = 3_ )

#### GGA_X_HERMAN(_ = 10_ )

#### GGA_X_HJS_B88(_ = 52_ )

#### GGA_X_HJS_B88_V2(_ = 4_ )

#### GGA_X_HJS_B97X(_ = 52_ )

#### GGA_X_HJS_PBE(_ = 52_ )

#### GGA_X_HJS_PBE_SOL(_ = 52_ )

#### GGA_X_HTBS(_ = 19_ )

#### GGA_X_ITYH(_ = 52_ )

#### GGA_X_KT1(_ = 14_ )

#### GGA_X_LAG(_ = 19_ )

#### GGA_X_LAMBDA_CH_N(_ = 4_ )

#### GGA_X_LAMBDA_LO_N(_ = 4_ )

#### GGA_X_LAMBDA_OC2_N(_ = 4_ )

#### GGA_X_LB(_ = 16_ )

#### GGA_X_LBM(_ = 18_ )

#### GGA_X_LG93(_ = 11_ )

#### GGA_X_LV_RPW86(_ = 5_ )

#### GGA_X_MB88(_ = 14_ )

#### GGA_X_MPBE(_ = 12_ )

#### GGA_X_MPW91(_ = 11_ )

#### GGA_X_N12(_ = 8_ )

#### GGA_X_OL2(_ = 18_ )

#### GGA_X_OPTB88_VDW(_ = 13_ )

#### GGA_X_OPTPBE_VDW(_ = 14_ )

#### GGA_X_OPTX(_ = 11_ )

#### GGA_X_PBE(_ = 10_ )

#### GGA_X_PBEA(_ = 12_ )

#### GGA_X_PBEFE(_ = 26_ )

#### GGA_X_PBEINT(_ = 6_ )

#### GGA_X_PBEK1_VDW(_ = 14_ )

#### GGA_X_PBE_JSJR(_ = 12_ )

#### GGA_X_PBE_MOL(_ = 4_ )

#### GGA_X_PBE_R(_ = 10_ )

#### GGA_X_PBE_SOL(_ = 11_ )

#### GGA_X_PBE_TCA(_ = 5_ )

#### GGA_X_PW86(_ = 10_ )

#### GGA_X_PW91(_ = 10_ )

#### GGA_X_Q2D(_ = 4_ )

#### GGA_X_RGE2(_ = 14_ )

#### GGA_X_RPBE(_ = 11_ )

#### GGA_X_RPW86(_ = 14_ )

#### GGA_X_SFAT(_ = 53_ )

#### GGA_X_SOGGA(_ = 15_ )

#### GGA_X_SOGGA11(_ = 15_ )

#### GGA_X_SSB(_ = 9_ )

#### GGA_X_SSB_D(_ = 9_ )

#### GGA_X_SSB_SW(_ = 9_ )

#### GGA_X_VMT84_GE(_ = 6_ )

#### GGA_X_VMT84_PBE(_ = 6_ )

#### GGA_X_VMT_GE(_ = 7_ )

#### GGA_X_VMT_PBE(_ = 7_ )

#### GGA_X_WC(_ = 11_ )

#### GGA_X_WPBEH(_ = 52_ )

#### GGA_X_XPBE(_ = 12_ )

#### HYB_GGA_XC_B1LYP(_ = 41_ )

#### HYB_GGA_XC_B1PW91(_ = 41_ )

#### HYB_GGA_XC_B1WC(_ = 41_ )

#### HYB_GGA_XC_B3LYP(_ = 40_ )

#### HYB_GGA_XC_B3LYP5(_ = 47_ )

#### HYB_GGA_XC_B3LYPs(_ = 45_ )

#### HYB_GGA_XC_B3P86(_ = 40_ )

#### HYB_GGA_XC_B3PW91(_ = 40_ )

#### HYB_GGA_XC_B97(_ = 40_ )

#### HYB_GGA_XC_B97_1(_ = 40_ )

#### HYB_GGA_XC_B97_1p(_ = 26_ )

#### HYB_GGA_XC_B97_2(_ = 41_ )

#### HYB_GGA_XC_B97_3(_ = 41_ )

#### HYB_GGA_XC_B97_K(_ = 41_ )

#### HYB_GGA_XC_BHANDH(_ = 43_ )

#### HYB_GGA_XC_BHANDHLYP(_ = 43_ )

#### HYB_GGA_XC_CAMY_B3LYP(_ = 47_ )

#### HYB_GGA_XC_CAMY_BLYP(_ = 45_ )

#### HYB_GGA_XC_CAM_B3LYP(_ = 43_ )

#### HYB_GGA_XC_CAP0(_ = 47_ )

#### HYB_GGA_XC_EDF2(_ = 47_ )

#### HYB_GGA_XC_HJS_B88(_ = 43_ )

#### HYB_GGA_XC_HJS_B97X(_ = 43_ )

#### HYB_GGA_XC_HJS_PBE(_ = 42_ )

#### HYB_GGA_XC_HJS_PBE_SOL(_ = 43_ )

#### HYB_GGA_XC_HPBEINT(_ = 47_ )

#### HYB_GGA_XC_HSE03(_ = 42_ )

#### HYB_GGA_XC_HSE06(_ = 42_ )

#### HYB_GGA_XC_LCY_BLYP(_ = 46_ )

#### HYB_GGA_XC_LCY_PBE(_ = 46_ )

#### HYB_GGA_XC_LC_VV10(_ = 46_ )

#### HYB_GGA_XC_LRC_WPBE(_ = 47_ )

#### HYB_GGA_XC_LRC_WPBEH(_ = 46_ )

#### HYB_GGA_XC_MB3LYP_RC04(_ = 43_ )

#### HYB_GGA_XC_MPW3LYP(_ = 41_ )

#### HYB_GGA_XC_MPW3PW(_ = 41_ )

#### HYB_GGA_XC_MPWLYP1M(_ = 45_ )

#### HYB_GGA_XC_O3LYP(_ = 40_ )

#### HYB_GGA_XC_PBE0_13(_ = 45_ )

#### HYB_GGA_XC_PBEH(_ = 40_ )

#### HYB_GGA_XC_REVB3LYP(_ = 45_ )

#### HYB_GGA_XC_SB98_1a(_ = 42_ )

#### HYB_GGA_XC_SB98_1b(_ = 42_ )

#### HYB_GGA_XC_SB98_1c(_ = 42_ )

#### HYB_GGA_XC_SB98_2a(_ = 42_ )

#### HYB_GGA_XC_SB98_2b(_ = 42_ )

#### HYB_GGA_XC_SB98_2c(_ = 42_ )

#### HYB_GGA_XC_TUNED_CAM_B3LYP(_ = 43_ )

#### HYB_GGA_XC_WB97(_ = 46_ )

#### HYB_GGA_XC_WB97X(_ = 46_ )

#### HYB_GGA_XC_WB97X_D(_ = 47_ )

#### HYB_GGA_XC_WB97X_V(_ = 46_ )

#### HYB_GGA_XC_X3LYP(_ = 41_ )

#### HYB_GGA_XC_mPW1K(_ = 40_ )

#### HYB_GGA_XC_mPW1PW(_ = 41_ )

#### HYB_GGA_X_N12_SX(_ = 8_ )

#### HYB_GGA_X_SOGGA11_X(_ = 42_ )

#### HYB_MGGA_XC_B86B95(_ = 44_ )

#### HYB_MGGA_XC_B88B95(_ = 44_ )

#### HYB_MGGA_XC_BB1K(_ = 44_ )

#### HYB_MGGA_XC_M05(_ = 43_ )

#### HYB_MGGA_XC_M05_2X(_ = 43_ )

#### HYB_MGGA_XC_M06(_ = 44_ )

#### HYB_MGGA_XC_M06_2X(_ = 45_ )

#### HYB_MGGA_XC_M06_HF(_ = 44_ )

#### HYB_MGGA_XC_M08_HX(_ = 46_ )

#### HYB_MGGA_XC_M08_SO(_ = 46_ )

#### HYB_MGGA_XC_M11(_ = 46_ )

#### HYB_MGGA_XC_MPW1B95(_ = 44_ )

#### HYB_MGGA_XC_MPWB1K(_ = 44_ )

#### HYB_MGGA_XC_PW6B95(_ = 45_ )

#### HYB_MGGA_XC_PW86B95(_ = 44_ )

#### HYB_MGGA_XC_PWB6K(_ = 45_ )

#### HYB_MGGA_XC_REVTPSSH(_ = 45_ )

#### HYB_MGGA_XC_TPSSH(_ = 45_ )

#### HYB_MGGA_XC_WB97M_V(_ = 53_ )

#### HYB_MGGA_XC_X1B95(_ = 44_ )

#### HYB_MGGA_XC_XB1K(_ = 44_ )

#### HYB_MGGA_X_DLDF(_ = 3_ )

#### HYB_MGGA_X_MN12_SX(_ = 24_ )

#### HYB_MGGA_X_MN15(_ = 26_ )

#### HYB_MGGA_X_MS2H(_ = 22_ )

#### HYB_MGGA_X_MVSH(_ = 47_ )

#### HYB_MGGA_X_SCAN0(_ = 26_ )

#### LDA_C_1D_CSC(_ = 1_ )

#### LDA_C_1D_LOOS(_ = 2_ )

#### LDA_C_2D_AMGB(_ = 1_ )

#### LDA_C_2D_PRM(_ = 1_ )

#### LDA_C_GL(_ = _ )

#### LDA_C_GOMBAS(_ = 2_ )

#### LDA_C_HL(_ = _ )

#### LDA_C_ML1(_ = 2_ )

#### LDA_C_ML2(_ = 2_ )

#### LDA_C_OB_PW(_ = 1_ )

#### LDA_C_OB_PZ(_ = 1_ )

#### LDA_C_PW(_ = 1_ )

#### LDA_C_PW_MOD(_ = 1_ )

#### LDA_C_PW_RPA(_ = 2_ )

#### LDA_C_PZ(_ = _ )

#### LDA_C_PZ_MOD(_ = 1_ )

#### LDA_C_RC04(_ = 2_ )

#### LDA_C_RPA(_ = _ )

#### LDA_C_VWN(_ = _ )

#### LDA_C_VWN_1(_ = 2_ )

#### LDA_C_VWN_2(_ = 2_ )

#### LDA_C_VWN_3(_ = 3_ )

#### LDA_C_VWN_4(_ = 3_ )

#### LDA_C_VWN_RPA(_ = _ )

#### LDA_C_WIGNER(_ = _ )

#### LDA_C_XALPHA(_ = _ )

#### LDA_C_vBH(_ = 1_ )

#### LDA_K_LP(_ = 5_ )

#### LDA_K_TF(_ = 5_ )

#### LDA_X(_ = _ )

#### LDA_XC_KSDT(_ = 25_ )

#### LDA_XC_TETER93(_ = 2_ )

#### LDA_XC_ZLP(_ = 4_ )

#### LDA_X_1D(_ = 2_ )

#### LDA_X_2D(_ = 1_ )

#### MGGA_C_BC95(_ = 24_ )

#### MGGA_C_CC06(_ = 22_ )

#### MGGA_C_CS(_ = 7_ )

#### MGGA_C_DLDF(_ = 3_ )

#### MGGA_C_M05(_ = 23_ )

#### MGGA_C_M05_2X(_ = 23_ )

#### MGGA_C_M06(_ = 23_ )

#### MGGA_C_M06_2X(_ = 23_ )

#### MGGA_C_M06_HF(_ = 23_ )

#### MGGA_C_M06_L(_ = 23_ )

#### MGGA_C_M08_HX(_ = 7_ )

#### MGGA_C_M08_SO(_ = 7_ )

#### MGGA_C_M11(_ = 7_ )

#### MGGA_C_M11_L(_ = 7_ )

#### MGGA_C_MN12_L(_ = 7_ )

#### MGGA_C_MN12_SX(_ = 7_ )

#### MGGA_C_MN15(_ = 26_ )

#### MGGA_C_MN15_L(_ = 26_ )

#### MGGA_C_PKZB(_ = 23_ )

#### MGGA_C_REVTPSS(_ = 24_ )

#### MGGA_C_SCAN(_ = 26_ )

#### MGGA_C_TPSS(_ = 23_ )

#### MGGA_C_TPSSLOC(_ = 24_ )

#### MGGA_C_VSXC(_ = 23_ )

#### MGGA_XC_B97M_V(_ = 25_ )

#### MGGA_XC_OTPSS_D(_ = 6_ )

#### MGGA_XC_TPSSLYP1W(_ = 24_ )

#### MGGA_XC_ZLP(_ = 4_ )

#### MGGA_X_2D_PRHG07(_ = 21_ )

#### MGGA_X_2D_PRHG07_PRP10(_ = 21_ )

#### MGGA_X_BJ06(_ = 20_ )

#### MGGA_X_BLOC(_ = 24_ )

#### MGGA_X_BR89(_ = 20_ )

#### MGGA_X_GVT4(_ = 20_ )

#### MGGA_X_LTA(_ = 20_ )

#### MGGA_X_M05(_ = 21_ )

#### MGGA_X_M05_2X(_ = 21_ )

#### MGGA_X_M06(_ = 21_ )

#### MGGA_X_M06_2X(_ = 21_ )

#### MGGA_X_M06_HF(_ = 21_ )

#### MGGA_X_M06_L(_ = 20_ )

#### MGGA_X_M08_HX(_ = 21_ )

#### MGGA_X_M08_SO(_ = 22_ )

#### MGGA_X_M11(_ = 22_ )

#### MGGA_X_M11_L(_ = 22_ )

#### MGGA_X_MBEEF(_ = 24_ )

#### MGGA_X_MBEEFVDW(_ = 25_ )

#### MGGA_X_MK00(_ = 23_ )

#### MGGA_X_MK00B(_ = 24_ )

#### MGGA_X_MN12_L(_ = 22_ )

#### MGGA_X_MN15_L(_ = 26_ )

#### MGGA_X_MODTPSS(_ = 24_ )

#### MGGA_X_MS0(_ = 22_ )

#### MGGA_X_MS1(_ = 22_ )

#### MGGA_X_MS2(_ = 22_ )

#### MGGA_X_MVS(_ = 25_ )

#### MGGA_X_PKZB(_ = 21_ )

#### MGGA_X_REVTPSS(_ = 21_ )

#### MGGA_X_RPP09(_ = 20_ )

#### MGGA_X_SCAN(_ = 26_ )

#### MGGA_X_TAU_HCTH(_ = 20_ )

#### MGGA_X_TB09(_ = 20_ )

#### MGGA_X_TPSS(_ = 20_ )

#### _static_ all_families()
List of strings with the libxc families.
Note that XC_FAMILY if removed from the string e.g. XC_FAMILY_LDA becomes LDA.


#### _static_ all_kinds()
List of strings with the libxc kinds.
Also in this case, the string is obtained by remove the

```
XC_
```

 prefix.
XC_CORRELATION –> CORRELATION.


#### as_dict()
Makes LibxcFunc obey the general json interface used in pymatgen for
easier serialization.


#### _static_ from_dict(d)
Makes LibxcFunc obey the general json interface used in pymatgen for
easier serialization.


#### _property_ info_dict()
Dictionary with metadata. see libxc_docs.json.


#### _property_ is_c_kind(_: boo_ )
True if this is a correlation-only functional.


#### _property_ is_gga_family(_: boo_ )
True if this functional belongs to the GGA family.


#### _property_ is_hyb_gga_family(_: boo_ )
True if this functional belongs to the hybrid + GGA family.


#### _property_ is_hyb_mgga_family(_: boo_ )
True if this functional belongs to the hybrid + meta-GGA family.


#### _property_ is_k_kind(_: boo_ )
True if this is a kinetic functional.


#### _property_ is_lda_family(_: boo_ )
True if this functional belongs to the LDA family.


#### _property_ is_mgga_family(_: boo_ )
True if this functional belongs to the meta-GGA family.


#### _property_ is_x_kind(_: boo_ )
True if this is an exchange-only functional.


#### _property_ is_xc_kind(_: boo_ )
True if this is a exchange+correlation functional.


#### to_json()
Returns a json string representation of the MSONable object.