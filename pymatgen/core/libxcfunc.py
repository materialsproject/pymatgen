# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Enumerator with the libxc identifiers.
This is a low level object, client code should not interact with LibxcFunc directly
but use the API provided by the Xcfunc object defined in core.xcfunc.py.
Part of this module is automatically generated so be careful when refactoring stuff.
Use the script ~pymatgen/dev_scripts/regen_libxcfunc.py to regenerate the enum values.
"""

import json
import os
from enum import Enum
from io import open

from monty.json import MontyEncoder

# The libxc version used to generate this file!
libxc_version = "3.0.0"

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = libxc_version
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo@gmail.com"
__status__ = "Production"
__date__ = "May 16, 2016"

# Loads libxc info from json file
with open(os.path.join(os.path.dirname(__file__), "libxc_docs.json"), "rt") as fh:
    _all_xcfuncs = {int(k): v for k, v in json.load(fh).items()}


# @unique
class LibxcFunc(Enum):
    """
    Enumerator with the identifiers. This object is used by Xcfunc
    declared in xcfunc.py to create an internal representation of the XC functional.
    This is a low level object, client code should not interact with LibxcFunc directly
    but use the API provided by Xcfunc.
    """

    # begin_include_dont_touch
    LDA_C_1D_CSC = 18
    LDA_C_1D_LOOS = 26
    LDA_C_2D_AMGB = 15
    LDA_C_2D_PRM = 16
    LDA_C_GOMBAS = 24
    LDA_C_HL = 4
    LDA_C_GL = 5
    LDA_C_vBH = 17
    LDA_C_ML1 = 22
    LDA_C_ML2 = 23
    LDA_C_PW = 12
    LDA_C_PW_MOD = 13
    LDA_C_OB_PW = 14
    LDA_C_PW_RPA = 25
    LDA_C_PZ = 9
    LDA_C_PZ_MOD = 10
    LDA_C_OB_PZ = 11
    LDA_C_RC04 = 27
    LDA_C_RPA = 3
    LDA_C_VWN = 7
    LDA_C_VWN_1 = 28
    LDA_C_VWN_2 = 29
    LDA_C_VWN_3 = 30
    LDA_C_VWN_4 = 31
    LDA_C_VWN_RPA = 8
    LDA_C_WIGNER = 2
    LDA_K_TF = 50
    LDA_K_LP = 51
    LDA_X = 1
    LDA_C_XALPHA = 6
    LDA_X_1D = 21
    LDA_X_2D = 19
    LDA_XC_KSDT = 259
    LDA_XC_TETER93 = 20
    LDA_XC_ZLP = 43
    GGA_C_AM05 = 135
    GGA_C_FT97 = 88
    GGA_C_LM = 137
    GGA_C_LYP = 131
    GGA_C_OP_B88 = 87
    GGA_C_OP_PBE = 86
    GGA_C_OP_G96 = 85
    GGA_C_OP_PW91 = 262
    GGA_C_OP_XALPHA = 84
    GGA_C_OPTC = 200
    GGA_C_P86 = 132
    GGA_C_PBE = 130
    GGA_C_PBE_SOL = 133
    GGA_C_XPBE = 136
    GGA_C_PBE_JRGX = 138
    GGA_C_RGE2 = 143
    GGA_C_APBE = 186
    GGA_C_SPBE = 89
    GGA_C_REGTPSS = 83
    GGA_C_ZPBESOL = 63
    GGA_C_PBEINT = 62
    GGA_C_ZPBEINT = 61
    GGA_C_PBELOC = 246
    GGA_C_BGCP = 39
    GGA_C_PBEFE = 258
    GGA_C_PW91 = 134
    GGA_C_Q2D = 47
    GGA_C_SOGGA11 = 152
    GGA_C_SOGGA11_X = 159
    GGA_C_TCA = 100
    GGA_C_REVTCA = 99
    GGA_C_WI0 = 153
    GGA_C_WI = 148
    GGA_C_WL = 147
    GGA_K_DK = 516
    GGA_K_PERDEW = 517
    GGA_K_VSK = 518
    GGA_K_VJKS = 519
    GGA_K_ERNZERHOF = 520
    GGA_K_MEYER = 57
    GGA_K_OL1 = 512
    GGA_X_OL2 = 183
    GGA_K_OL2 = 513
    GGA_K_PEARSON = 511
    GGA_K_TFVW = 52
    GGA_K_VW = 500
    GGA_K_GE2 = 501
    GGA_K_GOLDEN = 502
    GGA_K_YT65 = 503
    GGA_K_BALTIN = 504
    GGA_K_LIEB = 505
    GGA_K_ABSP1 = 506
    GGA_K_ABSP2 = 507
    GGA_K_GR = 508
    GGA_K_LUDENA = 509
    GGA_K_GP85 = 510
    GGA_X_2D_B86 = 128
    GGA_X_2D_B86_MGC = 124
    GGA_X_2D_B88 = 127
    GGA_X_2D_PBE = 129
    GGA_X_AIRY = 192
    GGA_X_LAG = 193
    GGA_X_AK13 = 56
    GGA_X_AM05 = 120
    GGA_X_B86 = 103
    GGA_X_B86_MGC = 105
    GGA_X_B86_R = 41
    GGA_X_B88 = 106
    GGA_X_OPTB88_VDW = 139
    GGA_X_MB88 = 149
    GGA_K_LLP = 522
    GGA_K_FR_B88 = 514
    GGA_K_THAKKAR = 523
    GGA_X_BAYESIAN = 125
    GGA_X_BPCCAC = 98
    GGA_X_C09X = 158
    GGA_X_CAP = 270
    GGA_X_DK87_R1 = 111
    GGA_X_DK87_R2 = 112
    GGA_X_EV93 = 35
    GGA_X_FT97_A = 114
    GGA_X_FT97_B = 115
    GGA_X_G96 = 107
    GGA_X_HCTH_A = 34
    GGA_X_HERMAN = 104
    GGA_X_HJS_PBE = 525
    GGA_X_HJS_PBE_SOL = 526
    GGA_X_HJS_B88 = 527
    GGA_X_HJS_B97X = 528
    GGA_X_HJS_B88_V2 = 46
    GGA_X_HTBS = 191
    GGA_X_ITYH = 529
    GGA_X_KT1 = 145
    GGA_XC_KT2 = 146
    GGA_X_LB = 160
    GGA_X_LBM = 182
    GGA_X_LG93 = 113
    GGA_X_LV_RPW86 = 58
    GGA_X_MPBE = 122
    GGA_X_N12 = 82
    GGA_X_GAM = 32
    GGA_X_OPTX = 110
    GGA_X_PBE = 101
    GGA_X_PBE_R = 102
    GGA_X_PBE_SOL = 116
    GGA_X_XPBE = 123
    GGA_X_PBE_JSJR = 126
    GGA_X_PBEK1_VDW = 140
    GGA_X_RGE2 = 142
    GGA_X_APBE = 184
    GGA_X_PBEINT = 60
    GGA_X_PBE_TCA = 59
    GGA_X_LAMBDA_LO_N = 45
    GGA_X_LAMBDA_CH_N = 44
    GGA_X_LAMBDA_OC2_N = 40
    GGA_X_PBE_MOL = 49
    GGA_X_BGCP = 38
    GGA_X_PBEFE = 265
    GGA_K_APBE = 185
    GGA_K_REVAPBE = 55
    GGA_K_TW1 = 187
    GGA_K_TW2 = 188
    GGA_K_TW3 = 189
    GGA_K_TW4 = 190
    GGA_K_APBEINT = 54
    GGA_K_REVAPBEINT = 53
    GGA_X_PBEA = 121
    GGA_X_PW86 = 108
    GGA_X_RPW86 = 144
    GGA_K_FR_PW86 = 515
    GGA_X_PW91 = 109
    GGA_X_MPW91 = 119
    GGA_K_LC94 = 521
    GGA_X_Q2D = 48
    GGA_X_RPBE = 117
    GGA_X_SFAT = 530
    GGA_X_SOGGA11 = 151
    GGA_X_SSB_SW = 90
    GGA_X_SSB = 91
    GGA_X_SSB_D = 92
    GGA_X_VMT_PBE = 71
    GGA_X_VMT_GE = 70
    GGA_X_VMT84_PBE = 69
    GGA_X_VMT84_GE = 68
    GGA_X_WC = 118
    GGA_X_WPBEH = 524
    GGA_XC_XLYP = 166
    GGA_XC_PBE1W = 173
    GGA_XC_MPWLYP1W = 174
    GGA_XC_PBELYP1W = 175
    GGA_XC_B97_D = 170
    GGA_XC_HCTH_93 = 161
    GGA_XC_HCTH_120 = 162
    GGA_XC_HCTH_147 = 163
    GGA_XC_HCTH_407 = 164
    GGA_C_HCTH_A = 97
    GGA_XC_B97_GGA1 = 96
    GGA_XC_HCTH_P14 = 95
    GGA_XC_HCTH_P76 = 94
    GGA_XC_HCTH_407P = 93
    GGA_C_N12 = 80
    GGA_C_N12_SX = 79
    GGA_C_GAM = 33
    GGA_XC_EDF1 = 165
    GGA_X_OPTPBE_VDW = 141
    GGA_XC_MOHLYP = 194
    GGA_XC_MOHLYP2 = 195
    GGA_X_SOGGA = 150
    GGA_XC_OBLYP_D = 67
    GGA_XC_OPWLYP_D = 66
    GGA_XC_OPBE_D = 65
    GGA_XC_TH_FL = 196
    GGA_XC_TH_FC = 197
    GGA_XC_TH_FCFO = 198
    GGA_XC_TH_FCO = 199
    GGA_XC_TH1 = 154
    GGA_XC_TH2 = 155
    GGA_XC_TH3 = 156
    GGA_XC_TH4 = 157
    GGA_XC_VV10 = 255
    HYB_GGA_XC_CAP0 = 477
    HYB_GGA_X_N12_SX = 81
    HYB_GGA_X_SOGGA11_X = 426
    HYB_GGA_XC_B97 = 407
    HYB_GGA_XC_B97_1 = 408
    HYB_GGA_XC_B97_2 = 410
    HYB_GGA_XC_B97_K = 413
    HYB_GGA_XC_B97_3 = 414
    HYB_GGA_XC_SB98_1a = 420
    HYB_GGA_XC_SB98_1b = 421
    HYB_GGA_XC_SB98_1c = 422
    HYB_GGA_XC_SB98_2a = 423
    HYB_GGA_XC_SB98_2b = 424
    HYB_GGA_XC_SB98_2c = 425
    HYB_GGA_XC_WB97 = 463
    HYB_GGA_XC_WB97X = 464
    HYB_GGA_XC_WB97X_V = 466
    HYB_GGA_XC_WB97X_D = 471
    HYB_GGA_XC_B97_1p = 266
    HYB_GGA_XC_LC_VV10 = 469
    HYB_GGA_XC_B1WC = 412
    HYB_GGA_XC_B1LYP = 416
    HYB_GGA_XC_B1PW91 = 417
    HYB_GGA_XC_mPW1PW = 418
    HYB_GGA_XC_mPW1K = 405
    HYB_GGA_XC_BHANDH = 435
    HYB_GGA_XC_BHANDHLYP = 436
    HYB_GGA_XC_MPWLYP1M = 453
    HYB_GGA_XC_B3PW91 = 401
    HYB_GGA_XC_B3LYP = 402
    HYB_GGA_XC_B3LYP5 = 475
    HYB_GGA_XC_B3P86 = 403
    HYB_GGA_XC_MPW3PW = 415
    HYB_GGA_XC_MPW3LYP = 419
    HYB_GGA_XC_MB3LYP_RC04 = 437
    HYB_GGA_XC_REVB3LYP = 454
    HYB_GGA_XC_B3LYPs = 459
    HYB_GGA_XC_CAM_B3LYP = 433
    HYB_GGA_XC_TUNED_CAM_B3LYP = 434
    HYB_GGA_XC_CAMY_B3LYP = 470
    HYB_GGA_XC_CAMY_BLYP = 455
    HYB_GGA_XC_EDF2 = 476
    HYB_GGA_XC_HSE03 = 427
    HYB_GGA_XC_HSE06 = 428
    HYB_GGA_XC_LRC_WPBEH = 465
    HYB_GGA_XC_LRC_WPBE = 473
    HYB_GGA_XC_HJS_PBE = 429
    HYB_GGA_XC_HJS_PBE_SOL = 430
    HYB_GGA_XC_HJS_B88 = 431
    HYB_GGA_XC_HJS_B97X = 432
    HYB_GGA_XC_LCY_BLYP = 468
    HYB_GGA_XC_LCY_PBE = 467
    HYB_GGA_XC_O3LYP = 404
    HYB_GGA_XC_X3LYP = 411
    HYB_GGA_XC_PBEH = 406
    HYB_GGA_XC_PBE0_13 = 456
    HYB_GGA_XC_HPBEINT = 472
    MGGA_XC_TPSSLYP1W = 242
    MGGA_C_BC95 = 240
    MGGA_C_CC06 = 229
    MGGA_C_CS = 72
    MGGA_C_M08_HX = 78
    MGGA_C_M08_SO = 77
    MGGA_C_M11 = 76
    MGGA_C_M11_L = 75
    MGGA_C_MN12_L = 74
    MGGA_C_MN12_SX = 73
    MGGA_C_MN15_L = 261
    MGGA_C_MN15 = 269
    MGGA_C_PKZB = 239
    MGGA_C_TPSS = 231
    MGGA_C_REVTPSS = 241
    MGGA_C_TPSSLOC = 247
    MGGA_C_SCAN = 267
    MGGA_C_M05 = 237
    MGGA_C_M05_2X = 238
    MGGA_C_VSXC = 232
    MGGA_C_M06_L = 233
    MGGA_C_M06_HF = 234
    MGGA_C_M06 = 235
    MGGA_C_M06_2X = 236
    MGGA_C_DLDF = 37
    MGGA_X_2D_PRHG07 = 210
    MGGA_X_2D_PRHG07_PRP10 = 211
    MGGA_X_BR89 = 206
    MGGA_X_BJ06 = 207
    MGGA_X_TB09 = 208
    MGGA_X_RPP09 = 209
    MGGA_X_GVT4 = 204
    MGGA_X_LTA = 201
    MGGA_X_M05 = 214
    MGGA_X_M05_2X = 215
    MGGA_X_M06_2X = 218
    MGGA_X_M06_L = 203
    MGGA_X_M06_HF = 216
    MGGA_X_M06 = 217
    MGGA_X_M08_HX = 219
    MGGA_X_M08_SO = 220
    MGGA_X_M11 = 225
    MGGA_X_M11_L = 226
    MGGA_X_MBEEF = 249
    MGGA_X_MBEEFVDW = 250
    MGGA_X_MK00 = 230
    MGGA_X_MK00B = 243
    MGGA_X_MN12_L = 227
    MGGA_X_MN15_L = 260
    MGGA_X_MS0 = 221
    MGGA_X_MS1 = 222
    MGGA_X_MS2 = 223
    MGGA_X_MVS = 257
    MGGA_X_PKZB = 213
    MGGA_X_SCAN = 263
    MGGA_X_TAU_HCTH = 205
    MGGA_X_TPSS = 202
    MGGA_X_MODTPSS = 245
    MGGA_X_REVTPSS = 212
    MGGA_X_BLOC = 244
    MGGA_XC_B97M_V = 254
    MGGA_XC_OTPSS_D = 64
    MGGA_XC_ZLP = 42
    HYB_MGGA_X_MVSH = 474
    HYB_MGGA_XC_M05 = 438
    HYB_MGGA_XC_M05_2X = 439
    HYB_MGGA_XC_B88B95 = 440
    HYB_MGGA_XC_B86B95 = 441
    HYB_MGGA_XC_PW86B95 = 442
    HYB_MGGA_XC_BB1K = 443
    HYB_MGGA_XC_MPW1B95 = 445
    HYB_MGGA_XC_MPWB1K = 446
    HYB_MGGA_XC_X1B95 = 447
    HYB_MGGA_XC_XB1K = 448
    HYB_MGGA_XC_M06_HF = 444
    HYB_MGGA_XC_M06 = 449
    HYB_MGGA_XC_M06_2X = 450
    HYB_MGGA_XC_PW6B95 = 451
    HYB_MGGA_XC_PWB6K = 452
    HYB_MGGA_XC_TPSSH = 457
    HYB_MGGA_XC_REVTPSSH = 458
    HYB_MGGA_X_DLDF = 36
    HYB_MGGA_XC_M08_HX = 460
    HYB_MGGA_XC_M08_SO = 461
    HYB_MGGA_XC_M11 = 462
    HYB_MGGA_X_MN12_SX = 248
    HYB_MGGA_X_MN15 = 268
    HYB_MGGA_X_MS2H = 224
    HYB_MGGA_X_SCAN0 = 264
    HYB_MGGA_XC_WB97M_V = 531

    # end_include_dont_touch

    def __init__(self, num):
        """
        Init.

        :param num: Number for the xc.
        """
        info = _all_xcfuncs[self.value]
        self.kind = info["Kind"]
        self.family = info["Family"]

    def __str__(self):
        return "name=%s, kind=%s, family=%s" % (self.name, self.kind, self.family)

    @staticmethod
    def all_families():
        """
        List of strings with the libxc families.
        Note that XC_FAMILY if removed from the string e.g. XC_FAMILY_LDA becomes LDA
        """
        return sorted(set(d["Family"] for d in _all_xcfuncs.values()))

    @staticmethod
    def all_kinds():
        """
        List of strings with the libxc kinds.
        Also in this case, the string is obtained by remove the XC_ prefix.
        XC_CORRELATION --> CORRELATION
        """
        return sorted(set(d["Kind"] for d in _all_xcfuncs.values()))

    @property
    def info_dict(self):
        """Dictionary with metadata. see libxc_docs.json"""
        return _all_xcfuncs[self.value]

    @property
    def is_x_kind(self):
        """True if this is an exchange-only functional"""
        return self.kind == "EXCHANGE"

    @property
    def is_c_kind(self):
        """True if this is a correlation-only functional"""
        return self.kind == "CORRELATION"

    @property
    def is_k_kind(self):
        """True if this is a kinetic functional"""
        return self.kind == "KINETIC"

    @property
    def is_xc_kind(self):
        """True if this is a exchange+correlation functional"""
        return self.kind == "EXCHANGE_CORRELATION"

    @property
    def is_lda_family(self):
        """True if this functional belongs to the LDA family."""
        return self.family == "LDA"

    @property
    def is_gga_family(self):
        """True if this functional belongs to the GGA family."""
        return self.family == "GGA"

    @property
    def is_mgga_family(self):
        """True if this functional belongs to the meta-GGA family."""
        return self.family == "MGGA"

    @property
    def is_hyb_gga_family(self):
        """True if this functional belongs to the hybrid + GGA family."""
        return self.family == "HYB_GGA"

    @property
    def is_hyb_mgga_family(self):
        """True if this functional belongs to the hybrid + meta-GGA family."""
        return self.family == "HYB_MGGA"

    def as_dict(self):
        """
        Makes LibxcFunc obey the general json interface used in pymatgen for
        easier serialization.
        """
        return {
            "name": self.name,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }

    @staticmethod
    def from_dict(d):
        """
        Makes LibxcFunc obey the general json interface used in pymatgen for
        easier serialization.
        """
        return LibxcFunc[d["name"]]

    def to_json(self):
        """
        Returns a json string representation of the MSONable object.
        """
        return json.dumps(self.as_dict(), cls=MontyEncoder)


if __name__ == "__main__":
    for xc in LibxcFunc:
        print(xc)
