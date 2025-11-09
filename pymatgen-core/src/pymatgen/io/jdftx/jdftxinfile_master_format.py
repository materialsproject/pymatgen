"""Master list of AbstractTag-type objects for JDFTx input file generation.

This module contains;
- MASTER_TAG_LIST: a dictionary mapping tag categories to dictionaries mapping
    tag names to AbstractTag-type objects.
- get_tag_object: a function that returns an AbstractTag-type object from
    MASTER_TAG_LIST given a tag name.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any

from pymatgen.core.periodic_table import Element
from pymatgen.io.jdftx.generic_tags import (
    AbstractTag,
    BoolTag,
    BoolTagContainer,
    DumpTagContainer,
    FloatTag,
    InitMagMomTag,
    IntTag,
    MultiformatTag,
    StrTag,
    TagContainer,
)
from pymatgen.io.jdftx.jdftxinfile_ref_options import (
    fluid_solvent_options,
    func_c_options,
    func_options,
    func_x_options,
    func_xc_options,
    jdftxdumpfreqoptions,
    jdftxdumpvaroptions,
    jdftxfluid_subtagdict,
    jdftxminimize_subtagdict,
    kinetic_functionals,
)

__author__ = "Jacob Clary, Ben Rich"

MASTER_TAG_LIST: dict[str, dict[str, Any]] = {
    "extrafiles": {
        "include": StrTag(can_repeat=True),
    },
    "structure": {
        "latt-scale": TagContainer(
            allow_list_representation=True,
            subtags={
                "s0": IntTag(write_tagname=False, optional=False),
                "s1": IntTag(write_tagname=False, optional=False),
                "s2": IntTag(write_tagname=False, optional=False),
            },
        ),
        "latt-move-scale": TagContainer(
            allow_list_representation=True,
            subtags={
                "s0": FloatTag(write_tagname=False, optional=False),
                "s1": FloatTag(write_tagname=False, optional=False),
                "s2": FloatTag(write_tagname=False, optional=False),
            },
        ),
        "coords-type": StrTag(options=["Cartesian", "Lattice"]),
        "lattice": MultiformatTag(
            can_repeat=False,
            optional=False,
            format_options=[
                # Explicitly define the lattice vectors
                TagContainer(
                    linebreak_nth_entry=3,
                    allow_list_representation=True,
                    subtags={
                        "R00": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R01": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R02": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R10": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R11": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R12": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R20": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R21": FloatTag(write_tagname=False, optional=False, prec=12),
                        "R22": FloatTag(write_tagname=False, optional=False, prec=12),
                    },
                ),
                TagContainer(
                    subtags={
                        "Triclinic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "b": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                        "alpha": FloatTag(write_tagname=False, optional=False),
                        "beta": FloatTag(write_tagname=False, optional=False),
                        "gamma": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "modification": StrTag(
                            options=["Base-Centered"],  # Single-option modification could be a boolean tag, but keeping
                            # as a string tagfor consistency with the other modified lattice tags
                            optional=False,
                            write_tagname=False,
                            write_value=True,
                        ),
                        "Monoclinic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "b": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                        "beta": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                # This is a duplicate of the above tag, but without the modification subtag.
                # This is not an elegant solution, but this approach only adds < 50 lines of code. The long term fix
                # will be to make major changes to the TagContainer `read` implementation in the way it expects
                # values to be ordered/present after parsing subtags with write_tagname=True.
                # Another option would be to make all modification options boolean tags, but this would allow mutually
                # exclusive tags to be present in the same tag, which is not ideal.
                TagContainer(
                    subtags={
                        "Monoclinic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "b": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                        "beta": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "modification": StrTag(
                            options=["Base-Centered", "Body-Centered", "Face-Centered"],
                            optional=False,
                            write_tagname=False,
                            write_value=True,
                        ),
                        "Orthorhombic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "b": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "Orthorhombic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "b": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "modification": StrTag(
                            options=["Body-Centered"],  # This could be a boolean tag, but keeping a string tag
                            # for consistency with the other modified lattice tags
                            optional=False,
                            write_tagname=False,
                            write_value=True,
                        ),
                        "Tetragonal": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "Tetragonal": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "Rhombohedral": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "alpha": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "Hexagonal": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                        "c": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "modification": StrTag(
                            options=["Body-Centered", "Face-Centered"],
                            optional=False,
                            write_tagname=False,
                            write_value=True,
                        ),
                        "Cubic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "Cubic": BoolTag(write_tagname=True, optional=False, write_value=False),
                        "a": FloatTag(write_tagname=False, optional=False),
                    }
                ),
            ],
        ),
        "ion": MultiformatTag(
            can_repeat=True,
            optional=False,
            allow_list_representation=False,
            format_options=[
                TagContainer(
                    allow_list_representation=False,
                    can_repeat=True,
                    subtags={
                        "species-id": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=[
                                value.symbol
                                for key, value in Element.__dict__.items()
                                if not key.startswith("_") and not callable(value)
                            ],  # Required in case bad species names gets fed
                        ),
                        "x0": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x1": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x2": FloatTag(write_tagname=False, optional=False, prec=12),
                        "v": TagContainer(
                            allow_list_representation=True,
                            subtags={
                                "vx0": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx1": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx2": FloatTag(write_tagname=False, optional=False, prec=12),
                            },
                        ),
                        "moveScale": IntTag(write_tagname=False, optional=False),
                    },
                ),
                TagContainer(
                    allow_list_representation=False,
                    can_repeat=True,
                    subtags={
                        "species-id": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=[
                                value.symbol
                                for key, value in Element.__dict__.items()
                                if not key.startswith("_") and not callable(value)
                            ],  # Required in case bad species names gets fed
                        ),
                        "x0": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x1": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x2": FloatTag(write_tagname=False, optional=False, prec=12),
                        "v": TagContainer(
                            allow_list_representation=True,
                            subtags={
                                "vx0": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx1": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx2": FloatTag(write_tagname=False, optional=False, prec=12),
                            },
                        ),
                        "moveScale": IntTag(write_tagname=False, optional=False),
                        "constraint type": StrTag(
                            options=["Linear", "Planar", "HyperPlane"],
                            write_tagname=False,
                            optional=False,
                        ),
                        "d0": FloatTag(write_tagname=False, optional=False, prec=12),
                        "d1": FloatTag(write_tagname=False, optional=False, prec=12),
                        "d2": FloatTag(write_tagname=False, optional=False, prec=12),
                    },
                ),
                TagContainer(
                    allow_list_representation=False,
                    can_repeat=True,
                    subtags={
                        "species-id": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=[
                                value.symbol
                                for key, value in Element.__dict__.items()
                                if not key.startswith("_") and not callable(value)
                            ],  # Required in case bad species names gets fed
                        ),
                        "x0": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x1": FloatTag(write_tagname=False, optional=False, prec=12),
                        "x2": FloatTag(write_tagname=False, optional=False, prec=12),
                        "v": TagContainer(
                            allow_list_representation=True,
                            subtags={
                                "vx0": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx1": FloatTag(write_tagname=False, optional=False, prec=12),
                                "vx2": FloatTag(write_tagname=False, optional=False, prec=12),
                            },
                        ),
                        "moveScale": IntTag(write_tagname=False, optional=False),
                        "HyperPlane": TagContainer(
                            write_tagname=True,
                            can_repeat=True,
                            subtags={
                                "d0": FloatTag(write_tagname=False, optional=False, prec=12),
                                "d1": FloatTag(write_tagname=False, optional=False, prec=12),
                                "d2": FloatTag(write_tagname=False, optional=False, prec=12),
                                "group": StrTag(write_tagname=False, optional=False),
                            },
                        ),
                    },
                ),
            ],
        ),
        "perturb-ion": TagContainer(
            subtags={
                "species": StrTag(write_tagname=False, optional=False),
                "atom": IntTag(write_tagname=False, optional=False),
                "dx0": FloatTag(write_tagname=False, optional=False),
                "dx1": FloatTag(write_tagname=False, optional=False),
                "dx2": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "core-overlap-check": StrTag(options=["additive", "vector", "none"]),
        "ion-species": StrTag(can_repeat=True, optional=False),
        "cache-projectors": BoolTag(),
        "ion-width": MultiformatTag(
            format_options=[
                StrTag(options=["Ecut", "fftbox"]),
                FloatTag(),
            ]
        ),
    },
    "symmetries": {
        "symmetries": StrTag(options=["automatic", "manual", "none"]),
        "symmetry-threshold": FloatTag(),
        "symmetry-matrix": TagContainer(
            linebreak_nth_entry=3,
            can_repeat=True,
            allow_list_representation=True,
            subtags={
                "s00": IntTag(write_tagname=False, optional=False),
                "s01": IntTag(write_tagname=False, optional=False),
                "s02": IntTag(write_tagname=False, optional=False),
                "s10": IntTag(write_tagname=False, optional=False),
                "s11": IntTag(write_tagname=False, optional=False),
                "s12": IntTag(write_tagname=False, optional=False),
                "s20": IntTag(write_tagname=False, optional=False),
                "s21": IntTag(write_tagname=False, optional=False),
                "s22": IntTag(write_tagname=False, optional=False),
                "a0": FloatTag(write_tagname=False, optional=False, prec=12),
                "a1": FloatTag(write_tagname=False, optional=False, prec=12),
                "a2": FloatTag(write_tagname=False, optional=False, prec=12),
            },
        ),
    },
    "k-mesh": {
        "kpoint": TagContainer(
            can_repeat=True,
            allow_list_representation=True,
            subtags={
                "k0": FloatTag(write_tagname=False, optional=False, prec=12),
                "k1": FloatTag(write_tagname=False, optional=False, prec=12),
                "k2": FloatTag(write_tagname=False, optional=False, prec=12),
                "weight": FloatTag(write_tagname=False, optional=False, prec=12),
            },
        ),
        "kpoint-folding": TagContainer(
            allow_list_representation=True,
            subtags={
                "n0": IntTag(write_tagname=False, optional=False),
                "n1": IntTag(write_tagname=False, optional=False),
                "n2": IntTag(write_tagname=False, optional=False),
            },
        ),
        "kpoint-reduce-inversion": BoolTag(),
    },
    "electronic": {
        "elec-ex-corr": MultiformatTag(
            format_options=[
                # note that hyb-HSE06 has a bug in JDFTx and should not be
                # used and is excluded here use the LibXC version instead
                # (hyb-gga-HSE06)
                StrTag(
                    write_tagname=True,
                    options=deepcopy(func_options),
                ),
                TagContainer(
                    subtags={
                        "funcX": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_x_options),
                        ),
                        "funcC": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_c_options),
                        ),
                    }
                ),
                TagContainer(
                    subtags={
                        "funcXC": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_xc_options),
                        )
                    }
                ),
            ]
        ),
        "elec-ex-corr-compare": MultiformatTag(
            can_repeat=True,
            format_options=[
                # note that hyb-HSE06 has a bug in JDFTx and should not be used
                # and is excluded here use the LibXC version instead
                # (hyb-gga-HSE06)
                StrTag(
                    write_tagname=True,
                    options=func_options,
                ),
                TagContainer(
                    subtags={
                        "funcX": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_x_options),
                        ),
                        "funcC": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_c_options),
                        ),
                    }
                ),
                TagContainer(
                    subtags={
                        "funcXC": StrTag(
                            write_tagname=False,
                            optional=False,
                            options=deepcopy(func_xc_options),
                        )
                    }
                ),
            ],
        ),
        "exchange-block-size": IntTag(),
        "exchange-outer-loop": IntTag(),
        "exchange-parameters": TagContainer(
            subtags={
                "exxScale": FloatTag(write_tagname=False, optional=False),
                "exxOmega": FloatTag(write_tagname=False),
            }
        ),
        "exchange-params": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "blockSize": IntTag(),
                "nOuterVxx": IntTag(),
            },
        ),
        "exchange-regularization": StrTag(
            options=[
                "AuxiliaryFunction",
                "None",
                "ProbeChargeEwald",
                "SphericalTruncated",
                "WignerSeitzTruncated",
            ]
        ),
        "tau-core": TagContainer(
            subtags={
                "species-id": StrTag(write_tagname=False, optional=False),
                "rCut": FloatTag(write_tagname=False),
                "plot": BoolTag(write_tagname=False),
            }
        ),
        "lj-override": FloatTag(),
        "van-der-waals": MultiformatTag(
            format_options=[
                StrTag(options=["D3"]),
                FloatTag(),
            ]
        ),
        "elec-cutoff": TagContainer(
            allow_list_representation=True,
            subtags={
                "Ecut": FloatTag(write_tagname=False, optional=False),
                "EcutRho": FloatTag(write_tagname=False),
            },
        ),
        "elec-smearing": TagContainer(
            allow_list_representation=True,
            subtags={
                "smearingType": StrTag(
                    options=["Cold", "Fermi", "Gauss", "MP1"],
                    write_tagname=False,
                    optional=False,
                ),
                "smearingWidth": FloatTag(write_tagname=False, optional=False),
            },
        ),
        "elec-n-bands": IntTag(),
        "spintype": StrTag(options=["no-spin", "spin-orbit", "vector-spin", "z-spin"]),
        "initial-magnetic-moments": InitMagMomTag(),
        "elec-initial-magnetization": TagContainer(
            subtags={
                "M": FloatTag(write_tagname=False, optional=False),
                "constrain": BoolTag(write_tagname=False, optional=False),
            }
        ),
        "target-Bz": FloatTag(),
        "elec-initial-charge": FloatTag(),
        "converge-empty-states": BoolTag(),
        "band-unfold": TagContainer(
            linebreak_nth_entry=3,
            allow_list_representation=True,
            subtags={
                "M00": IntTag(write_tagname=False, optional=False),
                "M01": IntTag(write_tagname=False, optional=False),
                "M02": IntTag(write_tagname=False, optional=False),
                "M10": IntTag(write_tagname=False, optional=False),
                "M11": IntTag(write_tagname=False, optional=False),
                "M12": IntTag(write_tagname=False, optional=False),
                "M20": IntTag(write_tagname=False, optional=False),
                "M21": IntTag(write_tagname=False, optional=False),
                "M22": IntTag(write_tagname=False, optional=False),
            },
        ),
        "basis": StrTag(options=["kpoint-dependent", "single"]),
        "fftbox": TagContainer(
            allow_list_representation=True,
            subtags={
                "S0": IntTag(write_tagname=False, optional=False),
                "S1": IntTag(write_tagname=False, optional=False),
                "S2": IntTag(write_tagname=False, optional=False),
            },
        ),
        "electric-field": TagContainer(
            allow_list_representation=True,
            subtags={
                "Ex": IntTag(write_tagname=False, optional=False),
                "Ey": IntTag(write_tagname=False, optional=False),
                "Ez": IntTag(write_tagname=False, optional=False),
            },
        ),
        "perturb-electric-field": TagContainer(
            allow_list_representation=True,
            subtags={
                "Ex": IntTag(write_tagname=False, optional=False),
                "Ey": IntTag(write_tagname=False, optional=False),
                "Ez": IntTag(write_tagname=False, optional=False),
            },
        ),
        "box-potential": TagContainer(
            can_repeat=True,
            subtags={
                "xmin": FloatTag(write_tagname=False, optional=False),
                "xmax": FloatTag(write_tagname=False, optional=False),
                "ymin": FloatTag(write_tagname=False, optional=False),
                "ymax": FloatTag(write_tagname=False, optional=False),
                "zmin": FloatTag(write_tagname=False, optional=False),
                "zmax": FloatTag(write_tagname=False, optional=False),
                "Vin": FloatTag(write_tagname=False, optional=False),
                "Vout": FloatTag(write_tagname=False, optional=False),
                "convolve_radius": FloatTag(write_tagname=False),
            },
        ),
        "ionic-gaussian-potential": TagContainer(
            can_repeat=True,
            subtags={
                "species": StrTag(write_tagname=False, optional=False),
                "U0": FloatTag(write_tagname=False, optional=False),
                "sigma": FloatTag(write_tagname=False, optional=False),
                "geometry": StrTag(
                    options=["Spherical", "Cylindrical", "Planar"],
                    write_tagname=False,
                    optional=False,
                ),
            },
        ),
        "bulk-epsilon": TagContainer(
            subtags={
                "DtotFile": StrTag(write_tagname=False, optional=False),
                "Ex": FloatTag(write_tagname=False),
                "Ey": FloatTag(write_tagname=False),
                "Ez": FloatTag(write_tagname=False),
            }
        ),
        "charged-defect": TagContainer(
            can_repeat=True,
            subtags={
                "x0": FloatTag(write_tagname=False, optional=False),
                "x1": FloatTag(write_tagname=False, optional=False),
                "x2": FloatTag(write_tagname=False, optional=False),
                "q": FloatTag(write_tagname=False, optional=False),
                "sigma": FloatTag(write_tagname=False, optional=False),
            },
        ),
        "charged-defect-correction": TagContainer(
            subtags={
                "Slab": TagContainer(
                    subtags={
                        "dir": StrTag(options=["100", "010", "001"], write_tagname=False),
                    }
                ),
                "DtotFile": StrTag(write_tagname=False, optional=False),
                "Eps": MultiformatTag(
                    format_options=[
                        FloatTag(write_tagname=False, optional=False),
                        StrTag(write_tagname=False, optional=False),
                    ]
                ),
                "rMin": FloatTag(write_tagname=False, optional=False),
                "rSigma": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "Cprime-params": TagContainer(
            subtags={
                "dk": FloatTag(write_tagname=False),
                "degeneracyThreshold": FloatTag(write_tagname=False),
                "vThreshold": FloatTag(write_tagname=False),
                "realSpaceTruncated": BoolTag(write_tagname=False),
            }
        ),
        "electron-scattering": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "eta": FloatTag(optional=False),
                "Ecut": FloatTag(),
                "fCut": FloatTag(),
                "omegaMax": FloatTag(),
                "RPA": BoolTag(),
                "dumpEpsilon": BoolTag(),
                "slabResponse": BoolTag(),
                "EcutTransverse": FloatTag(),
                "computeRange": TagContainer(
                    subtags={
                        "iqStart": FloatTag(write_tagname=False, optional=False),
                        "iqStop": FloatTag(write_tagname=False, optional=False),
                    }
                ),
            },
        ),
        "perturb-test": BoolTag(write_value=False),
        "perturb-wavevector": TagContainer(
            subtags={
                "q0": FloatTag(write_tagname=False, optional=False),
                "q1": FloatTag(write_tagname=False, optional=False),
                "q2": FloatTag(write_tagname=False, optional=False),
            }
        ),
    },
    "truncation": {
        "coulomb-interaction": MultiformatTag(
            format_options=[
                # note that the first 2 and last 2 TagContainers could be
                # combined, but keep separate so there is less ambiguity on
                # formatting
                TagContainer(
                    subtags={
                        "truncationType": StrTag(
                            options=["Periodic", "Isolated"],
                            write_tagname=False,
                            optional=False,
                        )
                    }
                ),
                TagContainer(
                    subtags={
                        "truncationType": StrTag(options=["Spherical"], write_tagname=False, optional=False),
                        "Rc": FloatTag(write_tagname=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "truncationType": StrTag(
                            options=["Slab", "Wire"],
                            write_tagname=False,
                            optional=False,
                        ),
                        "dir": StrTag(
                            options=["001", "010", "100"],
                            write_tagname=False,
                            optional=False,
                        ),
                    }
                ),
                TagContainer(
                    subtags={
                        "truncationType": StrTag(options=["Cylindrical"], write_tagname=False, optional=False),
                        "dir": StrTag(
                            options=["001", "010", "100"],
                            write_tagname=False,
                            optional=False,
                        ),
                        "Rc": FloatTag(write_tagname=False),
                    }
                ),
            ]
        ),
        "coulomb-truncation-embed": TagContainer(
            subtags={
                "c0": FloatTag(write_tagname=False, optional=False),
                "c1": FloatTag(write_tagname=False, optional=False),
                "c2": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "coulomb-truncation-ion-margin": FloatTag(),
    },
    "restart": {
        "initial-state": StrTag(),
        "elec-initial-eigenvals": StrTag(),
        "elec-initial-fillings": TagContainer(
            subtags={
                "read": BoolTag(write_value=False, optional=False),
                "filename": StrTag(write_tagname=False, optional=False),
                "nBandsOld": IntTag(write_tagname=False),
            }
        ),
        "wavefunction": MultiformatTag(
            format_options=[
                TagContainer(subtags={"lcao": BoolTag(write_value=False, optional=False)}),
                TagContainer(subtags={"random": BoolTag(write_value=False, optional=False)}),
                TagContainer(
                    subtags={
                        "read": BoolTag(write_value=False, optional=False),
                        "filename": StrTag(write_tagname=False, optional=False),
                        "nBandsOld": IntTag(write_tagname=False),
                        "EcutOld": FloatTag(write_tagname=False),
                    }
                ),
                TagContainer(
                    subtags={
                        "read-rs": BoolTag(write_value=False, optional=False),
                        "filename-pattern": StrTag(write_tagname=False, optional=False),
                        "nBandsOld": IntTag(write_tagname=False),
                        "NxOld": IntTag(write_tagname=False),
                        "NyOld": IntTag(write_tagname=False),
                        "NzOld": IntTag(write_tagname=False),
                    }
                ),
            ]
        ),
        "fluid-initial-state": StrTag(),
        "perturb-incommensurate-wavefunctions": TagContainer(
            subtags={
                "filename": StrTag(write_tagname=False, optional=False),
                "EcutOld": IntTag(write_tagname=False),
            }
        ),
        "perturb-rhoExternal": StrTag(),
        "perturb-Vexternal": StrTag(),
        "fix-electron-density": StrTag(),
        "fix-electron-potential": StrTag(),
        "Vexternal": MultiformatTag(
            format_options=[
                TagContainer(subtags={"filename": StrTag(write_value=False, optional=False)}),
                TagContainer(
                    subtags={
                        "filenameUp": StrTag(write_value=False, optional=False),
                        "filenameDn": StrTag(write_tagname=False, optional=False),
                    }
                ),
            ]
        ),
        "rhoExternal": TagContainer(
            subtags={
                "filename": StrTag(write_tagname=False, optional=False),
                "includeSelfEnergy": FloatTag(write_tagname=False),
            }
        ),
        "slab-epsilon": TagContainer(
            subtags={
                "DtotFile": StrTag(write_tagname=False, optional=False),
                "sigma": FloatTag(write_tagname=False, optional=False),
                "Ex": FloatTag(write_tagname=False),
                "Ey": FloatTag(write_tagname=False),
                "Ez": FloatTag(write_tagname=False),
            }
        ),
    },
    "minimization": {
        "lcao-params": TagContainer(
            subtags={
                "nIter": IntTag(write_tagname=False),
                "Ediff": FloatTag(write_tagname=False),
                "smearingWidth": FloatTag(write_tagname=False),
            }
        ),
        "elec-eigen-algo": StrTag(options=["CG", "Davidson"]),
        "ionic-minimize": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                **deepcopy(jdftxminimize_subtagdict),
            },
        ),
        "lattice-minimize": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                **deepcopy(jdftxminimize_subtagdict),
            },
        ),
        "electronic-minimize": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                **deepcopy(jdftxminimize_subtagdict),
            },
        ),
        "electronic-scf": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "energyDiffThreshold": FloatTag(),
                "history": IntTag(),
                "mixFraction": FloatTag(),
                "nIterations": IntTag(),
                "qMetric": FloatTag(),
                "residualThreshold": FloatTag(),
                "eigDiffThreshold": FloatTag(),
                "mixedVariable": StrTag(),
                "mixFractionMag": FloatTag(),
                "nEigSteps": IntTag(),
                "qKappa": FloatTag(),
                "qKerker": FloatTag(),
                "verbose": BoolTag(),
            },
        ),
        "fluid-minimize": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                **deepcopy(jdftxminimize_subtagdict),
            },
        ),
        "davidson-band-ratio": FloatTag(),
        "wavefunction-drag": BoolTag(),
        "subspace-rotation-factor": TagContainer(
            subtags={
                "factor": FloatTag(write_tagname=False, optional=False),
                "adjust": BoolTag(write_tagname=False, optional=False),
            }
        ),
        "perturb-minimize": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "algorithm": StrTag(options=["MINRES", "CGIMINRES"]),
                "CGBypass": BoolTag(),
                "nIterations": IntTag(),
                "recomputeResidual": BoolTag(),
                "residualDiffThreshold": FloatTag(),
                "residualTol": FloatTag(),
            },
        ),
    },
    "fluid": {
        "target-mu": TagContainer(
            allow_list_representation=True,
            subtags={
                "mu": FloatTag(write_tagname=False, optional=False),
                "outerLoop": BoolTag(write_tagname=False),
            },
        ),
        "fluid": TagContainer(
            subtags={
                "type": StrTag(
                    options=[
                        "None",
                        "LinearPCM",
                        "NonlinearPCM",
                        "SaLSA",
                        "ClassicalDFT",
                    ],
                    write_tagname=False,
                    optional=False,
                ),
                "Temperature": FloatTag(write_tagname=False),
                "Pressure": FloatTag(write_tagname=False),
            }
        ),
        "fluid-solvent": MultiformatTag(
            can_repeat=True,  # 11/27
            format_options=[
                TagContainer(
                    can_repeat=True,
                    subtags={
                        "name": StrTag(
                            options=fluid_solvent_options,
                            write_tagname=False,
                        ),
                        "concentration": FloatTag(write_tagname=False),
                        "functional": StrTag(
                            options=[
                                "BondedVoids",
                                "FittedCorrelations",
                                "MeanFieldLJ",
                                "ScalarEOS",
                            ],
                            write_tagname=False,
                        ),
                        **deepcopy(jdftxfluid_subtagdict),
                    },
                ),
                TagContainer(
                    can_repeat=True,
                    subtags={
                        "name": StrTag(
                            options=fluid_solvent_options,
                            write_tagname=False,
                        ),
                        "concentration": StrTag(options=["bulk"], write_tagname=False),
                        "functional": StrTag(
                            options=[
                                "BondedVoids",
                                "FittedCorrelations",
                                "MeanFieldLJ",
                                "ScalarEOS",
                            ],
                            write_tagname=False,
                        ),
                        **deepcopy(jdftxfluid_subtagdict),
                    },
                ),
            ],
        ),
        "fluid-anion": TagContainer(
            subtags={
                "name": StrTag(options=["Cl-", "ClO4-", "F-"], write_tagname=False, optional=False),
                "concentration": FloatTag(write_tagname=False, optional=False),
                "functional": StrTag(
                    options=[
                        "BondedVoids",
                        "FittedCorrelations",
                        "MeanFieldLJ",
                        "ScalarEOS",
                    ],
                    write_tagname=False,
                ),
                **deepcopy(jdftxfluid_subtagdict),
            }
        ),
        "fluid-cation": TagContainer(
            subtags={
                "name": StrTag(options=["K+", "Na+"], write_tagname=False, optional=False),
                "concentration": FloatTag(write_tagname=False, optional=False),
                "functional": StrTag(
                    options=[
                        "BondedVoids",
                        "FittedCorrelations",
                        "MeanFieldLJ",
                        "ScalarEOS",
                    ],
                    write_tagname=False,
                ),
                **deepcopy(jdftxfluid_subtagdict),
            }
        ),
        "fluid-dielectric-constant": TagContainer(
            subtags={
                "epsBulkOverride": FloatTag(write_tagname=False),
                "epsInfOverride": FloatTag(write_tagname=False),
            }
        ),
        "fluid-dielectric-tensor": TagContainer(
            subtags={
                "epsBulkXX": FloatTag(write_tagname=False, optional=False),
                "epsBulkYY": FloatTag(write_tagname=False, optional=False),
                "epsBulkZZ": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "fluid-ex-corr": TagContainer(
            subtags={
                "kinetic": StrTag(write_tagname=False, optional=False, options=deepcopy(kinetic_functionals)),
                "exchange-correlation": StrTag(write_tagname=False, options=deepcopy(func_options)),
            }
        ),
        "fluid-mixing-functional": TagContainer(
            can_repeat=True,
            subtags={
                "fluid1": StrTag(
                    options=[
                        "CCl4",
                        "CH3CN",
                        "CHCl3",
                        "Cl-",
                        "ClO4-",
                        "CustomAnion",
                        "CustomCation",
                        "F-",
                        "H2O",
                        "Na(H2O)4+",
                        "Na+",
                    ],
                    write_tagname=False,
                    optional=False,
                ),
                "fluid2": StrTag(
                    options=[
                        "CCl4",
                        "CH3CN",
                        "CHCl3",
                        "Cl-",
                        "ClO4-",
                        "CustomAnion",
                        "CustomCation",
                        "F-",
                        "H2O",
                        "Na(H2O)4+",
                        "Na+",
                    ],
                    write_tagname=False,
                    optional=False,
                ),
                "energyScale": FloatTag(write_tagname=False, optional=False),
                "lengthScale": FloatTag(write_tagname=False),
                "FMixType": StrTag(options=["LJPotential", "GaussianKernel"], write_tagname=False),
            },
        ),
        "fluid-vdwScale": FloatTag(),
        "fluid-gummel-loop": TagContainer(
            subtags={
                "maxIterations": IntTag(write_tagname=False, optional=False),
                "Atol": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "fluid-solve-frequency": StrTag(options=["Default", "Gummel", "Inner"]),
        "fluid-site-params": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            can_repeat=True,
            subtags={
                "component": StrTag(
                    options=[
                        "CCl4",
                        "CH3CN",
                        "CHCl3",
                        "Cl-",
                        "ClO4-",
                        "CustomAnion",
                        "CustomCation",
                        "F-",
                        "H2O",
                        "Na(H2O)4+",
                        "Na+",
                    ],
                    optional=False,
                ),
                "siteName": StrTag(optional=False),
                "aElec": FloatTag(),
                "alpha": FloatTag(),
                "aPol": FloatTag(),
                "elecFilename": StrTag(),
                "elecFilenameG": StrTag(),
                "rcElec": FloatTag(),
                "Rhs": FloatTag(),
                "sigmaElec": FloatTag(),
                "sigmaNuc": FloatTag(),
                "Zelec": FloatTag(),
                "Znuc": FloatTag(),
            },
        ),
        "pcm-variant": StrTag(
            options=[
                "CANDLE",
                "CANON",
                "FixedCavity",
                "GLSSA13",
                "LA12",
                "SCCS_anion",
                "SCCS_cation",
                "SCCS_g03",
                "SCCS_g03beta",
                "SCCS_g03p",
                "SCCS_g03pbeta",
                "SCCS_g09",
                "SCCS_g09beta",
                "SGA13",
                "SoftSphere",
            ]
        ),
        "pcm-nonlinear-scf": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "energyDiffThreshold": FloatTag(),
                "history": IntTag(),
                "mixFraction": FloatTag(),
                "nIterations": IntTag(),
                "qMetric": FloatTag(),
                "residualThreshold": FloatTag(),
            },
        ),
        "pcm-params": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "cavityFile": StrTag(),
                "cavityPressure": FloatTag(),
                "cavityScale": FloatTag(),
                "cavityTension": FloatTag(),
                "eta_wDiel": FloatTag(),
                "ionSpacing": FloatTag(),
                "lMax": FloatTag(),
                "nc": FloatTag(),
                "pCavity": FloatTag(),
                "rhoDelta": FloatTag(),
                "rhoMax": FloatTag(),
                "rhoMin": FloatTag(),
                "screenOverride": FloatTag(),
                "sigma": FloatTag(),
                "sqrtC6eff": FloatTag(),
                "Zcenter": FloatTag(),
                "zMask0": FloatTag(),
                "zMaskH": FloatTag(),
                "zMaskIonH": FloatTag(),
                "zMaskSigma": FloatTag(),
                "Ztot": FloatTag(),
            },
        ),
    },
    "dynamics": {
        "vibrations": TagContainer(
            subtags={
                "dr": FloatTag(),
                "centralDiff": BoolTag(),
                "useConstraints": BoolTag(),
                "translationSym": BoolTag(),
                "rotationSym": BoolTag(),
                "omegaMin": FloatTag(),
                "T": FloatTag(),
                "omegaResolution": FloatTag(),
            }
        ),
        "barostat-velocity": TagContainer(
            subtags={
                "v1": FloatTag(write_tagname=False, optional=False),
                "v2": FloatTag(write_tagname=False, optional=False),
                "v3": FloatTag(write_tagname=False, optional=False),
                "v4": FloatTag(write_tagname=False, optional=False),
                "v5": FloatTag(write_tagname=False, optional=False),
                "v6": FloatTag(write_tagname=False, optional=False),
                "v7": FloatTag(write_tagname=False, optional=False),
                "v8": FloatTag(write_tagname=False, optional=False),
                "v9": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "thermostat-velocity": TagContainer(
            subtags={
                "v1": FloatTag(write_tagname=False, optional=False),
                "v2": FloatTag(write_tagname=False, optional=False),
                "v3": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "ionic-dynamics": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "B0": FloatTag(),
                "chainLengthP": FloatTag(),
                "chainLengthT": FloatTag(),
                "dt": FloatTag(),
                "nSteps": IntTag(),
                "P0": FloatTag(),  # can accept numpy.nan
                "statMethod": StrTag(options=["Berendsen", "None", "NoseHoover"]),
                "stress0": TagContainer(  # can accept numpy.nan
                    subtags={
                        "xx": FloatTag(write_tagname=False, optional=False),
                        "yy": FloatTag(write_tagname=False, optional=False),
                        "zz": FloatTag(write_tagname=False, optional=False),
                        "yz": FloatTag(write_tagname=False, optional=False),
                        "zx": FloatTag(write_tagname=False, optional=False),
                        "xy": FloatTag(write_tagname=False, optional=False),
                    }
                ),
                "T0": FloatTag(),
                "tDampP": FloatTag(),
                "tDampT": FloatTag(),
            },
        ),
    },
    "export": {
        "dump-name": TagContainer(
            allow_list_representation=True,
            subtags={
                "format": StrTag(
                    write_tagname=False,
                    optional=False,
                ),
                "freq1": StrTag(
                    options=["Init", "Ionic", "Electronic", "Fluid", "Gummel", "End"],
                    write_tagname=False,
                    optional=True,
                ),
                "format1": StrTag(
                    write_tagname=False,
                    optional=True,
                ),
                "freq2": StrTag(
                    options=["Init", "Ionic", "Electronic", "Fluid", "Gummel", "End"],
                    write_tagname=False,
                    optional=True,
                ),
                "format2": StrTag(
                    write_tagname=False,
                    optional=True,
                ),
            },
        ),
        "dump-interval": TagContainer(
            can_repeat=True,
            subtags={
                "freq": StrTag(
                    options=["Init", "Ionic", "Electronic", "Fluid", "Gummel", "End"],
                    write_tagname=False,
                    optional=False,
                ),
                "var": IntTag(write_tagname=False, optional=False),
            },
        ),
        "dump-only": BoolTag(write_value=False),
        "band-projection-params": TagContainer(
            subtags={
                "ortho": BoolTag(write_tagname=False, optional=False),
                "norm": BoolTag(write_tagname=False, optional=False),
            }
        ),
        "density-of-states": TagContainer(
            linebreak_nth_entry=1,
            subtags={
                "Total": BoolTag(write_value=False),
                "Slice": TagContainer(
                    can_repeat=True,
                    subtags={
                        "c0": FloatTag(write_tagname=False, optional=False),
                        "c1": FloatTag(write_tagname=False, optional=False),
                        "c2": FloatTag(write_tagname=False, optional=False),
                        "r": FloatTag(write_tagname=False, optional=False),
                        "i0": FloatTag(write_tagname=False, optional=False),
                        "i1": FloatTag(write_tagname=False, optional=False),
                        "i2": FloatTag(write_tagname=False, optional=False),
                    },
                ),
                "Sphere": TagContainer(
                    can_repeat=True,
                    subtags={
                        "c0": FloatTag(write_tagname=False, optional=False),
                        "c1": FloatTag(write_tagname=False, optional=False),
                        "c2": FloatTag(write_tagname=False, optional=False),
                        "r": FloatTag(write_tagname=False, optional=False),
                    },
                ),
                "AtomSlice": TagContainer(
                    can_repeat=True,
                    subtags={
                        "species": StrTag(write_tagname=False, optional=False),
                        "atomIndex": IntTag(write_tagname=False, optional=False),
                        "r": FloatTag(write_tagname=False, optional=False),
                        "i0": FloatTag(write_tagname=False, optional=False),
                        "i1": FloatTag(write_tagname=False, optional=False),
                        "i2": FloatTag(write_tagname=False, optional=False),
                    },
                ),
                "AtomSphere": TagContainer(
                    can_repeat=True,
                    subtags={
                        "species": StrTag(write_tagname=False, optional=False),
                        "atomIndex": IntTag(write_tagname=False, optional=False),
                        "r": FloatTag(write_tagname=False, optional=False),
                    },
                ),
                "File": StrTag(),
                "Orbital": TagContainer(
                    can_repeat=True,
                    subtags={
                        "species": StrTag(write_tagname=False, optional=False),
                        "atomIndex": IntTag(write_tagname=False, optional=False),
                        "orbDesc": StrTag(write_tagname=False, optional=False),
                    },
                ),
                "OrthoOrbital": TagContainer(
                    can_repeat=True,
                    subtags={
                        "species": StrTag(write_tagname=False, optional=False),
                        "atomIndex": IntTag(write_tagname=False, optional=False),
                        "orbDesc": StrTag(write_tagname=False, optional=False),
                    },
                ),
                "Etol": FloatTag(),
                "Esigma": FloatTag(),
                "EigsOverride": StrTag(),
                "Occupied": BoolTag(write_value=False),
                "Complete": BoolTag(write_value=False),
                "SpinProjected": TagContainer(
                    can_repeat=True,
                    subtags={
                        "theta": FloatTag(write_tagname=False, optional=False),
                        "phi": FloatTag(write_tagname=False, optional=False),
                    },
                ),
                "SpinTotal": BoolTag(write_value=False),
            },
        ),
        "dump-Eresolved-density": TagContainer(
            subtags={
                "Emin": FloatTag(write_tagname=False, optional=False),
                "Emax": FloatTag(write_tagname=False, optional=False),
            }
        ),
        "dump-fermi-density": MultiformatTag(
            can_repeat=True,
            format_options=[
                BoolTag(write_value=False),
                FloatTag(),
            ],
        ),
        "bgw-params": TagContainer(
            linebreak_nth_entry=1,
            multiline_tag=False,
            subtags={
                "nBandsDense": IntTag(),
                "nBandsV": IntTag(),
                "blockSize": IntTag(),
                "clusterSize": IntTag(),
                "Ecut_rALDA": FloatTag(),
                "EcutChiFluid": FloatTag(),
                "rpaExx": BoolTag(),
                "saveVxc": BoolTag(),
                "saveVxx": BoolTag(),
                "offDiagV": BoolTag(),
                "elecOnly": BoolTag(),
                "freqBroaden_eV": FloatTag(),
                "freqNimag": IntTag(),
                "freqPlasma": FloatTag(),
                "freqReMax_eV": FloatTag(),
                "freqReStep_eV": FloatTag(),
                "kernelSym_rALDA": BoolTag(),
                "kFcut_rALDA": FloatTag(),
                "q0": TagContainer(
                    subtags={
                        "q0x": FloatTag(write_tagname=False, optional=False),
                        "q0y": FloatTag(write_tagname=False, optional=False),
                        "q0z": FloatTag(write_tagname=False, optional=False),
                    }
                ),
            },
        ),
        "forces-output-coords": StrTag(options=["Cartesian", "Contravariant", "Lattice", "Positions"]),
        "polarizability": TagContainer(
            subtags={
                "eigenBasis": StrTag(
                    options=["External", "NonInteracting", "Total"],
                    write_tagname=False,
                    optional=False,
                ),
                "Ecut": FloatTag(write_tagname=False),
                "nEigs": IntTag(write_tagname=False),
            }
        ),
        "polarizability-kdiff": TagContainer(
            subtags={
                "dk0": FloatTag(write_tagname=False, optional=False),
                "dk1": FloatTag(write_tagname=False, optional=False),
                "dk2": FloatTag(write_tagname=False, optional=False),
                "dkFilenamePattern": StrTag(write_tagname=False),
            }
        ),
        "potential-subtraction": BoolTag(),
    },
    "misc": {
        "debug": StrTag(
            options=[
                "Ecomponents",
                "EigsFillings",
                "Fluid",
                "Forces",
                "KpointsBasis",
                "MuSearch",
                "Symmetries",
            ],
            can_repeat=True,
        ),
        "pcm-nonlinear-debug": TagContainer(
            subtags={
                "linearDielectric": BoolTag(write_tagname=False, optional=False),
                "linearScreening": BoolTag(write_tagname=False, optional=False),
            }
        ),
    },
}


def get_dump_tag_container() -> DumpTagContainer:
    """
    Initialize a dump tag container.

    Returns:
        DumpTagContainer: The dump tag container.
    """
    subtags2: dict[str, AbstractTag] = {}  # Called "subtags2" to avoid name conflict with the "subtags" variable
    for freq in jdftxdumpfreqoptions:
        subsubtags: dict[str, AbstractTag] = {}
        for var in jdftxdumpvaroptions:
            subsubtags[var] = BoolTag(write_value=False)
        subtags2[freq] = BoolTagContainer(subtags=subsubtags, write_tagname=True)
    return DumpTagContainer(subtags=subtags2, write_tagname=True, can_repeat=True)


MASTER_TAG_LIST["export"]["dump"] = get_dump_tag_container()


__PHONON_TAGS__: list[str] = ["phonon"]
__WANNIER_TAGS__: list[str] = [
    "wannier",
    "wannier-center-pinned",
    "wannier-dump-name",
    "wannier-initial-state",
    "wannier-minimize",
    "defect-supercell",
]
__TAG_LIST__ = [tag for group in MASTER_TAG_LIST for tag in MASTER_TAG_LIST[group]]
__TAG_GROUPS__ = {tag: group for group in MASTER_TAG_LIST for tag in MASTER_TAG_LIST[group]}


def get_tag_object(tag: str) -> AbstractTag:
    """Get the tag object for a given tag name.

    Args:
        tag (str): The tag name.

    Returns:
        AbstractTag: The tag object.
    """
    return MASTER_TAG_LIST[__TAG_GROUPS__[tag]][tag]


def get_tag_object_on_val(tag: str, val: Any) -> AbstractTag:
    """Get the tag object for a given tag name.

    Args:
        tag (str): The tag name.

    Returns:
        AbstractTag: The tag object.
    """
    tag_object = get_tag_object(tag)
    if isinstance(tag_object, MultiformatTag):
        if isinstance(val, str):
            i = tag_object.get_format_index_for_str_value(tag, val)
        else:
            i, _ = tag_object._determine_format_option(tag, val)
        tag_object = tag_object.format_options[i]
    return tag_object
