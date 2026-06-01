from __future__ import annotations

import os

from monty.serialization import loadfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
AFLOW_PROTOTYPE_LIBRARY = loadfn(f"{MODULE_DIR}/aflow_prototypes.json.gz")
WYCKOFF_MULTIPLICITY_DICT = loadfn(f"{MODULE_DIR}/wyckoff-position-multiplicities.json.gz")
WYCKOFF_POSITION_PARAM_DICT = loadfn(f"{MODULE_DIR}/wyckoff-position-params.json.gz")
WYCKOFF_POSITION_RELAB_DICT = {
    spg_num: [{int(key): line for key, line in val.items()} for val in vals]
    for spg_num, vals in loadfn(f"{MODULE_DIR}/wyckoff-position-relabelings.json.gz").items()
}
