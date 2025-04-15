from __future__ import annotations

import json

from monty.serialization import loadfn

from pymatgen.core import PKG_DIR

PTABLE_YAML_PATH = "periodic_table.yaml"

data = loadfn(PTABLE_YAML_PATH)

with open(f"{PKG_DIR}/core/periodic_table.json", mode="w", encoding="utf-8") as file:
    json.dump(data, file)
