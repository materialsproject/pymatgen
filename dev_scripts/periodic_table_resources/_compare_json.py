#!/usr/bin/env python3

"""Temporary dev script to compare two periodic table JSONs to make sure they're exactly the same.

TODO: remove after finalizing JSON generator.

Known differences:
- New JSON has extra properties (present in current YAML, missing in JSON): "Ground level", "NMR Quadrupole Moment"
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from pymatgen.core import ROOT

NEW_JSON = f"{Path(__file__).parent}/_periodic_table.json"
OLD_JSON = f"{ROOT}/pymatgen/core/periodic_table.json"

known_diff_properties: tuple[str, ...] = (
    "Electrical resistivity",  # TODO: current recording would add an extra space, e.g.:
    # Old: &gt; 10<sup>23</sup>10<sup>-8</sup> &Omega; m (str)
    # New: &gt; 10<sup>23</sup> 10<sup>-8</sup> &Omega; m (str)
)

ABS_TOL: float = 0.001


def main():
    with open(OLD_JSON, encoding="utf-8") as f:
        old_data = json.load(f)

    with open(NEW_JSON, encoding="utf-8") as f:
        new_data = json.load(f)

    for element, old_props in old_data.items():
        print(f"🔍 Checking element: {element}")

        if element not in new_data:
            print(f"❌ Missing element in new JSON: {element}")
            continue

        new_props = new_data[element]
        for prop_name, old_value in old_props.items():
            if prop_name in known_diff_properties:
                continue

            # Old JSON use "no data" as placeholder, new JSON just dropped NaN
            if isinstance(old_value, str) and old_value.startswith("no data"):
                continue

            # Skip empty values
            if old_value is None or not old_value:
                continue

            if prop_name not in new_props:
                print(f"  ❌ Missing property '{prop_name}' in new JSON")
                print(f"       Old: {old_value}")
                continue

            new_value = new_props[prop_name]

            # Try coercion for numeric strings
            try:
                old_coerced = float(old_value.split()[0]) if isinstance(old_value, str) else old_value
                new_coerced = float(new_value.split()[0]) if isinstance(new_value, str) else new_value
            except (ValueError, TypeError):
                old_coerced = old_value
                new_coerced = new_value

            # Numeric tolerance comparison
            if (
                isinstance(old_coerced, float | int)
                and isinstance(new_coerced, float | int)
                and math.isclose(old_coerced, new_coerced, abs_tol=ABS_TOL)
            ):
                continue  # Close enough

            # Fallback strict comparison
            if old_coerced != new_coerced:
                print(f"  ❌ Mismatch in property '{prop_name}':")
                print(f"       Old: {old_value} ({type(old_value).__name__})")
                print(f"       New: {new_value} ({type(new_value).__name__})")


if __name__ == "__main__":
    main()
