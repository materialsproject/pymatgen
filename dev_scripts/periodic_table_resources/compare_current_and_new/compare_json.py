from __future__ import annotations

import json

KNOWN_DIFF = {
    "metallic_radius",  # https://github.com/materialsproject/pymatgen/pull/4380
    "mineral_hardness",  # comment dropped from value: 0.5 (graphite; diamond is 10.0)(no units) → 0.5
    "ionic_radii",  # rounding: {'2': 1.36, '3': 1.0979999999999999} → {'2': 1.36, '3': 1.098}
    "average_cationic_radius",  # rounding: 0.9674999999999999 → 0.9675
    "average_ionic_radius",  # rounding: 0.9674999999999999 → 0.9675
    "electron_affinity",  # in #4344 the tolerance was too high
    "X",  # X for some elements is NaN, and cannot be compared correctly
    "coefficient_of_linear_thermal_expansion",  # rounding: 1.42e-05 → 1.4199999999999998e-05
    "ionization_energies",  # default changed from [] to None
}


def compare_jsons(path1: str, path2: str, log_path: str = "diff.log") -> bool:
    """
    Compare two shallow JSON dictionaries and log the differences.

    Args:
        path1 (str): Path to first JSON file.
        path2 (str): Path to second JSON file.
        log_path (str): Path to output log file.

    Returns:
        bool: True if differences were found, False otherwise.
    """
    with open(path1, encoding="utf-8") as f1, open(path2, encoding="utf-8") as f2:
        d1 = json.load(f1)
        d2 = json.load(f2)

    differences = []

    # Compare all properties
    all_props = set(d1.keys()) | set(d2.keys())
    for prop in sorted(all_props):
        elements1 = d1.get(prop, {})
        elements2 = d2.get(prop, {})

        all_elements = set(elements1.keys()) | set(elements2.keys())
        for el in sorted(all_elements):
            v1 = elements1.get(el)
            v2 = elements2.get(el)

            if v1 != v2:
                differences.append((prop, el, v1, v2))

    with open(log_path, "w", encoding="utf-8") as log:
        if not differences:
            log.write("✅ No differences found.\n")
            print("✅ No differences found.")
            return False

        log.write("❌ Differences found:\n")
        for prop, el, v1, v2 in differences:
            if prop not in KNOWN_DIFF:
                log.write(f"[{prop}][{el}]: {v1} → {v2}\n")

    print(f"❌ Differences written to {log_path}")
    return True
