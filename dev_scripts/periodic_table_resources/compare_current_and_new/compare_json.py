# compare_json.py
from __future__ import annotations

import json

from deepdiff import DeepDiff


def compare_jsons(path1: str, path2: str, log_path: str = "diff.log") -> bool:
    """
    Compare two JSON files and write a human-readable diff to a log file.

    Args:
        path1 (str): Path to the first JSON file.
        path2 (str): Path to the second JSON file.
        log_path (str): Path to the output log file.

    Returns:
        bool: True if differences were found, False if identical.
    """
    with open(path1) as f1, open(path2) as f2:
        d1 = json.load(f1)
        d2 = json.load(f2)

    diff = DeepDiff(d1, d2, ignore_order=True, significant_digits=6)

    with open(log_path, "w") as log:
        if not diff:
            log.write("✅ No differences found.\n")
            print("✅ No differences found.")
            return False

        log.write("❌ Differences found:\n")

        if "dictionary_item_added" in diff:
            log.write("\n🔼 Added items:\n")
            for item in diff["dictionary_item_added"]:
                log.write(f"  + {item}\n")

        if "dictionary_item_removed" in diff:
            log.write("\n🔽 Removed items:\n")
            for item in diff["dictionary_item_removed"]:
                log.write(f"  - {item}\n")

        if "values_changed" in diff:
            log.write("\n🌀 Changed values:\n")
            for item, change in diff["values_changed"].items():
                log.write(f"  ~ {item}: {change['old_value']} → {change['new_value']}\n")

        if "type_changes" in diff:
            log.write("\n⚠️ Type changes:\n")
            for item, change in diff["type_changes"].items():
                log.write(f"  ! {item}: {change['old_type']} → {change['new_type']}\n")

    print(f"❌ Differences written to: {log_path}")
    return True
