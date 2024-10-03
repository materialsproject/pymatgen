#!/usr/bin/env python
"""
This script regenerates the enum values in pymatgen.core.libxc_func.py.
It requires in input the path of the `libxc_docs.txt` file contained in libxc/src
The script parses this file, creates a new JSON file inside pymatgen.core
and update the enum values declared in LibxcFunc.
The script must be executed inside pymatgen/dev_scripts.
"""

from __future__ import annotations

import json
import os
import sys
from copy import deepcopy


def parse_libxc_docs(path):
    """Parse libxc_docs.txt file, return dictionary {libxc_id: info_dict}."""

    def parse_section(section):
        dct = {}
        for line in section:
            key, value = line.split(":")
            dct[key.strip()] = value.strip()

        return int(dct["Number"]), dct

    dct = {}
    with open(path) as file:
        section = []
        for line in file:
            if not line.startswith("-"):
                section += [line]
            else:
                num, entry = parse_section(section)
                if num in dct:
                    raise RuntimeError(f"{num=} should not be present in {dct=}.")
                dct[num] = entry
                section = []
        if section:
            raise RuntimeError(f"Expected empty section, got {section=}")

    return dct


def write_libxc_docs_json(xc_funcs, json_path):
    """Write JSON file with libxc metadata to path jpath."""
    xc_funcs = deepcopy(xc_funcs)

    # Remove XC_FAMILY from Family and XC_ from Kind to make strings more human-readable.
    for dct in xc_funcs.values():
        dct["Family"] = dct["Family"].replace("XC_FAMILY_", "", 1)
        dct["Kind"] = dct["Kind"].replace("XC_", "", 1)

    # Build lightweight version with a subset of keys.
    for num, dct in xc_funcs.items():
        xc_funcs[num] = {key: dct[key] for key in ("Family", "Kind", "References")}
        # Descriptions are optional
        for opt in ("Description 1", "Description 2"):
            desc = dct.get(opt)
            if desc is not None:
                xc_funcs[num][opt] = desc

    with open(json_path, "w") as fh:
        json.dump(xc_funcs, fh)

    return xc_funcs


def main():
    """Main function."""
    if "-h" in sys.argv or "--help" in sys.argv:
        print(__doc__)
        print("Usage: regen_libxcfunc.py path_to_libxc_docs.txt")
        return 0

    try:
        path = sys.argv[1]
    except IndexError:
        print(__doc__)
        print("Usage: regen_libxcfunc.py path_to_libxc_docs.txt")
        return 1

    xc_funcs = parse_libxc_docs(path)

    # Generate new JSON file in pycore
    pmg_core = os.path.abspath("../pymatgen/core/")
    json_path = f"{pmg_core}/libxc_docs.json"
    write_libxc_docs_json(xc_funcs, json_path)

    # Build new enum list.
    enum_list = []
    for num, d in xc_funcs.items():
        # Remove XC_ from codename
        codename = d["Codename"][3:]
        enum_list += [f"    {codename} = {num}"]
    enum_list = "\n".join(enum_list) + "\n"

    # Re-generate enumerations.
    # [0] read py module.
    xc_funcpy_path = f"{pmg_core}/libxcfunc.py"
    with open(xc_funcpy_path) as file:
        lines = file.readlines()

    # [1] insert new enum values in list
    start = lines.index("#begin_include_dont_touch\n")
    stop = lines.index("#end_include_dont_touch\n")
    lines.insert(stop, enum_list)
    del lines[start + 1 : stop]

    # [2] write new py module
    with open(xc_funcpy_path, mode="w", encoding="utf-8") as file:
        file.writelines(lines)

    print("Files have been regenerated")
    print("Remember to update __version__ in libxcfuncs.py!")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
