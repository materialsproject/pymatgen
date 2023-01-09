#!/usr/bin/env python
"""
This script regenerates the enum values in pymatgen.core.libxc_func.py.
It requires in input the path of the `libxc_docs.txt` file contained in libxc/src
The script parses this file, creates a new json file inside pymatgen.core
and update the enum values declared in LibxcFunc.
The script must be executed inside pymatgen/dev_scripts.
"""

from __future__ import annotations

import json
import os
import sys


def parse_libxc_docs(path):
    """
    Parse libxc_docs.txt file, return dictionary with mapping:
    libxc_id --> info_dict
    """

    def parse_section(section):
        d = {}
        for l in section:
            key, value = l.split(":")
            key = key.strip()
            d[key] = value.strip()

        return int(d["Number"]), d

    d = {}
    with open(path) as fh:
        section = []
        for line in fh:
            if not line.startswith("-"):
                section.append(line)
            else:
                num, entry = parse_section(section)
                assert num not in d
                d[num] = entry
                section = []
        assert not section

    return d


def write_libxc_docs_json(xcfuncs, jpath):
    """Write json file with libxc metadata to path jpath."""
    from copy import deepcopy

    xcfuncs = deepcopy(xcfuncs)

    # Remove XC_FAMILY from Family and XC_ from Kind to make strings more human-readable.
    for d in xcfuncs.values():
        d["Family"] = d["Family"].replace("XC_FAMILY_", "", 1)
        d["Kind"] = d["Kind"].replace("XC_", "", 1)

    # Build lightweight version with a subset of keys.
    for num, d in xcfuncs.items():
        xcfuncs[num] = {k: d[k] for k in ("Family", "Kind", "References")}
        # Descriptions are optional
        for opt in ("Description 1", "Description 2"):
            desc = d.get(opt)
            if desc is not None:
                xcfuncs[num][opt] = desc

    with open(jpath, "w") as fh:
        json.dump(xcfuncs, fh)

    return xcfuncs


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

    xcfuncs = parse_libxc_docs(path)

    # Generate new json file in pycore
    pycore = os.path.abspath("../pymatgen/core/")
    jpath = os.path.join(pycore, "libxc_docs.json")
    write_libxc_docs_json(xcfuncs, jpath)

    # Build new enum list.
    enum_list = []
    for num, d in xcfuncs.items():
        # Remove XC_ from codename
        codename = d["Codename"][3:]
        enum_list.append(f"    {codename} = {num}")
    enum_list = "\n".join(enum_list) + "\n"

    # Re-generate enumerations.
    # [0] read py module.
    xcfuncpy_path = os.path.join(pycore, "libxcfunc.py")
    with open(xcfuncpy_path) as fh:
        lines = fh.readlines()

    # [1] insert new enum values in list
    start = lines.index("#begin_include_dont_touch\n")
    stop = lines.index("#end_include_dont_touch\n")
    lines.insert(stop, enum_list)
    del lines[start + 1 : stop]

    # [2] write new py module
    with open(xcfuncpy_path, "w") as fh:
        fh.writelines(lines)

    print("Files have been regenerated")
    print("Remember to update libxc_version in libxcfuncs.py!")

    return 0


if __name__ == "__main__":
    sys.exit(main())
