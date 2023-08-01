# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

"""
Implementation for `pmg query` CLI.
"""

from __future__ import annotations

import json
import re

from monty.serialization import dumpfn
from tabulate import tabulate

from pymatgen.ext.matproj import MPRester


def do_query(args):
    """
    Perform query to the Materials Project

    Args:
        args (dict): Args from argparse.
    """
    m = MPRester()
    try:
        criteria = json.loads(args.criteria)
    except json.decoder.JSONDecodeError:
        criteria = args.criteria
    if args.structure:
        count = 0
        for d in m.query(criteria, properties=["structure", "task_id"]):
            s = d["structure"]
            formula = re.sub(r"\s+", "", s.formula)
            if args.structure == "poscar":
                fname = f"POSCAR.{d['task_id']}_{formula}"
            else:
                fname = f"{d['task_id']}-{formula}.{args.structure}"
            s.to(filename=fname)
            count += 1
        print(f"{count} structures written!")
    elif args.entries:
        entries = m.get_entries(criteria)
        dumpfn(entries, args.entries)
        print(f"{len(entries)} entries written to {args.entries}!")
    else:
        props = ["e_above_hull", "spacegroup"]
        props += args.data
        entries = m.get_entries(criteria, property_data=props)
        t = []
        headers = ["mp-id", "Formula", "Spacegroup", "E/atom (eV)", "E above hull (eV)", *args.data]
        for e in entries:
            row = [
                e.entry_id,
                e.composition.reduced_formula,
                e.data["spacegroup"]["symbol"],
                e.energy_per_atom,
                e.data["e_above_hull"],
            ]
            row += [e.data[s] for s in args.data]

            t.append(row)

        t = sorted(t, key=lambda x: x[headers.index("E above hull (eV)")])
        print(tabulate(t, headers=headers, tablefmt="pipe", floatfmt=".3f"))
