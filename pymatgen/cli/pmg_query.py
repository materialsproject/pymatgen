# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

"""Implementation for `pmg query` CLI."""

from __future__ import annotations

import json
import re

from monty.serialization import dumpfn
from tabulate import tabulate


def do_query(args):
    """Perform query to the Materials Project.

    Args:
        args (dict): Args from argparse.
    """
    from pymatgen.ext.matproj import MPRester

    mpr = MPRester()
    try:
        criteria = json.loads(args.criteria)
    except json.decoder.JSONDecodeError:
        criteria = args.criteria
    if args.structure:
        count = 0
        for dct in mpr.query(criteria, properties=["structure", "task_id"]):
            struct = dct["structure"]
            formula = re.sub(r"\s+", "", struct.formula)
            if args.structure == "poscar":
                fname = f"POSCAR.{dct['task_id']}_{formula}"
            else:
                fname = f"{dct['task_id']}-{formula}.{args.structure}"
            struct.to(filename=fname)
            count += 1
        print(f"{count} structures written!")
    elif args.entries:
        entries = mpr.get_entries(criteria)
        dumpfn(entries, args.entries)
        print(f"{len(entries)} entries written to {args.entries}!")
    else:
        props = ["e_above_hull", "spacegroup"]
        props += args.data
        entries = mpr.get_entries(criteria, property_data=props)
        table = []
        headers = ["mp-id", "Formula", "Spacegroup", "E/atom (eV)", "E above hull (eV)", *args.data]
        for entry in entries:
            row = [
                entry.entry_id,
                entry.composition.reduced_formula,
                entry.data["spacegroup"]["symbol"],
                entry.energy_per_atom,
                entry.data["e_above_hull"],
            ]
            row += [entry.data[s] for s in args.data]

            table.append(row)

        table = sorted(table, key=lambda x: x[headers.index("E above hull (eV)")])
        print(tabulate(table, headers=headers, tablefmt="pipe", floatfmt=".3f"))
