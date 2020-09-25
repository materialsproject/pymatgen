# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.


"""
Implementation for `pmg query` CLI.
"""

import json
import re

from tabulate import tabulate

from monty.serialization import dumpfn

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
                fname = "POSCAR.%s_%s" % (d["task_id"], formula)
            else:
                fname = "%s-%s.%s" % (d["task_id"], formula, args.structure)
            s.to(filename=fname)
            count += 1
        print("%d structures written!" % count)
    elif args.entries:
        entries = m.get_entries(criteria)
        dumpfn(entries, args.entries)
        print("%d entries written to %s!" % (len(entries), args.entries))
    else:
        props = ["e_above_hull", "spacegroup"]
        props += args.data
        entries = m.get_entries(criteria, property_data=props)
        t = []
        headers = ["mp-id", "Formula", "Spacegroup", "E/atom (eV)",
                   "E above hull (eV)"] + args.data
        for e in entries:
            row = [e.entry_id, e.composition.reduced_formula,
                   e.data["spacegroup"]["symbol"],
                   e.energy_per_atom, e.data["e_above_hull"]]
            row += [e.data[s] for s in args.data]

            t.append(row)

        t = sorted(t, key=lambda x: x[headers.index("E above hull (eV)")])
        print(tabulate(t, headers=headers, tablefmt="pipe", floatfmt=".3f"))
