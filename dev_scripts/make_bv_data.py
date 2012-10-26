#!/usr/bin/env python

import json
import os

from pymatgen import Specie, PMGJSONEncoder
from pymatgen.analysis.bond_valence import BondValenceParam


param = []
species_set = []
with open('bvparm2011', 'r') as f:
    for l in f:
        toks = l.split()
        sp1 = Specie(toks[0], int(toks[1]))
        sp2 = Specie(toks[2], int(toks[3]))
        species = frozenset([sp1, sp2])
        if species not in species_set and toks[7] != "unchecked":
            param.append(BondValenceParam([sp1, sp2], float(toks[4]),
                                             float(toks[5])))
            species_set.append(species)

with open(os.path.join("..", "pymatgen", "analysis", "bvparm2011.json"),
          "w") as f:
    json.dump(param, f, cls=PMGJSONEncoder)
