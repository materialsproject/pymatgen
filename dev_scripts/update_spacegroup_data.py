"""Script to update symm_ops.json and symm_data.yaml in symmetry module due to issues #3845 and #3862.
symm_ops.json:
- adds Hermann_mauguin point group key and short Hermann Mauguin space group symbol
- converts screw axis notation to symm_data standard
symm_data.json
- removes mapping of rhombohedral space group types onto symbol + appended H
- replaces I/P2_12_121 key with I/P2_12_12_1
"""

from __future__ import annotations

import sys

from monty.serialization import dumpfn, loadfn

from pymatgen.symmetry.groups import PointGroup

__author__ = "Katharina Ueltzen @kaueltzen"
__date__ = "2024-06-06"

SYMM_OPS = loadfn("../src/pymatgen/symmetry/symm_ops.json")
SYMM_DATA = loadfn("../src/pymatgen/symmetry/symm_data.json")


def convert_symmops_to_sg_encoding(symbol: str) -> str:
    """
    Utility function to convert SYMMOPS space group type symbol notation
    into SYMM_DATA["space_group_encoding"] key notation with underscores before
    translational part of screw axes.
    Args:
        symbol (str): "hermann_mauguin" or "universal_h_m" key of symmops.json
    Returns:
        str: symbol in the format of SYMM_DATA["space_group_encoding"] keys
    """
    symbol_representation = symbol.split(":")
    representation = ":" + "".join(symbol_representation[1].split(" ")) if len(symbol_representation) > 1 else ""

    blickrichtungen = symbol_representation[0].split(" ")
    blickrichtungen_new = []
    for br in blickrichtungen:
        if len(br) > 1 and br[0].isdigit() and br[1].isdigit():
            blickrichtungen_new.append(br[0] + "_" + br[1:])
        else:
            blickrichtungen_new.append(br)
    return "".join(blickrichtungen_new) + representation


def remove_identity_from_full_hermann_mauguin(symbol: str) -> str:
    """
    Utility function to remove identity along blickrichtung (except in P1).
    Args:
        symbol (str): "hermann_mauguin" key of symmops.json
    Returns:
        str: short "hermann_mauguin" key
    """
    if symbol in ("P 1", "C 1", "P 1 "):
        return symbol
    blickrichtungen = symbol.split(" ")
    if blickrichtungen[1].startswith("3"):
        return symbol
    blickrichtungen_new = [br + " " for br in blickrichtungen if br != "1"]
    return "".join(blickrichtungen_new)


new_symm_data = {}
for k, v in SYMM_DATA["space_group_encoding"].items():
    if k.endswith("H"):
        new_symm_data[k.removesuffix("H")] = v
    elif k == "I2_12_121":
        new_symm_data["I2_12_12_1"] = v
    elif k == "P2_12_121":
        new_symm_data["P2_12_12_1"] = v
    else:
        new_symm_data[k] = v

SYMM_DATA["space_group_encoding"] = new_symm_data

for sg_symbol in SYMM_DATA["space_group_encoding"]:
    SYMM_DATA["space_group_encoding"][sg_symbol]["point_group"] = PointGroup.from_space_group(
        sg_symbol=sg_symbol
    ).symbol

for spg in SYMM_OPS:
    if "(" in spg["hermann_mauguin"]:
        spg["hermann_mauguin"] = spg["hermann_mauguin"].split("(")[0]

    short_h_m = remove_identity_from_full_hermann_mauguin(spg["hermann_mauguin"])
    spg["short_h_m"] = convert_symmops_to_sg_encoding(short_h_m)
    spg["hermann_mauguin_u"] = convert_symmops_to_sg_encoding(spg["hermann_mauguin"])

for spg_idx, spg in enumerate(SYMM_OPS):
    try:
        pg = PointGroup.from_space_group(spg["short_h_m"])
    except AssertionError as e:
        print(spg, str(e))
        sys.exit(1)
    SYMM_OPS[spg_idx]["point_group"] = pg.symbol

dumpfn(SYMM_DATA, "../src/pymatgen/symmetry/symm_data.json")
dumpfn(SYMM_OPS, "../src/pymatgen/symmetry/symm_ops.json")
