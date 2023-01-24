"""
Utility functions for assisting with cp2k IO
"""

from __future__ import annotations

import os
import re

import numpy as np
from monty.io import zopen

from pymatgen.core import Molecule, Structure


def postprocessor(data: str) -> str | int | float | bool | None:
    """
    Helper function to post process the results of the pattern matching functions in Cp2kOutput
    and turn them to Python types.

    Args:
        data (str): The data to be post processed.

    Raises:
        ValueError: If the data cannot be parsed.

    Returns:
        str | int | float | bool | None: The post processed data.
    """
    data = data.strip().replace(" ", "_")  # remove leading/trailing whitespace, replace spaces with _

    if data.lower() in ("false", "no", "f"):
        return False
    if data.lower() == "none":
        return None
    if data.lower() in ("true", "yes", "t"):
        return True
    if re.match(r"^-?\d+$", data):
        try:
            return int(data)
        except ValueError as exc:
            raise ValueError(f"Error parsing {data!r} as int in CP2K file.") from exc
    if re.match(r"^[+\-]?(?=.)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?$", data):
        try:
            return float(data)
        except ValueError as exc:
            raise ValueError(f"Error parsing {data!r} as float in CP2K file.") from exc
    if re.match(r"\*+", data):
        return np.NaN
    return data


def preprocessor(data: str, dir: str = ".") -> str:
    """
    Cp2k contains internal preprocessor flags that are evaluated before execution. This helper
    function recognizes those preprocessor flags and replaces them with an equivalent cp2k input
    (this way everything is contained neatly in the cp2k input structure, even if the user preferred
    to use the flags.

    CP2K preprocessor flags (with arguments) are:

        @INCLUDE FILENAME: Insert the contents of FILENAME into the file at
            this location.
        @SET VAR VALUE: set a variable, VAR, to have the value, VALUE.
        $VAR or ${VAR}: replace these with the value of the variable, as set
            by the @SET flag.
        @IF/@ELIF: Not implemented yet.

    Args:
        data (str): cp2k input to preprocess
        dir (str, optional): Path for include files. Default is '.' (current directory).

    Returns:
        Preprocessed string
    """
    includes = re.findall(r"(@include.+)", data, re.IGNORECASE)
    for incl in includes:
        inc = incl.split()
        assert len(inc) == 2  # @include filename
        inc = inc[1].strip("'")
        inc = inc.strip('"')
        with zopen(os.path.join(dir, inc)) as f:
            data = re.sub(rf"{incl}", f.read(), data)
    variable_sets = re.findall(r"(@SET.+)", data, re.IGNORECASE)
    for match in variable_sets:
        v = match.split()
        assert len(v) == 3  # @SET VAR value
        var, value = v[1:]
        data = re.sub(rf"{match}", "", data)
        data = re.sub(r"\${?" + var + "}?", value, data)

    c1 = re.findall(r"@IF", data, re.IGNORECASE)
    c2 = re.findall(r"@ELIF", data, re.IGNORECASE)
    if len(c1) > 0 or len(c2) > 0:
        raise NotImplementedError("This cp2k input processor does not currently support conditional blocks.")
    return data


def chunk(string: str):
    """
    Chunk the string from a cp2k basis or potential file
    """
    lines = iter(line for line in (l.strip() for l in string.split("\n")) if line and not line.startswith("#"))
    chunks: list = []
    for line in lines:
        if line.split()[0].isalpha():
            chunks.append([])
        chunks[-1].append(line)
    chunks = ["\n".join(c) for c in chunks]
    return chunks


def natural_keys(text: str):
    """
    Sort text by numbers coming after an underscore with natural number
    convention,
    Ex: [file_1, file_12, file_2] becomes [file_1, file_2, file_12]
    """

    def atoi(t):
        return int(t) if t.isdigit() else t

    return [atoi(c) for c in re.split(r"_(\d+)", text)]


def get_unique_site_indices(structure: Structure | Molecule):
    """
    Get unique site indices for a structure according to site properties. Whatever site-property
    has the most unique values is used for indexing.

    For example, if you have magnetic CoO with half Co atoms having a positive moment, and the
    other half having a negative moment. Then this function will create a dict of sites for
    Co_1, Co_2, O. This function also deals with "Species" properties like oxi_state and spin by
    pushing them to site properties.

    This creates unique sites, based on site properties, but does not have anything to do with
    turning those site properties into CP2K input parameters. This will only be done for properties
    which can be turned into CP2K input parameters, which are stored in parsable_site_properties.
    """
    spins = []
    oxi_states = []
    parsable_site_properties = {
        "magmom",
        "oxi_state",
        "spin",
        "u_minus_j",
        "basis",
        "potential",
        "ghost",
        "aux_basis",
    }

    for site in structure:
        for sp in site.species:
            oxi_states.append(getattr(sp, "oxi_state", 0))
            spins.append(getattr(sp, "_properties", {}).get("spin", 0))

    structure.add_site_property("oxi_state", oxi_states)
    structure.add_site_property("spin", spins)
    structure.remove_oxidation_states()
    items = [
        (
            site.species_string,
            *[
                structure.site_properties[k][i]
                for k in structure.site_properties
                if k.lower() in parsable_site_properties
            ],
        )
        for i, site in enumerate(structure)
    ]
    unique_itms = list(set(items))
    _sites: dict[tuple, list] = {u: [] for u in unique_itms}
    for i, itm in enumerate(items):
        _sites[itm].append(i)
    sites = {}
    nums = {s: 1 for s in structure.symbol_set}
    for s in _sites:
        sites[f"{s[0]}_{nums[s[0]]}"] = _sites[s]
        nums[s[0]] += 1

    return sites


def get_truncated_coulomb_cutoff(inp_struct: Structure):
    """
    Get the truncated Coulomb cutoff for a given structure.
    """
    m = inp_struct.lattice.matrix
    m = (abs(m) > 1e-5) * m
    a, b, c = m[0], m[1], m[2]
    x = abs(np.dot(a, np.cross(b, c)) / np.linalg.norm(np.cross(b, c)))
    y = abs(np.dot(b, np.cross(a, c)) / np.linalg.norm(np.cross(a, c)))
    z = abs(np.dot(c, np.cross(a, b)) / np.linalg.norm(np.cross(a, b)))
    return np.floor(100 * min([x, y, z]) / 2) / 100
