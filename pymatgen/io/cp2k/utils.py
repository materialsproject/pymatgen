"""
Utility functions for assisting with cp2k IO
"""

import os
import re
import warnings
from pathlib import Path

import numpy as np
from monty.io import zopen
from ruamel.yaml import YAML

from pymatgen.core import SETTINGS

MODULE_DIR = Path(__file__).resolve().parent


def _postprocessor(s):
    """
    Helper function to post process the results of the pattern matching functions in Cp2kOutput and turn them to
    python types.
    """
    s = s.rstrip()  # Remove leading/trailing whitespace
    s = s.replace(" ", "_")  # Remove whitespaces

    if s.lower() == "no":
        return False
    if s.lower() == "none":
        return None
    if s.lower() == "yes" or s.lower() == "true":
        return True
    if re.match(r"^-?\d+$", s):
        try:
            return int(s)
        except ValueError:
            raise OSError("Error in parsing CP2K file.")
    if re.match(r"^[+\-]?(?=.)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?$", s):
        try:
            return float(s)
        except ValueError:
            raise OSError("Error in parsing CP2K file.")
    if re.match(r"\*+", s):
        try:
            return np.NaN
        except ValueError:
            raise OSError("Error in parsing CP2K file.")
    return s


def _preprocessor(s, d="."):
    """
    Cp2k contains internal preprocessor flags that are evaluated before
    execution. This helper function recognizes those preprocessor flags
    and replacees them with an equivalent cp2k input (this way everything
    is contained neatly in the cp2k input structure, even if the user
    preferred to use the flags.

    CP2K preprocessor flags (with arguments) are:

        @INCLUDE FILENAME: Insert the contents of FILENAME into the file at
            this location.
        @SET VAR VALUE: set a variable, VAR, to have the value, VALUE.
        $VAR or ${VAR}: replace these with the value of the variable, as set
            by the @SET flag.
        @IF/@ELIF: Not implemented yet.

    Args:
        s (str): string representation of cp2k input to preprocess
    """
    includes = re.findall(r"(@include.+)", s, re.IGNORECASE)
    for incl in includes:
        inc = incl.split()
        assert len(inc) == 2  # @include filename
        inc = inc[1].strip("'")
        inc = inc.strip('"')
        with zopen(os.path.join(d, inc)) as f:
            s = re.sub(rf"{incl}", f.read(), s)
    variable_sets = re.findall(r"(@SET.+)", s, re.IGNORECASE)
    for match in variable_sets:
        v = match.split()
        assert len(v) == 3  # @SET VAR value
        var, value = v[1:]
        s = re.sub(rf"{match}", "", s)
        s = re.sub(r"\${?" + var + "}?", value, s)

    c1 = re.findall(r"@IF", s, re.IGNORECASE)
    c2 = re.findall(r"@ELIF", s, re.IGNORECASE)
    if len(c1) > 0 or len(c2) > 0:
        raise NotImplementedError("This cp2k input processor does not currently support conditional blocks.")
    return s


def natural_keys(text):
    """
    Sort text by numbers coming after an underscore with natural number
    convention,
    Ex: [file_1, file_12, file_2] becomes [file_1, file_2, file_12]
    """

    def atoi(t):
        return int(t) if t.isdigit() else t

    return [atoi(c) for c in re.split(r"_(\d+)", text)]


def get_basis_and_potential(species, basis_and_potential_map, functional="PBE"):
    """
    Retrieve basis/potential dictionary

    Args:
        species: (list) list of species for which to get the potential/basis strings
        basis_and_potential_map: (dict or str) a keyword string or a dictionary
            specifying how bases and/or potentials should be assigned
        functional: (str) functional type. Default: 'PBE'

    Returns:
        (dict) of the form {'specie': {'potential': potential, 'basis': basis}...}
    """

    potential_filename = SETTINGS.get("PMG_DEFAULT_CP2K_POTENTIAL_FILE", "GTH_POTENTIALS")
    basis_filenames = ["BASIS_MOLOPT", "BASIS_MOLOPT_UCL"]

    functional = functional or SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")

    basis_and_potential = {
        "basis_filenames": basis_filenames,
        "potential_filename": potential_filename,
    }

    with open(os.path.join(MODULE_DIR, "settings.yaml")) as f:
        yaml = YAML(typ="unsafe", pure=True)
        settings = yaml.load(f)

    if basis_and_potential_map == "best":
        basis_and_potential.update(
            {
                s: {
                    "basis": settings[s]["basis_sets"]["best_basis"],
                    "potential": [p for p in settings[s]["potentials"]["gth_potentials"] if functional in p][0],
                }
                for s in species
            }
        )
    elif basis_and_potential_map == "preferred":
        basis_and_potential.update(
            {
                s: {
                    "basis": settings[s]["basis_sets"]["preferred_basis"],
                    "potential": [p for p in settings[s]["potentials"]["gth_potentials"] if functional in p][0],
                }
                for s in species
            }
        )
    else:
        basis_and_potential.update(basis_and_potential_map)

    return basis_and_potential


def get_aux_basis(basis_type, default_basis_type="cpFIT"):
    """
    Get auxiliary basis info for a list of species.

    Args:
        basis_type (dict): dict of auxiliary basis sets to use. i.e:
            basis_type = {'Si': 'cFIT', 'O': 'cpFIT'}. Basis type needs to
            exist for that species.

            Basis types:
                FIT
                cFIT
                pFIT
                cpFIT
                GTH-def2
                aug-{FIT,cFIT,pFIT,cpFIT, GTH-def2}

        default_basis_type (str) default basis type if n

    """
    with open(os.path.join(MODULE_DIR, "settings.yaml")) as f:
        yaml = YAML(typ="unsafe", pure=True)
        settings = yaml.load(f)
        aux_bases = {
            s: [settings[s]["basis_sets"]["preferred_aux_basis"]]
            if "preferred_aux_basis" in settings[s]["basis_sets"]
            else settings[s]["basis_sets"]["aux_basis"]
            for s in settings
        }

    default_basis_type = default_basis_type or SETTINGS.get("PMG_CP2K_DEFAULT_AUX_BASIS_TYPE")
    basis_type = {k: basis_type[k] if basis_type[k] else default_basis_type for k in basis_type}
    basis = {k: {} for k in basis_type}
    for k in basis_type:
        if aux_bases.get(k) is None:
            basis[k] = None
            continue
        for i in aux_bases[k]:
            if i.startswith(basis_type[k]):
                basis[k] = i
                break
    for k in basis:
        if not basis[k]:
            if aux_bases[k]:
                basis[k] = aux_bases[k][0]
            else:
                raise LookupError(f"No basis of that type found for: {k}")
    return basis


def get_unique_site_indices(structure):
    """
    Get unique site indices for a structure according to site properties. Whatever site-property has the most
    unique values is used for indexing.

    For example, if you have magnetic CoO with half Co atoms having a positive moment, and the other
    half having a negative moment. Then this function will create a dict of sites for Co_1, Co_2, O. This
    function also deals with "Species" properties like oxi_state and spin by pushing them to site
    properties.

    This creates unique sites, based on site properties, but does not have anything to do with turning
    those site properties into CP2K input parameters. This will only be done for properties which can be
    turned into CP2K input parameters, which are stored in parsable_site_properties.
    """
    spins = []
    oxi_states = []
    parsable_site_properties = {"magmom", "oxi_state", "spin", "u_minus_j", "basis", "potential", "ghost", "aux_basis"}

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
    _sites = {u: [] for u in unique_itms}
    for i, itm in enumerate(items):
        _sites[itm].append(i)
    sites = {}
    nums = {s: 1 for s in structure.symbol_set}
    for s in _sites:
        sites[f"{s[0]}_{nums[s[0]]}"] = _sites[s]
        nums[s[0]] += 1

    return sites


def get_cutoff_from_basis(els, bases, rel_cutoff=50):
    """
    Gets the appropriate cutoff for the calculation given the elements/basis sets being used
    and the desired relative cutoff.

    Args:
        els: list of element symbols
        bases: corresponding basis set names for the elements
        rel_cutoff: The desired relative cutoff
    Returns:
        Ideal cutoff for calculation.
    """
    with open(os.path.join(MODULE_DIR, "settings.yaml")) as f:
        yaml = YAML(typ="unsafe", pure=True)
        _exponents = yaml.load(f)
        _exponents = {
            k: v["basis_sets"].get("basis_set_largest_exponents")
            for k, v in _exponents.items()
            if v["basis_sets"].get("basis_set_largest_exponents")
        }
        exponents = {el.upper(): {b.upper(): v for b, v in basis.items()} for el, basis in _exponents.items()}
        return max(np.ceil(exponents[el.upper()][basis.upper()]) * rel_cutoff for el, basis in zip(els, bases))


# TODO this is not comprehensive. There are so many libxc functionals (e.g. see r2scan)
# and their availability changes with libxc version. I see no easy way to deal with this
# So these are just a few of the most common ones.
# TODO whenever ADMM gets interfaced with LibXC, this should include hybrid functionals
# and the sets will get more streamlined.
def get_xc_functionals(name):
    """
    Get the XC functionals for a given functional name. This utility does
    not deal with defining hybrid functionals since those may or may not
    require ADMM, which is not supported by LibXC and so needs to be manually
    defined.

    Args:
        name: Name of the functional.
    """
    name = name.upper()
    if name == "PBE":
        return ["PBE"]
    if name in ("LDA", "PADE"):
        return ["PADE"]
    if name == "B3LYP":
        return ["B3LYP"]
    if name == "BLYP":
        return ["BLYP"]
    if name == "SCAN":
        return ["MGGA_X_SCAN", "MGGA_C_SCAN"]
    if name == "SCANL":
        return ["MGGA_X_SCANL", "MGGA_C_SCANL"]
    if name == "R2SCAN":
        return ["MGGA_X_R2SCAN", "MGGA_C_R2SCAN"]
    if name == "R2SCANL":
        return ["MGGA_X_R2SCANL", "MGGA_C_R2SCANL"]
    warnings.warn(f"Unknown XC functionals: {name}")
    return [name]


def get_truncated_coulomb_cutoff(inp_struct):
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
