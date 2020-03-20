"""
Utility functions for assisting with creating cp2k inputs
"""

import os
import re
import numpy as np
from monty.re import regrep
from pathlib import Path

from pymatgen import SETTINGS

MODULE_DIR = Path(__file__).resolve().parent


# TODO: Setting the default basis set to triple zeta double valence potential (highest accuracy). Check this.
def get_basis_and_potential(
    species, functional="PBE", basis_type="MOLOPT", cardinality="DZVP"
):

    """
    Given a specie and a potential/basis type, this function accesses the available basis sets and potentials in
    the available "*_POTENTIAL.yaml" files. These files should come with this module distribution, but can be
    updated or remade if needed (see utils.py). Generally, the GTH potentials are used with the GTH basis sets.

    Note: as with most cp2k inputs, the convention is to use all caps, so use type="GTH" instead of "gth"

    Args:
        species: (list) list of species for which to get the potential/basis strings
        functional: (str) functional type. Default: 'PBE'
        basis_type: (str) the basis set type. Default: 'MOLOPT'
        cardinality: (str) basis set cardinality. Default: 'DZVP'

            functionals available in CP2K:
                - BLYP
                - BP
                - HCTH120
                - HCTH407
                - PADE
                - PBE
                - PBEsol
                - OLYP

    Returns:
        (dict) of the form {'specie': {'potential': potential, 'basis': basis}...}
    """
    cp2k_data = SETTINGS.get("PMG_CP2K_DATA_DIR", "")
    basis_filename = SETTINGS.get("PMG_DEFAULT_CP2K_BASIS_FILE", "BASIS_MOLOPT")
    basis_path = os.path.join(cp2k_data, basis_filename)

    potential_filename = SETTINGS.get(
        "PMG_DEFAULT_CP2K_POTENTIAL_FILE", "GTH_POTENTIALS"
    )
    potential_path = os.path.join(cp2k_data, potential_filename)

    functional = functional or SETTINGS.get("PMG_DEFAULT_FUNCTIONAL", "PBE")
    cardinality = cardinality or SETTINGS.get(
        "PMG_DEFAULT_BASIS_CARDINALITY", "DZVP"
    )

    basis_and_potential = {
        "basis_filename": basis_filename,
        "potential_filename": potential_filename,
    }

    patterns = {
        specie: re.compile(r"^[\s+]?{}\s+(.+) ".format(specie))
        for specie in species
    }
    matches_basis = regrep(basis_path, patterns=patterns)
    matches_potential = regrep(potential_path, patterns=patterns)

    for k in patterns.keys():
        basis_and_potential[k] = {}
        for m in [i[0] for i in matches_basis.get(k, [])]:
            if m[0].__contains__(cardinality.upper()):
                if m[0].__contains__(basis_type.upper()):
                    basis_and_potential[k]["basis"] = m[0]

        for m in [i[0] for i in matches_potential.get(k, [])]:
            if m[0].__contains__(functional.upper()):
                basis_and_potential[k]["potential"] = m[0]

    return basis_and_potential


def get_aux_basis(species, basis_filename=[], basis_type="cpFIT"):
    """
    Get auxiliary basis info for a list of species

    Args:
        species (list): list of species to get info for
        basis_filename (list): list of basis set filenames to look in
        basis_type (str): default basis type to look for. Otherwise, follow defaults
    """

    _aux = [
        basis_type,
        "FIT",
        "cFIT",
        "pFIT",
        "cpFIT",
        "aug-FIT",
        "aug-cFIT",
        "aug-pFIT",
        "aug-cpFIT",
    ]

    cp2k_data = SETTINGS.get("PMG_CP2K_DATA_DIR", "")
    basis_filenames = basis_filename or ["BASIS_ADMM", "BASIS_ADMM_MOLOPT"]
    basis = {"basis_filename": basis_filenames}
    basis.update({k: {} for k in species})

    for basis_filename in basis_filenames:
        basis_path = os.path.join(cp2k_data, basis_filename)
        patterns = {
            specie: re.compile(r"^[\s+]?{}\s+(.+)[\s+]?$".format(specie))
            for specie in species
        }
        matches = regrep(basis_path, patterns=patterns)

        for k in patterns.keys():
            found = False
            for a in _aux:
                for m in [i[0] for i in matches.get(k, [])]:
                    if m[0].startswith(a):
                        basis[k]["basis"] = m[0]
                        found = True
                        break
                if found:
                    break
    return basis


def get_unique_site_indices(structure):
    """
    Get unique site indices for a structure according to site properties. Whatever site-property has the most
    unique values is used for indexing
    """
    sites = {}
    property = None
    for s in structure.symbol_set:
        s_ids = structure.indices_from_symbol(s)
        unique = [0]
        for site_prop, vals in structure.site_properties.items():
            _unique = np.unique([vals[i] for i in s_ids])
            if len(unique) < len(_unique):
                unique = _unique
                property = site_prop
        if property is None:
            sites[s] = s_ids
        else:
            for i, u in enumerate(unique):
                sites[s + "_" + str(i + 1)] = []
                for j, site in zip(
                    s_ids,
                    [structure.site_properties[property][ids] for ids in s_ids],
                ):
                    if site == u:
                        sites[s + "_" + str(i + 1)].append(j)
    return sites
