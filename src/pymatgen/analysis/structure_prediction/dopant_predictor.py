"""Predicting potential dopants."""

from __future__ import annotations

import warnings

import numpy as np

from pymatgen.analysis.structure_prediction.substitution_probability import SubstitutionPredictor
from pymatgen.core import Element, Species


def get_dopants_from_substitution_probabilities(
    structure, num_dopants=5, threshold=0.001, match_oxi_sign=False
) -> dict:
    """Get dopant suggestions based on substitution probabilities.

    Args:
        structure (Structure): A pymatgen structure decorated with
            oxidation states.
        num_dopants (int): The number of suggestions to return for
            n- and p-type dopants.
        threshold (float): Probability threshold for substitutions.
        match_oxi_sign (bool): Whether to force the dopant and original species
            to have the same sign of oxidation state. E.g. If the original site
            is in a negative charge state, then only negative dopants will be
            returned.

    Returns:
        dict: Dopant suggestions, given as a dictionary with keys "n_type" and
            "p_type". The suggestions for each doping type are given as a list of
            dictionaries, each with they keys:

            - "probability": The probability of substitution.
            - "dopant_species": The dopant species.
            - "original_species": The substituted species.
    """
    els_have_oxi_states = [hasattr(s, "oxi_state") for s in structure.species]

    if not all(els_have_oxi_states):
        raise ValueError("All sites in structure must have oxidation states to predict dopants.")

    sp = SubstitutionPredictor(threshold=threshold)

    subs = [sp.list_prediction([s]) for s in set(structure.species)]
    subs = [
        {
            "probability": pred["probability"],
            "dopant_species": next(iter(pred["substitutions"])),
            "original_species": next(iter(pred["substitutions"].values())),
        }
        for species_preds in subs
        for pred in species_preds
    ]
    subs.sort(key=lambda x: x["probability"], reverse=True)

    return _get_dopants(subs, num_dopants, match_oxi_sign)


def get_dopants_from_shannon_radii(bonded_structure, num_dopants=5, match_oxi_sign=False):
    """Get dopant suggestions based on Shannon radii differences.

    Args:
        bonded_structure (StructureGraph): A pymatgen structure graph
            decorated with oxidation states. For example, generated using the
            CrystalNN.get_bonded_structure() method.
        num_dopants (int): The number of suggestions to return for
            n- and p-type dopants.
        match_oxi_sign (bool): Whether to force the dopant and original species
            to have the same sign of oxidation state. E.g. If the original site
            is in a negative charge state, then only negative dopants will be
            returned.

    Returns:
        dict: Dopant suggestions, given as a dictionary with keys "n_type" and
            "p_type". The suggestions for each doping type are given as a list of
            dictionaries, each with they keys:

            - "radii_diff": The difference between the Shannon radii of the species.
            - "dopant_species": The dopant species.
            - "original_species": The substituted species.
    """
    # get a list of all Species for all elements in all their common oxidation states
    all_species = [Species(el.symbol, oxi) for el in Element for oxi in el.common_oxidation_states]

    # get a series of tuples with (coordination number, specie)
    cn_and_species = {
        (
            bonded_structure.get_coordination_of_site(idx),
            bonded_structure.structure[idx].specie,
        )
        for idx in range(len(bonded_structure))
    }

    cn_to_radii_map = {}
    possible_dopants = []

    for cn, species in cn_and_species:
        cn_roman = _int_to_roman(cn)

        try:
            species_radius = species.get_shannon_radius(cn_roman)
        except KeyError:
            warnings.warn(
                f"Shannon radius not found for {species} with coordination number {cn}.\nSkipping...", stacklevel=2
            )
            continue

        if cn not in cn_to_radii_map:
            cn_to_radii_map[cn] = _shannon_radii_from_cn(all_species, cn_roman, radius_to_compare=species_radius)

        shannon_radii = cn_to_radii_map[cn]

        possible_dopants += [
            {
                "radii_diff": p["radii_diff"],
                "dopant_species": p["species"],
                "original_species": species,
            }
            for p in shannon_radii
        ]

    possible_dopants.sort(key=lambda x: abs(x["radii_diff"]))

    return _get_dopants(possible_dopants, num_dopants, match_oxi_sign)


def _get_dopants(substitutions, num_dopants, match_oxi_sign) -> dict:
    """Utility method to get n- and p-type dopants from a list of substitutions."""
    n_type = [
        pred
        for pred in substitutions
        if pred["dopant_species"].oxi_state > pred["original_species"].oxi_state
        and (
            not match_oxi_sign
            or np.sign(pred["dopant_species"].oxi_state) == np.sign(pred["original_species"].oxi_state)
        )
    ]
    p_type = [
        pred
        for pred in substitutions
        if pred["dopant_species"].oxi_state < pred["original_species"].oxi_state
        and (
            not match_oxi_sign
            or np.sign(pred["dopant_species"].oxi_state) == np.sign(pred["original_species"].oxi_state)
        )
    ]

    return {"n_type": n_type[:num_dopants], "p_type": p_type[:num_dopants]}


def _shannon_radii_from_cn(species_list, cn_roman, radius_to_compare=0):
    """
    Utility func to get Shannon radii for a particular coordination number.

    As the Shannon radii depends on charge state and coordination number,
    species without an entry for a particular coordination number will
    be skipped.

    Args:
        species_list (list): A list of Species to get the Shannon radii for.
        cn_roman (str): The coordination number as a roman numeral. See
            Species.get_shannon_radius for more details.
        radius_to_compare (float, optional): If set, the data will be returned
            with a "radii_diff" key, containing the difference between the
            shannon radii and this radius.

    Returns:
        list[dict]: The Shannon radii for all Species in species. Formatted
            as a list of dictionaries, with the keys:

            - "species": The species with charge state.
            - "radius": The Shannon radius for the species.
            - "radius_diff": The difference between the Shannon radius and the
                radius_to_compare optional argument.
    """
    shannon_radii = []

    for s in species_list:
        try:
            radius = s.get_shannon_radius(cn_roman)
            shannon_radii.append(
                {
                    "species": s,
                    "radius": radius,
                    "radii_diff": radius - radius_to_compare,
                }
            )
        except KeyError:
            pass

    return shannon_radii


def _int_to_roman(number):
    """Utility method to convert an int (less than 20) to a roman numeral."""
    roman_conv = [(10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]

    result = []
    for arabic, roman in roman_conv:
        factor, number = divmod(number, arabic)
        result.append(roman * factor)
        if number == 0:
            break
    return "".join(result)
