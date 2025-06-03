"""
This module provides conversion between structure entries following the
OPTIMADE (https://optimade.org) standard and pymatgen Structure objects.

The code is adapted from the `optimade.adapters.structures.pymatgen` module in
optimade-python-tools (https://github.com/Materials-Consortia/optimade-python-tools),
and aims to work without requiring the explicit installation of the `optimade-python-tools`.

"""

from __future__ import annotations

import itertools
import json
import math
import re
from functools import reduce
from typing import TYPE_CHECKING

import orjson

from pymatgen.core.structure import Lattice, Structure

if TYPE_CHECKING:
    from collections.abc import Generator
    from typing import Any

    from pymatgen.core.structure import IStructure

__author__ = "Matthew Evans"


def _pymatgen_species(
    nsites: int,
    species_at_sites: list[str],
) -> list[dict[str, float]]:
    """Create list of {"symbol": "concentration"} per site for constructing pymatgen Species objects.
    Removes vacancies, if they are present.

    This function is adapted from the `optimade.adapters.structures.pymatgen` module in `optimade-python-tools`,
    with some of the generality removed (in terms of partial occupancy).

    """
    species = [{"name": _, "concentration": [1.0], "chemical_symbols": [_]} for _ in set(species_at_sites)]
    species_dict = {_["name"]: _ for _ in species}

    pymatgen_species = []
    for site_number in range(nsites):
        species_name = species_at_sites[site_number]
        current_species = species_dict[species_name]

        chemical_symbols = []
        concentration = []
        for index, symbol in enumerate(current_species["chemical_symbols"]):
            if symbol == "vacancy":
                # Skip. This is how pymatgen handles vacancies;
                # to not include them, while keeping the concentration in a site less than 1.
                continue
            chemical_symbols.append(symbol)
            concentration.append(current_species["concentration"][index])

        pymatgen_species.append(dict(zip(chemical_symbols, concentration, strict=True)))

    return pymatgen_species


def _optimade_anonymous_element_generator() -> Generator[str, None, None]:
    """Generator that yields the next symbol in the A, B, Aa, ... Az OPTIMADE anonymous
    element naming scheme.

    """
    from string import ascii_lowercase

    for size in itertools.count(1):
        for tuple_strings in itertools.product(ascii_lowercase, repeat=size):
            list_strings = list(tuple_strings)
            list_strings[0] = list_strings[0].upper()
            yield "".join(list_strings)


def _optimade_reduce_or_anonymize_formula(formula: str, alphabetize: bool = True, anonymize: bool = False) -> str:
    """Takes an input formula, reduces it and either alphabetizes or anonymizes it
    following the OPTIMADE standard.

    """

    numbers: list[int] = [int(n.strip() or 1) for n in re.split(r"[A-Z][a-z]*", formula)[1:]]
    # Need to remove leading 1 from split and convert to ints

    species: list[str] = re.findall("[A-Z][a-z]*", formula)

    gcd = reduce(math.gcd, numbers)

    if not len(species) == len(numbers):
        raise ValueError(f"Something is wrong with the input formula: {formula}")

    numbers = [n // gcd for n in numbers]

    if anonymize:
        numbers = sorted(numbers, reverse=True)
        species = [s for _, s in zip(numbers, _optimade_anonymous_element_generator(), strict=False)]

    elif alphabetize:
        species, numbers = zip(*sorted(zip(species, numbers, strict=True)), strict=True)  # type: ignore[assignment]

    return "".join(f"{s}{n if n != 1 else ''}" for n, s in zip(numbers, species, strict=True))


class OptimadeStructureAdapter:
    """Adapter serves as a bridge between OPTIMADE structures and pymatgen objects."""

    @staticmethod
    def get_optimade_structure(structure: Structure | IStructure, **kwargs) -> dict[str, str | dict[str, Any]]:
        """Get a dictionary in the OPTIMADE Structure format from a pymatgen structure or molecule.

        Args:
            structure (Structure): pymatgen Structure
            **kwargs: passed to the ASE Atoms constructor

        Returns:
            A dictionary serialization of the structure in the OPTIMADE format.

        """
        if not structure.is_ordered:
            raise ValueError("OPTIMADE Adapter currently only supports ordered structures")

        attributes: dict[str, Any] = {}
        attributes["cartesian_site_positions"] = structure.lattice.get_cartesian_coords(structure.frac_coords).tolist()
        attributes["lattice_vectors"] = structure.lattice.matrix.tolist()
        attributes["species_at_sites"] = [_.symbol for _ in structure.species]
        attributes["species"] = [
            {"name": _.symbol, "chemical_symbols": [_.symbol], "concentration": [1]}
            for _ in set(structure.composition.elements)
        ]
        attributes["dimension_types"] = [int(_) for _ in structure.lattice.pbc]
        attributes["nperiodic_dimensions"] = sum(attributes["dimension_types"])
        attributes["nelements"] = len(structure.composition.elements)
        attributes["chemical_formula_anonymous"] = _optimade_reduce_or_anonymize_formula(
            structure.composition.formula, anonymize=True
        )
        attributes["elements"] = sorted([_.symbol for _ in structure.composition.elements])
        attributes["chemical_formula_reduced"] = _optimade_reduce_or_anonymize_formula(
            structure.composition.formula, anonymize=False
        )
        attributes["chemical_formula_descriptive"] = structure.composition.formula
        attributes["elements_ratios"] = [structure.composition.get_atomic_fraction(e) for e in attributes["elements"]]
        attributes["nsites"] = len(attributes["species_at_sites"])

        attributes["last_modified"] = None
        attributes["immutable_id"] = None
        attributes["structure_features"] = []

        return {"attributes": attributes}

    @staticmethod
    def get_structure(resource: dict) -> Structure:
        """Get pymatgen structure from an OPTIMADE structure resource.

        Args:
            resource: OPTIMADE structure resource as a dictionary, JSON string, or the
                corresponding attributes dictionary (i.e., `resource["attributes"]`).

        Returns:
            Structure: Equivalent pymatgen Structure

        """
        if isinstance(resource, str):
            try:
                resource = orjson.loads(resource)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Could not decode the input OPTIMADE resource as JSON: {exc}")

        if "attributes" not in resource:
            resource = {"attributes": resource}

        _id = resource.get("id", None)
        attributes = resource["attributes"]
        properties: dict[str, Any] = {"optimade_id": _id}

        # Take any prefixed attributes and save them as properties
        if custom_properties := {k: v for k, v in attributes.items() if k.startswith("_")}:
            properties["optimade_attributes"] = custom_properties

        return Structure(
            lattice=Lattice(
                attributes["lattice_vectors"],
                [bool(d) for d in attributes["dimension_types"]],  # type: ignore[arg-type]
            ),
            species=_pymatgen_species(
                nsites=attributes["nsites"],
                species_at_sites=attributes["species_at_sites"],
            ),
            coords=attributes["cartesian_site_positions"],
            coords_are_cartesian=True,
            properties=properties,
        )
