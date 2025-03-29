#!/usr/bin/env python3

"""Create `core.periodic_table.json` from source files.

Parser interface requirement:
    In cases where a separate parser is needed for certain properties
    (e.g., when the data must be parsed from a non-YAML format such as HTML, CSV, etc.),
    please follow the required interface:

Each parser function must return:
    dict[PropStr, dict[Element, ElemPropertyValue]]

That is:
- The top-level dictionary maps **property names** (as strings) to their corresponding data.
- For each property, the value is another dictionary mapping:
    - `Element` to `ElemPropertyValue`(includes the actual value and optional unit)

This ensures that all parsers, regardless of data source, return a consistent format that
can be merged into the overall dataset using `generate_json`.

TODO:
    - would zipped JSON be more efficient (IO bound or not?)
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd
from ruamel.yaml import YAML

from pymatgen.core import Element

if TYPE_CHECKING:
    from typing import Any, TypeAlias

    from pymatgen.core.units import Unit
    from pymatgen.util.typing import PathLike


ELEMENTS: tuple[str, ...] = tuple([elem.name for elem in Element])
ISOTOPES: tuple[str, ...] = tuple(elem for elem in Element.__members__ if elem not in ELEMENTS)

DEFAULT_VALUE: str = "no data"  # The default value if not provided


@dataclass
class ElemPropertyValue:
    value: Any = DEFAULT_VALUE
    unit: Unit | None = None


PropStr: TypeAlias = str
Sources: TypeAlias = dict[PropStr, dict[Element, ElemPropertyValue]]


def parse_yaml(file: PathLike) -> Sources:
    """Parse a YAML file.

    Expected YAML format:
        We expect each YAML file to contain one or more properties.
        Each property should follow this structure:
            - `unit` (optional): The unit of measurement for the values.
            - `data`: Dict mapping each element symbol (e.g., "Fe") to its corresponding value.

    Args:
        working_dir (PathLike): directory containing all YAMLs.

    Returns:
        Sources: A property to {Element: ElemPropertyValue} mapping.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"YAML file {file} does not exist")

    print(f"Parsing YAML file: '{file}'")

    yaml = YAML()
    with open(file, encoding="utf-8") as f:
        raw = yaml.load(f)

    result: Sources = {}

    for prop_name, prop_info in raw.items():
        print(f"  - Found property: '{prop_name}' ")

        # TODO: convert to pymatgen Unit
        unit = prop_info.get("unit")
        data = prop_info.get("data")

        result[prop_name] = {Element(elem): ElemPropertyValue(value=value, unit=unit) for elem, value in data.items()}

    return result


def parse_csv(
    file: PathLike,
    dtype: Any | None = None,
) -> Sources:
    """Parse a CSV file.

    Expected CSV format:
        We expected each CSV file to contain one or more properties,
        where each property occupies a single column and the column name
        would be used as the property name.
        The index (first) column should be the elements.  TODO: clarify

    Args:
        file (PathLike): The CSV file to parse.

    Returns:
        Sources: A property to {Element: ElemPropertyValue} mapping.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"CSV file {file} does not exist")

    print(f"Parsing CSV file: '{file}'")

    # Expect the first column to be the element symbols
    data_df = pd.read_csv(file)
    data_df = data_df.set_index(data_df.columns[0])

    result: Sources = {}

    for prop in data_df.columns:
        print(f"  - Found property: '{prop}' ")

        prop_values = {}
        for symbol, value in data_df[prop].items():
            if pd.isna(value):
                value = DEFAULT_VALUE
            elif dtype is not None:
                try:
                    value = dtype(value)
                except (ValueError, TypeError):
                    value = str(value)

            prop_values[Element(symbol)] = ElemPropertyValue(value=value)
        result[prop] = prop_values

    return result


def generate_yaml_and_json(*sources: Sources) -> None:
    """Generate the intermediate YAML and final JSON from sources."""
    # Check for duplicate and combine all sources
    combined: Sources = {}
    seen_props: set[PropStr] = set()

    for source in sources:
        for prop_name, element_map in source.items():
            if prop_name in seen_props:
                raise ValueError(f"Duplicate property found: '{prop_name}' in source")
            seen_props.add(prop_name)
            combined[prop_name] = element_map

    # Save an intermediate YAML copy for manual inspection,
    # otherwise JSON is hard to read
    # TODO: WIP

    # Output to JSON (element->property->value format, and drop metadata)
    # TODO: WIP


def main():
    RESOURCES_DIR: str = "periodic_table_resources"

    generate_yaml_and_json(
        parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/ionization_energies.yaml"),
        parse_csv(f"{RESOURCES_DIR}/radii.csv", dtype=lambda x: float(x) / 100),
    )


if __name__ == "__main__":
    main()
