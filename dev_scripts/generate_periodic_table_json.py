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
    - use pymatgen Unit
    - convert Shannon Radii to CSV (better version control and no `openpyxl` needed)
    - gen_iupac_ordering
    - add_electron_affinities
"""

from __future__ import annotations

import csv
import os
import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import pandas as pd
from ruamel.yaml import YAML

from pymatgen.core import Element

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any

    from pymatgen.core.units import Unit
    from pymatgen.util.typing import PathLike


ELEMENTS: tuple[str, ...] = tuple([elem.name for elem in Element])
ISOTOPES: tuple[str, ...] = tuple(elem for elem in Element.__members__ if elem not in ELEMENTS)

DEFAULT_VALUE: str = "no data"  # The default value if not provided


@dataclass
class ElemPropertyValue:
    value: Any = DEFAULT_VALUE
    # unit: Unit | None = None  # Don't allow per-value unit for now
    # reference: str | None = None


@dataclass
class Property:
    name: str
    unit: Unit | None = None
    reference: str | None = None
    data: dict[Element, ElemPropertyValue] = field(default_factory=dict)


def parse_yaml(file: PathLike) -> list[Property]:
    """Parse a YAML file.

    Expected YAML format:
        We expect each YAML file to contain one or more properties.
        Each property should follow this structure:
            - `unit` (optional): The unit of measurement for the values.
            - `data`: Dict mapping each element symbol (e.g., "Fe") to its corresponding value.

    Args:
        working_dir (PathLike): directory containing all YAMLs.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"YAML file {file} does not exist")

    print(f"Parsing YAML file: '{file}'")

    yaml = YAML()
    with open(file, encoding="utf-8") as f:
        raw = yaml.load(f)

    result: list[Property] = []

    for prop_name, prop_info in raw.items():
        print(f"  - Found property: '{prop_name}' ")

        data_block = prop_info.get("data")
        data = {Element(elem): ElemPropertyValue(value=val) for elem, val in data_block.items()}

        result.append(
            Property(
                name=prop_name,
                unit=prop_info.get("unit"),
                reference=prop_info.get("reference"),
                data=data,
            )
        )

    return result


def parse_csv(
    file: PathLike,
    transform: Callable | None = None,
    unit: Unit | None = None,
    reference: str | None = None,
) -> list[Property]:
    """Parse a CSV file.

    Expected CSV format:
        - The CSV must contain a column named "element" (case-insensitive),
          which will be used as the index for mapping values to elements.
        - All other columns are treated as distinct property names, with each
          column header used as the property name (e.g., "Atomic mass").
        - Each row corresponds to a single element, and each cell represents
          the value of a property for that element.

    Args:
        file (PathLike): The CSV file to parse.
        transform (Callable): Optional function to convert each value.
            If provided, it will be applied to each non-null cell value.
        unit (Unit): Unit passed to Property.
        reference (str): Reference passed to Property.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"CSV file {file} does not exist")

    print(f"Parsing CSV file: '{file}'")

    # Try to locate the column named "element" (case-insensitive)
    data_df = pd.read_csv(file)
    try:
        index_col = data_df.columns[data_df.columns.str.lower().str.strip() == "element"][0]
    except IndexError as exc:
        raise ValueError(f"Could not find an 'element' column in CSV: {file}") from exc

    data_df = data_df.set_index(index_col)

    result: list[Property] = []

    for prop in data_df.columns:
        print(f"  - Found property: '{prop}' ")

        data: dict[Element, ElemPropertyValue] = {}

        for symbol, value in data_df[prop].items():
            if pd.isna(value):
                value = DEFAULT_VALUE
            elif transform is not None:
                try:
                    value = transform(value)
                except (ValueError, TypeError):
                    warnings.warn(f"Cannot transform {value=}, keep as string", stacklevel=2)
                    value = str(value)

            data[Element(symbol)] = ElemPropertyValue(value=value)

        result.append(Property(name=prop, unit=unit, reference=reference, data=data))

    return result


def parse_ionic_radii(
    file: PathLike,
    prop_base: str = "Ionic radii",
    unit: Unit = "nm",
) -> list[Property]:
    """Parse ionic radii from CSV.

    CSV Format:
        - Must include columns: "Element", "Spin", and oxidation states (-3 to +8).
        - "Spin" can be empty, "hs" (high spin) or "ls" (low spin).

    Behavior:
        - For each row:
            - If "Spin" is empty, the element is added under "Ionic radii".
            - If "Spin" is "hs" or "ls", the element is added under "Ionic radii hs"/"ls".
            - High spin ("hs") data is also copied into the base "Ionic radii" property.
        - Radii values are divided by 100 (converted from pm to nm).
    """
    print(f"Parsing {prop_base} CSV file: '{file}'")

    result: dict[str, dict[Element, ElemPropertyValue]] = {
        prop_base: {},
        f"{prop_base} hs": {},
        f"{prop_base} ls": {},
    }

    with open(file, encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            elem = Element(row["Element"].strip())
            spin = row.get("Spin", "").strip().lower()
            prop_name = f"{prop_base} {spin}".strip()

            # Collect non-empty fields
            ox_state_data = {
                ox: float(val) / 100  # NOTE: radii divided by 100
                for ox, val in row.items()
                if ox not in {"Element", "Spin"} and val.strip() != ""
            }

            result[prop_name][elem] = ElemPropertyValue(value=ox_state_data)
            # Copy high-spin radii to the base "Ionic radii"
            if spin == "hs":
                result[prop_base][elem] = ElemPropertyValue(value=ox_state_data)

    return [Property(name=name, unit=unit, data=data) for name, data in result.items()]


def parse_shannon_radii(file: PathLike):
    pass


def generate_yaml_and_json(*sources: list[Property]) -> None:
    """Generate the intermediate YAML and final JSON from Properties."""
    # Flatten all Property
    all_properties: list[Property] = [prop for group in sources for prop in group]

    # Check for duplicate
    prop_names: list[str] = [prop.name for prop in all_properties]
    if len(prop_names) != len(set(prop_names)):
        raise ValueError("Duplicate property name found in Property list.")

    # Save an intermediate YAML copy for manual inspection, otherwise JSON is hard to read
    # TODO: WIP

    # Output to JSON (element->property->value format, and drop metadata)
    # TODO: WIP


def main():
    RESOURCES_DIR: str = "periodic_table_resources"

    generate_yaml_and_json(
        parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/ionization_energies_nist.yaml"),  # Parsed from HTML
        parse_csv(f"{RESOURCES_DIR}/radii.csv", transform=lambda x: float(x) / 100, unit="nm"),
        parse_ionic_radii(f"{RESOURCES_DIR}/ionic_radii.csv"),
        # parse_shannon_radii(f"{RESOURCES_DIR}/Shannon_Radii.xlsx"),
    )


if __name__ == "__main__":
    main()
