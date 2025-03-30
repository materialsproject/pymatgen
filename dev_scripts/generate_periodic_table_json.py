#!/usr/bin/env python3

"""Create `core.periodic_table.json` from source files.

Each source file may be parsed using a common or custom parser. In cases where
a custom parser is required, it should return either a single `Property` or
a sequence of `Property`.

TODO:
    - use pymatgen Unit
"""

from __future__ import annotations

import csv
import os
import warnings
from dataclasses import dataclass
from io import StringIO
from itertools import product
from typing import TYPE_CHECKING

import pandas as pd
import requests
from ruamel.yaml import YAML

from pymatgen.core import Element

if TYPE_CHECKING:
    from collections.abc import Callable
    from typing import Any, Literal

    from pymatgen.core.units import Unit
    from pymatgen.util.typing import PathLike


ELEMENTS: tuple[str, ...] = tuple([elem.name for elem in Element])
ISOTOPES: tuple[str, ...] = tuple(elem for elem in Element.__members__ if elem not in ELEMENTS)

DEFAULT_VALUE: str = "no data"  # The default value if not provided


@dataclass
class ElemPropertyValue:
    value: Any = DEFAULT_VALUE
    # unit: Unit | None = None  # Don't allow per-value unit for now
    reference: str | None = None  # Parser not implemented


@dataclass
class Property:
    name: str
    data: dict[Element, ElemPropertyValue]
    unit: Unit | None = None
    reference: str | None = None


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
    unit: Unit,
    prop_base: str = "Ionic radii",
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
    print(f"Parsing {prop_base} file: '{file}'")
    print(f"  - Provide property: '{prop_base} (ls/hs)'")

    result: dict[str, dict[Element, ElemPropertyValue]] = {
        prop_base: {},
        f"{prop_base} hs": {},
        f"{prop_base} ls": {},
    }

    with open(file, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
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


def parse_shannon_radii(file: PathLike, unit: Unit) -> Property:
    """Parse Shannon radii from CSV.

    For each element, the ElemPropertyValue has the following structure:
        value[charge][coordination][spin_state] = {
            "crystal_radius": float,
            "ionic_radius": float,
        }

    Empty spin states are stored as empty strings.
    Charges and coordinations are kept as strings (instead of converting to int).
    """
    print(f"Parsing Shannon radii file: '{file}'")
    print("  - Provide property: 'Shannon radii'")

    nested_per_element: dict[Element, dict[str, dict[str, dict[str, dict[str, float]]]]] = {}

    with open(file, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            elem: Element = Element(row["Element"].strip())
            charge: str = row["Charge"].strip()
            coordination: str = row["Coordination"].strip()
            spin: Literal["Low Spin", "High Spin", ""] = row["Spin State"].strip()

            nested_per_element.setdefault(elem, {}).setdefault(charge, {}).setdefault(coordination, {})[spin] = {
                "crystal_radius": float(row["Crystal Radius"]),
                "ionic_radius": float(row["Ionic Radius"]),
            }

    # Flatten into {Element: ElemPropertyValue} mapping
    data: dict[Element, ElemPropertyValue] = {
        elem: ElemPropertyValue(value=nested_data) for elem, nested_data in nested_per_element.items()
    }
    return Property(name="Shannon radii", data=data, unit=unit)


def get_electron_affinities() -> Property:
    """Get electron affinities data from Wikipedia."""
    print("Getting electron affinities from Wikipedia:")
    print("  - Provide property: 'Electron affinity'")

    url: str = "https://en.wikipedia.org/wiki/Electron_affinity_(data_page)"
    tables = pd.read_html(StringIO(requests.get(url, timeout=5).text))

    # Get the "Elements Electron affinity" table (with unit eV)
    ea_df = next(
        table
        for table in tables
        if "Electron affinity (kJ/mol)" in table.columns and table["Name"].astype(str).str.contains("Hydrogen").any()
    )
    ea_df = ea_df.drop(columns=["References", "Electron affinity (kJ/mol)", "Element"])

    # Drop superheavy elements
    max_z: int = max(Element(element).Z for element in Element.__members__)
    ea_df = ea_df[pd.to_numeric(ea_df["Z"], errors="coerce") <= max_z]

    # # Drop heavy isotopes (by detecting duplicate Z)
    # ea_df = ea_df.drop_duplicates(subset="Z", keep="first")

    # Ensure we cover all elements up to Uranium (Z=92)
    if not (z_values := set(ea_df["Z"])).issuperset(range(1, 93)):
        raise ValueError(f"miss electron affinities: {set(range(1, 93)) - z_values}")

    # Clean up electron affinity data
    ea_df["Electron affinity (eV)"] = (
        ea_df["Electron affinity (eV)"]
        .str.replace("−", "-")  # Use the Unicode minus (-)  # noqa: RUF001
        .str.split("(")  # remove uncertainty ("0.754 195(19)")
        .str[0]
        .str.replace(" ", "")  # "0.754 195(19)"
    ).astype(float)

    data = {
        Element.from_name(name.strip()): value
        for name, value in zip(ea_df["Name"], ea_df["Electron affinity (eV)"], strict=True)
    }
    return Property(name="Electron affinity", data=data, unit="eV")


def generate_iupac_ordering() -> list[Property]:
    print("Generating IUPAC ordering:")
    print("  - Provide property: 'iupac_ordering'")  # TODO: duplicate
    print("  - Provide property: 'IUPAC ordering'")

    _order = [
        ([18], range(6, 0, -1)),  # noble gasses
        ([1], range(7, 1, -1)),  # alkali metals
        ([2], range(7, 1, -1)),  # alkali earth metals
        (range(17, 2, -1), [9]),  # actinides
        (range(17, 2, -1), [8]),  # lanthanides
        ([3], (5, 4)),  # Y, Sc
        ([4], (6, 5, 4)),  # Hf -> Ti
        ([5], (6, 5, 4)),  # Ta -> V
        ([6], (6, 5, 4)),  # W -> Cr
        ([7], (6, 5, 4)),  # Re -> Mn
        ([8], (6, 5, 4)),  # Os -> Fe
        ([9], (6, 5, 4)),  # Ir -> Co
        ([10], (6, 5, 4)),  # Pt -> Ni
        ([11], (6, 5, 4)),  # Au -> Cu
        ([12], (6, 5, 4)),  # Hg -> Zn
        ([13], range(6, 1, -1)),  # Tl -> B
        ([14], range(6, 1, -1)),  # Pb -> C
        ([15], range(6, 1, -1)),  # Bi -> N
        ([1], [1]),  # Hydrogen
        ([16], range(6, 1, -1)),  # Po -> O
        ([17], range(6, 1, -1)),  # At -> F
    ]
    order: list[tuple[int, int]] = [item for sublist in (list(product(x, y)) for x, y in _order) for item in sublist]

    iupac_ordering_dict = dict(
        zip(
            [Element.from_row_and_group(row, group) for group, row in order],
            range(len(order)),
            strict=True,
        )
    )

    return [
        Property(name="iupac_ordering", data=iupac_ordering_dict),
        Property(name="IUPAC ordering", data=iupac_ordering_dict),
    ]


def main():
    # Generate and gather Property
    RESOURCES_DIR: str = "periodic_table_resources"

    properties = (
        *parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        *parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
        *parse_yaml(f"{RESOURCES_DIR}/ionization_energies_nist.yaml"),  # Parsed from HTML
        *parse_csv(f"{RESOURCES_DIR}/radii.csv", transform=lambda x: float(x) / 100, unit="nm"),
        *parse_ionic_radii(f"{RESOURCES_DIR}/ionic_radii.csv", unit="nm"),
        parse_shannon_radii(f"{RESOURCES_DIR}/Shannon_Radii.csv", unit="nm"),
        get_electron_affinities(),
        *generate_iupac_ordering(),
    )

    # Check for duplicate
    prop_names: list[str] = [prop.name for prop in properties]
    if len(prop_names) != len(set(prop_names)):
        raise ValueError("Duplicate property name found in Property list.")

    # Save an intermediate YAML copy for manual inspection, otherwise JSON is hard to read
    # TODO: WIP

    # Output to JSON (element->property->value format, and drop metadata)
    # TODO: WIP


if __name__ == "__main__":
    main()
