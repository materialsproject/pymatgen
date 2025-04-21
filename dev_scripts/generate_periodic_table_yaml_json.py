#!/usr/bin/env python3

"""Create `core.periodic_table.json` from source files, and as such
you should NOT modify the JSON directly, but work on the data source
and then run this script to regenerate the JSON/YAML.

Each source file may be parsed using a common or custom parser. In cases where
a custom parser is required, it should return either a single `Property` or
a sequence of `Property`.

The YAML file is a readable aggregation of all properties in the following structure:
    Property name -> {
        unit: <unit string or null>,
        reference: <reference string or null>,
        data: {
            <element as str>: <value>
        }
    }

The JSON file is a compact, production-format structure (no metadata):
    <element symbol> -> {
        <property name>: <value unit>
    }

    Units are stored separately in a special top-level key:
            "_unit" -> {
                <property name>: <unit string>
            }

TODO: this script is in a transitional state, it would now:
    - Generate a duplicate JSON instead of overwriting the one in `core` dir
    - Append unit as string instead of as a separate field
    - Unit should be converted to pymatgen Unit
"""

from __future__ import annotations

import csv
import json
import os
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from io import StringIO
from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import requests
from ruamel.yaml import YAML

from pymatgen.core import Element

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from typing import Any, Literal

    from pymatgen.util.typing import PathLike


DEFAULT_VALUE: str = "no data"  # The default value if not provided

RESOURCES_DIR: str = f"{Path(__file__).parent}/periodic_table_resources"


@dataclass
class ElemPropertyValue:
    value: Any = DEFAULT_VALUE
    reference: str | None = None  # Parser not implemented


@dataclass
class Property:
    name: str
    data: dict[Element, ElemPropertyValue]
    unit: str | None = None
    reference: str | None = None

    # TODO: don't do unit conversion for now
    # def __post_init__(self):
    #     if self.unit is not None and isinstance(self.unit, str):
    #         self.unit = Unit(self.unit)


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

        if "data" not in prop_info:
            raise ValueError(f"Missing 'data' block for property '{prop_name}' in YAML {file}")

        data = {Element(elem): ElemPropertyValue(value=val) for elem, val in prop_info["data"].items()}

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
    *,
    transform: Callable | None = None,
    unit: str | None = None,
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
        unit (str): Unit passed to Property.
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
    *,
    unit: str,
    reference: str,
    prop_base: str = "Ionic radii",
) -> list[Property]:
    """Parse ionic radii from CSV.

    Behavior:
        - For each row:
            - If "Spin" is empty, the element is added under "Ionic radii".
            - If "Spin" is "hs" or "ls", the element is added under "Ionic radii hs"/"ls".
            - High spin ("hs") data is also copied into the base "Ionic radii" property.
        - Radii values are divided by 100 (converted from pm to nm).
    """
    print(f"Parsing {prop_base} file: '{file}'")
    print(f"  - Provide properties: '{prop_base} (ls/hs)'")

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
                ox: float(val) for ox, val in row.items() if ox not in {"Element", "Spin"} and val.strip() != ""
            }

            result[prop_name][elem] = ElemPropertyValue(value=ox_state_data)
            # Copy high-spin radii to the base "Ionic radii"
            if spin == "hs":
                result[prop_base][elem] = ElemPropertyValue(value=ox_state_data)

    return [Property(name=name, unit=unit, data=data, reference=reference) for name, data in result.items()]


def parse_shannon_radii(
    file: PathLike,
    *,
    unit: str,
    reference: str | None = None,
) -> Property:
    """Parse Shannon radii from CSV.

    For each element, the ElemPropertyValue has the following structure:
        value[charge][coordination][spin_state] = {
            "crystal_radius": float,
            "ionic_radius": float,
        }

    Empty spin states are stored as empty strings.
    Charges and coordinations are kept as strings (instead of converting to int).

    TODO: data source/reference is unknown
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
    return Property(name="Shannon radii", data=data, unit=unit, reference=reference)


def get_and_parse_electronic_affinities(prop_name: str = "Electron affinity") -> Property:
    """Get electronic affinities from Wikipedia and save a local YAML copy.

    TODO:
        current wikipedia crawler drops data for Deuterium.
    """
    # Get data table from Wikipedia
    url: str = "https://en.wikipedia.org/wiki/Electron_affinity_(data_page)"
    tables = pd.read_html(StringIO(requests.get(url, timeout=5).text))

    # Get the "Elements Electron affinity" table (with unit eV)
    ea_df = next(
        table
        for table in tables
        if f"{prop_name} (kJ/mol)" in table.columns and table["Name"].astype(str).str.contains("Hydrogen").any()
    )
    ea_df = ea_df.drop(columns=["References", f"{prop_name} (kJ/mol)", "Element"])

    # Drop superheavy elements
    max_z: int = max(Element(element).Z for element in Element.__members__)
    ea_df = ea_df[pd.to_numeric(ea_df["Z"], errors="coerce") <= max_z]

    # Drop heavy isotopes
    # TODO: correctly handle Deuterium
    ea_df = ea_df.drop_duplicates(subset="Z", keep="first")

    # Ensure we cover all elements up to Uranium (Z=92)
    if not (z_values := set(ea_df["Z"])).issuperset(range(1, 93)):
        raise ValueError(f"miss electron affinities: {set(range(1, 93)) - z_values}")

    # Clean up electron affinity data
    ea_df[f"{prop_name} (eV)"] = (
        ea_df[f"{prop_name} (eV)"]
        .str.replace("−", "-")  # Use the Unicode minus (-)  # noqa: RUF001
        .str.replace("(", "")  # "0.754 195(19)"
        .str.replace(")", "")
        .str.replace(" ", "")
        .str.replace("†", "")  # Remove dagger mark
    ).astype(float)

    # Save a local YAML copy
    yaml = YAML()
    yaml.default_flow_style = False

    data: dict[Element, float] = {
        Element.from_name(name.strip()): value
        for name, value in zip(ea_df["Name"], ea_df[f"{prop_name} (eV)"], strict=True)
    }

    data[Element.D] = 0.754674  # TODO: fix this

    output_data = {
        prop_name: {
            "unit": "eV",
            "reference": url,
            "data": {el.name: val for el, val in sorted(data.items(), key=lambda x: x[0].Z)},
        }
    }

    with open(f"{RESOURCES_DIR}/_electron_affinities.yaml", "w", encoding="utf-8") as f:
        yaml.dump(output_data, f)

    return parse_yaml(f"{RESOURCES_DIR}/_electron_affinitIES.yaml")[0]


def generate_iupac_ordering() -> Property:
    """Ordering according to Table VI of "Nomenclature of Inorganic Chemistry
    (IUPAC Recommendations 2005)". This ordering effectively follows the
    groups and rows of the periodic table, except the Lanthanides, Actinides
    and hydrogen.
    """
    print("Generating IUPAC ordering:")
    print("  - Provide property: 'IUPAC ordering'")

    _orders = [
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
    orders: list[tuple[int, int]] = [item for sublist in (list(product(x, y)) for x, y in _orders) for item in sublist]

    iupac_ordering_dict: dict[Element, ElemPropertyValue] = {
        Element.from_row_and_group(row, group): ElemPropertyValue(value=index)
        for index, (group, row) in enumerate(orders)
    }

    return Property(name="IUPAC ordering", data=iupac_ordering_dict)


def generate_yaml_and_json(
    properties: Iterable[Property],
    *,
    yaml_file: PathLike,
    json_file: PathLike,
) -> None:
    """
    Generate a human-readable YAML and a production-ready JSON file from properties.

    Args:
        properties (Iterable[Property]): A sequence of Property objects to serialize.
        yaml_file (PathLike): Path to output YAML file (for development).
        json_file (PathLike): Path to output JSON file (for production).

    Raises:
        ValueError: If duplicate property names are found in the input.
    """
    # Check for duplicate
    counter: Counter = Counter([prop.name for prop in properties])
    if duplicates := [name for name, count in counter.items() if count > 1]:
        raise ValueError(f"Duplicate property found: {', '.join(duplicates)}")

    # Sort properties by name alphabetically
    properties = sorted(properties, key=lambda prop: prop.name)

    # Save a YAML copy for development
    yaml_data: dict[str, dict[Literal["unit", "reference", "data"], Any]] = {}
    for prop in properties:
        # Sort elements by atomic number (Z)
        sorted_data: dict[str, Any] = dict(
            sorted(((elem.name, val.value) for elem, val in prop.data.items()), key=lambda pair: Element(pair[0]).Z)
        )

        yaml_data[prop.name] = {
            "unit": str(prop.unit) if prop.unit is not None else None,
            "reference": prop.reference,
            "data": sorted_data,
        }

    yaml = YAML()
    yaml.default_flow_style = None
    with open(yaml_file, "w", encoding="utf-8") as f:
        yaml.dump(yaml_data, f)

    print("=" * 50)
    print(f"Saved YAML to: {yaml_file}")

    # Output to JSON (element -> property -> value format, and drop metadata)
    element_to_props: dict[str, dict[str, Any]] = defaultdict(dict)

    # # Insert units under a special `_unit` key
    # element_to_props["_unit"] = {}

    for prop in properties:
        # # Store unit for this property if available
        # if prop.unit is not None:
        #     element_to_props["_unit"][prop.name] = str(prop.unit)

        for elem, prop_val in prop.data.items():
            unit = prop.unit
            if unit is None:
                element_to_props[elem.name][prop.name] = prop_val.value
            else:
                output = str(prop_val.value) + " " + str(unit)
                element_to_props[elem.name][prop.name] = output

    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(element_to_props, f, indent=4)

    print(f"Saved JSON to: {json_file}")


def main():
    # Generate and gather Property
    properties: tuple[Property, ...] = (
        *parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        *parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
        *parse_yaml(f"{RESOURCES_DIR}/ionization_energies_nist.yaml"),  # Parsed from HTML
        get_and_parse_electronic_affinities(),
        *parse_csv(f"{RESOURCES_DIR}/radii.csv", transform=float, unit=None),
        *parse_ionic_radii(
            f"{RESOURCES_DIR}/ionic_radii.csv", unit=None, reference="https://en.wikipedia.org/wiki/Ionic_radius"
        ),
        parse_shannon_radii(f"{RESOURCES_DIR}/Shannon_Radii.csv", unit=None),
        generate_iupac_ordering(),
    )

    generate_yaml_and_json(
        properties,
        yaml_file=f"{RESOURCES_DIR}/_periodic_table.yaml",
        json_file=f"{RESOURCES_DIR}/_periodic_table.json",
        # json_file=f"{PKG_DIR}/core/periodic_table.json",
    )


if __name__ == "__main__":
    main()
