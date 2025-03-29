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
    - `Element` to `ElemPropertyValue`(includes the actual value and optional unit/reference)

This ensures that all parsers, regardless of data source, return a consistent format that
can be merged into the overall dataset using `generate_json`.

TODO:
    - allow reference from either property or element
    - would zipped JSON be more efficient (IO bound or not?)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ruamel.yaml import YAML

from pymatgen.core import Element

if TYPE_CHECKING:
    from typing import Any, TypeAlias

    from pymatgen.core.units import Unit
    from pymatgen.util.typing import PathLike


ELEMENTS: tuple[str, ...] = tuple([elem.name for elem in Element])
ISOTOPES: tuple[str, ...] = tuple(elem for elem in Element.__members__ if elem not in ELEMENTS)

# The global default value if not provided
DEFAULT_VALUE: str = "no data"

RESOURCES_DIR: str = "periodic_table_resources"

PropStr: TypeAlias = str


@dataclass
class ElemPropertyValue:
    value: Any = DEFAULT_VALUE
    unit: Unit | None = None
    reference: str | None = None


def parse_yaml(file: PathLike) -> dict[PropStr, dict[Element, ElemPropertyValue]]:
    """Parse a YAML file.

    Expected YAML format:
        We expect each YAML file to contain one or more properties.
        Each property should follow this structure:
            - `unit` (optional): The unit of measurement for the values.
            - `reference` (optional): The reference or source from which the data is derived.
            - `data`: Dict mapping each element symbol (e.g., "Fe") to its corresponding value.

    Args:
        working_dir (PathLike): directory containing all YAMLs.

    Returns:
        dict[str, dict[Element, ElemPropertyValue]]: A property to
            {Element: ElemPropertyValue} mapping.
    """
    print(f"Parsing YAML file: '{file}'")

    yaml = YAML()
    with open(file, encoding="utf-8") as f:
        raw = yaml.load(f)

    result: dict[PropStr, dict[Element, ElemPropertyValue]] = {}

    for prop_name, prop_info in raw.items():
        print(f"  - Found property: '{prop_name}' ")

        # TODO: convert unit to pymatgen Unit?
        unit = prop_info.get("unit")
        reference = prop_info.get("reference")
        data = prop_info.get("data")

        result[prop_name] = {
            Element(elem): ElemPropertyValue(value=value, unit=unit, reference=reference)
            for elem, value in data.items()
        }

    return result


def generate_json(*sources: dict[PropStr, dict[Element, ElemPropertyValue]]) -> None:
    """Generate the final JSON from sources, each source may contain multiple properties."""
    # Check for duplicate and combine all sources
    combined: dict[PropStr, dict[Element, ElemPropertyValue]] = {}
    seen_props: set[PropStr] = set()

    for source in sources:
        for prop_name, element_map in source.items():
            if prop_name in seen_props:
                raise ValueError(f"Duplicate property found: '{prop_name}' in source")
            seen_props.add(prop_name)
            combined[prop_name] = element_map

    # Output to JSON (element->property->value format, and drop metadata)
    # TODO: WIP


def main():
    # TODO: save an intermediate YAML copy for manual inspection,
    # otherwise JSON is hard to read
    generate_json(
        parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/ionization_energies.yaml"),
    )


if __name__ == "__main__":
    main()
