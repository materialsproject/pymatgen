#!/usr/bin/env python3

"""Create `core.periodic_table.json` from source files.

Parser interface requirement:
    In cases where a separate parser is needed for certain property,
    please follow the interface requirement:  # TODO:
        - Should generate Mapping[Element, ElemPropertyValue]

TODO:
    - allow reference from either property or element
    - allow global default + per-property default
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

RESOURCES_DIR: str = "periodic_table_resources"

PropStr: TypeAlias = str


@dataclass
class ElemPropertyValue:
    value: Any = "no data"  # TODO: better missing value handling?
    unit: Unit | None = None
    reference: str | None = None


def parse_yaml(file: PathLike) -> dict[PropStr, dict[Element, ElemPropertyValue]]:
    """Parse a YAML file.

    Expected YAML format:
        TODO:

    TODO:
        - Allow a "per property" default value aside from the global
            default as "no data", for example "Is named isotope"
            should default to false instead of "no data"

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
        print(f"    - Found property: '{prop_name}' ")

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
    generate_json(
        # parse_yaml(f"{RESOURCES_DIR}/elemental_properties.yaml"),
        parse_yaml(f"{RESOURCES_DIR}/oxidation_states.yaml"),
    )


if __name__ == "__main__":
    main()
