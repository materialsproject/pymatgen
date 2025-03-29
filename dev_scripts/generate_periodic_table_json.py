#!/usr/bin/env python3

"""Create `core.periodic_table.json` from source files.

Parser interface requirement:
    In cases where a separate parser is needed for certain property,
    please follow the interface requirement:  # TODO:
        - Should generate Mapping[Element, ElemPropertyValue]

Workflow:
    - Parse basic properties: atomic number/mass and so on  # TODO: revise this
    - TODO: finish workflow

TODO:
    - allow reference from either property or element
    - Would zipped JSON be more efficient (IO bound or not?)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
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


def parse_yamls(working_dir: PathLike) -> dict[PropStr, dict[Element, ElemPropertyValue]]:
    """Parse all YAML files in given directory.

    Expected YAML format:
        TODO:

    Args:
        working_dir (PathLike): directory containing all YAMLs.

    Returns:
        dict[str, dict[Element, ElemPropertyValue]]: A property to
            {Element: ElemPropertyValue} mapping.
    """

    def parse_yaml(file: PathLike) -> dict[PropStr, dict[Element, ElemPropertyValue]]:
        print(f"Parsing YAML file: {yaml_file}")

        with open(file, encoding="utf-8") as f:
            raw = yaml.load(f)

        result: dict[PropStr, dict[Element, ElemPropertyValue]] = {}

        for prop_name, prop_info in raw.items():
            print(f"    - Found property '{prop_name}' ")

            unit = prop_info.get("unit")
            reference = prop_info.get("reference")
            data = prop_info.get("data")

            result[prop_name] = {
                Element(elem): ElemPropertyValue(value=value, unit=unit, reference=reference)
                for elem, value in data.items()
            }

        return result

    yaml = YAML()

    yaml_files = [file for file in Path(working_dir).glob("*.yaml") if not file.name.startswith(".")]

    combined: dict[PropStr, dict[str, ElemPropertyValue]] = {}
    seen_props: set[PropStr] = set()

    for yaml_file in yaml_files:
        prop_data = parse_yaml(yaml_file)

        for prop_name, element_map in prop_data.items():
            if prop_name in seen_props:
                raise ValueError(f"Duplicate property found: '{prop_name}' in file {yaml_file}")

            seen_props.add(prop_name)
            combined[prop_name] = element_map

    return combined


def main():
    parse_yamls(RESOURCES_DIR)


if __name__ == "__main__":
    main()
