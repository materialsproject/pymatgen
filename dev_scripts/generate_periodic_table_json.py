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
    - Would zipped JSON be more efficient (IO bound or not?)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pymatgen.core import Element

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core.units import Unit


ELEMENTS: tuple[str, ...] = tuple([elem.name for elem in Element])
ISOTOPES: tuple[str, ...] = tuple(elem for elem in Element.__members__ if elem not in ELEMENTS)

RESOURCES_DIR: str = "periodic_table_resources"


@dataclass
class ElemPropertyValue:
    value: Any = "no data"  # TODO: better missing value handling?
    unit: Unit | None = None
    reference: str | None = None
    comment: str | None = None


def main():
    parse_oxi_state()


if __name__ == "__main__":
    main()
