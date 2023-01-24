# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This package contains core modules and classes for representing structures and
operations on them.
"""

from __future__ import annotations

import os
import warnings
from typing import Any

from ruamel.yaml import YAML

from .composition import Composition as Composition
from .lattice import Lattice as Lattice
from .operations import SymmOp as SymmOp
from .periodic_table import DummySpecies as DummySpecies
from .periodic_table import Element as Element
from .periodic_table import Species as Species
from .sites import PeriodicSite as PeriodicSite
from .sites import Site as Site
from .structure import IMolecule as IMolecule
from .structure import IStructure as IStructure
from .structure import Molecule as Molecule
from .structure import Structure as Structure
from .units import ArrayWithUnit as ArrayWithUnit
from .units import FloatWithUnit as FloatWithUnit
from .units import Unit as Unit

__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ = "shyuep@gmail.com"
__version__ = "2023.1.20"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".config", ".pmgrc.yaml")
OLD_SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings() -> dict[str, Any]:
    settings: dict[str, Any] = {}

    # Load .pmgrc.yaml file
    yaml = YAML()
    for file_path in (SETTINGS_FILE, OLD_SETTINGS_FILE):
        try:
            with open(file_path) as yml_file:
                settings = yaml.load(yml_file) or {}
            break
        except FileNotFoundError:
            continue
        except Exception as exc:
            # If there are any errors, default to using environment variables
            # if present.
            warnings.warn(f"Error loading {file_path}: {exc}.\nYou may need to reconfigure your yaml file.")

    # Override .pmgrc.yaml with env vars (if present)
    for key, val in os.environ.items():
        if key.startswith("PMG_"):
            settings[key] = val
        elif key in ("VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"):
            settings[f"PMG_{key}"] = val

    return settings


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
