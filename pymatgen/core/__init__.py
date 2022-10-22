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

from .composition import Composition
from .lattice import Lattice
from .operations import SymmOp
from .periodic_table import DummySpecies, Element, Species
from .sites import PeriodicSite, Site
from .structure import IMolecule, IStructure, Molecule, Structure
from .units import ArrayWithUnit, FloatWithUnit, Unit

__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ = "shyuep@gmail.com"
__version__ = "2022.10.22"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".config", ".pmgrc.yaml")
OLD_SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings() -> dict[str, Any]:
    settings: dict[str, Any] = {}

    # Load .pmgrc.yaml file
    yaml = YAML()
    try:
        with open(SETTINGS_FILE) as yml_file:
            settings = yaml.load(yml_file) or {}
    except FileNotFoundError:
        try:
            with open(OLD_SETTINGS_FILE) as yml_file:
                settings = yaml.load(yml_file)
        except FileNotFoundError:
            pass
    except Exception as ex:
        # If there are any errors, default to using environment variables
        # if present.
        warnings.warn(f"Error loading .pmgrc.yaml: {ex}. You may need to reconfigure your yaml file.")

    # Override .pmgrc.yaml with env vars if present
    for key, val in os.environ.items():
        if key.startswith("PMG_"):
            settings[key] = val
        elif key in ("VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"):
            settings[f"PMG_{key}"] = val

    return settings


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
