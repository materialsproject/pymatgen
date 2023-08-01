# pylint: disable=C0414,W0718,C0301
# ruff: noqa: PLC0414
"""This package contains core modules and classes for representing structures and operations on them."""

from __future__ import annotations

import os
import warnings
from typing import Any

from ruamel.yaml import YAML

from pymatgen.core.composition import Composition as Composition
from pymatgen.core.lattice import Lattice as Lattice
from pymatgen.core.operations import SymmOp as SymmOp
from pymatgen.core.periodic_table import DummySpecie as DummySpecie
from pymatgen.core.periodic_table import DummySpecies as DummySpecies
from pymatgen.core.periodic_table import Element as Element
from pymatgen.core.periodic_table import Species as Species
from pymatgen.core.sites import PeriodicSite as PeriodicSite
from pymatgen.core.sites import Site as Site
from pymatgen.core.structure import IMolecule as IMolecule
from pymatgen.core.structure import IStructure as IStructure
from pymatgen.core.structure import Molecule as Molecule
from pymatgen.core.structure import SiteCollection as SiteCollection
from pymatgen.core.structure import Structure as Structure
from pymatgen.core.units import ArrayWithUnit as ArrayWithUnit
from pymatgen.core.units import FloatWithUnit as FloatWithUnit
from pymatgen.core.units import Unit as Unit

__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong, Matthew Horton, Janosh Riebesell"
__maintainer_email__ = "shyuep@gmail.com"
__version__ = "2023.7.17"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".config", ".pmgrc.yaml")
OLD_SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PKG_DIR = os.path.dirname(MODULE_DIR)
ROOT = os.path.dirname(PKG_DIR)


def _load_pmg_settings() -> dict[str, Any]:
    settings: dict[str, Any] = {}

    # Load .pmgrc.yaml file
    yaml = YAML()
    for file_path in (SETTINGS_FILE, OLD_SETTINGS_FILE):
        try:
            with open(file_path, encoding="utf-8") as yml_file:
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
