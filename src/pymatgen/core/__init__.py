"""This package contains core modules and classes for representing structures and operations on them."""

from __future__ import annotations

import os
import warnings
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import DummySpecie, DummySpecies, Element, Species, get_el_sp
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.structure import IMolecule, IStructure, Molecule, PeriodicNeighbor, SiteCollection, Structure
from pymatgen.core.units import ArrayWithUnit, FloatWithUnit, Unit
from ruamel.yaml import YAML

__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong, Matthew Horton, Janosh Riebesell"
__maintainer_email__ = "shyuep@gmail.com"
try:
    __version__ = version("pymatgen")
except PackageNotFoundError:  # pragma: no cover
    # package is not installed
    pass


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".config", ".pmgrc.yaml")
OLD_SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PKG_DIR = os.path.dirname(MODULE_DIR)
ROOT = os.path.dirname(PKG_DIR)


def _load_pmg_settings() -> dict[str, Any]:
    settings: dict[str, Any] = {}

    # PMG_CONFIG_FILE takes precedence over default settings location
    settings_file = os.getenv("PMG_CONFIG_FILE") or SETTINGS_FILE

    # Load .pmgrc.yaml file
    yaml = YAML()
    for file_path in (settings_file, OLD_SETTINGS_FILE):
        try:
            with open(file_path, encoding="utf-8") as yml_file:
                settings = yaml.load(yml_file) or {}
            break
        except FileNotFoundError:
            continue
        except Exception as exc:
            # If there are any errors, default to using environment variables
            # if present.
            warnings.warn(f"Error loading {file_path}: {exc}.\nYou may need to reconfigure your YAML file.")

    # Override .pmgrc.yaml with env vars (if present)
    for key, val in os.environ.items():
        if key.startswith("PMG_"):
            settings[key] = val
        elif key in ("VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"):
            settings[f"PMG_{key}"] = val

    for key, val in settings.items():
        if key.endswith("_DIR"):
            # Enables the use of $HOME and ~ in paths.
            val = os.path.expanduser(val)
            val = os.path.expandvars(val)
            settings[key] = val

    return settings


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
