"""This package contains core modules and classes for representing structures and operations on them."""

from __future__ import annotations

import os
import warnings
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from ruamel.yaml import YAML

from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import DummySpecie, DummySpecies, Element, Species, get_el_sp
from pymatgen.core.sites import PeriodicSite, Site
from pymatgen.core.structure import IMolecule, IStructure, Molecule, PeriodicNeighbor, SiteCollection, Structure
from pymatgen.core.units import ArrayWithUnit, FloatWithUnit, Unit

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


def _expand_strings(data):
    if isinstance(data, dict):
        return {k: _expand_strings(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [_expand_strings(i) for i in data]
    elif isinstance(data, str):
        return os.path.expandvars(os.path.expanduser(data))
    else:
        return data

def _load_yaml_file(file_path: str) -> Dict[str, Any]:
    yaml = YAML()
    try:
        expanded_path = os.path.expanduser(file_path)
        with open(expanded_path, encoding="utf-8") as yml_file:
            return _expand_strings(yaml.load(yml_file)) or {}
    except FileNotFoundError:
        warnings.warn(f"File not found: {file_path}")
        return {}
    except Exception as exc:
        warnings.warn(f"Error loading {file_path}: {exc}.\nYou may need to reconfigure your YAML file.")
        return {}

def _load_pmg_settings(file_path: str = None) -> Dict[str, Any]:
    settings: Dict[str, Any] = {}

    # Load settings from the specified file or default files
    if file_path:
        settings.update(_load_yaml_file(file_path))
    else:
        for fp in (SETTINGS_FILE, OLD_SETTINGS_FILE):
            settings.update(_load_yaml_file(fp))
            if settings:
                break

    # Override with environment variables (if present)
    for key, val in os.environ.items():
        if key.startswith("PMG_"):
            settings[key] = val
        elif key in ("VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"):
            settings[f"PMG_{key}"] = val

    return settings
    

SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
