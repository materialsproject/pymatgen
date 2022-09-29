# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This package contains core modules and classes for representing structures and
operations on them.
"""

import os
import warnings

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
__version__ = "2022.9.21"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    # Load environment variables by default as backup
    d = {}
    for k, v in os.environ.items():
        if k.startswith("PMG_"):
            d[k] = v
        elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
            d["PMG_" + k] = v

    # Override anything in env vars with that in yml file
    if os.path.exists(SETTINGS_FILE):
        try:
            yaml = YAML()
            with open(SETTINGS_FILE) as f:
                d_yml = yaml.load(f)
            d.update(d_yml)
        except Exception as ex:
            # If there are any errors, default to using environment variables
            # if present.
            warnings.warn(f"Error loading .pmgrc.yaml: {ex}. You may need to reconfigure your yaml file.")

    return d


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
