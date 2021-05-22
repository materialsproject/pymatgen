# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This package contains core modules and classes for representing structures and
operations on them.
"""

import os

try:
    import ruamel.yaml as yaml
except ImportError:
    try:
        import ruamel_yaml as yaml  # type: ignore  # noqa
    except ImportError:
        import yaml  # type: ignore # noqa

from .composition import Composition  # noqa
from .lattice import Lattice  # noqa
from .operations import SymmOp  # noqa
from .periodic_table import DummySpecies, Element, Species  # noqa
from .sites import PeriodicSite, Site  # noqa
from .structure import IMolecule, IStructure, Molecule, Structure  # noqa
from .units import ArrayWithUnit, FloatWithUnit, Unit  # noqa


__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ = "shyuep@gmail.com"
__version__ = "2022.0.8"


SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".pmgrc.yaml")


def _load_pmg_settings():
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        for k, v in os.environ.items():
            if k.startswith("PMG_"):
                d[k] = v
            elif k in ["VASP_PSP_DIR", "MAPI_KEY", "DEFAULT_FUNCTIONAL"]:
                d["PMG_" + k] = v
    d = d or {}
    return dict(d)


SETTINGS = _load_pmg_settings()
locals().update(SETTINGS)
