# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
# pylint: disable=C0413

"""
Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. This is the root package.
"""
__author__ = "Pymatgen Development Team"
__email__ = "pymatgen@googlegroups.com"
__maintainer__ = "Shyue Ping Ong"
__maintainer_email__ = "shyuep@gmail.com"
__version__ = "2021.3.5"

# Useful aliases for commonly used objects and modules.
# Allows from pymatgen import <class> for quick usage.
# Note that these have to come after the SETTINGS have been loaded. Otherwise, import does not work.

from .core.composition import Composition  # noqa
from .core.lattice import Lattice  # noqa
from .core.operations import SymmOp  # noqa
from .core.periodic_table import DummySpecie, DummySpecies, Element, Specie, Species  # noqa
from .core.sites import PeriodicSite, Site  # noqa
from .core.structure import IMolecule, IStructure, Molecule, Structure  # noqa
from .core.units import ArrayWithUnit, FloatWithUnit, Unit  # noqa
from .electronic_structure.core import Orbital, Spin  # noqa
from .ext.matproj import MPRester  # noqa
