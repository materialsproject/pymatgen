"""
This package contains core modules and classes for representing structures and
operations on them.
"""

__author__ = "Shyue Ping Ong"
__date__ = "Dec 15, 2010 7:21:29 PM"

from .periodic_table import Element, Specie, DummySpecie, \
    smart_element_or_specie
from .composition import Composition
from .structure import Structure, IStructure, Molecule, IMolecule
from .bonds import CovalentBond, get_bond_length
from .lattice import Lattice
from .sites import Site, PeriodicSite
from .operations import SymmOp
from .units import *