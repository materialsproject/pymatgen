from __future__ import annotations

import unittest

import pytest
from pytest import approx

from pymatgen.core import Element, Site
from pymatgen.core.bonds import CovalentBond, get_bond_length, get_bond_order, obtain_all_bond_lengths

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Jul 26, 2012"


class TestCovalentBond(unittest.TestCase):
    def test_length(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        assert approx(CovalentBond(site1, site2).length - 0.92195444572928864) == 0

    def test_get_bond_order(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1.08])
        assert approx(CovalentBond(site1, site2).get_bond_order() - 1) == 0
        bond = CovalentBond(Site("C", [0, 0, 0]), Site("Br", [0, 0, 2]))
        assert approx(bond.get_bond_order(0.5, 1.9) - 0.894736842105263) == 0

    def test_is_bonded(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1])
        assert CovalentBond.is_bonded(site1, site2)
        site2 = Site("H", [0, 0, 1.5])
        assert not CovalentBond.is_bonded(site1, site2)
        site1 = Site("U", [0, 0, 0])
        with pytest.raises(ValueError, match="No bond data for elements H - U"):
            CovalentBond.is_bonded(site1, site2)
        assert CovalentBond.is_bonded(site1, site2, default_bl=2)

    def test_str(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        assert CovalentBond(site1, site2) is not None


class TestFunc(unittest.TestCase):
    def test_get_bond_length(self):
        assert approx(get_bond_length("C", "C", 1) - 1.54) == 0
        assert approx(get_bond_length("C", "C", 2) - 1.34) == 0
        assert approx(get_bond_length("C", "H", 1) - 1.08) == 0
        assert get_bond_length("C", "H", 2) == 0.95
        assert approx(get_bond_length("C", "Br", 1) - 1.85) == 0

    def test_obtain_all_bond_lengths(self):
        assert obtain_all_bond_lengths("C", "C") == {1.0: 1.54, 2.0: 1.34, 3.0: 1.2}
        with pytest.raises(ValueError, match="No bond data for elements Br - C"):
            obtain_all_bond_lengths("Br", Element("C"))
        assert obtain_all_bond_lengths("C", Element("Br"), 1.76) == {1: 1.76}
        bond_lengths_dict = obtain_all_bond_lengths("C", "N")
        bond_lengths_dict[4] = 999
        assert obtain_all_bond_lengths("C", "N") == {1.0: 1.47, 2.0: 1.3, 3.0: 1.16}

    def test_get_bond_order(self):
        assert approx(get_bond_order("C", "C", 1) - 3) == 0
        assert approx(get_bond_order("C", "C", 1.2) - 3) == 0
        assert approx(get_bond_order("C", "C", 1.25) - 2.642857142857143) == 0
        assert approx(get_bond_order("C", "C", 1.34) - 2) == 0
        assert approx(get_bond_order("C", "C", 1.4) - 1.7) == 0  # bond length in benzene
        assert approx(get_bond_order("C", "C", 1.54) - 1) == 0
        assert approx(get_bond_order("C", "C", 2.5) - 0) == 0
        assert approx(get_bond_order("C", "C", 9999) - 0) == 0
        assert approx(get_bond_order("C", "Br", 1.9, default_bl=1.9) - 1) == 0
        assert approx(get_bond_order("C", "Br", 2, default_bl=1.9) - 0.7368421052631575) == 0
        assert approx(get_bond_order("C", "Br", 1.9, tol=0.5, default_bl=1.9) - 1) == 0
        assert approx(get_bond_order("C", "Br", 2, tol=0.5, default_bl=1.9) - 0.894736842105263) == 0
        with pytest.raises(ValueError, match="No bond data for elements Br - C"):
            get_bond_order("C", "Br", 1.9)
        assert approx(get_bond_order("N", "N", 1.25) - 2) == 0
