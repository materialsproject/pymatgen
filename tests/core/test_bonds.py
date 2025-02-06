from __future__ import annotations

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


class TestCovalentBond:
    def test_length(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0.7, 0.6])
        assert CovalentBond(site1, site2).length == approx(0.92195444572928864)

    def test_get_bond_order(self):
        site1 = Site("C", [0, 0, 0])
        site2 = Site("H", [0, 0, 1.08])
        assert CovalentBond(site1, site2).get_bond_order() == approx(1)
        bond = CovalentBond(Site("C", [0, 0, 0]), Site("Br", [0, 0, 2]))
        assert bond.get_bond_order(0.5, 1.9) == approx(0.894736842105263)

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


class TestFunc:
    def test_get_bond_length(self):
        assert get_bond_length("C", "C", 1) == approx(1.54)
        assert get_bond_length("C", "C", 2) == approx(1.34)
        assert get_bond_length("C", "H", 1) == approx(1.08)
        assert get_bond_length("C", "H", 2) == approx(0.95)
        assert get_bond_length("C", "Br", 1) == approx(1.85)

    def test_obtain_all_bond_lengths(self):
        assert obtain_all_bond_lengths("C", "C") == approx({1.0: 1.54, 2.0: 1.34, 3.0: 1.2})
        with pytest.raises(ValueError, match="No bond data for elements Br - C"):
            obtain_all_bond_lengths("Br", Element("C"))
        assert obtain_all_bond_lengths("C", Element("Br"), 1.76) == approx({1: 1.76})
        bond_lengths_dict = obtain_all_bond_lengths("C", "N")
        bond_lengths_dict[4] = 999
        assert obtain_all_bond_lengths("C", "N") == approx({1.0: 1.47, 2.0: 1.3, 3.0: 1.16})

    def test_get_bond_order(self):
        assert get_bond_order("C", "C", 1) == approx(3)
        assert get_bond_order("C", "C", 1.2) == approx(3)
        assert get_bond_order("C", "C", 1.25) == approx(2.642857142857143)
        assert get_bond_order("C", "C", 1.34) == approx(2)
        assert get_bond_order("C", "C", 1.4) == approx(1.7)  # bond length in benzene
        assert get_bond_order("C", "C", 1.54) == approx(1)
        assert get_bond_order("C", "C", 2.5) == approx(0)
        assert get_bond_order("C", "C", 9999) == approx(0)
        assert get_bond_order("C", "Br", 1.9, default_bl=1.9) == approx(1)
        assert get_bond_order("C", "Br", 2, default_bl=1.9) == approx(0.7368421052631575)
        assert get_bond_order("C", "Br", 1.9, tol=0.5, default_bl=1.9) == approx(1)
        assert get_bond_order("C", "Br", 2, tol=0.5, default_bl=1.9) == approx(0.894736842105263)
        with pytest.raises(ValueError, match="No bond data for elements Br - C"):
            get_bond_order("C", "Br", 1.9)
        assert get_bond_order("N", "N", 1.25) == approx(2)
