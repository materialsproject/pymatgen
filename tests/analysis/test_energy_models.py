from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.analysis.energy_models import EwaldElectrostaticModel, IsingModel, SymmetryModel
from pymatgen.core import Species
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR


class TestEwaldElectrostaticModel(unittest.TestCase):
    def test_get_energy(self):
        coords = [[0, 0, 0], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25]]
        lattice = Lattice([[3.0, 0.0, 0.0], [1.0, 3.0, 0], [0, -2.0, 3.0]])
        struct = Structure(
            lattice,
            [
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
                {"Si4+": 0.5, "O2-": 0.25, "P5+": 0.25},
            ],
            coords,
        )

        m = EwaldElectrostaticModel()
        # large tolerance because scipy constants changed between 0.16.1 and 0.17
        assert m.get_energy(struct) == approx(-264.66364858, abs=1e-2)  # Result from GULP
        s2 = Structure.from_file(f"{TEST_FILES_DIR}/Li2O.cif")
        assert m.get_energy(s2) == approx(-145.39050015844839, abs=1e-4)

    def test_as_from_dict(self):
        m = EwaldElectrostaticModel()
        d = m.as_dict()
        assert isinstance(EwaldElectrostaticModel.from_dict(d), EwaldElectrostaticModel)


class TestSymmetryModel(unittest.TestCase):
    def test_get_energy(self):
        m = SymmetryModel()
        s2 = Structure.from_file(f"{TEST_FILES_DIR}/Li2O.cif")
        assert m.get_energy(s2) == approx(-225)

    def test_as_from_dict(self):
        m = SymmetryModel(symprec=0.2)
        d = m.as_dict()
        o = SymmetryModel.from_dict(d)
        assert isinstance(o, SymmetryModel)
        assert o.symprec == approx(0.2)


class TestIsingModel(unittest.TestCase):
    def test_get_energy(self):
        m = IsingModel(5, 6)

        struct = Structure.from_file(f"{TEST_FILES_DIR}/LiFePO4.cif")
        struct.replace_species({"Fe": Species("Fe", 2, spin=4)})
        assert m.get_energy(struct) == approx(172.81260515787977)
        struct[4] = Species("Fe", 2, spin=-4)
        struct[5] = Species("Fe", 2, spin=-4)
        assert m.get_energy(struct) == approx(51.97424405382921)

    def test_as_from_dict(self):
        m = IsingModel(5, 4)
        d = m.as_dict()
        o = IsingModel.from_dict(d)
        assert isinstance(o, IsingModel)
        assert o.j == approx(5)
