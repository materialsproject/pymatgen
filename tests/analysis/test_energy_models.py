from __future__ import annotations

from pytest import approx

from pymatgen.analysis.energy_models import EwaldElectrostaticModel, IsingModel, SymmetryModel
from pymatgen.core import Species
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR


class TestEwaldElectrostaticModel:
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

        model = EwaldElectrostaticModel()
        # large tolerance because scipy constants changed between 0.16.1 and 0.17
        assert model.get_energy(struct) == approx(-264.66364858, abs=1e-2)  # Result from GULP
        s2 = Structure.from_file(f"{TEST_FILES_DIR}/cif/Li2O.cif")
        assert model.get_energy(s2) == approx(-145.39050015844839, abs=1e-4)

    def test_as_from_dict(self):
        model = EwaldElectrostaticModel()
        dct = model.as_dict()
        restored = EwaldElectrostaticModel.from_dict(dct)
        assert isinstance(restored, EwaldElectrostaticModel)


class TestSymmetryModel:
    def test_get_energy(self):
        model = SymmetryModel()
        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/Li2O.cif")
        assert model.get_energy(struct) == approx(-225)

    def test_as_from_dict(self):
        model = SymmetryModel(symprec=0.2)
        restored = SymmetryModel.from_dict(model.as_dict())
        assert isinstance(restored, SymmetryModel)
        assert restored.symprec == approx(0.2)


class TestIsingModel:
    def test_get_energy(self):
        model = IsingModel(5, 6)

        struct = Structure.from_file(f"{TEST_FILES_DIR}/cif/LiFePO4.cif")
        struct.replace_species({"Fe": Species("Fe", 2, spin=4)})
        assert model.get_energy(struct) == approx(172.81260515787977)
        struct[4] = Species("Fe", 2, spin=-4)
        struct[5] = Species("Fe", 2, spin=-4)
        assert model.get_energy(struct) == approx(51.97424405382921)

    def test_as_from_dict(self):
        model = IsingModel(5, 4)
        restored = IsingModel.from_dict(model.as_dict())
        assert isinstance(restored, IsingModel)
        assert restored.j == approx(5)
