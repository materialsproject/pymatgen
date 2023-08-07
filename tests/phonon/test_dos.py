from __future__ import annotations

import json

from pytest import approx

from pymatgen.core.periodic_table import Element
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/NaCl_ph_dos.json") as f:
            self.dos = PhononDos.from_dict(json.load(f))
        with open(f"{TEST_FILES_DIR}/NaCl_complete_ph_dos.json") as f:
            self.structure = CompletePhononDos.from_dict(json.load(f)).structure

    def test_properties(self):
        assert self.dos.densities[15] == approx(0.0001665998)
        assert self.dos.frequencies[20] == approx(0.0894965119)
        assert self.dos.get_interpolated_value(3.0) == approx(1.2915532670115628)
        assert len(self.dos.frequencies) == 201
        assert len(self.dos.densities) == 201

    def test_get_smeared_densities(self):
        smeared = self.dos.get_smeared_densities(0.01)
        assert smeared[20] == approx(0.00084171007635058825)
        dens = self.dos.densities
        assert sum(dens) == approx(sum(smeared))

    def test_dict_methods(self):
        s = json.dumps(self.dos.as_dict())
        assert s is not None
        self.assert_msonable(self.dos)

    def test_thermodynamic_functions(self):
        assert self.dos.cv(300, structure=self.structure) == approx(48.049366665412485, abs=1e-4)
        assert self.dos.internal_energy(300, structure=self.structure) == approx(15527.596956593827, abs=1e-4)
        assert self.dos.helmholtz_free_energy(300, structure=self.structure) == approx(-6998.034212172695, abs=1e-4)
        assert self.dos.entropy(300, structure=self.structure) == approx(75.08543723748751, abs=1e-4)
        assert self.dos.zero_point_energy(structure=self.structure) == approx(4847.462485708741, abs=1e-4)


class TestCompleteDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/NaCl_complete_ph_dos.json") as f:
            self.cdos = CompletePhononDos.from_dict(json.load(f))

    def test_properties(self):
        site_Na = self.cdos.structure[0]
        site_Cl = self.cdos.structure[1]

        assert len(self.cdos.frequencies) == 201
        assert self.cdos.pdos[site_Na][30] == approx(0.008058208)
        assert self.cdos.get_site_dos(site_Na).densities[30] == approx(0.008058208)
        assert self.cdos.pdos[site_Cl][30] == approx(0.0119040783)

        assert Element.Na in self.cdos.get_element_dos()
        assert Element.Cl in self.cdos.get_element_dos()

        sum_dos = self.cdos.get_element_dos()[Element.Na] + self.cdos.get_element_dos()[Element.Cl]
        assert sum_dos.frequencies == approx(self.cdos.frequencies)
        assert sum_dos.densities == approx(self.cdos.densities)

    def test_dict_methods(self):
        s = json.dumps(self.cdos.as_dict())
        assert s is not None
        self.assert_msonable(self.cdos)

    def test_str(self):
        assert str(self.cdos) is not None
