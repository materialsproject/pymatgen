from __future__ import annotations

import json
import re

from pytest import approx

from pymatgen.core import Element
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestPhononDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/NaCl_ph_dos.json") as file:
            self.dos = PhononDos.from_dict(json.load(file))
        with open(f"{TEST_FILES_DIR}/NaCl_complete_ph_dos.json") as file:
            self.structure = CompletePhononDos.from_dict(json.load(file)).structure

    def test_repr(self):
        assert repr(self.dos) == "PhononDos(frequencies=(201,), densities=(201,), n_positive_freqs=183)"

    def test_str(self):
        assert re.match(
            r"#Frequency\s+Density\s+\n-0.66954\s+0.00000\n-0.63158\s+0.00000\n-0.59363\s+0.00000", str(self.dos)
        )

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
        json_str = json.dumps(self.dos.as_dict())
        assert json_str is not None
        self.assert_msonable(self.dos)

    def test_thermodynamic_functions(self):
        assert self.dos.cv(300, structure=self.structure) == approx(48.049366665412485, abs=1e-4)
        assert self.dos.internal_energy(300, structure=self.structure) == approx(15527.596956593827, abs=1e-4)
        assert self.dos.helmholtz_free_energy(300, structure=self.structure) == approx(-6998.034212172695, abs=1e-4)
        assert self.dos.entropy(300, structure=self.structure) == approx(75.08543723748751, abs=1e-4)
        assert self.dos.zero_point_energy(structure=self.structure) == approx(4847.462485708741, abs=1e-4)

    def test_add(self):
        dos_2x = self.dos + self.dos
        assert dos_2x.frequencies == approx(self.dos.frequencies)
        assert dos_2x.densities == approx(2 * self.dos.densities)

        dos_3x = self.dos + dos_2x
        assert dos_3x.frequencies == approx(self.dos.frequencies)
        assert dos_3x.densities == approx(3 * self.dos.densities)

        # test commutativity
        assert dos_2x + 42 == 42 + dos_2x

    def test_sub(self):
        dos_0 = self.dos - self.dos
        assert dos_0.frequencies == approx(self.dos.frequencies)
        assert dos_0.densities == approx(self.dos.densities * 0)

        dos_1 = self.dos - dos_0
        assert dos_1.frequencies == approx(self.dos.frequencies)
        assert dos_1.densities == approx(self.dos.densities)

    def test_mul(self):
        dos_2x = self.dos * 2
        assert dos_2x.frequencies == approx(self.dos.frequencies)
        assert dos_2x.densities == approx(2 * self.dos.densities)

        # test commutativity
        assert dos_2x * 1.234 == 1.234 * dos_2x

    def test_eq(self):
        assert self.dos == self.dos
        assert self.dos != 42
        assert self.dos != 2 * self.dos
        assert 2 * self.dos == self.dos + self.dos


class TestCompletePhononDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/NaCl_complete_ph_dos.json") as file:
            self.cdos = CompletePhononDos.from_dict(json.load(file))

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
        json_str = json.dumps(self.cdos.as_dict())
        assert json_str is not None
        self.assert_msonable(self.cdos)

    def test_str(self):
        assert str(self.cdos).startswith(
            "Complete phonon DOS for Full Formula (Na1 Cl1)\nReduced Formula: NaCl\n"
            "abc   :   4.023651   4.023651   4.023651\nangles"
        )
