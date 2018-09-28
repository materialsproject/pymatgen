# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest

from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution, DefectEntry
from pymatgen.util.testing import PymatgenTest


class DefectsCoreTest(PymatgenTest):
    def test_vacancy(self):
        struc = PymatgenTest.get_structure("VO2")
        V_index = struc.indices_from_symbol("V")[0]
        vac = Vacancy(struc, struc[V_index])

        # test generation and super cell
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 1, "O": 4})

        vac_struc = vac.generate_defect_structure(2)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 15, "O": 32})

        vac_struc = vac.generate_defect_structure(3)
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 53, "O": 108})

        vac_struc = vac.generate_defect_structure([[2., 0, 0], [0, 0, -3.], [0, 2., 0]])
        self.assertEqual(vac_struc.composition.as_dict(), {"V": 23, "O": 48})

        # test charge
        vac = Vacancy(struc, struc[V_index])
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, 0.0)

        vac = Vacancy(struc, struc[V_index], charge=1.0)
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, 1.0)

        vac = Vacancy(struc, struc[V_index], charge=-1.0)
        vac_struc = vac.generate_defect_structure(1)
        self.assertEqual(vac_struc.charge, -1.0)

        # test multiplicity
        vac = Vacancy(struc, struc[V_index])
        self.assertEqual(vac.multiplicity, 2)

        O_index = struc.indices_from_symbol("O")[0]
        vac = Vacancy(struc, struc[O_index])
        self.assertEqual(vac.multiplicity, 4)

        # Test composition
        self.assertEqual(dict(vac.defect_composition.as_dict()), {"V": 2, "O": 3})

        # test lattice value error occurs for differnet lattices
        sc_scaled_struc = struc.copy()
        sc_scaled_struc.make_supercell(2)
        self.assertRaises( ValueError, Vacancy, struc, sc_scaled_struc[V_index])
        self.assertRaises( ValueError, Vacancy, sc_scaled_struc, struc[V_index])

    def test_interstitial(self):
        struc = PymatgenTest.get_structure("VO2")
        V_index = struc.indices_from_symbol("V")[0]

        int_site = PeriodicSite("V", struc[V_index].coords + [0.1, 0.1, 0.1], struc.lattice)
        interstitial = Interstitial(struc, int_site)

        # test generation and super cell
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 3, "O": 4})
        # Ensure the site is in the right place
        self.assertEqual(int_site, int_struc.get_sites_in_sphere(int_site.coords, 0.1)[0][0])

        int_struc = interstitial.generate_defect_structure(2)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 17, "O": 32})

        int_struc = interstitial.generate_defect_structure(3)
        self.assertEqual(int_struc.composition.as_dict(), {"V": 55, "O": 108})

        int_struc = interstitial.generate_defect_structure([[2., 0, 0], [0, 0, -3.], [0, 2., 0]])
        self.assertEqual(int_struc.composition.as_dict(), {"V": 25, "O": 48})

        # test charge
        interstitial = Interstitial(struc, int_site)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, 0.0)

        interstitial = Interstitial(struc, int_site, charge=1.0)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, 1.0)

        interstitial = Interstitial(struc, int_site, charge=-1.0)
        int_struc = interstitial.generate_defect_structure(1)
        self.assertEqual(int_struc.charge, -1.0)

        # test multiplicity
        interstitial = Interstitial(struc, int_site)
        self.assertEqual(interstitial.multiplicity, 1.0)

        interstitial = Interstitial(struc, int_site, multiplicity=4.0)
        self.assertEqual(interstitial.multiplicity, 4.0)

        # Test composoition
        self.assertEqual(dict(interstitial.defect_composition.as_dict()), {"V": 3, "O": 4})

    def test_substitution(self):
        struc = PymatgenTest.get_structure("VO2")
        V_index = struc.indices_from_symbol("V")[0]

        sub_site = PeriodicSite("Sr", struc[V_index].coords, struc.lattice, coords_are_cartesian=True)
        substitution = Substitution(struc, sub_site)

        # test generation and super cell
        sub_struc = substitution.generate_defect_structure(1)
        self.assertEqual(sub_struc.composition.as_dict(), {"V": 1, "Sr": 1, "O": 4})

        sub_struc = substitution.generate_defect_structure(2)
        self.assertEqual(sub_struc.composition.as_dict(), {"V": 15, "Sr": 1, "O": 32})

        sub_struc = substitution.generate_defect_structure(3)
        self.assertEqual(sub_struc.composition.as_dict(), {"V": 53, "Sr": 1, "O": 108})

        sub_struc = substitution.generate_defect_structure([[2., 0, 0], [0, 0, -3.], [0, 2., 0]])
        self.assertEqual(sub_struc.composition.as_dict(), {"V": 23, "O": 48, "Sr": 1})

        # test charge
        substitution = Substitution(struc, sub_site)
        sub_struc = substitution.generate_defect_structure(1)
        self.assertEqual(sub_struc.charge, 0.0)

        substitution = Substitution(struc, sub_site, charge=1.0)
        sub_struc = substitution.generate_defect_structure(1)
        self.assertEqual(sub_struc.charge, 1.0)

        substitution = Substitution(struc, sub_site, charge=-1.0)
        sub_struc = substitution.generate_defect_structure(1)
        self.assertEqual(sub_struc.charge, -1.0)

        # test multiplicity
        substitution = Substitution(struc, sub_site)
        self.assertEqual(substitution.multiplicity, 2.0)

        O_index = struc.indices_from_symbol("O")[0]
        sub_site = PeriodicSite("Sr", struc[O_index].coords, struc.lattice, coords_are_cartesian=True)
        substitution = Substitution(struc, sub_site)
        self.assertEqual(substitution.multiplicity, 4)

        # Test composoition
        self.assertEqual(dict(substitution.defect_composition.as_dict()), {"V": 2, "Sr": 1, "O": 3})


class DefectEntryTest(PymatgenTest):
    def setUp(self):
        self.struc = PymatgenTest.get_structure("VO2")
        V_index = self.struc.indices_from_symbol("V")[0]
        sub_site = PeriodicSite("Sr", self.struc[V_index].coords, self.struc.lattice, coords_are_cartesian=True)
        self.substitution = Substitution(self.struc, sub_site)

    def test_init(self):
        entry = DefectEntry(self.substitution, 2.5)
        entry_doc = entry.as_dict()
        re_entry = DefectEntry.from_dict(entry_doc)
        self.assertNotEqual(re_entry, None)

    def test_corrections(self):
        entry = DefectEntry(self.substitution, 2.5)

        self.assertAlmostEqual(entry.energy, 2.5)

        entry.corrections["pot_corr"] = -0.3
        self.assertAlmostEqual(entry.energy, 2.2)

    def test_formation_energy(self):
        entry = DefectEntry(self.substitution, 2.5, corrections={"pot_corr": -0.3})

        # Test chemical potentials on formation energy
        self.assertAlmostEqual(entry.formation_energy(), 2.2)
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2}), 2.0)
        self.assertAlmostEqual(entry.formation_energy({"V": 0.2}), 2.4)
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2, "V": 0.2}), 2.2)
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2, "V": 0.2, "O": 2}), 2.2)

        # Test Fermi level on formation energy
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2, "V": 0.2}, fermi_level=0.2), 2.2)
        entry.parameters["vbm"] = 0
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2, "V": 0.2}, fermi_level=0.2), 2.2)
        entry.defect._charge = 1
        self.assertAlmostEqual(entry.formation_energy({"Sr": 0.2, "V": 0.2}, fermi_level=0.2), 2.4)

    def test_defect_concentration(self):
        entry = DefectEntry(self.substitution, .5, corrections={})
        entry.defect._charge = -1

        chem_pots = {"Sr": 0., "V": 0., "O": 0.}
        self.assertAlmostEqual(entry.defect_concentration(chem_pots), 1.2878309944593931e14)

        # #test temperature dependence
        self.assertAlmostEqual(entry.defect_concentration(chem_pots, temperature=600), 2.040208007417593e+18)

        # test fermi level dependence
        self.assertAlmostEqual(entry.defect_concentration(chem_pots, fermi_level=.3), 1.4113592133771723e+19)


if __name__ == "__main__":
    unittest.main()
