# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import os
import unittest
import numpy as np

from pymatgen.core import Structure
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.defects.core import Vacancy, Interstitial, Substitution, \
    DefectEntry, create_saturated_interstitial_structure
from pymatgen.util.testing import PymatgenTest

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

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

        # test that structure has all velocities equal to [0., 0., 0.](previously caused failures for structure printing)
        # self.assertTrue( (np.array(sc_scaled_struc.site_properties['velocities']) == 0.).all())
        # self.assertEqual( len(sc_scaled_struc.site_properties['velocities']), len(sc_scaled_struc))

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
        self.assertEqual(interstitial.multiplicity, 8.0)

        # test manual setting of multiplicity
        interstitial = Interstitial(struc, int_site, multiplicity=4.0)
        self.assertEqual(interstitial.multiplicity, 4.0)

        # Test composoition
        self.assertEqual(dict(interstitial.defect_composition.as_dict()), {"V": 3, "O": 4})

        # test that structure has all velocities equal if velocities previously existed
        # (previously caused failures for structure printing)
        vel_struc = Structure( struc.lattice, struc.species, struc.frac_coords,
                               site_properties= {'velocities': [[0., 0., 0.]]*len(struc) } )
        interstitial = Interstitial(vel_struc, int_site, charge=-1.0)
        int_struc = interstitial.generate_defect_structure(1)

        self.assertTrue( (np.array(int_struc.site_properties['velocities']) == 0.).all())
        self.assertEqual( len(int_struc.site_properties['velocities']), len(int_struc))

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

        # Test composition
        self.assertEqual(dict(substitution.defect_composition.as_dict()), {"V": 2, "Sr": 1, "O": 3})

        # test that structure has all velocities equal if velocities previously existed
        # (previously caused failures for structure printing)
        vel_struc = Structure( struc.lattice, struc.species, struc.frac_coords,
                               site_properties= {'velocities': [[0., 0., 0.]]*len(struc) } )
        substitution = Substitution(vel_struc, sub_site)
        sub_struc = substitution.generate_defect_structure(1)

        self.assertTrue( (np.array(sub_struc.site_properties['velocities']) == 0.).all())
        self.assertEqual( len(sub_struc.site_properties['velocities']), len(sub_struc))


class create_saturated_interstitial_structureTest(PymatgenTest):

    def test_sublattice_generation(self):
        struc = PymatgenTest.get_structure("CsCl")
        sc_struc = struc.copy()
        sc_struc.make_supercell(3)

        # test for vacancy and sub (should not change structure)
        Cs_index = sc_struc.indices_from_symbol("Cs")[0]
        cs_vac = Vacancy(sc_struc, sc_struc[Cs_index])
        decorated_cs_vac = create_saturated_interstitial_structure(cs_vac)
        self.assertEqual(len(decorated_cs_vac), len(sc_struc))

        Cl_index = sc_struc.indices_from_symbol("Cl")[0]

        cl_vac = Vacancy(sc_struc, sc_struc[Cl_index])
        decorated_cl_vac = create_saturated_interstitial_structure(cl_vac)
        self.assertEqual(len(decorated_cl_vac), len(sc_struc))

        sub_site = PeriodicSite("Sr", sc_struc[Cs_index].coords, sc_struc.lattice,
                                coords_are_cartesian=True)

        sub = Substitution(sc_struc, sub_site)
        decorated_sub = create_saturated_interstitial_structure(sub)
        self.assertEqual(len(decorated_sub), len(sc_struc))

        # test interstitial in symmorphic structure type
        inter_site = PeriodicSite("H", [0., 1.05225, 2.1045], struc.lattice,
                                  coords_are_cartesian=True)  # voronoi type
        interstitial = Interstitial(struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 14)

        inter_site = PeriodicSite("H", [0.10021429, 0.10021429, 2.1045], struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 14)

        inter_site = PeriodicSite("H", [4.10878571, 1.10235714, 2.1045], struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 26)

        inter_site = PeriodicSite("H", [0., 0., 0.5], struc.lattice,
                                  coords_are_cartesian=False)  # a reasonable guess type
        interstitial = Interstitial(struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 5)

        # test interstitial in non-symmorphic structure type
        # (voronoi and InFit generator of different types...)
        ns_struc = Structure.from_file(os.path.join(test_dir, "CuCl.cif"))

        inter_site = PeriodicSite("H", [0.45173594, 0.41157895, 5.6604067], ns_struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(ns_struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 40)

        inter_site = PeriodicSite("H", [0.47279906, 0.82845998, 5.62015285], ns_struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(ns_struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 40)

        inter_site = PeriodicSite("H", [0.70845255, 6.50298148, 5.16979425], ns_struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(ns_struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 40)

        inter_site = PeriodicSite("H", [0.98191329, 0.36460337, 4.64718203], ns_struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(ns_struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 40)

        inter_site = PeriodicSite("H", [0.39286561, 3.92702149, 1.05802631], ns_struc.lattice,
                                  coords_are_cartesian=True)  # InFit type
        interstitial = Interstitial(ns_struc, inter_site)
        decorated_inter = create_saturated_interstitial_structure(interstitial)
        self.assertEqual(len(decorated_inter), 40)


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
