# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.defects.generators import VacancyGenerator, \
    SubstitutionGenerator, InterstitialGenerator, VoronoiInterstitialGenerator,\
    SimpleChargeGenerator


class VacancyGeneratorTest(PymatgenTest):
    def test_vacancy_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc)

        vacs = list(vac_gen)
        self.assertEqual(len(vacs), 2)

        multiplicities = {str(v.site.specie): v.multiplicity for v in vacs}
        self.assertEqual(multiplicities, {"O": 4, "V": 2})

    def test_vacancy_gen_charges(self):
        # Ensure correct BV charges are assigned
        struc = PymatgenTest.get_structure("VO2")
        vac_gen = VacancyGenerator(struc, include_bv_charge=True)
        for vac in vac_gen:
            if str(vac.site.specie) == "V":
                self.assertEqual(vac.charge, -4)
            if str(vac.site.specie) == "O":
                self.assertEqual(vac.charge, 2)


class SubstitutionGeneratorTest(PymatgenTest):
    def test_substitution_gen(self):
        struc = PymatgenTest.get_structure("VO2")

        #test antisite
        sub_gen = SubstitutionGenerator(struc, "V")
        sub = list(sub_gen)
        self.assertEqual(len(sub), 1)
        self.assertEqual(sub[0].site.specie.symbol, 'V')
        self.assertEqual(sub[0].multiplicity, 4)
        #   test vacant site symbol
        defindex = sorted(
            struc.get_sites_in_sphere(sub[0].site.coords, 2,
                                       include_index=True),
            key=lambda x: x[1])[0][2]
        self.assertEqual(struc[defindex].specie.symbol, 'O')

        #test substitutional
        sub_gen = SubstitutionGenerator(struc, "S")
        subs = list(sub_gen)
        self.assertEqual(len(subs), 2)
        name_sets = set([s.name for s in subs])
        true_name_sets = set(['Sub_S_on_O_mult4', 'Sub_S_on_V_mult2'])
        self.assertEqual(true_name_sets, name_sets)

class InterstitialGeneratorTest(PymatgenTest):
    def test_int_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        int_gen = InterstitialGenerator(struc, "Li")

        ints = list(int_gen)
        self.assertEqual(len(ints), 4)

        multiplicities = [i.multiplicity for i in ints]
        self.assertEqual(multiplicities, [8, 8, 4, 4])

        self.assertEqual(str(ints[0].site.specie), "Li")
        self.assertEqual(str(ints[1].site.specie), "Li")

        self.assertArrayAlmostEqual(ints[0].site.coords, (0.9106, 0.3078, 0.3078), decimal=4)
        self.assertArrayAlmostEqual(ints[1].site.coords, (1.5177, 1.7444, 0.3078,), decimal=4)


class VoronoiInterstitialGeneratorTest(PymatgenTest):
    def test_int_gen(self):
        struc = PymatgenTest.get_structure("VO2")
        int_gen = VoronoiInterstitialGenerator(struc, "Li")

        ints = list(int_gen)
        self.assertEqual(len(ints), 3)

        multiplicities = [i.multiplicity for i in ints]
        self.assertEqual(multiplicities, [8, 8, 4])

        self.assertEqual(str(ints[0].site.specie), "Li")
        self.assertEqual(str(ints[1].site.specie), "Li")
        self.assertEqual(str(ints[2].site.specie), "Li")
        # self.assertEqual(str(ints[3].site.specie), "Li")

        self.assertArrayAlmostEqual(ints[0].site.coords, (1.5177146, 2.6784354, 3.9481299))
        self.assertArrayAlmostEqual(ints[1].site.coords, (1.7357692, 3.8392513, 3.8392513))
        # self.assertArrayAlmostEqual(ints[2].site.coords, (1.5177146, 3.7168193, 3.7168193))
        self.assertArrayAlmostEqual(ints[2].site.coords, (2.2765713, 2.2575138, 4.5150233))

class SimpleChargeGeneratorTest(PymatgenTest):

    def test_charge_gen(self):
        struc = PymatgenTest.get_structure("VO2")

        #assemble set of defects to get charges for
        vac_gen = VacancyGenerator(struc)
        vacs = list(vac_gen)
        full_subs = []
        for sub_elt in ['V', 'O', 'S']:
            sub_gen = SubstitutionGenerator(struc, sub_elt)
            full_subs.extend( list(sub_gen))
        int_gen = VoronoiInterstitialGenerator(struc, "H")
        inters = list(int_gen)
        defect_list = list(set().union( vacs, full_subs, inters))

        #test simple charges
        true_charges = {'Vac_O_mult4': 2, 'Int_H_Voronoi1_mult8': 0,
                        'Int_H_Voronoi2_mult8': 0, 'Vac_V_mult2': -4,
                        'Sub_S_on_V_mult2': 0, 'Int_H_Voronoi3_mult4': 0,
                        'Int_H_Voronoi4_mult4': 0, 'Sub_O_on_V_mult2': -2,
                        'Sub_S_on_O_mult4': 0, 'Sub_V_on_O_mult4': 1}
        for defect in defect_list:
            scg = SimpleChargeGenerator(defect)
            charged_defects_list = list(scg)
            def_name = charged_defects_list[0].name
            charge = charged_defects_list[0].charge
            self.assertEqual(len(charged_defects_list), 1)
            self.assertEqual( true_charges[def_name], charge)




if __name__ == "__main__":
    unittest.main()
