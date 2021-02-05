import os
import unittest

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.lobster_env import LobsterNeighbors
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__email__ = "janine.george@uclouvain.be"
__date__ = "Jan 14, 2021"

test_dir_env = os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp/environments")
this_dir = os.path.dirname(os.path.abspath(__file__))


class TestLobsterNeighbors(unittest.TestCase):
    def setUp(self):
        # test additional conditions first
        # only consider cation anion bonds

        self.chemenvlobster1 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=1)

        # all bonds
        self.chemenvlobster0 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=0)

        # only cation cation, anion anion bonds
        self.chemenvlobster5 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=5)

        # only cation cation bonds
        self.chemenvlobster6 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=6)

        # 2,3,4 are not tested so far
        self.chemenvlobster2 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=2)

        self.chemenvlobster3 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=3)

        self.chemenvlobster4 = LobsterNeighbors(are_coops=False,
                                                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190"),
                                                structure=Structure.from_file(
                                                    os.path.join(test_dir_env, "POSCAR.mp_190")),
                                                additional_condition=4)

        # search for other testcase where 2,3,4 arrive at different results
        self.chemenvlobster0_second = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                    "ICOHPLIST.lobster.mp_353"),
                                                       structure=Structure.from_file(
                                                           os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                       additional_condition=0)
        self.chemenvlobster1_second = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                    "ICOHPLIST.lobster.mp_353"),
                                                       structure=Structure.from_file(
                                                           os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                       additional_condition=1)

        self.chemenvlobster2_second = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                    "ICOHPLIST.lobster.mp_353"),
                                                       structure=Structure.from_file(
                                                           os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                       additional_condition=2)

        self.chemenvlobster5_second = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                    "ICOHPLIST.lobster.mp_353"),
                                                       structure=Structure.from_file(
                                                           os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                       additional_condition=5)

        self.chemenvlobster5_second_percentage = LobsterNeighbors(are_coops=False,
                                                                  filename_ICOHP=os.path.join(test_dir_env,
                                                                                              "ICOHPLIST.lobster.mp_353"),
                                                                  structure=Structure.from_file(
                                                                      os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                                  additional_condition=5, perc_strength_ICOHP=1.0)

        self.chemenvlobster6_second = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                    "ICOHPLIST.lobster.mp_353"),
                                                       structure=Structure.from_file(
                                                           os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                       additional_condition=6)

        # TODO: use charge instead of valence
        self.chemenvlobster1_charges = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                     "ICOHPLIST.lobster.mp_353"),
                                                        structure=Structure.from_file(
                                                            os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                        valences_from_charges=True,
                                                        filename_CHARGE=os.path.join(test_dir_env,
                                                                                     "CHARGE.lobster.mp-353.gz"),
                                                        additional_condition=1)

    def test_use_of_coop(self):
        #TODO: check if ValueError is raised
        with self.assertRaises(ValueError):
            test=LobsterNeighbors(are_coops=True, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                         "ICOHPLIST.lobster.mp_353"),
                                                            structure=Structure.from_file(
                                                                os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                            valences_from_charges=True,
                                                            filename_CHARGE=os.path.join(test_dir_env,
                                                                                         "CHARGE.lobster.mp-353.gz"),
                                                            additional_condition=1)

    def test_wrong_additional_correction(self):
        with self.assertRaises(ValueError):
            test=LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                                                         "ICOHPLIST.lobster.mp_353"),
                                                            structure=Structure.from_file(
                                                                os.path.join(test_dir_env, "POSCAR.mp_353")),
                                                            valences_from_charges=True,
                                                            filename_CHARGE=os.path.join(test_dir_env,
                                                                                         "CHARGE.lobster.mp-353.gz"),
                                                            additional_condition=10)
    def test_no_bva_possible(self):
        pass
        #TODO: write something where not BondValenceAnalysis can be done (provide a certain structure where it cannot
        # be done)


    def test_set_limits(self):
        test = LobsterNeighbors(are_coops=False, filename_ICOHP=os.path.join(test_dir_env,
                                                                             "ICOHPLIST.lobster.mp_353"),
                                structure=Structure.from_file(
                                    os.path.join(test_dir_env, "POSCAR.mp_353")),
                                valences_from_charges=True,
                                filename_CHARGE=os.path.join(test_dir_env,
                                                             "CHARGE.lobster.mp-353.gz"),
                                additional_condition=1,limits=[-100000,0])

    def test_molecules_allowed(self):
        self.chemenvlobster1.molecules_allowed

    def test_get_nn_info(self):
        # NO_ADDITIONAL_CONDITION = 0
        # ONLY_ANION_CATION_BONDS = 1
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        # ONLY_CATION_CATION_BONDS=6

        # All bonds
        # ReO3
        self.assertEqual(len(self.chemenvlobster0.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 6)
        self.assertEqual(len(self.chemenvlobster0.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 2)
        # ONLY_ANION_CATION_BONDS = 1
        self.assertEqual(len(self.chemenvlobster1.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 6)
        self.assertEqual(len(self.chemenvlobster1.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 2)
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        self.assertEqual(len(self.chemenvlobster2.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 6)
        self.assertEqual(len(self.chemenvlobster2.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 2)
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        self.assertEqual(len(self.chemenvlobster3.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 6)
        self.assertEqual(len(self.chemenvlobster3.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 2)
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        self.assertEqual(len(self.chemenvlobster4.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 6)
        self.assertEqual(len(self.chemenvlobster4.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 2)
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        self.assertEqual(len(self.chemenvlobster5.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 0)
        self.assertEqual(len(self.chemenvlobster5.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=1)), 0)
        # ONLY_CATION_CATION_BONDS=6
        self.assertEqual(len(self.chemenvlobster6.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_190")), n=0)), 0)
        self.assertEqual(len(self.chemenvlobster6.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                    "POSCAR.mp_190")),
                                                         n=1)), 0)

        # All bonds
        # mp-353, Ag2O
        # all bonds
        self.assertEqual(len(self.chemenvlobster0_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=0)), 8)

        # ONLY_ANION_CATION_BONDS = 1
        self.assertEqual(len(self.chemenvlobster1_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=0)), 2)

        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        self.assertEqual(len(self.chemenvlobster2_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=0)), 2)
        self.assertEqual(len(self.chemenvlobster2_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=4)), 4)

        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        self.assertEqual(len(self.chemenvlobster5_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=0)), 6)
        self.assertEqual(len(self.chemenvlobster5_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=4)), 0)
        # ONLY_CATION_CATION_BONDS=6
        self.assertEqual(len(self.chemenvlobster6_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=0)), 6)
        self.assertEqual(len(self.chemenvlobster6_second.get_nn(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                           "POSCAR.mp_353")),
                                                                n=4)), 0)

        self.assertEqual(len(self.chemenvlobster5_second_percentage.get_nn(structure=Structure.from_file(
            os.path.join(test_dir_env,
                         "POSCAR.mp_353")),
            n=0)), 0)

    def test_structure_graph(self):
        sg = self.chemenvlobster1_second.get_bonded_structure(structure=Structure.from_file(os.path.join(test_dir_env,
                                                                                                         "POSCAR.mp_353")))
        self.assertEqual(type(sg), StructureGraph)

    def test_order_parameter(self):
        self.assertAlmostEqual(
            self.chemenvlobster1_second.get_local_order_parameters(structure=Structure.from_file(os.path.join(
                test_dir_env,
                "POSCAR.mp_353")),
                n=0)["linear"], 1.0)

    def test_get_structure_environments(self):
        # TODO: make sure Davids bugfix is already implemented
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        self.assertEqual(lse.coordination_environments[0][0]['ce_symbol'], 'L:2')
        self.assertEqual(lse.coordination_environments[5][0]['ce_symbol'], 'T:4')

        lse2 = self.chemenvlobster1.get_light_structure_environment()
        self.assertEqual(lse2.coordination_environments[0][0]['ce_symbol'], 'O:6')

    def test_get_strucuture_environments_further_tests(self):
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        lse.as_dict()
        lse.get_statistics()

    def test_get_info_icohps_neighbors(self):
        results = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])
        self.assertAlmostEqual(results[0], -33.26058)
        for bond in results[1]:
            self.assertAlmostEqual(bond, -5.54345, 3)
        self.assertAlmostEqual(results[2], 6)
        self.assertAlmostEqual(results[3], ['27', '30', '48', '49', '64', '73'])

        results2 = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[])
        self.assertAlmostEqual(results2[0], -33.26058)
        for bond in results2[1]:
            self.assertAlmostEqual(bond, -5.54345, 3)
        self.assertAlmostEqual(results2[2], 6)
        self.assertAlmostEqual(results2[3], ['27', '30', '48', '49', '64', '73'])

    def test_get_sum_icohps_between_neighbors_of_atom(self):
        # will only look at icohps between cations or anions
        self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[1])
        self.assertEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[2], 1)
        self.assertAlmostEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[0], -0.05507)
        self.assertEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[0])[2], 15)



if __name__ == "__main__":
    unittest.main()
