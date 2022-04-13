import os
import unittest

import numpy as np

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.util.testing import PymatgenTest

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__email__ = "janine.george@uclouvain.be"
__date__ = "Jan 14, 2021"

test_dir_env = os.path.join(PymatgenTest.TEST_FILES_DIR, "cohp", "environments")
this_dir = os.path.dirname(os.path.abspath(__file__))


class TestLobsterNeighbors(unittest.TestCase):
    def setUp(self):
        # test additional conditions first
        # only consider cation anion bonds

        self.chemenvlobster1 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=1,
        )

        # all bonds
        self.chemenvlobster0 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=0,
        )

        # only cation cation, anion anion bonds
        self.chemenvlobster5 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=5,
        )

        # only cation cation bonds
        self.chemenvlobster6 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=6,
        )

        # 2,3,4 are not tested so far
        self.chemenvlobster2 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=2,
        )

        self.chemenvlobster3 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=3,
        )

        self.chemenvlobster4 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=4,
        )

        # search for other testcase where 2,3,4 arrive at different results
        self.chemenvlobster0_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=0,
        )
        self.chemenvlobster1_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=1,
        )

        self.chemenvlobster2_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=2,
        )

        self.chemenvlobster5_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=5,
        )

        self.chemenvlobster5_second_percentage = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=5,
            perc_strength_ICOHP=1.0,
        )

        self.chemenvlobster6_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            additional_condition=6,
        )

        # TODO: use charge instead of valence
        self.chemenvlobster1_charges = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=1,
        )
        self.chemenvlobster1_charges_loewdin = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=1,
            which_charge="Loewdin",
        )
        self.chemenvlobster6_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=6,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster5_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=5,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster4_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=4,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster3_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=3,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster2_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=2,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster1_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=1,
            adapt_extremum_to_add_cond=True,
        )

        self.chemenvlobster0_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster0_NaSi = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.NaSi.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.NaSi.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.NaSi.gz"),
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )

    def test_use_of_coop(self):
        with self.assertRaises(ValueError):
            test = LobsterNeighbors(
                are_coops=True,
                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
                structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
                valences_from_charges=True,
                filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
                additional_condition=1,
            )

    def test_cation_anion_mode_without_ions(self):
        with self.assertRaises(ValueError) as err:
            test = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=os.path.join(test_dir_env, "../ICOHPLIST.lobster"),
                structure=Structure.from_file(os.path.join(test_dir_env, "../POSCAR")),
                valences_from_charges=False,
                additional_condition=1,
            )
        self.assertEqual(
            str(err.exception), "Valences cannot be assigned, additional_conditions 1 and 3 and 5 and 6 will not work"
        )
        with self.assertRaises(ValueError) as err:
            test = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=os.path.join(test_dir_env, "../ICOHPLIST.lobster"),
                structure=Structure.from_file(os.path.join(test_dir_env, "../POSCAR")),
                valences_from_charges=False,
                additional_condition=1,
                valences=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            )
        self.assertEqual(
            str(err.exception), "All valences are equal to 0, additional_conditions 1 and 3 and 5 and 6 will not work"
        )

    def test_wrong_additional_correction(self):
        with self.assertRaises(ValueError):
            test = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
                structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
                valences_from_charges=True,
                filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
                additional_condition=10,
            )

    def test_set_limits(self):
        test = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_353.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")),
            valences_from_charges=True,
            filename_CHARGE=os.path.join(test_dir_env, "CHARGE.lobster.mp-353.gz"),
            additional_condition=1,
            limits=[-100000, 0],
        )
        self.assertListEqual(test.limits, [-100000, 0])

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
        self.assertEqual(
            len(
                self.chemenvlobster0.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster0.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            2,
        )
        # ONLY_ANION_CATION_BONDS = 1
        self.assertEqual(
            len(
                self.chemenvlobster1.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster1.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            2,
        )
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        self.assertEqual(
            len(
                self.chemenvlobster2.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster2.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            2,
        )
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        self.assertEqual(
            len(
                self.chemenvlobster3.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster3.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            2,
        )
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        self.assertEqual(
            len(
                self.chemenvlobster4.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster4.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            2,
        )
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        self.assertEqual(
            len(
                self.chemenvlobster5.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            0,
        )
        self.assertEqual(
            len(
                self.chemenvlobster5.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            0,
        )
        # ONLY_CATION_CATION_BONDS=6
        self.assertEqual(
            len(
                self.chemenvlobster6.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=0
                )
            ),
            0,
        )

        self.assertEqual(
            len(
                self.chemenvlobster6.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")), n=1
                )
            ),
            0,
        )

        # All bonds
        # mp-353, Ag2O
        # all bonds
        self.assertEqual(
            len(
                self.chemenvlobster0_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            8,
        )

        # ONLY_ANION_CATION_BONDS = 1
        self.assertEqual(
            len(
                self.chemenvlobster1_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            2,
        )

        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        self.assertEqual(
            len(
                self.chemenvlobster2_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            2,
        )
        self.assertEqual(
            len(
                self.chemenvlobster2_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=4
                )
            ),
            4,
        )

        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        self.assertEqual(
            len(
                self.chemenvlobster5_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster5_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=4
                )
            ),
            0,
        )
        # ONLY_CATION_CATION_BONDS=6
        self.assertEqual(
            len(
                self.chemenvlobster6_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            6,
        )
        self.assertEqual(
            len(
                self.chemenvlobster6_second.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=4
                )
            ),
            0,
        )

        self.assertEqual(
            len(
                self.chemenvlobster5_second_percentage.get_nn(
                    structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
                )
            ),
            0,
        )

    def test_structure_graph(self):
        sg = self.chemenvlobster1_second.get_bonded_structure(
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz"))
        )
        self.assertEqual(type(sg), StructureGraph)

    def test_order_parameter(self):
        self.assertAlmostEqual(
            self.chemenvlobster1_second.get_local_order_parameters(
                structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_353.gz")), n=0
            )["linear"],
            1.0,
        )

    def test_get_structure_environments(self):
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        self.assertEqual(lse.coordination_environments[0][0]["ce_symbol"], "L:2")
        self.assertEqual(lse.coordination_environments[5][0]["ce_symbol"], "T:4")

        lse2 = self.chemenvlobster1.get_light_structure_environment()
        self.assertEqual(lse2.coordination_environments[0][0]["ce_symbol"], "O:6")

    def test_get_strucuture_environments_further_tests(self):
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        lse.as_dict()
        lse.get_statistics()
        self.assertTrue(lse.uniquely_determines_coordination_environments)

    def test_get_info_icohps_neighbors(self):
        results = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])
        self.assertAlmostEqual(results[0], -33.26058)
        for bond in results[1]:
            self.assertAlmostEqual(bond, -5.54345, 3)
        self.assertAlmostEqual(results[2], 6)
        self.assertAlmostEqual(results[3], ["27", "30", "48", "49", "64", "73"])

        results2 = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=None)
        self.assertAlmostEqual(results2[0], -33.26058)
        for bond in results2[1]:
            self.assertAlmostEqual(bond, -5.54345, 3)
        self.assertAlmostEqual(results2[2], 6)
        self.assertAlmostEqual(results2[3], ["27", "30", "48", "49", "64", "73"])
        self.assertAlmostEqual(
            results2[4], [["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]]
        )

    def test_get_sum_icohps_between_neighbors_of_atom(self):
        # will only look at icohps between cations or anions
        self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[1])
        self.assertEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[2], 1)
        self.assertAlmostEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[0], -0.05507)
        self.assertEqual(self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[0])[2], 15)
        # use an example where this is easier to test (e.g., linear environment?)

        chemenv_here = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp-7000.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp-7000.gz")),
            additional_condition=1,
        )
        self.assertEqual(len(chemenv_here.get_info_icohps_between_neighbors(isites=[0])[4]), 6)

    def test_get_plot_label(self):
        self.assertEqual(
            self.chemenvlobster1._get_plot_label(
                atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
                per_bond=False,
            ),
            "6 x O-Re",
        )
        self.assertEqual(
            self.chemenvlobster1._get_plot_label(
                atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
                per_bond=False,
            ),
            "5 x O-Re, 1 x O-Si",
        )

        self.assertEqual(
            self.chemenvlobster1._get_plot_label(
                atoms=[["Si1", "O2"], ["Si1", "O2"], ["Si1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
                per_bond=False,
            ),
            "4 x O-Si, 2 x O-Re",
        )

        self.assertEqual(
            self.chemenvlobster1._get_plot_label(
                atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
                per_bond=True,
            ),
            "6 x O-Re (per bond)",
        )

    def test_get_info_cohps_to_neighbors(self):
        chemenvlobster1 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=os.path.join(test_dir_env, "ICOHPLIST.lobster.mp_190_2.gz"),
            structure=Structure.from_file(os.path.join(test_dir_env, "POSCAR.mp_190.gz")),
            additional_condition=1,
        )
        self.assertEqual(
            chemenvlobster1.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"), isites=[0], only_bonds_to=["O"]
            )[0],
            "6 x O-Re (per bond)",
        )
        self.assertEqual(
            type(
                chemenvlobster1.get_info_cohps_to_neighbors(
                    path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
                    isites=[0],
                    only_bonds_to=["O"],
                )[1]
            ),
            Cohp,
        )

        cophthing = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
        )[1]
        self.assertAlmostEqual(
            np.sum([cophthing.icohp[Spin.up], cophthing.icohp[Spin.down]], axis=0)[300],
            chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])[0],
        )

        # summed_spin_channel
        cophthing = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
            summed_spin_channels=True,
        )[1]
        self.assertAlmostEqual(
            cophthing.icohp[Spin.up][300], chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])[0]
        )

        self.assertEqual(
            chemenvlobster1.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
                isites=[0],
                only_bonds_to=["Te"],
            )[0],
            None,
        )

        self.assertEqual(
            chemenvlobster1.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
                isites=[0],
                only_bonds_to=["Te"],
            )[1],
            None,
        )

        self.assertEqual(
            self.chemenvlobster0_NaSi.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.NaSi.gz"),
                isites=[8],
                onlycation_isites=False,
                only_bonds_to=["Na"],
            )[0],
            "1 x Na-Si (per bond)",
        )
        self.assertEqual(
            self.chemenvlobster0_NaSi.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.NaSi.gz"),
                isites=[8],
                onlycation_isites=False,
                only_bonds_to=["Si"],
            )[0],
            "3 x Si-Si (per bond)",
        )

        chemenvlobster1.plot_cohps_of_neighbors(
            path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
        )

        chemenvlobster1.plot_cohps_of_neighbors(
            path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
            xlim=[-10, 10],
            ylim=None,
        )

        with self.assertRaises(ValueError):
            # icohplist and cohpcar do not fit together
            self.chemenvlobster1.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

        with self.assertRaises(ValueError):
            # icohplist and cohpcar do not fit together
            self.chemenvlobster2.get_info_cohps_to_neighbors(
                path_to_COHPCAR=os.path.join(test_dir_env, "COHPCAR.lobster.mp-190.gz"),
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )


if __name__ == "__main__":
    unittest.main()
