#!/usr/bin/env python


__author__ = "waroquiers"

import os
import unittest

import numpy as np

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
    AllCoordinationGeometries,
)
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import (
    AbstractGeometry,
    LocalGeometryFinder,
    symmetry_measure,
)
from pymatgen.core.structure import Lattice, Structure
from pymatgen.util.testing import PymatgenTest

json_files_dir = os.path.join(
    PymatgenTest.TEST_FILES_DIR,
    "chemenv",
    "json_test_files",
)


class CoordinationGeometryFinderTest(PymatgenTest):
    def setUp(self):
        self.lgf = LocalGeometryFinder()
        self.lgf.setup_parameters(
            centering_type="standard",
            structure_refinement=self.lgf.STRUCTURE_REFINEMENT_NONE,
        )
        self.lgf2 = LocalGeometryFinder(print_citation=True)

    #     self.strategies = [SimplestChemenvStrategy(), SimpleAbundanceChemenvStrategy()]

    def test_abstract_geometry(self):
        cg_ts3 = self.lgf.allcg["TS:3"]
        cg_tet = self.lgf.allcg["T:4"]
        abstract_geom = AbstractGeometry.from_cg(cg=cg_ts3, centering_type="central_site")
        self.assertArrayAlmostEqual(abstract_geom.centre, [0.0, 0.0, 0.0])
        abstract_geom = AbstractGeometry.from_cg(cg=cg_ts3, centering_type="centroid")
        self.assertArrayAlmostEqual(abstract_geom.centre, [0.0, 0.0, 0.33333333333])
        with self.assertRaises(ValueError) as cm:
            AbstractGeometry.from_cg(
                cg=cg_ts3,
                centering_type="central_site",
                include_central_site_in_centroid=True,
            )
        self.assertEqual(
            str(cm.exception),
            "The center is the central site, no calculation of the centroid, "
            "variable include_central_site_in_centroid should be set to False",
        )
        abstract_geom = AbstractGeometry.from_cg(
            cg=cg_ts3, centering_type="centroid", include_central_site_in_centroid=True
        )
        self.assertArrayAlmostEqual(abstract_geom.centre, [0.0, 0.0, 0.25])

        # WHY ARE WE TESTING STRINGS????
        # self.assertEqual(abstract_geom.__str__(),
        #                  '\nAbstract Geometry with 3 points :\n'
        #                  '  [-1.    0.   -0.25]\n'
        #                  '  [ 1.    0.   -0.25]\n'
        #                  '  [ 0.   0.   0.75]\n'
        #                  'Points are referenced to the centroid (calculated with the central site) :\n'
        #                  '  [ 0.   0.   0.25]\n')

        symm_dict = symmetry_measure([[0.0, 0.0, 0.0]], [1.1, 2.2, 3.3])
        self.assertAlmostEqual(symm_dict["symmetry_measure"], 0.0)
        self.assertEqual(symm_dict["scaling_factor"], None)
        self.assertEqual(symm_dict["rotation_matrix"], None)

        tio2_struct = self.get_structure("TiO2")

        envs = self.lgf.compute_coordination_environments(structure=tio2_struct, indices=[0])
        self.assertAlmostEqual(envs[0][0]["csm"], 1.5309987846957258)
        self.assertAlmostEqual(envs[0][0]["ce_fraction"], 1.0)
        self.assertEqual(envs[0][0]["ce_symbol"], "O:6")
        self.assertEqual(sorted(envs[0][0]["permutation"]), sorted([0, 4, 1, 5, 2, 3]))

        self.lgf.setup_random_structure(coordination=5)
        self.assertEqual(len(self.lgf.structure), 6)

        self.lgf.setup_random_indices_local_geometry(coordination=5)
        self.assertEqual(self.lgf.icentral_site, 0)
        self.assertEqual(len(self.lgf.indices), 5)

        self.lgf.setup_ordered_indices_local_geometry(coordination=5)
        self.assertEqual(self.lgf.icentral_site, 0)
        self.assertEqual(self.lgf.indices, list(range(1, 6)))

        self.lgf.setup_explicit_indices_local_geometry(explicit_indices=[3, 5, 2, 0, 1, 4])
        self.assertEqual(self.lgf.icentral_site, 0)
        self.assertEqual(self.lgf.indices, [4, 6, 3, 1, 2, 5])

        LiFePO4_struct = self.get_structure("LiFePO4")
        isite = 10
        envs_LiFePO4 = self.lgf.compute_coordination_environments(structure=LiFePO4_struct, indices=[isite])
        self.assertAlmostEqual(envs_LiFePO4[isite][0]["csm"], 0.140355832317)
        nbs_coords = [
            np.array([6.16700437, -4.55194317, -5.89031356]),
            np.array([4.71588167, -4.54248093, -3.75553856]),
            np.array([6.88012571, -5.79877503, -3.73177541]),
            np.array([6.90041188, -3.32797839, -3.71812416]),
        ]

        # test to check that one can pass voronoi_distance_cutoff
        struct = Structure(Lattice.cubic(25), ["O", "C", "O"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.17], [0.0, 0.0, 2.34]])
        self.lgf.setup_structure(structure=struct)
        self.lgf.compute_structure_environments(voronoi_distance_cutoff=25)

        self.lgf.setup_structure(LiFePO4_struct)
        self.lgf.setup_local_geometry(isite, coords=nbs_coords)

        perfect_tet = AbstractGeometry.from_cg(
            cg=cg_tet, centering_type="centroid", include_central_site_in_centroid=False
        )
        points_perfect_tet = perfect_tet.points_wcs_ctwcc()
        res = self.lgf.coordination_geometry_symmetry_measures_fallback_random(
            coordination_geometry=cg_tet, NRANDOM=5, points_perfect=points_perfect_tet
        )
        (
            permutations_symmetry_measures,
            permutations,
            algos,
            local2perfect_maps,
            perfect2local_maps,
        ) = res
        for perm_csm_dict in permutations_symmetry_measures:
            self.assertAlmostEqual(perm_csm_dict["symmetry_measure"], 0.140355832317)

    #
    # def _strategy_test(self, strategy):
    #     files = []
    #     for (dirpath, dirnames, filenames) in os.walk(json_files_dir):
    #         files.extend(filenames)
    #         break
    #
    #     for ifile, json_file in enumerate(files):
    #         with self.subTest(json_file=json_file):
    #             with open("{}/{}".format(json_files_dir, json_file), 'r') as f:
    #                 dd = json.load(f)
    #
    #             atom_indices = dd['atom_indices']
    #             expected_geoms = dd['expected_geoms']
    #
    #             struct = Structure.from_dict(dd['structure'])
    #
    #             struct = self.lgf.setup_structure(struct)
    #             se = self.lgf.compute_structure_environments_detailed_voronoi(only_indices=atom_indices,
    #                                                                           maximum_distance_factor=1.5)
    #
    #             #All strategies should get the correct environment with their default parameters
    #             strategy.set_structure_environments(se)
    #             for ienv, isite in enumerate(atom_indices):
    #                 ce = strategy.get_site_coordination_environment(struct[isite])
    #                 try:
    #                     coord_env = ce[0]
    #                 except TypeError:
    #                     coord_env = ce
    #                 #Check that the environment found is the expected one
    #                 self.assertEqual(coord_env, expected_geoms[ienv])
    #
    # def test_simplest_chemenv_strategy(self):
    #     strategy = SimplestChemenvStrategy()
    #     self._strategy_test(strategy)
    #
    # def test_simple_abundance_chemenv_strategy(self):
    #     strategy = SimpleAbundanceChemenvStrategy()
    #     self._strategy_test(strategy)

    def test_perfect_environments(self):
        allcg = AllCoordinationGeometries()
        indices_CN = {
            1: [0],
            2: [1, 0],
            3: [1, 0, 2],
            4: [2, 0, 3, 1],
            5: [2, 3, 1, 0, 4],
            6: [0, 2, 3, 1, 5, 4],
            7: [2, 6, 0, 3, 4, 5, 1],
            8: [1, 2, 6, 3, 7, 0, 4, 5],
            9: [5, 2, 6, 0, 4, 7, 3, 8, 1],
            10: [8, 5, 6, 3, 0, 7, 2, 4, 9, 1],
            11: [7, 6, 4, 1, 2, 5, 0, 8, 9, 10, 3],
            12: [5, 8, 9, 0, 3, 1, 4, 2, 6, 11, 10, 7],
            13: [4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10],
            20: [8, 12, 11, 0, 14, 10, 13, 6, 18, 1, 9, 17, 3, 19, 5, 7, 15, 2, 16, 4],
        }

        for coordination in range(1, 21):
            for mp_symbol in allcg.get_implemented_geometries(coordination=coordination, returned="mp_symbol"):
                cg = allcg.get_geometry_from_mp_symbol(mp_symbol=mp_symbol)
                self.lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
                self.lgf.setup_test_perfect_environment(
                    mp_symbol,
                    randomness=False,
                    indices=indices_CN[coordination],
                    random_translation="NONE",
                    random_rotation="NONE",
                    random_scale="NONE",
                )
                se = self.lgf.compute_structure_environments(
                    only_indices=[0],
                    maximum_distance_factor=1.01 * cg.distfactor_max,
                    min_cn=cg.coordination_number,
                    max_cn=cg.coordination_number,
                    only_symbols=[mp_symbol],
                )
                self.assertAlmostEqual(
                    se.get_csm(0, mp_symbol)["symmetry_measure"],
                    0.0,
                    delta=1e-8,
                    msg=f"Failed to get perfect environment with mp_symbol {mp_symbol}",
                )

    def test_disable_hints(self):
        allcg = AllCoordinationGeometries()
        mp_symbol = "SH:13"
        mp_symbols = ["SH:13", "HP:12"]
        cg = allcg.get_geometry_from_mp_symbol(mp_symbol=mp_symbol)
        mypoints = cg.points
        mypoints[-1] = [0.9 * cc for cc in mypoints[-1]]
        self.lgf.allcg = AllCoordinationGeometries(only_symbols=[mp_symbol])
        self.lgf.setup_test_perfect_environment(
            mp_symbol,
            randomness=False,
            indices=[4, 11, 5, 12, 1, 2, 8, 3, 0, 6, 9, 7, 10],
            random_translation="NONE",
            random_rotation="NONE",
            random_scale="NONE",
            points=mypoints,
        )
        se_nohints = self.lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.02 * cg.distfactor_max,
            min_cn=12,
            max_cn=13,
            only_symbols=mp_symbols,
            get_from_hints=False,
        )
        se_hints = self.lgf.compute_structure_environments(
            only_indices=[0],
            maximum_distance_factor=1.02 * cg.distfactor_max,
            min_cn=12,
            max_cn=13,
            only_symbols=mp_symbols,
            get_from_hints=True,
        )
        with self.assertRaises(KeyError):
            abc = se_nohints.ce_list[0][12]
            abc.minimum_geometries()
        self.assertAlmostEqual(se_hints.ce_list[0][13][0], se_nohints.ce_list[0][13][0])
        self.assertTrue(set(se_nohints.ce_list[0].keys()).issubset(set(se_hints.ce_list[0].keys())))


if __name__ == "__main__":
    unittest.main()
