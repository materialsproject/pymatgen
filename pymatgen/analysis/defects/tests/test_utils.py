# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import unittest
import os
import numpy as np
import random

from pymatgen.analysis.defects.utils import QModel, eV_to_k, \
    generate_reciprocal_vectors_squared, genrecip, \
    closestsites, StructureMotifInterstitial, TopographyAnalyzer, \
    ChargeDensityAnalyzer, converge, calculate_vol, \
    tune_for_gamma, generate_R_and_G_vecs

from pymatgen.util.testing import PymatgenTest

from pymatgen.core import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Chgcar

try:
    from skimage.feature import peak_local_max
except ImportError:
    peak_local_max = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', 'chgden')


class DefectsUtilsTest(PymatgenTest):
    def test_qmodel(self):
        qm = QModel()
        modqm = QModel(beta=2., expnorm=0.5, gamma=0.1)

        # test rho_rec
        self.assertEqual(qm.rho_rec(1.), 0.77880078307140488)
        self.assertEqual(modqm.rho_rec(1.), 0.6814583156907158)

        # test rho_rec_limit0
        self.assertEqual(qm.rho_rec_limit0, -0.25)
        self.assertEqual(modqm.rho_rec_limit0, -0.51)

    def test_eV_to_k(self):
        self.assertAlmostEqual(eV_to_k(1.), 0.9681404248678961)

    def test_genrecip(self):
        a = 6.
        lattconsts = [a, a / 2., 3. * a]
        lattvectors = [[lattconsts[i] if i == j else 0. for j in range(3)] for i
                       in range(3)]
        recip_list = list( genrecip(lattvectors[0], lattvectors[1], lattvectors[2], 300))
        self.assertEqual( len(recip_list), 25620)

    def test_generate_reciprocal_vectors_squared(self):
        # test cubic case
        a = 6.
        lattvectors = [[a if i == j else 0. for j in range(3)] for i in
                       range(3)]
        brecip = [1.0966227112321507 for i in range(6)]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared(lattvectors[0],
                                                     lattvectors[1],
                                                     lattvectors[2], 1.3)),
            brecip)

        # test orthorhombic case
        lattconsts = [a, a / 2., 3. * a]
        lattvectors = [[lattconsts[i] if i == j else 0. for j in range(3)] for i
                       in range(3)]
        brval = 0.4873878716587337
        brecip = [brval, brval / 4., brval / 4., brval]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared(lattvectors[0],
                                                     lattvectors[1],
                                                     lattvectors[2], 1.)),
            brecip)

        # test triclinic case
        lattvectors = [[1.5, 0.2, 0.3], [0.3, 1.2, .2], [0.5, 0.4, 1.3]]
        brval = 24.28330561545568
        brecip = [brval, brval]
        self.assertAlmostEqual(
            list(generate_reciprocal_vectors_squared(lattvectors[0],
                                                     lattvectors[1],
                                                     lattvectors[2], 30.)),
            brecip)

    def test_closest_sites(self):
        struct = PymatgenTest.get_structure("VO2")

        # test O vacancy
        dstruct = struct.copy()
        dstruct.remove_sites([0])
        pos = struct.sites[0].coords
        bsite, dsite = closestsites(struct, dstruct, pos)
        self.assertEqual(bsite[2], 0)  # test against index
        self.assertEqual(dsite[2], 4)

        # test V vacancy
        dstruct = struct.copy()
        dstruct.remove_sites([4])
        pos = struct.sites[4].coords
        bsite, dsite = closestsites(struct, dstruct, pos)
        self.assertEqual(bsite[2], 4)  # test against index
        self.assertEqual(dsite[2], 1)

    def test_converges(self):
        self.assertAlmostEqual(converge(np.sqrt, 0.1, 0.1, 1.0),
                               0.6324555320336759)

    def test_tune_for_gamma(self):
        lattice = Lattice( [[ 4.692882, -8.12831 ,  0.],
                            [ 4.692882,  8.12831 ,  0.],
                            [ 0.,  0., 10.03391 ]])
        epsilon = 10. * np.identity(3)
        gamma = tune_for_gamma( lattice, epsilon)
        self.assertAlmostEqual(gamma, 0.19357221)

    def test_generate_R_and_G_vecs(self):
        gamma = 0.19357221
        prec = 28
        lattice = Lattice( [[ 4.692882, -8.12831 ,  0.],
                            [ 4.692882,  8.12831 ,  0.],
                            [ 0.,  0., 10.03391 ]])
        epsilon = 10. * np.identity(3)
        g_vecs, recip_summation, r_vecs, real_summation = generate_R_and_G_vecs( gamma, prec,
                                                                                 lattice, epsilon)
        self.assertEqual(len(g_vecs[0]), 16418)
        self.assertAlmostEqual(recip_summation[0], 2.8946556e-15)
        self.assertEqual(len(r_vecs[0]), 16299)
        self.assertAlmostEqual(real_summation[0], 0.00679361)

class StructureMotifInterstitialTest(PymatgenTest):
    def setUp(self):
        self.silicon = Structure(
            Lattice.from_lengths_and_angles([5.47, 5.47, 5.47],
                                            [90.0, 90.0, 90.0]),
            ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"],
            [[0.000000, 0.000000, 0.500000], [0.750000, 0.750000, 0.750000],
             [0.000000, 0.500000, 1.000000],
             [0.750000, 0.250000, 0.250000], [0.500000, 0.000000, 1.000000],
             [0.250000, 0.750000, 0.250000],
             [0.500000, 0.500000, 0.500000], [0.250000, 0.250000, 0.750000]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=False,
            site_properties=None)
        self.smi = StructureMotifInterstitial(
            self.silicon,
            "Si",
            motif_types=["tetrahedral", "octahedral"],
            op_threshs=[0.3, 0.5],
            dl=0.4,
            doverlap=1.0,
            facmaxdl=1.51)
        self.diamond = Structure(
            Lattice([[2.189, 0, 1.264], [0.73, 2.064, 1.264], [0, 0, 2.528]]),
            ["C0+", "C0+"],
            [[2.554, 1.806, 4.423], [0.365, 0.258, 0.632]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None)
        self.nacl = Structure(
            Lattice([[3.485, 0, 2.012], [1.162, 3.286, 2.012], [0, 0, 4.025]]),
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [2.324, 1.643, 4.025]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None)
        self.cscl = Structure(
            Lattice([[4.209, 0, 0], [0, 4.209, 0], [0, 0, 4.209]]),
            ["Cl1-", "Cs1+"],
            [[2.105, 2.105, 2.105], [0, 0, 0]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None)
        self.square_pyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["C", "C", "C", "C", "C", "C"],
            [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
             [0, 0, 1]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None)
        self.trigonal_bipyramid = Structure(
            Lattice([[100, 0, 0], [0, 100, 0], [0, 0, 100]]),
            ["P", "Cl", "Cl", "Cl", "Cl", "Cl"],
            [[0, 0, 0], [0, 0, 2.14], [0, 2.02, 0], [1.74937, -1.01, 0],
             [-1.74937, -1.01, 0], [0, 0, -2.14]],
            validate_proximity=False,
            to_unit_cell=False,
            coords_are_cartesian=True,
            site_properties=None)

    def test_all(self):
        self.assertIsInstance(self.smi, StructureMotifInterstitial)
        self.assertEqual(len(self.smi.enumerate_defectsites()), 1)
        self.assertIsInstance(self.smi.enumerate_defectsites()[0], PeriodicSite)
        self.assertEqual("Si",
                         self.smi.enumerate_defectsites()[0].species_string)
        self.assertEqual("tetrahedral", self.smi.get_motif_type(0))

        elem_cn_dict = self.smi.get_coordinating_elements_cns(0)
        self.assertEqual(len(list(elem_cn_dict.keys())), 1)
        self.assertEqual(list(elem_cn_dict.keys())[0], "Si")
        self.assertEqual(elem_cn_dict["Si"], 4)

        structs = self.smi.make_supercells_with_defects(np.array([1, 1, 1]))
        self.assertEqual(len(structs), 2)
        self.assertIsInstance(structs[0], Structure)

    def tearDown(self):
        del self.smi
        del self.silicon
        del self.diamond
        del self.nacl
        del self.cscl


class TopographyAnalyzerTest(unittest.TestCase):
    def setUp(self):
        feo4 = Structure.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        feo4.remove_species(["Li"])
        feo4.remove_oxidation_states()
        self.feo4 = feo4

    def test_topography_analyzer(self):
        # check interstitial sites for FePO4 using Voronoi Tessellation
        vor_feo4 = TopographyAnalyzer(self.feo4, framework_ions=["O"],
                                      cations=["P", "Fe"], check_volume=False)
        vor_feo4.cluster_nodes(tol=1.2)
        vor_feo4.remove_collisions(1.2)
        s_feo4 = vor_feo4.get_structure_with_nodes()
        sites_feo4 = np.array(
            [s_feo4[i].frac_coords for i in range(len(s_feo4)) if
             s_feo4[i].species_string == "X0+"])

        # check total number of vnodes
        self.assertAlmostEqual(len(vor_feo4.vnodes), 24)

        # check four sites that match Li sites in LiFePO4(mp-19017)
        site_predicted = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0, 0.5],
                          [0, 0.5, 0]]
        for i in range(0, 4):
            is_site_matched = False
            for site in sites_feo4:
                distance = s_feo4.lattice. \
                    get_distance_and_image(site, site_predicted[i])
                if distance[0] < 0.01:
                    is_site_matched = True
                else:
                    continue
            self.assertTrue(is_site_matched)

    def test_calculate_vol(self):
        s = Structure.from_file(os.path.join(test_dir, "LiFePO4.cif"))
        a = TopographyAnalyzer(s, framework_ions=["O"],
                               cations=["P", "Fe"], check_volume=False)
        coords = [s[i].coords for i in [20, 23, 25, 17, 24, 19]]
        vol = calculate_vol(coords=coords)
        vol_expected = 12.8884  # LiO6 volume calculated by VESTA
        self.assertAlmostEqual(vol, vol_expected, 4)


@unittest.skipIf(not peak_local_max,
                 "skimage.feature.peak_local_max module not present.")
class ChgDenAnalyzerTest(unittest.TestCase):
    def setUp(self):
        # This is a CHGCAR_sum file with reduced grid size
        chgcar_path = os.path.join(test_dir, "CHGCAR.FePO4")
        chg_FePO4 = Chgcar.from_file(chgcar_path)
        self.chgcar_path = chgcar_path
        self.chg_FePO4 = chg_FePO4
        self.ca_FePO4 = ChargeDensityAnalyzer(chg_FePO4)
        self.s_LiFePO4 = Structure.from_file(
            os.path.join(test_dir, "LiFePO4.cif"))

    def test_get_local_extrema(self):
        ca = ChargeDensityAnalyzer.from_file(self.chgcar_path)
        threshold_frac = random.random()
        threshold_abs_min = random.randrange(2, 14)
        threshold_abs_max = random.randrange(27e2, 28e4)

        # Minima test
        full_list_min = self.ca_FePO4.get_local_extrema(find_min=True,
                                                        threshold_frac=1.0)
        frac_list_min_frac = self.ca_FePO4.get_local_extrema(find_min=True,
                                                             threshold_frac=threshold_frac)
        frac_list_min_abs = self.ca_FePO4.get_local_extrema(find_min=True,
                                                            threshold_abs=threshold_abs_min)

        self.assertAlmostEqual(len(full_list_min) * threshold_frac,
                               len(frac_list_min_frac), delta=1)

        ca.get_local_extrema(find_min=True)
        df_expected = ca.extrema_df[
            ca.extrema_df["Charge Density"] <= threshold_abs_min]
        self.assertEqual(len(frac_list_min_abs), len(df_expected))

        # Maxima test
        full_list_max = self.ca_FePO4.get_local_extrema(find_min=False,
                                                        threshold_frac=1.0)
        frac_list_max = self.ca_FePO4.get_local_extrema(find_min=False,
                                                        threshold_frac=threshold_frac)
        frac_list_max_abs = self.ca_FePO4.get_local_extrema(find_min=False,
                                                            threshold_abs=threshold_abs_max)

        self.assertAlmostEqual(len(full_list_max) * threshold_frac,
                               len(frac_list_max), delta=1)

        # Local maxima should finds all center of atoms
        self.assertEqual(len(self.ca_FePO4.structure), len(full_list_max))

        ca.get_local_extrema(find_min=False)
        df_expected = ca.extrema_df[
            ca.extrema_df["Charge Density"] >= threshold_abs_max]
        self.assertEqual(len(frac_list_max_abs), len(df_expected))

    def test_remove_collisions(self):
        ca = ChargeDensityAnalyzer(self.chg_FePO4)
        ca.get_local_extrema(threshold_frac=0)
        ca.remove_collisions()  # should not trigger error
        self.assertEqual(len(ca.extrema_df), 0)

        self.ca_FePO4.get_local_extrema(find_min=False, threshold_frac=1.0)
        self.ca_FePO4.remove_collisions(min_dist=0.5)
        self.assertEqual(len(self.ca_FePO4.extrema_df), 0)

    def test_cluster_nodes(self):
        ca = ChargeDensityAnalyzer(self.chg_FePO4)
        ca.get_local_extrema()
        ca.cluster_nodes(tol=20)
        self.assertEqual(len(ca.extrema_df), 1)

    def test_get_structure_with_nodes(self):
        s_FePO4 = self.ca_FePO4.get_structure_with_nodes(find_min=True)

        sites_predicted = np.array([
            self.s_LiFePO4[i].frac_coords
            for i in range(len(self.s_LiFePO4))
            if self.s_LiFePO4[i].species_string == "Li"
        ])
        sites_guess = np.array(
            [s_FePO4[i].frac_coords for i in range(len(s_FePO4)) if
             s_FePO4[i].species_string == "X0+"])
        distances = s_FePO4.lattice.get_all_distances(sites_predicted,
                                                      sites_guess).flatten()
        distances = [d for d in distances if d < 0.1]
        self.assertEqual(len(distances), len(sites_predicted))

    def test_from_file(self):
        ca = ChargeDensityAnalyzer.from_file(self.chgcar_path)
        self.assertTrue(isinstance(ca, ChargeDensityAnalyzer))

    def test_sort_sites_by_integrated_chg(self):
        print(self.chgcar_path)
        ca = ChargeDensityAnalyzer.from_file(self.chgcar_path)
        ca.get_local_extrema()
        ca.sort_sites_by_integrated_chg()
        print(ca._extrema_df.iloc[0], 0.5)
        print(ca._extrema_df.iloc[0]['avg_charge_den'])
        self.assertAlmostEqual(ca._extrema_df.iloc[0]['a'], 0.0)
        self.assertAlmostEqual(ca._extrema_df.iloc[0]['b'], 0.5)
        self.assertAlmostEqual(ca._extrema_df.iloc[0]['c'], 0.0)
        self.assertAlmostEqual(ca._extrema_df.iloc[0]['Charge Density'],
                               1.65288944124)
        self.assertAlmostEqual(ca._extrema_df.iloc[0]['avg_charge_den'],
                               0.006831484178753711)


if __name__ == "__main__":
    unittest.main()
