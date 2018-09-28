# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import json
from io import open
import warnings

from pymatgen.electronic_structure.bandstructure import Kpoint
from pymatgen import Lattice
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.bandstructure import (BandStructureSymmLine,
                                                         get_reconstructed_band_structure)
from pymatgen.util.testing import PymatgenTest

from monty.serialization import loadfn

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class KpointTest(unittest.TestCase):

    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.kpoint = Kpoint([0.1, 0.4, -0.5], self.lattice, label="X")

    def test_properties(self):
        self.assertEqual(self.kpoint.frac_coords[0], 0.1)
        self.assertEqual(self.kpoint.frac_coords[1], 0.4)
        self.assertEqual(self.kpoint.frac_coords[2], -0.5)
        self.assertEqual(self.kpoint.a, 0.1)
        self.assertEqual(self.kpoint.b, 0.4)
        self.assertEqual(self.kpoint.c, -0.5)
        self.assertEqual(self.lattice, Lattice.cubic(10.0))
        self.assertEqual(self.kpoint.cart_coords[0], 1.0)
        self.assertEqual(self.kpoint.cart_coords[1], 4.0)
        self.assertEqual(self.kpoint.cart_coords[2], -5.0)
        self.assertEqual(self.kpoint.label, "X")


class BandStructureSymmLine_test(PymatgenTest):

    def setUp(self):
        self.bs = loadfn(os.path.join(test_dir, "Cu2O_361_bandstructure.json"))
        self.bs2 = loadfn(os.path.join(test_dir, "CaO_2605_bandstructure.json"))
        self.bs_spin = loadfn(os.path.join(test_dir, "NiO_19009_bandstructure.json"))
        self.bs_cbm0 = loadfn(os.path.join(test_dir, "InN_22205_bandstructure.json"))
        self.bs_cu = loadfn(os.path.join(test_dir, "Cu_30_bandstructure.json"))
        self.bs_diff_spins = loadfn(os.path.join(test_dir, "VBr2_971787_bandstructure.json"))
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_basic(self):
        self.assertArrayAlmostEqual(self.bs.projections[Spin.up][10][12][0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        self.assertArrayAlmostEqual(self.bs.projections[Spin.up][25][0][
                                 Orbital.dyz.value],
                             [0.0, 0.0, 0.0011, 0.0219, 0.0219, 0.069])
        self.assertAlmostEqual(
            self.bs.get_projection_on_elements()[Spin.up][25][10]['O'], 0.0328)
        self.assertAlmostEqual(
            self.bs.get_projection_on_elements()[Spin.up][22][25]['Cu'], 0.8327)
        proj = self.bs.get_projections_on_elements_and_orbitals({'Cu': ['s',
                                                                        'd']})
        self.assertAlmostEqual(
            proj[Spin.up][25][0]['Cu']['s'], 0.0027)
        self.assertAlmostEqual(
            proj[Spin.up][25][0]['Cu']['d'], 0.8495999999999999)

        self.assertEqual(self.bs2.nb_bands, 16)
        self.assertAlmostEqual(self.bs2.bands[Spin.up][5][10], 0.5608)
        self.assertAlmostEqual(self.bs2.bands[Spin.up][5][10], 0.5608)
        self.assertEqual(self.bs2.branches[5]['name'], "L-U")
        self.assertEqual(self.bs2.branches[5]['start_index'], 80)
        self.assertEqual(self.bs2.branches[5]['end_index'], 95)
        self.assertAlmostEqual(self.bs2.distance[70], 4.2335127528765737)

        self.assertEqual(self.bs_spin.nb_bands, 27)
        self.assertAlmostEqual(self.bs_spin.bands[Spin.up][5][10], 0.262)
        self.assertAlmostEqual(self.bs_spin.bands[Spin.down][5][10],
                               1.6156)

    def test_properties(self):
        self.one_kpoint = self.bs2.kpoints[31]
        self.assertEqual(self.one_kpoint.frac_coords[0], 0.5)
        self.assertEqual(self.one_kpoint.frac_coords[1], 0.25)
        self.assertEqual(self.one_kpoint.frac_coords[2], 0.75)
        self.assertAlmostEqual(self.one_kpoint.cart_coords[0], 0.64918757)
        self.assertAlmostEqual(self.one_kpoint.cart_coords[1], 1.29837513)
        self.assertAlmostEqual(self.one_kpoint.cart_coords[2], 0.0)
        self.assertEqual(self.one_kpoint.label, "W")

        self.assertAlmostEqual(self.bs2.efermi, 2.6211967, "wrong fermi energy")

    def test_get_branch(self):
        self.assertAlmostEqual(self.bs2.get_branch(110)[0]['name'], "U-W")

    def test_get_direct_band_gap_dict(self):
        direct_dict = self.bs_diff_spins.get_direct_band_gap_dict()
        self.assertEqual(direct_dict[Spin.down]['value'], 4.5365)

        for bs in [self.bs2, self.bs_spin]:
            dg_dict = bs.get_direct_band_gap_dict()
            for spin, v in bs.bands.items():
                kpt = dg_dict[spin]['kpoint_index']
                vb, cb = dg_dict[spin]['band_indices']
                gap = v[cb][kpt] - v[vb][kpt]
                self.assertEqual(gap, dg_dict[spin]['value'])
        self.assertRaises(ValueError, self.bs_cu.get_direct_band_gap_dict)

    def test_get_direct_band_gap(self):
        self.assertAlmostEqual(self.bs2.get_direct_band_gap(),
                               4.0125999999999999)
        self.assertTrue(self.bs_diff_spins.get_direct_band_gap() > 0)
        self.assertEqual(self.bs_cu.get_direct_band_gap(), 0)

    def test_is_metal(self):
        self.assertFalse(self.bs2.is_metal(), "wrong metal assignment")
        self.assertFalse(self.bs_spin.is_metal(), "wrong metal assignment")
        self.assertTrue(self.bs_cu.is_metal(), "wrong metal assignment")

    def test_get_cbm(self):
        cbm = self.bs2.get_cbm()
        self.assertAlmostEqual(cbm['energy'], 5.8709, "wrong CBM energy")
        self.assertEqual(cbm['band_index'][Spin.up][0], 8, "wrong CBM band index")
        self.assertEqual(cbm['kpoint_index'][0], 15, "wrong CBM kpoint index")
        self.assertEqual(cbm['kpoint'].frac_coords[0], 0.5, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm['kpoint'].frac_coords[1], 0.0, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm['kpoint'].frac_coords[2], 0.5, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm['kpoint'].label, "X", "wrong CBM kpoint label")
        cbm_spin = self.bs_spin.get_cbm()
        self.assertAlmostEqual(cbm_spin['energy'], 8.0458, "wrong CBM energy")
        self.assertEqual(cbm_spin['band_index'][Spin.up][0], 12, "wrong CBM band index")
        self.assertEqual(len(cbm_spin['band_index'][Spin.down]), 0, "wrong CBM band index")
        self.assertEqual(cbm_spin['kpoint_index'][0], 0, "wrong CBM kpoint index")
        self.assertEqual(cbm_spin['kpoint'].frac_coords[0], 0.0, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm_spin['kpoint'].frac_coords[1], 0.0, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm_spin['kpoint'].frac_coords[2], 0.0, "wrong CBM kpoint frac coords")
        self.assertEqual(cbm_spin['kpoint'].label, "\\Gamma", "wrong CBM kpoint label")

    def test_get_vbm(self):
        vbm = self.bs2.get_vbm()
        self.assertAlmostEqual(vbm['energy'], 2.2361, "wrong VBM energy")
        self.assertEqual(len(vbm['band_index'][Spin.up]), 3, "wrong VBM number of bands")
        self.assertEqual(vbm['band_index'][Spin.up][0], 5, "wrong VBM band index")
        self.assertEqual(vbm['kpoint_index'][0], 0, "wrong VBM kpoint index")
        self.assertEqual(vbm['kpoint'].frac_coords[0], 0.0, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm['kpoint'].frac_coords[1], 0.0, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm['kpoint'].frac_coords[2], 0.0, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm['kpoint'].label, "\\Gamma", "wrong VBM kpoint label")
        vbm_spin = self.bs_spin.get_vbm()
        self.assertAlmostEqual(vbm_spin['energy'], 5.731, "wrong VBM energy")
        self.assertEqual(len(vbm_spin['band_index'][Spin.up]), 2, "wrong VBM number of bands")
        self.assertEqual(len(vbm_spin['band_index'][Spin.down]), 0, "wrong VBM number of bands")
        self.assertEqual(vbm_spin['band_index'][Spin.up][0], 10, "wrong VBM band index")
        self.assertEqual(vbm_spin['kpoint_index'][0], 79, "wrong VBM kpoint index")
        self.assertEqual(vbm_spin['kpoint'].frac_coords[0], 0.5, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm_spin['kpoint'].frac_coords[1], 0.5, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm_spin['kpoint'].frac_coords[2], 0.5, "wrong VBM kpoint frac coords")
        self.assertEqual(vbm_spin['kpoint'].label, "L", "wrong VBM kpoint label")

    def test_get_band_gap(self):
        bg = self.bs2.get_band_gap()
        self.assertAlmostEqual(bg['energy'], 3.6348, "wrong gap energy")
        self.assertEqual(bg['transition'], "\\Gamma-X", "wrong kpoint transition")
        self.assertFalse(bg['direct'], "wrong nature of the gap")
        bg_spin = self.bs_spin.get_band_gap()
        self.assertAlmostEqual(bg_spin['energy'], 2.3148, "wrong gap energy")
        self.assertEqual(bg_spin['transition'], "L-\\Gamma", "wrong kpoint transition")
        self.assertFalse(bg_spin['direct'], "wrong nature of the gap")
        bg_cbm0 = self.bs_cbm0.get_band_gap()
        self.assertAlmostEqual(bg_cbm0['energy'], 0, places=3, msg="wrong gap energy")

    def test_get_sym_eq_kpoints_and_degeneracy(self):
        bs = self.bs2
        cbm_k = bs.get_cbm()['kpoint'].frac_coords
        vbm_k = bs.get_vbm()['kpoint'].frac_coords
        self.assertEqual(bs.get_kpoint_degeneracy(cbm_k), None)
        bs.structure = loadfn(os.path.join(test_dir, "CaO_2605_structure.json"))
        self.assertEqual(bs.get_kpoint_degeneracy(cbm_k), 3)
        self.assertEqual(bs.get_kpoint_degeneracy(vbm_k), 1)
        cbm_eqs = bs.get_sym_eq_kpoints(cbm_k)
        self.assertTrue([0.5, 0., 0.5] in cbm_eqs)
        self.assertTrue([0., 0.5, 0.5] in cbm_eqs)
        self.assertTrue([0.5, 0.5, 0.] in cbm_eqs)
        vbm_eqs = bs.get_sym_eq_kpoints(vbm_k)
        self.assertTrue([0., 0., 0.] in vbm_eqs)

    def test_as_dict(self):
        s = json.dumps(self.bs.as_dict())
        self.assertIsNotNone(s)
        s = json.dumps(self.bs2.as_dict())
        self.assertIsNotNone(s)
        s = json.dumps(self.bs_spin.as_dict())
        self.assertIsNotNone(s)

    def test_old_format_load(self):
        with open(os.path.join(test_dir, "bs_ZnS_old.json"),
                  "r", encoding='utf-8') as f:
            d = json.load(f)
            bs_old = BandStructureSymmLine.from_dict(d)
            self.assertEqual(bs_old.get_projection_on_elements()[
                                 Spin.up][0][0]['Zn'], 0.0971)


class ReconstructBandStructureTest(PymatgenTest):

    def setUp(self):
        self.bs_cu = loadfn(os.path.join(test_dir, "Cu_30_bandstructure.json"))
        self.bs_cu2 = loadfn(os.path.join(test_dir, "Cu_30_bandstructure.json"))
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_reconstruct_band_structure(self):
        bs = get_reconstructed_band_structure([self.bs_cu, self.bs_cu2])
        self.assertEqual(bs.bands[Spin.up].shape, (20, 700), "wrong number of bands or kpoints")

    def test_vasprun_bs(self):

        bsv = BSVasprun(os.path.join(test_dir, "vasprun.xml"),
                        parse_projected_eigen=True,
                        parse_potcar_file=True)
        bs = bsv.get_band_structure(kpoints_filename=os.path.join(test_dir, "KPOINTS.band"),
                                    line_mode=True)
        bs.get_projection_on_elements()


if __name__ == '__main__':
    unittest.main()
