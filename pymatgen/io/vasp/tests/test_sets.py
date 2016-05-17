# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import os
import shutil
import tempfile

import numpy as np
from monty.json import MontyDecoder

from pymatgen.io.vasp.sets import *
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints
from pymatgen import Specie, Lattice, Structure
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.outputs import Vasprun

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

dec = MontyDecoder()


class MITMPVaspInputSetTest(unittest.TestCase):

    def setUp(self):
        if "VASP_PSP_DIR" not in os.environ:
            os.environ["VASP_PSP_DIR"] = test_dir
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.structure

        self.mitparamset = MITVaspInputSet()
        self.mitparamset_unsorted = MITVaspInputSet(sort_structure=False)
        self.mithseparamset = MITHSEVaspInputSet()
        self.paramset = MPVaspInputSet()
        self.userparamset = MPVaspInputSet(
            user_incar_settings={'MAGMOM': {"Fe": 10, "S": -5, "Mn3+": 100}}
        )
        self.mitggaparam = MITGGAVaspInputSet()
        self.mpstaticparamset = MPStaticVaspInputSet()
        self.mpnscfparamsetu = MPNonSCFVaspInputSet(
            {"NBANDS": 50}, mode="Uniform")
        self.mpnscfparamsetl = MPNonSCFVaspInputSet(
            {"NBANDS": 60}, mode="Line")
        self.mphseparamset = MPHSEVaspInputSet()
        self.mpbshseparamsetl = MPBSHSEVaspInputSet(mode="Line")
        self.mpbshseparamsetu = MPBSHSEVaspInputSet(
            mode="Uniform", added_kpoints=[[0.5, 0.5, 0.0]])
        self.mpdielparamset = MPStaticDielectricDFPTVaspInputSet()

    def test_get_poscar(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        s_unsorted = self.mitparamset_unsorted.get_poscar(struct).structure
        s_sorted = self.mitparamset.get_poscar(struct).structure

        self.assertEqual(s_unsorted[0].specie.symbol, 'Fe')
        self.assertEqual(s_sorted[0].specie.symbol, 'Mn')

    def test_get_potcar_symbols(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.25, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["P", "Fe", "O"], coords)

        syms = self.paramset.get_potcar_symbols(struct)
        self.assertEqual(syms, ['Fe_pv', 'P', 'O'])

        syms = MPVaspInputSet(sort_structure=False).get_potcar_symbols(struct)
        self.assertEqual(syms, ['P', 'Fe_pv', 'O'])

    def test_false_potcar_hash(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.25, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["P", "Fe", "O"], coords)

        self.mitparamset.potcar_settings['Fe']['symbol'] = 'Fe_pv'
        self.assertRaises(ValueError, self.mitparamset.get_potcar, struct, check_hash=True)
        self.mitparamset.potcar_settings['Fe']['symbol'] = 'Fe'

    def test_lda_potcar(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["P", "Fe"], coords)
        p = MITVaspInputSet(potcar_functional="LDA").get_potcar(struct)
        self.assertEqual(p.functional, 'LDA')

    def test_get_nelect(self):
        coords = [[0]*3, [0.5]*3, [0.75]*3]
        lattice = Lattice.cubic(4)
        s = Structure(lattice, ['Si', 'Si', 'Fe'], coords)
        self.assertAlmostEqual(MITVaspInputSet().get_nelect(s), 16)

    def test_get_incar(self):
        incar = self.paramset.get_incar(self.struct)

        self.assertEqual(incar['LDAUU'], [5.3, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)

        incar = self.mitparamset.get_incar(self.struct)
        self.assertEqual(incar['LDAUU'], [4.0, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)

        incar_gga = self.mitggaparam.get_incar(self.struct)
        self.assertNotIn("LDAU", incar_gga)

        incar_static = self.mpstaticparamset.get_incar(self.struct)
        self.assertEqual(incar_static["NSW"], 0)

        incar_nscfl = self.mpnscfparamsetl.get_incar(self.struct)
        self.assertEqual(incar_nscfl["NBANDS"], 60)

        incar_nscfu = self.mpnscfparamsetu.get_incar(self.struct)
        self.assertEqual(incar_nscfu["ISYM"], 0)

        incar_hse = self.mphseparamset.get_incar(self.struct)
        self.assertEqual(incar_hse['LHFCALC'], True)
        self.assertEqual(incar_hse['HFSCREEN'], 0.2)

        incar_hse_bsl = self.mpbshseparamsetl.get_incar(self.struct)
        self.assertEqual(incar_hse_bsl['LHFCALC'], True)
        self.assertEqual(incar_hse_bsl['HFSCREEN'], 0.2)
        self.assertEqual(incar_hse_bsl['NSW'], 0)

        incar_hse_bsu = self.mpbshseparamsetu.get_incar(self.struct)
        self.assertEqual(incar_hse_bsu['LHFCALC'], True)
        self.assertEqual(incar_hse_bsu['HFSCREEN'], 0.2)
        self.assertEqual(incar_hse_bsu['NSW'], 0)

        incar_diel = self.mpdielparamset.get_incar(self.struct)
        self.assertEqual(incar_diel['IBRION'], 8)
        self.assertEqual(incar_diel['LEPSILON'], True)

        si = 14
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        #Silicon structure for testing.
        latt = Lattice(np.array([[3.8401979337, 0.00, 0.00],
                              [1.9200989668, 3.3257101909, 0.00],
                              [0.00, -2.2171384943, 3.1355090603]]))
        struct = Structure(latt, [si, si], coords)
        incar = self.paramset.get_incar(struct)
        self.assertNotIn("LDAU", incar)

        incar = self.mithseparamset.get_incar(self.struct)
        self.assertTrue(incar['LHFCALC'])

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        incar = self.paramset.get_incar(struct)
        self.assertNotIn('LDAU', incar)

        #check fluorides
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = self.paramset.get_incar(struct)
        self.assertEqual(incar['LDAUU'], [5.3, 0])
        self.assertEqual(incar['MAGMOM'], [5, 0.6])

        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = self.mitparamset.get_incar(struct)
        self.assertEqual(incar['LDAUU'], [4.0, 0])

        #Make sure this works with species.
        struct = Structure(lattice, ["Fe2+", "O2-"], coords)
        incar = self.paramset.get_incar(struct)
        self.assertEqual(incar['LDAUU'], [5.3, 0])

        struct = Structure(lattice, ["Fe", "Mn"], coords,
                           site_properties={'magmom': (5.2, -4.5)})
        incar = self.paramset.get_incar(struct)
        self.assertEqual(incar['MAGMOM'], [-4.5, 5.2])
        incar = self.mpstaticparamset.get_incar(struct)
        self.assertEqual(incar['MAGMOM'], [-4.5, 5.2])
        incar = self.mitparamset_unsorted.get_incar(struct)
        self.assertEqual(incar['MAGMOM'], [5.2, -4.5])

        struct = Structure(lattice, [Specie("Fe", 2, {'spin': 4.1}), "Mn"],
                           coords)
        incar = self.paramset.get_incar(struct)
        self.assertEqual(incar['MAGMOM'], [5, 4.1])
        incar = self.mpnscfparamsetl.get_incar(struct)
        self.assertEqual(incar.get('MAGMOM', None), None)

        struct = Structure(lattice, ["Mn3+", "Mn4+"], coords)
        incar = self.mitparamset.get_incar(struct)
        self.assertEqual(incar['MAGMOM'], [4, 3])
        incar = self.mpnscfparamsetu.get_incar(struct)
        self.assertEqual(incar.get('MAGMOM', None), None)

        self.assertEqual(self.userparamset.get_incar(struct)['MAGMOM'],
                         [100, 0.6])

        #sulfide vs sulfate test

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.25, 0.5, 0])

        struct = Structure(lattice, ["Fe", "Fe", "S"], coords)
        incar = self.mitparamset.get_incar(struct)
        self.assertEqual(incar['LDAUU'], [1.9, 0])

        #Make sure Matproject sulfides are ok.
        self.assertNotIn('LDAUU', self.paramset.get_incar(struct))
        self.assertNotIn('LDAUU', self.mpstaticparamset.get_incar(struct))

        struct = Structure(lattice, ["Fe", "S", "O"], coords)
        incar = self.mitparamset.get_incar(struct)
        self.assertEqual(incar['LDAUU'], [4.0, 0, 0])

        #Make sure Matproject sulfates are ok.
        self.assertEqual(self.paramset.get_incar(struct)['LDAUU'], [5.3, 0, 0])
        self.assertEqual(self.mpnscfparamsetl.get_incar(struct)['LDAUU'],
                         [5.3, 0, 0])

        self.assertEqual(self.userparamset.get_incar(struct)['MAGMOM'],
                         [10, -5, 0.6])

    def test_optics(self):
        self.mpopticsparamset = MPOpticsNonSCFVaspInputSet.from_previous_vasp_run(
            '{}/static_silicon'.format(test_dir), output_dir='optics_test_dir',
            nedos=1145)
        self.assertTrue(os.path.exists('optics_test_dir/CHGCAR'))
        incar = Incar.from_file('optics_test_dir/INCAR')
        self.assertTrue(incar['LOPTICS'])
        self.assertEqual(incar['NEDOS'], 1145)

        #Remove the directory in which the inputs have been created
        shutil.rmtree('optics_test_dir')

    def test_get_kpoints(self):
        kpoints = self.paramset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

        kpoints = self.mitparamset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

        kpoints = self.mpstaticparamset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[6, 6, 4]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

        kpoints = self.mpnscfparamsetl.get_kpoints(self.struct)
        self.assertEqual(kpoints.num_kpts, 140)
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Reciprocal)

        kpoints = self.mpnscfparamsetu.get_kpoints(self.struct)
        self.assertEqual(kpoints.num_kpts, 168)

        kpoints = self.mpbshseparamsetl.get_kpoints(self.struct)
        self.assertAlmostEqual(kpoints.num_kpts, 164)
        self.assertAlmostEqual(kpoints.kpts[10][0], 0.0)
        self.assertAlmostEqual(kpoints.kpts[10][1], 0.5)
        self.assertAlmostEqual(kpoints.kpts[10][2], 0.16666667)
        self.assertAlmostEqual(kpoints.kpts[26][0], 0.0714285714286)
        self.assertAlmostEqual(kpoints.kpts[26][1], 0.0)
        self.assertAlmostEqual(kpoints.kpts[26][2], 0.0)
        self.assertAlmostEqual(kpoints.kpts[-1][0], 0.5)
        self.assertAlmostEqual(kpoints.kpts[-1][1], 0.5)
        self.assertAlmostEqual(kpoints.kpts[-1][2], 0.5)

        kpoints = self.mpbshseparamsetu.get_kpoints(self.struct)
        self.assertAlmostEqual(kpoints.num_kpts, 25)
        self.assertAlmostEqual(kpoints.kpts[10][0], 0.0)
        self.assertAlmostEqual(kpoints.kpts[10][1], 0.5)
        self.assertAlmostEqual(kpoints.kpts[10][2], 0.16666667)
        self.assertAlmostEqual(kpoints.kpts[-1][0], 0.5)
        self.assertAlmostEqual(kpoints.kpts[-1][1], 0.5)
        self.assertAlmostEqual(kpoints.kpts[-1][2], 0.0)

        recip_paramset = MPVaspInputSet(force_gamma=True)
        recip_paramset.kpoints_settings = {"reciprocal_density": 40}
        kpoints = recip_paramset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_get_all_vasp_input(self):
        d = self.mitparamset.get_all_vasp_input(self.struct)
        self.assertEqual(d["INCAR"]["ISMEAR"], -5)
        self.struct.make_supercell(4)
        d = self.mitparamset.get_all_vasp_input(self.struct)
        self.assertEqual(d["INCAR"]["ISMEAR"], 0)

    def test_to_from_dict(self):
        self.mitparamset = MITVaspInputSet()
        self.mithseparamset = MITHSEVaspInputSet()
        self.paramset = MPVaspInputSet()
        self.userparamset = MPVaspInputSet(
            user_incar_settings={'MAGMOM': {"Fe": 10, "S": -5, "Mn3+": 100}}
        )

        d = self.mitparamset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.incar_settings["LDAUU"]["O"]["Fe"], 4)

        d = self.mitggaparam.as_dict()
        v = dec.process_decoded(d)
        self.assertNotIn("LDAUU", v.incar_settings)

        d = self.mithseparamset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.incar_settings["LHFCALC"], True)

        d = self.mphseparamset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.incar_settings["LHFCALC"], True)

        d = self.paramset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.incar_settings["LDAUU"]["O"]["Fe"], 5.3)

        d = self.userparamset.as_dict()
        v = dec.process_decoded(d)
        #self.assertEqual(type(v), MPVaspInputSet)
        self.assertEqual(v.incar_settings["MAGMOM"],
                         {"Fe": 10, "S": -5, "Mn3+": 100})


class MITMDVaspInputSetTest(unittest.TestCase):

    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.structure
        self.mitmdparam = MITMDVaspInputSet(300, 1200, 10000)

    def test_get_potcar_symbols(self):
        syms = self.mitmdparam.get_potcar_symbols(self.struct)
        self.assertEqual(syms, ['Fe', 'P', 'O'])

    def test_get_incar(self):
        incar = self.mitmdparam.get_incar(self.struct)
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar['EDIFF'], 2.4e-5)

    def test_get_kpoints(self):
        kpoints = self.mitmdparam.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [(1, 1, 1)])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_to_from_dict(self):
        d = self.mitmdparam.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MITMDVaspInputSet)
        self.assertEqual(v.incar_settings["TEBEG"], 300)


class MITNEBVaspInputSetTest(unittest.TestCase):

    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.structure
        self.vis = MITNEBVaspInputSet(nimages=10, hubbard_off=True)

    def test_get_potcar_symbols(self):
        syms = self.vis.get_potcar_symbols(self.struct)
        self.assertEqual(syms, ['Fe', 'P', 'O'])

    def test_get_incar(self):
        incar = self.vis.get_incar(self.struct)
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar['EDIFF'], 0.00005)

    def test_get_kpoints(self):
        kpoints = self.vis.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

    def test_to_from_dict(self):
        d = self.vis.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.incar_settings["IMAGES"], 10)

    def test_write_inputs(self):
        c1 = [[0.5] * 3, [0.9] * 3]
        c2 = [[0.5] * 3, [0.9, 0.1, 0.1]]
        s1 = Structure(Lattice.cubic(5), ['Si', 'Si'], c1)
        s2 = Structure(Lattice.cubic(5), ['Si', 'Si'], c2)
        structs = []
        for s in s1.interpolate(s2, 3, pbc=True):
            structs.append(Structure.from_sites(s.sites, to_unit_cell=True))

        fc = self.vis._process_structures(structs)[2].frac_coords
        self.assertTrue(np.allclose(fc, [[0.5]*3,[0.9, 1.033333, 1.0333333]]))


class MVLVaspInputSetTest(PymatgenTest):

    def setUp(self):
        self.mvlparam = MVLElasticInputSet()

    def test_get_incar(self):
        incar = self.mvlparam.get_incar(self.get_structure("Graphite"))
        self.assertEqual(incar["IBRION"], 6)
        self.assertEqual(incar["NFREE"], 2)
        self.assertEqual(incar["POTIM"], 0.015)
        self.assertNotIn("NPAR", incar)


class MPStaticSetTest(PymatgenTest):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_init(self):
        prev_run = os.path.join(test_dir, "relaxation")

        vis = MPStaticSet.from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Monkhorst)

        # Check as from dict.
        vis = MPStaticSet.from_dict(vis.as_dict())
        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Monkhorst)

        non_prev_vis = MPStaticSet(vis.structure,
                                   user_incar_settings={"LORBIT": 12, "LWAVE": True})
        self.assertEqual(non_prev_vis.incar["NSW"], 0)
        # Check that the ENCUT and Kpoints style has NOT been inherited.
        self.assertEqual(non_prev_vis.incar["ENCUT"], 520)
        # Check that user incar settings are applied.
        self.assertEqual(non_prev_vis.incar["LORBIT"], 12)
        self.assertTrue(non_prev_vis.incar["LWAVE"])

        self.assertEqual(non_prev_vis.kpoints.style,
                         Kpoints.supported_modes.Gamma)
        v2 = MPStaticSet.from_dict(non_prev_vis.as_dict())
        self.assertEqual(v2.incar["ENCUT"], 520)
        # Check that user incar settings are applied.
        self.assertEqual(v2.incar["LORBIT"], 12)
        leps_vis = MPStaticSet.from_prev_calc(prev_calc_dir=prev_run,
                                              lepsilon=True)
        self.assertTrue(leps_vis.incar["LEPSILON"])
        self.assertEqual(leps_vis.incar["IBRION"], 8)
        self.assertNotIn("NPAR", leps_vis.incar)
        self.assertNotIn("NSW", leps_vis.incar)

    def tearDown(self):
        shutil.rmtree(self.tmp)


class MPNonSCFDerivedSetTest(PymatgenTest):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_init(self):
        prev_run = os.path.join(test_dir, "relaxation")
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, mode="Line", copy_chgcar=False)
        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

        # Check as from dict.
        vis = MPNonSCFSet.from_dict(vis.as_dict())
        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

        vis.write_input(self.tmp)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "CHGCAR")))

        vis = MPNonSCFSet.from_prev_calc(prev_calc_dir=prev_run,
                                                mode="Line", copy_chgcar=True)
        vis.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "CHGCAR")))

        # Code below is just to make sure that the parameters are the same
        # between the old MPStaticVaspInputSet and the new MPStaticSet.
        # TODO: Delete code below in future.
        MPNonSCFVaspInputSet.from_previous_vasp_run(
            previous_vasp_dir=prev_run, output_dir=self.tmp, mode="Line")

        incar = Incar.from_file(os.path.join(self.tmp, "INCAR"))

        for k, v1 in vis.incar.items():
            v2 = incar.get(k)
            try:
                v1 = v1.upper()
                v2 = v2.upper()
            except:
                # Convert strings to upper case for comparison. Ignore other
                # types.
                pass
            self.assertEqual(v1, v2, str(v1)+str(v2))
        kpoints = Kpoints.from_file(os.path.join(self.tmp, "KPOINTS"))
        self.assertEqual(kpoints.style, vis.kpoints.style)
        self.assertArrayAlmostEqual(kpoints.kpts, vis.kpoints.kpts)

    def test_optics(self):
        prev_run = os.path.join(test_dir, "relaxation")
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, copy_chgcar=False, optics=True,
            mode="Uniform", nedos=2001)

        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertTrue(vis.incar["LOPTICS"])
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

    def tearDown(self):
        shutil.rmtree(self.tmp)


class MagmomLdauTest(PymatgenTest):

    def test_structure_from_prev_run(self):
        vrun = Vasprun(os.path.join(test_dir, "vasprun.xml.magmom_ldau"))
        structure = vrun.final_structure
        poscar = Poscar(structure)
        structure_decorated = get_structure_from_prev_run(vrun, sym_prec=0)
        ldau_ans = {'LDAUU': [5.3, 0.0], 'LDAUL': [2, 0], 'LDAUJ': [0.0, 0.0]}
        magmom_ans = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
        ldau_dict = {}
        for key in ('LDAUU', 'LDAUJ', 'LDAUL'):
            if hasattr(structure_decorated[0], key.lower()):
                m = dict(
                    [(site.specie.symbol, getattr(site, key.lower()))
                     for site in structure_decorated])
                ldau_dict[key] = [m[sym] for sym in poscar.site_symbols]
        magmom = [site.magmom for site in structure_decorated]
        self.assertEqual(ldau_dict, ldau_ans)
        self.assertEqual(magmom, magmom_ans)


if __name__ == '__main__':
    unittest.main()
