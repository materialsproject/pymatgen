# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import os
import shutil

import numpy as np
from monty.json import MontyDecoder

from pymatgen.io.vasp.sets import MITVaspInputSet, MITHSEVaspInputSet, \
    MPVaspInputSet, MITGGAVaspInputSet, MITNEBVaspInputSet,\
    MPStaticVaspInputSet, MPNonSCFVaspInputSet, MITMDVaspInputSet,\
    MPHSEVaspInputSet, MPBSHSEVaspInputSet, MPStaticDielectricDFPTVaspInputSet,\
    MPOpticsNonSCFVaspInputSet
from pymatgen.io.vasp.inputs import Poscar, Incar
from pymatgen import Specie, Lattice, Structure

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
        self.assertEqual(kpoints.style, 'Monkhorst')

        kpoints = self.mitparamset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, 'Monkhorst')

        kpoints = self.mpstaticparamset.get_kpoints(self.struct)
        self.assertEqual(kpoints.kpts, [[6, 6, 4]])
        self.assertEqual(kpoints.style, 'Monkhorst')

        kpoints = self.mpnscfparamsetl.get_kpoints(self.struct)
        self.assertEqual(kpoints.num_kpts, 140)
        self.assertEqual(kpoints.style, 'Reciprocal')

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
        self.assertEqual(kpoints.style, 'Gamma')

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
        self.assertEqual(kpoints.style, 'Monkhorst')

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
            structs.append(Structure.from_sites(s.sites,
                                        to_unit_cell=True))

        fc = self.vis._process_structures(structs)[2].frac_coords
        self.assertTrue(np.allclose(fc, [[0.5]*3,[0.9, 1.033333, 1.0333333]]))


if __name__ == '__main__':
    unittest.main()
