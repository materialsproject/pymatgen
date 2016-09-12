# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest2 as unittest
import tempfile
from monty.tempfile import ScratchDir
from monty.json import MontyDecoder

from pymatgen.io.vasp.sets import *
from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints
from pymatgen import Specie, Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.outputs import Vasprun

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

dec = MontyDecoder()


class MITMPRelaxSetTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        if "VASP_PSP_DIR" not in os.environ:
            os.environ["VASP_PSP_DIR"] = test_dir
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.structure = poscar.structure
        self.coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        self.lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])

        self.mitset = MITRelaxSet(self.structure)
        self.mitset_unsorted = MITRelaxSet(self.structure, sort_structure=False)
        self.mpset = MPRelaxSet(self.structure)

    def test_poscar(self):

        structure = Structure(self.lattice, ["Fe", "Mn"], self.coords)
        mitparamset = MITRelaxSet(structure, sort_structure=False)
        s_unsorted = mitparamset.poscar.structure
        mitparamset = MITRelaxSet(structure, sort_structure=True)
        s_sorted = mitparamset.poscar.structure
        self.assertEqual(s_unsorted[0].specie.symbol, 'Fe')
        self.assertEqual(s_sorted[0].specie.symbol, 'Mn')

    def test_potcar_symbols(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.25, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        structure = Structure(lattice, ["P", "Fe", "O"], coords)
        mitparamset = MITRelaxSet(structure)
        syms = mitparamset.potcar_symbols
        self.assertEqual(syms, ['Fe', 'P', 'O'])
        paramset = MPRelaxSet(structure, sort_structure=False)
        syms = paramset.potcar_symbols
        self.assertEqual(syms, ['P', 'Fe_pv', 'O'])

    def test_lda_potcar(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        p = MITRelaxSet(structure, potcar_functional="LDA").potcar
        self.assertEqual(p.functional, 'LDA')

    def test_nelect(self):
        coords = [[0]*3, [0.5]*3, [0.75]*3]
        lattice = Lattice.cubic(4)
        s = Structure(lattice, ['Si', 'Si', 'Fe'], coords)
        self.assertAlmostEqual(MITRelaxSet(s).nelect, 16)

    def test_get_incar(self):

        incar = self.mpset.incar

        self.assertEqual(incar['LDAUU'], [5.3, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)

        incar = self.mitset.incar
        self.assertEqual(incar['LDAUU'], [4.0, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)

        si = 14
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        #Silicon structure for testing.
        latt = Lattice(np.array([[3.8401979337, 0.00, 0.00],
                              [1.9200989668, 3.3257101909, 0.00],
                              [0.00, -2.2171384943, 3.1355090603]]))
        struct = Structure(latt, [si, si], coords)
        incar = MPRelaxSet(struct).incar
        self.assertNotIn("LDAU", incar)

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice([[3.8401979337, 0.00, 0.00],
                           [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        incar = MPRelaxSet(struct).incar
        self.assertNotIn('LDAU', incar)

        #check fluorides
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar['LDAUU'], [5.3, 0])
        self.assertEqual(incar['MAGMOM'], [5, 0.6])

        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar['LDAUU'], [4.0, 0])

        #Make sure this works with species.
        struct = Structure(lattice, ["Fe2+", "O2-"], coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar['LDAUU'], [5.3, 0])

        struct = Structure(lattice, ["Fe", "Mn"], coords,
                           site_properties={'magmom': (5.2, -4.5)})
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar['MAGMOM'], [-4.5, 5.2])

        incar = MITRelaxSet(struct, sort_structure=False).incar
        self.assertEqual(incar['MAGMOM'], [5.2, -4.5])

        struct = Structure(lattice, [Specie("Fe", 2, {'spin': 4.1}), "Mn"],
                           coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar['MAGMOM'], [5, 4.1])


        struct = Structure(lattice, ["Mn3+", "Mn4+"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar['MAGMOM'], [4, 3])

        userset = MPRelaxSet(struct,
            user_incar_settings={'MAGMOM': {"Fe": 10, "S": -5, "Mn3+": 100}}
        )
        self.assertEqual(userset.incar['MAGMOM'], [100, 0.6])

        #sulfide vs sulfate test

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.25, 0.5, 0])

        struct = Structure(lattice, ["Fe", "Fe", "S"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar['LDAUU'], [1.9, 0])

        #Make sure Matproject sulfides are ok.
        self.assertNotIn('LDAUU', MPRelaxSet(struct).incar)

        struct = Structure(lattice, ["Fe", "S", "O"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar['LDAUU'], [4.0, 0, 0])

        #Make sure Matproject sulfates are ok.
        self.assertEqual(MPRelaxSet(struct).incar['LDAUU'], [5.3, 0, 0])

    def test_get_kpoints(self):
        kpoints = MPRelaxSet(self.structure).kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

        kpoints = MPRelaxSet(self.structure, user_kpoints_settings={
            "reciprocal_density": 1000}).kpoints
        self.assertEqual(kpoints.kpts, [[6, 11, 13]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

        kpoints = self.mitset.kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

        recip_paramset = MPRelaxSet(self.structure, force_gamma=True)
        recip_paramset.kpoints_settings = {"reciprocal_density": 40}
        kpoints = recip_paramset.kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_all_input(self):
        d = self.mitset.all_input
        self.assertEqual(d["INCAR"]["ISMEAR"], -5)
        s = self.structure.copy()
        s.make_supercell(4)
        paramset = MITRelaxSet(s)
        d = paramset.all_input
        self.assertEqual(d["INCAR"]["ISMEAR"], 0)

    def test_as_from_dict(self):
        mitset = MITRelaxSet(self.structure)
        mpset = MPRelaxSet(self.structure)
        mpuserset = MPRelaxSet(self.structure,
            user_incar_settings={'MAGMOM': {"Fe": 10, "S": -5, "Mn3+": 100}}
        )

        d = mitset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.config_dict["INCAR"]["LDAUU"]["O"]["Fe"], 4)

        d = mpset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.config_dict["INCAR"]["LDAUU"]["O"]["Fe"], 5.3)

        d = mpuserset.as_dict()
        v = dec.process_decoded(d)
        #self.assertEqual(type(v), MPVaspInputSet)
        self.assertEqual(v.user_incar_settings["MAGMOM"],
                         {"Fe": 10, "S": -5, "Mn3+": 100})

    def test_hubbard_off_and_ediff_override(self):
        p = MPRelaxSet(self.structure, user_incar_settings={"LDAU": False,
                                                            "EDIFF": 1e-10})
        self.assertNotIn("LDAUU", p.incar)
        self.assertEqual(p.incar["EDIFF"], 1e-10)

    def test_write_input(self):
        with ScratchDir(".") as d:
            self.mitset.write_input(d, make_dir_if_not_present=True)
            for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
                self.assertTrue(os.path.exists(f))
            self.assertFalse(os.path.exists("Fe4P4O16.cif"))
            self.mitset.write_input(d, make_dir_if_not_present=True,
                                    include_cif=True)
            self.assertTrue(os.path.exists("Fe4P4O16.cif"))


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
        self.assertEqual(non_prev_vis.kpoints.kpts, [[13, 11, 11]])
        non_prev_vis = MPStaticSet(vis.structure, reciprocal_density=200)
        self.assertEqual(non_prev_vis.kpoints.kpts, [[15, 13, 13]])

    def tearDown(self):
        shutil.rmtree(self.tmp)


class MPNonSCFSetTest(PymatgenTest):

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


class MITMDSetTest(unittest.TestCase):

    def setUp(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.structure
        self.mitmdparam = MITMDSet(self.struct, 300, 1200, 10000)

    def test_params(self):
        param = self.mitmdparam
        syms = param.potcar_symbols
        self.assertEqual(syms, ['Fe', 'P', 'O'])
        incar = param.incar
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar['EDIFF'], 2.4e-5)
        kpoints = param.kpoints
        self.assertEqual(kpoints.kpts, [(1, 1, 1)])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_as_from_dict(self):
        d = self.mitmdparam.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MITMDSet)
        self.assertEqual(v.config_dict["INCAR"]["TEBEG"], 300)


class MITNEBSetTest(unittest.TestCase):

    def setUp(self):
        c1 = [[0.5] * 3, [0.9] * 3]
        c2 = [[0.5] * 3, [0.9, 0.1, 0.1]]
        s1 = Structure(Lattice.cubic(5), ['Si', 'Si'], c1)
        s2 = Structure(Lattice.cubic(5), ['Si', 'Si'], c2)
        structs = []
        for s in s1.interpolate(s2, 3, pbc=True):
            structs.append(Structure.from_sites(s.sites, to_unit_cell=True))
        self.structures = structs
        self.vis = MITNEBSet(self.structures)

    def test_potcar_symbols(self):
        syms = self.vis.potcar_symbols
        self.assertEqual(syms, ['Si'])

    def test_incar(self):
        incar = self.vis.incar
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar['EDIFF'], 0.00005)

    def test_kpoints(self):
        kpoints = self.vis.kpoints
        self.assertEqual(kpoints.kpts, [[8, 8, 8]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)

    def test_as_from_dict(self):
        d = self.vis.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v.config_dict["INCAR"]["IMAGES"], 2)

    def test_write_input(self):
        with ScratchDir(".") as d:
            self.vis.write_input(d, write_cif=True,
                                 write_endpoint_inputs=True,
                                 write_path_cif=True)
            self.assertTrue(os.path.exists("INCAR"))
            self.assertTrue(os.path.exists("KPOINTS"))
            self.assertTrue(os.path.exists("POTCAR"))
            self.assertTrue(os.path.exists("00/POSCAR"))
            self.assertTrue(os.path.exists("01/POSCAR"))
            self.assertTrue(os.path.exists("02/POSCAR"))
            self.assertTrue(os.path.exists("03/POSCAR"))
            self.assertFalse(os.path.exists("04/POSCAR"))
            self.assertTrue(os.path.exists("00/INCAR"))
            self.assertTrue(os.path.exists("path.cif"))


class MPSOCSetTest(PymatgenTest):

    def test_from_prev_calc(self):
        prev_run = os.path.join(test_dir, "fe_monomer")
        vis = MPSOCSet.from_prev_calc(prev_calc_dir=prev_run, magmom=[3],
                                      saxis=(1, 0, 0))
        self.assertEqual(vis.incar["ISYM"], -1)
        self.assertTrue(vis.incar["LSORBIT"])
        self.assertEqual(vis.incar["ICHARG"], 11)
        self.assertEqual(vis.incar["SAXIS"], [1, 0, 0])
        self.assertEqual(vis.incar["MAGMOM"], [[0, 0, 3]])


class MVLSlabSetTest(PymatgenTest):

    def setUp(self):

        if "VASP_PSP_DIR" not in os.environ:
            os.environ["VASP_PSP_DIR"] = test_dir
        s = PymatgenTest.get_structure("Li2O")
        gen = SlabGenerator(s, (1, 0, 0), 10, 10)
        self.slab = gen.get_slab()
        self.bulk = self.slab.oriented_unit_cell

        vis_bulk = MVLSlabSet(self.bulk, bulk=True)
        vis = MVLSlabSet(self.slab)
        vis_bulk_gpu = MVLSlabSet(self.bulk, bulk=True, gpu=True)

        self.d_bulk = vis_bulk.all_input
        self.d_slab = vis.all_input
        self.d_bulk_gpu = vis_bulk_gpu.all_input

    def test_bulk(self):

        incar_bulk = self.d_bulk["INCAR"]
        incar_bulk_gpu = self.d_bulk_gpu["INCAR"]
        poscar_bulk = self.d_bulk["POSCAR"]

        self.assertEqual(incar_bulk["ISIF"], 3)
        self.assertEqual(poscar_bulk.structure.formula,
                         self.bulk.formula)
        # Test VASP-gpu compatibility
        self.assertEqual(incar_bulk_gpu["KPAR"], 1)
        self.assertTrue("NPAR" not in incar_bulk_gpu.keys())

    def test_slab(self):

        incar_slab = self.d_slab["INCAR"]
        poscar_slab = self.d_slab["POSCAR"]
        potcar_slab = self.d_slab["POTCAR"]

        self.assertEqual(incar_slab["AMIN"], 0.01)
        self.assertEqual(incar_slab["AMIX"], 0.2)
        self.assertEqual(incar_slab["BMIX"], 0.001)
        self.assertEqual(incar_slab["NELMIN"], 8)
        # No volume relaxation during slab calculations
        self.assertEqual(incar_slab["ISIF"], 2)
        self.assertEqual(potcar_slab.functional, 'PBE')
        self.assertEqual(potcar_slab.symbols[0], u'Li_sv')
        self.assertEqual(potcar_slab.symbols[1], u'O')
        self.assertEqual(poscar_slab.structure.formula,
                         self.slab.formula)

    def test_kpoints(self):

        kpoints_slab = self.d_slab["KPOINTS"].kpts[0]
        kpoints_bulk = self.d_bulk["KPOINTS"].kpts[0]

        self.assertEqual(kpoints_bulk[0], kpoints_slab[0])
        self.assertEqual(kpoints_bulk[1], kpoints_slab[1])
        self.assertEqual(kpoints_bulk[0], 15)
        self.assertEqual(kpoints_bulk[1], 15)
        self.assertEqual(kpoints_bulk[2], 15)
        # The last kpoint in a slab should always be 1
        self.assertEqual(kpoints_slab[2], 1)


class MVLElasticSetTest(PymatgenTest):

    def test_incar(self):
        mvlparam = MVLElasticSet(self.get_structure("Graphite"))
        incar = mvlparam.incar
        self.assertEqual(incar["IBRION"], 6)
        self.assertEqual(incar["NFREE"], 2)
        self.assertEqual(incar["POTIM"], 0.015)
        self.assertNotIn("NPAR", incar)


class MPHSEBSTest(PymatgenTest):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_init(self):
        prev_run = os.path.join(test_dir, "static_silicon")
        vis = MPHSEBSSet.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform")
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(len(vis.kpoints.kpts), 31)

        vis = MPHSEBSSet.from_prev_calc(prev_calc_dir=prev_run, mode="Line")
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(vis.incar['HFSCREEN'], 0.2)
        self.assertEqual(vis.incar['NSW'], 0)
        self.assertEqual(len(vis.kpoints.kpts), 195)


class FuncTest(PymatgenTest):

    def test_batch_write_input(self):
        with ScratchDir("."):
            structures = [PymatgenTest.get_structure("Li2O"),
                          PymatgenTest.get_structure("LiFePO4")]
            batch_write_input(structures)
            for d in ['Li4Fe4P4O16_1', 'Li2O1_0']:
                for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
                    self.assertTrue(os.path.exists(os.path.join(d, f)))

if __name__ == '__main__':
    unittest.main()
