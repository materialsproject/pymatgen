# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import hashlib
import tempfile
import unittest

import pytest  # type: ignore
from _pytest.monkeypatch import MonkeyPatch  # type: ignore
from monty.json import MontyDecoder

from pymatgen import SETTINGS
from pymatgen.core import Species, Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import *
from pymatgen.util.testing import PymatgenTest

MODULE_DIR = Path(__file__).resolve().parent

dec = MontyDecoder()


class SetChangeCheckTest(PymatgenTest):
    def test_sets_changed(self):
        # WARNING!
        # These tests will fail when you change an input set.
        # They are included as a sanity check: if you want to change
        # an input set, please make sure to notify the users for that set.
        # For sets starting with "MVL" this is @shyuep, for sets starting
        # with "MP" this is @shyuep and @mkhorton.
        os.chdir(MODULE_DIR / "..")
        input_sets = glob.glob("*.yaml")
        hashes = {}
        for input_set in input_sets:
            with open(input_set, "r") as f:
                hashes[input_set] = hashlib.sha1(f.read().encode("utf-8")).hexdigest()
        known_hashes = {
            'MVLGWSet.yaml': 'f4df9516cf7dd923b37281172c662a70fa32bebc',
            'MVLRelax52Set.yaml': 'eb538ffb45c0cd13f13df48afc1e71c44d2e34b2',
            'MPHSERelaxSet.yaml': '2bb969e64b57ff049077c8ec10e64f94c9c97f42',
            'VASPIncarBase.yaml': '19762515f8deefb970f2968fca48a0d67f7964d4',
            'MPSCANRelaxSet.yaml': '2604952d387f6531bfc37641ac3b1ffcce9f1bc1',
            'MPRelaxSet.yaml': '4ea97d776fbdc7e168036f73e9176012a56c0a45',
            'MITRelaxSet.yaml': '1a0970f8cad9417ec810f7ab349dc854eaa67010',
            'vdW_parameters.yaml': '66541f58b221c8966109156f4f651b2ca8aa76da'}

        self.assertDictEqual(
            hashes,
            known_hashes,
            "These tests will fail when you change an input set. "
            "They are included as a sanity check: if you want to "
            "change an input set, please make sure to notify the "
            "users for that set. "
            'For sets starting with "MVL" this is @shyuep, '
            'for sets starting with "MP" this is @shyuep and @mkhorton.',
        )


class MITMPRelaxSetTest(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.monkeypatch = MonkeyPatch()

        filepath = cls.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        cls.structure = poscar.structure
        cls.coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        cls.lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )

        cls.mitset = MITRelaxSet(cls.structure)
        cls.mitset_unsorted = MITRelaxSet(cls.structure, sort_structure=False)
        cls.mpset = MPRelaxSet(cls.structure)

    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_metal_check(self):
        structure = Structure.from_spacegroup(
            "Fm-3m", Lattice.cubic(3), ["Cu"], [[0, 0, 0]]
        )

        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            vis = MITRelaxSet(structure)
            incar = vis.incar
            # Verify some things
            self.assertIn("ISMEAR", str(w[-1].message))

    def test_poscar(self):
        structure = Structure(self.lattice, ["Fe", "Mn"], self.coords)
        mitparamset = MITRelaxSet(structure, sort_structure=False)
        s_unsorted = mitparamset.poscar.structure
        mitparamset = MITRelaxSet(structure, sort_structure=True)
        s_sorted = mitparamset.poscar.structure
        self.assertEqual(s_unsorted[0].specie.symbol, "Fe")
        self.assertEqual(s_sorted[0].specie.symbol, "Mn")

    def test_potcar_symbols(self):
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.75, 0.25, 0.75])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        structure = Structure(lattice, ["P", "Fe", "O"], coords)
        mitparamset = MITRelaxSet(structure)
        syms = mitparamset.potcar_symbols
        self.assertEqual(syms, ["Fe", "P", "O"])
        paramset = MPRelaxSet(structure, sort_structure=False)
        syms = paramset.potcar_symbols
        self.assertEqual(syms, ["P", "Fe_pv", "O"])

    def test_potcar_validation(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        # Use pytest's monkeypatch to temporarily point pymatgen to a directory
        # containing the wrong POTCARs (LDA potcars in a PBE directory)
        with self.monkeypatch.context() as m:
            m.setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(self.TEST_FILES_DIR / "wrong_potcars"))
            with pytest.warns(BadInputSetWarning, match="not known by pymatgen"):
                MITRelaxSet(structure).potcar

    def test_lda_potcar(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        p = MITRelaxSet(structure, potcar_functional="LDA").potcar
        self.assertEqual(p.functional, "LDA")

    def test_nelect(self):
        coords = [[0] * 3, [0.5] * 3, [0.75] * 3]
        lattice = Lattice.cubic(4)
        s = Structure(lattice, ["Si", "Si", "Fe"], coords)
        self.assertAlmostEqual(MITRelaxSet(s).nelect, 16)

        # Check that it works even when oxidation states are present. Was a bug
        # previously.
        s = Structure(lattice, ["Si4+", "Si4+", "Fe2+"], coords)
        self.assertAlmostEqual(MITRelaxSet(s).nelect, 16)
        self.assertAlmostEqual(MPRelaxSet(s).nelect, 22)

        # Check that it works for disordered structure. Was a bug previously
        s = Structure(lattice, ["Si4+", "Fe2+", "Si4+"], coords)
        self.assertAlmostEqual(MITRelaxSet(s).nelect, 16)
        self.assertAlmostEqual(MPRelaxSet(s).nelect, 22)

    def test_get_incar(self):

        incar = self.mpset.incar

        self.assertEqual(incar["LDAUU"], [5.3, 0, 0])
        self.assertAlmostEqual(incar["EDIFF"], 0.0012)

        incar = self.mitset.incar
        self.assertEqual(incar["LDAUU"], [4.0, 0, 0])
        self.assertAlmostEqual(incar["EDIFF"], 1e-5)

        si = 14
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        # Silicon structure for testing.
        latt = Lattice(
            np.array(
                [
                    [3.8401979337, 0.00, 0.00],
                    [1.9200989668, 3.3257101909, 0.00],
                    [0.00, -2.2171384943, 3.1355090603],
                ]
            )
        )
        struct = Structure(latt, [si, si], coords)
        incar = MPRelaxSet(struct).incar
        self.assertNotIn("LDAU", incar)

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        lattice = Lattice(
            [
                [3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603],
            ]
        )
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        incar = MPRelaxSet(struct).incar
        self.assertNotIn("LDAU", incar)

        # check fluorides
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar["LDAUU"], [5.3, 0])
        self.assertEqual(incar["MAGMOM"], [5, 0.6])

        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar["LDAUU"], [4.0, 0])

        # Make sure this works with species.
        struct = Structure(lattice, ["Fe2+", "O2-"], coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar["LDAUU"], [5.3, 0])

        struct = Structure(
            lattice, ["Fe", "Mn"], coords, site_properties={"magmom": (5.2, -4.5)}
        )
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar["MAGMOM"], [-4.5, 5.2])

        incar = MITRelaxSet(struct, sort_structure=False).incar
        self.assertEqual(incar["MAGMOM"], [5.2, -4.5])

        struct = Structure(lattice, [Species("Fe", 2, {"spin": 4.1}), "Mn"], coords)
        incar = MPRelaxSet(struct).incar
        self.assertEqual(incar["MAGMOM"], [5, 4.1])

        struct = Structure(lattice, ["Mn3+", "Mn4+"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar["MAGMOM"], [4, 3])

        userset = MPRelaxSet(
            struct, user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}}
        )
        self.assertEqual(userset.incar["MAGMOM"], [100, 0.6])

        noencutset = MPRelaxSet(struct, user_incar_settings={"ENCUT": None})
        self.assertNotIn("ENCUT", noencutset.incar)

        # sulfide vs sulfate test

        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])
        coords.append([0.25, 0.5, 0])

        struct = Structure(lattice, ["Fe", "Fe", "S"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar["LDAUU"], [1.9, 0])

        # Make sure Matproject sulfides are ok.
        self.assertNotIn("LDAUU", MPRelaxSet(struct).incar)

        struct = Structure(lattice, ["Fe", "S", "O"], coords)
        incar = MITRelaxSet(struct).incar
        self.assertEqual(incar["LDAUU"], [4.0, 0, 0])

        # Make sure Matproject sulfates are ok.
        self.assertEqual(MPRelaxSet(struct).incar["LDAUU"], [5.3, 0, 0])

        # test for default LDAUU value
        userset_ldauu_fallback = MPRelaxSet(
            struct, user_incar_settings={"LDAUU": {"Fe": 5.0, "S": 0}}
        )
        self.assertEqual(userset_ldauu_fallback.incar["LDAUU"], [5.0, 0, 0])

        # Expected to be oxide (O is the most electronegative atom)
        s = Structure(lattice, ["Fe", "O", "S"], coords)
        incar = MITRelaxSet(s).incar
        self.assertEqual(incar["LDAUU"], [4.0, 0, 0])

        # Expected to be chloride (Cl is the most electronegative atom)
        s = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = MITRelaxSet(s, user_incar_settings={"LDAU": True}).incar
        self.assertFalse("LDAUU" in incar)  # LDAU = False

        # User set a compound to be sulfide by specifing values of "LDAUL" etc.
        s = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = MITRelaxSet(
            s,
            user_incar_settings={
                "LDAU": True,
                "LDAUL": {"Fe": 3},
                "LDAUU": {"Fe": 1.8},
            },
        ).incar
        self.assertEqual(incar["LDAUL"], [3.0, 0, 0])
        self.assertEqual(incar["LDAUU"], [1.8, 0, 0])

        # test that van-der-Waals parameters are parsed correctly
        incar = MITRelaxSet(struct, vdw="optB86b").incar
        self.assertEqual(incar["GGA"], "Mk")
        self.assertEqual(incar["LUSE_VDW"], True)
        self.assertEqual(incar["PARAM1"], 0.1234)

        # Test that NELECT is updated when a charge is present
        si = 14
        coords = list()
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        # Silicon structure for testing.
        latt = Lattice(
            np.array(
                [
                    [3.8401979337, 0.00, 0.00],
                    [1.9200989668, 3.3257101909, 0.00],
                    [0.00, -2.2171384943, 3.1355090603],
                ]
            )
        )
        struct = Structure(latt, [si, si], coords, charge=1)
        mpr = MPRelaxSet(struct, use_structure_charge=True)
        self.assertEqual(
            mpr.incar["NELECT"], 7, "NELECT not properly set for nonzero charge"
        )

        # test that NELECT does not get set when use_structure_charge = False
        mpr = MPRelaxSet(struct, use_structure_charge=False)
        self.assertFalse(
            "NELECT" in mpr.incar.keys(),
            "NELECT should not be set when " "use_structure_charge is False",
        )

        struct = Structure(latt, ["Co", "O"], coords)
        mpr = MPRelaxSet(struct)
        self.assertEqual(mpr.incar["MAGMOM"], [0.6, 0.6])
        struct = Structure(latt, ["Co4+", "O"], coords)
        mpr = MPRelaxSet(struct)
        self.assertEqual(mpr.incar["MAGMOM"], [1, 0.6])

        # test passing user_incar_settings and user_kpoint_settings of None
        sets = [
            MPRelaxSet(struct, user_incar_settings=None, user_kpoints_settings=None),
            MPStaticSet(struct, user_incar_settings=None, user_kpoints_settings=None),
            MPNonSCFSet(struct, user_incar_settings=None, user_kpoints_settings=None)
        ]
        for mp_set in sets:
            self.assertNotEqual(mp_set.kpoints, None)
            self.assertNotEqual(mp_set.incar, None)

    def test_get_kpoints(self):
        kpoints = MPRelaxSet(self.structure).kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 5]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

        kpoints = MPRelaxSet(
            self.structure, user_kpoints_settings={"reciprocal_density": 1000}
        ).kpoints
        self.assertEqual(kpoints.kpts, [[6, 10, 13]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

        kpoints_obj = Kpoints(kpts=[[3, 3, 3]])
        kpoints_return = MPRelaxSet(
            self.structure, user_kpoints_settings=kpoints_obj
        ).kpoints
        self.assertEqual(kpoints_return.kpts, [[3, 3, 3]])

        kpoints = self.mitset.kpoints
        self.assertEqual(kpoints.kpts, [[25]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Automatic)

        recip_paramset = MPRelaxSet(self.structure, force_gamma=True)
        recip_paramset.kpoints_settings = {"reciprocal_density": 40}
        kpoints = recip_paramset.kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 5]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_get_vasp_input(self):
        d = self.mitset.get_vasp_input()
        self.assertEqual(d["INCAR"]["ISMEAR"], -5)
        s = self.structure.copy()
        s.make_supercell(4)
        paramset = MPRelaxSet(s)
        d = paramset.get_vasp_input()
        self.assertEqual(d["INCAR"]["ISMEAR"], 0)

    def test_MPMetalRelaxSet(self):
        mpmetalset = MPMetalRelaxSet(self.get_structure("Sn"))
        incar = mpmetalset.incar
        self.assertEqual(incar["ISMEAR"], 1)
        self.assertEqual(incar["SIGMA"], 0.2)
        kpoints = mpmetalset.kpoints
        self.assertArrayAlmostEqual(kpoints.kpts[0], [5, 5, 5])

    def test_as_from_dict(self):
        mitset = MITRelaxSet(self.structure)
        mpset = MPRelaxSet(self.structure)
        mpuserset = MPRelaxSet(
            self.structure,
            user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}},
        )

        d = mitset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v._config_dict["INCAR"]["LDAUU"]["O"]["Fe"], 4)

        d = mpset.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v._config_dict["INCAR"]["LDAUU"]["O"]["Fe"], 5.3)

        d = mpuserset.as_dict()
        v = dec.process_decoded(d)
        # self.assertEqual(type(v), MPVaspInputSet)
        self.assertEqual(
            v.user_incar_settings["MAGMOM"], {"Fe": 10, "S": -5, "Mn3+": 100}
        )

    def test_hubbard_off_and_ediff_override(self):
        p = MPRelaxSet(
            self.structure, user_incar_settings={"LDAU": False, "EDIFF": 1e-10}
        )
        self.assertNotIn("LDAUU", p.incar)
        self.assertEqual(p.incar["EDIFF"], 1e-10)

    def test_write_input(self):
        self.mitset.write_input(".", make_dir_if_not_present=True)
        for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
            self.assertTrue(os.path.exists(f))
        self.assertFalse(os.path.exists("Fe4P4O16.cif"))

        self.mitset.write_input(".", make_dir_if_not_present=True, include_cif=True)
        self.assertTrue(os.path.exists("Fe4P4O16.cif"))
        for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR", "Fe4P4O16.cif"]:
            os.remove(f)

        self.mitset.write_input(".", make_dir_if_not_present=True, potcar_spec=True)

        for f in ["INCAR", "KPOINTS", "POSCAR"]:
            self.assertTrue(os.path.exists(f))
        self.assertFalse(os.path.exists("POTCAR"))
        self.assertTrue(os.path.exists("POTCAR.spec"))
        for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"]:
            os.remove(f)

    def test_user_potcar_settings(self):
        vis = MPRelaxSet(self.structure, user_potcar_settings={"Fe": "Fe"})
        potcar = vis.potcar
        self.assertEqual(potcar.symbols, ["Fe", "P", "O"])


class MPStaticSetTest(PymatgenTest):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        warnings.simplefilter("ignore")

    def test_init(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"

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

        non_prev_vis = MPStaticSet(
            vis.structure, user_incar_settings={"LORBIT": 12, "LWAVE": True}
        )
        self.assertEqual(non_prev_vis.incar["NSW"], 0)
        # Check that the ENCUT and Kpoints style has NOT been inherited.
        self.assertEqual(non_prev_vis.incar["ENCUT"], 520)
        # Check that user incar settings are applied.
        self.assertEqual(non_prev_vis.incar["LORBIT"], 12)
        self.assertTrue(non_prev_vis.incar["LWAVE"])

        self.assertEqual(non_prev_vis.kpoints.style, Kpoints.supported_modes.Gamma)
        v2 = MPStaticSet.from_dict(non_prev_vis.as_dict())
        self.assertEqual(v2.incar["ENCUT"], 520)
        # Check that user incar settings are applied.
        self.assertEqual(v2.incar["LORBIT"], 12)
        leps_vis = MPStaticSet.from_prev_calc(prev_calc_dir=prev_run, lepsilon=True)
        self.assertTrue(leps_vis.incar["LEPSILON"])
        self.assertEqual(leps_vis.incar["IBRION"], 8)
        self.assertNotIn("NPAR", leps_vis.incar)
        self.assertNotIn("NSW", leps_vis.incar)
        self.assertEqual(non_prev_vis.kpoints.kpts, [[11, 10, 10]])
        non_prev_vis = MPStaticSet(vis.structure, reciprocal_density=200)
        self.assertEqual(non_prev_vis.kpoints.kpts, [[14, 12, 12]])
        # Check LCALCPOL flag
        lcalcpol_vis = MPStaticSet.from_prev_calc(prev_calc_dir=prev_run, lcalcpol=True)
        self.assertTrue(lcalcpol_vis.incar["LCALCPOL"])

    def test_user_incar_kspacing(self):
        # Make sure user KSPACING settings properly overrides KPOINTS.
        si = self.get_structure("Si")
        vis = MPRelaxSet(si, user_incar_settings={"KSPACING": 0.22})
        self.assertEqual(vis.incar["KSPACING"], 0.22)
        self.assertEqual(vis.kpoints, None)

    def test_kspacing_override(self):
        # If KSPACING is set and user_kpoints_settings are given,
        # make sure the user_kpoints_settings override KSPACING
        si = self.get_structure("Si")
        vis = MPRelaxSet(
            si,
            user_incar_settings={"KSPACING": 0.22},
            user_kpoints_settings={"reciprocal_density": 1000},
        )
        self.assertEqual(vis.incar.get("KSPACING"), None)
        self.assertIsInstance(vis.kpoints, Kpoints)

    def test_override_from_prev_calc(self):
        # test override_from_prev
        prev_run = self.TEST_FILES_DIR / "relaxation"

        vis = MPStaticSet(_dummy_structure)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Monkhorst)

        # Check LCALCPOL flag
        lcalcpol_vis = MPStaticSet(_dummy_structure, lcalcpol=True)
        lcalcpol_vis = lcalcpol_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(lcalcpol_vis.incar["LCALCPOL"])

    def test_standardize_structure(self):
        sga = SpacegroupAnalyzer(self.get_structure("Si"))
        original_structure = sga.get_conventional_standard_structure()
        sm = StructureMatcher(primitive_cell=False, scale=False)

        vis = MPStaticSet(original_structure)
        self.assertTrue(sm.fit(vis.structure, original_structure))

        vis = MPStaticSet(original_structure, standardize=True)
        self.assertFalse(sm.fit(vis.structure, original_structure))

    def test_write_input_zipped(self):
        vis = MPStaticSet(self.get_structure("Si"))
        vis.write_input(output_dir=".", potcar_spec=True, zip_output=True)

        self.assertTrue(os.path.exists("MPStaticSet.zip"))
        with ZipFile("MPStaticSet.zip", "r") as zip:
            contents = zip.namelist()
            self.assertSetEqual(
                set(contents), {"INCAR", "POSCAR", "POTCAR.spec", "KPOINTS"}
            )
            spec = zip.open("POTCAR.spec", "r").read().decode()
            self.assertEqual(spec, "Si")

        os.remove("MPStaticSet.zip")

    def test_conflicting_arguments(self):
        with pytest.raises(ValueError, match="deprecated"):
            si = self.get_structure("Si")
            vis = MPStaticSet(si, potcar_functional="PBE", user_potcar_functional="PBE")

    def test_grid_size_from_struct(self):
        # TODO grab a bunch_of_calculations store as a list of tuples
        # (structure, ngx, ..., ngxf, ...) where all the grid size values are generated by vasp
        # check that the code produces the same grid sizes
        fname = self.TEST_FILES_DIR / "grid_data_files" / "vasp_inputs_for_grid_check.json"
        parsed_vasp_data = loadfn(fname)
        for tt in parsed_vasp_data:
            ng = [tt['input']['parameters'][ik] for ik in ['NGX', 'NGY', "NGZ"]]
            ngf = [tt['input']['parameters'][ik] for ik in ['NGXF', 'NGYF', "NGZF"]]
            struct = tt['input']['structure']
            static_set = MPStaticSet(struct)
            matched = static_set.calculate_ng() == (ng, ngf)
            self.assertTrue(matched)

    def tearDown(self):
        shutil.rmtree(self.tmp)
        warnings.simplefilter("default")


class MPNonSCFSetTest(PymatgenTest):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        warnings.simplefilter("ignore")

    def test_init(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"
        # check boltztrap mode
        vis = MPNonSCFSet.from_prev_calc(prev_calc_dir=prev_run, mode="Boltztrap")
        self.assertEqual(vis.incar["ISMEAR"], 0)

        # check uniform mode
        vis = MPNonSCFSet.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform")
        self.assertEqual(vis.incar["ISMEAR"], -5)
        self.assertEqual(vis.incar["ISYM"], 2)

        # check uniform mode with automatic nedos
        vis = MPNonSCFSet.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform",
                                         nedos=0)
        self.assertEqual(vis.incar["NEDOS"], 12217)

        # test line mode
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run,
            mode="Line",
            copy_chgcar=False,
            user_incar_settings={"SIGMA": 0.025},
        )

        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        # Check that the user_incar_settings works
        self.assertEqual(vis.incar["SIGMA"], 0.025)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

        # Check as from dict.
        vis = MPNonSCFSet.from_dict(vis.as_dict())
        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

        vis.write_input(self.tmp)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "CHGCAR")))

        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, mode="Line", copy_chgcar=True
        )
        # check ISMEAR set correctly for line mode
        self.assertEqual(vis.incar["ISMEAR"], 0)
        vis.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "CHGCAR")))
        os.remove(os.path.join(self.tmp, "CHGCAR"))

        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, standardize=True, mode="Line", copy_chgcar=True
        )
        vis.write_input(self.tmp)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "CHGCAR")))

    def test_override_from_prev(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"

        # test override_from_prev
        vis = MPNonSCFSet(_dummy_structure, mode="Boltztrap")
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["ISMEAR"], 0)

        vis = MPNonSCFSet(_dummy_structure, mode="Uniform")
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["ISMEAR"], -5)
        self.assertEqual(vis.incar["ISYM"], 2)

        vis = MPNonSCFSet(_dummy_structure, mode="Uniform", nedos=0)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["NEDOS"], 12217)

        # test line mode
        vis = MPNonSCFSet(
            _dummy_structure,
            mode="Line",
            copy_chgcar=False,
            user_incar_settings={"SIGMA": 0.025},
        )
        vis.override_from_prev_calc(prev_calc_dir=prev_run)

        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["ENCUT"], 600)
        self.assertEqual(vis.incar["SIGMA"], 0.025)
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Reciprocal)

        vis = MPNonSCFSet(_dummy_structure, mode="Line", copy_chgcar=True)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["ISMEAR"], 0)
        vis.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "CHGCAR")))
        os.remove(os.path.join(self.tmp, "CHGCAR"))

        vis = MPNonSCFSet(
            _dummy_structure, standardize=True, mode="Line", copy_chgcar=True
        )
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        vis.write_input(self.tmp)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "CHGCAR")))

    def test_kpoints(self):
        # test k-points are generated in the correct format
        prev_run = self.TEST_FILES_DIR / "relaxation"
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, mode="Uniform", copy_chgcar=False
        )
        self.assertEqual(np.array(vis.kpoints.kpts).shape, (1, 3))

        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, mode="Line", copy_chgcar=False
        )
        self.assertNotEqual(np.array(vis.kpoints.kpts).shape, (1, 3))

        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run, mode="Boltztrap", copy_chgcar=False
        )
        self.assertNotEqual(np.array(vis.kpoints.kpts).shape, (1, 3))

    def test_optics(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run,
            copy_chgcar=False,
            optics=True,
            mode="Uniform",
            nedos=2001,
        )

        self.assertEqual(vis.incar["NSW"], 0)
        # Check that the ENCUT has been inherited.
        self.assertEqual(vis.incar["ENCUT"], 600)

        # check NEDOS and ISMEAR set correctly
        self.assertEqual(vis.incar["NEDOS"], 2001)
        self.assertEqual(vis.incar["ISMEAR"], -5)
        self.assertEqual(vis.incar["ISYM"], 2)

        self.assertTrue(vis.incar["LOPTICS"])
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Gamma)

    def test_user_kpoint_override(self):
        user_kpoints_override = Kpoints(
            style=Kpoints.supported_modes.Gamma, kpts=((1, 1, 1),)
        )  # the default kpoints style is reciprocal

        prev_run = self.TEST_FILES_DIR / "relaxation"
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=prev_run,
            copy_chgcar=False,
            optics=True,
            mode="Uniform",
            nedos=2001,
            user_kpoints_settings=user_kpoints_override,
        )
        self.assertEqual(vis.kpoints.style, Kpoints.supported_modes.Gamma)

    def tearDown(self):
        shutil.rmtree(self.tmp)
        warnings.simplefilter("default")


class MagmomLdauTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_structure_from_prev_run(self):
        vrun = Vasprun(self.TEST_FILES_DIR / "vasprun.xml.magmom_ldau")
        structure = vrun.final_structure
        poscar = Poscar(structure)
        structure_decorated = get_structure_from_prev_run(vrun)
        ldau_ans = {"LDAUU": [5.3, 0.0], "LDAUL": [2, 0], "LDAUJ": [0.0, 0.0]}
        magmom_ans = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
        ldau_dict = {}
        for key in ("LDAUU", "LDAUJ", "LDAUL"):
            if hasattr(structure_decorated[0], key.lower()):
                m = dict(
                    [
                        (site.specie.symbol, getattr(site, key.lower()))
                        for site in structure_decorated
                    ]
                )
                ldau_dict[key] = [m[sym] for sym in poscar.site_symbols]
        magmom = [site.magmom for site in structure_decorated]
        self.assertEqual(ldau_dict, ldau_ans)
        self.assertEqual(magmom, magmom_ans)

    def test_ln_magmom(self):
        YAML_PATH = os.path.join(os.path.dirname(__file__), "../VASPIncarBase.yaml")
        MAGMOM_SETTING = loadfn(YAML_PATH)["INCAR"]["MAGMOM"]
        structure = Structure.from_file(self.TEST_FILES_DIR / "La4Fe4O12.cif")
        structure.add_oxidation_state_by_element({"La": +3, "Fe": +3, "O": -2})
        for ion in MAGMOM_SETTING:
            s = structure.copy()
            s.replace_species({"La3+": ion})
            vis = MPRelaxSet(s)
            fe_pos = vis.poscar.comment.index("Fe")
            if fe_pos == 0:
                magmom_ans = [5] * 4 + [MAGMOM_SETTING[ion]] * 4 + [0.6] * 12
            else:
                magmom_ans = [MAGMOM_SETTING[ion]] * 4 + [5] * 4 + [0.6] * 12

            self.assertEqual(vis.incar["MAGMOM"], magmom_ans)


class MITMDSetTest(PymatgenTest):
    def setUp(self):
        filepath = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.structure
        self.mitmdparam = MITMDSet(self.struct, 300, 1200, 10000)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_params(self):
        param = self.mitmdparam
        syms = param.potcar_symbols
        self.assertEqual(syms, ["Fe", "P", "O"])
        incar = param.incar
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar["EDIFF"], 1e-5)
        kpoints = param.kpoints
        self.assertEqual(kpoints.kpts, [(1, 1, 1)])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_as_from_dict(self):
        d = self.mitmdparam.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MITMDSet)
        self.assertEqual(v._config_dict["INCAR"]["TEBEG"], 300)
        self.assertEqual(v._config_dict["INCAR"]["PREC"], "Low")


class MVLNPTMDSetTest(PymatgenTest):
    def setUp(self):
        file_path = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(file_path)
        self.struct = poscar.structure
        self.mvl_npt_set = MVLNPTMDSet(
            self.struct, start_temp=0, end_temp=300, nsteps=1000
        )
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        npt_set = self.mvl_npt_set

        syms = npt_set.potcar_symbols
        self.assertEqual(syms, ["Fe", "P", "O"])

        incar = npt_set.incar
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar["EDIFF"], 1e-5)
        self.assertEqual(incar["LANGEVIN_GAMMA_L"], 1)
        self.assertEqual(incar["LANGEVIN_GAMMA"], [10, 10, 10])
        enmax = max(
            [npt_set.potcar[i].keywords["ENMAX"] for i in range(self.struct.ntypesp)]
        )
        self.assertAlmostEqual(incar["ENCUT"], 1.5 * enmax)
        self.assertEqual(incar["IALGO"], 48)
        self.assertEqual(incar["ISIF"], 3)
        self.assertEqual(incar["MDALGO"], 3)
        self.assertEqual(incar["SMASS"], 0)
        self.assertEqual(incar["PREC"], "Low")

        kpoints = npt_set.kpoints
        self.assertEqual(kpoints.kpts, [(1, 1, 1)])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_as_from_dict(self):
        d = self.mvl_npt_set.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MVLNPTMDSet)
        self.assertEqual(v._config_dict["INCAR"]["NSW"], 1000)


class MITNEBSetTest(PymatgenTest):
    def setUp(self):
        c1 = [[0.5] * 3, [0.9] * 3]
        c2 = [[0.5] * 3, [0.9, 0.1, 0.1]]
        s1 = Structure(Lattice.cubic(5), ["Si", "Si"], c1)
        s2 = Structure(Lattice.cubic(5), ["Si", "Si"], c2)
        structs = []
        for s in s1.interpolate(s2, 3, pbc=True):
            structs.append(Structure.from_sites(s.sites, to_unit_cell=True))
        self.structures = structs
        self.vis = MITNEBSet(self.structures)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_potcar_symbols(self):
        syms = self.vis.potcar_symbols
        self.assertEqual(syms, ["Si"])

    def test_incar(self):
        incar = self.vis.incar
        self.assertNotIn("LDAUU", incar)
        self.assertAlmostEqual(incar["EDIFF"], 0.00001)

    def test_kpoints(self):
        kpoints = self.vis.kpoints
        self.assertEqual(kpoints.kpts, [[25]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Automatic)

    def test_as_from_dict(self):
        d = self.vis.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(v._config_dict["INCAR"]["IMAGES"], 2)

    def test_write_input(self):
        self.vis.write_input(
            ".", write_cif=True, write_endpoint_inputs=True, write_path_cif=True
        )
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
        for d in ["00", "01", "02", "03"]:
            shutil.rmtree(d)
        for f in ["INCAR", "KPOINTS", "POTCAR", "path.cif"]:
            os.remove(f)


class MPSOCSetTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_from_prev_calc(self):
        prev_run = self.TEST_FILES_DIR / "fe_monomer"
        vis = MPSOCSet.from_prev_calc(
            prev_calc_dir=prev_run,
            magmom=[3],
            saxis=(1, 0, 0),
            user_incar_settings={"SIGMA": 0.025},
        )
        self.assertEqual(vis.incar["ISYM"], -1)
        self.assertTrue(vis.incar["LSORBIT"])
        self.assertEqual(vis.incar["ICHARG"], 11)
        self.assertEqual(vis.incar["SAXIS"], [1, 0, 0])
        self.assertEqual(vis.incar["MAGMOM"], [[0, 0, 3]])
        self.assertEqual(vis.incar["SIGMA"], 0.025)

    def test_override_from_prev_calc(self):
        # test override_from_prev_calc
        prev_run = self.TEST_FILES_DIR / "fe_monomer"
        vis = MPSOCSet(
            _dummy_structure,
            magmom=[3],
            saxis=(1, 0, 0),
            user_incar_settings={"SIGMA": 0.025},
        )
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(vis.incar["ISYM"], -1)
        self.assertTrue(vis.incar["LSORBIT"])
        self.assertEqual(vis.incar["ICHARG"], 11)
        self.assertEqual(vis.incar["SAXIS"], [1, 0, 0])
        self.assertEqual(vis.incar["MAGMOM"], [[0, 0, 3]])
        self.assertEqual(vis.incar["SIGMA"], 0.025)


class MPNMRSetTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        filepath = self.TEST_FILES_DIR / "Li.cif"
        structure = Structure.from_file(filepath)

        vis = MPNMRSet(structure)
        self.assertTrue(vis.incar.get("LCHIMAG", None))
        self.assertEqual(vis.incar.get("QUAD_EFG", None), None)

        vis = MPNMRSet(structure, mode="efg")
        self.assertFalse(vis.incar.get("LCHIMAG", None))
        self.assertEqual(vis.incar.get("QUAD_EFG", None), [-0.808])

        vis = MPNMRSet(structure, mode="efg", isotopes=["Li-7"])
        self.assertFalse(vis.incar.get("LCHIMAG", None))
        self.assertEqual(vis.incar.get("QUAD_EFG", None), [-40.1])


class MVLSlabSetTest(PymatgenTest):
    def setUp(self):
        s = self.get_structure("Li2O")
        gen = SlabGenerator(s, (1, 0, 0), 10, 10)
        self.slab = gen.get_slab()
        self.bulk = self.slab.oriented_unit_cell

        vis_bulk = MVLSlabSet(self.bulk, bulk=True)
        vis = MVLSlabSet(self.slab)
        vis_dipole = MVLSlabSet(self.slab, auto_dipole=True)

        self.d_bulk = vis_bulk.get_vasp_input()
        self.d_slab = vis.get_vasp_input()
        self.d_dipole = vis_dipole.get_vasp_input()
        self.vis = vis
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_user_incar_settings(self):
        # Make sure user incar settings properly override AMIX.
        si = self.get_structure("Si")
        vis = MVLSlabSet(si, user_incar_settings={"AMIX": 0.1})
        self.assertEqual(vis.incar["AMIX"], 0.1)

    def test_bulk(self):
        incar_bulk = self.d_bulk["INCAR"]
        poscar_bulk = self.d_bulk["POSCAR"]

        self.assertEqual(incar_bulk["ISIF"], 3)
        self.assertEqual(incar_bulk["EDIFF"], 1e-4)
        self.assertEqual(incar_bulk["EDIFFG"], -0.02)
        self.assertEqual(poscar_bulk.structure.formula, self.bulk.formula)

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
        self.assertEqual(potcar_slab.functional, "PBE")
        self.assertEqual(potcar_slab.symbols[1], u"O")
        self.assertEqual(potcar_slab.symbols[0], u"Li_sv")
        self.assertEqual(poscar_slab.structure.formula, self.slab.formula)
        # Test auto-dipole
        dipole_incar = self.d_dipole["INCAR"]
        self.assertTrue(dipole_incar["LDIPOL"])
        self.assertArrayAlmostEqual(
            dipole_incar["DIPOL"], [0.2323, 0.2323, 0.2165], decimal=4
        )
        self.assertEqual(dipole_incar["IDIPOL"], 3)

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

    def test_as_dict(self):
        vis_dict = self.vis.as_dict()
        new = MVLSlabSet.from_dict(vis_dict)


class MVLElasticSetTest(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        mvlparam = MVLElasticSet(self.get_structure("Graphite"))
        incar = mvlparam.incar
        self.assertEqual(incar["IBRION"], 6)
        self.assertEqual(incar["NFREE"], 2)
        self.assertEqual(incar["POTIM"], 0.015)
        self.assertNotIn("NPAR", incar)


class MVLGWSetTest(PymatgenTest):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.s = PymatgenTest.get_structure("Li2O")
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")
        shutil.rmtree(self.tmp)

    def test_static(self):
        mvlgwsc = MVLGWSet(self.s)
        incar = mvlgwsc.incar
        self.assertEqual(incar["SIGMA"], 0.01)
        kpoints = mvlgwsc.kpoints
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        symbols = mvlgwsc.potcar.symbols
        self.assertEqual(symbols, ["Li_sv_GW", "O_GW"])

    def test_diag(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"
        mvlgwdiag = MVLGWSet.from_prev_calc(prev_run, copy_wavecar=True, mode="diag")
        mvlgwdiag.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVECAR")))
        self.assertEqual(mvlgwdiag.incar["NBANDS"], 32)
        self.assertEqual(mvlgwdiag.incar["ALGO"], "Exact")
        self.assertTrue(mvlgwdiag.incar["LOPTICS"])

        # test override_from_prev_calc
        mvlgwdiag = MVLGWSet(_dummy_structure, copy_wavecar=True, mode="diag")
        mvlgwdiag.override_from_prev_calc(prev_calc_dir=prev_run)
        mvlgwdiag.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVECAR")))
        self.assertEqual(mvlgwdiag.incar["NBANDS"], 32)
        self.assertEqual(mvlgwdiag.incar["ALGO"], "Exact")
        self.assertTrue(mvlgwdiag.incar["LOPTICS"])

    def test_bse(self):
        prev_run = self.TEST_FILES_DIR / "relaxation"
        mvlgwgbse = MVLGWSet.from_prev_calc(prev_run, copy_wavecar=True, mode="BSE")
        mvlgwgbse.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVECAR")))
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVEDER")))

        prev_run = self.TEST_FILES_DIR / "relaxation"
        mvlgwgbse = MVLGWSet.from_prev_calc(prev_run, copy_wavecar=False, mode="GW")
        self.assertEqual(mvlgwgbse.incar["NOMEGA"], 80)
        self.assertEqual(mvlgwgbse.incar["ENCUTGW"], 250)
        self.assertEqual(mvlgwgbse.incar["ALGO"], "GW0")
        mvlgwgbse1 = MVLGWSet.from_prev_calc(prev_run, copy_wavecar=False, mode="BSE")
        self.assertEqual(mvlgwgbse1.incar["ANTIRES"], 0)
        self.assertEqual(mvlgwgbse1.incar["NBANDSO"], 20)
        self.assertEqual(mvlgwgbse1.incar["ALGO"], "BSE")

        # test override_from_prev_calc
        prev_run = self.TEST_FILES_DIR / "relaxation"
        mvlgwgbse = MVLGWSet(_dummy_structure, copy_wavecar=True, mode="BSE")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        mvlgwgbse.write_input(self.tmp)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVECAR")))
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "WAVEDER")))

        prev_run = self.TEST_FILES_DIR / "relaxation"
        mvlgwgbse = MVLGWSet(_dummy_structure, copy_wavecar=True, mode="GW")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(mvlgwgbse.incar["NOMEGA"], 80)
        self.assertEqual(mvlgwgbse.incar["ENCUTGW"], 250)
        self.assertEqual(mvlgwgbse.incar["ALGO"], "GW0")

        mvlgwgbse1 = MVLGWSet(_dummy_structure, copy_wavecar=False, mode="BSE")
        mvlgwgbse1.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertEqual(mvlgwgbse1.incar["ANTIRES"], 0)
        self.assertEqual(mvlgwgbse1.incar["NBANDSO"], 20)
        self.assertEqual(mvlgwgbse1.incar["ALGO"], "BSE")


class MPHSEBSTest(PymatgenTest):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_init(self):
        prev_run = self.TEST_FILES_DIR / "static_silicon"
        vis = MPHSEBSSet.from_prev_calc(prev_calc_dir=prev_run, mode="uniform")
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(len(vis.kpoints.kpts), 16)

        vis = MPHSEBSSet.from_prev_calc(prev_calc_dir=prev_run, mode="gap")
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(len(vis.kpoints.kpts), 18)

        vis = MPHSEBSSet.from_prev_calc(prev_calc_dir=prev_run, mode="line")
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(vis.incar["HFSCREEN"], 0.2)
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["ISYM"], 3)
        self.assertEqual(len(vis.kpoints.kpts), 180)

    def test_override_from_prev_calc(self):
        prev_run = self.TEST_FILES_DIR / "static_silicon"
        vis = MPHSEBSSet(_dummy_structure, mode="uniform")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(len(vis.kpoints.kpts), 16)

        vis = MPHSEBSSet(_dummy_structure, mode="gap")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(len(vis.kpoints.kpts), 18)

        vis = MPHSEBSSet(_dummy_structure, mode="line")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(vis.incar["LHFCALC"])
        self.assertEqual(vis.incar["HFSCREEN"], 0.2)
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["ISYM"], 3)
        self.assertEqual(len(vis.kpoints.kpts), 180)


class MVLScanRelaxSetTest(PymatgenTest):
    def setUp(self):
        file_path = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(file_path)
        self.struct = poscar.structure
        self.mvl_scan_set = MVLScanRelaxSet(
            self.struct, potcar_functional="PBE_52", user_incar_settings={"NSW": 500}
        )
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        incar = self.mvl_scan_set.incar
        self.assertIn("METAGGA", incar)
        self.assertIn("LASPH", incar)
        self.assertIn("ADDGRID", incar)
        self.assertEqual(incar["NSW"], 500)

        # Test SCAN+rVV10
        scan_rvv10_set = MVLScanRelaxSet(self.struct, vdw="rVV10")
        self.assertEqual(scan_rvv10_set.incar["BPARAM"], 15.7)

    def test_potcar(self):
        self.assertEqual(self.mvl_scan_set.potcar.functional, "PBE_52")

        test_potcar_set_1 = MVLScanRelaxSet(self.struct, potcar_functional="PBE_54")
        self.assertEqual(test_potcar_set_1.potcar.functional, "PBE_54")

        self.assertRaises(
            ValueError, MVLScanRelaxSet, self.struct, potcar_functional="PBE"
        )

    def test_as_from_dict(self):
        d = self.mvl_scan_set.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MVLScanRelaxSet)
        self.assertEqual(v._config_dict["INCAR"]["METAGGA"], "SCAN")
        self.assertEqual(v.user_incar_settings["NSW"], 500)


class MPScanRelaxSetTest(PymatgenTest):
    def setUp(self):
        file_path = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(file_path)
        self.struct = poscar.structure
        self.mp_scan_set = MPScanRelaxSet(
            self.struct, potcar_functional="PBE_52", user_incar_settings={"NSW": 500}
        )
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        incar = self.mp_scan_set.incar
        self.assertEqual(incar["METAGGA"], "R2scan")
        self.assertEqual(incar["LASPH"], True)
        self.assertEqual(incar["ENAUG"], 1360)
        self.assertEqual(incar["ENCUT"], 680)
        self.assertEqual(incar["NSW"], 500)
        # the default POTCAR contains metals
        self.assertEqual(incar["KSPACING"], 0.22)
        self.assertEqual(incar["ISMEAR"], 2)
        self.assertEqual(incar["SIGMA"], 0.2)

    def test_scan_substitute(self):
        mp_scan_sub = MPScanRelaxSet(
            self.struct, potcar_functional="PBE_52", user_incar_settings={"METAGGA": "SCAN"})
        incar = mp_scan_sub.incar
        self.assertEqual(incar["METAGGA"], "Scan")

    def test_nonmetal(self):
        # Test that KSPACING and ISMEAR change with a nonmetal structure
        file_path = self.TEST_FILES_DIR / "POSCAR.O2"
        struct = Poscar.from_file(file_path, check_for_POTCAR=False).structure
        scan_nonmetal_set = MPScanRelaxSet(struct, bandgap=1.1)
        incar = scan_nonmetal_set.incar
        self.assertAlmostEqual(incar["KSPACING"], 0.29125, places=5)
        self.assertEqual(incar["ISMEAR"], -5)
        self.assertEqual(incar["SIGMA"], 0.05)

    def test_incar_overrides(self):
        # use 'user_incar_settings' to override the KSPACING, ISMEAR, and SIGMA
        # parameters that MPScanSet normally determines
        mp_scan_set2 = MPScanRelaxSet(
            self.struct,
            user_incar_settings={"KSPACING": 0.5, "ISMEAR": 0, "SIGMA": 0.05},
        )
        incar = mp_scan_set2.incar
        self.assertEqual(incar["KSPACING"], 0.5)
        self.assertEqual(incar["ISMEAR"], 0)
        self.assertEqual(incar["SIGMA"], 0.05)

    # Test SCAN+rVV10
    def test_rvv10(self):
        scan_rvv10_set = MPScanRelaxSet(self.struct, vdw="rVV10")
        self.assertIn("LUSE_VDW", scan_rvv10_set.incar)
        self.assertEqual(scan_rvv10_set.incar["BPARAM"], 15.7)

    def test_other_vdw(self):
        # should raise a warning.
        # IVDW key should not be present in the incar
        with pytest.warns(UserWarning, match=r"not supported at this time"):
            scan_vdw_set = MPScanRelaxSet(self.struct, vdw="DFTD3")
            self.assertNotIn("LUSE_VDW", scan_vdw_set.incar)
            self.assertNotIn("IVDW", scan_vdw_set.incar)

    def test_potcar(self):
        self.assertEqual(self.mp_scan_set.potcar.functional, "PBE_52")

        # the default functional should be PBE_54
        test_potcar_set_1 = MPScanRelaxSet(self.struct)
        self.assertEqual(test_potcar_set_1.potcar.functional, "PBE_54")

        self.assertRaises(
            ValueError, MPScanRelaxSet, self.struct, potcar_functional="PBE"
        )

    def test_as_from_dict(self):
        d = self.mp_scan_set.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MPScanRelaxSet)
        self.assertEqual(v._config_dict["INCAR"]["METAGGA"], "R2SCAN")
        self.assertEqual(v.user_incar_settings["NSW"], 500)

    def test_write_input(self):
        self.mp_scan_set.write_input(
            "."
        )
        self.assertTrue(os.path.exists("INCAR"))
        self.assertFalse(os.path.exists("KPOINTS"))
        self.assertTrue(os.path.exists("POTCAR"))
        self.assertTrue(os.path.exists("POSCAR"))

        for f in ["INCAR", "POSCAR", "POTCAR"]:
            os.remove(f)


class MPScanStaticSetTest(PymatgenTest):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        warnings.simplefilter("ignore")

    def test_init(self):
        # test inheriting from a previous SCAN relaxation
        prev_run = self.TEST_FILES_DIR / "scan_relaxation"

        vis = MPScanStaticSet.from_prev_calc(prev_calc_dir=prev_run)
        # check that StaticSet settings were applied
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["LREAL"], False)
        self.assertEqual(vis.incar["LORBIT"], 11)
        self.assertEqual(vis.incar["LVHAR"], True)
        self.assertEqual(vis.incar["ISMEAR"], -5)
        # Check that ENCUT and other INCAR settings were inherited.
        self.assertEqual(vis.incar["ENCUT"], 680)
        self.assertEqual(vis.incar["METAGGA"], "R2scan")
        self.assertEqual(vis.incar["KSPACING"], 0.34292842)

        # Check as from dict.
        # check that StaticSet settings were applied
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["LREAL"], False)
        self.assertEqual(vis.incar["LORBIT"], 11)
        self.assertEqual(vis.incar["LVHAR"], True)
        # Check that ENCUT and KSPACING were inherited.
        self.assertEqual(vis.incar["ENCUT"], 680)
        self.assertEqual(vis.incar["METAGGA"], "R2scan")
        self.assertEqual(vis.incar["KSPACING"], 0.34292842)

        non_prev_vis = MPScanStaticSet(
            vis.structure, user_incar_settings={"ENCUT": 800, "LORBIT": 12, "LWAVE": True}
        )
        # check that StaticSet settings were applied
        self.assertEqual(non_prev_vis.incar["NSW"], 0)
        self.assertEqual(non_prev_vis.incar["LREAL"], False)
        self.assertEqual(non_prev_vis.incar["LVHAR"], True)
        self.assertEqual(vis.incar["ISMEAR"], -5)
        # Check that ENCUT and other INCAR settings were inherited.
        self.assertEqual(non_prev_vis.incar["METAGGA"], "R2scan")
        # the KSPACING will have the default value here, since no previous calc
        self.assertEqual(non_prev_vis.incar["KSPACING"], 0.22)
        # Check that user incar settings are applied.
        self.assertEqual(non_prev_vis.incar["ENCUT"], 800)
        self.assertEqual(non_prev_vis.incar["LORBIT"], 12)
        self.assertTrue(non_prev_vis.incar["LWAVE"])

        v2 = MPScanStaticSet.from_dict(non_prev_vis.as_dict())
        # Check that user incar settings are applied.
        self.assertEqual(v2.incar["ENCUT"], 800)
        self.assertEqual(v2.incar["LORBIT"], 12)
        self.assertTrue(non_prev_vis.incar["LWAVE"])

        # Check LCALCPOL flag
        lcalcpol_vis = MPScanStaticSet.from_prev_calc(prev_calc_dir=prev_run, lcalcpol=True)
        self.assertTrue(lcalcpol_vis.incar["LCALCPOL"])

        # Check LEPSILON flag
        lepsilon_vis = MPScanStaticSet.from_prev_calc(prev_calc_dir=prev_run, lepsilon=True)
        self.assertTrue(lepsilon_vis.incar["LEPSILON"])
        self.assertTrue(lepsilon_vis.incar["LPEAD"])
        self.assertEqual(lepsilon_vis.incar["IBRION"], 8)
        self.assertIsNone(lepsilon_vis.incar.get("NSW"))
        self.assertIsNone(lepsilon_vis.incar.get("NPAR"))

    def test_override_from_prev_calc(self):
        # test override_from_prev
        prev_run = self.TEST_FILES_DIR / "scan_relaxation"

        vis = MPScanStaticSet(_dummy_structure)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        # check that StaticSet settings were applied
        self.assertEqual(vis.incar["NSW"], 0)
        self.assertEqual(vis.incar["LREAL"], False)
        self.assertEqual(vis.incar["LORBIT"], 11)
        self.assertEqual(vis.incar["LVHAR"], True)
        self.assertEqual(vis.incar["ISMEAR"], -5)
        # Check that ENCUT and other INCAR settings were inherited.
        self.assertEqual(vis.incar["ENCUT"], 680)
        self.assertEqual(vis.incar["METAGGA"], "R2scan")
        self.assertEqual(vis.incar["KSPACING"], 0.34292842)

        # Check LCALCPOL flag
        lcalcpol_vis = MPScanStaticSet(_dummy_structure, lcalcpol=True)
        lcalcpol_vis = lcalcpol_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(lcalcpol_vis.incar["LCALCPOL"])

        # Check LEPSILON flag
        lepsilon_vis = MPScanStaticSet(_dummy_structure, lepsilon=True)
        lepsilon_vis = lepsilon_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        self.assertTrue(lepsilon_vis.incar["LEPSILON"])
        self.assertTrue(lepsilon_vis.incar["LPEAD"])
        self.assertEqual(lepsilon_vis.incar["IBRION"], 8)
        self.assertIsNone(lepsilon_vis.incar.get("NSW"))
        self.assertIsNone(lepsilon_vis.incar.get("NPAR"))

    def test_conflicting_arguments(self):
        with pytest.raises(ValueError, match="deprecated"):
            si = self.get_structure("Si")
            vis = MPScanStaticSet(si, potcar_functional="PBE", user_potcar_functional="PBE")

    def tearDown(self):
        shutil.rmtree(self.tmp)
        warnings.simplefilter("default")


class FuncTest(PymatgenTest):
    def test_batch_write_input(self):
        structures = [
            PymatgenTest.get_structure("Li2O"),
            PymatgenTest.get_structure("LiFePO4"),
        ]
        batch_write_input(structures)
        for d in ["Li4Fe4P4O16_1", "Li2O1_0"]:
            for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
                self.assertTrue(os.path.exists(os.path.join(d, f)))
        for d in ["Li4Fe4P4O16_1", "Li2O1_0"]:
            shutil.rmtree(d)


class MVLGBSetTest(PymatgenTest):
    def setUp(self):
        filepath = self.TEST_FILES_DIR / "Li.cif"
        self.s = Structure.from_file(filepath)

        self.bulk = MVLGBSet(self.s)
        self.slab = MVLGBSet(self.s, slab_mode=True)

        self.d_bulk = self.bulk.get_vasp_input()
        self.d_slab = self.slab.get_vasp_input()
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_bulk(self):
        incar_bulk = self.d_bulk["INCAR"]
        self.assertEqual(incar_bulk["ISIF"], 3)

    def test_slab(self):
        incar_slab = self.d_slab["INCAR"]
        self.assertEqual(incar_slab["ISIF"], 2)

    def test_kpoints(self):
        kpoints = self.d_slab["KPOINTS"]
        k_a = int(40 / (self.s.lattice.abc[0]) + 0.5)
        k_b = int(40 / (self.s.lattice.abc[1]) + 0.5)
        self.assertEqual(kpoints.kpts, [[k_a, k_b, 1]])


class MVLRelax52SetTest(PymatgenTest):
    def setUp(self):
        file_path = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(file_path)
        self.struct = poscar.structure
        self.mvl_rlx_set = MVLRelax52Set(
            self.struct, potcar_functional="PBE_54", user_incar_settings={"NSW": 500}
        )
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_incar(self):
        incar = self.mvl_rlx_set.incar
        self.assertIn("NSW", incar)
        self.assertEqual(incar["LREAL"], "Auto")

    def test_potcar(self):
        self.assertEqual(self.mvl_rlx_set.potcar.functional, "PBE_54")
        self.assertIn("Fe", self.mvl_rlx_set.potcar.symbols)

        self.struct.remove_species(["Fe"])
        test_potcar_set_1 = MVLRelax52Set(self.struct, potcar_functional="PBE_52")
        self.assertEqual(test_potcar_set_1.potcar.functional, "PBE_52")

        self.assertRaises(
            ValueError, MVLRelax52Set, self.struct, potcar_functional="PBE"
        )

    def test_potcar_functional_warning(self):
        with pytest.warns(FutureWarning, match="argument is deprecated"):
            test_potcar_set_1 = MVLRelax52Set(self.struct, potcar_functional="PBE_52")

    def test_as_from_dict(self):
        d = self.mvl_rlx_set.as_dict()
        v = dec.process_decoded(d)
        self.assertEqual(type(v), MVLRelax52Set)
        self.assertEqual(v.incar["NSW"], 500)


class LobsterSetTest(PymatgenTest):

    def setUp(self):
        file_path = self.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(file_path)
        self.struct = poscar.structure
        # test for different parameters!
        self.lobsterset1 = LobsterSet(self.struct, isym=-1, ismear=-5)
        self.lobsterset2 = LobsterSet(self.struct, isym=0, ismear=0)
        # only allow isym=-1 and isym=0
        with self.assertRaises(ValueError):
            self.lobsterset_new = LobsterSet(self.struct, isym=2, ismear=0)
        # only allow ismear=-5 and ismear=0
        with self.assertRaises(ValueError):
            self.lobsterset_new = LobsterSet(self.struct, isym=-1, ismear=2)
        # test if one can still hand over grid density of kpoints
        self.lobsterset3 = LobsterSet(
            self.struct, isym=0, ismear=0, user_kpoints_settings={"grid_density": 6000}
        )
        # check if users can overwrite settings in this class with the help of user_incar_settings
        self.lobsterset4 = LobsterSet(self.struct, user_incar_settings={"ALGO": "Fast"})
        # use basis functions supplied by user
        self.lobsterset5 = LobsterSet(
            self.struct,
            user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s", "O": "2p 2s"},
        )
        with self.assertRaises(ValueError):
            self.lobsterset6 = LobsterSet(
                self.struct, user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s"}
            )
        self.lobsterset7 = LobsterSet(
            self.struct,
            address_basis_file=os.path.join(MODULE_DIR, "../../lobster/lobster_basis/BASIS_PBE_54_standard.yaml"),
        )
        with pytest.warns(BadInputSetWarning, match="Overriding the POTCAR"):
            self.lobsterset6 = LobsterSet(self.struct)

    def test_incar(self):
        incar1 = self.lobsterset1.incar
        self.assertIn("NBANDS", incar1)
        self.assertEqual(incar1["NBANDS"], 116)
        self.assertEqual(incar1["NSW"], 0)
        self.assertEqual(incar1["ISMEAR"], -5)
        self.assertEqual(incar1["ISYM"], -1)
        self.assertEqual(incar1["ALGO"], "Normal")
        self.assertEqual(incar1["EDIFF"], 1e-6)
        incar2 = self.lobsterset2.incar
        self.assertEqual(incar2["ISYM"], 0)
        self.assertEqual(incar2["ISMEAR"], 0)
        incar4 = self.lobsterset4.incar
        self.assertEqual(incar4["ALGO"], "Fast")

    def test_kpoints(self):
        kpoints1 = self.lobsterset1.kpoints
        self.assertTrue(kpoints1.comment.split(" ")[6], 6138)
        kpoints2 = self.lobsterset2.kpoints
        self.assertTrue(kpoints2.comment.split(" ")[6], 6138)
        kpoints3 = self.lobsterset3.kpoints
        self.assertTrue(kpoints3.comment.split(" ")[6], 6000)

    def test_potcar(self):
        # PBE_54 is preferred at the moment
        self.assertEqual(self.lobsterset1.potcar_functional, "PBE_54")

    def test_as_from_dict(self):
        dict_here = self.lobsterset1.as_dict()

        lobsterset_new = LobsterSet.from_dict(dict_here)
        # test relevant parts again
        incar1 = lobsterset_new.incar
        self.assertIn("NBANDS", incar1)
        self.assertEqual(incar1["NBANDS"], 116)
        self.assertEqual(incar1["NSW"], 0)
        self.assertEqual(incar1["NSW"], 0)
        self.assertEqual(incar1["ISMEAR"], -5)
        self.assertEqual(incar1["ISYM"], -1)
        self.assertEqual(incar1["ALGO"], "Normal")
        kpoints1 = lobsterset_new.kpoints
        self.assertTrue(kpoints1.comment.split(" ")[6], 6138)
        self.assertEqual(lobsterset_new.potcar_functional, "PBE_54")


_dummy_structure = Structure(
    [1, 0, 0, 0, 1, 0, 0, 0, 1],
    ["I"],
    [[0, 0, 0]],
    site_properties={"magmom": [[0, 0, 1]]},
)

if __name__ == "__main__":
    unittest.main()
