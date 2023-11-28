from __future__ import annotations

import hashlib
import os
from glob import glob
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pytest
from _pytest.monkeypatch import MonkeyPatch
from monty.json import MontyDecoder
from monty.serialization import loadfn
from numpy.testing import assert_allclose
from pytest import approx, mark

import pymatgen
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import SETTINGS, Lattice, Species, Structure
from pymatgen.core.composition import Composition
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.units import FloatWithUnit
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, PotcarSingle
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import (
    BadInputSetWarning,
    DictSet,
    LobsterSet,
    MatPESStaticSet,
    MITMDSet,
    MITNEBSet,
    MITRelaxSet,
    MPAbsorptionSet,
    MPHSEBSSet,
    MPHSERelaxSet,
    MPMDSet,
    MPMetalRelaxSet,
    MPNMRSet,
    MPNonSCFSet,
    MPRelaxSet,
    MPScanRelaxSet,
    MPScanStaticSet,
    MPSOCSet,
    MPStaticSet,
    MVLElasticSet,
    MVLGBSet,
    MVLGWSet,
    MVLNPTMDSet,
    MVLRelax52Set,
    MVLScanRelaxSet,
    MVLSlabSet,
    VaspInputSet,
    batch_write_input,
    get_structure_from_prev_run,
    get_valid_magmom_struct,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

MODULE_DIR = Path(pymatgen.io.vasp.__file__).parent

dec = MontyDecoder()

NO_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR") is None
skip_if_no_psp_dir = mark.skipif(NO_PSP_DIR, reason="PMG_VASP_PSP_DIR is not set.")

dummy_structure = Structure(
    [1, 0, 0, 0, 1, 0, 0, 0, 1],
    ["I"],
    [[0, 0, 0]],
    site_properties={"magmom": [[0, 0, 1]]},
)


@pytest.mark.parametrize(
    "input_set",
    [MPRelaxSet, MPHSERelaxSet, MVLRelax52Set, MPAbsorptionSet],
)
def test_yb_2_warning(input_set: VaspInputSet) -> None:
    # https://github.com/materialsproject/pymatgen/pull/2972

    structure = Structure(
        lattice=Lattice.cubic(5),
        species=("Yb", "O"),
        coords=((0, 0, 0), (0.5, 0.5, 0.5)),
    )

    with pytest.warns(BadInputSetWarning) as record:
        input_set(structure)

    expected = "The structure contains Ytterbium (Yb) and this InputSet uses the Yb_2"
    assert expected in str(record[0].message)


class TestSetChangeCheck(PymatgenTest):
    def test_sets_changed(self):
        msg = (
            "WARNING! These tests will fail when you change an input set. They are included "
            "as a sanity check. When changing an input set, please make sure to "
            "notify the users for that set. For sets starting with 'MVL' this is @shyuep, for "
            "sets starting with 'MP' this is @shyuep and @mkhorton. For sets starting with 'MatPES' "
            "this is @shyuep and @janosh."
        )

        input_sets = glob(f"{MODULE_DIR}/*.yaml")
        hashes = {}
        for input_set in input_sets:
            with open(input_set) as file:
                text = file.read().encode("utf-8")
                name = os.path.basename(input_set)
                hashes[name] = hashlib.sha1(text).hexdigest()

        known_hashes = {
            "MatPESStaticSet.yaml": "8edecff2bbd1932c53159f56a8e6340e900aaa2f",
            "MITRelaxSet.yaml": "1a0970f8cad9417ec810f7ab349dc854eaa67010",
            "MPAbsorptionSet.yaml": "5931e1cb3cf8ba809b3d4f4a5960d728c682adf1",
            "MPHSERelaxSet.yaml": "0d0d96a620461071cfd416ec9d5d6a8d2dfd0855",
            "MPRelaxSet.yaml": "f2949cdc5dc8cd0bee6d39a5df0d6a6b7c144821",
            "MPSCANRelaxSet.yaml": "2d31ee637cb5d4d96f2e0aba3772a52cbcceb348",
            "MVLGWSet.yaml": "104ae93c3b3be19a13b0ee46ebdd0f40ceb96597",
            "MVLRelax52Set.yaml": "4cfc6b1bd0548e45da3bde4a9c65b3249da13ecd",
            "PBE54Base.yaml": "ec317781a7f344beb54c17a228db790c0eb49282",
            "PBE64Base.yaml": "480c41c2448cb25706181de268090618e282b264",
            "VASPIncarBase.yaml": "19762515f8deefb970f2968fca48a0d67f7964d4",
            "vdW_parameters.yaml": "04bb09bb563d159565bcceac6a11e8bdf0152b79",
        }

        for input_set in hashes:
            assert hashes[input_set] == known_hashes[input_set], f"{input_set=}\n{msg}"


class TestDictSet(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        cls.structure = Structure.from_file(filepath)

    def test_as_dict(self):
        # https://github.com/materialsproject/pymatgen/pull/3031
        dict_set = DictSet(self.structure, config_dict={"INCAR": {}}, user_potcar_functional="PBE_54")
        assert {*dict_set.as_dict()} >= {
            "structure",
            "config_dict",
            "user_incar_settings",
            "user_kpoints_settings",
            "user_potcar_settings",
        }
        assert dict_set.potcar_functional == dict_set.user_potcar_functional


class TestMITMPRelaxSet(PymatgenTest):
    @classmethod
    def setUpClass(cls):
        cls.set = MITRelaxSet
        cls.mp_set = MPRelaxSet
        cls.monkeypatch = MonkeyPatch()

        filepath = f"{TEST_FILES_DIR}/POSCAR"
        cls.structure = Structure.from_file(filepath)
        cls.coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        cls.lattice = Lattice([[3.8401979337, 0, 0], [1.9200989668, 3.3257101909, 0], [0, -2.2171384943, 3.1355090603]])

        cls.mit_set = cls.set(cls.structure)
        cls.mit_set_unsorted = cls.set(cls.structure, sort_structure=False)
        cls.mp_set = MPRelaxSet(cls.structure)

    def test_no_structure_init(self):
        # basic test of initialization with no structure.
        vis = MPRelaxSet()
        assert vis.as_dict()["structure"] is None
        with pytest.raises(RuntimeError, match="No structure is associated with the input set!"):
            _ = vis.incar
        vis.structure = self.structure
        assert vis.incar["LDAUU"] == [5.3, 0, 0]
        assert vis.as_dict()["structure"] is not None

    def test_metal_check(self):
        structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Cu"], [[0, 0, 0]])

        with pytest.warns(
            BadInputSetWarning,
            match="Relaxation of likely metal with ISMEAR < 1 detected. See VASP recommendations on ISMEAR for metals.",
        ) as warns:
            vis = self.set(structure)
            _ = vis.incar
        assert len(warns) == 1

    def test_poscar(self):
        structure = Structure(self.lattice, ["Fe", "Mn"], self.coords)
        mit_param_set = self.set(structure, sort_structure=False)
        s_unsorted = mit_param_set.poscar.structure
        mit_param_set = self.set(structure, sort_structure=True)
        s_sorted = mit_param_set.poscar.structure
        assert s_unsorted[0].specie.symbol == "Fe"
        assert s_sorted[0].specie.symbol == "Mn"

    def test_potcar_symbols(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75], [0.75, 0.25, 0.75]]
        lattice = Lattice([[3.8401979337, 0, 0], [1.9200989668, 3.3257101909, 0], [0, -2.2171384943, 3.1355090603]])
        structure = Structure(lattice, ["P", "Fe", "O"], coords)
        mit_param_set = self.set(structure)
        syms = mit_param_set.potcar_symbols
        assert syms == ["Fe", "P", "O"]
        param_set = MPRelaxSet(structure, sort_structure=False)
        syms = param_set.potcar_symbols
        assert syms == ["P", "Fe_pv", "O"]

    def test_potcar_validation(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        # Use pytest's monkeypatch to temporarily point pymatgen to a directory
        # containing the wrong POTCARs (LDA potcars in a PBE directory)
        with self.monkeypatch.context() as m:
            m.setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(f"{TEST_FILES_DIR}/wrong_potcars"))
            with pytest.warns(BadInputSetWarning, match="not known by pymatgen"):
                _ = self.set(structure).potcar

    def test_potcar_special_defaults(self):
        # https://github.com/materialsproject/pymatgen/pull/3022
        for user_potcar_settings in [{"Fe": "Fe_pv"}, {"W": "W_pv"}, None]:
            for species in [("W", "W"), ("Fe", "W"), ("Fe", "Fe")]:
                struct = Structure(lattice=Lattice.cubic(3), species=species, coords=[[0, 0, 0], [0.5, 0.5, 0.5]])
                relax_set = MPRelaxSet(
                    structure=struct, user_potcar_functional="PBE_54", user_potcar_settings=user_potcar_settings
                )
                expected = {  # noqa: SIM222
                    **({"W": "W_sv"} if "W" in struct.symbol_set else {}),
                    **(user_potcar_settings or {}),
                } or None
                assert relax_set.user_potcar_settings == expected

    @skip_if_no_psp_dir
    def test_lda_potcar(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        p = self.set(structure, user_potcar_functional="LDA").potcar
        assert p.functional == "LDA"

    @skip_if_no_psp_dir
    def test_nelect(self):
        coords = [[0] * 3, [0.5] * 3, [0.75] * 3]
        lattice = Lattice.cubic(4)
        struct = Structure(lattice, ["Si", "Si", "Fe"], coords)
        assert self.set(struct).nelect == 16
        assert MPRelaxSet(struct).nelect == 22

        # Expect same answer when oxidation states are present. Was a bug previously.
        oxi_struct = Structure(lattice, ["Si4+", "Si4+", "Fe2+"], coords)
        assert self.set(oxi_struct).nelect == 16
        assert MPRelaxSet(oxi_struct).nelect == 22

        # disordered structure are not supported
        disordered = Structure.from_spacegroup("Im-3m", Lattice.cubic(3), [Composition("Fe0.5Mn0.5")], [[0, 0, 0]])
        with pytest.raises(
            ValueError, match="Disordered structure with partial occupancies cannot be converted into POSCAR"
        ):
            _ = self.set(disordered).nelect

    @skip_if_no_psp_dir
    def test_estimate_nbands(self):
        # estimate_nbands is a function of n_elect, n_ions, magmom, non-collinearity of magnetism, and n_par
        coords = [[0] * 3, [0.5] * 3, [0.75] * 3]
        lattice = Lattice.cubic(4)

        # pure Si
        struct = Structure(lattice, ["Si", "Si", "Si"], coords)
        assert self.set(struct).estimate_nbands() == 11
        assert MPRelaxSet(struct).estimate_nbands() == 11

        # Si + Fe
        struct = Structure(lattice, ["Si", "Si", "Fe"], coords)
        assert self.set(struct).estimate_nbands() == 15
        assert MPRelaxSet(struct).estimate_nbands() == 18

        # Si + Fe with NPAR = 4
        uis = {"NPAR": 4}
        assert self.set(struct, user_incar_settings=uis).estimate_nbands() == approx(16)
        assert MPRelaxSet(struct, user_incar_settings=uis).estimate_nbands() == approx(20)

        # Si + Fe with noncollinear magnetism turned on
        uis = {"LNONCOLLINEAR": True}
        assert self.set(struct, user_incar_settings=uis).estimate_nbands() == approx(30)
        assert MPRelaxSet(struct, user_incar_settings=uis).estimate_nbands() == approx(36)

    @skip_if_no_psp_dir
    def test_get_incar(self):
        incar = self.mp_set.incar

        assert incar["LDAUU"] == [5.3, 0, 0]
        assert incar["EDIFF"] == approx(0.0012)

        incar = self.mit_set.incar
        assert incar["LDAUU"] == [4.0, 0, 0]
        assert incar["EDIFF"] == approx(1e-5)

        si = 14
        coords = []
        coords.extend((np.array([0, 0, 0]), np.array([0.75, 0.5, 0.75])))

        # Silicon structure for testing.
        lattice = Lattice(
            [[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        )
        struct = Structure(lattice, [si, si], coords)
        incar = MPRelaxSet(struct).incar
        assert "LDAU" not in incar

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice(
            [[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        )
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        incar = MPRelaxSet(struct).incar
        assert "LDAU" not in incar

        # check fluorides
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPRelaxSet(struct).incar
        assert incar["LDAUU"] == [5.3, 0]
        assert incar["MAGMOM"] == [5, 0.6]

        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = self.set(struct).incar
        assert incar["LDAUU"] == [4.0, 0]

        # This seems counterintuitive at first, but even if the prior INCAR has a MAGMOM flag,
        # because the structure has no site properties, the default MAGMOM is assigned from the
        # config dictionary.
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPStaticSet(struct, prev_incar=f"{TEST_FILES_DIR}/INCAR").incar
        assert incar["MAGMOM"] == [5, 0.6]

        # Make sure this works with species.
        struct = Structure(lattice, ["Fe2+", "O2-"], coords)
        incar = MPRelaxSet(struct).incar
        assert incar["LDAUU"] == [5.3, 0]

        struct = Structure(lattice, ["Fe", "Mn"], coords, site_properties={"magmom": (5.2, -4.5)})
        incar = MPRelaxSet(struct).incar
        assert incar["MAGMOM"] == [-4.5, 5.2]

        incar = self.set(struct, sort_structure=False).incar
        assert incar["MAGMOM"] == [5.2, -4.5]

        struct = Structure(lattice, [Species("Fe2+", spin=4.1), "Mn"], coords)
        incar = MPRelaxSet(struct).incar
        assert incar["MAGMOM"] == [5, 4.1]

        struct = Structure(lattice, ["Mn3+", "Mn4+"], coords)
        incar = self.set(struct).incar
        assert incar["MAGMOM"] == [4, 3]

        user_set = MPRelaxSet(struct, user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}})
        assert user_set.incar["MAGMOM"] == [100, 0.6]

        no_encut_set = MPRelaxSet(struct, user_incar_settings={"ENCUT": None})
        assert "ENCUT" not in no_encut_set.incar

        # sulfide vs sulfate test
        coords = [[0, 0, 0], [0.75, 0.5, 0.75], [0.25, 0.5, 0]]

        struct = Structure(lattice, ["Fe", "Fe", "S"], coords)
        incar = self.set(struct).incar
        assert incar["LDAUU"] == [1.9, 0]

        # Make sure MP sulfides are ok.
        assert "LDAUU" not in MPRelaxSet(struct).incar

        struct = Structure(lattice, ["Fe", "S", "O"], coords)
        incar = self.set(struct).incar
        assert incar["LDAUU"] == [4.0, 0, 0]

        # Make sure MP sulfates are ok.
        assert MPRelaxSet(struct).incar["LDAUU"] == [5.3, 0, 0]

        # test for default LDAUU value
        user_set_ldauu_fallback = MPRelaxSet(struct, user_incar_settings={"LDAUU": {"Fe": 5.0, "S": 0}})
        assert user_set_ldauu_fallback.incar["LDAUU"] == [5.0, 0, 0]

        # Expected to be oxide (O is the most electronegative atom)
        struct = Structure(lattice, ["Fe", "O", "S"], coords)
        incar = self.set(struct).incar
        assert incar["LDAUU"] == [4.0, 0, 0]

        # Expected to be chloride (Cl is the most electronegative atom)
        struct = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = self.set(struct, user_incar_settings={"LDAU": True}).incar
        assert "LDAUU" not in incar  # LDAU = False

        # User set a compound to be sulfide by specifying values of "LDAUL" etc.
        struct = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = self.set(
            struct,
            user_incar_settings={"LDAU": True, "LDAUL": {"Fe": 3}, "LDAUU": {"Fe": 1.8}},
        ).incar
        assert incar["LDAUL"] == [3.0, 0, 0]
        assert incar["LDAUU"] == [1.8, 0, 0]

        # test that van-der-Waals parameters are parsed correctly
        incar = self.set(struct, vdw="optB86b").incar
        assert incar["GGA"] == "Mk"
        assert incar["LUSE_VDW"]
        assert incar["PARAM1"] == 0.1234

        # Test that NELECT is updated when a charge is present
        si = 14
        coords = []
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        # Silicon structure for testing.
        lattice = Lattice(
            [[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        )
        struct = Structure(lattice, [si, si], coords, charge=1)
        mpr = MPRelaxSet(struct, use_structure_charge=True)
        assert mpr.incar["NELECT"] == 7, "NELECT not properly set for nonzero charge"

        # test that NELECT does not get set when use_structure_charge = False
        mpr = MPRelaxSet(struct, use_structure_charge=False)
        assert "NELECT" not in mpr.incar, "NELECT should not be set when use_structure_charge is False"

        struct = Structure(lattice, ["Co", "O"], coords)
        mpr = MPRelaxSet(struct)
        assert mpr.incar["MAGMOM"] == [0.6, 0.6]
        struct = Structure(lattice, ["Co4+", "O"], coords)
        mpr = MPRelaxSet(struct)
        assert mpr.incar["MAGMOM"] == [1, 0.6]

        # test passing user_incar_settings and user_kpoint_settings of None
        for set_cls in [MPRelaxSet, MPStaticSet, MPNonSCFSet]:
            mp_set = set_cls(struct, user_incar_settings=None, user_kpoints_settings=None)
            assert mp_set.kpoints is not None
            assert mp_set.incar is not None

    def test_get_kpoints(self):
        kpoints = MPRelaxSet(self.structure).kpoints
        assert kpoints.kpts == [[2, 4, 5]]
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints = MPRelaxSet(self.structure, user_kpoints_settings={"reciprocal_density": 1000}).kpoints
        assert kpoints.kpts == [[6, 10, 13]]
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints_obj = Kpoints(kpts=[[3, 3, 3]])
        kpoints_return = MPRelaxSet(self.structure, user_kpoints_settings=kpoints_obj).kpoints
        assert kpoints_return.kpts == [[3, 3, 3]]

        kpoints = self.mit_set.kpoints
        assert kpoints.kpts == [[25]]
        assert kpoints.style == Kpoints.supported_modes.Automatic

        recip_param_set = MPRelaxSet(self.structure, force_gamma=True)
        recip_param_set.kpoints_settings = {"reciprocal_density": 40}
        kpoints = recip_param_set.kpoints
        assert kpoints.kpts == [[2, 4, 5]]
        assert kpoints.style == Kpoints.supported_modes.Gamma

    @skip_if_no_psp_dir
    def test_get_vasp_input(self):
        d = self.mit_set.get_vasp_input()
        assert d["INCAR"]["ISMEAR"] == -5
        struct = self.structure.copy()
        struct.make_supercell(4)
        paramset = MPRelaxSet(struct)
        d = paramset.get_vasp_input()
        assert d["INCAR"]["ISMEAR"] == 0

    def test_mp_metal_relax_set(self):
        mp_metal_set = MPMetalRelaxSet(self.get_structure("Sn"))
        incar = mp_metal_set.incar
        assert incar["ISMEAR"] == 1
        assert incar["SIGMA"] == 0.2
        kpoints = mp_metal_set.kpoints
        assert_allclose(kpoints.kpts[0], [5, 5, 5])

    def test_as_from_dict(self):
        mit_set = self.set(self.structure)
        mp_set = MPRelaxSet(self.structure)
        mp_user_set = MPRelaxSet(
            self.structure,
            user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}},
        )

        dct = mit_set.as_dict()
        val = dec.process_decoded(dct)
        assert val._config_dict["INCAR"]["LDAUU"]["O"]["Fe"] == 4

        dct = mp_set.as_dict()
        val = dec.process_decoded(dct)
        assert val._config_dict["INCAR"]["LDAUU"]["O"]["Fe"] == 5.3

        dct = mp_user_set.as_dict()
        val = dec.process_decoded(dct)
        # assert isinstance(val, MPVaspInputSet)
        assert val.user_incar_settings["MAGMOM"] == {"Fe": 10, "S": -5, "Mn3+": 100}

    def test_hubbard_off_and_ediff_override(self):
        input_set = MPRelaxSet(self.structure, user_incar_settings={"LDAU": False, "EDIFF": 1e-10})
        assert "LDAUU" not in input_set.incar
        assert input_set.incar["EDIFF"] == 1e-10
        # after testing, we have determined LMAXMIX should still be 4 for d-block
        # even if U is turned off (thanks Andrew Rosen for reporting)
        assert input_set.incar["LMAXMIX"] == 4

    def test_incar_lmaxmix(self):
        # https://github.com/materialsproject/pymatgen/issues/3040

        # structure containing neither f- nor d-electrons
        structure_f = self.get_structure("Si")
        assert "LMAXMIX" not in MPRelaxSet(structure_f).incar

        # structure containing d-electrons but no f-electrons
        structure_d = self.get_structure("LiFePO4")
        assert MPRelaxSet(structure_d).incar["LMAXMIX"] == 4

        # structure containing f-electrons but no d-electrons
        structure_f = structure_d.copy()
        structure_f.replace_species({"Fe": "La"})
        assert MPRelaxSet(structure_f).incar["LMAXMIX"] == 6

        # structure containing f- and d-electrons
        structure_f_and_d = structure_d.copy()
        structure_f_and_d.replace_species({"P": "Ce"})
        assert MPRelaxSet(structure_f_and_d).incar["LMAXMIX"] == 6

        # explicit LMAXMIX in settings overrides automatic selection
        structure_override = Structure(
            lattice=Lattice.cubic(3),
            species=("Fe", "Fe"),
            coords=((0, 0, 0), (0.5, 0.5, 0.5)),
        )
        set_override = MPRelaxSet(structure_override, user_incar_settings={"LMAXMIX": 3})
        assert set_override.incar["LMAXMIX"] == 3

    @skip_if_no_psp_dir
    def test_write_input(self):
        vasp_files = {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}
        self.mit_set.write_input(self.tmp_path)
        assert {*os.listdir(self.tmp_path)} == vasp_files

        self.mit_set.write_input(self.tmp_path, include_cif=True)
        assert {*os.listdir(self.tmp_path)} == {*vasp_files, "Fe4P4O16.cif"}

    @skip_if_no_psp_dir
    def test_write_input_potcar_spec(self):
        self.mit_set.write_input(self.tmp_path, potcar_spec=True)
        assert {*os.listdir(self.tmp_path)} == {"INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"}

    @skip_if_no_psp_dir
    def test_user_potcar_settings(self):
        vis = MPRelaxSet(self.structure, user_potcar_settings={"Fe": "Fe"})
        potcar = vis.potcar
        assert potcar.symbols == ["Fe", "P", "O"]

    @skip_if_no_psp_dir
    def test_valid_magmom_struct(self):
        # First test the helper function
        struct = self.structure.copy()
        get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="v")
        props = [site.properties for site in struct]
        assert props == [{"magmom": [1.0, 1.0, 1.0]}] * len(props)

        struct = self.structure.copy()
        get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="s")
        props = [site.properties for site in struct]
        assert props == [{"magmom": 1.0}] * len(props)
        struct.insert(0, "Li", [0, 0, 0])
        get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="a")
        props = [site.properties for site in struct]
        assert props == [{"magmom": 1.0}] * len(props)

        struct = self.structure.copy()
        get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="v")
        struct.insert(0, "Li", [0, 0, 0], properties={"magmom": 10.0})
        with pytest.raises(TypeError) as exc:
            get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="a")
        assert "Magmom type conflict" in str(exc.value)

        # Test the behavior of MPRelaxSet to automatically fill in the missing magmom
        struct = self.structure.copy()
        get_valid_magmom_struct(structure=struct, inplace=True, spin_mode="s")
        struct.insert(0, "Li", [0, 0, 0])

        # vis = MPRelaxSet(struct, user_potcar_settings={"Fe": "Fe"}, validate_magmom=False)
        # with pytest.raises(TypeError, match=r"float\(\) argument must be a string or a (real )?number,
        #                    not 'NoneType'"):
        #     vis.get_vasp_input()

        vis = MPRelaxSet(struct, user_potcar_settings={"Fe": "Fe"}, validate_magmom=True)
        assert vis.get_vasp_input()["INCAR"]["MAGMOM"] == [1.0] * len(struct)

        # Test the behavior of constraining the net magnetic moment with a non-integer
        struct = self.structure.copy()
        with pytest.warns(UserWarning, match="constrain_total_magmom"):
            vis = MPRelaxSet(
                struct,
                user_incar_settings={"MAGMOM": {"Fe": 5.1}},
                user_potcar_settings={"Fe": "Fe"},
                constrain_total_magmom=True,
            )
            vis.incar.items()


class TestMPStaticSet(PymatgenTest):
    def setUp(self):
        self.set = MPStaticSet

    def test_init(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["NSW"] == 0
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600
        assert vis.kpoints.style == Kpoints.supported_modes.Monkhorst

        # Check as from dict.
        vis = self.set.from_dict(vis.as_dict())
        assert vis.incar["NSW"] == 0
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600
        assert vis.kpoints.style == Kpoints.supported_modes.Monkhorst

        non_prev_vis = self.set(vis.structure, user_incar_settings={"LORBIT": 12, "LWAVE": True})
        assert non_prev_vis.incar["NSW"] == 0
        # Check that the ENCUT and Kpoints style has NOT been inherited.
        assert non_prev_vis.incar["ENCUT"] == 520
        # Check that user incar settings are applied.
        assert non_prev_vis.incar["LORBIT"] == 12
        assert non_prev_vis.incar["LWAVE"]

        assert non_prev_vis.kpoints.style == Kpoints.supported_modes.Gamma
        v2 = self.set.from_dict(non_prev_vis.as_dict())
        assert v2.incar["ENCUT"] == 520
        # Check that user incar settings are applied.
        assert v2.incar["LORBIT"] == 12
        leps_vis = self.set.from_prev_calc(prev_calc_dir=prev_run, lepsilon=True)
        assert leps_vis.incar["LEPSILON"]
        assert leps_vis.incar["IBRION"] == 8
        assert leps_vis.incar["EDIFF"] == 1e-5
        assert "NPAR" not in leps_vis.incar
        assert "NSW" not in leps_vis.incar
        assert non_prev_vis.kpoints.kpts == [[11, 10, 10]]
        non_prev_vis = self.set(vis.structure, reciprocal_density=200)
        assert non_prev_vis.kpoints.kpts == [[14, 12, 12]]
        # Check LCALCPOL flag
        lcalc_pol_vis = self.set.from_prev_calc(prev_calc_dir=prev_run, lcalcpol=True)
        assert lcalc_pol_vis.incar["LCALCPOL"]

        # Check warning if LASPH is set to False for meta-GGAs/hybrids/+U/vdW
        with pytest.warns(BadInputSetWarning, match=r"LASPH"):
            vis = self.set(vis.structure, user_incar_settings={"METAGGA": "M06L", "LASPH": False})
            vis.incar.items()
        with pytest.warns(BadInputSetWarning, match=r"LASPH"):
            vis = self.set(vis.structure, user_incar_settings={"LHFCALC": True, "LASPH": False})
            vis.incar.items()
        with pytest.warns(BadInputSetWarning, match=r"LASPH"):
            vis = self.set(vis.structure, user_incar_settings={"LUSE_VDW": True, "LASPH": False})
            vis.incar.items()
        with pytest.warns(BadInputSetWarning, match=r"LASPH"):
            dummy_struc = Structure(
                lattice=[[0, 2, 2], [2, 0, 2], [2, 2, 0]],
                species=["Fe", "O"],
                coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
            vis = self.set(dummy_struc, user_incar_settings={"LDAU": True, "LASPH": False})
            vis.incar.items()

    def test_user_incar_kspacing(self):
        # Make sure user KSPACING settings properly overrides KPOINTS.
        si = self.get_structure("Si")
        vis = self.set(si, user_incar_settings={"KSPACING": 0.22})
        assert vis.incar["KSPACING"] == 0.22
        assert vis.kpoints is None

    def test_kspacing_override(self):
        # If KSPACING is set and user_kpoints_settings are given,
        # make sure the user_kpoints_settings override KSPACING
        si = self.get_structure("Si")
        vis = self.set(
            si,
            user_incar_settings={"KSPACING": 0.22},
            user_kpoints_settings={"reciprocal_density": 1000},
        )
        assert vis.incar.get("KSPACING") is None
        assert isinstance(vis.kpoints, Kpoints)

    def test_override_from_prev_calc(self):
        # test override_from_prev
        prev_run = f"{TEST_FILES_DIR}/relaxation"

        vis = self.set(dummy_structure)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["NSW"] == 0
        assert vis.incar["ENCUT"] == 600
        assert vis.kpoints.style == Kpoints.supported_modes.Monkhorst

        # Check LCALCPOL flag
        lcalcpol_vis = self.set(dummy_structure, lcalcpol=True)
        lcalcpol_vis = lcalcpol_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert lcalcpol_vis.incar["LCALCPOL"]

    def test_standardize_structure(self):
        sga = SpacegroupAnalyzer(self.get_structure("Si"))
        original_structure = sga.get_conventional_standard_structure()
        sm = StructureMatcher(primitive_cell=False, scale=False)

        vis = self.set(original_structure)
        assert sm.fit(vis.structure, original_structure)

        vis = self.set(original_structure, standardize=True)
        assert not sm.fit(vis.structure, original_structure)

    def test_write_input_zipped(self):
        vis = self.set(self.get_structure("Si"))
        vis.write_input(output_dir=self.tmp_path, potcar_spec=True, zip_output=True)

        assert os.path.isfile(f"{self.tmp_path}/MPStaticSet.zip")
        with ZipFile(f"{self.tmp_path}/MPStaticSet.zip", "r") as zip_file:
            contents = zip_file.namelist()
            assert set(contents).issuperset({"INCAR", "POSCAR", "POTCAR.spec", "KPOINTS"})
            spec = zip_file.open("POTCAR.spec", "r").read().decode()
            assert spec == "Si"

    def test_grid_size_from_struct(self):
        # TODO grab a bunch_of_calculations store as a list of tuples
        # (structure, ngx, ..., ngxf, ...) where all the grid size values are generated by vasp
        # check that the code produces the same grid sizes
        fname = f"{TEST_FILES_DIR}/grid_data_files/vasp_inputs_for_grid_check.json"
        parsed_vasp_data = loadfn(fname)
        for tt in parsed_vasp_data:
            ng = [tt["input"]["parameters"][ik] for ik in ["NGX", "NGY", "NGZ"]]
            ngf = [tt["input"]["parameters"][ik] for ik in ["NGXF", "NGYF", "NGZF"]]
            struct = tt["input"]["structure"]
            static_set = self.set(struct)
            matched = static_set.calculate_ng() == (ng, ngf)
            assert matched

        assert static_set.calculate_ng() == ([30, 48, 50], [60, 96, 100])
        # test `custom_encut` kwarg for final structure in above test using
        # an (obviously fictitious) custom encut.
        assert static_set.calculate_ng(custom_encut=2000) == ([56, 96, 96], [112, 192, 192])

        assert static_set.calculate_ng() == ([30, 48, 50], [60, 96, 100])
        # test `custom_prec` kwarg for final structure in above test using "NORMAL".
        assert static_set.calculate_ng(custom_prec="NORMAL") == ([24, 36, 40], [48, 72, 80])


class TestMatPESStaticSet(PymatgenTest):
    def setUp(self):
        self.struct = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR")
        self.prev_incar = Incar.from_file(f"{TEST_FILES_DIR}/INCAR")

    def test_default(self):
        input_set = MatPESStaticSet(self.struct)
        incar = input_set.incar
        assert incar["ALGO"] == "Normal"
        assert incar["EDIFF"] == 1.0e-05
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["GGA"] == "Pe"
        assert incar["ISMEAR"] == 0
        assert incar["ISPIN"] == 2
        assert incar["KSPACING"] == 0.22
        assert incar["LAECHG"]
        assert incar["LASPH"]
        assert incar["LCHARG"]
        assert incar["LMIXTAU"]
        assert incar.get("LDAU") is None
        assert incar["LORBIT"] == 11
        assert not incar["LREAL"]
        assert not incar["LWAVE"]
        assert incar["NELM"] == 200
        assert incar["NSW"] == 0
        assert incar["PREC"] == "Accurate"
        assert incar["SIGMA"] == 0.05
        assert incar["LMAXMIX"] == 6
        assert input_set.potcar_symbols == ["Fe_pv", "P", "O"]
        assert input_set.potcar_symbols == ["Fe_pv", "P", "O"]
        assert input_set.potcar_functional == "PBE_64"
        assert input_set.kpoints is None
        # test POTCAR files are default PBE_64 PSPs and functional
        # only runs if POTCAR files to compare against are available
        if os.path.isdir(f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_64"):
            assert str(input_set.potcar[0]) == str(PotcarSingle.from_symbol_and_functional("Fe_pv", "PBE_64"))

    def test_with_prev_incar(self):
        default_prev = MatPESStaticSet(structure=self.struct, prev_incar=self.prev_incar)
        incar = default_prev.incar
        # check that prev_incar is used.
        assert incar["NPAR"] == 8
        assert incar["LMAXMIX"] == 4
        # test some incar parameters from prev_incar are not inherited
        for key in "ISPIND LSCALU ISIF SYSTEM HFSCREEN NSIM ENCUTFOCK NKRED LPLANE TIME LHFCALC".split():
            assert key not in incar, f"{key=} should be absent, got {incar[key]=}"
        # test if default in MatPESStaticSet is prioritized.
        assert incar["ALGO"] == "Normal"
        assert incar["EDIFF"] == 1.0e-05
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["GGA"] == "Pe"
        assert incar["ISMEAR"] == 0
        assert incar["ISPIN"] == 2
        assert incar["KSPACING"] == 0.22
        assert incar["LAECHG"]
        assert incar["LASPH"]
        assert incar["LCHARG"]
        assert incar["LMIXTAU"]
        assert incar.get("LDAU") is None
        assert incar["LORBIT"] == 11
        assert incar["LREAL"] is False
        assert incar["LWAVE"] is False
        assert incar["NELM"] == 200
        assert incar["NSW"] == 0
        assert incar["PREC"] == "Accurate"
        assert incar["SIGMA"] == 0.05
        # test POTCAR files are default PBE_64 PSPs and functional
        assert default_prev.potcar_symbols == ["Fe_pv", "P", "O"]
        assert default_prev.potcar_functional == "PBE_64"
        assert default_prev.kpoints is None

    def test_r2scan(self):
        scan = MatPESStaticSet(self.struct, xc_functional="R2SCAN")
        incar_scan = scan.incar
        assert incar_scan["METAGGA"] == "R2scan"
        assert incar_scan.get("GGA") is None
        assert incar_scan["ALGO"] == "All"
        assert incar_scan.get("LDAU") is None
        # test POTCAR files are default PBE_64 PSPs and functional
        assert scan.potcar_symbols == ["Fe_pv", "P", "O"]
        assert scan.potcar_functional == "PBE_64"
        assert scan.kpoints is None

    def test_default_u(self):
        default_u = MatPESStaticSet(self.struct, xc_functional="PBE+U")
        incar_u = default_u.incar
        assert incar_u["LDAU"] is True
        assert incar_u["GGA"] == "Pe"
        assert incar_u["ALGO"] == "Normal"
        # test POTCAR files are default PBE_64 PSPs and functional
        assert incar_u["LDAUU"] == [5.3, 0, 0]
        assert default_u.potcar_symbols == ["Fe_pv", "P", "O"]
        assert default_u.potcar_functional == "PBE_64"
        assert default_u.kpoints is None

    def test_functionals(self):
        xc_functional = "LDA"
        with pytest.raises(
            ValueError, match=f"Unrecognized {xc_functional=}. Supported exchange-correlation functionals are "
        ):
            MatPESStaticSet(self.struct, xc_functional=xc_functional)

        with pytest.warns(UserWarning, match="inconsistent with the recommended PBE_64"):
            diff_potcar = MatPESStaticSet(self.struct, user_potcar_functional="PBE")
            assert str(diff_potcar.potcar[0]) == str(PotcarSingle.from_symbol_and_functional("Fe_pv", "PBE"))

    def test_from_prev_calc(self):
        vis = MatPESStaticSet.from_prev_calc(f"{TEST_FILES_DIR}/relaxation")
        incar = vis.incar
        assert incar["GGA"] == "Pe"
        assert incar["ALGO"] == "Normal"
        assert vis.potcar_symbols == ["Li_sv"]
        assert vis.kpoints is None


class TestMPNonSCFSet(PymatgenTest):
    def setUp(self):
        self.set = MPNonSCFSet

    @skip_if_no_psp_dir
    def test_init(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        # check boltztrap mode

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Boltztrap")
        assert vis.incar["ISMEAR"] == 0

        # check uniform mode
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform")
        assert vis.incar["ISMEAR"] == -5
        assert vis.incar["ISYM"] == 2

        # check uniform mode with automatic nedos
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform", nedos=0)
        assert vis.incar["NEDOS"] == 12217

        # test line mode
        vis = self.set.from_prev_calc(
            prev_calc_dir=prev_run,
            mode="Line",
            copy_chgcar=False,
            user_incar_settings={"SIGMA": 0.025},
        )

        assert vis.incar["NSW"] == 0
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600
        # Check that the user_incar_settings works
        assert vis.incar["SIGMA"] == 0.025
        assert vis.kpoints.style == Kpoints.supported_modes.Reciprocal

        # Check as from dict.
        vis = self.set.from_dict(vis.as_dict())
        assert vis.incar["NSW"] == 0
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600
        assert vis.kpoints.style == Kpoints.supported_modes.Reciprocal

        vis.write_input(self.tmp_path)
        assert not os.path.isfile(f"{self.tmp_path}/CHGCAR")

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Line", copy_chgcar=True)
        # check ISMEAR set correctly for line mode
        assert vis.incar["ISMEAR"] == 0
        vis.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/CHGCAR")
        os.remove(f"{self.tmp_path}/CHGCAR")  # needed for next assert

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, standardize=True, mode="Line", copy_chgcar=True)
        vis.write_input(self.tmp_path)
        assert not os.path.isfile(f"{self.tmp_path}/CHGCAR")

    @skip_if_no_psp_dir
    def test_override_from_prev(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"

        # test override_from_prev
        vis = self.set(dummy_structure, mode="Boltztrap")
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["ISMEAR"] == 0

        vis = self.set(dummy_structure, mode="Uniform")
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["ISMEAR"] == -5
        assert vis.incar["ISYM"] == 2

        vis = self.set(dummy_structure, mode="Uniform", nedos=0)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["NEDOS"] == 12217

        # test line mode
        vis = self.set(
            dummy_structure,
            mode="Line",
            copy_chgcar=False,
            user_incar_settings={"SIGMA": 0.025},
        )
        vis.override_from_prev_calc(prev_calc_dir=prev_run)

        assert vis.incar["NSW"] == 0
        assert vis.incar["ENCUT"] == 600
        assert vis.incar["SIGMA"] == 0.025
        assert vis.kpoints.style == Kpoints.supported_modes.Reciprocal

        vis = self.set(dummy_structure, mode="Line", copy_chgcar=True)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["ISMEAR"] == 0
        vis.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/CHGCAR")
        os.remove(f"{self.tmp_path}/CHGCAR")  # needed for next assert

        vis = self.set(dummy_structure, standardize=True, mode="Line", copy_chgcar=True)
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        vis.write_input(self.tmp_path)
        assert not os.path.isfile(f"{self.tmp_path}/CHGCAR")

    def test_kpoints(self):
        # test k-points are generated in the correct format
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape == (1, 3)

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Line", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape != (1, 3)

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Boltztrap", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape != (1, 3)

    def test_optics(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        vis = self.set.from_prev_calc(
            prev_calc_dir=prev_run,
            copy_chgcar=False,
            optics=True,
            mode="Uniform",
            nedos=2001,
        )

        assert vis.incar["NSW"] == 0
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600

        # check NEDOS and ISMEAR set correctly
        assert vis.incar["NEDOS"] == 2001
        assert vis.incar["ISMEAR"] == -5
        assert vis.incar["ISYM"] == 2

        assert vis.incar["LOPTICS"]
        assert vis.kpoints.style == Kpoints.supported_modes.Gamma

    def test_user_kpoint_override(self):
        # default kpoints style is reciprocal, try setting to gamma
        user_kpoints_override = Kpoints(style=Kpoints.supported_modes.Gamma, kpts=((1, 1, 1),))

        prev_run = f"{TEST_FILES_DIR}/relaxation"
        vis = self.set.from_prev_calc(
            prev_calc_dir=prev_run,
            copy_chgcar=False,
            optics=True,
            mode="Uniform",
            nedos=2001,
            user_kpoints_settings=user_kpoints_override,
        )
        assert vis.kpoints.style == Kpoints.supported_modes.Gamma


class TestMagmomLdau(PymatgenTest):
    def test_structure_from_prev_run(self):
        vrun = Vasprun(f"{TEST_FILES_DIR}/vasprun.xml.magmom_ldau")
        structure = vrun.final_structure
        poscar = Poscar(structure)
        structure_decorated = get_structure_from_prev_run(vrun)
        ldau_ans = {"LDAUU": [5.3, 0.0], "LDAUL": [2, 0], "LDAUJ": [0.0, 0.0]}
        magmom_ans = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
        ldau_dict = {}
        for key in ("LDAUU", "LDAUJ", "LDAUL"):
            if hasattr(structure_decorated[0], key.lower()):
                m = {site.specie.symbol: getattr(site, key.lower()) for site in structure_decorated}
                ldau_dict[key] = [m[sym] for sym in poscar.site_symbols]
        magmom = [site.magmom for site in structure_decorated]
        assert ldau_dict == ldau_ans
        assert magmom == magmom_ans

    def test_ln_magmom(self):
        yaml_path = MODULE_DIR / "VASPIncarBase.yaml"
        magmom_setting = loadfn(yaml_path)["INCAR"]["MAGMOM"]
        structure = Structure.from_file(f"{TEST_FILES_DIR}/La4Fe4O12.cif")
        structure.add_oxidation_state_by_element({"La": +3, "Fe": +3, "O": -2})
        for ion in magmom_setting:
            struct = structure.copy()
            struct.replace_species({"La3+": ion})
            vis = MPRelaxSet(struct)
            fe_pos = vis.poscar.comment.index("Fe")
            if fe_pos == 0:
                magmom_ans = [5] * 4 + [magmom_setting[ion]] * 4 + [0.6] * 12
            else:
                magmom_ans = [magmom_setting[ion]] * 4 + [5] * 4 + [0.6] * 12

            assert vis.incar["MAGMOM"] == magmom_ans


class TestMITMDSet(PymatgenTest):
    def setUp(self):
        self.set = MITMDSet
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(filepath)
        self.mit_md_param = self.set(self.struct, 300, 1200, 10000)

    def test_params(self):
        param = self.mit_md_param
        syms = param.potcar_symbols
        assert syms == ["Fe", "P", "O"]
        incar = param.incar
        assert "LDAUU" not in incar
        assert incar["EDIFF"] == approx(1e-5)
        kpoints = param.kpoints
        assert kpoints.kpts == [(1, 1, 1)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_as_from_dict(self):
        dct = self.mit_md_param.as_dict()
        input_set = dec.process_decoded(dct)
        assert isinstance(input_set, self.set)
        assert input_set._config_dict["INCAR"]["TEBEG"] == 300
        assert input_set._config_dict["INCAR"]["PREC"] == "Low"


@skip_if_no_psp_dir
class TestMVLNPTMDSet(PymatgenTest):
    def setUp(self):
        file_path = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mvl_npt_set = MVLNPTMDSet(self.struct, start_temp=0, end_temp=300, nsteps=1000)

    def test_incar(self):
        npt_set = self.mvl_npt_set

        syms = npt_set.potcar_symbols
        assert syms == ["Fe", "P", "O"]

        incar = npt_set.incar
        assert "LDAUU" not in incar
        assert incar["EDIFF"] == approx(1e-5)
        assert incar["LANGEVIN_GAMMA_L"] == 1
        assert incar["LANGEVIN_GAMMA"] == [10, 10, 10]
        enmax = max(npt_set.potcar[i].keywords["ENMAX"] for i in range(self.struct.ntypesp))
        assert incar["ENCUT"] == approx(1.5 * enmax)
        assert incar["ALGO"] == "Fast"
        assert incar["ISIF"] == 3
        assert incar["MDALGO"] == 3
        assert incar["SMASS"] == 0
        assert incar["PREC"] == "Low"

        kpoints = npt_set.kpoints
        assert kpoints.kpts == [(1, 1, 1)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_as_from_dict(self):
        dct = self.mvl_npt_set.as_dict()
        input_set = dec.process_decoded(dct)
        assert isinstance(input_set, MVLNPTMDSet)
        assert input_set._config_dict["INCAR"]["NSW"] == 1000


class TestMPMDSet(PymatgenTest):
    def setUp(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(filepath)
        self.struct_with_H = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR_hcp")
        self.mp_md_set_noTS = MPMDSet(self.struct, start_temp=0, end_temp=300, nsteps=1000)
        self.mp_md_set_noTS_with_H = MPMDSet(self.struct_with_H, start_temp=0, end_temp=300, nsteps=1000)
        self.mp_md_set_TS1 = MPMDSet(self.struct, start_temp=0, end_temp=300, nsteps=1000, time_step=1.0)

    def test_incar_no_ts(self):
        mpmdset = self.mp_md_set_noTS
        incar = mpmdset.incar

        assert incar["TEBEG"] == 0
        assert incar["TEEND"] == 300
        assert incar["NSW"] == 1000
        assert incar["IBRION"] == 0
        assert incar["ISIF"] == 0
        assert incar["ISYM"] == 0
        assert incar["POTIM"] == 2.0

    def test_incar_no_ts_with_h(self):
        mpmdset = self.mp_md_set_noTS_with_H
        incar = mpmdset.incar

        assert incar["POTIM"] == 0.5
        assert incar["NSW"] == 4000

    def test_incar_ts1(self):
        mpmdset = self.mp_md_set_TS1
        incar = mpmdset.incar

        assert incar["POTIM"] == 1.0
        assert incar["NSW"] == 1000

    def test_as_from_dict(self):
        d = self.mp_md_set_noTS.as_dict()
        v = dec.process_decoded(d)
        assert isinstance(v, MPMDSet)
        assert v._config_dict["INCAR"]["NSW"] == 1000


class TestMITNEBSet(PymatgenTest):
    def setUp(self):
        c1 = [[0.5] * 3, [0.9] * 3]
        c2 = [[0.5] * 3, [0.9, 0.1, 0.1]]
        s1 = Structure(Lattice.cubic(5), ["Si", "Si"], c1)
        s2 = Structure(Lattice.cubic(5), ["Si", "Si"], c2)
        self.structures = [Structure.from_sites(s.sites, to_unit_cell=True) for s in s1.interpolate(s2, 3, pbc=True)]
        self.vis = MITNEBSet(self.structures)

    def test_potcar_symbols(self):
        syms = self.vis.potcar_symbols
        assert syms == ["Si"]

    def test_incar(self):
        incar = self.vis.incar
        assert "LDAUU" not in incar
        assert incar["EDIFF"] == approx(0.00001)

    def test_kpoints(self):
        kpoints = self.vis.kpoints
        assert kpoints.kpts == [[25]]
        assert kpoints.style == Kpoints.supported_modes.Automatic

    def test_as_from_dict(self):
        d = self.vis.as_dict()
        v = dec.process_decoded(d)
        assert v._config_dict["INCAR"]["IMAGES"] == 2

    @skip_if_no_psp_dir
    def test_write_input(self):
        self.vis.write_input(".", write_cif=True, write_endpoint_inputs=True, write_path_cif=True)
        for file in "INCAR KPOINTS POTCAR 00/POSCAR 01/POSCAR 02/POSCAR 03/POSCAR 00/INCAR path.cif".split():
            assert os.path.isfile(file), f"{file=} not written"
        assert not os.path.isfile("04/POSCAR")


class TestMPSOCSet(PymatgenTest):
    def setUp(self):
        self.set = MPSOCSet

    def test_from_prev_calc(self):
        prev_run = f"{TEST_FILES_DIR}/fe_monomer"
        vis = self.set.from_prev_calc(
            prev_calc_dir=prev_run,
            magmom=[3],
            saxis=(1, 0, 0),
            user_incar_settings={"SIGMA": 0.025},
        )
        assert vis.incar["ISYM"] == -1
        assert vis.incar["LSORBIT"]
        assert vis.incar["ICHARG"] == 11
        assert vis.incar["SAXIS"] == [1, 0, 0]
        assert vis.incar["MAGMOM"] == [[0, 0, 3]]
        assert vis.incar["SIGMA"] == 0.025

    def test_override_from_prev_calc(self):
        # test override_from_prev_calc
        prev_run = f"{TEST_FILES_DIR}/fe_monomer"
        vis = self.set(
            dummy_structure,
            magmom=[3],
            saxis=(1, 0, 0),
            user_incar_settings={"SIGMA": 0.025},
        )
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["ISYM"] == -1
        assert vis.incar["LSORBIT"]
        assert vis.incar["ICHARG"] == 11
        assert vis.incar["SAXIS"] == [1, 0, 0]
        assert vis.incar["MAGMOM"] == [[0, 0, 3]]
        assert vis.incar["SIGMA"] == 0.025


class TestMPNMRSet(PymatgenTest):
    def test_incar(self):
        filepath = f"{TEST_FILES_DIR}/Li.cif"
        structure = Structure.from_file(filepath)

        vis = MPNMRSet(structure)
        assert vis.incar.get("LCHIMAG")
        assert vis.incar.get("QUAD_EFG") is None

        vis = MPNMRSet(structure, mode="efg")
        assert not vis.incar.get("LCHIMAG")
        assert vis.incar.get("QUAD_EFG") == [-0.808]
        for q in vis.incar["QUAD_EFG"]:
            assert isinstance(q, float)
            assert not isinstance(q, FloatWithUnit)

        vis = MPNMRSet(structure, mode="efg", isotopes=["Li-7"])
        assert not vis.incar.get("LCHIMAG")
        assert vis.incar.get("QUAD_EFG") == [-40.1]


@skip_if_no_psp_dir
class TestMVLSlabSet(PymatgenTest):
    def setUp(self):
        self.set = MVLSlabSet
        struct = self.get_structure("Li2O")
        gen = SlabGenerator(struct, (1, 0, 0), 10, 10)
        self.slab = gen.get_slab()
        self.bulk = self.slab.oriented_unit_cell

        vis_bulk = self.set(self.bulk, bulk=True)
        vis = self.set(self.slab)
        vis_dipole = self.set(self.slab, auto_dipole=True)

        self.d_bulk = vis_bulk.get_vasp_input()
        self.d_slab = vis.get_vasp_input()
        self.d_dipole = vis_dipole.get_vasp_input()
        self.vis = vis

    def test_user_incar_settings(self):
        # Make sure user incar settings properly override AMIX.
        si = self.get_structure("Si")
        vis = self.set(si, user_incar_settings={"AMIX": 0.1})
        assert vis.incar["AMIX"] == 0.1

    def test_bulk(self):
        incar_bulk = self.d_bulk["INCAR"]
        poscar_bulk = self.d_bulk["POSCAR"]

        assert incar_bulk["ISIF"] == 3
        assert incar_bulk["EDIFF"] == 1e-4
        assert incar_bulk["EDIFFG"] == -0.02
        assert poscar_bulk.structure.formula == self.bulk.formula

    def test_slab(self):
        incar_slab = self.d_slab["INCAR"]
        poscar_slab = self.d_slab["POSCAR"]
        potcar_slab = self.d_slab["POTCAR"]

        assert incar_slab["AMIN"] == 0.01
        assert incar_slab["AMIX"] == 0.2
        assert incar_slab["BMIX"] == 0.001
        assert incar_slab["NELMIN"] == 8
        # No volume relaxation during slab calculations
        assert incar_slab["ISIF"] == 2
        assert potcar_slab.functional == "PBE"
        assert potcar_slab.symbols[1] == "O"
        assert potcar_slab.symbols[0] == "Li_sv"
        assert poscar_slab.structure.formula == self.slab.formula
        # Test auto-dipole
        dipole_incar = self.d_dipole["INCAR"]
        assert dipole_incar["LDIPOL"]
        assert_allclose(dipole_incar["DIPOL"], [0.2323, 0.2323, 0.2165], atol=1e-4)
        assert dipole_incar["IDIPOL"] == 3

    def test_kpoints(self):
        kpoints_slab = self.d_slab["KPOINTS"].kpts[0]
        kpoints_bulk = self.d_bulk["KPOINTS"].kpts[0]

        assert kpoints_bulk[0] == kpoints_slab[0]
        assert kpoints_bulk[1] == kpoints_slab[1]
        assert kpoints_bulk[0] == 15
        assert kpoints_bulk[1] == 15
        assert kpoints_bulk[2] == 15
        # The last kpoint in a slab should always be 1
        assert kpoints_slab[2] == 1

    def test_as_from_dict(self):
        vis_dict = self.vis.as_dict()
        vis = self.set.from_dict(vis_dict)
        assert {*vis.as_dict()} == {*vis_dict}


class TestMVLElasticSet(PymatgenTest):
    def test_incar(self):
        mvlparam = MVLElasticSet(self.get_structure("Graphite"))
        incar = mvlparam.incar
        assert incar["IBRION"] == 6
        assert incar["NFREE"] == 2
        assert incar["POTIM"] == 0.015
        assert "NPAR" not in incar


@skip_if_no_psp_dir
class TestMVLGWSet(PymatgenTest):
    def setUp(self):
        self.set = MVLGWSet
        self.struct = PymatgenTest.get_structure("Li2O")

    def test_static(self):
        mvlgwsc = self.set(self.struct)
        incar = mvlgwsc.incar
        assert incar["SIGMA"] == 0.01
        kpoints = mvlgwsc.kpoints
        assert kpoints.style == Kpoints.supported_modes.Gamma
        symbols = mvlgwsc.potcar.symbols
        assert symbols == ["Li_sv_GW", "O_GW"]

    def test_diag(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        mvlgwdiag = self.set.from_prev_calc(prev_run, copy_wavecar=True, mode="diag")
        mvlgwdiag.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert mvlgwdiag.incar["NBANDS"] == 32
        assert mvlgwdiag.incar["ALGO"] == "Exact"
        assert mvlgwdiag.incar["LOPTICS"]

        # test override_from_prev_calc
        mvlgwdiag = self.set(dummy_structure, copy_wavecar=True, mode="diag")
        mvlgwdiag.override_from_prev_calc(prev_calc_dir=prev_run)
        mvlgwdiag.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert mvlgwdiag.incar["NBANDS"] == 32
        assert mvlgwdiag.incar["ALGO"] == "Exact"
        assert mvlgwdiag.incar["LOPTICS"]

    def test_bse(self):
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        mvlgwgbse = self.set.from_prev_calc(prev_run, copy_wavecar=True, mode="BSE")
        mvlgwgbse.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")

        prev_run = f"{TEST_FILES_DIR}/relaxation"
        mvlgwgbse = self.set.from_prev_calc(prev_run, copy_wavecar=False, mode="GW")
        assert mvlgwgbse.incar["NOMEGA"] == 80
        assert mvlgwgbse.incar["ENCUTGW"] == 250
        assert mvlgwgbse.incar["ALGO"] == "GW0"
        mvlgwgbse1 = self.set.from_prev_calc(prev_run, copy_wavecar=False, mode="BSE")
        assert mvlgwgbse1.incar["ANTIRES"] == 0
        assert mvlgwgbse1.incar["NBANDSO"] == 20
        assert mvlgwgbse1.incar["ALGO"] == "BSE"

        # test override_from_prev_calc
        prev_run = f"{TEST_FILES_DIR}/relaxation"
        mvlgwgbse = self.set(dummy_structure, copy_wavecar=True, mode="BSE")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        mvlgwgbse.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")

        prev_run = f"{TEST_FILES_DIR}/relaxation"
        mvlgwgbse = self.set(dummy_structure, copy_wavecar=True, mode="GW")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        assert mvlgwgbse.incar["NOMEGA"] == 80
        assert mvlgwgbse.incar["ENCUTGW"] == 250
        assert mvlgwgbse.incar["ALGO"] == "GW0"

        mvlgwgbse1 = self.set(dummy_structure, copy_wavecar=False, mode="BSE")
        mvlgwgbse1.override_from_prev_calc(prev_calc_dir=prev_run)
        assert mvlgwgbse1.incar["ANTIRES"] == 0
        assert mvlgwgbse1.incar["NBANDSO"] == 20
        assert mvlgwgbse1.incar["ALGO"] == "BSE"


class TestMPHSEBS(PymatgenTest):
    def setUp(self):
        self.set = MPHSEBSSet

    def test_init(self):
        prev_run = f"{TEST_FILES_DIR}/static_silicon"
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="uniform")
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 16

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="gap")
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 18

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="line")
        assert vis.incar["LHFCALC"]
        assert vis.incar["HFSCREEN"] == 0.2
        assert vis.incar["NSW"] == 0
        assert vis.incar["ISYM"] == 3
        assert len(vis.kpoints.kpts) == 180

        with pytest.warns(BadInputSetWarning, match=r"Hybrid functionals"):
            vis = self.set(PymatgenTest.get_structure("Li2O"), user_incar_settings={"ALGO": "Fast"})
            vis.incar.items()

    def test_override_from_prev_calc(self):
        prev_run = f"{TEST_FILES_DIR}/static_silicon"
        vis = self.set(dummy_structure, mode="uniform")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 16

        vis = self.set(dummy_structure, mode="gap")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 18

        vis = self.set(dummy_structure, mode="line")
        vis = vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis.incar["LHFCALC"]
        assert vis.incar["HFSCREEN"] == 0.2
        assert vis.incar["NSW"] == 0
        assert vis.incar["ISYM"] == 3
        assert len(vis.kpoints.kpts) == 180


class TestMVLScanRelaxSet(PymatgenTest):
    def setUp(self):
        self.set = MVLScanRelaxSet
        file_path = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mvl_scan_set = self.set(self.struct, user_potcar_functional="PBE_52", user_incar_settings={"NSW": 500})

    def test_incar(self):
        incar = self.mvl_scan_set.incar
        assert "METAGGA" in incar
        assert "LASPH" in incar
        assert "ADDGRID" in incar
        assert incar["NSW"] == 500

        # Test SCAN+rVV10
        scan_rvv10_set = self.set(self.struct, vdw="rVV10")
        assert scan_rvv10_set.incar["BPARAM"] == 15.7

    # @skip_if_no_psp_dir
    # def test_potcar(self):
    #
    #     test_potcar_set_1 = self.set(self.struct, user_potcar_functional="PBE_54")
    #     assert test_potcar_set_1.potcar.functional == "PBE_54"
    #
    #     with pytest.raises(
    #         ValueError, match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54'\)"
    #     ):
    #         self.set(self.struct, user_potcar_functional="PBE")
    #
    #     # https://github.com/materialsproject/pymatgen/pull/3022
    #     # same test also in MITMPRelaxSetTest above (for redundancy,
    #     # should apply to all classes inheriting from DictSet)
    #     for user_potcar_settings in [{"Fe": "Fe_pv"}, {"W": "W_pv"}, None]:
    #         for species in [("W", "W"), ("Fe", "W"), ("Fe", "Fe")]:
    #             struct = Structure(lattice=Lattice.cubic(3), species=species, coords=[[0, 0, 0], [0.5, 0.5, 0.5]])
    #             relax_set = MPRelaxSet(
    #                 structure=struct, user_potcar_functional="PBE_54", user_potcar_settings=user_potcar_settings
    #             )
    #             expected = {
    #                 **({"W": "W_sv"} if "W" in struct.symbol_set else {}),
    #                 **(user_potcar_settings or {}),
    #             } or None
    #             assert relax_set.user_potcar_settings == expected

    def test_as_from_dict(self):
        d = self.mvl_scan_set.as_dict()
        v = dec.process_decoded(d)
        assert isinstance(v, self.set)
        assert v._config_dict["INCAR"]["METAGGA"] == "SCAN"
        assert v.user_incar_settings["NSW"] == 500


class TestMPScanRelaxSet(PymatgenTest):
    def setUp(self):
        file_path = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mp_scan_set = MPScanRelaxSet(
            self.struct, user_potcar_functional="PBE_52", user_incar_settings={"NSW": 500}
        )

    def test_incar(self):
        incar = self.mp_scan_set.incar
        assert incar["METAGGA"] == "R2scan"
        assert incar["LASPH"]
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["NSW"] == 500
        # the default POTCAR contains metals
        assert incar["KSPACING"] == 0.22
        assert incar["ISMEAR"] == 2
        assert incar["SIGMA"] == 0.2

        # https://github.com/materialsproject/pymatgen/pull/3036
        for bandgap in (-1e-12, 1e-5, 1e-3):
            set_near_0_bandgap = MPScanRelaxSet(self.struct, bandgap=bandgap)
            expected = (
                {"KSPACING": 0.22, "SIGMA": 0.2, "ISMEAR": 2} if bandgap < 1e-4 else {"ISMEAR": -5, "SIGMA": 0.05}
            )
            for key, val in expected.items():
                assert set_near_0_bandgap.incar.get(key) == val

    def test_scan_substitute(self):
        mp_scan_sub = MPScanRelaxSet(
            self.struct,
            user_potcar_functional="PBE_52",
            user_incar_settings={"METAGGA": "SCAN"},
        )
        incar = mp_scan_sub.incar
        assert incar["METAGGA"] == "Scan"

    def test_bandgap_tol(self):
        # Test that the bandgap tolerance is applied correctly
        bandgap = 0.01
        for bandgap_tol, expected_kspacing in ((0.001, 0.2668137888), (0.02, 0.22)):
            incar = MPScanRelaxSet(self.struct, bandgap=0.01, bandgap_tol=bandgap_tol).incar
            assert incar["KSPACING"] == approx(expected_kspacing, abs=1e-5), f"{bandgap_tol=}, {bandgap=}"
            assert incar["ISMEAR"] == -5 if bandgap > bandgap_tol else 2
            assert incar["SIGMA"] == 0.05 if bandgap > bandgap_tol else 0.2

    def test_kspacing(self):
        # Test that KSPACING is capped at 0.44 for insulators
        file_path = f"{TEST_FILES_DIR}/POSCAR.O2"
        struct = Structure.from_file(file_path)
        for bandgap, expected in ((10, 0.44), (3, 0.4136617), (1.1, 0.3064757), (0.5, 0.2832948), (0, 0.22)):
            incar = MPScanRelaxSet(struct, bandgap=bandgap).incar
            assert incar["KSPACING"] == approx(expected, abs=1e-5)
            assert incar["ISMEAR"] == -5 if bandgap > 1e-4 else 2
            assert incar["SIGMA"] == 0.05 if bandgap > 1e-4 else 0.2

    def test_incar_overrides(self):
        # use 'user_incar_settings' to override the KSPACING, ISMEAR, and SIGMA
        # parameters that MPScanSet normally determines
        mp_scan_set2 = MPScanRelaxSet(
            self.struct,
            user_incar_settings={"KSPACING": 0.5, "ISMEAR": 0, "SIGMA": 0.05},
        )
        incar = mp_scan_set2.incar
        assert incar["KSPACING"] == 0.5
        assert incar["ISMEAR"] == 0
        assert incar["SIGMA"] == 0.05

    # Test SCAN+rVV10
    def test_rvv10(self):
        scan_rvv10_set = MPScanRelaxSet(self.struct, vdw="rVV10")
        assert "LUSE_VDW" in scan_rvv10_set.incar
        assert scan_rvv10_set.incar["BPARAM"] == 15.7

    def test_other_vdw(self):
        # should raise a warning.
        # IVDW key should not be present in the incar
        with pytest.warns(UserWarning, match=r"not supported at this time"):
            scan_vdw_set = MPScanRelaxSet(self.struct, vdw="DFTD3")
            assert "LUSE_VDW" not in scan_vdw_set.incar
            assert "IVDW" not in scan_vdw_set.incar

    @skip_if_no_psp_dir
    def test_potcar(self):
        assert self.mp_scan_set.potcar.functional == "PBE_52"

        # the default functional should be PBE_54
        input_set = MPScanRelaxSet(self.struct)
        assert input_set.potcar.functional == "PBE_54"

        with pytest.raises(
            ValueError, match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54'\)"
        ):
            MPScanRelaxSet(self.struct, user_potcar_functional="PBE")

    def test_as_from_dict(self):
        d = self.mp_scan_set.as_dict()
        input_set = dec.process_decoded(d)
        assert isinstance(input_set, MPScanRelaxSet)
        assert input_set._config_dict["INCAR"]["METAGGA"] == "R2SCAN"
        assert input_set.user_incar_settings["NSW"] == 500

    @skip_if_no_psp_dir
    def test_write_input(self):
        self.mp_scan_set.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/INCAR")
        assert not os.path.isfile(f"{self.tmp_path}/KPOINTS")
        assert os.path.isfile(f"{self.tmp_path}/POTCAR")
        assert os.path.isfile(f"{self.tmp_path}/POSCAR")


class TestMPScanStaticSet(PymatgenTest):
    def setUp(self):
        self.set = MPScanStaticSet
        self.prev_run = f"{TEST_FILES_DIR}/scan_relaxation"
        # test inheriting from a previous SCAN relaxation
        self.vis = self.set.from_prev_calc(prev_calc_dir=self.prev_run)

    def test_init(self):
        vis, prev_run = self.vis, self.prev_run
        # check that StaticSet settings were applied
        assert vis.incar["NSW"] == 0
        assert vis.incar["LREAL"] is False
        assert vis.incar["LORBIT"] == 11
        assert vis.incar["LVHAR"]
        assert vis.incar["ISMEAR"] == -5
        # Check that ENCUT and other INCAR settings were inherited.
        assert vis.incar["ENCUT"] == 680
        assert vis.incar["METAGGA"] == "R2scan"
        assert vis.incar["KSPACING"] == 0.34292842

        non_prev_vis = self.set(
            vis.structure,
            user_incar_settings={"ENCUT": 800, "LORBIT": 12, "LWAVE": True},
        )
        # check that StaticSet settings were applied
        assert non_prev_vis.incar["NSW"] == 0
        assert non_prev_vis.incar["LREAL"] is False
        assert non_prev_vis.incar["LVHAR"]
        assert vis.incar["ISMEAR"] == -5
        # Check that ENCUT and other INCAR settings were inherited.
        assert non_prev_vis.incar["METAGGA"] == "R2scan"
        # the KSPACING will have the default value here, since no previous calc
        assert non_prev_vis.incar["KSPACING"] == 0.22
        # Check that user incar settings are applied.
        assert non_prev_vis.incar["ENCUT"] == 800
        assert non_prev_vis.incar["LORBIT"] == 12
        assert non_prev_vis.incar["LWAVE"]

        v2 = self.set.from_dict(non_prev_vis.as_dict())
        # Check that user incar settings are applied.
        assert v2.incar["ENCUT"] == 800
        assert v2.incar["LORBIT"] == 12
        assert non_prev_vis.incar["LWAVE"]

        # Check LCALCPOL flag
        lcalcpol_vis = self.set.from_prev_calc(prev_calc_dir=prev_run, lcalcpol=True)
        assert lcalcpol_vis.incar["LCALCPOL"]

        # Check LEPSILON flag
        lepsilon_vis = self.set.from_prev_calc(prev_calc_dir=prev_run, lepsilon=True)
        assert lepsilon_vis.incar["LEPSILON"]
        assert lepsilon_vis.incar["LPEAD"]
        assert lepsilon_vis.incar["IBRION"] == 8
        assert lepsilon_vis.incar.get("NSW") is None
        assert lepsilon_vis.incar.get("NPAR") is None

    def test_as_from_dict(self):
        vis_dict = self.vis.as_dict()
        assert vis_dict == self.set.from_dict(vis_dict).as_dict()

    def test_override_from_prev_calc(self):
        vis, prev_run = self.vis, self.prev_run
        vis.override_from_prev_calc(prev_calc_dir=prev_run)
        # check that StaticSet settings were applied
        assert vis.incar["NSW"] == 0
        assert vis.incar["LREAL"] is False
        assert vis.incar["LORBIT"] == 11
        assert vis.incar["LVHAR"]
        assert vis.incar["ISMEAR"] == -5
        # Check that ENCUT and other INCAR settings were inherited.
        assert vis.incar["ENCUT"] == 680
        assert vis.incar["METAGGA"] == "R2scan"
        assert vis.incar["KSPACING"] == 0.34292842

        # Check LCALCPOL flag
        lcalcpol_vis = self.set(dummy_structure, lcalcpol=True)
        lcalcpol_vis = lcalcpol_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert lcalcpol_vis.incar["LCALCPOL"]

        # Check LEPSILON flag
        lepsilon_vis = self.set(dummy_structure, lepsilon=True)
        lepsilon_vis = lepsilon_vis.override_from_prev_calc(prev_calc_dir=prev_run)
        assert lepsilon_vis.incar["LEPSILON"]
        assert lepsilon_vis.incar["LPEAD"]
        assert lepsilon_vis.incar["IBRION"] == 8
        assert lepsilon_vis.incar.get("NSW") is None
        assert lepsilon_vis.incar.get("NPAR") is None


class TestFunc(PymatgenTest):
    @skip_if_no_psp_dir
    def test_batch_write_input(self):
        structs = list(map(PymatgenTest.get_structure, ["Li2O", "LiFePO4"]))

        batch_write_input(structs)
        for d in ["Li4Fe4P4O16_1", "Li2O1_0"]:
            for f in ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]:
                assert os.path.isfile(os.path.join(d, f))


@skip_if_no_psp_dir
class TestMVLGBSet(PymatgenTest):
    def setUp(self):
        filepath = f"{TEST_FILES_DIR}/Li.cif"
        self.struct = Structure.from_file(filepath)

        self.bulk = MVLGBSet(self.struct)
        self.slab = MVLGBSet(self.struct, slab_mode=True)

        self.d_bulk = self.bulk.get_vasp_input()
        self.d_slab = self.slab.get_vasp_input()

    def test_bulk(self):
        incar_bulk = self.d_bulk["INCAR"]
        assert incar_bulk["ISIF"] == 3

    def test_slab(self):
        incar_slab = self.d_slab["INCAR"]
        assert incar_slab["ISIF"] == 2

    def test_kpoints(self):
        kpoints = self.d_slab["KPOINTS"]
        k_a = int(40 / (self.struct.lattice.abc[0]) + 0.5)
        k_b = int(40 / (self.struct.lattice.abc[1]) + 0.5)
        assert kpoints.kpts == [[k_a, k_b, 1]]


class TestMVLRelax52Set(PymatgenTest):
    def setUp(self):
        self.set = MVLRelax52Set
        file_path = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mvl_rlx_set = self.set(self.struct, user_potcar_functional="PBE_54", user_incar_settings={"NSW": 500})

    def test_incar(self):
        incar = self.mvl_rlx_set.incar
        assert "NSW" in incar
        assert incar["LREAL"] == "Auto"

    @skip_if_no_psp_dir
    def test_potcar(self):
        assert self.mvl_rlx_set.potcar.functional == "PBE_54"
        assert "Fe" in self.mvl_rlx_set.potcar.symbols

        self.struct.remove_species(["Fe"])
        test_potcar_set_1 = self.set(self.struct, user_potcar_functional="PBE_52")
        assert test_potcar_set_1.potcar.functional == "PBE_52"

        with pytest.raises(
            ValueError, match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54'\)"
        ):
            self.set(self.struct, user_potcar_functional="PBE")

    def test_as_from_dict(self):
        dct = self.mvl_rlx_set.as_dict()
        vasp_input = dec.process_decoded(dct)
        assert isinstance(vasp_input, self.set)
        assert vasp_input.incar["NSW"] == 500


class TestLobsterSet(PymatgenTest):
    def setUp(self):
        self.set = LobsterSet
        file_path = f"{TEST_FILES_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        # test for different parameters!
        self.lobsterset1 = self.set(self.struct, isym=-1, ismear=-5)
        self.lobsterset2 = self.set(self.struct, isym=0, ismear=0)
        # only allow isym=-1 and isym=0
        with pytest.raises(ValueError, match="Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0"):
            self.lobsterset_new = self.set(self.struct, isym=2, ismear=0)
        with pytest.raises(ValueError, match="Lobster usually works with ismear=-5 or ismear=0"):
            self.lobsterset_new = self.set(self.struct, isym=-1, ismear=2)
        # test if one can still hand over grid density of kpoints
        self.lobsterset3 = self.set(self.struct, isym=0, ismear=0, user_kpoints_settings={"grid_density": 6000})
        # check if users can overwrite settings in this class with the help of user_incar_settings
        self.lobsterset4 = self.set(self.struct, user_incar_settings={"ALGO": "Fast"})
        # use basis functions supplied by user
        self.lobsterset5 = self.set(
            self.struct,
            user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s", "O": "2p 2s"},
        )
        with pytest.raises(ValueError, match="There are no basis functions for the atom type O"):
            self.lobsterset6 = self.set(self.struct, user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s"}).incar
        self.lobsterset7 = self.set(
            self.struct,
            address_basis_file=f"{MODULE_DIR}/../lobster/lobster_basis/BASIS_PBE_54_standard.yaml",
        )
        with pytest.warns(BadInputSetWarning, match="Overriding the POTCAR"):
            self.lobsterset6 = self.set(self.struct)

        # test W_sw
        self.lobsterset8 = self.set(Structure.from_file(f"{TEST_FILES_DIR}/cohp/POSCAR.W"))

    def test_incar(self):
        incar1 = self.lobsterset1.incar
        assert "NBANDS" in incar1
        assert incar1["NBANDS"] == 116
        assert incar1["NSW"] == 0
        assert incar1["ISMEAR"] == -5
        assert incar1["ISYM"] == -1
        assert incar1["ALGO"] == "Normal"
        assert incar1["EDIFF"] == 1e-6
        incar2 = self.lobsterset2.incar
        assert incar2["ISYM"] == 0
        assert incar2["ISMEAR"] == 0
        incar4 = self.lobsterset4.incar
        assert incar4["ALGO"] == "Fast"

    def test_kpoints(self):
        kpoints1 = self.lobsterset1.kpoints
        assert kpoints1.comment.split()[6], 6138
        kpoints2 = self.lobsterset2.kpoints
        assert kpoints2.comment.split()[6], 6138
        kpoints3 = self.lobsterset3.kpoints
        assert kpoints3.comment.split()[6], 6000

    @skip_if_no_psp_dir
    def test_potcar(self):
        # PBE_54 is preferred at the moment
        assert self.lobsterset1.user_potcar_functional == "PBE_54"

    def test_as_from_dict(self):
        dict_here = self.lobsterset1.as_dict()

        lobsterset_new = self.set.from_dict(dict_here)
        # test relevant parts again
        incar1 = lobsterset_new.incar
        assert "NBANDS" in incar1
        assert incar1["NBANDS"] == 116
        assert incar1["NSW"] == 0
        assert incar1["NSW"] == 0
        assert incar1["ISMEAR"] == -5
        assert incar1["ISYM"] == -1
        assert incar1["ALGO"] == "Normal"
        kpoints1 = lobsterset_new.kpoints
        assert kpoints1.comment.split()[6], 6138
        assert lobsterset_new.user_potcar_functional == "PBE_54"


@skip_if_no_psp_dir
class TestMPAbsorptionSet(PymatgenTest):
    def setUp(self):
        file_path = f"{TEST_FILES_DIR}/absorption/static/POSCAR"
        self.structure = Structure.from_file(file_path)

    def test_ipa(self):
        prev_run = f"{TEST_FILES_DIR}/absorption/static"
        absorption_ipa = MPAbsorptionSet.from_prev_calc(prev_calc_dir=prev_run, copy_wavecar=True, mode="IPA")
        absorption_ipa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert absorption_ipa.incar["NBANDS"] == 32
        assert absorption_ipa.incar["ALGO"] == "Exact"
        assert absorption_ipa.incar["LOPTICS"]

        # test override_from_prev_calc
        absorption_ipa = MPAbsorptionSet(dummy_structure, copy_wavecar=True, mode="IPA")
        absorption_ipa.override_from_prev_calc(prev_calc_dir=prev_run)
        absorption_ipa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert absorption_ipa.incar["NBANDS"] == 32
        assert absorption_ipa.incar["ALGO"] == "Exact"
        assert absorption_ipa.incar["LOPTICS"]

    def test_rpa(self):
        prev_run = f"{TEST_FILES_DIR}/absorption/ipa"
        absorption_rpa = MPAbsorptionSet.from_prev_calc(prev_run, copy_wavecar=True, mode="RPA")
        absorption_rpa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")
        assert absorption_rpa.incar["NOMEGA"] == 1000
        assert absorption_rpa.incar["NBANDS"] == 48
        assert absorption_rpa.incar["ALGO"] == "CHI"

        # test override_from_prev_calc
        prev_run = f"{TEST_FILES_DIR}/absorption/ipa"
        absorption_rpa = MPAbsorptionSet(dummy_structure, copy_wavecar=True, mode="RPA")
        absorption_rpa.override_from_prev_calc(prev_calc_dir=prev_run)
        absorption_rpa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")
        assert absorption_rpa.incar["NOMEGA"] == 1000
        assert absorption_rpa.incar["NBANDS"] == 48
        assert absorption_rpa.incar["ALGO"] == "CHI"
