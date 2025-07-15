from __future__ import annotations

import hashlib
import os
from glob import glob
from zipfile import ZipFile

import numpy as np
import pytest
from monty.json import MontyDecoder
from monty.serialization import loadfn
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import SETTINGS, Lattice, Species, Structure
from pymatgen.core.composition import Composition
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.units import FloatWithUnit
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, PotcarSingle
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.sets import (
    MODULE_DIR,
    BadInputSetWarning,
    CINEBSet,
    DictSet,
    LobsterSet,
    MatPESStaticSet,
    MITMDSet,
    MITNEBSet,
    MITRelaxSet,
    MP24RelaxSet,
    MP24StaticSet,
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
    NEBSet,
    VaspInputGenerator,
    VaspInputSet,
    batch_write_input,
    get_structure_from_prev_run,
    get_valid_magmom_struct,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/io/vasp"

pytest.MonkeyPatch().setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(FAKE_POTCAR_DIR))

NO_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR") is None
skip_if_no_psp_dir = pytest.mark.skipif(NO_PSP_DIR, reason="PMG_VASP_PSP_DIR is not set")

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


class TestSetChangeCheck(MatSciTest):
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
            with open(input_set, encoding="utf-8") as file:
                text = file.read().encode("utf-8")
                name = os.path.basename(input_set)
                hashes[name] = hashlib.sha256(text).hexdigest()

        known_hashes = {
            "MVLGWSet.yaml": "eba4564a18b99494a08ab6fdbe5364e7212b5992c7a9ef109001ce314a5b33db",
            "MVLRelax52Set.yaml": "3660879566a9ee2ab289e81d7916335b2f33ab24dcb3c16ba7aaca9ff22dfbad",
            "MPHSERelaxSet.yaml": "1779cb6a6af43ad54a12aec22882b9b8aa3469b764e29ac4ab486960d067b811",
            "VASPIncarBase.yaml": "8c1ce90d6697e45b650e1881e2b3d82a733dba17fb1bd73747a38261ec65a4c4",
            "MPSCANRelaxSet.yaml": "ad652ea740d06f9edd979494f31e25074b82b9fffdaaf7eff2ae5541fb0e6288",
            "PBE64Base.yaml": "40e7e42159f59543b17f512666916001045f7644f422ccc45b8466d6a1cf0c48",
            "MPRelaxSet.yaml": "c9b0a519588fb3709509a9f9964632692584905e2961a0fe2e5f657561913083",
            "MITRelaxSet.yaml": "0b4bec619fa860dac648584853c3b3d5407e4148a85d0e95024fbd1dc315669d",
            "vdW_parameters.yaml": "7d2599a855533865335a313c043b6f89e03fc2633c88b6bc721723d94cc862bd",
            "MatPESStaticSet.yaml": "4ec60ad4bbbb9a756f1b3fea8ca4eab8fc767d8f6a67332e7af3908c910fd7c5",
            "MPAbsorptionSet.yaml": "e49cd0ab87864f1c244e9b5ceb4703243116ec1fbb8958a374ddff07f7a5625c",
            "PBE54Base.yaml": "cdffe123eca8b19354554b60a7f8de9b8776caac9e1da2bd2a0516b7bfac8634",
            "MP24RelaxSet.yaml": "35a5d4456f01d644cf41218725c5e0896c59e1658045ecd1544579cbb1ed7b85",
        }

        for input_set, hash_str in hashes.items():
            assert hash_str == known_hashes[input_set], f"{input_set=}\n{msg}"


class TestVaspInputSet(MatSciTest):
    @classmethod
    def setup_class(cls):
        filepath = f"{VASP_IN_DIR}/POSCAR"
        cls.structure = Structure.from_file(filepath)

    def test_as_dict(self):
        # https://github.com/materialsproject/pymatgen/pull/3031
        dict_set = VaspInputSet(self.structure, config_dict={"INCAR": {}}, user_potcar_functional="PBE_54")
        assert {*dict_set.as_dict()} >= {
            "@class",
            "@module",
            "@version",
            "auto_ismear",
            "bandgap_tol",
            "config_dict",
            "constrain_total_magmom",
            "files_to_transfer",
            "force_gamma",
            "inherit_incar",
            "international_monoclinic",
            "reduce_structure",
            "sort_structure",
            "standardize",
            "structure",
            "sym_prec",
            "use_structure_charge",
            "user_incar_settings",
            "user_kpoints_settings",
            "user_potcar_functional",
            "user_potcar_settings",
            "validate_magmom",
            "vdw",
        }
        assert dict_set.potcar_functional == dict_set.user_potcar_functional


class TestMITMPRelaxSet(MatSciTest):
    @classmethod
    def setup_class(cls):
        cls.set = MITRelaxSet
        cls.mp_set = MPRelaxSet

        filepath = f"{VASP_IN_DIR}/POSCAR"
        cls.structure = Structure.from_file(filepath)
        cls.coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        cls.lattice = Lattice(
            [
                [3.8401979337, 0, 0],
                [1.9200989668, 3.3257101909, 0],
                [0, -2.2171384943, 3.1355090603],
            ]
        )

        cls.mit_set = cls.set(cls.structure)
        cls.mit_set_unsorted = cls.set(cls.structure, sort_structure=False)
        cls.mp_set = MPRelaxSet(cls.structure)

    def test_pbe64(self):
        vis = MPRelaxSet(self.structure, user_potcar_functional="PBE_64")
        assert vis.potcar[0].keywords["TITEL"] == "PAW_PBE Fe_pv 02Aug2007"
        assert vis.potcar[1].keywords["TITEL"] == "PAW_PBE P 06Sep2000"
        assert vis.potcar[2].keywords["TITEL"] == "PAW_PBE O 08Apr2002"

    def test_no_structure_init(self):
        # basic test of initialization with no structure.
        vis = MPRelaxSet()
        assert vis.as_dict()["structure"] is None
        assert vis.inherit_incar is False
        props_to_test = ("incar", "kpoints", "poscar")
        for prop in props_to_test:
            with pytest.raises(RuntimeError, match="No structure is associated with the input set!"):
                _ = getattr(vis, prop)

        assert vis.structure is None
        inputs = vis.get_input_set(structure=self.structure)
        assert {*inputs} == {"INCAR", "KPOINTS", "POSCAR", "POTCAR"}
        assert_allclose(vis.incar["LDAUU"], [5.3, 0, 0])
        assert vis.as_dict()["structure"] is not None
        assert "structure" not in vis.as_dict(verbosity=1)

    def test_warnings(self):
        structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Cu"], [[0, 0, 0]])

        vis = self.set(structure)
        with pytest.warns(BadInputSetWarning) as warns_metal:
            _ = vis.incar
        assert len(warns_metal) == 1
        vasp_docs_link = "See VASP recommendations on ISMEAR for metals (https://www.vasp.at/wiki/index.php/ISMEAR)."
        assert (
            str(warns_metal[0].message)
            == f"Relaxation of likely metal with ISMEAR < 0 ({vis.incar['ISMEAR']}). {vasp_docs_link}"
        )

        # test different warning for ismear == 0 and sigma > 0.05
        vis = self.set(structure, user_incar_settings={"ISMEAR": 0, "SIGMA": 0.1})
        with pytest.warns(BadInputSetWarning) as warns_metal:
            _ = vis.incar
        assert len(warns_metal) == 1
        assert (
            str(warns_metal[0].message)
            == f"ISMEAR = 0 with a small SIGMA ({vis.incar['SIGMA']}) detected. {vasp_docs_link}"
        )

        with pytest.warns(
            BadInputSetWarning,
            match="Large KSPACING value detected with ISMEAR = -5. Ensure that VASP "
            "generates an adequate number of KPOINTS, lower KSPACING, or set ISMEAR = 0",
        ) as warns_kspacing:
            vis = self.set(structure, user_incar_settings={"KSPACING": 1, "ISMEAR": -5})
            _ = vis.incar
            for warn in warns_kspacing:
                print("scoots:", warn.message)
        assert len(warns_kspacing) == 2

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
        lattice = Lattice(
            [
                [3.8401979337, 0, 0],
                [1.9200989668, 3.3257101909, 0],
                [0, -2.2171384943, 3.1355090603],
            ]
        )
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
        with pytest.MonkeyPatch().context() as monkeypatch:
            monkeypatch.setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(f"{VASP_IN_DIR}/wrong_potcars"))
            with pytest.warns(BadInputSetWarning, match="not known by pymatgen"):
                _ = self.set(structure).potcar

    def test_potcar_special_defaults(self):
        # https://github.com/materialsproject/pymatgen/pull/3022
        for user_potcar_settings in [{"Fe": "Fe_pv"}, {"W": "W_pv"}, None]:
            for species in [("W", "W"), ("Fe", "W"), ("Fe", "Fe")]:
                struct = Structure(
                    lattice=Lattice.cubic(3),
                    species=species,
                    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
                )
                relax_set = MPRelaxSet(
                    structure=struct,
                    user_potcar_functional="PBE_54",
                    user_potcar_settings=user_potcar_settings,
                )
                expected = {  # noqa: SIM222
                    **({"W": "W_sv"} if "W" in struct.symbol_set else {}),
                    **(user_potcar_settings or {}),
                } or None
                assert relax_set.user_potcar_settings == expected

    @skip_if_no_psp_dir
    def test_lda_potcar(self):
        structure = Structure(self.lattice, ["P", "Fe"], self.coords)
        potcar = self.set(structure, user_potcar_functional="LDA").potcar
        assert potcar.functional == "LDA"

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
            ValueError,
            match="Disordered structure with partial occupancies cannot be converted into POSCAR",
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

        assert_allclose(incar["LDAUU"], [5.3, 0, 0])
        assert incar["EDIFF"] == approx(0.0012)

        incar = self.mit_set.incar
        assert_allclose(incar["LDAUU"], [4.0, 0, 0])
        assert incar["EDIFF"] == approx(1e-5)

        si = 14
        coords = []
        coords.extend((np.array([0, 0, 0]), np.array([0.75, 0.5, 0.75])))

        # Silicon structure for testing.
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, [si, si], coords)
        incar = MPRelaxSet(struct).incar
        assert "LDAU" not in incar

        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, ["Fe", "Mn"], coords)

        incar = MPRelaxSet(struct).incar
        assert "LDAU" not in incar

        # check fluorides
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPRelaxSet(struct).incar
        assert_allclose(incar["LDAUU"], [5.3, 0])
        assert_allclose(incar["MAGMOM"], [5, 0.6])

        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = self.set(struct).incar
        assert_allclose(incar["LDAUU"], [4.0, 0])

        # This seems counterintuitive at first, but even if the prior INCAR has a MAGMOM flag,
        # because the structure has no site properties, the default MAGMOM is assigned from the
        # config dictionary.
        struct = Structure(lattice, ["Fe", "F"], coords)
        incar = MPStaticSet(struct, prev_incar=f"{VASP_IN_DIR}/INCAR").incar
        assert_allclose(incar["MAGMOM"], [5, 0.6])

        # Make sure this works with species.
        struct = Structure(lattice, ["Fe2+", "O2-"], coords)
        incar = MPRelaxSet(struct).incar
        assert_allclose(incar["LDAUU"], [5.3, 0])

        struct = Structure(lattice, ["Fe", "Mn"], coords, site_properties={"magmom": (5.2, -4.5)})
        incar = MPRelaxSet(struct).incar
        assert_allclose(incar["MAGMOM"], [-4.5, 5.2])

        incar = self.set(struct, sort_structure=False).incar
        assert_allclose(incar["MAGMOM"], [5.2, -4.5])

        struct = Structure(lattice, [Species("Fe2+", spin=4.1), "Mn"], coords)
        incar = MPRelaxSet(struct).incar
        assert_allclose(incar["MAGMOM"], [5, 4.1])

        struct = Structure(lattice, ["Mn3+", "Mn4+"], coords)
        incar = self.set(struct).incar
        assert incar["MAGMOM"] == [4, 3]

        user_set = MPRelaxSet(struct, user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}})
        assert_allclose(user_set.incar["MAGMOM"], [100, 0.6])

        no_encut_set = MPRelaxSet(struct, user_incar_settings={"ENCUT": None})
        assert "ENCUT" not in no_encut_set.incar

        # sulfide vs sulfate test
        coords = [[0, 0, 0], [0.75, 0.5, 0.75], [0.25, 0.5, 0]]

        struct = Structure(lattice, ["Fe", "Fe", "S"], coords)
        incar = self.set(struct).incar
        assert_allclose(incar["LDAUU"], [1.9, 0])

        # Make sure MP sulfides are ok.
        assert "LDAUU" not in MPRelaxSet(struct).incar

        struct = Structure(lattice, ["Fe", "S", "O"], coords)
        incar = self.set(struct).incar
        assert_allclose(incar["LDAUU"], [4.0, 0, 0])

        # Make sure MP sulfates are ok.
        assert_allclose(MPRelaxSet(struct).incar["LDAUU"], [5.3, 0, 0])

        # test for default LDAUU value
        user_set_ldauu_fallback = MPRelaxSet(struct, user_incar_settings={"LDAUU": {"Fe": 5.0, "S": 0}})
        assert_allclose(user_set_ldauu_fallback.incar["LDAUU"], [5.0, 0, 0])

        # Expected to be oxide (O is the most electronegative atom)
        struct = Structure(lattice, ["Fe", "O", "S"], coords)
        incar = self.set(struct).incar
        assert_allclose(incar["LDAUU"], [4.0, 0, 0])

        # Expected to be chloride (Cl is the most electronegative atom)
        struct = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = self.set(struct, user_incar_settings={"LDAU": True}).incar
        assert "LDAUU" not in incar  # LDAU = False

        # User set a compound to be sulfide by specifying values of "LDAUL" etc.
        struct = Structure(lattice, ["Fe", "Cl", "S"], coords)
        incar = self.set(
            struct,
            user_incar_settings={
                "LDAU": True,
                "LDAUL": {"Fe": 3},
                "LDAUU": {"Fe": 1.8},
            },
        ).incar
        assert_allclose(incar["LDAUL"], [3.0, 0, 0])
        assert_allclose(incar["LDAUU"], [1.8, 0, 0])

        # test that van-der-Waals parameters are parsed correctly

        vdw_par = loadfn(f"{MODULE_DIR}/vdW_parameters.yaml")
        with pytest.raises(
            KeyError,
            match=f"Invalid or unsupported van-der-Waals functional. Supported functionals are {', '.join(vdw_par)}.",
        ):
            self.set(struct, vdw="optB86")
        incar = self.set(struct, vdw="optB86b").incar
        assert incar["GGA"] == "Mk"
        assert incar["LUSE_VDW"]
        assert incar["PARAM1"] == approx(0.1234)

        # Test that NELECT is updated when a charge is present
        si = 14
        coords = []
        coords.append(np.array([0, 0, 0]))
        coords.append(np.array([0.75, 0.5, 0.75]))

        # Silicon structure for testing.
        lattice = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(lattice, [si, si], coords, charge=1)
        mpr = MPRelaxSet(struct, use_structure_charge=True)
        assert mpr.incar["NELECT"] == 7, "NELECT not properly set for nonzero charge"

        # test that NELECT does not get set when use_structure_charge = False
        mpr = MPRelaxSet(struct, use_structure_charge=False)
        assert "NELECT" not in mpr.incar, "NELECT should not be set when use_structure_charge is False"

        struct = Structure(lattice, ["Co", "O"], coords)
        mpr = MPRelaxSet(struct)
        assert_allclose(mpr.incar["MAGMOM"], [0.6, 0.6])
        struct = Structure(lattice, ["Co4+", "O"], coords)
        mpr = MPRelaxSet(struct)
        assert_allclose(mpr.incar["MAGMOM"], [1, 0.6])

        # test passing user_incar_settings and user_kpoint_settings of None
        for set_cls in [MPRelaxSet, MPStaticSet, MPNonSCFSet]:
            mp_set = set_cls(struct, user_incar_settings=None, user_kpoints_settings=None)
            assert mp_set.kpoints is not None
            assert mp_set.incar is not None

    def test_get_kpoints(self):
        kpoints = MPRelaxSet(self.structure).kpoints
        assert kpoints.kpts == [(2, 4, 5)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints = MPRelaxSet(self.structure, user_kpoints_settings={"reciprocal_density": 1000}).kpoints
        assert kpoints.kpts == [(6, 10, 13)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

        kpoints_obj = Kpoints(kpts=[(3, 3, 3)])
        kpoints_return = MPRelaxSet(self.structure, user_kpoints_settings=kpoints_obj).kpoints
        assert kpoints_return.kpts == [(3, 3, 3)]

        kpoints = self.mit_set.kpoints
        assert kpoints.kpts == [(25,)]
        assert kpoints.style == Kpoints.supported_modes.Automatic

        recip_param_set = MPRelaxSet(self.structure, force_gamma=True)
        recip_param_set.kpoints_settings = {"reciprocal_density": 40}
        kpoints = recip_param_set.kpoints
        assert kpoints.kpts == [(2, 4, 5)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

    @skip_if_no_psp_dir
    @pytest.mark.filterwarnings("ignore:get_vasp_input is deprecated")
    def test_get_vasp_input(self):
        dct = self.mit_set.get_vasp_input()
        assert dct["INCAR"]["ISMEAR"] == -5
        struct = self.structure.copy()
        struct.make_supercell(4)
        relax_set = MPRelaxSet(struct)
        dct = relax_set.get_vasp_input()
        assert dct["INCAR"]["ISMEAR"] == 0

    def test_mp_metal_relax_set(self):
        mp_metal_set = MPMetalRelaxSet(self.get_structure("Sn"))
        incar = mp_metal_set.incar
        assert incar["ISMEAR"] == 1
        assert incar["SIGMA"] == approx(0.2)
        kpoints = mp_metal_set.kpoints
        assert_allclose(kpoints.kpts[0], (5, 5, 5))

    def test_as_from_dict(self):
        mit_set = self.set(self.structure)
        mp_set = MPRelaxSet(self.structure)
        mp_user_set = MPRelaxSet(
            self.structure,
            user_incar_settings={"MAGMOM": {"Fe": 10, "S": -5, "Mn3+": 100}},
        )

        dct = mit_set.as_dict()
        val = MontyDecoder().process_decoded(dct)
        assert val._config_dict["INCAR"]["LDAUU"]["O"]["Fe"] == 4

        dct = mp_set.as_dict()
        val = MontyDecoder().process_decoded(dct)
        assert val._config_dict["INCAR"]["LDAUU"]["O"]["Fe"] == approx(5.3)

        dct = mp_user_set.as_dict()
        val = MontyDecoder().process_decoded(dct)
        assert isinstance(val, VaspInputSet)
        assert val.user_incar_settings["MAGMOM"] == {"Fe": 10, "S": -5, "Mn3+": 100}

    def test_hubbard_off_and_ediff_override(self):
        input_set = MPRelaxSet(self.structure, user_incar_settings={"LDAU": False, "EDIFF": 1e-10})
        assert "LDAUU" not in input_set.incar
        assert input_set.incar["EDIFF"] == approx(1e-10)
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
        assert {*os.listdir(self.tmp_path)} == {
            "INCAR",
            "KPOINTS",
            "POSCAR",
            "POTCAR.spec",
        }

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
        #     vis.get_input_set()

        vis = MPRelaxSet(struct, user_potcar_settings={"Fe": "Fe"}, validate_magmom=True)
        assert_allclose(vis.get_input_set()["INCAR"]["MAGMOM"], [1.0] * len(struct))

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

    def test_write_input_and_from_directory(self):
        structure = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.0), ["Fe"], [[0.0, 0.0, 0.0]])

        vis = self.set(structure=structure)
        input_set = vis.get_input_set()

        vis.write_input(output_dir=".")
        assert all(os.path.isfile(file) for file in ("INCAR", "KPOINTS", "POSCAR", "POTCAR"))
        input_set_from_dir = self.set().from_directory(".")

        assert all(input_set_from_dir[k] == input_set[k] for k in ("INCAR", "KPOINTS", "POTCAR"))
        # for some reason the POSCARs are not identical, but their structures and as_dict()'s are
        assert input_set_from_dir["POSCAR"].structure == input_set["POSCAR"].structure
        assert input_set_from_dir["POSCAR"].as_dict() == input_set["POSCAR"].as_dict()

    def test_get_nedos(self):
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.pbesol.xml.gz")
        vis = self.set(structure=vrun.structures[-1])
        # no `prev_vasprun` --> default value of NEDOS
        assert vis._get_nedos(0.1) == 2000
        vis.prev_vasprun = vrun
        assert vis._get_nedos(0.1) == approx(741, abs=1)


class TestMPStaticSet(MatSciTest):
    def setup_method(self):
        self.set = MPStaticSet

    def test_init(self):
        prev_run = f"{TEST_DIR}/fixtures/relaxation"

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run)
        assert vis.inherit_incar is True
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
        assert leps_vis.incar["EDIFF"] == approx(1e-5)
        assert "NPAR" not in leps_vis.incar
        assert leps_vis.incar["NSW"] == 1
        assert non_prev_vis.kpoints.kpts == [(11, 10, 10)]
        non_prev_vis = self.set(vis.structure, reciprocal_density=200)
        assert non_prev_vis.kpoints.kpts == [(14, 12, 12)]
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
            dummy_struct = Structure(
                lattice=[[0, 2, 2], [2, 0, 2], [2, 2, 0]],
                species=["Fe", "O"],
                coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
            )
            vis = self.set(dummy_struct, user_incar_settings={"LDAU": True, "LASPH": False})
            vis.incar.items()

    def test_user_incar_kspacing(self):
        # Make sure user KSPACING settings properly overrides KPOINTS.
        si = self.get_structure("Si")
        vis = self.set(si, user_incar_settings={"KSPACING": 0.22})
        assert vis.incar["KSPACING"] == approx(0.22)
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
        prev_run = f"{TEST_DIR}/fixtures/relaxation"

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
        with ZipFile(f"{self.tmp_path}/MPStaticSet.zip", mode="r") as zip_file:
            contents = zip_file.namelist()
            assert set(contents).issuperset({"INCAR", "POSCAR", "POTCAR.spec", "KPOINTS"})
            spec = zip_file.open("POTCAR.spec", mode="r").read().decode()
            assert spec == "Si"

    def test_grid_size_from_struct(self):
        # TODO grab a bunch_of_calculations store as a list of tuples
        # (structure, ngx, ..., ngxf, ...) where all the grid size values are generated by vasp
        # check that the code produces the same grid sizes
        fname = f"{TEST_DIR}/fixtures/grid_data_files/vasp_inputs_for_grid_check.json"
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
        assert static_set.calculate_ng(custom_encut=2000) == (
            [56, 96, 96],
            [112, 192, 192],
        )

        assert static_set.calculate_ng() == ([30, 48, 50], [60, 96, 100])
        # test `custom_prec` kwarg for final structure in above test using "NORMAL".
        assert static_set.calculate_ng(custom_prec="NORMAL") == (
            [24, 36, 40],
            [48, 72, 80],
        )


class TestMatPESStaticSet(MatSciTest):
    def setup_method(self):
        self.struct = Structure.from_file(f"{VASP_IN_DIR}/POSCAR")
        self.prev_incar = Incar.from_file(f"{VASP_IN_DIR}/INCAR")

    def test_default(self):
        input_set = MatPESStaticSet(self.struct)
        incar = input_set.incar
        assert incar["ALGO"] == "Normal"
        assert incar["EDIFF"] == approx(1.0e-05)
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["GGA"] == "Pe"
        assert incar["ISMEAR"] == 0
        assert incar["ISPIN"] == 2
        assert incar["KSPACING"] == approx(0.22)
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
        assert incar["SIGMA"] == approx(0.05)
        assert incar["LMAXMIX"] == 6
        assert input_set.potcar_symbols == ["Fe_pv", "P", "O"]

        assert input_set.potcar_functional == "PBE_64"  # test POTCARs default to PBE_64
        assert input_set.kpoints is None
        if os.path.isdir(f"{FAKE_POTCAR_DIR}/POT_PAW_PBE_64"):
            # this part only runs if POTCAR files are available
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
        assert incar["EDIFF"] == approx(1.0e-05)
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["GGA"] == "Pe"
        assert incar["ISMEAR"] == 0
        assert incar["ISPIN"] == 2
        assert incar["KSPACING"] == approx(0.22)
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
        assert incar["SIGMA"] == approx(0.05)
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
        assert_allclose(incar_u["LDAUU"], [5.3, 0, 0])
        assert default_u.potcar_symbols == ["Fe_pv", "P", "O"]
        assert default_u.potcar_functional == "PBE_64"
        assert default_u.kpoints is None

    def test_functionals(self):
        xc_functional = "LDA"
        with pytest.raises(
            ValueError,
            match=f"Unrecognized {xc_functional=}. Supported exchange-correlation functionals are ",
        ):
            MatPESStaticSet(self.struct, xc_functional=xc_functional)

        with pytest.warns(UserWarning, match="inconsistent with the recommended PBE_64"):
            diff_potcar = MatPESStaticSet(self.struct, user_potcar_functional="PBE")
            assert str(diff_potcar.potcar[0]) == str(PotcarSingle.from_symbol_and_functional("Fe_pv", "PBE"))

    def test_from_prev_calc(self):
        vis = MatPESStaticSet.from_prev_calc(f"{TEST_DIR}/fixtures/relaxation")
        incar = vis.incar
        assert incar["GGA"] == "Pe"
        assert incar["ALGO"] == "Normal"
        assert vis.potcar_symbols == ["Li_sv"]
        assert vis.kpoints is None


class TestMPNonSCFSet(MatSciTest):
    def setup_method(self):
        self.set = MPNonSCFSet

    @skip_if_no_psp_dir
    def test_init(self):
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        # check mode belong to ["line", "uniform", "boltztrap"]
        valid_modes = ("line", "uniform", "boltztrap")
        mode = "none"
        with pytest.raises(
            ValueError,
            match=f"Invalid {mode=}. Supported modes for NonSCF runs are {', '.join(map(repr, valid_modes))}",
        ):
            vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode=mode)

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
        assert vis.incar["LORBIT"] == 11
        # Check that the ENCUT has been inherited.
        assert vis.incar["ENCUT"] == 600
        # Check that the user_incar_settings works
        assert vis.incar["SIGMA"] == approx(0.025)
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
        prev_run = f"{TEST_DIR}/fixtures/relaxation"

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
        assert vis.incar["SIGMA"] == approx(0.025)
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
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Uniform", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape == (1, 3)

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Line", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape != (1, 3)

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="Boltztrap", copy_chgcar=False)
        assert np.array(vis.kpoints.kpts).shape != (1, 3)

    def test_optics(self):
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
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

        # Check NEDOS and ISMEAR set correctly
        assert vis.incar["NEDOS"] == 2001
        assert vis.incar["ISMEAR"] == -5
        assert vis.incar["ISYM"] == 2

        # Check OPTICS=True INCAR settings, LREAL needs to be False when LOPTICS=True
        assert vis.incar["CSHIFT"] == approx(1e-5)
        assert vis.incar["LREAL"] is False
        assert vis.incar["LOPTICS"]

        # Check spin off for magmom <0.02 for each site
        assert vis.incar["ISPIN"] == 1
        assert vis.kpoints.style == Kpoints.supported_modes.Gamma

    def test_user_kpoint_override(self):
        # default kpoints style is reciprocal, try setting to gamma
        user_kpoints_override = Kpoints(style=Kpoints.supported_modes.Gamma, kpts=((1, 1, 1),))

        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        vis = self.set.from_prev_calc(
            prev_calc_dir=prev_run,
            copy_chgcar=False,
            optics=True,
            mode="Uniform",
            nedos=2001,
            user_kpoints_settings=user_kpoints_override,
        )
        assert vis.kpoints.style == Kpoints.supported_modes.Gamma


class TestMagmomLdau(MatSciTest):
    def test_structure_from_prev_run(self):
        vrun = Vasprun(f"{VASP_OUT_DIR}/vasprun.magmom_ldau.xml.gz")
        structure = vrun.final_structure
        poscar = Poscar(structure)
        struct_magmom_decorated = get_structure_from_prev_run(vrun)
        ldau_ans = {"LDAUU": [5.3, 0.0], "LDAUL": [2, 0], "LDAUJ": [0.0, 0.0]}
        magmom_ans = [5.0, 5.0, 5.0, 5.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6]
        ldau_dict = {}
        for key in ("LDAUU", "LDAUJ", "LDAUL"):
            if hasattr(struct_magmom_decorated[0], key.lower()):
                magmoms = {site.specie.symbol: getattr(site, key.lower()) for site in struct_magmom_decorated}
                ldau_dict[key] = [magmoms[sym] for sym in poscar.site_symbols]
        magmom = [site.magmom for site in struct_magmom_decorated]
        assert ldau_dict == ldau_ans
        assert magmom == magmom_ans

    def test_ln_magmom(self):
        yaml_path = f"{MODULE_DIR}/VASPIncarBase.yaml"
        magmom_setting = loadfn(yaml_path)["INCAR"]["MAGMOM"]
        structure = Structure.from_file(f"{TEST_FILES_DIR}/cif/La4Fe4O12.cif")
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


class TestMITMDSet(MatSciTest):
    def setup_method(self):
        self.set = MITMDSet
        filepath = f"{VASP_IN_DIR}/POSCAR"
        self.struct = Structure.from_file(filepath)
        self.mit_md_param = self.set(structure=self.struct, start_temp=300, end_temp=1200, nsteps=10000)

    def test_params(self):
        param = self.mit_md_param
        syms = param.potcar_symbols
        assert syms == ["Fe", "P", "O"]
        incar = param.incar
        assert "LDAUU" not in incar
        assert incar["EDIFF"] == approx(1e-5)
        assert incar["ALGO"] == "Fast"
        assert incar["ISMEAR"] == 0
        assert incar["IBRION"] == 0
        assert incar["ISYM"] == 0
        kpoints = param.kpoints
        assert kpoints.kpts == [(1, 1, 1)]
        assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_as_from_dict(self):
        dct = self.mit_md_param.as_dict()
        input_set = MontyDecoder().process_decoded(dct)
        assert isinstance(input_set, self.set)
        assert input_set.incar["TEBEG"] == 300
        assert input_set.incar["TEEND"] == 1200
        assert input_set.incar["PREC"] == "Low"
        assert input_set.incar["POTIM"] == 2
        assert input_set.incar["NSW"] == 10000

    def test_user_heat_speed(self):
        vis = self.set(self.struct, start_temp=0, end_temp=1000, time_step=0.5)
        vis.nsteps = (vis.end_temp - vis.start_temp) / vis.time_step
        assert vis.incar["TEBEG"] == 0
        assert vis.incar["TEEND"] == 1000
        assert vis.incar["POTIM"] == approx(0.5)
        assert vis.incar["NSW"] == 2000


@skip_if_no_psp_dir
class TestMVLNPTMDSet(MatSciTest):
    def setup_method(self):
        file_path = f"{VASP_IN_DIR}/POSCAR"
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
        enmax = max(npt_set.potcar[idx].keywords["ENMAX"] for idx in range(self.struct.n_elems))
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
        input_set = MontyDecoder().process_decoded(dct)
        assert isinstance(input_set, MVLNPTMDSet)
        assert input_set.incar["NSW"] == 1000


class TestMPMDSet(MatSciTest):
    def setup_method(self):
        filepath = f"{VASP_IN_DIR}/POSCAR"
        self.struct = Structure.from_file(filepath)
        self.struct_with_H = Structure.from_file(f"{VASP_IN_DIR}/POSCAR_hcp")
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
        assert incar["POTIM"] == approx(2.0)

    def test_incar_no_ts_with_h(self):
        mpmdset = self.mp_md_set_noTS_with_H
        incar = mpmdset.incar

        assert incar["POTIM"] == approx(0.5)
        assert incar["NSW"] == 4000

    def test_incar_ts1(self):
        mpmdset = self.mp_md_set_TS1
        incar = mpmdset.incar

        assert incar["POTIM"] == approx(1.0)
        assert incar["NSW"] == 1000

    def test_as_from_dict(self):
        dct = self.mp_md_set_noTS.as_dict()
        v = MontyDecoder().process_decoded(dct)
        assert isinstance(v, MPMDSet)
        assert v.incar["NSW"] == 1000


class TestNEBSet(MatSciTest):
    def setup_method(self):
        c1 = [[0.5] * 3, [0.9] * 3]
        c2 = [[0.5] * 3, [0.9, 0.1, 0.1]]
        s1 = Structure(Lattice.cubic(5), ["Si", "Si"], c1)
        s2 = Structure(Lattice.cubic(5), ["Si", "Si"], c2)
        self.structures = [Structure.from_sites(s.sites, to_unit_cell=True) for s in s1.interpolate(s2, 3, pbc=True)]
        self.vis = NEBSet(self.structures)
        self.vis_MIT = MITNEBSet(self.structures)
        self.vis_cineb = CINEBSet(self.structures)

    def test_potcar_symbols(self):
        syms = self.vis.potcar_symbols
        assert syms == ["Si"]

    def test_unset_encut(self):
        vis = NEBSet(self.structures, unset_encut=True)
        assert "ENCUT" not in vis.incar
        n_structs = 2
        with pytest.raises(
            ValueError,
            match=f"You need at least 3 structures for an NEB, got {n_structs}",
        ):
            _ = NEBSet(self.structures[:n_structs], unset_encut=True)

    def test_incar(self):
        incar = self.vis.incar
        assert "LDAUU" not in incar
        assert incar["EDIFF"] == approx(0.00005)
        assert self.vis_MIT.incar["EDIFF"] == approx(0.00005)
        assert incar["NSW"] == 200
        assert incar["IBRION"] == 3
        assert "LCLIMB" in self.vis_cineb.incar

    def test_kpoints(self):
        kpoints = self.vis.kpoints
        assert kpoints.kpts == [(5, 5, 5)]
        assert kpoints.style == Kpoints.supported_modes.Gamma
        assert self.vis_MIT.kpoints.kpts == [(25,)]
        assert self.vis_MIT.kpoints.style == Kpoints.supported_modes.Automatic

    def test_as_from_dict(self):
        dct = self.vis_MIT.as_dict()
        v = MontyDecoder().process_decoded(dct)
        assert v.incar["IMAGES"] == 2
        dct = self.vis.as_dict()
        v = MontyDecoder().process_decoded(dct)
        assert v.incar["IMAGES"] == 2

    @skip_if_no_psp_dir
    def test_write_input(self):
        self.vis.write_input(".", write_cif=True, write_endpoint_inputs=True, write_path_cif=True)
        for file in "INCAR KPOINTS POTCAR 00/POSCAR 01/POSCAR 02/POSCAR 03/POSCAR 00/INCAR path.cif".split():
            assert os.path.isfile(file), f"{file=} not written"
        # check structures match
        assert len(self.vis.structures[0]) + 3 == len(Structure.from_file("path.cif"))
        assert not os.path.isfile("04/POSCAR")


class TestMPSOCSet(MatSciTest):
    def setup_method(self):
        self.set = MPSOCSet

    def test_from_prev_calc(self):
        prev_run = f"{TEST_DIR}/fixtures/fe_monomer"
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
        assert vis.incar["SIGMA"] == approx(0.025)

    def test_override_from_prev_calc(self):
        # test override_from_prev_calc
        prev_run = f"{TEST_DIR}/fixtures/fe_monomer"
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
        assert vis.incar["SIGMA"] == approx(0.025)


class TestMPNMRSet(MatSciTest):
    def test_incar(self):
        filepath = f"{TEST_FILES_DIR}/cif/Li.cif"
        structure = Structure.from_file(filepath)

        vis = MPNMRSet(structure)
        assert vis.incar.get("LCHIMAG")
        assert vis.incar.get("QUAD_EFG") is None

        vis = MPNMRSet(structure, mode="efg")
        assert not vis.incar.get("LCHIMAG")
        assert_allclose(vis.incar.get("QUAD_EFG"), [-0.808])
        for q in vis.incar["QUAD_EFG"]:
            assert isinstance(q, float)
            assert not isinstance(q, FloatWithUnit)

        vis = MPNMRSet(structure, mode="efg", isotopes=["Li-7"])
        assert not vis.incar.get("LCHIMAG")
        assert_allclose(vis.incar.get("QUAD_EFG"), [-40.1])


@skip_if_no_psp_dir
class TestMVLSlabSet(MatSciTest):
    def setup_method(self):
        self.set = MVLSlabSet
        struct = self.get_structure("Li2O")
        gen = SlabGenerator(struct, (1, 0, 0), 10, 10)
        self.slab = gen.get_slab()
        self.bulk = self.slab.oriented_unit_cell

        vis_bulk = self.set(self.bulk, bulk=True)
        vis = self.set(self.slab)
        vis_dipole = self.set(self.slab, auto_dipole=True)

        self.d_bulk = vis_bulk.get_input_set()
        self.d_slab = vis.get_input_set()
        self.d_dipole = vis_dipole.get_input_set()
        self.vis = vis

    def test_user_incar_settings(self):
        # Make sure user incar settings properly override AMIX.
        si = self.get_structure("Si")
        vis = self.set(si, user_incar_settings={"AMIX": 0.1})
        assert vis.incar["AMIX"] == approx(0.1)

    def test_bulk(self):
        incar_bulk = self.d_bulk["INCAR"]
        poscar_bulk = self.d_bulk["POSCAR"]

        assert incar_bulk["ISIF"] == 3
        assert incar_bulk["EDIFF"] == approx(1e-5)
        assert incar_bulk["EDIFFG"] == approx(-0.05)
        assert poscar_bulk.structure.formula == self.bulk.formula

    def test_slab(self):
        incar_slab = self.d_slab["INCAR"]
        poscar_slab = self.d_slab["POSCAR"]
        potcar_slab = self.d_slab["POTCAR"]

        assert incar_slab["AMIN"] == approx(0.01)
        assert incar_slab["AMIX"] == approx(0.2)
        assert incar_slab["BMIX"] == approx(0.001)
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
        # DIPOL must be a plain Python list of floats
        assert isinstance(dipole_incar["DIPOL"], list)
        for comp in dipole_incar["DIPOL"]:
            assert isinstance(comp, float)

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
        assert "structure" not in self.vis.as_dict(verbosity=1)


class TestMVLElasticSet(MatSciTest):
    def test_incar(self):
        mvlparam = MVLElasticSet(self.get_structure("Graphite"))
        incar = mvlparam.incar
        assert incar["IBRION"] == 6
        assert incar["NFREE"] == 2
        assert incar["POTIM"] == approx(0.015)
        assert "NPAR" not in incar


@skip_if_no_psp_dir
class TestMVLGWSet(MatSciTest):
    def setup_method(self):
        self.set = MVLGWSet
        self.struct = MatSciTest.get_structure("Li2O")

    def test_static(self):
        assert self.set.mode == "STATIC"
        mode = "EVGW"
        with pytest.raises(
            ValueError,
            match=f"Invalid {mode=}, supported modes are {', '.join(map(repr, MVLGWSet.SUPPORTED_MODES))}",
        ):
            _ = self.set(mode=mode)
        mvlgwsc = self.set(self.struct)
        incar = mvlgwsc.incar
        assert incar["SIGMA"] == approx(0.01)
        kpoints = mvlgwsc.kpoints
        assert kpoints.style == Kpoints.supported_modes.Gamma
        symbols = mvlgwsc.potcar.symbols
        assert symbols == ["Li_sv_GW", "O_GW"]

    def test_diag(self):
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
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
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        mvlgwgbse = self.set.from_prev_calc(prev_run, copy_wavecar=True, mode="BSE")
        mvlgwgbse.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")

        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        mvlgwgbse = self.set.from_prev_calc(prev_run, copy_wavecar=False, mode="GW")
        assert mvlgwgbse.incar["NOMEGA"] == 80
        assert mvlgwgbse.incar["ENCUTGW"] == 250
        assert mvlgwgbse.incar["ALGO"] == "Gw0"
        mvlgwgbse1 = self.set.from_prev_calc(prev_run, copy_wavecar=False, mode="BSE")
        assert mvlgwgbse1.incar["ANTIRES"] == 0
        assert mvlgwgbse1.incar["NBANDSO"] == 20
        assert mvlgwgbse1.incar["ALGO"] == "Bse"

        # test override_from_prev_calc
        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        mvlgwgbse = self.set(dummy_structure, copy_wavecar=True, mode="BSE")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        mvlgwgbse.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")

        prev_run = f"{TEST_DIR}/fixtures/relaxation"
        mvlgwgbse = self.set(dummy_structure, copy_wavecar=True, mode="GW")
        mvlgwgbse.override_from_prev_calc(prev_calc_dir=prev_run)
        assert mvlgwgbse.incar["NOMEGA"] == 80
        assert mvlgwgbse.incar["ENCUTGW"] == 250
        assert mvlgwgbse.incar["ALGO"] == "Gw0"

        mvlgwgbse1 = self.set(dummy_structure, copy_wavecar=False, mode="BSE")
        mvlgwgbse1.override_from_prev_calc(prev_calc_dir=prev_run)
        assert mvlgwgbse1.incar["ANTIRES"] == 0
        assert mvlgwgbse1.incar["NBANDSO"] == 20
        assert mvlgwgbse1.incar["ALGO"] == "Bse"


class TestMPHSERelaxSet(MatSciTest):
    def setup_method(self):
        self.structure = dummy_structure
        self.set = MPHSERelaxSet

    def test_vdw_and_lasph_none(self):
        vis = self.set(self.structure, vdw=None)
        assert vis.incar["LASPH"], "LASPH is not set to True"
        vdw_keys = {"VDW_SR", "VDW_S8", "VDW_A1", "VDW_A2"}
        assert all(key not in vis.incar for key in vdw_keys), "Unexpected vdW parameters are set"

    def test_vdw_and_lasph_dftd3(self):
        vis = self.set(self.structure, vdw="dftd3")
        assert vis.incar["LASPH"], "LASPH is not set to True"
        assert vis.incar["VDW_SR"] == approx(1.129), "VDW_SR is not set correctly"
        assert vis.incar["VDW_S8"] == approx(0.109), "VDW_S8 is not set correctly"

    def test_vdw_and_lasph_dftd3_bj(self):
        vis = self.set(self.structure, vdw="dftd3-bj")
        assert vis.incar["LASPH"], "LASPH is not set to True"
        assert vis.incar["VDW_A1"] == approx(0.383), "VDW_A1 is not set correctly"
        assert vis.incar["VDW_S8"] == approx(2.310), "VDW_S8 is not set correctly"
        assert vis.incar["VDW_A2"] == approx(5.685), "VDW_A2 is not set correctly"

    def test_user_incar_settings(self):
        user_incar_settings = {"LASPH": False, "VDW_SR": 1.5}
        vis = self.set(self.structure, vdw="dftd3", user_incar_settings=user_incar_settings)
        assert not vis.incar["LASPH"], "LASPH user setting not applied"
        assert vis.incar["VDW_SR"] == approx(1.5), "VDW_SR user setting not applied"

    def test_from_prev_calc(self):
        prev_run = os.path.join(TEST_DIR, "fixtures", "relaxation")

        # Test for dftd3
        vis_d3 = self.set.from_prev_calc(prev_calc_dir=prev_run, vdw="dftd3")
        assert vis_d3.incar["LASPH"]
        assert "VDW_SR" in vis_d3.incar
        assert "VDW_S8" in vis_d3.incar

        # Test for dftd3-bj
        vis_bj = self.set.from_prev_calc(prev_calc_dir=prev_run, vdw="dftd3-bj")
        assert vis_bj.incar["LASPH"]
        assert "VDW_A1" in vis_bj.incar
        assert "VDW_A2" in vis_bj.incar
        assert "VDW_S8" in vis_bj.incar

    def test_override_from_prev_calc(self):
        prev_run = os.path.join(TEST_DIR, "fixtures", "relaxation")

        # Test for dftd3
        vis_d3 = self.set(self.structure, vdw="dftd3")
        vis_d3 = vis_d3.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis_d3.incar["LASPH"]
        assert "VDW_SR" in vis_d3.incar
        assert "VDW_S8" in vis_d3.incar

        # Test for dftd3-bj
        vis_bj = self.set(self.structure, vdw="dftd3-bj")
        vis_bj = vis_bj.override_from_prev_calc(prev_calc_dir=prev_run)
        assert vis_bj.incar["LASPH"]
        assert "VDW_A1" in vis_bj.incar
        assert "VDW_A2" in vis_bj.incar
        assert "VDW_S8" in vis_bj.incar


class TestMPHSEBS(MatSciTest):
    def setup_method(self):
        self.set = MPHSEBSSet

    def test_init(self):
        prev_run = f"{TEST_DIR}/fixtures/static_silicon"
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="uniform")
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 16

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="gap")
        assert vis.incar["LHFCALC"]
        assert len(vis.kpoints.kpts) == 18

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="line")
        assert vis.incar["LHFCALC"]
        assert vis.incar["HFSCREEN"] == approx(0.2)
        assert vis.incar["NSW"] == 0
        assert vis.incar["ISYM"] == 3
        assert len(vis.kpoints.kpts) == 180

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="uniform", reciprocal_density=50)
        assert vis.reciprocal_density == 50

        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="uniform", reciprocal_density=100)
        assert vis.reciprocal_density == 100

        uks = {"reciprocal_density": 100}
        vis = self.set.from_prev_calc(prev_calc_dir=prev_run, mode="uniform", user_kpoints_settings=uks)
        assert vis.reciprocal_density == 100

        with pytest.warns(BadInputSetWarning, match=r"Hybrid functionals"):
            vis = self.set(MatSciTest.get_structure("Li2O"), user_incar_settings={"ALGO": "Fast"})
            vis.incar.items()

    def test_override_from_prev_calc(self):
        prev_run = f"{TEST_DIR}/fixtures/static_silicon"
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
        assert vis.incar["HFSCREEN"] == approx(0.2)
        assert vis.incar["NSW"] == 0
        assert vis.incar["ISYM"] == 3
        assert len(vis.kpoints.kpts) == 180


class TestMVLScanRelaxSet(MatSciTest):
    def setup_method(self):
        self.set = MVLScanRelaxSet
        file_path = f"{VASP_IN_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mvl_scan_set = self.set(
            self.struct,
            user_potcar_functional="PBE_54",
            user_incar_settings={"NSW": 500},
        )

    def test_incar(self):
        incar = self.mvl_scan_set.incar
        assert incar["METAGGA"] == "Scan"
        assert incar["ADDGRID"] is True
        assert incar["LASPH"] is True
        assert incar["SIGMA"] == approx(0.05)
        assert incar["ISMEAR"] == -5
        assert incar["NSW"] == 500
        # Test SCAN+rVV10
        scan_rvv10_set = self.set(self.struct, vdw="rVV10")
        assert scan_rvv10_set.incar["BPARAM"] == approx(15.7)

    @skip_if_no_psp_dir
    def test_potcar(self):
        assert self.mvl_scan_set.potcar.functional == "PBE_54"

        # the default functional of MVLScanRelex is PBE_52
        input_set = MVLScanRelaxSet(self.struct)
        assert input_set.potcar.functional == "PBE_52"

        with pytest.raises(
            ValueError,
            match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54', 'PBE_64'\)",
        ):
            MVLScanRelaxSet(self.struct, user_potcar_functional="PBE")

    @pytest.mark.xfail(reason="TODO: need someone to fix this")
    @skip_if_no_psp_dir
    def test_potcar_need_fix(self):
        test_potcar_set_1 = self.set(self.struct, user_potcar_functional="PBE_54")
        assert test_potcar_set_1.potcar.functional == "PBE_54"

        with pytest.raises(
            ValueError,
            match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54', 'PBE_64'\)",
        ):
            self.set(self.struct, user_potcar_functional="PBE")

        # https://github.com/materialsproject/pymatgen/pull/3022
        # same test also in MITMPRelaxSetTest above (for redundancy,
        # should apply to all classes inheriting from VaspInputSet)
        for user_potcar_settings in [{"Fe": "Fe_pv"}, {"W": "W_pv"}, None]:
            for species in [("W", "W"), ("Fe", "W"), ("Fe", "Fe")]:
                struct = Structure(
                    lattice=Lattice.cubic(3),
                    species=species,
                    coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
                )
                relax_set = MPRelaxSet(
                    structure=struct,
                    user_potcar_functional="PBE_54",
                    user_potcar_settings=user_potcar_settings,
                )
                expected = {
                    **({"W": "W_sv"} if "W" in struct.symbol_set else {}),
                    **(user_potcar_settings or {}),
                }
                assert relax_set.user_potcar_settings == expected

    def test_as_from_dict(self):
        dct = self.mvl_scan_set.as_dict()
        v = MontyDecoder().process_decoded(dct)
        assert isinstance(v, self.set)
        assert v.incar["METAGGA"] == "Scan"
        assert v.user_incar_settings["NSW"] == 500


class TestMPScanRelaxSet(MatSciTest):
    def setup_method(self):
        file_path = f"{VASP_IN_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mp_scan_set = MPScanRelaxSet(
            self.struct,
            user_potcar_functional="PBE_52",
            user_incar_settings={"NSW": 500},
        )

    def test_incar(self):
        incar = self.mp_scan_set.incar
        assert incar["METAGGA"] == "R2scan"
        assert incar["LASPH"]
        assert incar["ENAUG"] == 1360
        assert incar["ENCUT"] == 680
        assert incar["NSW"] == 500
        # the default POTCAR contains metals, but no prev calc set --> bandgap unknown
        assert incar["KSPACING"] == approx(0.22)
        assert incar["ISMEAR"] == 0
        assert incar["SIGMA"] == approx(0.2)

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
            assert incar["SIGMA"] == approx(0.05 if bandgap > bandgap_tol else 0.2)

    def test_kspacing(self):
        # Test that KSPACING is capped at 0.44 for insulators
        file_path = f"{VASP_IN_DIR}/POSCAR_O2"
        struct = Structure.from_file(file_path)
        for bandgap, expected in (
            (10, 0.44),
            (3, 0.4136617),
            (1.1, 0.3064757),
            (0.5, 0.2832948),
            (0, 0.22),
        ):
            incar = MPScanRelaxSet(struct, bandgap=bandgap).incar
            assert incar["KSPACING"] == approx(expected, abs=1e-5)
            assert incar["ISMEAR"] == -5 if bandgap > 1e-4 else 2
            assert incar["SIGMA"] == approx(0.05 if bandgap > 1e-4 else 0.2)

    def test_incar_overrides(self):
        # use 'user_incar_settings' to override the KSPACING, ISMEAR, and SIGMA
        # parameters that MPScanSet normally determines
        mp_scan_set2 = MPScanRelaxSet(
            self.struct,
            user_incar_settings={"KSPACING": 0.5, "ISMEAR": 0, "SIGMA": 0.05},
        )
        incar = mp_scan_set2.incar
        assert incar["KSPACING"] == approx(0.5)
        assert incar["ISMEAR"] == 0
        assert incar["SIGMA"] == approx(0.05)

    # Test SCAN+rVV10
    def test_rvv10(self):
        scan_rvv10_set = MPScanRelaxSet(self.struct, vdw="rVV10")
        assert "LUSE_VDW" in scan_rvv10_set.incar
        assert scan_rvv10_set.incar["BPARAM"] == approx(15.7)

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
            ValueError,
            match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54', 'PBE_64'\)",
        ):
            MPScanRelaxSet(self.struct, user_potcar_functional="PBE")

    def test_as_from_dict(self):
        dct = self.mp_scan_set.as_dict()
        input_set = MontyDecoder().process_decoded(dct)
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


class TestMPScanStaticSet(MatSciTest):
    def setup_method(self):
        self.set = MPScanStaticSet
        self.prev_run = f"{TEST_DIR}/fixtures/scan_relaxation"
        # test inheriting from a previous SCAN relaxation
        self.vis = self.set.from_prev_calc(prev_calc_dir=self.prev_run)

    def test_init(self):
        vis, prev_run = self.vis, self.prev_run
        # check that StaticSet settings were applied
        assert vis.inherit_incar is True
        assert vis.incar["NSW"] == 0
        assert vis.incar["LREAL"] is False
        assert vis.incar["LORBIT"] == 11
        assert vis.incar["LVHAR"]
        assert vis.incar["ISMEAR"] == -5
        # Check that ENCUT and other INCAR settings were inherited.
        assert vis.incar["ENCUT"] == 680
        assert vis.incar["METAGGA"] == "R2scan"
        assert vis.incar["KSPACING"] == approx(0.34292842)

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
        assert non_prev_vis.incar["KSPACING"] == approx(0.22)
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
        assert lepsilon_vis.incar["NSW"] == 1
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
        assert vis.incar["KSPACING"] == approx(0.34292842)

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
        assert lepsilon_vis.incar["NSW"] == 1
        assert lepsilon_vis.incar.get("NPAR") is None


class TestFunc(MatSciTest):
    @skip_if_no_psp_dir
    def test_batch_write_input(self):
        structs = list(map(MatSciTest.get_structure, ("Li2O", "LiFePO4")))

        batch_write_input(structs, sanitize=True)
        for formula in ("Li4Fe4P4O16_1", "Li2O1_0"):
            for file in ("INCAR", "KPOINTS", "POSCAR", "POTCAR"):
                assert os.path.isfile(f"{formula}/{file}")


@skip_if_no_psp_dir
class TestMVLGBSet(MatSciTest):
    def setup_method(self):
        filepath = f"{TEST_FILES_DIR}/cif/Li.cif"
        self.struct = Structure.from_file(filepath)

        self.bulk = MVLGBSet(self.struct)
        self.slab = MVLGBSet(self.struct, slab_mode=True)

        self.d_bulk = self.bulk.get_input_set()
        self.d_slab = self.slab.get_input_set()

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
        assert kpoints.kpts == [(k_a, k_b, 1)]


class TestMVLRelax52Set(MatSciTest):
    def setup_method(self):
        self.set = MVLRelax52Set
        file_path = f"{VASP_IN_DIR}/POSCAR"
        self.struct = Structure.from_file(file_path)
        self.mvl_rlx_set = self.set(
            self.struct,
            user_potcar_functional="PBE_54",
            user_incar_settings={"NSW": 500},
        )

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
            ValueError,
            match=r"Invalid user_potcar_functional='PBE', must be one of \('PBE_52', 'PBE_54', 'PBE_64'\)",
        ):
            self.set(self.struct, user_potcar_functional="PBE")

    def test_as_from_dict(self):
        dct = self.mvl_rlx_set.as_dict()
        vasp_input = MontyDecoder().process_decoded(dct)
        assert isinstance(vasp_input, self.set)
        assert vasp_input.incar["NSW"] == 500


class TestLobsterSet(MatSciTest):
    def setup_method(self):
        self.set = LobsterSet
        file_path = f"{VASP_IN_DIR}/POSCAR"
        file_path2 = f"{VASP_IN_DIR}/POSCAR.lobster.spin_DOS"
        self.struct = Structure.from_file(file_path)
        self.struct2 = Structure.from_file(file_path2)

        # test for different parameters!
        self.lobsterset1 = self.set(self.struct, isym=-1, ismear=-5)
        self.lobsterset2 = self.set(self.struct, isym=0, ismear=0)
        # only allow isym=-1 and isym=0
        with pytest.raises(
            ValueError,
            match="Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0",
        ):
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
        self.lobsterset8 = self.set(Structure.from_file(f"{TEST_FILES_DIR}/electronic_structure/cohp/POSCAR.W"))

        # test if potcar selection is consistent with PBE_54
        self.lobsterset9 = self.set(self.struct2)

    def test_incar(self):
        incar1 = self.lobsterset1.incar
        assert "NBANDS" in incar1
        assert incar1["NBANDS"] == 116
        assert incar1["NSW"] == 0
        assert incar1["ISMEAR"] == -5
        assert incar1["ISYM"] == -1
        assert incar1["ALGO"] == "Normal"
        assert incar1["EDIFF"] == approx(1e-6)
        incar2 = self.lobsterset2.incar
        assert incar2["ISYM"] == 0
        assert incar2["ISMEAR"] == 0
        incar4 = self.lobsterset4.incar
        assert incar4["ALGO"] == "Fast"

    def test_kpoints(self):
        kpoints1 = self.lobsterset1.kpoints
        assert kpoints1.comment.split()[5] == "6138"
        kpoints2 = self.lobsterset2.kpoints
        assert kpoints2.comment.split()[5] == "6138"
        kpoints3 = self.lobsterset3.kpoints
        assert kpoints3.comment.split()[5] == "6000"

    @skip_if_no_psp_dir
    def test_potcar(self):
        # PBE_54 is preferred at the moment
        functional, symbol = "PBE_54", "K_sv"
        assert self.lobsterset1.user_potcar_functional == functional
        # test if potcars selected are consistent with PBE_54
        assert self.lobsterset2.potcar.symbols == ["Fe_pv", "P", "O"]
        # test if error raised contains correct potcar symbol for K element as PBE_54 set
        with pytest.raises(
            FileNotFoundError,
            match=f"You do not have the right POTCAR with {functional=} and {symbol=}",
        ):
            _ = self.lobsterset9.potcar.symbols

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
        assert kpoints1.comment.split()[5] == "6138"
        assert lobsterset_new.user_potcar_functional == "PBE_54"


@skip_if_no_psp_dir
class TestMPAbsorptionSet(MatSciTest):
    def setup_method(self):
        file_path = f"{TEST_DIR}/fixtures/absorption/static/POSCAR"
        self.structure = Structure.from_file(file_path)
        self.set = MPAbsorptionSet
        with pytest.raises(ValueError, match=r"STATIC not one of the support modes : \('IPA', 'RPA'\)"):
            self.set = MPAbsorptionSet(self.structure, mode="STATIC")

    def test_ipa(self):
        prev_run = f"{TEST_DIR}/fixtures/absorption/static"
        absorption_ipa = MPAbsorptionSet.from_prev_calc(
            prev_calc_dir=prev_run,
            user_incar_settings={"NEDOS": 3000},
            copy_wavecar=True,
            mode="IPA",
        )
        absorption_ipa.write_input(self.tmp_path)
        assert absorption_ipa.inherit_incar is True
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert absorption_ipa.incar["ENCUT"] == 680
        assert absorption_ipa.incar["NEDOS"] == 3000
        assert absorption_ipa.incar["NBANDS"] == 32
        assert absorption_ipa.incar["ALGO"] == "Exact"
        assert absorption_ipa.incar["LREAL"] is False
        assert absorption_ipa.incar["LWAVE"] is True
        assert absorption_ipa.incar["LOPTICS"]

        # test override_from_prev_calc
        absorption_ipa = MPAbsorptionSet(dummy_structure, copy_wavecar=True, mode="IPA")
        absorption_ipa.override_from_prev_calc(prev_calc_dir=prev_run)
        absorption_ipa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert absorption_ipa.incar["ENCUT"] == 680
        assert absorption_ipa.incar["NEDOS"] == 2001
        assert absorption_ipa.incar["NBANDS"] == 32
        assert absorption_ipa.incar["ALGO"] == "Exact"
        assert absorption_ipa.incar["LREAL"] is False
        assert absorption_ipa.incar["LWAVE"] is True
        assert absorption_ipa.incar["LOPTICS"]

    def test_rpa(self):
        prev_run = f"{TEST_DIR}/fixtures/absorption/ipa"
        absorption_rpa = MPAbsorptionSet.from_prev_calc(
            prev_run, user_incar_settings={"NEDOS": 3000}, copy_wavecar=True, mode="RPA"
        )
        absorption_rpa.write_input(self.tmp_path)
        assert absorption_rpa.inherit_incar
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")
        assert absorption_rpa.incar["ENCUT"] == 680
        assert absorption_rpa.incar["NEDOS"] == 3000
        assert absorption_rpa.incar["NOMEGA"] == 1000
        assert absorption_rpa.incar["NBANDS"] == 48
        assert absorption_rpa.incar["NKREDX"] == 13
        assert absorption_rpa.incar["ALGO"] == "Chi"
        assert absorption_rpa.incar["LREAL"] is False
        assert "LOPTICS" not in absorption_rpa.incar
        assert "LWAVE" not in absorption_rpa.incar

        # test override_from_prev_calc
        prev_run = f"{TEST_DIR}/fixtures/absorption/ipa"
        absorption_rpa = MPAbsorptionSet(dummy_structure, copy_wavecar=True, mode="RPA")
        absorption_rpa.override_from_prev_calc(prev_calc_dir=prev_run)
        absorption_rpa.write_input(self.tmp_path)
        assert os.path.isfile(f"{self.tmp_path}/WAVECAR")
        assert os.path.isfile(f"{self.tmp_path}/WAVEDER")
        assert absorption_rpa.incar["ENCUT"] == 680
        assert absorption_rpa.incar["NEDOS"] == 2001
        assert absorption_rpa.incar["NOMEGA"] == 1000
        assert absorption_rpa.incar["NBANDS"] == 48
        assert absorption_rpa.incar["NKREDX"] == 13
        assert absorption_rpa.incar["ALGO"] == "Chi"

        assert absorption_rpa.incar["LREAL"] is False
        assert "LOPTICS" not in absorption_rpa.incar
        assert "LWAVE" not in absorption_rpa.incar

    def test_kpoints(self):
        # Check IPA kpoints
        prev_run = f"{TEST_DIR}/fixtures/absorption/static"
        absorption_ipa = MPAbsorptionSet.from_prev_calc(prev_calc_dir=prev_run, mode="IPA")
        kpoints1 = absorption_ipa.kpoints
        assert kpoints1.kpts == [(13, 13, 13)]
        assert kpoints1.style == Kpoints.supported_modes.Gamma
        # Check RPA kpoints
        prev_run = f"{TEST_DIR}/fixtures/absorption/ipa"
        absorption_rpa = MPAbsorptionSet.from_prev_calc(prev_run, mode="RPA")
        kpoints2 = absorption_rpa.kpoints
        assert kpoints2.kpts == [(13, 13, 13)]
        assert kpoints2.style == Kpoints.supported_modes.Gamma

    def test_as_from_dict(self):
        # IPA_as_dict
        prev_run = f"{TEST_DIR}/fixtures/absorption/static"
        absorption_ipa = MPAbsorptionSet.from_prev_calc(prev_calc_dir=prev_run, mode="IPA")
        dct = absorption_ipa.as_dict()
        vasp_input = MontyDecoder().process_decoded(dct)
        assert vasp_input.incar["ALGO"] == "Exact"
        assert vasp_input.incar["LOPTICS"]
        assert vasp_input.incar["GGA"] == "Ps"
        # RPA_as_dict
        prev_run = f"{TEST_DIR}/fixtures/absorption/ipa"
        absorption_rpa = MPAbsorptionSet.from_prev_calc(prev_run, mode="RPA")
        dct = absorption_rpa.as_dict()
        vasp_input = MontyDecoder().process_decoded(dct)
        assert vasp_input.incar["ALGO"] == "Chi"
        assert vasp_input.incar["NBANDS"] == 48
        assert vasp_input.incar["GGA"] == "Ps"


def test_vasp_input_set_alias():
    assert VaspInputSet is VaspInputGenerator


def test_dict_set_alias():
    with pytest.warns(
        FutureWarning,
        match="DictSet is deprecated, and will be removed on 2025-12-31\nUse VaspInputSet",
    ):
        DictSet()
        assert isinstance(DictSet(), VaspInputSet)


class TestMP24Sets:
    def setup_method(self):
        self.relax_set = MP24RelaxSet
        self.static_set = MP24StaticSet

        filepath = f"{VASP_IN_DIR}/POSCAR"
        self.structure = Structure.from_file(filepath)

    @staticmethod
    def matches_ref(test_val, ref_val) -> bool:
        if isinstance(ref_val, float):
            return test_val == approx(ref_val)
        return test_val == ref_val

    def test_kspacing(self):
        bandgaps = [0.0, 0.1, 0.5, 1.0, 2.0, 3, 5, 10, 1e4]
        expected_kspacing = [
            0.22,
            0.22000976174867864,
            0.2200466614246148,
            0.22056799311325073,
            0.2876525546497567,
            0.40800309817134106,
            0.44000114629141485,
            0.4999999999540808,
            0.5,
        ]
        for i, bandgap in enumerate(bandgaps):
            assert self.relax_set()._multi_sigmoid_interp(bandgap) == approx(expected_kspacing[i])

    def test_default(self):
        vis = self.relax_set(structure=self.structure)
        expected_incar_relax = {
            "ALGO": "Normal",
            "EDIFF": 1.0e-05,
            "EDIFFG": -0.02,
            "ENAUG": 1360,
            "ENCUT": 680,
            "GGA_COMPAT": False,
            "KSPACING": 0.22,
            "ISMEAR": 0,
            "SIGMA": 0.05,
            "METAGGA": "R2scan",
            "LMAXMIX": 6,
            "LREAL": False,
        }

        assert all(self.matches_ref(vis.incar[k], v) for k, v in expected_incar_relax.items())

        assert self.relax_set(self.structure, xc_functional="r2SCAN")._config_dict == vis._config_dict
        assert vis.inherit_incar is False
        assert vis.dispersion is None
        assert vis.potcar_functional == "PBE_64"
        assert vis.potcar_symbols == ["Fe_pv", "P", "O"]
        assert vis.kpoints is None

        vis = self.static_set(structure=self.structure, dispersion="rVV10")

        expected_incar_static = {
            "ISMEAR": -5,
            "NSW": 0,
            "LORBIT": 11,
            "METAGGA": "R2scan",
            "LUSE_VDW": True,
            "BPARAM": 11.95,
            "CPARAM": 0.0093,
        }

        assert all(self.matches_ref(vis.incar[k], v) for k, v in expected_incar_static.items())
        assert (
            self.static_set(self.structure, xc_functional="r2SCAN", dispersion="rVV10")._config_dict == vis._config_dict
        )

    def test_non_default_xc_func(self):
        for xc_functional, vasp_name in {"PBE": "Pe", "PBEsol": "Ps"}.items():
            vis = self.relax_set(structure=self.structure, xc_functional=xc_functional)
            assert vis.incar.get("METAGGA") is None
            assert vis.incar["GGA"] == vasp_name

            vis = self.static_set(structure=self.structure, xc_functional=xc_functional, dispersion="D4")
            assert vis.incar.get("METAGGA") is None
            print(self.relax_set(structure=self.structure, xc_functional=xc_functional).incar)
            assert vis.incar["GGA"] == vasp_name
            assert all(
                isinstance(vis.incar[k], float | int)
                for k in (
                    "IVDW",
                    "VDW_A1",
                    "VDW_A2",
                    "VDW_S6",
                    "VDW_S8",
                )
            )
