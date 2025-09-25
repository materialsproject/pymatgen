from __future__ import annotations

import orjson
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.electronic_structure.cohp import (
    Cohp,
    CompleteCohp,
    IcohpCollection,
    IcohpValue,
    get_integrated_cohp_in_energy_range,
)
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp"


class TestCohp:
    def setup_method(self):
        with open(f"{TEST_DIR}/cohp.json", "rb") as file:
            self.cohp = Cohp.from_dict(orjson.loads(file.read()))
        self.cohp_only = Cohp(self.cohp.efermi, self.cohp.energies, self.cohp.cohp)
        with open(f"{TEST_DIR}/coop.json", "rb") as file:
            self.coop = Cohp.from_dict(orjson.loads(file.read()))
        with open(f"{TEST_DIR}/cobi.json", "rb") as file:
            self.cobi = Cohp.from_dict(orjson.loads(file.read()))

    def test_as_from_dict(self):
        with open(f"{TEST_DIR}/cohp.json", "rb") as file:
            cohp_dict = orjson.loads(file.read())
        assert self.cohp.as_dict() == cohp_dict

        with open(f"{TEST_DIR}/cobi.json", "rb") as file:
            cobi_dict = orjson.loads(file.read())
        assert self.cobi.as_dict() == cobi_dict

    def test_attributes(self):
        assert len(self.cohp.energies) == 301
        assert self.cohp.efermi == approx(9.75576)
        assert self.coop.efermi == approx(5.90043)
        assert not self.cohp.are_coops
        assert self.coop.are_coops
        assert not self.coop.are_cobis
        assert not self.cobi.are_coops
        assert self.cobi.are_cobis

    def test_get_icohp(self):
        assert self.cohp.get_icohp() == self.cohp.get_cohp(integrated=True)
        assert self.cohp_only.get_icohp() is None

    def test_get_interpolated_value(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {Spin.up: -0.10218, Spin.down: -0.19701}
        icoop_ef_dict = {Spin.up: 0.24714}
        icohp_ef = self.cohp.get_interpolated_value(self.cohp.efermi, integrated=True)
        icoop_ef = self.coop.get_interpolated_value(self.coop.efermi, integrated=True)
        assert icohp_ef_dict == approx(icohp_ef)
        assert icoop_ef_dict == approx(icoop_ef)
        with pytest.raises(ValueError, match="ICOHP is empty"):
            self.cohp_only.get_interpolated_value(5.0, integrated=True)

    def test_str(self):
        header = "#Energy          COOPUp          ICOOPUp        \n"

        with open(f"{TEST_DIR}/cohp.str", encoding="utf-8") as file:
            str_cohp = file.read()
        assert str(self.cohp) == str_cohp
        assert str(self.coop).strip().startswith(header)

        with open(f"{TEST_DIR}/coop.str", encoding="utf-8") as file:
            str_coop = file.read()
        assert str(self.coop) == str_coop
        assert str(self.coop).strip().startswith(header)

    def test_antibnd_states_below_efermi(self):
        assert self.cohp.has_antibnd_states_below_efermi(spin=None) == {
            Spin.up: True,
            Spin.down: True,
        }
        assert self.cohp.has_antibnd_states_below_efermi(spin=None, limit=0.5) == {
            Spin.up: False,
            Spin.down: False,
        }
        assert self.cohp.has_antibnd_states_below_efermi(spin=Spin.up, limit=0.5) == {Spin.up: False}


class TestIcohpValue:
    def setup_method(self):
        # without spin polarization
        label = "1"
        atom1 = "K1"
        atom2 = "F2"
        length = "2.3"
        translation = [-1, 0, 0]
        num = 1
        icohp = {Spin.up: -2.0}
        are_coops = False
        self.icohpvalue = IcohpValue(
            label=label,
            atom1=atom1,
            atom2=atom2,
            length=length,
            translation=translation,
            num=num,
            icohp=icohp,
            are_coops=are_coops,
        )

        label_sp = "1"
        atom1_sp = "K1"
        atom2_sp = "F2"
        length_sp = "2.3"
        translation_sp = [-1, 0, 0]
        num_sp = 1
        icohp_sp = {Spin.up: -1.1, Spin.down: -1.0}
        are_coops_sp = False
        self.icohpvalue_sp = IcohpValue(
            label=label_sp,
            atom1=atom1_sp,
            atom2=atom2_sp,
            length=length_sp,
            translation=translation_sp,
            num=num_sp,
            icohp=icohp_sp,
            are_coops=are_coops_sp,
        )

    def test_attributes(self):
        # without spin polarization
        assert self.icohpvalue_sp.num_bonds == 1
        assert self.icohpvalue_sp.are_coops is False
        assert self.icohpvalue_sp.is_spin_polarized
        assert self.icohpvalue.icohp == approx({Spin.up: -2.0})

        # with spin polarization
        assert self.icohpvalue_sp.num_bonds == 1
        assert self.icohpvalue_sp.are_coops is False
        assert self.icohpvalue_sp.is_spin_polarized
        assert self.icohpvalue_sp.icohp == approx({Spin.up: -1.1, Spin.down: -1.0})

    def test_icohpvalue(self):
        # without spin polarization
        assert self.icohpvalue.icohpvalue(spin=Spin.up) == approx(-2.0)

        # with spin polarization
        assert self.icohpvalue_sp.icohpvalue(spin=Spin.up) == approx(-1.1)
        assert self.icohpvalue_sp.icohpvalue(spin=Spin.down) == approx(-1.0)

    def test_summed_icohp(self):
        # without spin polarization
        assert self.icohpvalue.summed_icohp == approx(-2.0)

        # with spin polarization
        assert self.icohpvalue_sp.summed_icohp == approx(-2.1)

    def test_str(self):
        # without spin polarization
        assert str(self.icohpvalue) == "ICOHP 1 between K1 and F2 ([-1, 0, 0]): -2.0 eV (Spin up)"

        # with spin polarization
        expected = "ICOHP 1 between K1 and F2 ([-1, 0, 0]): -1.1 eV (Spin up) and -1.0 eV (Spin down)"
        assert str(self.icohpvalue_sp) == expected


class TestCombinedIcohp:
    def setup_method(self):
        # without spin polarization:
        are_coops = are_cobis = is_spin_polarized = False
        list_atom2 = ["K2", "K2", "K2", "K2", "K2", "K2"]
        list_icohp = [
            {Spin.up: -0.40075},
            {Spin.up: -0.40074},
            {Spin.up: -0.40079},
            {Spin.up: -0.40079},
            {Spin.up: -0.40074},
            {Spin.up: -0.40075},
        ]
        list_icoop = [
            {Spin.up: 0.02342},
            {Spin.up: 0.02342},
            {Spin.up: 0.02343},
            {Spin.up: 0.02343},
            {Spin.up: 0.02342},
            {Spin.up: 0.02342},
        ]
        list_labels = ["1", "2", "3", "4", "5", "6"]
        list_length = [2.71199, 2.71199, 2.71199, 2.71199, 2.71199, 2.71199]
        list_num = [1, 1, 1, 1, 1, 1]
        list_atom1 = ["F1", "F1", "F1", "F1", "F1", "F1"]
        list_translation = [
            [0, -1, -1],
            [-1, 0, -1],
            [0, 0, -1],
            [-1, -1, 0],
            [0, -1, 0],
            [-1, 0, 0],
        ]
        self.icohpcollection_KF = IcohpCollection(
            is_spin_polarized=is_spin_polarized,
            are_coops=are_coops,
            are_cobis=are_cobis,
            list_labels=list_labels,
            list_atom1=list_atom1,
            list_atom2=list_atom2,
            list_length=list_length,
            list_translation=list_translation,
            list_num=list_num,
            list_icohp=list_icohp,
        )

        self.icoopcollection_KF = IcohpCollection(
            is_spin_polarized=is_spin_polarized,
            are_coops=True,
            list_labels=list_labels,
            list_atom1=list_atom1,
            list_atom2=list_atom2,
            list_length=list_length,
            list_translation=list_translation,
            list_num=list_num,
            list_icohp=list_icoop,
        )
        self.icohpcollection_orbitalwise = IcohpCollection.from_dict(
            {
                "@module": "pymatgen.electronic_structure.cohp",
                "@class": "IcohpCollection",
                "@version": None,
                "list_labels": ["1", "2"],
                "list_atom1": ["O5", "O5"],
                "list_atom2": ["Ta2", "Ta2"],
                "list_length": [1.99474, 1.99474],
                "list_translation": [[0, 0, -1], [0, 0, 0]],
                "list_num": [1, 1],
                "list_icohp": [
                    {"1": 0.29324, "-1": 0.29324},
                    {"1": 0.29324, "-1": 0.29324},
                ],
                "is_spin_polarized": True,
                "list_orb_icohp": [
                    {
                        "2s-6s": {
                            "icohp": {"1": 0.0247, "-1": 0.0247},
                            "orbitals": [[2, 0], [6, 0]],
                        },
                        "2s-5py": {
                            "icohp": {"1": 8e-05, "-1": 8e-05},
                            "orbitals": [[2, 0], [5, 1]],
                        },
                    },
                    {
                        "2s-6s": {
                            "icohp": {"1": 0.0247, "-1": 0.0247},
                            "orbitals": [[2, 0], [6, 0]],
                        },
                        "2s-5py": {
                            "icohp": {"1": 0.5, "-1": 0},
                            "orbitals": [[2, 0], [5, 1]],
                        },
                    },
                ],
                "are_coops": False,
                "are_cobis": True,
            }
        )
        # with spin polarization:
        list_atom2_sp = ["Fe7", "Fe9"]
        list_labels_sp = ["1", "2"]
        list_translation_sp = [[0, 0, 0], [0, 0, 0]]
        list_length_sp = [2.83189, 2.45249]
        list_atom1_sp = ["Fe8", "Fe8"]
        is_spin_polarized_sp = True
        are_coops_sp = False
        list_num_sp = [2, 1]
        list_icohp_sp = [
            {Spin.up: -0.10218, Spin.down: -0.19701},
            {Spin.up: -0.28485, Spin.down: -0.58279},
        ]
        list_icoop_sp = [
            {Spin.up: -0.11389, Spin.down: -0.20828},
            {Spin.up: -0.04087, Spin.down: -0.05756},
        ]

        self.icohpcollection_Fe = IcohpCollection(
            is_spin_polarized=is_spin_polarized_sp,
            are_coops=are_coops_sp,
            are_cobis=False,
            list_labels=list_labels_sp,
            list_atom1=list_atom1_sp,
            list_atom2=list_atom2_sp,
            list_length=list_length_sp,
            list_translation=list_translation_sp,
            list_num=list_num_sp,
            list_icohp=list_icohp_sp,
        )
        self.icoopcollection_Fe = IcohpCollection(
            is_spin_polarized=is_spin_polarized_sp,
            are_coops=True,
            list_labels=list_labels_sp,
            list_atom1=list_atom1_sp,
            list_atom2=list_atom2_sp,
            list_length=list_length_sp,
            list_translation=list_translation_sp,
            list_num=list_num_sp,
            list_icohp=list_icoop_sp,
        )

    def test_get_icohp_by_label(self):
        # without spin polarization

        # ICOHPs
        assert self.icohpcollection_KF.get_icohp_by_label("1") == approx(-0.40075)
        assert self.icohpcollection_KF.get_icohp_by_label("2") == approx(-0.40074)
        assert self.icohpcollection_KF.get_icohp_by_label("3") == approx(-0.40079)
        assert self.icohpcollection_KF.get_icohp_by_label("4") == approx(-0.40079)
        assert self.icohpcollection_KF.get_icohp_by_label("5") == approx(-0.40074)
        assert self.icohpcollection_KF.get_icohp_by_label("6") == approx(-0.40075)

        # with spin polarization
        # summed spin
        # ICOHPs
        assert self.icohpcollection_Fe.get_icohp_by_label("1") == approx(-0.10218 - 0.19701)
        assert self.icohpcollection_Fe.get_icohp_by_label("2") == approx(-0.28485 - 0.58279)

        # Spin up
        # ICOHPs
        assert self.icohpcollection_Fe.get_icohp_by_label("1", summed_spin_channels=False) == approx(-0.10218)
        assert self.icohpcollection_Fe.get_icohp_by_label("2", summed_spin_channels=False) == approx(-0.28485)

        # Spin down
        # ICOHPs
        assert self.icohpcollection_Fe.get_icohp_by_label("1", summed_spin_channels=False, spin=Spin.down) == approx(
            -0.19701
        )
        assert self.icohpcollection_Fe.get_icohp_by_label("2", summed_spin_channels=False, spin=Spin.down) == approx(
            -0.58279
        )

        # orbitalwise
        assert self.icohpcollection_orbitalwise.get_icohp_by_label("1", orbitals="2s-6s") == approx(0.0494)
        assert self.icohpcollection_orbitalwise.get_icohp_by_label(
            "1", orbitals="2s-6s", spin=Spin.up, summed_spin_channels=False
        ) == approx(0.0247)
        assert self.icohpcollection_orbitalwise.get_icohp_by_label(
            "1", orbitals="2s-6s", spin=Spin.down, summed_spin_channels=False
        ) == approx(0.0247)
        assert self.icohpcollection_orbitalwise.get_icohp_by_label(
            "2", orbitals="2s-5py", spin=Spin.up, summed_spin_channels=False
        ) == approx(0.5)

    def test_get_summed_icohp_by_label_list(self):
        # without spin polarization
        assert self.icohpcollection_KF.get_summed_icohp_by_label_list(
            ["1", "2", "3", "4", "5", "6"], divisor=6.0
        ) == approx(-0.40076)

        # with spin polarization
        sum1 = (-0.10218 - 0.19701 - 0.28485 - 0.58279) / 2.0
        sum2 = (-0.10218 - 0.28485) / 2.0
        sum3 = (-0.19701 - 0.58279) / 2.0
        assert self.icohpcollection_Fe.get_summed_icohp_by_label_list(["1", "2"], divisor=2.0) == approx(sum1)
        assert self.icohpcollection_Fe.get_summed_icohp_by_label_list(
            ["1", "2"], summed_spin_channels=False, divisor=2.0
        ) == approx(sum2)
        assert self.icohpcollection_Fe.get_summed_icohp_by_label_list(
            ["1", "2"], summed_spin_channels=False, spin=Spin.down, divisor=2.0
        ) == approx(sum3)

    def test_get_icohp_dict_by_bondlengths(self):
        # without spin polarization
        icohpvalue = {}
        icohpvalue["1"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40075},
            "are_coops": False,
            "are_cobis": False,
            "label": "1",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [0, -1, -1],
            "orbitals": None,
        }
        icohpvalue["2"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40074},
            "are_coops": False,
            "are_cobis": False,
            "label": "2",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [-1, 0, -1],
            "orbitals": None,
        }
        icohpvalue["3"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40079},
            "are_coops": False,
            "are_cobis": False,
            "label": "3",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [0, 0, -1],
            "orbitals": None,
        }
        icohpvalue["4"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40079},
            "are_coops": False,
            "are_cobis": False,
            "label": "4",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [-1, -1, 0],
            "orbitals": None,
        }
        icohpvalue["5"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40074},
            "are_coops": False,
            "are_cobis": False,
            "label": "5",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [0, -1, 0],
            "orbitals": None,
        }
        icohpvalue["6"] = {
            "@module": "pymatgen.electronic_structure.cohp",
            "num": 1,
            "length": 2.71199,
            "icohp": {Spin.up: -0.40075},
            "are_coops": False,
            "are_cobis": False,
            "label": "6",
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "translation": [-1, 0, 0],
            "orbitals": None,
        }

        dict_KF = self.icohpcollection_KF.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=8.0)
        for key, value in sorted(dict_KF.items()):
            v = value.as_dict()
            if "@version" in v:
                v.pop("@version")
            assert v == icohpvalue[key]

        assert self.icohpcollection_KF.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=1.0) == {}

        # with spin polarization
        icohpvalue_spin = {}
        icohpvalue_spin["1"] = {
            "num": 2,
            "atom2": "Fe7",
            "translation": [0, 0, 0],
            "@module": "pymatgen.electronic_structure.cohp",
            "are_coops": False,
            "are_cobis": False,
            "atom1": "Fe8",
            "label": "1",
            "length": 2.83189,
            "@class": "IcohpValue",
            "icohp": {Spin.up: -0.10218, Spin.down: -0.19701},
            "orbitals": None,
        }
        icohpvalue_spin["2"] = {
            "num": 1,
            "atom2": "Fe9",
            "translation": [0, 0, 0],
            "@module": "pymatgen.electronic_structure.cohp",
            "are_coops": False,
            "are_cobis": False,
            "atom1": "Fe8",
            "label": "2",
            "length": 2.45249,
            "@class": "IcohpValue",
            "icohp": {Spin.up: -0.28485, Spin.down: -0.58279},
            "orbitals": None,
        }

        dict_Fe = self.icohpcollection_Fe.get_icohp_dict_by_bondlengths(minbondlength=0.0, maxbondlength=8.0)
        for key, value in sorted(dict_Fe.items()):
            v = value.as_dict()
            if "@version" in v:
                v.pop("@version")
            assert v == icohpvalue_spin[key]

        dict_Fe2 = self.icohpcollection_Fe.get_icohp_dict_by_bondlengths(minbondlength=2.5, maxbondlength=2.9)
        assert len(dict_Fe2) == 1
        for key, value in sorted(dict_Fe2.items()):
            v = value.as_dict()
            if "@version" in v:
                v.pop("@version")
            assert v == icohpvalue_spin[key]

    def test_get_icohp_dict_of_site(self):
        # without spin polarization
        icohpvalue = {}
        icohpvalue["1"] = {
            "translation": [0, -1, -1],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "1",
            "icohp": {Spin.up: -0.40075},
            "orbitals": None,
        }
        icohpvalue["2"] = {
            "translation": [-1, 0, -1],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "2",
            "icohp": {Spin.up: -0.40074},
            "orbitals": None,
        }
        icohpvalue["3"] = {
            "translation": [0, 0, -1],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "3",
            "icohp": {Spin.up: -0.40079},
            "orbitals": None,
        }
        icohpvalue["4"] = {
            "translation": [-1, -1, 0],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "4",
            "icohp": {Spin.up: -0.40079},
            "orbitals": None,
        }
        icohpvalue["5"] = {
            "translation": [0, -1, 0],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "5",
            "icohp": {Spin.up: -0.40074},
            "orbitals": None,
        }
        icohpvalue["6"] = {
            "translation": [-1, 0, 0],
            "are_coops": False,
            "are_cobis": False,
            "@module": "pymatgen.electronic_structure.cohp",
            "length": 2.71199,
            "atom2": "K2",
            "@class": "IcohpValue",
            "atom1": "F1",
            "num": 1,
            "label": "6",
            "icohp": {Spin.up: -0.40075},
            "orbitals": None,
        }

        dict_KF = self.icohpcollection_KF.get_icohp_dict_of_site(site=0)

        for key, value in sorted(dict_KF.items()):
            v = value.as_dict()
            if "@version" in v:
                v.pop("@version")
            assert v == icohpvalue[key]

        # compare number of results dependent on minsummedicohp, maxsummedicohp,minbondlength, maxbondlength, and
        # only_bonds_to
        dict_KF_2 = self.icohpcollection_KF.get_icohp_dict_of_site(
            site=0,
            minsummedicohp=None,
            maxsummedicohp=-0.0,
            minbondlength=0.0,
            maxbondlength=8.0,
        )
        dict_KF_3 = self.icohpcollection_KF.get_icohp_dict_of_site(
            site=0,
            minsummedicohp=None,
            maxsummedicohp=-0.5,
            minbondlength=0.0,
            maxbondlength=8.0,
        )
        dict_KF_4 = self.icohpcollection_KF.get_icohp_dict_of_site(
            site=0,
            minsummedicohp=0.0,
            maxsummedicohp=None,
            minbondlength=0.0,
            maxbondlength=8.0,
        )
        dict_KF_5 = self.icohpcollection_KF.get_icohp_dict_of_site(
            site=0,
            minsummedicohp=None,
            maxsummedicohp=None,
            minbondlength=0.0,
            maxbondlength=2.0,
        )
        dict_KF_6 = self.icohpcollection_KF.get_icohp_dict_of_site(
            site=0,
            minsummedicohp=None,
            maxsummedicohp=None,
            minbondlength=3.0,
            maxbondlength=8.0,
        )
        dict_KF_7 = self.icohpcollection_KF.get_icohp_dict_of_site(site=0, only_bonds_to=["K"])
        dict_KF_8 = self.icohpcollection_KF.get_icohp_dict_of_site(site=1, only_bonds_to=["K"])
        dict_KF_9 = self.icohpcollection_KF.get_icohp_dict_of_site(site=1, only_bonds_to=["F"])

        assert len(dict_KF_2) == 6
        assert len(dict_KF_3) == 0
        assert len(dict_KF_4) == 0
        assert len(dict_KF_5) == 0
        assert len(dict_KF_6) == 0
        assert len(dict_KF_7) == 6
        assert len(dict_KF_8) == 0
        assert len(dict_KF_9) == 6

        # spin polarization

        dict_Fe = self.icohpcollection_Fe.get_icohp_dict_of_site(site=0)
        assert len(dict_Fe) == 0

        # Fe8
        dict_Fe2 = self.icohpcollection_Fe.get_icohp_dict_of_site(site=7)
        assert len(dict_Fe2) == 2
        # Test the values

        icohplist_Fe = {}
        icohplist_Fe["1"] = {
            "are_coops": False,
            "are_cobis": False,
            "translation": [0, 0, 0],
            "icohp": {Spin.down: -0.19701, Spin.up: -0.10218},
            "length": 2.83189,
            "@module": "pymatgen.electronic_structure.cohp",
            "atom1": "Fe8",
            "atom2": "Fe7",
            "label": "1",
            "orbitals": None,
            "@class": "IcohpValue",
            "num": 2,
        }
        icohplist_Fe["2"] = {
            "are_coops": False,
            "are_cobis": False,
            "translation": [0, 0, 0],
            "icohp": {Spin.down: -0.58279, Spin.up: -0.28485},
            "length": 2.45249,
            "@module": "pymatgen.electronic_structure.cohp",
            "atom1": "Fe8",
            "atom2": "Fe9",
            "label": "2",
            "orbitals": None,
            "@class": "IcohpValue",
            "num": 1,
        }

        for key, value in sorted(dict_Fe2.items()):
            v = value.as_dict()
            if "@version" in v:
                v.pop("@version")
            assert v == icohplist_Fe[key]

        # Fe9
        dict_Fe3 = self.icohpcollection_Fe.get_icohp_dict_of_site(site=8)
        assert len(dict_Fe3) == 1

        # compare number of results dependent on minsummedicohp, maxsummedicohp,minbondlength, maxbondlength
        # Fe8
        dict_Fe4 = self.icohpcollection_Fe.get_icohp_dict_of_site(
            site=7,
            minsummedicohp=-0.3,
            maxsummedicohp=None,
            minbondlength=0.0,
            maxbondlength=8.0,
        )
        assert len(dict_Fe4) == 1
        values = list(dict_Fe4.values())
        v = values[0].as_dict()
        if "@version" in v:
            v.pop("@version")
        assert v == icohplist_Fe["1"]

        dict_Fe5 = self.icohpcollection_Fe.get_icohp_dict_of_site(
            site=7,
            minsummedicohp=None,
            maxsummedicohp=-0.3,
            minbondlength=0.0,
            maxbondlength=8.0,
        )
        assert len(dict_Fe5) == 1
        values = list(dict_Fe5.values())
        v = values[0].as_dict()
        if "@version" in v:
            v.pop("@version")
        assert v == icohplist_Fe["2"]

        dict_Fe6 = self.icohpcollection_Fe.get_icohp_dict_of_site(
            site=7,
            minsummedicohp=None,
            maxsummedicohp=None,
            minbondlength=0.0,
            maxbondlength=2.5,
        )

        assert len(dict_Fe6) == 1
        values = list(dict_Fe6.values())
        v = values[0].as_dict()
        if "@version" in v:
            v.pop("@version")
        assert v == icohplist_Fe["2"]

        dict_Fe7 = self.icohpcollection_Fe.get_icohp_dict_of_site(
            site=7,
            minsummedicohp=None,
            maxsummedicohp=None,
            minbondlength=2.5,
            maxbondlength=8.0,
        )
        assert len(dict_Fe7) == 1
        values = list(dict_Fe7.values())
        v = values[0].as_dict()
        if "@version" in v:
            v.pop("@version")
        assert v == icohplist_Fe["1"]

    def test_extremum_icohpvalue(self):
        # without spin polarization
        # ICOHPs
        assert self.icohpcollection_KF.extremum_icohpvalue() == approx(-0.40079)
        # ICOOPs
        assert self.icoopcollection_KF.extremum_icohpvalue() == approx(0.02343)
        # with spin polarization
        # summed spin
        # ICOHPs
        assert self.icohpcollection_Fe.extremum_icohpvalue() == approx(-0.86764)
        assert self.icoopcollection_Fe.extremum_icohpvalue() == approx(-0.09842999999999999)
        # ICOOPs
        # spin up
        # ICOHPs
        assert self.icohpcollection_Fe.extremum_icohpvalue(summed_spin_channels=False) == approx(-0.28485)
        # ICOOPs
        assert self.icoopcollection_Fe.extremum_icohpvalue(summed_spin_channels=False) == approx(-0.04087)
        # spin down
        # ICOHPs
        assert self.icohpcollection_Fe.extremum_icohpvalue(summed_spin_channels=False, spin=Spin.down) == approx(
            -0.58279
        )
        # ICOOPs
        assert self.icoopcollection_Fe.extremum_icohpvalue(summed_spin_channels=False, spin=Spin.down) == approx(
            -0.05756
        )


class TestCompleteCohp(MatSciTest):
    def setup_method(self):
        filepath = f"{TEST_DIR}/complete_cohp_lobster.json"
        with open(filepath, "rb") as file:
            self.cohp_lobster_dict = CompleteCohp.from_dict(orjson.loads(file.read()))
        filepath = f"{TEST_DIR}/complete_coop_lobster.json"
        with open(filepath, "rb") as file:
            self.coop_lobster_dict = CompleteCohp.from_dict(orjson.loads(file.read()))
        filepath = f"{TEST_DIR}/complete_cohp_lmto.json"
        with open(filepath, "rb") as file:
            self.cohp_lmto_dict = CompleteCohp.from_dict(orjson.loads(file.read()))
        filepath = f"{TEST_DIR}/complete_cohp_orbitalwise.json"
        with open(filepath, "rb") as file:
            self.cohp_orb_dict = CompleteCohp.from_dict(orjson.loads(file.read()))
        # Lobster 3.0
        filepath = f"{TEST_DIR}/complete_cohp_forb.json"
        with open(filepath, "rb") as file:
            self.cohp_lobster_forb_dict = CompleteCohp.from_dict(orjson.loads(file.read()))

            # Lobster 2.0
        filepath = f"{TEST_DIR}/COPL.BiSe"
        structure = f"{TEST_DIR}/CTRL.BiSe"
        self.cohp_lmto = CompleteCohp.from_file("lmto", filename=filepath, structure_file=structure)
        filepath = f"{TEST_DIR}/COHPCAR.lobster.gz"
        structure = f"{TEST_DIR}/POSCAR"
        self.cohp_lobster = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)
        # with open(f"{TEST_DIR}/complete_cohp_lobster.json", "w", encoding="utf-8") as file:
        #     json.dump(self.cohp_lobster.as_dict(), file)
        filepath = f"{TEST_DIR}/COOPCAR.lobster.BiSe.gz"
        structure = f"{TEST_DIR}/POSCAR.BiSe"
        self.coop_lobster = CompleteCohp.from_file(
            "lobster", filename=filepath, structure_file=structure, are_coops=True
        )
        filepath = f"{TEST_DIR}/COHPCAR.lobster.orbitalwise.gz"
        structure = f"{TEST_DIR}/POSCAR.orbitalwise"
        self.cohp_orb = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)
        # with open(f"{TEST_DIR}/complete_cohp_orbitalwise.json", "w", encoding="utf-8") as file:
        #     json.dump(self.cohp_orb.as_dict(), file)
        filepath = f"{TEST_DIR}/COHPCAR.lobster.notot.orbitalwise.gz"
        self.cohp_notot = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)
        # Lobster 3.0
        filepath = f"{TEST_DIR}/COHPCAR.lobster.Na2UO4.gz"
        structure = f"{TEST_DIR}/POSCAR.Na2UO4"
        self.cohp_lobster_forb = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)

        # spinpolarized case:
        filepath = f"{TEST_DIR}/environments/COHPCAR.lobster.mp-190.gz"
        structure = f"{TEST_DIR}/environments/POSCAR.mp_190.gz"
        self.cohp_lobster_spin_polarized = CompleteCohp.from_file(
            "lobster", filename=filepath, structure_file=structure
        )
        # COBI
        filepath = f"{TEST_DIR}/COBICAR.lobster.gz"
        structure = f"{TEST_DIR}/POSCAR.COBI"

        self.cobi = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure, are_cobis=True)

        # COBI multi-center
        filepath = f"{TEST_DIR}/COBICAR.lobster.GeTe.multi.orbitalwise.full"
        structure = f"{TEST_DIR}/POSCAR.GeTe"
        self.cobi_multi = CompleteCohp.from_file(
            "lobster",
            filename=filepath,
            structure_file=structure,
            are_multi_center_cobis=True,
        )

        # COBI multi-center
        filepath = f"{TEST_DIR}/COBICAR.lobster.B2H6.spin"
        structure = f"{TEST_DIR}/POSCAR.B2H6"
        self.cobi_multi_B2H6 = CompleteCohp.from_file(
            "lobster",
            filename=filepath,
            structure_file=structure,
            are_multi_center_cobis=True,
        )

        # COBI multi-center
        filepath = f"{TEST_DIR}/COBICAR.lobster.B2H6.spin.average.2"
        structure = f"{TEST_DIR}/POSCAR.B2H6"
        self.cobi_multi_B2H6_average2 = CompleteCohp.from_file(
            "lobster", filename=filepath, structure_file=structure, are_cobis=True
        )

    def test_attributes(self):
        assert not self.cohp_lobster.are_coops
        assert not self.cohp_lobster.are_cobis
        assert not self.cohp_lobster_dict.are_coops
        assert not self.cohp_lmto.are_coops
        assert not self.cohp_lmto_dict.are_coops
        assert self.coop_lobster.are_coops
        assert self.coop_lobster_dict.are_coops
        assert not self.cohp_lobster_forb.are_coops
        assert not self.cohp_lobster_forb_dict.are_coops

        assert len(self.cohp_lobster.energies) == 301
        assert len(self.cohp_lmto.energies) == 801
        assert len(self.coop_lobster.energies) == 241
        assert len(self.cohp_lobster_forb.energies) == 7

        assert self.cohp_lobster.efermi == approx(9.75576)
        assert self.cohp_lmto.efermi == approx(-2.3433)
        assert self.coop_lobster.efermi == approx(5.90043)
        assert self.cohp_lobster_forb.efermi == approx(4.12875)

        assert self.cobi.are_cobis
        assert not self.cobi.are_coops

        assert self.cohp_lobster_forb.cohp[Spin.up][0] == approx(0.00000)
        assert self.cohp_lobster_forb.icohp[Spin.up][0] == approx(-0.09040)

    def test_average_multi_center_cobi(self):
        # tests if the averages for a mult-center cobi are computed in the same way as in Lobster
        for cohp1, cohp2 in zip(
            self.cobi_multi_B2H6.get_cohp_by_label("average").cohp[Spin.up],
            self.cobi_multi_B2H6_average2.get_cohp_by_label("average").cohp[Spin.up],
            strict=True,
        ):
            assert cohp1 == approx(cohp2, abs=1e-4)

        for cohp1, cohp2 in zip(
            self.cobi_multi_B2H6.get_cohp_by_label("average").cohp[Spin.down],
            self.cobi_multi_B2H6_average2.get_cohp_by_label("average").cohp[Spin.down],
            strict=True,
        ):
            assert cohp1 == approx(cohp2, abs=1e-4)

        for icohp1, icohp2 in zip(
            self.cobi_multi_B2H6.get_cohp_by_label("average").icohp[Spin.up],
            self.cobi_multi_B2H6_average2.get_cohp_by_label("average").icohp[Spin.up],
            strict=True,
        ):
            assert icohp1 == approx(icohp2, abs=1e-4)

        for icohp1, icohp2 in zip(
            self.cobi_multi_B2H6.get_cohp_by_label("average").icohp[Spin.down],
            self.cobi_multi_B2H6_average2.get_cohp_by_label("average").icohp[Spin.down],
            strict=True,
        ):
            assert icohp1 == approx(icohp2, abs=1e-4)

    def test_dict(self):
        # The JSON files are dict representations of the COHPs from the LMTO
        # and LOBSTER calculations and should thus be the same.

        def is_equal(a, b):
            a_dict = a.as_dict()
            b_dict = b.as_dict()
            del a_dict["structure"]
            del b_dict["structure"]
            return a_dict == b_dict and a.structure == b.structure

        assert is_equal(self.cohp_lobster, self.cohp_lobster_dict)
        assert is_equal(self.cohp_orb, self.cohp_orb_dict)
        # Lobster 3.0, including f orbitals

        assert is_equal(self.cohp_lobster_forb, self.cohp_lobster_forb_dict)

        # Testing the LMTO dicts will be more involved. Since the average
        # is calculated and not read, there may be differences in rounding
        # with a very small number of matrix elements, which would cause the
        # test to fail
        cohp_lmto_dict = self.cohp_lmto.as_dict()
        for key in ["COHP", "ICOHP"]:
            assert_allclose(
                cohp_lmto_dict[key]["average"]["1"],
                self.cohp_lmto_dict.as_dict()[key]["average"]["1"],
                5,
            )
            # check if the same dicts are generated
            cobi_new = CompleteCohp.from_dict(self.cobi_multi.as_dict())
            assert is_equal(self.cobi_multi, cobi_new)
        # for key in cohp_lmto_dict:
        #     if key not in ["COHP", "ICOHP"]:
        #         assert cohp_lmto_dict[key] == self.cohp_lmto_dict.as_dict()[key]
        #     else:
        #         for bond in cohp_lmto_dict[key]:
        #             if bond != "average":
        #                 assert cohp_lmto_dict[key][bond] == self.cohp_lmto_dict.as_dict()[key][bond]

    def test_icohp_values(self):
        # icohp_ef are the ICHOP(Ef) values taken from
        # the ICOHPLIST.lobster file.
        icohp_ef_dict = {
            "1": {Spin.up: -0.10218, Spin.down: -0.19701},
            "2": {Spin.up: -0.28485, Spin.down: -0.58279},
        }
        all_cohps_lobster = self.cohp_lobster.all_cohps
        for bond, val in icohp_ef_dict.items():
            icohp_ef = all_cohps_lobster[bond].get_interpolated_value(self.cohp_lobster.efermi, integrated=True)
            assert val == icohp_ef

        icoop_ef_dict = {
            "1": {Spin.up: 0.14245},
            "2": {Spin.up: -0.04118},
            "3": {Spin.up: 0.14245},
            "4": {Spin.up: -0.04118},
            "5": {Spin.up: -0.03516},
            "6": {Spin.up: 0.10745},
            "7": {Spin.up: -0.03516},
            "8": {Spin.up: 0.10745},
            "9": {Spin.up: -0.12395},
            "10": {Spin.up: 0.24714},
            "11": {Spin.up: -0.12395},
        }
        all_coops_lobster = self.coop_lobster.all_cohps
        for bond, val in icoop_ef_dict.items():
            icoop_ef = all_coops_lobster[bond].get_interpolated_value(self.coop_lobster.efermi, integrated=True)
            assert val == icoop_ef

    def test_get_cohp_by_label(self):
        assert self.cohp_orb.get_cohp_by_label("1").energies[0] == approx(-11.7225)
        assert self.cohp_orb.get_cohp_by_label("1").energies[5] == approx(-11.47187)
        assert not self.cohp_orb.get_cohp_by_label("1").are_coops
        assert self.cohp_orb.get_cohp_by_label("1").cohp[Spin.up][0] == approx(0.0)
        assert self.cohp_orb.get_cohp_by_label("1").cohp[Spin.up][300] == approx(0.03392)
        assert self.cohp_orb.get_cohp_by_label("average").cohp[Spin.up][230] == approx(-0.08792)
        assert self.cohp_orb.get_cohp_by_label("average").energies[230] == approx(-0.19368000000000007)
        assert not self.cohp_orb.get_cohp_by_label("average").are_coops
        # test methods from super class that could be overwritten
        assert self.cohp_orb.get_icohp()[Spin.up][3] == approx(0.0)
        assert self.cohp_orb.get_cohp()[Spin.up][3] == approx(0.0)

    def test_get_cohp_by_label_summed_spin(self):
        # files without spin polarization
        assert self.cohp_orb.get_cohp_by_label("1", summed_spin_channels=True).energies[0] == approx(-11.7225)
        assert self.cohp_orb.get_cohp_by_label("1", summed_spin_channels=True).energies[5] == approx(-11.47187)
        assert not self.cohp_orb.get_cohp_by_label("1", summed_spin_channels=True).are_coops
        assert self.cohp_orb.get_cohp_by_label("1", summed_spin_channels=True).cohp[Spin.up][0] == approx(0.0)
        assert self.cohp_orb.get_cohp_by_label("1", summed_spin_channels=True).cohp[Spin.up][300] == approx(0.03392)
        assert self.cohp_orb.get_cohp_by_label("average", summed_spin_channels=True).cohp[Spin.up][230] == approx(
            -0.08792
        )
        assert self.cohp_orb.get_cohp_by_label("average", summed_spin_channels=True).energies[230] == approx(
            -0.19368000000000007
        )
        assert not self.cohp_orb.get_cohp_by_label("average", summed_spin_channels=True).are_coops

        # file with spin polarization
        assert self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=False).cohp[Spin.up][
            300
        ] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=True).cohp[Spin.up][300]
        )
        assert self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=False).cohp[Spin.down][
            300
        ] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=True).cohp[Spin.up][300]
        )
        assert self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=True).energies[0] == approx(
            -15.03759 + 1.96204
        )
        assert self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=True).energies[5] == approx(
            -14.78697 + 1.96204
        )
        assert not self.cohp_lobster_spin_polarized.get_cohp_by_label("1", summed_spin_channels=True).are_coops

    def test_get_summed_cohp_by_label_list(self):
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"]).energies[0] == approx(-11.7225)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"]).energies[0] == approx(-11.7225)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"]).energies[5] == approx(-11.47187)
        assert not self.cohp_orb.get_summed_cohp_by_label_list(["1"]).are_coops
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"]).cohp[Spin.up][0] == approx(0.0)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"]).cohp[Spin.up][0] == approx(0.0)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"]).cohp[Spin.up][300] == approx(0.03392 * 2.0)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"], divisor=2).cohp[Spin.up][300] == approx(0.03392)

    def test_get_summed_cohp_by_label_list_summed_spin(self):
        # files without spin polarization
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).energies[0] == approx(
            -11.7225
        )
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"], summed_spin_channels=True).energies[0] == approx(
            -11.7225
        )
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).energies[5] == approx(
            -11.47187
        )
        assert not self.cohp_orb.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).are_coops
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).cohp[Spin.up][0] == approx(
            0.0
        )
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"], summed_spin_channels=True).cohp[Spin.up][
            0
        ] == approx(0.0)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"], summed_spin_channels=True).cohp[Spin.up][
            300
        ] == approx(0.03392 * 2.0)
        assert self.cohp_orb.get_summed_cohp_by_label_list(["1", "1"], summed_spin_channels=True, divisor=2).cohp[
            Spin.up
        ][300] == approx(0.03392)

        # file with spin polarization
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(["1"], summed_spin_channels=False).cohp[
            Spin.up
        ][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).cohp[
                Spin.up
            ][300]
        )
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(["1"], summed_spin_channels=False).cohp[
            Spin.down
        ][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(["1"], summed_spin_channels=True).cohp[
                Spin.up
            ][300]
        )
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(
            ["1", "1"], summed_spin_channels=True
        ).energies[0] == approx(-15.03759 + 1.96204)
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(
            ["1"], summed_spin_channels=True
        ).energies[5] == approx(-14.78697 + 1.96204)
        assert not self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_list(
            ["1"], summed_spin_channels=True
        ).are_coops

    def test_get_summed_cohp_by_label_and_orbital_list(self):
        ref = self.cohp_orb.orb_res_cohp["1"]["4s-4px"]
        ref2 = self.cohp_orb.orb_res_cohp["1"]["4px-4pz"]
        cohp_label = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1"], ["4s-4px"])
        cohp_label2 = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1", "1"], ["4s-4px", "4s-4px"])
        cohp_label2x = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
            ["1", "1"], ["4s-4px", "4s-4px"], divisor=2
        )
        cohp_label3 = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1", "1"], ["4px-4pz", "4s-4px"])

        assert_allclose(cohp_label.cohp[Spin.up], ref["COHP"][Spin.up])
        assert_allclose(cohp_label2.cohp[Spin.up], ref["COHP"][Spin.up] * 2.0)
        assert_allclose(cohp_label3.cohp[Spin.up], ref["COHP"][Spin.up] + ref2["COHP"][Spin.up])
        assert_allclose(cohp_label.icohp[Spin.up], ref["ICOHP"][Spin.up])
        assert_allclose(cohp_label2.icohp[Spin.up], ref["ICOHP"][Spin.up] * 2.0)
        assert_allclose(cohp_label2x.icohp[Spin.up], ref["ICOHP"][Spin.up])
        assert_allclose(cohp_label3.icohp[Spin.up], ref["ICOHP"][Spin.up] + ref2["ICOHP"][Spin.up])
        expected_msg = "label_list and orbital_list don't have the same length"
        with pytest.raises(ValueError, match=expected_msg):
            self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1"], ["4px-4pz", "4s-4px"])
        with pytest.raises(ValueError, match=expected_msg):
            self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1", "2"], ["4s-4px"])

    def test_get_summed_cohp_by_label_and_orbital_list_summed_spin_channels(self):
        ref = self.cohp_orb.orb_res_cohp["1"]["4s-4px"]
        ref2 = self.cohp_orb.orb_res_cohp["1"]["4px-4pz"]
        cohp_label = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["4s-4px"], summed_spin_channels=True
        )
        cohp_label2 = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
            ["1", "1"], ["4s-4px", "4s-4px"], summed_spin_channels=True
        )
        cohp_label2x = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
            ["1", "1"], ["4s-4px", "4s-4px"], divisor=2, summed_spin_channels=True
        )
        cohp_label3 = self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
            ["1", "1"], ["4px-4pz", "4s-4px"], summed_spin_channels=True
        )

        assert_allclose(cohp_label.cohp[Spin.up], ref["COHP"][Spin.up])
        assert_allclose(cohp_label2.cohp[Spin.up], ref["COHP"][Spin.up] * 2.0)
        assert_allclose(cohp_label3.cohp[Spin.up], ref["COHP"][Spin.up] + ref2["COHP"][Spin.up])
        assert_allclose(cohp_label.icohp[Spin.up], ref["ICOHP"][Spin.up])
        assert_allclose(cohp_label2.icohp[Spin.up], ref["ICOHP"][Spin.up] * 2.0)
        assert_allclose(cohp_label2x.icohp[Spin.up], ref["ICOHP"][Spin.up])
        assert_allclose(cohp_label3.icohp[Spin.up], ref["ICOHP"][Spin.up] + ref2["ICOHP"][Spin.up])
        expected_msg = "label_list and orbital_list don't have the same length"
        with pytest.raises(ValueError, match=expected_msg):
            self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(
                ["1"], ["4px-4pz", "4s-4px"], summed_spin_channels=True
            )
        with pytest.raises(ValueError, match=expected_msg):
            self.cohp_orb.get_summed_cohp_by_label_and_orbital_list(["1", "2"], ["4s-4px"], summed_spin_channels=True)

        # files with spin polarization
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["6s-6s"], summed_spin_channels=False
        ).cohp[Spin.up][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
                ["1"], ["6s-6s"], summed_spin_channels=True
            ).cohp[Spin.up][300]
        )
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["6s-6s"], summed_spin_channels=False
        ).cohp[Spin.down][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
                ["1"], ["6s-6s"], summed_spin_channels=True
            ).cohp[Spin.up][300]
        )
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["6s-6s"], summed_spin_channels=True
        ).energies[0] == approx(-15.03759 + 1.96204)
        assert self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["6s-6s"], summed_spin_channels=True
        ).energies[5] == approx(-14.78697 + 1.96204)
        assert not self.cohp_lobster_spin_polarized.get_summed_cohp_by_label_and_orbital_list(
            ["1"], ["6s-6s"], summed_spin_channels=True
        ).are_coops

    def test_orbital_resolved_cohp(self):
        # When read from a COHPCAR file, total COHPs are calculated from
        # the orbital-resolved COHPs if the total is missing. This may be
        # case for LOBSTER version 2.2.0 and earlier due to a bug with the
        # cohpgenerator keyword. The calculated total should be approximately
        # the total COHP calculated by LOBSTER. Due to numerical errors in
        # the LOBSTER calculation, the precision is not very high though.

        assert_allclose(
            self.cohp_orb.all_cohps["1"].cohp[Spin.up],
            self.cohp_notot.all_cohps["1"].cohp[Spin.up],
            atol=1e-3,
        )
        assert_allclose(
            self.cohp_orb.all_cohps["1"].icohp[Spin.up],
            self.cohp_notot.all_cohps["1"].icohp[Spin.up],
            atol=1e-3,
        )

        # Tests different methods for getting orbital-resolved COHPs
        ref = self.cohp_orb.orb_res_cohp["1"]["4s-4px"]
        cohp_label = self.cohp_orb.get_orbital_resolved_cohp("1", "4s-4px")
        assert cohp_label.cohp == ref["COHP"]
        assert cohp_label.icohp == ref["ICOHP"]
        orbitals = [[Orbital.s, Orbital.px], ["s", "px"], [0, 3]]
        cohps = [self.cohp_orb.get_orbital_resolved_cohp("1", [[4, orb[0]], [4, orb[1]]]) for orb in orbitals]
        for cohp in cohps:
            assert cohp.as_dict() == cohp_label.as_dict()

    def test_orbital_resolved_cohp_summed_spin_channels(self):
        ref = self.cohp_orb.orb_res_cohp["1"]["4s-4px"]
        cohp_label = self.cohp_orb.get_orbital_resolved_cohp("1", "4s-4px", summed_spin_channels=True)
        assert cohp_label.cohp == ref["COHP"]
        assert cohp_label.icohp == ref["ICOHP"]
        orbitals = [[Orbital.s, Orbital.px], ["s", "px"], [0, 3]]
        cohps = [
            self.cohp_orb.get_orbital_resolved_cohp("1", [[4, orb[0]], [4, orb[1]]], summed_spin_channels=True)
            for orb in orbitals
        ]

        for cohp in cohps:
            assert cohp.as_dict() == cohp_label.as_dict()

        # spin polarization
        assert self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp(
            "1", "6s-6s", summed_spin_channels=False
        ).cohp[Spin.up][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp("1", "6s-6s", summed_spin_channels=True).cohp[
                Spin.up
            ][300]
        )
        assert self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp(
            "1", "6s-6s", summed_spin_channels=False
        ).cohp[Spin.down][300] * 2 == approx(
            self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp("1", "6s-6s", summed_spin_channels=True).cohp[
                Spin.up
            ][300]
        )
        assert self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp(
            "1", "6s-6s", summed_spin_channels=True
        ).energies[0] == approx(-15.03759 + 1.96204)
        assert self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp(
            "1", "6s-6s", summed_spin_channels=True
        ).energies[5] == approx(-14.78697 + 1.96204)
        assert not self.cohp_lobster_spin_polarized.get_orbital_resolved_cohp(
            "1", "6s-6s", summed_spin_channels=True
        ).are_coops


class TestMethod:
    def setup_method(self):
        filepath = f"{TEST_DIR}/COHPCAR.lobster.gz"
        structure = f"{TEST_DIR}/POSCAR"
        self.cohp_lobster = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)

        filepath = f"{TEST_DIR}/COHPCAR.lobster.orbitalwise.gz"
        structure = f"{TEST_DIR}/POSCAR.orbitalwise"
        self.cohp_orb = CompleteCohp.from_file("lobster", filename=filepath, structure_file=structure)

        filepath = f"{TEST_DIR}/environments/COHPCAR.lobster.mp-190.gz"
        structure = f"{TEST_DIR}/environments/POSCAR.mp_190.gz"
        self.cohp_lobster_spin_polarized = CompleteCohp.from_file(
            "lobster", filename=filepath, structure_file=structure
        )

    def test_get_integrated_cohp_in_energy_range_full(self):
        # integration up to Fermi level

        cohp = self.cohp_lobster
        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )
        assert result == approx(-0.10218 - 0.19701)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )
        assert result[Spin.up] == approx(-0.10218)
        assert result[Spin.down] == approx(-0.19701)

        # One without spin polarization

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=None,
            relative_E_Fermi=False,
            summed_spin_channels=False,
        )
        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )

        assert result == approx(-4.36062)
        # something else for orbital resolved version
        # self.cohp_lobster_spin_polarized

        result = get_integrated_cohp_in_energy_range(
            self.cohp_lobster_spin_polarized,
            label="1",
            orbital="6s-6s",
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-0.00006)
        assert result[Spin.down] == approx(-0.00006)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_lobster_spin_polarized,
            label="1",
            orbital="6s-6s",
            energy_range=None,
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )

        assert result == approx(-0.00006 * 2)

    def test_get_integrated_cohp_in_energy_range_onefloat(self):
        # only one float is given for energy range
        cohp = self.cohp_lobster
        fermi = cohp.efermi
        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=-0.60201,
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )

        assert result == approx(-0.10218 - 0.19701 + 0.14894 + 0.21889)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=-0.60201,
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-0.10218 + 0.14894)
        assert result[Spin.down] == approx(-0.19701 + 0.21889)
        # only one float is given for energy range (relative to E-fermi)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=-0.60201 + fermi,
            relative_E_Fermi=False,
            summed_spin_channels=True,
        )
        assert result == approx(-0.10218 - 0.19701 + 0.14894 + 0.21889)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=-0.60201 + fermi,
            relative_E_Fermi=False,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-0.10218 + 0.14894)
        assert result[Spin.down] == approx(-0.19701 + 0.21889)

        # without spin
        fermi = self.cohp_orb.efermi
        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=-14.0350 + fermi,
            relative_E_Fermi=False,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=-14.03509,
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )
        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=-14.03509,
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )

        assert result == approx(-4.36062)

    def test_get_integrated_cohp_in_energy_range_whole_range(self):
        cohp = self.cohp_lobster
        fermi = cohp.efermi
        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=[-0.60201, 0],
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )
        assert result == approx(-0.10218 - 0.19701 + 0.14894 + 0.21889)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=[-0.60201, 0],
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-0.10218 + 0.14894)
        assert result[Spin.down] == approx(-0.19701 + 0.21889)
        # whole energy range

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=[-0.60201 + fermi, 0 + fermi],
            relative_E_Fermi=False,
            summed_spin_channels=True,
        )
        assert result == approx(-0.10218 - 0.19701 + 0.14894 + 0.21889)

        result = get_integrated_cohp_in_energy_range(
            cohp,
            label="1",
            orbital=None,
            energy_range=[-0.60201 + fermi, 0 + fermi],
            relative_E_Fermi=False,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-0.10218 + 0.14894)
        assert result[Spin.down] == approx(-0.19701 + 0.21889)

        # without spin
        fermi = self.cohp_orb.efermi
        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=[-14.0350 + fermi, fermi],
            relative_E_Fermi=False,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=[-14.0350, 0],
            relative_E_Fermi=True,
            summed_spin_channels=False,
        )

        assert result[Spin.up] == approx(-4.36062)

        result = get_integrated_cohp_in_energy_range(
            self.cohp_orb,
            label="1",
            orbital=None,
            energy_range=[-14.0350, 0],
            relative_E_Fermi=True,
            summed_spin_channels=True,
        )

        assert result == approx(-4.36062)
