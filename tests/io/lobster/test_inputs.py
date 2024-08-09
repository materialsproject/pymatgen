from __future__ import annotations

import numpy as np
import pytest
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.io.lobster import Lobsterin
from pymatgen.io.lobster.inputs import get_all_possible_basis_combinations
from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp"

__author__ = "Janine George, Marco Esters"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "0.2"
__email__ = "janine.george@uclouvain.be, esters@uoregon.edu"
__date__ = "Dec 10, 2017"


class TestLobsterin(PymatgenTest):
    def setUp(self):
        self.Lobsterin = Lobsterin.from_file(f"{TEST_DIR}/lobsterin.1")
        self.Lobsterin2 = Lobsterin.from_file(f"{TEST_DIR}/lobsterin.2")
        self.Lobsterin3 = Lobsterin.from_file(f"{TEST_DIR}/lobsterin.3")
        self.Lobsterin4 = Lobsterin.from_file(f"{TEST_DIR}/lobsterin.4.gz")

    def test_from_file(self):
        # Test reading from file
        assert self.Lobsterin["cohpstartenergy"] == approx(-15.0)
        assert self.Lobsterin["cohpendenergy"] == approx(5.0)
        assert self.Lobsterin["basisset"] == "pbeVaspFit2015"
        assert self.Lobsterin["gaussiansmearingwidth"] == approx(0.1)
        assert self.Lobsterin["basisfunctions"][0] == "Fe 3d 4p 4s"
        assert self.Lobsterin["basisfunctions"][1] == "Co 3d 4p 4s"
        assert self.Lobsterin["skipdos"]
        assert self.Lobsterin["skipcohp"]
        assert self.Lobsterin["skipcoop"]
        assert self.Lobsterin["skippopulationanalysis"]
        assert self.Lobsterin["skipgrosspopulation"]

        # Test if comments are correctly removed
        assert self.Lobsterin == self.Lobsterin2

    def test_duplicates_from_file(self):
        with open(f"{TEST_DIR}/lobsterin.1") as file:
            original_file = file.readlines()

        # String and float keywords does not allow duplicates
        float_dup_file = original_file.copy()
        float_dup_file.append("cohpstartenergy -15.0")

        float_tmp_file = self.tmp_path / "tmp_lobster_in_float"
        float_tmp_file.write_text("\n".join(float_dup_file))

        with pytest.raises(ValueError, match="Same keyword cohpstartenergy twice!"):
            _ = Lobsterin.from_file(float_tmp_file)

        # Boolean and list keywords allow duplicates
        bool_dup_file = original_file.copy()
        bool_dup_file.append("skipdos")

        bool_tmp_file = self.tmp_path / "tmp_lobster_in_bool"
        bool_tmp_file.write_text("\n".join(bool_dup_file))

        _ = Lobsterin.from_file(bool_tmp_file)  # no error should be raised

    def test_magic_methods(self):
        """Test __getitem__, __setitem__ and __contains__,
        should be case independent.
        """
        # Test __setitem__
        assert self.Lobsterin["COHPSTARTENERGY"] == approx(-15.0)

        with pytest.raises(KeyError, match="Key hello is currently not available"):
            self.Lobsterin["HELLO"] = True

        # Test __getitem__
        self.Lobsterin["skipCOHP"] = False
        assert self.Lobsterin["skipcohp"] is False
        assert self.Lobsterin.get("skipCOHP") is False

        with pytest.raises(KeyError, match="key='World' is not available"):
            _ = self.Lobsterin["World"]

        # Test __contains__
        assert "COHPSTARTENERGY" in self.Lobsterin
        assert "cohpstartenergy" in self.Lobsterin

        assert "helloworld" not in self.Lobsterin

    def test_initialize_from_dict(self):
        # initialize from dict
        lobsterin = Lobsterin(
            {
                "cohpstartenergy": -15.0,
                "cohpendenergy": 5.0,
                "basisset": "pbeVaspFit2015",
                "gaussiansmearingwidth": 0.1,
                "basisfunctions": ["Fe 3d 4p 4s", "Co 3d 4p 4s"],
                "skipdos": True,
                "skipcohp": True,
                "skipcoop": True,
                "skippopulationanalysis": True,
                "skipgrosspopulation": True,
            }
        )
        assert lobsterin["cohpstartenergy"] == approx(-15.0)
        assert lobsterin["cohpendenergy"] == approx(5.0)
        assert lobsterin["basisset"] == "pbeVaspFit2015"
        assert lobsterin["gaussiansmearingwidth"] == approx(0.1)
        assert lobsterin["basisfunctions"][0] == "Fe 3d 4p 4s"
        assert lobsterin["basisfunctions"][1] == "Co 3d 4p 4s"
        assert {*lobsterin} >= {"skipdos", "skipcohp", "skipcoop", "skippopulationanalysis", "skipgrosspopulation"}
        with pytest.raises(KeyError, match="There are duplicates for the keywords!"):
            lobsterin2 = Lobsterin({"cohpstartenergy": -15.0, "cohpstartEnergy": -20.0})
        lobsterin2 = Lobsterin({"cohpstartenergy": -15.0})
        # can only calculate nbands if basis functions are provided
        with pytest.raises(ValueError, match="No basis functions are provided. The program cannot calculate nbands."):
            lobsterin2._get_nbands(structure=Structure.from_file(f"{VASP_IN_DIR}/POSCAR_Fe3O4"))

    def test_standard_settings(self):
        # test standard settings
        for option in [
            "standard",
            "standard_from_projection",
            "standard_with_fatband",
            "onlyprojection",
            "onlydos",
            "onlycohp",
            "onlycoop",
            "onlycobi",
            "onlycohpcoop",
            "onlycohpcoopcobi",
        ]:
            lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
                f"{VASP_IN_DIR}/POSCAR_Fe3O4",
                f"{VASP_IN_DIR}/INCAR.lobster",
                f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz",
                option=option,
            )
            assert lobsterin1["cohpstartenergy"] == approx(-35.0)
            assert lobsterin1["cohpendenergy"] == approx(5.0)
            assert lobsterin1["basisset"] == "pbeVaspFit2015"
            assert lobsterin1["gaussiansmearingwidth"] == approx(0.1)
            assert lobsterin1["basisfunctions"][0] == "Fe 3d 4p 4s "
            assert lobsterin1["basisfunctions"][1] == "O 2p 2s "

            if option in [
                "standard",
                "standard_with_fatband",
                "onlyprojection",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                assert lobsterin1["saveProjectiontoFile"]
            if option in [
                "standard",
                "standard_with_fatband",
                "onlycohp",
                "onlycoop",
                "onlycohpcoop",
            ]:
                assert lobsterin1["cohpGenerator"] == "from 0.1 to 6.0 orbitalwise"
            if option == "standard":
                assert "skipdos" not in lobsterin1
                assert "skipcohp" not in lobsterin1
                assert "skipcoop" not in lobsterin1
            if option == "standard_with_fatband":
                assert lobsterin1["createFatband"] == ["Fe 3d 4p 4s ", "O 2p 2s "]
                assert "skipdos" not in lobsterin1
                assert "skipcohp" not in lobsterin1
                assert "skipcoop" not in lobsterin1
            if option == "standard_from_projection":
                assert lobsterin1["loadProjectionFromFile"], True
            if option in [
                "onlyprojection",
                "onlycohp",
                "onlycoop",
                "onlycobi",
                "onlycohpcoop",
                "onlycohpcoopcobi",
            ]:
                assert lobsterin1["skipdos"], True
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipMadelungEnergy"], True

            if option == "onlydos":
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
                assert lobsterin1["skipMadelungEnergy"], True
            if option == "onlycohp":
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
            if option == "onlycoop":
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcobi"], True
            if option == "onlyprojection":
                assert lobsterin1["skipdos"], True
            if option == "onlymadelung":
                assert lobsterin1["skipPopulationAnalysis"], True
                assert lobsterin1["skipGrossPopulation"], True
                assert lobsterin1["skipcohp"], True
                assert lobsterin1["skipcoop"], True
                assert lobsterin1["skipcobi"], True
                assert lobsterin1["skipdos"], True
        # test basis functions by dict
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            f"{VASP_IN_DIR}/INCAR.lobster",
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        assert lobsterin_new["basisfunctions"] == ["Fe 3d 4p 4s", "O 2s 2p"]

        # test gaussian smearing
        lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            f"{VASP_IN_DIR}/INCAR.lobster2",
            dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
            option="standard",
        )
        assert "gaussiansmearingwidth" not in lobsterin_new

        # fatband and ISMEAR=-5 does not work together
        with pytest.raises(ValueError, match="ISMEAR has to be 0 for a fatband calculation with Lobster"):
            lobsterin_new = Lobsterin.standard_calculations_from_vasp_files(
                f"{VASP_IN_DIR}/POSCAR_Fe3O4",
                f"{VASP_IN_DIR}/INCAR.lobster2",
                dict_for_basis={"Fe": "3d 4p 4s", "O": "2s 2p"},
                option="standard_with_fatband",
            )

    def test_standard_with_energy_range_from_vasprun(self):
        # test standard_with_energy_range_from_vasprun
        lobsterin_comp = Lobsterin.standard_calculations_from_vasp_files(
            f"{VASP_IN_DIR}/POSCAR_C2",
            f"{VASP_IN_DIR}/INCAR_C2",
            f"{VASP_IN_DIR}/POTCAR_C2.gz",
            f"{VASP_OUT_DIR}/vasprun.C2.xml.gz",
            option="standard_with_energy_range_from_vasprun",
        )
        assert lobsterin_comp["COHPstartEnergy"] == -28.3679
        assert lobsterin_comp["COHPendEnergy"] == 32.8968
        assert lobsterin_comp["COHPSteps"] == 301

    def test_diff(self):
        # test diff
        assert self.Lobsterin.diff(self.Lobsterin2)["Different"] == {}
        assert self.Lobsterin.diff(self.Lobsterin2)["Same"]["cohpstartenergy"] == approx(-15.0)

        # test diff in both directions
        for entry in self.Lobsterin.diff(self.Lobsterin3)["Same"]:
            assert entry in self.Lobsterin3.diff(self.Lobsterin)["Same"]
        for entry in self.Lobsterin3.diff(self.Lobsterin)["Same"]:
            assert entry in self.Lobsterin.diff(self.Lobsterin3)["Same"]
        for entry in self.Lobsterin.diff(self.Lobsterin3)["Different"]:
            assert entry in self.Lobsterin3.diff(self.Lobsterin)["Different"]
        for entry in self.Lobsterin3.diff(self.Lobsterin)["Different"]:
            assert entry in self.Lobsterin.diff(self.Lobsterin3)["Different"]

        assert (
            self.Lobsterin.diff(self.Lobsterin3)["Different"]["skipcohp"]["lobsterin1"]
            == self.Lobsterin3.diff(self.Lobsterin)["Different"]["skipcohp"]["lobsterin2"]
        )

    def test_diff_case_insensitivity(self):
        """Test case-insensitivity of diff method."""
        with open(f"{TEST_DIR}/lobsterin.1", encoding="utf-8") as file:
            lobsterin_content = file.read()

        lobsterin_content.replace("COHPstartEnergy -15.0", "cohpSTARTEnergy -15.0")
        lobsterin_content.replace("skipcohp", "skipCOHP")

        tmp_file_path = self.tmp_path / "tmp_lobster_in"
        tmp_file_path.write_text(lobsterin_content)

        lobsterin_diff_case = Lobsterin.from_file(tmp_file_path)
        assert self.Lobsterin.diff(lobsterin_diff_case)["Different"] == {}

    def test_dict_functionality(self):
        for key in ("COHPstartEnergy", "COHPstartEnergy", "COhPstartenergy"):
            start_energy = self.Lobsterin.get(key)
            assert start_energy == -15.0, f"{start_energy=}, {key=}"

        lobsterin_copy = self.Lobsterin.copy()
        lobsterin_copy.update({"cohpstarteNergy": -10.00})
        assert lobsterin_copy["cohpstartenergy"] == -10.0
        lobsterin_copy.pop("cohpstarteNergy")
        assert "cohpstartenergy" not in lobsterin_copy
        lobsterin_copy.pop("cohpendenergY")
        lobsterin_copy["cohpsteps"] = 100
        assert lobsterin_copy["cohpsteps"] == 100
        len_before = len(lobsterin_copy.items())
        assert len_before == 9, f"{len_before=}"

        lobsterin_copy.popitem()
        len_after = len(lobsterin_copy.items())
        assert len_after == len_before - 1

        # Test case sensitivity of |= operator
        self.Lobsterin |= {"skipCOHP": True}  # Camel case
        assert self.Lobsterin["skipcohp"] is True

        self.Lobsterin |= {"skipcohp": False}  # lower case
        assert self.Lobsterin["skipcohp"] is False

    def test_read_write_lobsterin(self):
        outfile_path = self.tmp_path / "lobsterin_test"
        lobsterin1 = Lobsterin.from_file(f"{TEST_DIR}/lobsterin.1")
        lobsterin1.write_lobsterin(outfile_path)
        lobsterin2 = Lobsterin.from_file(outfile_path)
        assert lobsterin1.diff(lobsterin2)["Different"] == {}

    def test_get_basis(self):
        # get basis functions
        lobsterin1 = Lobsterin({})
        potcar = Potcar.from_file(f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz")
        potcar_names = [name["symbol"] for name in potcar.spec]

        assert lobsterin1.get_basis(
            Structure.from_file(f"{TEST_FILES_DIR}/cif/Fe3O4.cif"),
            potcar_symbols=potcar_names,
        ) == ["Fe 3d 4p 4s ", "O 2p 2s "]
        potcar = Potcar.from_file(f"{TEST_DIR}/POTCAR.GaAs")
        potcar_names = [name["symbol"] for name in potcar.spec]
        assert lobsterin1.get_basis(
            Structure.from_file(f"{TEST_DIR}/POSCAR.GaAs"),
            potcar_symbols=potcar_names,
        ) == ["Ga 3d 4p 4s ", "As 4p 4s "]

    def test_get_all_possible_basis_functions(self):
        potcar = Potcar.from_file(f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz")
        potcar_names = [name["symbol"] for name in potcar.spec]
        result = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(f"{TEST_FILES_DIR}/cif/Fe3O4.cif"),
            potcar_symbols=potcar_names,
        )
        assert result[0] == {"Fe": "3d 4s", "O": "2p 2s"}
        assert result[1] == {"Fe": "3d 4s 4p", "O": "2p 2s"}

        potcar2 = Potcar.from_file(f"{FAKE_POTCAR_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz")
        Potcar_names2 = [name["symbol"] for name in potcar2.spec]
        result2 = Lobsterin.get_all_possible_basis_functions(
            Structure.from_file(f"{TEST_FILES_DIR}/cif/Fe.cif"),
            potcar_symbols=Potcar_names2,
        )
        assert result2[0] == {"Fe": "3d 4s"}

    def test_get_potcar_symbols(self):
        lobsterin1 = Lobsterin({})
        assert lobsterin1._get_potcar_symbols(f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz") == ["Fe", "O"]
        assert lobsterin1._get_potcar_symbols(f"{TEST_DIR}/POTCAR.GaAs") == ["Ga_d", "As"]

    def test_write_lobsterin(self):
        # write lobsterin, read it and compare it
        outfile_path = self.tmp_path / "lobsterin_test"
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            f"{VASP_IN_DIR}/INCAR.lobster",
            f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz",
            option="standard",
        )
        lobsterin1.write_lobsterin(outfile_path)
        lobsterin2 = Lobsterin.from_file(outfile_path)
        assert lobsterin1.diff(lobsterin2)["Different"] == {}

    def test_write_incar(self):
        # write INCAR and compare
        outfile_path = self.tmp_path / "INCAR_test"
        lobsterin1 = Lobsterin.standard_calculations_from_vasp_files(
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            f"{VASP_IN_DIR}/INCAR.lobster",
            f"{VASP_IN_DIR}/POTCAR_Fe3O4.gz",
            option="standard",
        )
        lobsterin1.write_INCAR(
            f"{VASP_IN_DIR}/INCAR.lobster3",
            outfile_path,
            f"{VASP_IN_DIR}/POSCAR_Fe3O4",
            isym=-1,
        )

        incar1 = Incar.from_file(f"{VASP_IN_DIR}/INCAR.lobster3")
        incar2 = Incar.from_file(outfile_path)

        assert incar1.diff(incar2)["Different"] == {
            "ISYM": {"INCAR1": 2, "INCAR2": -1},
            "NBANDS": {"INCAR1": None, "INCAR2": 86},
            "NSW": {"INCAR1": 500, "INCAR2": 0},
            "LWAVE": {"INCAR1": False, "INCAR2": True},
        }

    def test_write_kpoints(self):
        # line mode
        outfile_path = self.tmp_path / "KPOINTS_test"
        outfile_path2 = self.tmp_path / "POSCAR_test"
        lobsterin1 = Lobsterin({})
        # test writing primitive cell
        lobsterin1.write_POSCAR_with_standard_primitive(
            POSCAR_input=f"{VASP_IN_DIR}/POSCAR_Fe3O4", POSCAR_output=outfile_path2
        )

        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            kpoints_line_density=58,
            isym=-1,
        )
        kpoint = Kpoints.from_file(outfile_path)
        assert kpoint.num_kpts == 562
        assert kpoint.kpts[-1][0] == approx(-0.5)
        assert kpoint.kpts[-1][1] == approx(0.5)
        assert kpoint.kpts[-1][2] == approx(0.5)
        assert kpoint.labels[-1] == "T"
        kpoint2 = Kpoints.from_file(f"{VASP_IN_DIR}/KPOINTS_band.lobster")

        labels = []
        number = 0
        for label in kpoint.labels:
            if label is not None:
                if number != 0:
                    if label != labels[number - 1]:
                        labels.append(label)
                        number += 1
                else:
                    labels.append(label)
                    number += 1

        labels2 = []
        number2 = 0
        for label in kpoint2.labels:
            if label is not None:
                if number2 != 0:
                    if label != labels2[number2 - 1]:
                        labels2.append(label)
                        number2 += 1
                else:
                    labels2.append(label)
                    number2 += 1
        assert labels == labels2

        # without line mode
        lobsterin1.write_KPOINTS(POSCAR_input=outfile_path2, KPOINTS_output=outfile_path, line_mode=False, isym=-1)
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{VASP_OUT_DIR}/IBZKPT.lobster")

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            assert list_kpoint[0] == approx(kpoint2.kpts[num_kpt][0])
            assert list_kpoint[1] == approx(kpoint2.kpts[num_kpt][1])
            assert list_kpoint[2] == approx(kpoint2.kpts[num_kpt][2])

        assert kpoint.num_kpts == 108

        # without line mode, use grid instead of reciprocal density
        lobsterin1.write_KPOINTS(
            POSCAR_input=outfile_path2,
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[6, 6, 3],
            isym=-1,
        )
        kpoint = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{VASP_OUT_DIR}/IBZKPT.lobster")

        for num_kpt, list_kpoint in enumerate(kpoint.kpts):
            assert list_kpoint[0] == approx(kpoint2.kpts[num_kpt][0])
            assert list_kpoint[1] == approx(kpoint2.kpts[num_kpt][1])
            assert list_kpoint[2] == approx(kpoint2.kpts[num_kpt][2])

        assert kpoint.num_kpts == 108

        #
        # #without line mode, using a certain grid, isym=0 instead of -1
        lobsterin1.write_KPOINTS(
            POSCAR_input=f"{TEST_DIR}/POSCAR.Li",
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[3, 3, 3],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_DIR}/IBZKPT_3_3_3_Li")
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint2.kpts,
                kpoint1.kpts_weights[ikpoint],
                kpoint2.kpts_weights,
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint1.kpts,
                kpoint2.kpts_weights[ikpoint],
                kpoint1.kpts_weights,
            )

        lobsterin1.write_KPOINTS(
            POSCAR_input=f"{TEST_DIR}/POSCAR.Li",
            KPOINTS_output=outfile_path,
            line_mode=False,
            from_grid=True,
            input_grid=[2, 2, 2],
            isym=0,
        )

        kpoint1 = Kpoints.from_file(outfile_path)
        kpoint2 = Kpoints.from_file(f"{TEST_DIR}/IBZKPT_2_2_2_Li")
        for ikpoint, kpoint in enumerate(kpoint1.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint2.kpts,
                kpoint1.kpts_weights[ikpoint],
                kpoint2.kpts_weights,
            )
        for ikpoint, kpoint in enumerate(kpoint2.kpts):
            assert self.is_kpoint_in_list(
                kpoint,
                kpoint1.kpts,
                kpoint2.kpts_weights[ikpoint],
                kpoint1.kpts_weights,
            )

    def is_kpoint_in_list(self, kpoint, kpointlist, weight, weightlist) -> bool:
        found = 0
        for ikpoint2, kpoint2 in enumerate(kpointlist):
            if (
                np.isclose(kpoint[0], kpoint2[0])
                and np.isclose(kpoint[1], kpoint2[1])
                and np.isclose(kpoint[2], kpoint2[2])
            ):
                if weight == weightlist[ikpoint2]:
                    found += 1
            elif (
                np.isclose(-kpoint[0], kpoint2[0])
                and np.isclose(-kpoint[1], kpoint2[1])
                and np.isclose(-kpoint[2], kpoint2[2])
            ) and weight == weightlist[ikpoint2]:
                found += 1
        return found == 1

    def test_as_from_dict(self):
        # tests as dict and from dict methods
        new_lobsterin = Lobsterin.from_dict(self.Lobsterin.as_dict())
        assert new_lobsterin == self.Lobsterin
        new_lobsterin.to_json()


class TestUtils(PymatgenTest):
    def test_get_all_possible_basis_combinations(self):
        # this basis is just for testing (not correct)
        min_basis = ["Li 1s 2s ", "Na 1s 2s", "Si 1s 2s"]
        max_basis = ["Li 1s 2p 2s ", "Na 1s 2p 2s", "Si 1s 2s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Li 1s 2s", "Na 1s 2s", "Si 1s 2s"],
            ["Li 1s 2s", "Na 1s 2s 2p", "Si 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s", "Si 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p", "Si 1s 2s"],
        ]

        min_basis = ["Li 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [["Li 1s 2s"], ["Li 1s 2s 2p"], ["Li 1s 2s 3s"], ["Li 1s 2s 2p 3s"]]

        min_basis = ["Li 1s 2s", "Na 1s 2s"]
        max_basis = ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Li 1s 2s", "Na 1s 2s"],
            ["Li 1s 2s", "Na 1s 2s 2p"],
            ["Li 1s 2s", "Na 1s 2s 3s"],
            ["Li 1s 2s", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 2p", "Na 1s 2s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p"],
            ["Li 1s 2s 2p", "Na 1s 2s 3s"],
            ["Li 1s 2s 2p", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 3s", "Na 1s 2s"],
            ["Li 1s 2s 3s", "Na 1s 2s 2p"],
            ["Li 1s 2s 3s", "Na 1s 2s 3s"],
            ["Li 1s 2s 3s", "Na 1s 2s 2p 3s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 2p"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 3s"],
            ["Li 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
        ]

        min_basis = ["Si 1s 2s 2p", "Na 1s 2s"]
        max_basis = ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"]
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        assert combinations_basis == [
            ["Si 1s 2s 2p", "Na 1s 2s"],
            ["Si 1s 2s 2p", "Na 1s 2s 2p"],
            ["Si 1s 2s 2p", "Na 1s 2s 3s"],
            ["Si 1s 2s 2p", "Na 1s 2s 2p 3s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 2p"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 3s"],
            ["Si 1s 2s 2p 3s", "Na 1s 2s 2p 3s"],
        ]
