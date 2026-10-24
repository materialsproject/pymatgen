from __future__ import annotations

import warnings

import matplotlib.pyplot as plt
import pytest
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.lobster_env import LobsterLightStructureEnvironments, LobsterNeighbors, check_ICOHPs
from pymatgen.core import Element, PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import Cohp, CompleteCohp, IcohpValue
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster import Charge, Icohplist
from tests.testing import TEST_FILES_DIR

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__email__ = "janine.george@uclouvain.be"
__date__ = "Jan 14, 2021"

TEST_DIR = f"{TEST_FILES_DIR}/analysis/lobster_env"


class TestLobsterNeighbors:
    def setup_method(self):
        # test additional conditions first
        # only consider cation anion bonds

        self.chem_env_lobster1 = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            additional_condition=1,
            perc_strength_icohp=0.3,
            noise_cutoff=0.0,
        )

        # all bonds
        self.chem_env_lobster0 = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=0,
        )
        # test __init_new__ using from_files
        self.chem_env_lobster0_from_file = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=0,
        )

        # only cation-cation, anion-anion bonds
        self.chem_env_lobster5 = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=5,
        )

        self.chem_env_lobster5_from_file = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=5,
        )

        # only cation-cation bonds
        self.chem_env_lobster6 = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            additional_condition=6,
        )

        self.chem_env_lobster6_from_file = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=6,
        )

        # 2,3,4 are not tested so far
        self.chem_env_lobster2 = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=2,
        )

        self.chem_env_lobster3 = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=3,
        )

        self.chem_env_lobster4 = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            additional_condition=4,
        )

        # search for other testcase where 2,3,4 arrive at different results
        self.chem_env_lobster0_second = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=0,
            perc_strength_icohp=0.05,
        )
        self.chem_env_lobster1_second = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=1,
        )

        self.chem_env_lobster2_second = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=2,
        )

        self.chem_env_lobster5_second = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=5,
            perc_strength_icohp=0.05,
        )

        self.chem_env_lobster5_second_percentage = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=5,
            perc_strength_icohp=1.0,
        )

        self.chem_env_lobster6_second = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            additional_condition=6,
            perc_strength_icohp=0.05,
        )
        # coop / cobi
        self.chem_env_lobster1_coop_NaCl = LobsterNeighbors.from_files(
            are_coops=True,
            structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_coop_NaCl_from_file = LobsterNeighbors.from_files(
            are_coops=True,
            structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_cobi_NaCl = LobsterNeighbors.from_files(
            are_coops=True,
            structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_cobi_NaCl_from_file = LobsterNeighbors.from_files(
            are_coops=True,
            structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_cobi_mp470 = LobsterNeighbors.from_files(
            are_coops=True,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-470.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOBILIST.lobster.mp-470.gz",
            additional_condition=1,
        )

        # TODO: use charge instead of valence
        self.chem_env_lobster1_charges = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
        )
        self.chem_env_lobster1_charges_noisecutoff = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-632319.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-632319.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-632319.gz",
            additional_condition=1,
            perc_strength_icohp=0.05,
            noise_cutoff=0.1,
        )
        self.chem_env_lobster1_charges_wo_noisecutoff = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-632319.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-632319.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-632319.gz",
            additional_condition=1,
            perc_strength_icohp=0.05,
            noise_cutoff=None,
        )
        self.chem_env_lobster1_charges_loewdin = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            which_charge="Loewdin",
        )
        self.chem_env_lobster6_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=6,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster5_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=5,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster4_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=4,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster3_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=3,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster2_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=2,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster1_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            adapt_extremum_to_add_cond=True,
        )

        self.chem_env_lobster0_charges_additional_condition = LobsterNeighbors.from_files(
            are_coops=False,
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster0_NaSi = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.NaSi.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.NaSi.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster_NaSi_wo_charges = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.NaSi.gz",
            valences_from_charges=False,
            charge_path=None,
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        # Test LobsterNeighbors using pymatgen objects
        self.obj_icohp = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz")
        self.obj_charge = Charge(filename=f"{TEST_DIR}/CHARGE.lobster.NaSi.gz")
        self.chem_env_w_obj = LobsterNeighbors(
            are_coops=False,
            icoxxlist_obj=self.obj_icohp,
            charge_obj=self.obj_charge,
            structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaSi.gz"),
            valences_from_charges=True,
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        # test LSE on_error
        self.chem_env_on_error = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-1018096.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-1018096.gz",
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-1018096.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
            valences_from_charges=True,
        )

    def test_cation_anion_mode_without_ions(self):
        with pytest.raises(
            ValueError,
            match="Valences cannot be assigned, additional_conditions 1, 3, 5 and 6 will not work",
        ):
            _ = LobsterNeighbors.from_files(
                are_coops=False,
                icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.gz",
                structure_path=f"{TEST_DIR}/POSCAR.gz",
                valences_from_charges=False,
                additional_condition=1,
            )
        with pytest.raises(
            ValueError,
            match="All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work",
        ):
            _ = LobsterNeighbors.from_files(
                are_coops=False,
                icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.gz",
                structure_path=f"{TEST_DIR}/POSCAR.gz",
                valences_from_charges=False,
                additional_condition=1,
                valences=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            )

    def test_wrong_additional_correction(self):
        with pytest.raises(
            ValueError,
            match=r"Unexpected additional_condition=10, must be one of \[0, 1, 2, 3, 4, 5, 6\]",
        ):
            LobsterNeighbors.from_files(
                are_coops=False,
                icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
                structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
                valences_from_charges=True,
                charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
                additional_condition=10,
            )

    def test_set_limits(self):
        test = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            limits=[-100000, 0],
        )
        assert test.limits == [-100000, 0]

    def test_molecules_allowed(self):
        assert not self.chem_env_lobster1.molecules_allowed

    @pytest.mark.filterwarnings("ignore:get_anion_types is deprecated")
    def test_get_anion_types(self):
        assert self.chem_env_lobster0_second.get_anion_types() == {Element("O")}
        assert self.chem_env_lobster0_second.anion_types == {Element("O")}

    def test_get_nn_info(self):
        # 0: NO_ADDITIONAL_CONDITION
        # 1: ONLY_ANION_CATION_BONDS
        # 2: NO_ELEMENT_TO_SAME_ELEMENT_BONDS
        # 3: ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
        # 4: ONLY_ELEMENT_TO_OXYGEN_BONDS
        # 5: DO_NOT_CONSIDER_ANION_CATION_BONDS
        # 6: ONLY_CATION_CATION_BONDS

        # All bonds
        # ReO3
        assert (
            len(
                self.chem_env_lobster0.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster0.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # 1: ONLY_ANION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster1.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster1.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 2
        )
        assert (
            len(
                self.chem_env_lobster1_charges_noisecutoff.get_nn(
                    structure=self.chem_env_lobster1_charges_noisecutoff.structure,
                    n=1,
                )
            )
            == 0
        )
        assert (
            len(
                self.chem_env_lobster1_charges_wo_noisecutoff.get_nn(
                    structure=self.chem_env_lobster1_charges_wo_noisecutoff.structure,
                    n=1,
                )
            )
            == 8
        )
        # 2: NO_ELEMENT_TO_SAME_ELEMENT_BONDS
        assert (
            len(
                self.chem_env_lobster2.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster2.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # 3: ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS
        assert (
            len(
                self.chem_env_lobster3.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster3.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # 4: ONLY_ELEMENT_TO_OXYGEN_BONDS
        assert (
            len(
                self.chem_env_lobster4.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster4.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # 5: DO_NOT_CONSIDER_ANION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster5.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 0
        )
        assert (
            len(
                self.chem_env_lobster5.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 0
        )
        # 6: ONLY_CATION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster6.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=0,
                )
            )
            == 0
        )

        assert (
            len(
                self.chem_env_lobster6.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-190.gz"),
                    n=1,
                )
            )
            == 0
        )

        # All bonds
        # mp-353, Ag2O
        # all bonds
        assert (
            len(
                self.chem_env_lobster0_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 8
        )

        # 1: ONLY_ANION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster1_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 2
        )

        assert (
            len(
                self.chem_env_lobster1_coop_NaCl.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chem_env_lobster1_cobi_NaCl.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chem_env_lobster1_cobi_mp470.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-470.gz"),
                    n=3,
                )
            )
            == 3
        )

        # 2: NO_ELEMENT_TO_SAME_ELEMENT_BONDS
        assert (
            len(
                self.chem_env_lobster2_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 2
        )
        assert (
            len(
                self.chem_env_lobster2_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=4,
                )
            )
            == 4
        )

        # 5: DO_NOT_CONSIDER_ANION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster5_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster5_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=4,
                )
            )
            == 4
        )

        # 6: ONLY_CATION_CATION_BONDS
        assert (
            len(
                self.chem_env_lobster6_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster6_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=4,
                )
            )
            == 0
        )

        assert (
            len(
                self.chem_env_lobster5_second_percentage.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"),
                    n=0,
                )
            )
            == 0
        )

        # NaSi
        # all bonds
        # passes after changes introduced in PR #4148 (correct ordering of all edge attributes)
        nn0 = self.chem_env_lobster0_NaSi.get_nn_info(structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaSi.gz"), n=0)
        nn11 = self.chem_env_lobster0_NaSi.get_nn_info(
            structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaSi.gz"), n=11
        )
        assert next(d["edge_properties"] for d in nn0 if d.get("site_index") == 11) == next(
            d["edge_properties"] for d in nn11 if d.get("site_index") == 0
        )

    def test_structure_graph(self):
        sg = self.chem_env_lobster1_second.get_bonded_structure(
            structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz")
        )
        assert isinstance(sg, StructureGraph)

    def test_extended_structure_graph(self):
        self.chem_env_lobsterNaCl = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.NaCl.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
            valences_from_charges=True,
            charge_path=f"{TEST_DIR}/CHARGE.lobster.NaCl.gz",
            blist_sg1_path=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
            blist_sg2_path=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
            add_additional_data_sg=True,
            id_blist_sg1="icobi",
            id_blist_sg2="icoop",
            additional_condition=1,
        )
        sg = self.chem_env_lobsterNaCl.get_bonded_structure(
            structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.NaCl.gz"),
            decorate=True,
            edge_properties=True,
            weights=True,
        )
        assert sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == approx(-0.59138)
        assert sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == approx(0.08682)
        assert sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == approx(0.03007)
        assert sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"
        assert sg.graph.get_edge_data(0, 1)[5]["bond_label"] == "30"
        assert isinstance(sg, StructureGraph)

    def test_raises_extended_structure_graph(self):
        with pytest.raises(ValueError, match="Algorithm can only work with ICOOPs, ICOBIs"):
            self.chem_env_lobsterNaCl = LobsterNeighbors.from_files(
                are_coops=False,
                icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.NaCl.gz",
                structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
                valences_from_charges=True,
                charge_path=f"{TEST_DIR}/CHARGE.lobster.NaCl.gz",
                blist_sg1_path=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
                blist_sg2_path=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
                add_additional_data_sg=True,
                id_blist_sg1="icopppp",
                id_blist_sg2="icoop",
                additional_condition=1,
            )

    def test_order_parameter(self):
        assert self.chem_env_lobster1_second.get_local_order_parameters(
            structure=Structure.from_file(f"{TEST_DIR}/CONTCAR.mp-353.gz"), n=0
        )["linear"] == approx(1.0)

    def test_get_structure_environments(self):
        lse = self.chem_env_lobster1_second.get_light_structure_environment()
        assert lse.coordination_environments[0][0]["ce_symbol"] == "L:2"
        assert lse.coordination_environments[5][0]["ce_symbol"] == "T:4"

        lse2 = self.chem_env_lobster1.get_light_structure_environment()
        assert lse2.coordination_environments[0][0]["ce_symbol"] == "O:6"

    def test_get_structure_environments_further_tests(self):
        lse = self.chem_env_lobster1_second.get_light_structure_environment()
        lse.as_dict()
        lse.get_statistics()
        assert lse.uniquely_determines_coordination_environments

    def test_get_info_icohps_neighbors(self):
        # Bond labels in results[3] are a set of identifiers — their iteration order is
        # not part of the API contract and shifts between platforms / dependency versions
        # (e.g. ASE neighbor list ordering). Compare as a sorted list.
        expected_bond_labels = sorted(["48", "27", "64", "73", "49", "30"])

        results = self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])
        assert results[0] == approx(-33.87452)
        for bond in results[1]:
            assert bond == approx(-5.64612, abs=1e-2)
        assert results[2] == 6
        assert sorted(results[3]) == expected_bond_labels

        results2 = self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=None)
        assert results2[0] == approx(-33.87452)
        for bond in results2[1]:
            assert bond == approx(-5.64455, abs=1e-2)
        assert results2[2] == 6

        assert sorted(results2[3]) == expected_bond_labels
        from collections import Counter

        assert Counter(map(tuple, results2[4])) == Counter(
            [("Re1", "O2"), ("Re1", "O2"), ("Re1", "O3"), ("Re1", "O3"), ("Re1", "O4"), ("Re1", "O4")]
        )

    def test_get_sum_icohps_between_neighbors_of_atom(self):
        # will only look at icohps between cations or anions
        self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=[1])
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[1])[2] == 1
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[1])[0] == approx(-0.04535)

        # confirmed by looking at the VESTA output
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[0])[2] == 12
        # use an example where this is easier to test (e.g., linear environment?)

        chemenv_here = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-7000.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-7000.gz",
            additional_condition=1,
            perc_strength_icohp=0.05,
            noise_cutoff=0.0,
        )

        assert len(chemenv_here.get_info_icohps_between_neighbors(isites=[0])[4]) == 2

    def test_get_plot_label(self):
        label = self.chem_env_lobster1._get_plot_label(
            atoms=[
                ["Re1", "O2"],
                ["Re1", "O2"],
                ["Re1", "O3"],
                ["Re1", "O3"],
                ["Re1", "O4"],
                ["Re1", "O4"],
            ],
            per_bond=False,
        )
        assert label == "6 x O-Re"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[
                ["Re1", "O2"],
                ["Re1", "O2"],
                ["Re1", "O3"],
                ["Re1", "O3"],
                ["Re1", "O4"],
                ["Si1", "O4"],
            ],
            per_bond=False,
        )
        assert label == "5 x O-Re, 1 x O-Si"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[
                ["Si1", "O2"],
                ["Si1", "O2"],
                ["Si1", "O3"],
                ["Re1", "O3"],
                ["Re1", "O4"],
                ["Si1", "O4"],
            ],
            per_bond=False,
        )
        assert label == "4 x O-Si, 2 x O-Re"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[
                ["Re1", "O2"],
                ["Re1", "O2"],
                ["Re1", "O3"],
                ["Re1", "O3"],
                ["Re1", "O4"],
                ["Re1", "O4"],
            ],
            per_bond=True,
        )
        assert label == "6 x O-Re (per bond)"

    def test_get_info_cohps_to_neighbors(self):
        chem_env_lobster1 = LobsterNeighbors.from_files(
            are_coops=False,
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
            structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
            additional_condition=1,
        )
        cohpcar_lobster_mp_190 = f"{TEST_DIR}/COHPCAR.lobster.mp-190.gz"
        plot_label, summed_cohpcar_mp_190 = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
        )
        assert plot_label == "6 x O-Re (per bond)"
        assert isinstance(summed_cohpcar_mp_190, Cohp)

        cohp_result = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
        )[1]
        print(cohp_result.icohp[Spin.up][700])
        print(chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])[0])
        assert cohp_result.icohp[Spin.up][700] == approx(chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])[0])

        # summed_spin_channel
        coph_thing = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
            summed_spin_channels=True,
        )[1]
        assert coph_thing.icohp[Spin.up][700] == approx(chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])[0])

        plot_label, summed_cohpcar_mp_190_Te = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["Te"],
        )

        assert plot_label is None
        assert summed_cohpcar_mp_190_Te is None

        plot_label, _summed_cohpcar_NaSi = self.chem_env_lobster0_NaSi.get_info_cohps_to_neighbors(
            path_to_cohpcar=f"{TEST_DIR}/COHPCAR.lobster.NaSi.gz",
            isites=[8],
            onlycation_isites=False,
            only_bonds_to=["Na"],
        )
        assert plot_label == "1 x Na-Si (per bond)"

        obj_cohpcar = CompleteCohp.from_file(
            filename=f"{TEST_DIR}/COHPCAR.lobster.NaSi.gz",
            fmt="LOBSTER",
            structure_file=f"{TEST_DIR}/CONTCAR.NaSi.gz",
        )
        plot_label_obj, _summed_cohpcar_NaSi_obj = self.chem_env_w_obj.get_info_cohps_to_neighbors(
            coxxcar_obj=obj_cohpcar,
            isites=[8],
            onlycation_isites=False,
            only_bonds_to=["Na"],
        )
        assert plot_label_obj == "1 x Na-Si (per bond)"

        info = self.chem_env_lobster0_NaSi.get_info_cohps_to_neighbors(
            path_to_cohpcar=f"{TEST_DIR}/COHPCAR.lobster.NaSi.gz",
            isites=[8],
            onlycation_isites=False,
            only_bonds_to=["Si"],
        )[0]
        assert info == "3 x Si-Si (per bond)"

        chem_env_lobster1.plot_cohps_of_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
        )

        chem_env_lobster1.plot_cohps_of_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
            xlim=[-10, 10],
            ylim=None,
        )

        expected_msg = "COHPCAR and ICOHPLIST do not fit together"
        with pytest.raises(ValueError, match=expected_msg):
            # icohplist and cohpcar do not fit together
            self.chem_env_lobster1.get_info_cohps_to_neighbors(
                path_to_cohpcar=f"{TEST_DIR}/COHPCAR.lobster.mp-190_2.gz",
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

        with pytest.raises(ValueError, match=expected_msg):
            # icohplist and cohpcar do not fit together
            self.chem_env_lobster2.get_info_cohps_to_neighbors(
                path_to_cohpcar=f"{TEST_DIR}/COHPCAR.lobster.mp-190_2.gz",
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

    def test_valences(self):
        assert self.chem_env_lobster1_charges_noisecutoff.valences == [
            0.75,
            -0.75,
        ]  # Mulliken
        assert self.chem_env_lobster1_charges_loewdin.valences == [
            0.51,
            0.51,
            0.51,
            0.51,
            -1.01,
            -1.01,
        ]

        assert_allclose(
            self.chem_env_w_obj.valences, [0.66] * 4 + [0.69] * 4 + [-0.68] * 3 + [-0.69] + [-0.67] * 4
        )  # charge_obj
        assert self.chem_env_lobster_NaSi_wo_charges.valences == [1] * 8 + [-1] * 8  # BVA

    def test_backward_compatibility_warning(self):
        with pytest.warns(UserWarning, match="You are using an older version for neighbor detection"):
            _ = LobsterNeighbors.from_files(
                icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
                structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
                backward_compatibility=True,
            )

    def test_on_error_case(self):
        with pytest.raises(
            ValueError,
            match=r"Environment cannot be determined for site 0\. Number of neighbors \(18\) is larger than 13\.",
        ):
            _ = self.chem_env_on_error.get_light_structure_environment(on_error="raise")

        with pytest.warns(
            UserWarning,
            match=r"Site 0 has 18 neighbors \(>13\)\. Using coordination number instead of geometry\.",
        ) as warnings_record:
            lse = self.chem_env_on_error.get_light_structure_environment(on_error="warn")

        assert lse.coordination_environments[0][0]["ce_symbol"] == "18"
        assert [str(warn.message) for warn in warnings_record] == [
            "Site 0 has 18 neighbors (>13). Using coordination number instead of geometry.",
            "Site 1 has 18 neighbors (>13). Using coordination number instead of geometry.",
        ]

        # if on_error is "ignore", no error should be raised and no warning should be issued
        with warnings.catch_warnings(record=True) as warnings_record:
            lse = self.chem_env_on_error.get_light_structure_environment(on_error="ignore")

        assert lse.coordination_environments[0][0]["ce_symbol"] == "18"
        assert warnings_record == []


# The tests below use module-level helpers and fixtures (instead of the expensive
# `setup_method` above) to cover remaining options and error branches.


def _neighbors_mp190(**kwargs) -> LobsterNeighbors:
    return LobsterNeighbors.from_files(
        icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-190.gz",
        structure_path=f"{TEST_DIR}/CONTCAR.mp-190.gz",
        **kwargs,
    )


def _neighbors_without_valences(additional_condition: int = 0) -> LobsterNeighbors:
    # Fe16 cell for which the bond valence analysis cannot assign valences
    with pytest.warns(UserWarning, match="Please consider using a newer LOBSTER version"):
        return LobsterNeighbors.from_files(
            icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.gz",
            structure_path=f"{TEST_DIR}/POSCAR.gz",
            additional_condition=additional_condition,
        )


@pytest.fixture(scope="module")
def ln_mp190() -> LobsterNeighbors:
    """ReO3 (mp-190), only cation-anion bonds. Only used for read-only queries."""
    return _neighbors_mp190(additional_condition=1)


@pytest.fixture(scope="module")
def ln_mp353() -> LobsterNeighbors:
    """Ag2O (mp-353) with Mulliken charges as valences. Only used for read-only queries."""
    return LobsterNeighbors.from_files(
        icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.mp-353.gz",
        structure_path=f"{TEST_DIR}/CONTCAR.mp-353.gz",
        charge_path=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
        valences_from_charges=True,
    )


@pytest.mark.parametrize(
    ("additional_condition", "expected_n_neighbors"),
    [
        (0, [0] * 6 + [1, 2, 1] + [0] * 7),
        # Fe-Fe bonds only: conditions 2 and 4 exclude all of them
        (2, [0] * 16),
        (4, [0] * 16),
    ],
)
def test_conditions_without_valences(additional_condition, expected_n_neighbors):
    """Conditions that do not need valences must work when BVAnalyzer fails."""
    lobster_nb = _neighbors_without_valences(additional_condition)
    assert lobster_nb.valences is None
    assert [len(neighbors) for neighbors in lobster_nb.list_neighisite] == expected_n_neighbors


def test_no_valences_errors_and_queries():
    lobster_nb = _neighbors_without_valences()

    with pytest.raises(ValueError, match="No cations and anions defined"):
        _ = lobster_nb.anion_types
    with pytest.raises(ValueError, match="No valences are provided"):
        lobster_nb.get_info_icohps_to_neighbors()
    with pytest.raises(ValueError, match="No valences are provided"):
        lobster_nb.get_info_icohps_between_neighbors()
    with pytest.raises(ValueError, match=r"self\.valences=None"):
        lobster_nb.get_light_structure_environment(only_cation_environments=True)

    # all sites are used if onlycation_isites=False
    info = lobster_nb.get_info_icohps_to_neighbors(onlycation_isites=False)
    assert info.total_icohp == approx(-2.33366)
    assert info.n_bonds == 4
    assert sorted(info.central_isites) == [6, 7, 7, 8]
    assert sum(info.list_icohps) == approx(info.total_icohp)

    info_between = lobster_nb.get_info_icohps_between_neighbors(onlycation_isites=False)
    assert info_between.central_isites is None
    assert info_between.n_bonds == len(info_between.labels) == len(info_between.atoms)


def test_get_info_icohps_between_neighbors_all_sites():
    lobster_nb = _neighbors_mp190(additional_condition=1)
    only_cations = lobster_nb.get_info_icohps_between_neighbors()
    all_sites = lobster_nb.get_info_icohps_between_neighbors(onlycation_isites=False)

    # 12 O-O interactions in the ReO6 octahedron, plus 3 Re-Re interactions across the O sites
    assert only_cations.n_bonds == 12
    assert only_cations.total_icohp == approx(-1.59574)
    assert all_sites.n_bonds == 15
    assert all_sites.total_icohp == approx(-1.73182)
    assert all_sites.list_icohps[-3:] == approx([-0.04535, -0.04535, -0.04538])
    assert all(len(atoms) == 2 for atoms in all_sites.atoms)


def test_get_nn_info_errors(ln_mp190):
    with pytest.raises(ValueError, match="LobsterEnv cannot use weights"):
        ln_mp190.get_nn_info(ln_mp190.structure, n=0, use_weights=True)

    supercell = ln_mp190.structure * (2, 1, 1)
    with pytest.raises(ValueError, match=r"Length of structure \(8\) and LobsterNeighbors \(4\) differ"):
        ln_mp190.get_nn_info(supercell, n=0)


def test_only_one_limit_given():
    with pytest.raises(ValueError, match="Please give two limits or leave them both at None"):
        _neighbors_mp190(limits=(-100.0, None))


def test_from_files_default_bonding_list_ids():
    """id_blist_sg1="icoop" and id_blist_sg2="icobi" (the defaults)."""
    lobster_nb = LobsterNeighbors.from_files(
        icoxxlist_path=f"{TEST_DIR}/ICOHPLIST.lobster.NaCl.gz",
        structure_path=f"{TEST_DIR}/CONTCAR.NaCl.gz",
        valences_from_charges=True,
        charge_path=f"{TEST_DIR}/CHARGE.lobster.NaCl.gz",
        blist_sg1_path=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
        blist_sg2_path=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
        add_additional_data_sg=True,
        additional_condition=1,
    )
    assert lobster_nb.bonding_list_1.are_coops
    assert lobster_nb.bonding_list_2.are_cobis
    edge = next(nb for nb in lobster_nb.get_nn_info(lobster_nb.structure, n=0) if nb["site_index"] == 1)
    assert edge["edge_properties"]["ICOOP"] == approx(0.03007)
    assert edge["edge_properties"]["ICOBI"] == approx(0.08682)


def test_light_structure_environment_without_neighbors():
    # ReO3 has no cation-cation bonds, hence no neighbors for any site
    lobster_nb = _neighbors_mp190(additional_condition=6)
    assert lobster_nb.valences == [6, -2, -2, -2]
    assert all(neighbors == [] for neighbors in lobster_nb.list_neighisite)

    lse = lobster_nb.get_light_structure_environment()
    for ce in lse.coordination_environments:
        assert ce[0]["ce_symbol"] is None
        assert ce[0]["csm"] is None
        assert ce[0]["permutation"] is None


def test_light_structure_environment_only_cations(ln_mp190):
    lse = ln_mp190.get_light_structure_environment(only_cation_environments=True)
    assert [ce[0]["ce_symbol"] for ce in lse.coordination_environments] == ["O:6", None, None, None]
    assert [len(nb_sets[0]) for nb_sets in lse.neighbors_sets] == [6, 0, 0, 0]
    assert lse.coordination_environments[1][0]["permutation"] == []
    assert lse.valences == [6, -2, -2, -2]


def test_light_structure_environment_only_indices(ln_mp190):
    lse = ln_mp190.get_light_structure_environment(only_indices=[1])
    assert [ce[0]["ce_symbol"] for ce in lse.coordination_environments] == [None, "L:2", None, None]
    assert [len(nb_sets[0]) for nb_sets in lse.neighbors_sets] == [0, 2, 0, 0]

    dct = lse.as_dict()
    assert dct["@class"] == "LobsterLightStructureEnvironments"
    assert dct["valences"] == [6, -2, -2, -2]
    # only the two Re neighbors of O2 are stored
    assert [nb["index"] for nb in dct["all_nbs_sites"]] == [0, 0]


def test_from_lobster_without_ce_and_neighbors(ln_mp190):
    structure = ln_mp190.structure
    n_sites = len(structure)
    lse = LobsterLightStructureEnvironments.from_Lobster(
        list_ce_symbol=None,
        list_csm=[None] * n_sites,
        list_permutation=[None] * n_sites,
        list_neighsite=[None] * n_sites,
        list_neighisite=[None] * n_sites,
        structure=structure,
    )
    assert lse.coordination_environments == [[None]] * n_sites
    assert all(len(nb_sets[0]) == 0 for nb_sets in lse.neighbors_sets)
    assert lse.valences is None


def test_from_lobster_non_integer_image(ln_mp190):
    structure = ln_mp190.structure
    shifted = PeriodicSite(structure[1].species, structure[1].frac_coords + 0.25, structure.lattice)
    with pytest.raises(ValueError, match="differences between one site in a periodic image cell is not integer"):
        LobsterLightStructureEnvironments.from_Lobster(
            list_ce_symbol=["S:1", None, None, None],
            list_csm=[0.0, None, None, None],
            list_permutation=[[0], None, None, None],
            list_neighsite=[[shifted], [], [], []],
            list_neighisite=[[1], [], [], []],
            structure=structure,
        )


@pytest.mark.parametrize(
    ("atom1", "expected_conditions"),
    [
        # O5-Ag1: anion-cation bond between different elements involving O
        ("O5", {0, 1, 2, 3, 4}),
        # Ag2-Ag1: cation-cation bond between the same element
        ("Ag2", {0, 5, 6}),
    ],
)
def test_find_relevant_atoms_site_as_second_atom(ln_mp353, atom1, expected_conditions):
    """The central site (Ag1, index 0) can also be the second atom of a bond."""
    icohp = IcohpValue(
        label="1",
        atom1=atom1,
        atom2="Ag1",
        length=2.5,
        translation=(0, 0, 1),
        num=1,
        icohp={Spin.up: -0.5},
    )
    expected_neighbor = LobsterNeighbors._get_atomnumber(atom1)
    for additional_condition in range(7):
        result = ln_mp353._find_relevant_atoms_additional_condition(
            site_idx=0, icohps={"1": icohp}, additional_condition=additional_condition
        )
        if additional_condition in expected_conditions:
            assert result == (["1"], [2.5], [expected_neighbor], [-0.5], [(0, 0, 1)])
        else:
            assert result == ([], [], [], [], [])


def test_determine_unit_cell(ln_mp190):
    site = PeriodicSite("Re", [-0.00001, 1.2, 2.99999], ln_mp190.structure.lattice)
    # coordinates are rounded to 4 decimals before flooring
    assert LobsterNeighbors._determine_unit_cell(site) == [0, 1, 3]


def test_check_icohps_warns():
    with pytest.warns(UserWarning, match="Our neighbor detection might fail"):
        check_ICOHPs(
            lengths_from_ICOHPs=[2.0, 2.005],
            selected_ICOHPs=[-1.0, -0.5],
            translation=[(1, 0, 0), (-1, 0, 0)],
        )

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        # same ICOHPs, or translations that are not opposite: no warning
        check_ICOHPs([2.0, 2.005], [-1.0, -1.0], [(1, 0, 0), (-1, 0, 0)])
        check_ICOHPs([2.0, 2.005], [-1.0, -0.5], [(1, 0, 0), (1, 0, 0)])


def test_get_info_cohps_to_neighbors_errors():
    lobster_nb = _neighbors_mp190(additional_condition=1)
    with pytest.raises(ValueError, match="Please provide either path_to_cohpcar or coxxcar_obj"):
        lobster_nb.get_info_cohps_to_neighbors(path_to_cohpcar=None, isites=[0])

    # NaSi COHPCAR has a different number of bonds than the ReO3 ICOHPLIST
    nasi_cohp = CompleteCohp.from_file(
        filename=f"{TEST_DIR}/COHPCAR.lobster.NaSi.gz",
        fmt="LOBSTER",
        structure_file=f"{TEST_DIR}/CONTCAR.NaSi.gz",
    )
    with pytest.raises(ValueError, match="COHPCAR and ICOHPLIST do not fit together"):
        lobster_nb.get_info_cohps_to_neighbors(coxxcar_obj=nasi_cohp, isites=[0])


def test_plot_cohps_of_neighbors_per_bond_integrated():
    lobster_nb = _neighbors_mp190(additional_condition=1)
    ax = lobster_nb.plot_cohps_of_neighbors(
        path_to_cohpcar=f"{TEST_DIR}/COHPCAR.lobster.mp-190.gz",
        isites=[0],
        per_bond=True,
        integrated=True,
        xlim=(-5, 5),
    )
    try:
        assert ax.get_xlim() == approx((-5, 5))
        assert ax.get_ylim() == approx((-10, 6))
        labels = [line.get_label() for line in ax.get_lines()]
        assert "6 x O-Re (per bond)" in labels
        # the COHPCAR is cached and the per-bond integrated COHP equals the mean ICOHP at E_F
        _label, summed = lobster_nb.get_info_cohps_to_neighbors(path_to_cohpcar=None, isites=[0], per_bond=True)
        per_bond_icohp = lobster_nb.get_info_icohps_to_neighbors(isites=[0]).total_icohp / 6
        assert summed.icohp[Spin.up][700] == approx(per_bond_icohp)
        assert per_bond_icohp == approx(-5.64575, abs=1e-5)
        # CohpPlotter plots -ICOHP on the x axis
        line = ax.get_lines()[labels.index("6 x O-Re (per bond)")]
        assert line.get_xdata()[700] == approx(-per_bond_icohp)
    finally:
        plt.close("all")
