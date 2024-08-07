from __future__ import annotations

import os
from unittest import TestCase

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core import Element
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import Cohp, CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster import Charge, Icohplist
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.util.testing import TEST_FILES_DIR

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__email__ = "janine.george@uclouvain.be"
__date__ = "Jan 14, 2021"

TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/cohp/environments"
module_dir = os.path.dirname(os.path.abspath(__file__))


class TestLobsterNeighbors(TestCase):
    def setUp(self):
        # test additional conditions first
        # only consider cation anion bonds

        self.chem_env_lobster1 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=1,
        )

        # all bonds
        self.chem_env_lobster0 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=0,
        )

        # only cation-cation, anion-anion bonds
        self.chem_env_lobster5 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=5,
        )

        # only cation-cation bonds
        self.chem_env_lobster6 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=6,
        )

        # 2,3,4 are not tested so far
        self.chem_env_lobster2 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=2,
        )

        self.chem_env_lobster3 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=3,
        )

        self.chem_env_lobster4 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
            additional_condition=4,
        )

        # search for other testcase where 2,3,4 arrive at different results
        self.chem_env_lobster0_second = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=0,
        )
        self.chem_env_lobster1_second = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=1,
        )

        self.chem_env_lobster2_second = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=2,
        )

        self.chem_env_lobster5_second = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=5,
        )

        self.chem_env_lobster5_second_percentage = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=5,
            perc_strength_icohp=1.0,
        )

        self.chem_env_lobster6_second = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            additional_condition=6,
        )
        # coop / cobi
        self.chem_env_lobster1_coop_NaCl = LobsterNeighbors(
            are_coops=True,
            filename_icohp=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_cobi_NaCl = LobsterNeighbors(
            are_coops=True,
            filename_icohp=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chem_env_lobster1_cobi_mp470 = LobsterNeighbors(
            are_coops=True,
            filename_icohp=f"{TEST_DIR}/ICOBILIST.lobster.mp_470.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_470.gz"),
            additional_condition=1,
        )

        # TODO: use charge instead of valence
        self.chem_env_lobster1_charges = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
        )
        self.chem_env_lobster1_charges_noisecutoff = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_632319.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_632319.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp_632319.gz",
            additional_condition=1,
            perc_strength_icohp=0.05,
            noise_cutoff=0.1,
        )
        self.chem_env_lobster1_charges_wo_noisecutoff = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_632319.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_632319.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp_632319.gz",
            additional_condition=1,
            perc_strength_icohp=0.05,
            noise_cutoff=None,
        )
        self.chem_env_lobster1_charges_loewdin = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            which_charge="Loewdin",
        )
        self.chem_env_lobster6_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=6,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster5_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=5,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster4_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=4,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster3_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=3,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster2_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=2,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster1_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            adapt_extremum_to_add_cond=True,
        )

        self.chem_env_lobster0_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster0_NaSi = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaSi.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.NaSi.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chem_env_lobster_NaSi_wo_charges = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaSi.gz"),
            valences_from_charges=False,
            filename_charge=None,
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        # Test LobsterNeighbors using pymatgen objects
        self.obj_icohp = Icohplist(filename=f"{TEST_DIR}/ICOHPLIST.lobster.NaSi.gz")
        self.obj_charge = Charge(filename=f"{TEST_DIR}/CHARGE.lobster.NaSi.gz")
        self.chem_env_w_obj = LobsterNeighbors(
            filename_icohp=None,
            are_coops=False,
            obj_icohp=self.obj_icohp,
            obj_charge=self.obj_charge,
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaSi.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.NaSi.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )

    def test_cation_anion_mode_without_ions(self):
        with pytest.raises(
            ValueError, match="Valences cannot be assigned, additional_conditions 1, 3, 5 and 6 will not work"
        ):
            _ = LobsterNeighbors(
                are_coops=False,
                filename_icohp=f"{TEST_DIR}/../ICOHPLIST.lobster",
                structure=Structure.from_file(f"{TEST_DIR}/../POSCAR"),
                valences_from_charges=False,
                additional_condition=1,
            )
        with pytest.raises(
            ValueError, match="All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work"
        ):
            _ = LobsterNeighbors(
                are_coops=False,
                filename_icohp=f"{TEST_DIR}/../ICOHPLIST.lobster",
                structure=Structure.from_file(f"{TEST_DIR}/../POSCAR"),
                valences_from_charges=False,
                additional_condition=1,
                valences=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            )

    def test_wrong_additional_correction(self):
        with pytest.raises(
            ValueError, match=r"Unexpected additional_condition=10, must be one of \[0, 1, 2, 3, 4, 5, 6\]"
        ):
            LobsterNeighbors(
                are_coops=False,
                filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
                structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                valences_from_charges=True,
                filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
                additional_condition=10,
            )

    def test_set_limits(self):
        test = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            limits=[-100000, 0],
        )
        assert test.limits == [-100000, 0]

    def test_molecules_allowed(self):
        assert not self.chem_env_lobster1.molecules_allowed

    def test_get_anion_types(self):
        assert self.chem_env_lobster0_second.get_anion_types() == {Element("O")}
        assert self.chem_env_lobster0_second.anion_types == {Element("O")}

    def test_get_nn_info(self):
        # NO_ADDITIONAL_CONDITION = 0
        # ONLY_ANION_CATION_BONDS = 1
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        # ONLY_CATION_CATION_BONDS=6

        # All bonds
        # ReO3
        assert (
            len(
                self.chem_env_lobster0.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster0.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ANION_CATION_BONDS = 1
        assert (
            len(
                self.chem_env_lobster1.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster1.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
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
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        assert (
            len(
                self.chem_env_lobster2.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster2.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        assert (
            len(
                self.chem_env_lobster3.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster3.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        assert (
            len(
                self.chem_env_lobster4.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster4.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        assert (
            len(
                self.chem_env_lobster5.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 0
        )
        assert (
            len(
                self.chem_env_lobster5.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 0
        )
        # ONLY_CATION_CATION_BONDS=6
        assert (
            len(
                self.chem_env_lobster6.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 0
        )

        assert (
            len(
                self.chem_env_lobster6.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
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
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 8
        )

        # ONLY_ANION_CATION_BONDS = 1
        assert (
            len(
                self.chem_env_lobster1_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 2
        )

        assert (
            len(
                self.chem_env_lobster1_coop_NaCl.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chem_env_lobster1_cobi_NaCl.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chem_env_lobster1_cobi_mp470.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_470.gz"),
                    n=3,
                )
            )
            == 3
        )

        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        assert (
            len(
                self.chem_env_lobster2_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 2
        )
        assert (
            len(
                self.chem_env_lobster2_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 4
        )

        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        assert (
            len(
                self.chem_env_lobster5_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster5_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 0
        )
        # ONLY_CATION_CATION_BONDS=6
        assert (
            len(
                self.chem_env_lobster6_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chem_env_lobster6_second.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 0
        )

        assert (
            len(
                self.chem_env_lobster5_second_percentage.get_nn(
                    structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 0
        )

    def test_structure_graph(self):
        sg = self.chem_env_lobster1_second.get_bonded_structure(
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz")
        )
        assert isinstance(sg, StructureGraph)

    def test_extended_structure_graph(self):
        self.chem_env_lobsterNaCl = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
            valences_from_charges=True,
            filename_charge=f"{TEST_DIR}/CHARGE.lobster.NaCl.gz",
            filename_blist_sg1=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
            filename_blist_sg2=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
            add_additional_data_sg=True,
            id_blist_sg1="icobi",
            id_blist_sg2="icoop",
            additional_condition=1,
        )
        sg = self.chem_env_lobsterNaCl.get_bonded_structure(
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
            decorate=True,
            edge_properties=True,
            weights=True,
        )
        assert sg.graph.get_edge_data(0, 1)[0]["ICOHP"] == approx(-0.56541)
        assert sg.graph.get_edge_data(0, 1)[0]["ICOBI"] == approx(0.08484)
        assert sg.graph.get_edge_data(0, 1)[0]["ICOOP"] == approx(0.02826)
        assert sg.graph.get_edge_data(0, 1)[0]["bond_label"] == "21"
        assert sg.graph.get_edge_data(0, 1)[5]["bond_label"] == "30"
        assert isinstance(sg, StructureGraph)

    def test_raises_extended_structure_graph(self):
        with pytest.raises(ValueError, match="Algorithm can only work with ICOOPs, ICOBIs"):
            self.chem_env_lobsterNaCl = LobsterNeighbors(
                are_coops=False,
                filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.NaCl.gz",
                structure=Structure.from_file(f"{TEST_DIR}/POSCAR.NaCl.gz"),
                valences_from_charges=True,
                filename_charge=f"{TEST_DIR}/CHARGE.lobster.NaCl.gz",
                filename_blist_sg1=f"{TEST_DIR}/ICOBILIST.lobster.NaCl.gz",
                filename_blist_sg2=f"{TEST_DIR}/ICOOPLIST.lobster.NaCl.gz",
                add_additional_data_sg=True,
                id_blist_sg1="icopppp",
                id_blist_sg2="icoop",
                additional_condition=1,
            )

    def test_order_parameter(self):
        assert self.chem_env_lobster1_second.get_local_order_parameters(
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_353.gz"), n=0
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
        results = self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])
        assert results[0] == approx(-33.26058)
        for bond in results[1]:
            assert bond == approx(-5.54345, abs=1e-3)
        assert results[2] == 6
        assert results[3] == ["27", "30", "48", "49", "64", "73"]

        results2 = self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=None)
        assert results2[0] == approx(-33.26058)
        for bond in results2[1]:
            assert bond == approx(-5.54345, abs=1e-3)
        assert results2[2] == 6
        assert results2[3] == ["27", "30", "48", "49", "64", "73"]
        assert results2[4] == [
            ["Re1", "O2"],
            ["Re1", "O2"],
            ["Re1", "O3"],
            ["Re1", "O3"],
            ["Re1", "O4"],
            ["Re1", "O4"],
        ]

    def test_get_sum_icohps_between_neighbors_of_atom(self):
        # will only look at icohps between cations or anions
        self.chem_env_lobster1.get_info_icohps_to_neighbors(isites=[1])
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[1])[2] == 1
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[1])[0] == approx(-0.05507)
        assert self.chem_env_lobster1.get_info_icohps_between_neighbors(isites=[0])[2] == 15
        # use an example where this is easier to test (e.g., linear environment?)

        chemenv_here = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp-7000.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp-7000.gz"),
            additional_condition=1,
        )
        assert len(chemenv_here.get_info_icohps_between_neighbors(isites=[0])[4]) == 6

    def test_get_plot_label(self):
        label = self.chem_env_lobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
            per_bond=False,
        )
        assert label == "6 x O-Re"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
            per_bond=False,
        )
        assert label == "5 x O-Re, 1 x O-Si"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[["Si1", "O2"], ["Si1", "O2"], ["Si1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
            per_bond=False,
        )
        assert label == "4 x O-Si, 2 x O-Re"

        label = self.chem_env_lobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
            per_bond=True,
        )
        assert label == "6 x O-Re (per bond)"

    def test_get_info_cohps_to_neighbors(self):
        chem_env_lobster1 = LobsterNeighbors(
            are_coops=False,
            filename_icohp=f"{TEST_DIR}/ICOHPLIST.lobster.mp_190_2.gz",
            structure=Structure.from_file(f"{TEST_DIR}/POSCAR.mp_190.gz"),
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

        coph_thing = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
        )[1]
        assert np.sum([coph_thing.icohp[Spin.up], coph_thing.icohp[Spin.down]], axis=0)[300] == approx(
            chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])[0]
        )

        # summed_spin_channel
        coph_thing = chem_env_lobster1.get_info_cohps_to_neighbors(
            path_to_cohpcar=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
            summed_spin_channels=True,
        )[1]
        assert coph_thing.icohp[Spin.up][300] == approx(chem_env_lobster1.get_info_icohps_to_neighbors(isites=[0])[0])

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
            filename=f"{TEST_DIR}/COHPCAR.lobster.NaSi.gz", fmt="LOBSTER", structure_file=f"{TEST_DIR}/POSCAR.NaSi.gz"
        )
        plot_label_obj, _summed_cohpcar_NaSi_obj = self.chem_env_w_obj.get_info_cohps_to_neighbors(
            obj_cohpcar=obj_cohpcar,
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
                path_to_cohpcar=cohpcar_lobster_mp_190,
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

        with pytest.raises(ValueError, match=expected_msg):
            # icohplist and cohpcar do not fit together
            self.chem_env_lobster2.get_info_cohps_to_neighbors(
                path_to_cohpcar=cohpcar_lobster_mp_190,
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

    def test_valences(self):
        assert self.chem_env_lobster1_charges_noisecutoff.valences == [0.75, -0.75]  # Mulliken
        assert self.chem_env_lobster1_charges_loewdin.valences == [0.27, 0.27, 0.27, 0.27, -0.54, -0.54]
        assert self.chem_env_w_obj.valences == [0.67] * 4 + [0.7] * 4 + [-0.7] * 4 + [-0.68] * 4  # charge_obj
        assert self.chem_env_lobster_NaSi_wo_charges.valences == [1] * 8 + [-1] * 8  # BVA
