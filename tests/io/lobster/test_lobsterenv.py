from __future__ import annotations

import os
import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core import Element
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.cohp import Cohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster.lobsterenv import LobsterNeighbors
from pymatgen.util.testing import TEST_FILES_DIR

__author__ = "Janine George"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.1"
__email__ = "janine.george@uclouvain.be"
__date__ = "Jan 14, 2021"

test_dir_env = f"{TEST_FILES_DIR}/cohp/environments"
this_dir = os.path.dirname(os.path.abspath(__file__))


class TestLobsterNeighbors(unittest.TestCase):
    def setUp(self):
        # test additional conditions first
        # only consider cation anion bonds

        self.chemenvlobster1 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=1,
        )

        # all bonds
        self.chemenvlobster0 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=0,
        )

        # only cation cation, anion anion bonds
        self.chemenvlobster5 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=5,
        )

        # only cation cation bonds
        self.chemenvlobster6 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=6,
        )

        # 2,3,4 are not tested so far
        self.chemenvlobster2 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=2,
        )

        self.chemenvlobster3 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=3,
        )

        self.chemenvlobster4 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=4,
        )

        # search for other testcase where 2,3,4 arrive at different results
        self.chemenvlobster0_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=0,
        )
        self.chemenvlobster1_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=1,
        )

        self.chemenvlobster2_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=2,
        )

        self.chemenvlobster5_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=5,
        )

        self.chemenvlobster5_second_percentage = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=5,
            perc_strength_ICOHP=1.0,
        )

        self.chemenvlobster6_second = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            additional_condition=6,
        )
        # coop / cobi
        self.chemenvlobster1_coop_NaCl = LobsterNeighbors(
            are_coops=True,
            filename_ICOHP=f"{test_dir_env}/ICOOPLIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chemenvlobster1_cobi_NaCl = LobsterNeighbors(
            are_coops=True,
            filename_ICOHP=f"{test_dir_env}/ICOBILIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
            additional_condition=1,
            noise_cutoff=None,
        )

        self.chemenvlobster1_cobi_mp470 = LobsterNeighbors(
            are_coops=True,
            filename_ICOHP=f"{test_dir_env}/ICOBILIST.lobster.mp_470.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_470.gz"),
            additional_condition=1,
        )

        # TODO: use charge instead of valence
        self.chemenvlobster1_charges = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
        )
        self.chemenvlobster1_charges_noisecutoff = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_632319.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_632319.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp_632319.gz",
            additional_condition=1,
            perc_strength_ICOHP=0.05,
            noise_cutoff=0.1,
        )
        self.chemenvlobster1_charges_wo_noisecutoff = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_632319.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_632319.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp_632319.gz",
            additional_condition=1,
            perc_strength_ICOHP=0.05,
            noise_cutoff=None,
        )
        self.chemenvlobster1_charges_loewdin = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            which_charge="Loewdin",
        )
        self.chemenvlobster6_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=6,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster5_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=5,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster4_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=4,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster3_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=3,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster2_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=2,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster1_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            adapt_extremum_to_add_cond=True,
        )

        self.chemenvlobster0_charges_additional_condition = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )
        self.chemenvlobster0_NaSi = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.NaSi.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaSi.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.NaSi.gz",
            additional_condition=0,
            adapt_extremum_to_add_cond=True,
        )

    def test_cation_anion_mode_without_ions(self):
        with pytest.raises(
            ValueError, match="Valences cannot be assigned, additional_conditions 1, 3, 5 and 6 will not work"
        ):
            _ = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=f"{test_dir_env}/../ICOHPLIST.lobster",
                structure=Structure.from_file(f"{test_dir_env}/../POSCAR"),
                valences_from_charges=False,
                additional_condition=1,
            )
        with pytest.raises(
            ValueError, match="All valences are equal to 0, additional_conditions 1, 3, 5 and 6 will not work"
        ):
            _ = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=f"{test_dir_env}/../ICOHPLIST.lobster",
                structure=Structure.from_file(f"{test_dir_env}/../POSCAR"),
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
                filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
                structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                valences_from_charges=True,
                filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
                additional_condition=10,
            )

    def test_set_limits(self):
        test = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_353.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.mp-353.gz",
            additional_condition=1,
            limits=[-100000, 0],
        )
        assert test.limits == [-100000, 0]

    def test_molecules_allowed(self):
        assert not self.chemenvlobster1.molecules_allowed

    def test_get_anion_types(self):
        assert self.chemenvlobster0_second.get_anion_types() == {Element("O")}
        assert self.chemenvlobster0_second.anion_types == {Element("O")}

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
                self.chemenvlobster0.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster0.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ANION_CATION_BONDS = 1
        assert (
            len(
                self.chemenvlobster1.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster1.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        assert (
            len(
                self.chemenvlobster1_charges_noisecutoff.get_nn(
                    structure=self.chemenvlobster1_charges_noisecutoff.structure,
                    n=1,
                )
            )
            == 0
        )
        assert (
            len(
                self.chemenvlobster1_charges_wo_noisecutoff.get_nn(
                    structure=self.chemenvlobster1_charges_wo_noisecutoff.structure,
                    n=1,
                )
            )
            == 8
        )
        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        assert (
            len(
                self.chemenvlobster2.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster2.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ANION_CATION_BONDS_AND_NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 3
        assert (
            len(
                self.chemenvlobster3.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster3.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # ONLY_ELEMENT_TO_OXYGEN_BONDS = 4
        assert (
            len(
                self.chemenvlobster4.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster4.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 2
        )
        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        assert (
            len(
                self.chemenvlobster5.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 0
        )
        assert (
            len(
                self.chemenvlobster5.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=1,
                )
            )
            == 0
        )
        # ONLY_CATION_CATION_BONDS=6
        assert (
            len(
                self.chemenvlobster6.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
                    n=0,
                )
            )
            == 0
        )

        assert (
            len(
                self.chemenvlobster6.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
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
                self.chemenvlobster0_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 8
        )

        # ONLY_ANION_CATION_BONDS = 1
        assert (
            len(
                self.chemenvlobster1_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 2
        )

        assert (
            len(
                self.chemenvlobster1_coop_NaCl.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chemenvlobster1_cobi_NaCl.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
                    n=0,
                )
            )
            == 6
        )

        assert (
            len(
                self.chemenvlobster1_cobi_mp470.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_470.gz"),
                    n=3,
                )
            )
            == 3
        )

        # NO_ELEMENT_TO_SAME_ELEMENT_BONDS = 2
        assert (
            len(
                self.chemenvlobster2_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 2
        )
        assert (
            len(
                self.chemenvlobster2_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 4
        )

        # DO_NOT_CONSIDER_ANION_CATION_BONDS=5
        assert (
            len(
                self.chemenvlobster5_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster5_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 0
        )
        # ONLY_CATION_CATION_BONDS=6
        assert (
            len(
                self.chemenvlobster6_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 6
        )
        assert (
            len(
                self.chemenvlobster6_second.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=4,
                )
            )
            == 0
        )

        assert (
            len(
                self.chemenvlobster5_second_percentage.get_nn(
                    structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"),
                    n=0,
                )
            )
            == 0
        )

    def test_structure_graph(self):
        sg = self.chemenvlobster1_second.get_bonded_structure(
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz")
        )
        assert isinstance(sg, StructureGraph)

    def test_extended_structure_graph(self):
        self.chemenvlobsterNaCl = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.NaCl.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
            valences_from_charges=True,
            filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.NaCl.gz",
            filename_blist_sg1=f"{test_dir_env}/ICOBILIST.lobster.NaCl.gz",
            filename_blist_sg2=f"{test_dir_env}/ICOOPLIST.lobster.NaCl.gz",
            add_additional_data_sg=True,
            id_blist_sg1="icobi",
            id_blist_sg2="icoop",
            additional_condition=1,
        )
        sg = self.chemenvlobsterNaCl.get_bonded_structure(
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
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
            self.chemenvlobsterNaCl = LobsterNeighbors(
                are_coops=False,
                filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.NaCl.gz",
                structure=Structure.from_file(f"{test_dir_env}/POSCAR.NaCl.gz"),
                valences_from_charges=True,
                filename_CHARGE=f"{test_dir_env}/CHARGE.lobster.NaCl.gz",
                filename_blist_sg1=f"{test_dir_env}/ICOBILIST.lobster.NaCl.gz",
                filename_blist_sg2=f"{test_dir_env}/ICOOPLIST.lobster.NaCl.gz",
                add_additional_data_sg=True,
                id_blist_sg1="icopppp",
                id_blist_sg2="icoop",
                additional_condition=1,
            )

    def test_order_parameter(self):
        assert self.chemenvlobster1_second.get_local_order_parameters(
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_353.gz"), n=0
        )["linear"] == approx(1.0)

    def test_get_structure_environments(self):
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        assert lse.coordination_environments[0][0]["ce_symbol"] == "L:2"
        assert lse.coordination_environments[5][0]["ce_symbol"] == "T:4"

        lse2 = self.chemenvlobster1.get_light_structure_environment()
        assert lse2.coordination_environments[0][0]["ce_symbol"] == "O:6"

    def test_get_strucuture_environments_further_tests(self):
        lse = self.chemenvlobster1_second.get_light_structure_environment()
        lse.as_dict()
        lse.get_statistics()
        assert lse.uniquely_determines_coordination_environments

    def test_get_info_icohps_neighbors(self):
        results = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])
        assert results[0] == approx(-33.26058)
        for bond in results[1]:
            assert bond == approx(-5.54345, abs=1e-3)
        assert results[2] == 6
        assert results[3] == ["27", "30", "48", "49", "64", "73"]

        results2 = self.chemenvlobster1.get_info_icohps_to_neighbors(isites=None)
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
        self.chemenvlobster1.get_info_icohps_to_neighbors(isites=[1])
        assert self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[2] == 1
        assert self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[1])[0] == approx(-0.05507)
        assert self.chemenvlobster1.get_info_icohps_between_neighbors(isites=[0])[2] == 15
        # use an example where this is easier to test (e.g., linear environment?)

        chemenv_here = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp-7000.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp-7000.gz"),
            additional_condition=1,
        )
        assert len(chemenv_here.get_info_icohps_between_neighbors(isites=[0])[4]) == 6

    def test_get_plot_label(self):
        label = self.chemenvlobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
            per_bond=False,
        )
        assert label == "6 x O-Re"

        label = self.chemenvlobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
            per_bond=False,
        )
        assert label == "5 x O-Re, 1 x O-Si"

        label = self.chemenvlobster1._get_plot_label(
            atoms=[["Si1", "O2"], ["Si1", "O2"], ["Si1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Si1", "O4"]],
            per_bond=False,
        )
        assert label == "4 x O-Si, 2 x O-Re"

        label = self.chemenvlobster1._get_plot_label(
            atoms=[["Re1", "O2"], ["Re1", "O2"], ["Re1", "O3"], ["Re1", "O3"], ["Re1", "O4"], ["Re1", "O4"]],
            per_bond=True,
        )
        assert label == "6 x O-Re (per bond)"

    def test_get_info_cohps_to_neighbors(self):
        chemenvlobster1 = LobsterNeighbors(
            are_coops=False,
            filename_ICOHP=f"{test_dir_env}/ICOHPLIST.lobster.mp_190_2.gz",
            structure=Structure.from_file(f"{test_dir_env}/POSCAR.mp_190.gz"),
            additional_condition=1,
        )
        cohpcar_lobster_mp_190 = f"{test_dir_env}/COHPCAR.lobster.mp-190.gz"
        plot_label, summed_cohpcar_mp_190 = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
        )
        assert plot_label == "6 x O-Re (per bond)"
        assert isinstance(summed_cohpcar_mp_190, Cohp)

        coph_thing = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
        )[1]
        assert np.sum([coph_thing.icohp[Spin.up], coph_thing.icohp[Spin.down]], axis=0)[300] == approx(
            chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])[0]
        )

        # summed_spin_channel
        coph_thing = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=None,
            per_bond=False,
            summed_spin_channels=True,
        )[1]
        assert coph_thing.icohp[Spin.up][300] == approx(chemenvlobster1.get_info_icohps_to_neighbors(isites=[0])[0])

        plot_label, summed_cohpcar_mp_190_Te = chemenvlobster1.get_info_cohps_to_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["Te"],
        )

        assert plot_label is None
        assert summed_cohpcar_mp_190_Te is None

        plot_label, _summed_cohpcar_NaSi = self.chemenvlobster0_NaSi.get_info_cohps_to_neighbors(
            path_to_COHPCAR=f"{test_dir_env}/COHPCAR.lobster.NaSi.gz",
            isites=[8],
            onlycation_isites=False,
            only_bonds_to=["Na"],
        )
        assert plot_label == "1 x Na-Si (per bond)"
        assert (
            self.chemenvlobster0_NaSi.get_info_cohps_to_neighbors(
                path_to_COHPCAR=f"{test_dir_env}/COHPCAR.lobster.NaSi.gz",
                isites=[8],
                onlycation_isites=False,
                only_bonds_to=["Si"],
            )[0]
            == "3 x Si-Si (per bond)"
        )

        chemenvlobster1.plot_cohps_of_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
        )

        chemenvlobster1.plot_cohps_of_neighbors(
            path_to_COHPCAR=cohpcar_lobster_mp_190,
            isites=[0],
            only_bonds_to=["O"],
            summed_spin_channels=True,
            xlim=[-10, 10],
            ylim=None,
        )

        expected_msg = "COHPCAR and ICOHPLIST do not fit together"
        with pytest.raises(ValueError, match=expected_msg):
            # icohplist and cohpcar do not fit together
            self.chemenvlobster1.get_info_cohps_to_neighbors(
                path_to_COHPCAR=cohpcar_lobster_mp_190,
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )

        with pytest.raises(ValueError, match=expected_msg):
            # icohplist and cohpcar do not fit together
            self.chemenvlobster2.get_info_cohps_to_neighbors(
                path_to_COHPCAR=cohpcar_lobster_mp_190,
                isites=[0],
                only_bonds_to=None,
                per_bond=False,
            )
