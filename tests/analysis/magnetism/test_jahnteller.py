from __future__ import annotations

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.magnetism.jahnteller import JahnTellerAnalyzer, Species
from pymatgen.core import Lattice, Structure
from tests.testing import TEST_FILES_DIR


class TestJahnTeller:
    def setup_method(self):
        self.jt = JahnTellerAnalyzer()

    def test_jahn_teller_species_analysis(self):
        # 1 d-shell electron
        magnitude = self.jt.get_magnitude_of_effect_from_species("Ti3+", "", "oct")
        assert magnitude == "weak"

        # 2 d-shell electrons
        magnitude = self.jt.get_magnitude_of_effect_from_species("Ti2+", "", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("V3+", "", "oct")
        assert magnitude == "weak"

        # 3
        magnitude = self.jt.get_magnitude_of_effect_from_species("V2+", "", "oct")
        assert magnitude == "none"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Cr3+", "", "oct")
        assert magnitude == "none"

        # 4
        magnitude = self.jt.get_magnitude_of_effect_from_species("Cr2+", "high", "oct")
        assert magnitude == "strong"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Cr2+", "low", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Mn3+", "high", "oct")
        assert magnitude == "strong"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Mn3+", "low", "oct")
        assert magnitude == "weak"

        # 5
        magnitude = self.jt.get_magnitude_of_effect_from_species("Mn2+", "high", "oct")
        assert magnitude == "none"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Mn2+", "low", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Fe3+", "high", "oct")
        assert magnitude == "none"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Fe3+", "low", "oct")
        assert magnitude == "weak"

        # 6
        magnitude = self.jt.get_magnitude_of_effect_from_species("Fe2+", "high", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Fe2+", "low", "oct")
        assert magnitude == "none"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Co3+", "high", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Co3+", "low", "oct")
        assert magnitude == "none"

        # 7
        magnitude = self.jt.get_magnitude_of_effect_from_species("Co2+", "high", "oct")
        assert magnitude == "weak"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Co2+", "low", "oct")
        assert magnitude == "strong"

        # 8
        magnitude = self.jt.get_magnitude_of_effect_from_species("Ni2+", "", "oct")
        assert magnitude == "none"

        # 9
        magnitude = self.jt.get_magnitude_of_effect_from_species("Cu2+", "", "oct")
        assert magnitude == "strong"

        # 10
        magnitude = self.jt.get_magnitude_of_effect_from_species("Cu+", "", "oct")
        assert magnitude == "none"
        magnitude = self.jt.get_magnitude_of_effect_from_species("Zn2+", "", "oct")
        assert magnitude == "none"

    def test_jahn_teller_structure_analysis(self):
        LiFePO4 = Structure.from_file(f"{TEST_FILES_DIR}/cif/LiFePO4.cif", primitive=True)

        Fe3O4 = Structure.from_file(f"{TEST_FILES_DIR}/cif/Fe3O4.cif", primitive=True)

        assert self.jt.is_jahn_teller_active(LiFePO4)
        assert self.jt.is_jahn_teller_active(Fe3O4)

        LiFePO4_analysis = {
            "active": True,
            "strength": "weak",
            "sites": [
                {
                    "ligand": "O2-",
                    "ligand_bond_length_spread": 0.2111,
                    "ligand_bond_lengths": {
                        2.2951,
                        2.2215,
                        2.2383,
                        2.1382,
                        2.084,
                        2.0863,
                    },
                    "strength": "weak",
                    "motif": "oct",
                    "motif_order_parameter": 0.1441,
                    "site_indices": [4, 5, 6, 7],
                    "species": "Fe2+",
                    "spin_state": "unknown",
                }
            ],
        }
        jt_predicted = self.jt.get_analysis(LiFePO4)
        # order does not matter
        jt_predicted["sites"][0]["ligand_bond_lengths"] = set(jt_predicted["sites"][0]["ligand_bond_lengths"])
        assert LiFePO4_analysis == jt_predicted

    def test_mu_so(self):
        SpeciesCo = Species(symbol="Co", oxidation_state=4)
        assert np.sqrt(3) == approx(JahnTellerAnalyzer.mu_so(SpeciesCo, "oct", "low"))
        assert np.sqrt(35) == approx(JahnTellerAnalyzer.mu_so(SpeciesCo, "oct", "high"))
        SpeciesNa = Species(symbol="Na", oxidation_state=1)
        assert None is JahnTellerAnalyzer.mu_so(SpeciesNa, "oct", "high")


MAG_ORDERINGS_DIR = f"{TEST_FILES_DIR}/analysis/magnetic_orderings"


def isolated_tetrahedron(metal: Species, bond_length: float = 2.2) -> Structure:
    """An isolated [MCl4] tetrahedron in a 12 A box, e.g. [FeCl4]2-."""
    offset = bond_length / np.sqrt(3)
    vertices = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    coords = np.array([[0, 0, 0], *(np.array(vertices) * offset)]) + 6
    return Structure(Lattice.cubic(12), [metal, *[Species("Cl", -1)] * 4], coords, coords_are_cartesian=True)


@pytest.fixture
def jt() -> JahnTellerAnalyzer:
    return JahnTellerAnalyzer()


@pytest.fixture
def la_mn_o3() -> Structure:
    return Structure.from_file(f"{MAG_ORDERINGS_DIR}/LaMnO3.json")


class TestJahnTellerStructureAnalysis:
    def test_strong_octahedral(self, jt, la_mn_o3):
        # high-spin Mn3+ (d4, t2g^3 eg^1) is the textbook strong J-T ion
        analysis, structure = jt.get_analysis_and_structure(la_mn_o3)
        assert analysis["active"]
        assert analysis["strength"] == "strong"
        (site,) = analysis["sites"]
        assert site["species"] == "Mn3+"
        assert site["motif"] == "oct"
        assert site["ligand"] == "O2-"
        assert site["site_indices"] == [4, 5, 6, 7]
        assert site["ligand_bond_length_spread"] == approx(0.2675)
        # returned structure is the oxidation-state-decorated primitive cell
        assert sorted(map(str, structure.composition.elements)) == ["La3+", "Mn3+", "O2-"]
        assert jt.is_jahn_teller_active(la_mn_o3)

    def test_op_threshold(self, jt, la_mn_o3):
        # demanding a near-ideal motif rules out every site
        analysis = jt.get_analysis(la_mn_o3, op_threshold=0.99)
        assert analysis == {
            "active": False,
            "sites": [
                {"site_indices": [0, 1, 2, 3], "strength": "none", "reason": "motif='unknown'"},
                {"site_indices": [4, 5, 6, 7], "strength": "none", "reason": "motif='unknown'"},
            ],
        }
        assert not jt.is_jahn_teller_active(la_mn_o3, op_threshold=0.99)

    def test_inactive(self, jt):
        # Ni2+ (d8, t2g^6 eg^2) is not J-T active
        nio = Structure.from_spacegroup(225, Lattice.cubic(4.17), ["Ni", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        assert jt.get_analysis(nio) == {
            "active": False,
            "sites": [
                {
                    "site_indices": [0],
                    "strength": "none",
                    "reason": "Not Jahn-Teller active for this electronic configuration.",
                }
            ],
        }
        assert not jt.is_jahn_teller_active(nio)

    def test_guesstimate_spin(self, jt):
        # octahedral Co3+ in Ca3Co2O6 carries ~1 muB in this calculation, so
        # it is estimated as low spin (t2g^6), which is not J-T active
        structure = Structure.from_file(f"{MAG_ORDERINGS_DIR}/Ca3Co2O6.json")

        analysis = jt.get_analysis(structure)
        assert analysis["strength"] == "weak"
        assert [(site["species"], site["spin_state"]) for site in analysis["sites"]] == [("Co3+", "unknown")]

        analysis = jt.get_analysis(structure, guesstimate_spin=True)
        assert not analysis["active"]
        assert [site["reason"] for site in analysis["sites"]] == [
            "Not Jahn-Teller active for this electronic configuration.",
            "motif='unknown'",  # trigonal prismatic Co
        ]

    def test_guesstimate_spin_single_spin_state(self, jt):
        # Cu2+ (d9) has a single spin state in octahedral coordination
        structure = Structure.from_file(f"{MAG_ORDERINGS_DIR}/CuO.json")
        analysis = jt.get_analysis(structure, guesstimate_spin=True)
        assert analysis["strength"] == "strong"
        assert [(site["species"], site["spin_state"]) for site in analysis["sites"]] == [("Cu2+", "undefined")]

    def test_tetrahedral(self, jt):
        analysis = jt.get_analysis(isolated_tetrahedron(Species("Fe", 2)), calculate_valences=False)
        assert analysis == {
            "active": True,
            "strength": "weak",
            "sites": [
                {
                    "strength": "weak",
                    "motif": "tet",
                    "motif_order_parameter": 1.0,
                    "spin_state": "unknown",
                    "species": "Fe2+",
                    "ligand": "Cl-",
                    "ligand_bond_lengths": [2.2] * 4,
                    "ligand_bond_length_spread": 0.0,
                    "site_indices": [0],
                }
            ],
        }

    def test_analysis_error(self, jt):
        # Pd has no 5s electrons, so its d count cannot be derived
        structure = isolated_tetrahedron(Species("Pd", 2))
        with pytest.raises(AttributeError, match="Invalid element Pd for crystal field calculation"):
            jt.get_analysis(structure, calculate_valences=False)

        with pytest.warns(UserWarning, match="Error analyzing PdCl4: Invalid element Pd"):
            assert not jt.is_jahn_teller_active(structure, calculate_valences=False)

        with pytest.warns(UserWarning, match="Error analyzing PdCl4: Invalid element Pd"):
            tagged = jt.tag_structure(structure, calculate_valences=False)
        assert "possible_jt_active" not in tagged.site_properties

    def test_tag_structure(self, jt, la_mn_o3):
        tagged = jt.tag_structure(la_mn_o3)
        assert tagged.site_properties["possible_jt_active"] == [site.specie.symbol == "Mn" for site in tagged]
        assert tagged.site_properties["possible_jt_active"].count(True) == 4


class TestJahnTellerMagnitude:
    @pytest.mark.parametrize(
        ("species", "expected"),
        [
            ("Ti4+", "none"),  # d0
            ("V4+", "weak"),  # d1, e^1
            ("Fe2+", "weak"),  # d6, e^3 t2^3
            (Species("Cr", 2), "weak"),  # d4, e^2 t2^2
        ],
    )
    def test_from_species_tet(self, jt, species, expected):
        assert jt.get_magnitude_of_effect_from_species(species, "", "tet") == expected

    def test_from_species_object(self, jt):
        assert jt.get_magnitude_of_effect_from_species(Species("Cu", 2), "high", "oct") == "strong"

    def test_from_species_unknown_motif(self, jt):
        assert jt.get_magnitude_of_effect_from_species("Cu2+", "", "square_planar") == "none"

    @pytest.mark.parametrize("species", ["Fe", "O2-"])
    def test_from_species_no_data(self, jt, species):
        # need a transition metal Species with an oxidation state
        with pytest.warns(UserWarning, match="No data for this species"):
            assert jt.get_magnitude_of_effect_from_species(species, "", "oct") == "none"

    @pytest.mark.parametrize(
        ("species", "match"),
        [
            ("Pd2+", "Invalid element Pd for crystal field calculation"),
            ("Fe9+", r"Invalid oxidation state 9\.0 for element Fe"),
        ],
    )
    def test_from_species_invalid(self, jt, species, match):
        with pytest.raises(AttributeError, match=match):
            jt.get_magnitude_of_effect_from_species(species, "", "oct")

    @pytest.mark.parametrize(
        ("motif", "spin_config", "expected"),
        [
            ("oct", {"e_g": 0, "t_2g": 3}, "none"),
            ("oct", {"e_g": 2, "t_2g": 6}, "none"),
            ("oct", {"e_g": 0, "t_2g": 4}, "weak"),
            ("oct", {"e_g": 2, "t_2g": 5}, "weak"),
            ("oct", {"e_g": 1, "t_2g": 3}, "strong"),
            ("oct", {"e_g": 3, "t_2g": 6}, "strong"),
            ("tet", {"e": 0, "t_2": 0}, "none"),
            ("tet", {"e": 1, "t_2": 0}, "weak"),
            ("tet", {"e": 3, "t_2": 3}, "weak"),
            ("square_planar", {}, "none"),
        ],
    )
    def test_from_spin_config(self, motif, spin_config, expected):
        assert JahnTellerAnalyzer.get_magnitude_of_effect_from_spin_config(motif, spin_config) == expected

    @pytest.mark.xfail(
        reason="tetrahedral rule tests e % 3 and t_2 % 2 instead of e % 2 and t_2 % 3, "
        "so filled/half-filled shells (d5 high spin, d10) are reported as 'weak'",
        strict=True,
    )
    @pytest.mark.parametrize("species", ["Mn2+", "Zn2+"])
    def test_from_species_tet_closed_shell(self, jt, species):
        assert jt.get_magnitude_of_effect_from_species(species, "", "tet") == "none"


class TestEstimateSpinState:
    @pytest.mark.parametrize(
        ("species", "motif", "magmom", "expected"),
        [
            ("Fe2+", "oct", 4.9, "high"),  # mu_so(high) = sqrt(24)
            ("Fe2+", "oct", 6.0, "high"),  # above high-spin value
            ("Fe2+", "oct", 0.0, "low"),  # mu_so(low) = 0
            ("Fe2+", "oct", 2.5, "unknown"),  # too far from either
            ("Mn3+", "oct", 5.0, "high"),
            ("Mn3+", "oct", 2.8, "low"),
            ("Fe3+", "tet", 5.9, "high"),
            ("Ni2+", "oct", 2.8, "undefined"),  # only one spin state
            ("O2-", "oct", 0.0, "undefined"),  # no crystal field data
        ],
    )
    def test_estimate_spin_state(self, species, motif, magmom, expected):
        assert JahnTellerAnalyzer._estimate_spin_state(species, motif, magmom) == expected

    @pytest.mark.xfail(
        reason="the sign of the moment is not discarded, so down-spin sites in AFM "
        "structures are always classified as low spin",
        strict=True,
    )
    def test_estimate_spin_state_negative_moment(self):
        assert JahnTellerAnalyzer._estimate_spin_state("Mn3+", "oct", -5.0) == "high"
