from __future__ import annotations

from unittest import TestCase

import numpy as np
from pytest import approx

from pymatgen.analysis.magnetism.jahnteller import JahnTellerAnalyzer, Species
from pymatgen.core import Structure
from pymatgen.util.testing import TEST_FILES_DIR


class TestJahnTeller(TestCase):
    def setUp(self):
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
                    "ligand_bond_lengths": {2.2951, 2.2215, 2.2383, 2.1382, 2.084, 2.0863},
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
