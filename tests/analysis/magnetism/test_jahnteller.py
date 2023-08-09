from __future__ import annotations

import unittest

import numpy as np
from pytest import approx

from pymatgen.analysis.magnetism.jahnteller import JahnTellerAnalyzer, Species
from pymatgen.io.cif import CifParser
from pymatgen.util.testing import TEST_FILES_DIR


class TestJahnTeller(unittest.TestCase):
    def setUp(self):
        self.jt = JahnTellerAnalyzer()

    def test_jahn_teller_species_analysis(self):
        # 1 d-shell electron
        m = self.jt.get_magnitude_of_effect_from_species("Ti3+", "", "oct")
        assert m == "weak"

        # 2 d-shell electrons
        m = self.jt.get_magnitude_of_effect_from_species("Ti2+", "", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("V3+", "", "oct")
        assert m == "weak"

        # 3
        m = self.jt.get_magnitude_of_effect_from_species("V2+", "", "oct")
        assert m == "none"
        m = self.jt.get_magnitude_of_effect_from_species("Cr3+", "", "oct")
        assert m == "none"

        # 4
        m = self.jt.get_magnitude_of_effect_from_species("Cr2+", "high", "oct")
        assert m == "strong"
        m = self.jt.get_magnitude_of_effect_from_species("Cr2+", "low", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("Mn3+", "high", "oct")
        assert m == "strong"
        m = self.jt.get_magnitude_of_effect_from_species("Mn3+", "low", "oct")
        assert m == "weak"

        # 5
        m = self.jt.get_magnitude_of_effect_from_species("Mn2+", "high", "oct")
        assert m == "none"
        m = self.jt.get_magnitude_of_effect_from_species("Mn2+", "low", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("Fe3+", "high", "oct")
        assert m == "none"
        m = self.jt.get_magnitude_of_effect_from_species("Fe3+", "low", "oct")
        assert m == "weak"

        # 6
        m = self.jt.get_magnitude_of_effect_from_species("Fe2+", "high", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("Fe2+", "low", "oct")
        assert m == "none"
        m = self.jt.get_magnitude_of_effect_from_species("Co3+", "high", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("Co3+", "low", "oct")
        assert m == "none"

        # 7
        m = self.jt.get_magnitude_of_effect_from_species("Co2+", "high", "oct")
        assert m == "weak"
        m = self.jt.get_magnitude_of_effect_from_species("Co2+", "low", "oct")
        assert m == "strong"

        # 8
        m = self.jt.get_magnitude_of_effect_from_species("Ni2+", "", "oct")
        assert m == "none"

        # 9
        m = self.jt.get_magnitude_of_effect_from_species("Cu2+", "", "oct")
        assert m == "strong"

        # 10
        m = self.jt.get_magnitude_of_effect_from_species("Cu+", "", "oct")
        assert m == "none"
        m = self.jt.get_magnitude_of_effect_from_species("Zn2+", "", "oct")
        assert m == "none"

    def test_jahn_teller_structure_analysis(self):
        parser = CifParser(f"{TEST_FILES_DIR}/LiFePO4.cif")
        LiFePO4 = parser.get_structures()[0]

        parser = CifParser(f"{TEST_FILES_DIR}/Fe3O4.cif")
        Fe3O4 = parser.get_structures()[0]

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
