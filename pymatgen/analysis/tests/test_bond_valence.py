from __future__ import annotations

import os

from pytest import approx

from pymatgen.analysis.bond_valence import BVAnalyzer, calculate_bv_sum, calculate_bv_sum_unordered
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class BVAnalyzerTest(PymatgenTest):
    def setUp(self):
        self.analyzer = BVAnalyzer()

    def test_get_valence(self):
        struct = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiMn2O4.json"))
        ans = [1, 1, 3, 3, 4, 4, -2, -2, -2, -2, -2, -2, -2, -2]
        assert self.analyzer.get_valences(struct) == ans
        struct = self.get_structure("LiFePO4")
        ans = [1] * 4 + [2] * 4 + [5] * 4 + [-2] * 16
        assert self.analyzer.get_valences(struct) == ans
        struct = self.get_structure("Li3V2(PO4)3")
        ans = [1] * 6 + [3] * 4 + [5] * 6 + [-2] * 24
        assert self.analyzer.get_valences(struct) == ans
        struct = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "Li4Fe3Mn1(PO4)4.json"))
        ans = [1] * 4 + [2] * 4 + [5] * 4 + [-2] * 16
        assert self.analyzer.get_valences(struct) == ans
        struct = self.get_structure("NaFePO4")
        assert self.analyzer.get_valences(struct) == ans

    def test_get_oxi_state_structure(self):
        struct = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiMn2O4.json"))
        news = self.analyzer.get_oxi_state_decorated_structure(struct)
        assert Species("Mn", 3) in news.composition.elements
        assert Species("Mn", 4) in news.composition.elements


class BondValenceSumTest(PymatgenTest):
    def test_calculate_bv_sum(self):
        struct = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiMn2O4.json"))
        neighbors = struct.get_neighbors(struct[0], 3.0)
        bv_sum = calculate_bv_sum(struct[0], neighbors)
        assert bv_sum == approx(0.7723402182087497)

    def test_calculate_bv_sum_unordered(self):
        struct = Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "LiMn2O4.json"))
        struct[0].species = Composition("Li0.5Na0.5")
        neighbors = struct.get_neighbors(struct[0], 3.0)
        bv_sum = calculate_bv_sum_unordered(struct[0], neighbors)
        assert bv_sum == approx(1.5494662306918852)
