from __future__ import annotations

import unittest

from pytest import approx

from pymatgen.analysis.disorder import get_warren_cowley_parameters
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest


class OrderParameterTest(PymatgenTest):
    def test_compute_warren_cowley_parameters(self):
        struct = Structure.from_prototype("CsCl", ["Mo", "W"], a=4)
        aij = get_warren_cowley_parameters(struct, r=3.4, dr=0.3)
        assert aij[(Element.Mo, Element.W)] == approx(-1.0)
        aij = get_warren_cowley_parameters(struct, r=4, dr=0.2)
        assert aij[(Element.Mo, Element.Mo)] == approx(1.0)
        struct = Structure.from_prototype("CsCl", ["Mo", "W"], a=4)
        struct = struct * 4

        # Swap the first and last sites to cause disorder
        struct[0] = "W"
        struct[len(struct) - 1] = "Mo"

        aij = get_warren_cowley_parameters(struct, r=3.4, dr=0.3)
        assert aij[(Element.Mo, Element.W)] == approx(-0.9453125)
        assert aij[(Element.Mo, Element.W)] == aij[(Element.W, Element.Mo)]


if __name__ == "__main__":
    unittest.main()
