from __future__ import annotations

from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.io.atat import Mcsqs
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/mcsqs"


class TestAtat(PymatgenTest):
    def test_mcsqs_import(self):
        test_string = """1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
0.000000 -1.000000 -2.000000
2.000000 -1.000000 0.000000
-1.000000 -1.000000 1.000000
0.000000 -2.000000 -1.000000 Mn
1.000000 -2.000000 -1.000000 Mn
0.000000 -1.000000 -1.000000 Mn
-0.000000 -2.000000 0.000000 Mn
1.000000 -2.000000 0.000000 Mn
0.000000 -1.000000 0.000000 Mn
1.000000 -1.000000 0.000000 Fe
1.000000 -3.000000 -1.000000 Mn
0.500000 -1.500000 -0.500000 Sr
1.500000 -1.500000 -0.500000 Ca
-0.500000 -1.500000 0.500000 Ca
0.500000 -1.500000 0.500000 Ca
1.500000 -2.500000 -1.500000 Ca
0.500000 -1.500000 -1.500000 Sr
0.500000 -2.500000 -0.500000 Sr
-0.500000 -1.500000 -0.500000 Ca
0.000000 -1.500000 -1.000000 O
1.000000 -1.500000 -1.000000 O
1.000000 -2.500000 0.000000 O
-0.000000 -1.500000 0.000000 O
1.000000 -1.500000 0.000000 O
0.000000 -0.500000 0.000000 O
0.000000 -2.500000 -1.000000 O
1.000000 -2.500000 -1.000000 O
0.500000 -2.000000 -1.000000 O
1.500000 -2.000000 -1.000000 O
0.500000 -1.000000 -1.000000 O
0.500000 -2.000000 0.000000 O
-0.500000 -1.000000 0.000000 O
0.500000 -1.000000 0.000000 O
1.500000 -1.000000 0.000000 O
-0.500000 -2.000000 -1.000000 O
0.000000 -2.000000 -0.500000 O
1.000000 -2.000000 -0.500000 O
0.000000 -1.000000 -0.500000 O
1.000000 -1.000000 -0.500000 O
1.000000 -2.000000 0.500000 O
0.000000 -1.000000 0.500000 O
1.000000 -2.000000 -1.500000 O
0.000000 -1.000000 -1.500000 O
"""

        mcsqs = Mcsqs.structure_from_str(test_string)

        assert mcsqs.composition.formula == "Sr3 Ca5 Mn7 Fe1 O24"
        assert mcsqs.lattice.a == approx(2.2360679775)
        assert mcsqs.lattice.b == approx(2.2360679775)
        assert mcsqs.lattice.c == approx(1.73205080757)

    def test_mcsqs_export(self):
        struct = self.get_structure("SrTiO3")
        struct.replace_species({"Sr2+": {"Sr2+": 0.5, "Ca2+": 0.5}})

        ref_string = """3.905000 0.000000 0.000000
-0.000000 3.905000 0.000000
0.000000 0.000000 3.905000
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
0.500000 0.500000 0.500000 Sr2+=0.5,Ca2+=0.5
0.000000 0.000000 0.000000 Ti4+=1.0
0.000000 0.000000 0.500000 O2-=1.0
0.000000 0.500000 0.000000 O2-=1.0
0.500000 0.000000 0.000000 O2-=1.0"""

        assert Mcsqs(struct).to_str() == ref_string

    def test_mcsqs_cif_nacl(self):
        # CIF file from str2cif (utility distributed with atat)
        struct_from_cif = Structure.from_file(f"{test_dir}/bestsqs_nacl.cif")

        # output file directly from mcsqs
        struct_from_out = Structure.from_file(f"{test_dir}/bestsqs_nacl.out")

        assert struct_from_cif.matches(struct_from_out)
        assert_allclose(
            struct_from_out.lattice.parameters,
            struct_from_cif.lattice.parameters,
            atol=1e-4,
        )

    def test_mcsqs_cif_pzt(self):
        # CIF file from str2cif (utility distributed with atat)
        struct_from_cif = Structure.from_file(f"{test_dir}/bestsqs_pzt.cif")

        # output file directly from mcsqs
        struct_from_out = Structure.from_file(f"{test_dir}/bestsqs_pzt.out")

        assert struct_from_cif.matches(struct_from_out)
        assert_allclose(
            struct_from_out.lattice.parameters,
            struct_from_cif.lattice.parameters,
            atol=1e-4,
        )
