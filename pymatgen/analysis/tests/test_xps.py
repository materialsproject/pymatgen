import unittest

from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.xps import XPS


class XPSTestCase(PymatgenTest):
    def test_from_dos(self):
        v = Vasprun(PymatgenTest.TEST_FILES_DIR / "vasprun.xml.LiF")
        dos = v.complete_dos
        xps = XPS.from_dos(dos)
        self.assertEqual(len(xps), 301)
        xps.smear(0.3)
        self.assertEqual(len(xps), 301)


if __name__ == "__main__":
    unittest.main()
