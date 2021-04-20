import unittest

from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.xps import XPS


class XPSTestCase(PymatgenTest):
    def test_from_dos(self):
        # mpr = MPRester()
        # dos_SnO2 = mpr.get_dos_by_material_id("mp-856")
        v = Vasprun(PymatgenTest.TEST_FILES_DIR / "vasprun.xml.LiF")
        dos = v.complete_dos
        xps = XPS.from_dos(dos, 0.5)
        self.assertEqual(len(xps), 301)


if __name__ == "__main__":
    unittest.main()
