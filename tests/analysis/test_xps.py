from __future__ import annotations

from pymatgen.analysis.xps import XPS
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class XPSTestCase(PymatgenTest):
    def test_from_dos(self):
        vasp_run = Vasprun(f"{TEST_FILES_DIR}/vasprun.xml.LiF")
        dos = vasp_run.complete_dos
        xps = XPS.from_dos(dos)
        assert len(xps) == 301
        xps.smear(0.3)
        assert len(xps) == 301
