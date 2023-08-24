from __future__ import annotations

import os

from monty.serialization import loadfn

from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestIRDielectricTensor(PymatgenTest):
    def setUp(self):
        self.ir_spectra = loadfn(f"{TEST_FILES_DIR}/ir_spectra_mp-991652_DDB.json")

    def test_basic(self):
        self.ir_spectra.write_json("test.json")
        ir_spectra = loadfn("test.json")
        irdict = ir_spectra.as_dict()
        ir_spectra.from_dict(irdict)
        ir_spectra.plot(show=False)

    def tearDown(self):
        if os.path.isfile("test.json"):
            os.remove("test.json")
