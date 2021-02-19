import os
import unittest

from monty.serialization import loadfn

from pymatgen.util.testing import PymatgenTest


class IRDielectricTensorTest(PymatgenTest):
    def setUp(self):
        self.ir_spectra = loadfn(os.path.join(PymatgenTest.TEST_FILES_DIR, "ir_spectra_mp-991652_DDB.json"))

    def test_basic(self):
        self.ir_spectra.write_json("test.json")
        ir_spectra = loadfn("test.json")
        irdict = ir_spectra.as_dict()
        ir_spectra.from_dict(irdict)
        ir_spectra.plot(show=False)

    def tearDown(self):
        if os.path.isfile("test.json"):
            os.remove("test.json")


if __name__ == "__main__":
    unittest.main()
