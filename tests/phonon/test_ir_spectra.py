from __future__ import annotations

from monty.serialization import loadfn

from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestIRDielectricTensor(PymatgenTest):
    def setUp(self):
        self.ir_spectra = loadfn(f"{TEST_FILES_DIR}/phonon/dos/ir_spectra_mp-991652_DDB.json")

    def test_basic(self):
        self.ir_spectra.write_json(f"{self.tmp_path}/test.json")
        ir_spectra = loadfn(f"{self.tmp_path}/test.json")
        irdict = ir_spectra.as_dict()
        ir_spectra.from_dict(irdict)
        ir_spectra.plot(show=False)
