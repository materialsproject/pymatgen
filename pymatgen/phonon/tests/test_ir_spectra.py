from __future__ import unicode_literals

import unittest
import os
import json
from io import open

from pymatgen.phonon.bandstructure import PhononBandStructure, PhononBandStructureSymmLine
from pymatgen.util.testing import PymatgenTest
from monty.serialization import loadfn


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')

class IRDielectricTensorGeneratorTest(PymatgenTest):

    def setUp(self):
        self.ir_spectra_generator = loadfn(os.path.join(test_dir,'ir_spectra_generator_mp-991652_DDB.json'))

    def test_basic(self):

        self.ir_spectra_generator.write_json('test.json')
        ir_spectra_generator = loadfn('test.json')
        irdict = ir_spectra_generator.as_dict()
        ir_spectra_generator.from_dict(irdict)
        ir_spectra = ir_spectra_generator.get_ir_spectra()
        ir_spectra.plot(show=False)

    def tearDown(self):
        if os.path.isfile('test.json'):
            os.remove('test.json')

if __name__ == '__main__':
    unittest.main()
