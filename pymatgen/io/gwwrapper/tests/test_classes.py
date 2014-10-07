from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import os
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.io.gwwrapper.datastructures import GWSpecs, GWConvergenceData, get_spec
from pymatgen.io.gwwrapper.codeinterfaces import AbinitInterface, AbstractCodeInterface, VaspInterface, get_code_interface

from pymatgen.io.gwwrapper.tests.test_helpers import structure

#test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",'test_files')


class GWSpecTest(PymatgenTest):
    def test_GWspect(self):
        spec = GWSpecs()
        self.assertIsInstance(spec, GWSpecs)
        self.assertEqual(spec.get_code(), 'ABINIT')
        self.assertIsInstance(spec.code_interface, AbinitInterface)

    def test_GWspect_test(self):
        spec = GWSpecs()
        spec.test()
        self.assertEqual(len(spec.warnings), 0)
        self.assertEqual(len(spec.errors), 0)

    def test_GWget_spec(self):
        spec = get_spec('GW')
        self.assertIsInstance(spec, GWSpecs)


class GWConvergenceDataTest(PymatgenTest):
    def test_GWConvergenceData(self):
        spec = GWSpecs()
        self.assertIsInstance(structure, Structure)
        structure.item = 'mp-149'
        conv_data = GWConvergenceData(spec=spec, structure=structure)
        self.assertIsInstance(conv_data, GWConvergenceData)


class GWTestCodeInterfaces(PymatgenTest):
    def test_VaspInterface(self):
        interface = get_code_interface('VASP')
        self.assertIsInstance(interface, VaspInterface)
        self.assertEqual(len(interface.conv_pars), 3)

    def test_AbinitInterface(self):
        interface = get_code_interface('ABINIT')
        self.assertIsInstance(interface, AbinitInterface)
        self.assertEqual(len(interface.conv_pars), 3)


