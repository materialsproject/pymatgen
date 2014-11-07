from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import os
import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.io.gwwrapper.datastructures import GWSpecs, GWConvergenceData, get_spec
from pymatgen.io.gwwrapper.codeinterfaces import AbinitInterface, AbstractCodeInterface, VaspInterface, get_code_interface
from pymatgen.io.gwwrapper.GWvaspinputsets import GWDFTDiagVaspInputSet, GWG0W0VaspInputSet, GWscDFTPrepVaspInputSet
from pymatgen.io.gwwrapper.GWvaspinputsets import SingleVaspGWWork
from pymatgen.io.gwwrapper.tests.test_helpers import structure

#test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",'test_files')

from pymatgen.io.vaspio.vasp_input import get_potcar_dir
POTCAR_DIR = get_potcar_dir()

class GWSpecTest(PymatgenTest):
    def test_GWspect(self):
        spec = GWSpecs()
        self.assertIsInstance(spec, GWSpecs)
        self.assertEqual(spec.get_code(), 'ABINIT')
        self.assertIsInstance(spec.code_interface, AbinitInterface)
        spec.data['code'] = 'VASP'
        spec.update_code_interface()
        self.assertEqual(spec.get_code(), 'VASP')
        self.assertIsInstance(spec.code_interface, VaspInterface)

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
        spec = GWSpecs()
        self.assertIsInstance(interface, VaspInterface)
        self.assertEqual(len(interface.conv_pars), 3)
        self.assertEqual(len(interface.supported_methods), 4)
        #self.assertEqual(interface.get_conv_res_test(spec_data=spec.data, structure=structure), {})

    def test_AbinitInterface(self):
        interface = get_code_interface('ABINIT')
        self.assertIsInstance(interface, AbinitInterface)
        self.assertEqual(len(interface.conv_pars), 3)
        self.assertEqual(len(interface.supported_methods), 2)


class GWVaspInputSetTests(PymatgenTest):

    def setUp(self):
        self.structure = structure
        self.spec = GWSpecs()
        self.spec.data['code'] = 'VASP'
        self.spec.update_code_interface()

    def test_GWscDFTPrepVaspInputSet(self):
        inpset = GWscDFTPrepVaspInputSet(structure=self.structure, spec=self.spec)
        self.assertIsInstance(inpset, GWscDFTPrepVaspInputSet)
        self.assertEqual(inpset.convs, {})

    @unittest.skipIf(POTCAR_DIR is None, "POTCAR dir is None")
    def test_GWDFTDiagVaspInputSet(self):
        self.maxDiff = None
        inpset = GWDFTDiagVaspInputSet(structure=self.structure, spec=self.spec)
        self.assertIsInstance(inpset, GWDFTDiagVaspInputSet)
        self.assertEqual(inpset.convs,
                         {u'NBANDS': {u'test_range': (10, 20, 30, 40, 50, 60, 70), u'control': u'gap', u'method': u'set_nbands'}})
        self.assertEqual(inpset.incar_settings, {u'ALGO': u'Exact', u'EDIFF': 1e-10, u'IBRION': -1, u'ICHARG': 1,
                                                 u'ISMEAR': -5, u'ISPIN': 1, u'LOPTICS': u'TRUE', u'LORBIT': 11,
                                                 u'LREAL': u'AUTO', u'LWAVE': True, u'MAGMOM': {u'Co': 5, u'Cr': 5,
                                                 u'Fe': 5, u'Mn': 5, u'Mn3+': 4, u'Mn4+': 3, u'Mo': 5, u'Ni': 5,
                                                 u'V': 5, u'W': 5}, u'NBANDS': 221, u'NELM': 1, u'NPAR': 13,
                                                 u'PREC': u'Medium', u'SIGMA': 0.01})

    @unittest.skipIf(POTCAR_DIR is None, "POTCAR dir is None")
    def test_GWG0W0VaspInputSet(self):
        inpset = GWG0W0VaspInputSet(structure=self.structure, spec=self.spec)
        self.assertIsInstance(inpset, GWG0W0VaspInputSet)
        self.assertEqual(inpset.convs,
                         {u'ENCUTGW': {u'test_range': (200, 400, 600, 800), u'control': u'gap', u'method': u'incar_settings'}})

    def test_SingleVaspGWWork(self):
        work = SingleVaspGWWork(structure=self.structure, spec=self.spec, job='prep')
        self.assertIsInstance(work, SingleVaspGWWork)
