from __future__ import division, print_function, unicode_literals

__author__ = 'setten'

import os
from pymatgen.util.testing import PymatgenTest
from pymatgen.core.structure import Structure
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.io.gwwrapper.datastructures import GWSpecs, GWConvergenceData, get_spec
from pymatgen.io.gwwrapper.codeinterfaces import AbinitInterface, NewCodeInterface, VaspInterface, get_code_interface
from pymatgen.io.gwwrapper.GWworkflows import GWG0W0VaspInputSet, SingleAbinitGWWorkFlow
from pymatgen.io.gwwrapper.helpers import refine_structure, clean, load_ps, read_extra_abivars, read_grid_from_file
from pymatgen.io.gwwrapper.helpers import expand
from pymatgen.io.gwwrapper.codeinterfaces import get_all_ecuteps, get_all_nbands, CODE_CLASSES


class GWFunctionsTest(PymatgenTest):

    def test_get_all_ecuteps(self):
        self.assertEqual(set(get_all_ecuteps()), set(['ecuteps', 'ENCUTGW', 'new_code_nbands']))

    def test_get_all_nbands(self):
        self.assertEqual(set(get_all_nbands()), set(['nscf_nbands', 'NBANDS', 'new_code_nbands']))


class GWConstantsTest(PymatgenTest):

    def test_code_classes(self):
        self.assertEqual(CODE_CLASSES, {'VASP': VaspInterface, 'ABINIT': AbinitInterface, 'NEW_CODE': NewCodeInterface})
