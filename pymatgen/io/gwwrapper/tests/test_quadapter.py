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
from pymatgen.io.abinitio.tasks import TaskManager

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",'test_files')

from pymatgen.io.abinitio.flows import AbinitFlow


class GWSpecTest(PymatgenTest):
    def test_fixes(self):
        flow = AbinitFlow(workdir=test_dir, manager=TaskManager.from_file(os.path.join(test_dir, "taskmanager.yml")))
        inp = {}
        flow.register_task(input=inp)
        flow.allocate()



