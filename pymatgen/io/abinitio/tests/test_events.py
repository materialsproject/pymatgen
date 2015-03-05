# coding: utf-8
from __future__ import unicode_literals, division, print_function

#import os
#import tempfile
#import shutil

from pymatgen.util.testing import PymatgenTest
#from monty.functools import lazy_property
#from pymatgen.core.lattice import Lattice
#from pymatgen.core.structure import Structure
from pymatgen.io.abinitio import events
#from pymatgen.io.abinitio.flows import *
#from pymatgen.io.abinitio.tasks import *
#from pymatgen.io.abinitio.pseudos import Pseudo

#_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", 'test_files')
#
#
#def ref_file(filename):
#    return os.path.join(_test_dir, filename)


class EventHandlersTest(PymatgenTest):
    def test_api(self):
        events.autodoc_event_handlers()


if __name__ == '__main__':
    import unittest
    unittest.main()
