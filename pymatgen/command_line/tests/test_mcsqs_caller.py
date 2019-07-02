import unittest

from pymatgen.core.structure import Structure
from pymatgen.command_line.critic2_caller import *
from monty.os.path import which

__author__ = "Handong ling
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__date__ = "June 2019"


@unittest.skipIf(not which('mcsqs'), "mcsqs executable not present")
class Mcsqs_CallerTest(unittest.TestCase):
	def setUp(self):

        self.pzt = Critic2Output(structure, reference_stdout)
        self.pztstrings = np.load(os.path.join(test_dir, "pztist.npy"), allow_pickle=True)
		self.pzt_struc = self.get_structure('Pb2TiZrO6')