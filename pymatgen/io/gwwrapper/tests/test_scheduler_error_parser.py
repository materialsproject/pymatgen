#!/usr/bin/env python

from __future__ import division
from __future__ import unicode_literals, print_function

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Mar 24, 2014"

import unittest

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.qadapters import *
from pymatgen.io.abinitio.qadapters import AbstractQueueAdapter


class QueueErrorParseTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_something(self):
        pass

if __name__ == '__main__':
    unittest.main()
