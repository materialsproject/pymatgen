# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import warnings

import os
import unittest

from monty.serialization import loadfn

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

class HeisenbergMapperTest(unittest.TestCase):
    """Todo:
        Put json files of magnetic orderings and energies in
        pymatgen/test_files with variable number of unique magnetic
        sites, sublattices, FM, FiM, and AFM ground state, etc.
    """

    def setUp(self):
        pass

    def tearDown(self):
        warnings.simplefilter('default')

    def test_graphs(self):
        pass

    def test_sites(self):
        pass

    def test_nn_interactions(self):
        pass

    def test_exchange_matrix(self):
        pass

    def test_exchange_params(self):
        pass

    def test_mean_field(self):
        pass





if __name__ == '__main__':
    unittest.main()