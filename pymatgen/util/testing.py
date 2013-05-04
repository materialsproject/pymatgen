"""
Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.
"""
import unittest
import numpy.testing.utils as nptu


class PymatgenTest(unittest.TestCase):
    """
    Extends unittest.TestCase with functions (taken from numpy.testing.utils)
    that support the comparison of ndarrays.
    """

    @staticmethod
    def assert_almost_equal(actual, desired, decimal=7, err_msg='',
                            verbose=True):
        return nptu.assert_almost_equal(actual, desired, decimal, err_msg,
                                        verbose)

    @staticmethod
    def assert_equal(actual, desired, err_msg='', verbose=True):
        return nptu.assert_equal(actual, desired, err_msg=err_msg,
                               verbose=verbose)
