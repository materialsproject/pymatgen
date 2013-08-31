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
    that support the comparison of arrays.
    """

    @staticmethod
    def assert_almost_equal(actual, desired, decimal=7, err_msg='',
                            verbose=True):
        """
        Alternative naming for assertArrayAlmostEqual.
        """
        return PymatgenTest.assertArrayAlmostEqual(
            actual, desired, decimal, err_msg, verbose)

    @staticmethod
    def assert_equal(actual, desired, err_msg='', verbose=True):
        """
        Alternative naming for assertArrayEqual.
        """
        return PymatgenTest.assertArrayEqual(actual, desired,
                                             err_msg=err_msg, verbose=verbose)

    @staticmethod
    def assertArrayAlmostEqual(actual, desired, decimal=7, err_msg='',
                               verbose=True):
        """
        Tests if two arrays are almost equal to a tolerance. The CamelCase
        naming is so that it is consistent with standard unittest methods.
        """
        return nptu.assert_almost_equal(actual, desired, decimal, err_msg,
                                        verbose)

    @staticmethod
    def assertArrayEqual(actual, desired, err_msg='', verbose=True):
        """
        Tests if two arrays are equal. The CamelCase naming is so that it is
         consistent with standard unittest methods.
        """
        return nptu.assert_equal(actual, desired, err_msg=err_msg,
                                 verbose=verbose)
