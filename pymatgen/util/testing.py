"""
Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.
"""
import os
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


    def serialize_with_pickle(self, objects, protocols=None, test_eq=True):
        """
        Test whether the object(s) can be serialized and deserialized with pickle.
        This method tries to serialize the objects with pickle and the protocols
        specified in input. Then it deserializes the pickle format and compares
        the two objects with the __eq__ operator if test_eq == True.

        Args:
            objects:
                Object or list of objects. 
            protocols:
                List of pickle protocols to test.

        Returns:
            Nested list with the objects deserialized with the specified protocols.
        """
        import tempfile
        # Use the python version so that we get the traceback in case of errors 
        import pickle as pickle  
        #import cPickle as pickle

        # Build a list even when we receive a single object.
        got_single_object = False
        if not isinstance(objects, (list, tuple)):
            got_single_object = True
            objects = [objects]

        # By default, all pickle protocols are tested.
        if protocols is None:
            protocols = set([0, 1, 2] + [pickle.HIGHEST_PROTOCOL])

        # This list will contains the object deserialized with the different protocols.
        objects_by_protocol = []

        for protocol in protocols:
            # Serialize and deserialize the object.
            mode = "w" if protocol == 0 else "wb"
            fd, tmpfile = tempfile.mkstemp(text="b" not in mode)

            with open(tmpfile, mode) as fh:
                pickle.dump(objects, fh, protocol=protocol)

            with open(tmpfile, "r") as fh:
                new_objects = pickle.load(fh)

            # Test for equality
            if test_eq:
                for old_obj, new_obj in zip(objects, new_objects):
                    self.assert_equal(old_obj, new_obj)

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(new_objects)

        # Return nested list so that client code can perform additional tests.
        if got_single_object:
            return [o[0] for o in objects_by_protocol]
        else:
            return objects_by_protocol
