# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import tempfile
import numpy.testing.utils as nptu
from io import open
from pathlib import Path
import json

from monty.json import MontyDecoder
from monty.serialization import loadfn
from monty.json import MSONable
from monty.dev import requires

from pymatgen import SETTINGS, MPRester

"""
Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.
"""


class PymatgenTest(unittest.TestCase):
    """
    Extends unittest.TestCase with functions (taken from numpy.testing.utils)
    that support the comparison of arrays.
    """
    _multiprocess_shared_ = True
    MODULE_DIR = Path(__file__).absolute().parent
    STRUCTURES_DIR = MODULE_DIR / "structures"

    """
    Dict for test structures to aid testing.
    """
    TEST_STRUCTURES = {}
    for fn in STRUCTURES_DIR.iterdir():
        TEST_STRUCTURES[fn.name.rsplit(".", 1)[0]] = loadfn(str(fn))

    @classmethod
    def get_structure(cls, name):
        return cls.TEST_STRUCTURES[name].copy()

    @classmethod
    @requires(SETTINGS.get("PMG_MAPI_KEY"), "PMG_MAPI_KEY needs to be set.")
    def get_mp_structure(cls, mpid):
        m = MPRester()
        return m.get_structure_by_material_id(mpid)

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
        Test whether the object(s) can be serialized and deserialized with
        pickle. This method tries to serialize the objects with pickle and the
        protocols specified in input. Then it deserializes the pickle format
        and compares the two objects with the __eq__ operator if
        test_eq == True.

        Args:
            objects: Object or list of objects.
            protocols: List of pickle protocols to test. If protocols is None,
                HIGHEST_PROTOCOL is tested.

        Returns:
            Nested list with the objects deserialized with the specified
            protocols.
        """
        # Use the python version so that we get the traceback in case of errors
        import pickle as pickle
        from pymatgen.util.serialization import pmg_pickle_load, \
            pmg_pickle_dump

        # Build a list even when we receive a single object.
        got_single_object = False
        if not isinstance(objects, (list, tuple)):
            got_single_object = True
            objects = [objects]

        if protocols is None:
            # protocols = set([0, 1, 2] + [pickle.HIGHEST_PROTOCOL])
            protocols = [pickle.HIGHEST_PROTOCOL]

        # This list will contains the object deserialized with the different
        # protocols.
        objects_by_protocol, errors = [], []

        for protocol in protocols:
            # Serialize and deserialize the object.
            mode = "wb"
            fd, tmpfile = tempfile.mkstemp(text="b" not in mode)

            try:
                with open(tmpfile, mode) as fh:
                    pmg_pickle_dump(objects, fh, protocol=protocol)
            except Exception as exc:
                errors.append("pickle.dump with protocol %s raised:\n%s" %
                              (protocol, str(exc)))
                continue

            try:
                with open(tmpfile, "rb") as fh:
                    new_objects = pmg_pickle_load(fh)
            except Exception as exc:
                errors.append("pickle.load with protocol %s raised:\n%s" %
                              (protocol, str(exc)))
                continue

            # Test for equality
            if test_eq:
                for old_obj, new_obj in zip(objects, new_objects):
                    self.assertEqual(old_obj, new_obj)

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(new_objects)

        if errors:
            raise ValueError("\n".join(errors))

        # Return nested list so that client code can perform additional tests.
        if got_single_object:
            return [o[0] for o in objects_by_protocol]
        else:
            return objects_by_protocol

    def tmpfile_write(self, string):
        """
        Write string to a temporary file. Returns the name of the temporary
        file.
        """
        fd, tmpfile = tempfile.mkstemp(text=True)

        with open(tmpfile, "w") as fh:
            fh.write(string)

        return tmpfile

    def assertMSONable(self, obj, test_if_subclass=True):
        """
        Tests if obj is MSONable and tries to verify whether the contract is
        fulfilled.

        By default, the method tests whether obj is an instance of MSONable.
        This check can be deactivated by setting test_if_subclass to False.
        """
        if test_if_subclass:
            self.assertIsInstance(obj, MSONable)
        self.assertDictEqual(obj.as_dict(), obj.__class__.from_dict(
            obj.as_dict()).as_dict())
        json.loads(obj.to_json(), cls=MontyDecoder)
