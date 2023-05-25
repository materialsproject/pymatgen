"""
Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from monty.json import MontyDecoder, MSONable
from monty.serialization import loadfn
from numpy.testing import assert_allclose

from pymatgen.core import SETTINGS, Structure


class PymatgenTest(unittest.TestCase):
    """
    Extends unittest.TestCase with functions (taken from numpy.testing.utils)
    that support the comparison of arrays.
    """

    _multiprocess_shared_ = True
    MODULE_DIR = Path(__file__).absolute().parent
    STRUCTURES_DIR = MODULE_DIR / "structures"
    try:
        TEST_FILES_DIR = Path(SETTINGS["PMG_TEST_FILES_DIR"])
    except KeyError:
        import warnings

        warnings.warn(
            "It is recommended that you set the PMG_TEST_FILES_DIR environment variable explicitly. "
            "Now using a fallback location based on relative path from this module."
        )
        TEST_FILES_DIR = MODULE_DIR / ".." / ".." / "test_files"

    TEST_STRUCTURES = {}  # Dict for test structures to aid testing.
    for fn in STRUCTURES_DIR.iterdir():
        TEST_STRUCTURES[fn.name.rsplit(".", 1)[0]] = loadfn(str(fn))

    @classmethod
    def get_structure(cls, name: str) -> Structure:
        """
        Get a structure from the template directories.

        Args:
            name (str): Name of structure file.

        Returns:
            Structure
        """
        return cls.TEST_STRUCTURES[name].copy()

    @staticmethod
    def assert_all_close(actual, desired, decimal=7, err_msg="", verbose=True):
        """
        Tests if two arrays are almost equal up to some relative or absolute tolerance.
        """
        # TODO (janosh): replace the decimal kwarg with assert_allclose() atol and rtol kwargs
        return assert_allclose(actual, desired, atol=10**-decimal, err_msg=err_msg, verbose=verbose)

    @staticmethod
    def assert_str_content_equal(actual, desired, err_msg="", verbose=True):
        """
        Tests if two strings are equal, ignoring things like trailing spaces, etc.
        """
        lines1 = actual.split("\n")
        lines2 = desired.split("\n")
        if len(lines1) != len(lines2):
            return False
        failed = []
        for l1, l2 in zip(lines1, lines2):
            if l1.strip() != l2.strip():
                failed.append(f"{l1} != {l2}")
        return len(failed) == 0

    def serialize_with_pickle(self, objects, protocols=None, test_eq=True):
        """
        Test whether the object(s) can be serialized and deserialized with
        pickle. This method tries to serialize the objects with pickle and the
        protocols specified in input. Then it deserializes the pickle format
        and compares the two objects with the __eq__ operator if
        test_eq is True.

        Args:
            objects: Object or list of objects.
            protocols: List of pickle protocols to test. If protocols is None,
                HIGHEST_PROTOCOL is tested.
            test_eq: If True, the deserialized object is compared with the
                original object using the __eq__ method.

        Returns:
            Nested list with the objects deserialized with the specified
            protocols.
        """
        # Use the python version so that we get the traceback in case of errors
        import pickle

        # Build a list even when we receive a single object.
        got_single_object = False
        if not isinstance(objects, (list, tuple)):
            got_single_object = True
            objects = [objects]

        if protocols is None:
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
                    pickle.dump(objects, fh, protocol=protocol)
            except Exception as exc:
                errors.append(f"pickle.dump with protocol {protocol} raised:\n{exc}")
                continue

            try:
                with open(tmpfile, "rb") as fh:
                    new_objects = pickle.load(fh)
            except Exception as exc:
                errors.append(f"pickle.load with protocol {protocol} raised:\n{exc}")
                continue

            # Test for equality
            if test_eq:
                for old_obj, new_obj in zip(objects, new_objects):
                    assert old_obj == new_obj

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(new_objects)

        if errors:
            raise ValueError("\n".join(errors))

        # Return nested list so that client code can perform additional tests.
        if got_single_object:
            return [o[0] for o in objects_by_protocol]
        return objects_by_protocol

    def assert_msonable(self, obj, test_if_subclass=True):
        """
        Test if obj is MSONable and verify the contract is fulfilled.

        By default, the method tests whether obj is an instance of MSONable.
        This check can be deactivated by setting test_if_subclass=False.
        """
        if test_if_subclass:
            assert isinstance(obj, MSONable)
        assert obj.as_dict() == obj.__class__.from_dict(obj.as_dict()).as_dict()
        json.loads(obj.to_json(), cls=MontyDecoder)
