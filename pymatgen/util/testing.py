"""Common test support for pymatgen test scripts.

This single module should provide all the common functionality for pymatgen
tests in a single location, so that test scripts can just import it and work
right away.
"""

from __future__ import annotations

import json
import string
import unittest
from pathlib import Path
from typing import TYPE_CHECKING, Any, ClassVar

import pytest
from monty.json import MontyDecoder, MSONable
from monty.serialization import loadfn

from pymatgen.core import SETTINGS, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence

MODULE_DIR = Path(__file__).absolute().parent
STRUCTURES_DIR = MODULE_DIR / "structures"
TEST_FILES_DIR = Path(SETTINGS.get("PMG_TEST_FILES_DIR", MODULE_DIR / ".." / ".." / "tests" / "files"))


class PymatgenTest(unittest.TestCase):
    """Extends unittest.TestCase with several assert methods for array and str comparison."""

    _multiprocess_shared_ = True

    # dict of lazily-loaded test structures (initialized to None)
    TEST_STRUCTURES: ClassVar[dict[str | Path, Structure | None]] = {key: None for key in STRUCTURES_DIR.glob("*")}

    @pytest.fixture(autouse=True)  # make all tests run a in a temporary directory accessible via self.tmp_path
    def _tmp_dir(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        # https://pytest.org/en/latest/how-to/unittest.html#using-autouse-fixtures-and-accessing-other-fixtures
        monkeypatch.chdir(tmp_path)  # change to pytest-provided temporary directory
        self.tmp_path = tmp_path

    @classmethod
    def get_structure(cls, name: str) -> Structure:
        """
        Lazily load a structure from pymatgen/util/testing/structures.

        Args:
            name (str): Name of structure file.

        Returns:
            Structure
        """
        struct = cls.TEST_STRUCTURES.get(name) or loadfn(f"{STRUCTURES_DIR}/{name}.json")
        cls.TEST_STRUCTURES[name] = struct
        return struct.copy()

    @staticmethod
    def assert_str_content_equal(actual, expected):
        """Tests if two strings are equal, ignoring things like trailing spaces, etc."""
        strip_whitespace = {ord(c): None for c in string.whitespace}
        return actual.translate(strip_whitespace) == expected.translate(strip_whitespace)

    def serialize_with_pickle(self, objects: Any, protocols: Sequence[int] | None = None, test_eq: bool = True):
        """Test whether the object(s) can be serialized and deserialized with
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
        # use pickle, not cPickle so that we get the traceback in case of errors
        import pickle

        # Build a list even when we receive a single object.
        got_single_object = False
        if not isinstance(objects, (list, tuple)):
            got_single_object = True
            objects = [objects]

        protocols = protocols or [pickle.HIGHEST_PROTOCOL]

        # This list will contain the objects deserialized with the different protocols.
        objects_by_protocol, errors = [], []

        for protocol in protocols:
            # Serialize and deserialize the object.
            tmpfile = self.tmp_path / f"tempfile_{protocol}.pkl"

            try:
                with open(tmpfile, "wb") as fh:
                    pickle.dump(objects, fh, protocol=protocol)
            except Exception as exc:
                errors.append(f"pickle.dump with {protocol=} raised:\n{exc}")
                continue

            try:
                with open(tmpfile, "rb") as fh:
                    unpickled_objs = pickle.load(fh)
            except Exception as exc:
                errors.append(f"pickle.load with {protocol=} raised:\n{exc}")
                continue

            # Test for equality
            if test_eq:
                for orig, unpickled in zip(objects, unpickled_objs):
                    assert (
                        orig == unpickled
                    ), f"Unpickled and original objects are unequal for {protocol=}\n{orig=}\n{unpickled=}"

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(unpickled_objs)

        if errors:
            raise ValueError("\n".join(errors))

        # Return nested list so that client code can perform additional tests.
        if got_single_object:
            return [o[0] for o in objects_by_protocol]
        return objects_by_protocol

    def assert_msonable(self, obj, test_is_subclass=True):
        """Test if obj is MSONable and verify the contract is fulfilled.

        By default, the method tests whether obj is an instance of MSONable.
        This check can be deactivated by setting test_is_subclass=False.
        """
        if test_is_subclass:
            assert isinstance(obj, MSONable)
        assert obj.as_dict() == obj.__class__.from_dict(obj.as_dict()).as_dict()
        json.loads(obj.to_json(), cls=MontyDecoder)
