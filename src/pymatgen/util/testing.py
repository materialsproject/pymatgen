"""This module implements testing utilities for materials science codes.

While the primary use is within pymatgen, the functionality is meant to be useful for external materials science
codes as well. For instance, obtaining example crystal structures to perform tests, specialized assert methods for
materials science, etc.
"""

from __future__ import annotations

import json
import pickle  # use pickle over cPickle to get traceback in case of errors
import string
from pathlib import Path
from typing import TYPE_CHECKING
from unittest import TestCase

import pytest
from monty.dev import deprecated
from monty.json import MontyDecoder, MontyEncoder, MSONable
from monty.serialization import loadfn

from pymatgen.core import ROOT, SETTINGS

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, ClassVar

    from pymatgen.core import Structure
    from pymatgen.util.typing import PathLike

_MODULE_DIR: Path = Path(__file__).absolute().parent

STRUCTURES_DIR: Path = _MODULE_DIR / "structures"

TEST_FILES_DIR: Path = Path(SETTINGS.get("PMG_TEST_FILES_DIR", f"{ROOT}/../tests/files"))
VASP_IN_DIR: str = f"{TEST_FILES_DIR}/io/vasp/inputs"
VASP_OUT_DIR: str = f"{TEST_FILES_DIR}/io/vasp/outputs"

# Fake POTCARs have original header information, meaning properties like number of electrons,
# nuclear charge, core radii, etc. are unchanged (important for testing) while values of the and
# pseudopotential kinetic energy corrections are scrambled to avoid VASP copyright infringement
FAKE_POTCAR_DIR: str = f"{VASP_IN_DIR}/fake_potcars"


class MatSciTest:
    """`pytest` based test framework extended to facilitate testing with
    the following methods:
    - tmp_path (attribute): Temporary directory.
    - get_structure: Load a Structure from `util.structures` with its name.
    - assert_str_content_equal: Check if two strings are equal (ignore whitespaces).
    - serialize_with_pickle: Test whether object(s) can be (de)serialized with pickle.
    - assert_msonable: Test if obj is MSONable and return its serialized string.
    """

    # dict of lazily-loaded test structures (initialized to None)
    TEST_STRUCTURES: ClassVar[dict[PathLike, Structure | None]] = dict.fromkeys(STRUCTURES_DIR.glob("*"))

    @pytest.fixture(autouse=True)
    def _tmp_dir(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Make all tests run a in a temporary directory accessible via self.tmp_path.

        References:
            https://docs.pytest.org/en/stable/how-to/tmp_path.html
        """
        monkeypatch.chdir(tmp_path)  # change to temporary directory
        self.tmp_path = tmp_path

    @staticmethod
    def assert_msonable(obj: Any, test_is_subclass: bool = True) -> str:
        """Test if an object is MSONable and verify the contract is fulfilled,
        and return the serialized object.

        By default, the method tests whether obj is an instance of MSONable.
        This check can be deactivated by setting `test_is_subclass` to False.

        Args:
            obj (Any): The object to be checked.
            test_is_subclass (bool): Check if object is an instance of MSONable
                or its subclasses.

        Returns:
            str: Serialized object.
        """
        obj_name = obj.__class__.__name__

        # Check if is an instance of MONable (or its subclasses)
        if test_is_subclass and not isinstance(obj, MSONable):
            raise TypeError(f"{obj_name} object is not MSONable")

        # Check if the object can be accurately reconstructed from its dict representation
        if obj.as_dict() != type(obj).from_dict(obj.as_dict()).as_dict():
            raise ValueError(f"{obj_name} object could not be reconstructed accurately from its dict representation.")

        # Verify that the deserialized object's class is a subclass of the original object's class
        json_str = json.dumps(obj.as_dict(), cls=MontyEncoder)
        round_trip = json.loads(json_str, cls=MontyDecoder)
        if not issubclass(type(round_trip), type(obj)):
            raise TypeError(f"The reconstructed {round_trip.__class__.__name__} object is not a subclass of {obj_name}")
        return json_str

    @staticmethod
    def assert_str_content_equal(actual: str, expected: str) -> None:
        """Test if two strings are equal, ignoring whitespaces.

        Args:
            actual (str): The string to be checked.
            expected (str): The reference string.

        Raises:
            AssertionError: When two strings are not equal.
        """
        strip_whitespace = {ord(c): None for c in string.whitespace}
        if actual.translate(strip_whitespace) != expected.translate(strip_whitespace):
            raise AssertionError(
                "Strings are not equal (whitespaces ignored):\n"
                f"{' Actual '.center(50, '=')}\n"
                f"{actual}\n"
                f"{' Expected '.center(50, '=')}\n"
                f"{expected}\n"
            )

    @classmethod
    def get_structure(cls, name: str) -> Structure:
        """
        Load a structure from `pymatgen.util.structures`.

        Args:
            name (str): Name of the structure file, for example "LiFePO4".

        Returns:
            Structure
        """
        try:
            struct = cls.TEST_STRUCTURES.get(name) or loadfn(f"{STRUCTURES_DIR}/{name}.json")
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"structure for {name} doesn't exist") from exc

        cls.TEST_STRUCTURES[name] = struct

        return struct.copy()

    def serialize_with_pickle(
        self,
        objects: Any,
        protocols: Sequence[int] | None = None,
        test_eq: bool = True,
    ) -> list:
        """Test whether the object(s) can be serialized and deserialized with
        `pickle`. This method tries to serialize the objects with `pickle` and the
        protocols specified in input. Then it deserializes the pickled format
        and compares the two objects with the `==` operator if `test_eq`.

        Args:
            objects (Any): Object or list of objects.
            protocols (Sequence[int]): List of pickle protocols to test.
                If protocols is None, HIGHEST_PROTOCOL is tested.
            test_eq (bool): If True, the deserialized object is compared
                with the original object using the `__eq__` method.

        Returns:
            list[Any]: Objects deserialized with the specified protocols.
        """
        # Build a list even when we receive a single object.
        got_single_object = False
        if not isinstance(objects, list | tuple):
            got_single_object = True
            objects = [objects]

        protocols = protocols or [pickle.HIGHEST_PROTOCOL]

        # This list will contain the objects deserialized with the different protocols.
        objects_by_protocol, errors = [], []

        for protocol in protocols:
            # Serialize and deserialize the object.
            tmpfile = self.tmp_path / f"tempfile_{protocol}.pkl"

            try:
                with open(tmpfile, "wb") as file:
                    pickle.dump(objects, file, protocol=protocol)
            except Exception as exc:
                errors.append(f"pickle.dump with {protocol=} raised:\n{exc}")
                continue

            try:
                with open(tmpfile, "rb") as file:
                    unpickled_objs = pickle.load(file)  # noqa: S301
            except Exception as exc:
                errors.append(f"pickle.load with {protocol=} raised:\n{exc}")
                continue

            # Test for equality
            if test_eq:
                for orig, unpickled in zip(objects, unpickled_objs, strict=True):
                    if orig != unpickled:
                        raise ValueError(
                            f"Unpickled and original objects are unequal for {protocol=}\n{orig=}\n{unpickled=}"
                        )

            # Save the deserialized objects and test for equality.
            objects_by_protocol.append(unpickled_objs)

        if errors:
            raise ValueError("\n".join(errors))

        # Return nested list so that client code can perform additional tests.
        if got_single_object:
            return [o[0] for o in objects_by_protocol]
        return objects_by_protocol


@deprecated(MatSciTest, deadline=(2026, 1, 1))
class PymatgenTest(TestCase, MatSciTest):
    """Extends unittest.TestCase with several assert methods for array and str comparison.

    Deprecated: please use `MatSciTest` instead (migrate from `unittest` to `pytest`).
    """
