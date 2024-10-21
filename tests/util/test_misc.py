from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pymatgen.util.misc import is_np_dict_equal


class TestIsNpDictEqual:
    def test_different_keys(self):
        """Test two dicts with different keys."""
        dict1 = {"a": np.array([1, 2, 3])}
        dict2 = {"a": np.array([1, 2, 3]), "b": "hello"}
        equal = is_np_dict_equal(dict1, dict2)
        # make sure it's not a np.bool
        assert isinstance(equal, bool)
        assert not equal

    def test_both_list(self):
        """Test two dicts where both have lists as values."""
        dict1 = {"a": [1, 2, 3]}
        dict2 = {"a": [1, 2, 3]}
        assert is_np_dict_equal(dict1, dict2)

    def test_both_np_array(self):
        """Test two dicts where both have NumPy arrays as values."""
        dict1 = {"a": np.array([1, 2, 3])}
        dict2 = {"a": np.array([1, 2, 3])}
        assert is_np_dict_equal(dict1, dict2)

    def test_one_np_one_list(self):
        """Test two dicts where one has a NumPy array and the other has a list."""
        dict1 = {"a": np.array([1, 2, 3])}
        dict2 = {"a": [1, 2, 3]}
        assert is_np_dict_equal(dict1, dict2)

    def test_nested_arrays(self):
        """Test two dicts with deeper nested arrays."""
        dict1 = {"a": np.array([[1, 2], [3, 4]])}
        dict2 = {"a": np.array([[1, 2], [3, 4]])}
        assert is_np_dict_equal(dict1, dict2)

        dict3 = {"a": np.array([[1, 2], [3, 5]])}
        assert not is_np_dict_equal(dict1, dict3)

    def test_diff_dtype(self):
        """Make sure it also works for other data types as value."""

        @dataclass
        class CustomClass:
            name: str
            value: int

        # Test with bool values
        dict1 = {"a": True}
        dict2 = {"a": True}
        assert is_np_dict_equal(dict1, dict2)

        dict3 = {"a": False}
        assert not is_np_dict_equal(dict1, dict3)

        # Test with string values
        dict4 = {"a": "hello"}
        dict5 = {"a": "hello"}
        assert is_np_dict_equal(dict4, dict5)

        dict6 = {"a": "world"}
        assert not is_np_dict_equal(dict4, dict6)

        # Test with a custom data class
        dict7 = {"a": CustomClass(name="test", value=1)}
        dict8 = {"a": CustomClass(name="test", value=1)}
        assert is_np_dict_equal(dict7, dict8)

        dict9 = {"a": CustomClass(name="test", value=2)}
        assert not is_np_dict_equal(dict7, dict9)

        # Test with None
        dict10 = {"a": None}
        dict11 = {"a": None}
        assert is_np_dict_equal(dict10, dict11)

        dict12 = {"a": None}
        dict13 = {"a": "non-none"}
        assert not is_np_dict_equal(dict12, dict13)

        # Test with nested complex lists
        dict14 = {"a": [[1, 2], ["hello", 3.0]]}
        dict15 = {"a": [np.array([1, 2]), ["hello", 3.0]]}
        assert is_np_dict_equal(dict14, dict15)

        dict16 = {"a": [[1, 2], ["world", 3.0]]}
        assert not is_np_dict_equal(dict14, dict16)

        # Test with unhashable dicts
        dict17 = {"a": {"key1": "val1", "key2": "val2"}}
        dict18 = {"a": {"key1": "val1", "key2": "val2"}}
        assert is_np_dict_equal(dict17, dict18)

        # Test with unhashable sets
        dict19 = {"a": {1, 2, 3}}
        dict20 = {"a": {1, 2, 3}}
        assert is_np_dict_equal(dict19, dict20)
