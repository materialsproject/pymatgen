from __future__ import annotations

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
