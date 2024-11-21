from __future__ import annotations

import pytest

from pymatgen.io.jdftx._input_utils import flatten_list


def test_flatten_list():
    assert flatten_list("", [1, 2, 3]) == [1, 2, 3]
    assert flatten_list("", [1, [2, 3]]) == [1, 2, 3]
    assert flatten_list("", [1, [2, [3]]]) == [1, 2, 3]
    assert flatten_list("", [1, [2, [3, [4]]]]) == [1, 2, 3, 4]
    with pytest.raises(TypeError):
        flatten_list("", 1)
