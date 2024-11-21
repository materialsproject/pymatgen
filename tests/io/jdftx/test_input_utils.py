from __future__ import annotations

import pytest

from pymatgen.io.jdftx.generic_tags import _flatten_list


def test_flatten_list():
    assert _flatten_list("", [1, 2, 3]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, 3]]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, [3]]]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, [3, [4]]]]) == [1, 2, 3, 4]
    with pytest.raises(TypeError):
        _flatten_list("", 1)
