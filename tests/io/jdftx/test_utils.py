from __future__ import annotations

import pytest

from pymatgen.io.jdftx.generic_tags import flatten_list
from pymatgen.io.jdftx.jdftxoutfileslice_helpers import get_start_lines


def test_flatten_list():
    assert flatten_list("", [1, 2, 3]) == [1, 2, 3]
    assert flatten_list("", [1, [2, 3]]) == [1, 2, 3]
    assert flatten_list("", [1, [2, [3]]]) == [1, 2, 3]
    assert flatten_list("", [1, [2, [3, [4]]]]) == [1, 2, 3, 4]
    with pytest.raises(TypeError):
        flatten_list("", 1)


def test_get_start_lines():
    with pytest.raises(ValueError, match="Outfile parser fed an empty file."):
        get_start_lines([])
    with pytest.raises(ValueError, match="No JDFTx calculations found in file."):
        get_start_lines(["\n", "\n"])
