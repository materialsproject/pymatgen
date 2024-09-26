from __future__ import annotations

import pytest

from pymatgen.io.jdftx.generic_tags import flatten_list
from pymatgen.io.jdftx.jdftxoutfileslice_helpers import find_first_range_key, get_start_lines


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


def test_find_first_range_key():
    out1 = find_first_range_key("barbie", ["barbie"])
    assert len(out1) == 1
    assert out1[0] == 0
    out1 = find_first_range_key("barbie", ["# barbie"])
    assert len(out1) == 0
    out1 = find_first_range_key("barbie", ["# barbie"], skip_pound=True)
    assert len(out1) == 1
    assert out1[0] == 0
    out1 = find_first_range_key("barbie", ["barbie", "barbie", "barbie"])
    assert len(out1) == 3
    out1 = find_first_range_key("barbie", ["barbie", "ken", "barbie"])
    assert len(out1) == 2
    out1 = find_first_range_key("barbie", ["barbie", "ken", "barbie"], startline=1)
    assert len(out1) == 1
