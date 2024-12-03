from __future__ import annotations

import pytest

from pymatgen.io.jdftx._output_utils import find_first_range_key, get_start_lines
from pymatgen.io.jdftx.joutstructures import _get_joutstructures_start_idx


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


def test_get_joutstructures_start_idx():
    start_flag = "barbie"
    assert _get_joutstructures_start_idx(["ken", "barbie"], out_slice_start_flag=start_flag) == 1
    assert _get_joutstructures_start_idx(["barbie", "ken"], out_slice_start_flag=start_flag) == 0
    assert _get_joutstructures_start_idx(["ken", "ken"], out_slice_start_flag=start_flag) is None