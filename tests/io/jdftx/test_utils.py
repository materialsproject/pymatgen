from __future__ import annotations

import pytest

from pymatgen.io.jdftx._utils import find_first_range_key, flatten_list, get_start_lines, multi_getattr, multi_hasattr
from pymatgen.io.jdftx.joutstructures import _get_joutstructures_start_idx


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


def test_get_joutstructures_start_idx():
    start_flag = "barbie"
    assert _get_joutstructures_start_idx(["ken", "barbie"], out_slice_start_flag=start_flag) == 1
    assert _get_joutstructures_start_idx(["barbie", "ken"], out_slice_start_flag=start_flag) == 0
    assert _get_joutstructures_start_idx(["ken", "ken"], out_slice_start_flag=start_flag) is None


def test_multihasattr():
    class A:
        def __init__(self):
            self.v1: int = 1

    class B:
        def __init__(self):
            self.a = A()
            self.v2: int = 2

    a = A()
    b = B()
    assert multi_hasattr(a, "v1")
    assert multi_hasattr(b, "a")
    assert multi_hasattr(b, "a.v1")
    assert not multi_hasattr(b, "a.v2")
    assert not multi_hasattr(b, "v1")


def test_multigetattr():
    class A:
        def __init__(self):
            self.v1: int = 1

    class B:
        def __init__(self):
            self.a = A()
            self.v2: int = 2

    a = A()
    b = B()
    assert multi_getattr(a, "v1") == 1
    assert multi_getattr(b, "v2") == 2
    assert multi_getattr(b, "a.v1") == 1
    with pytest.raises(AttributeError):
        multi_getattr(b, "v1")
