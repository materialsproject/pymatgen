from __future__ import annotations

import pytest

from pymatgen.io.jdftx.generic_tags import _flatten_list
from pymatgen.io.jdftx.inputs import _multi_getattr as multi_getattr
from pymatgen.io.jdftx.inputs import _multi_hasattr as multi_hasattr


def test_flatten_list():
    assert _flatten_list("", [1, 2, 3]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, 3]]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, [3]]]) == [1, 2, 3]
    assert _flatten_list("", [1, [2, [3, [4]]]]) == [1, 2, 3, 4]
    with pytest.raises(TypeError):
        _flatten_list("", 1)


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
