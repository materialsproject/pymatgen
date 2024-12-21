from __future__ import annotations


def test_JDFTXInfile():
    from pymatgen.io import jdftx

    assert hasattr(jdftx, "JDFTXInfile")


def test_JDFTXOutfile():
    from pymatgen.io import jdftx

    assert hasattr(jdftx, "JDFTXOutfile")
