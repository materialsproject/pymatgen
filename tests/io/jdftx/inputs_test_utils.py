"""Shared test utilities for JDFTx input files.

This module contains shared testing functions and example data + known values for JDFTx input files.
This module will be combined with shared_test_utils.py upon final implementation.
"""

from __future__ import annotations

from pathlib import Path

from pymatgen.io.jdftx.inputs import JDFTXInfile, JDFTXStructure
from pymatgen.util.testing import TEST_FILES_DIR

from .shared_test_utils import assert_same_value


def assert_idential_jif(jif1: JDFTXInfile | dict, jif2: JDFTXInfile | dict):
    djif1 = jif1.as_dict() if isinstance(jif1, JDFTXInfile) else jif1
    djif2 = jif2.as_dict() if isinstance(jif2, JDFTXInfile) else jif2
    assert_same_value(djif1, djif2)


def assert_equiv_jdftxstructure(struc1: JDFTXStructure, struc2: JDFTXStructure) -> None:
    """Check if two JDFTXStructure objects are equivalent.

    Check if two JDFTXStructure objects are equivalent.

    Parameters:
    ----------
    struc1: JDFTXStructure
        The first JDFTXStructure object.
    struc2: JDFTXStructure
        The second JDFTXStructure object.
    """
    d1 = struc1.as_dict()
    d2 = struc2.as_dict()
    assert_idential_jif(d1, d2)


ex_in_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "test_jdftx_in_files"

ex_infile1_fname = ex_in_files_dir / "CO.in"
ex_infile1_knowns = {
    "dump-name": {"format": "$VAR"},
    "initial-state": "$VAR",
    "elec-ex-corr": "gga",
    "van-der-waals": "D3",
    "elec-cutoff": {"Ecut": 20.0, "EcutRho": 100.0},
    "elec-n-bands": 15,
    "kpoint-folding": {"n0": 1, "n1": 1, "n2": 1},
    "spintype": "z-spin",
    "core-overlap-check": "none",
    "converge-empty-states": True,
    "latt-move-scale": {"s0": 0.0, "s1": 0.0, "s2": 0.0},
    "symmetries": "none",
    "fluid": {"type": "LinearPCM"},
    "pcm-variant": "CANDLE",
    "fluid-solvent": [{"name": "H2O"}],
    "fluid-cation": {"name": "Na+", "concentration": 0.5},
    "fluid-anion": {"name": "F-", "concentration": 0.5},
    "initial-magnetic-moments": "C 1 O 1",
}

ex_infile2_fname = ex_in_files_dir / "example_sp.in"
ex_infile3_fname = ex_in_files_dir / "ct_slab_001.in"
