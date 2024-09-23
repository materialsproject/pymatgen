from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
from pytest import approx

from pymatgen.io.jdftx.jdftxinfile import JDFTXInfile
from pymatgen.util.testing import TEST_FILES_DIR

if TYPE_CHECKING:
    from collections.abc import Callable

    from pymatgen.util.typing import PathLike


ex_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "example_files"
dump_files_dir = Path(TEST_FILES_DIR) / "io" / "jdftx" / "new_files"

ex_infile1_fname = ex_files_dir / "CO.in"
ex_infile1_knowns = {
    "dump-name": "$VAR",
    "initial-state": "$VAR",
    "elec-ex-corr": "gga",
    "van-der-waals": "D3",
    "elec-cutoff": {"Ecut": 20.0, "EcutRho": 100.0},
    "elec-n-bands": 15,
    # "kpoint-folding": [1,1,1],
    "kpoint-folding": {"n0": 1, "n1": 1, "n2": 1},
    "spintype": "z-spin",
    "core-overlap-check": "none",
    "converge-empty-states": True,
    "latt-move-scale": {"s0": 0.0, "s1": 0.0, "s2": 0.0},
    # "latt-move-scale": [0,0,0],
    "symmetries": "none",
    "fluid": {"type": "LinearPCM"},
    "pcm-variant": "CANDLE",
    "fluid-solvent": [{"name": "H2O"}],
    "fluid-cation": {"name": "Na+", "concentration": 0.5},
    "fluid-anion": {"name": "F-", "concentration": 0.5},
}


@pytest.mark.parametrize(
    ("infile_fname", "bool_func"),
    [
        (ex_infile1_fname, lambda jif: all(jif["kpoint-folding"][x] == 1 for x in jif["kpoint-folding"])),
        (ex_infile1_fname, lambda jif: jif["elec-n-bands"] == 15),
    ],
)
def test_JDFTXInfile_known_lambda(infile_fname: str, bool_func: Callable[[JDFTXInfile], bool]):
    jif = JDFTXInfile.from_file(infile_fname)
    assert bool_func(jif)


@pytest.mark.parametrize(
    ("val_key", "val"),
    [
        ("lattice", np.eye(3)),
        ("fluid-solvent", "H2O 0.5"),
        ("fluid-solvent", "H2O"),
        ("latt-scale", "1 1 1"),
        ("latt-scale", ["1", "1", "1"]),
        ("latt-scale", [1, 1, 1]),
        ("latt-scale", {"s0": 1, "s1": 1, "s2": 1}),
    ],
)
def test_JDFTXInfile_set_values(val_key: str, val: Any):
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif[val_key] = val
    JDFTXInfile_self_consistency_tester(jif)


@pytest.mark.parametrize(
    ("val_key", "val"),
    [
        ("fluid-solvent", "H2O 0.5"),
        ("fluid-solvent", "H2O"),
        ("dump", "End DOS"),
        ("dump", "End DOS BandEigs"),
        ("dump-interval", "Electronic 1"),
        ("ion", "Fe 1 1 1 0"),
    ],
)
def test_JDFTXInfile_append_values(val_key: str, val: Any):
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif.append_tag(val_key, val)
    JDFTXInfile_self_consistency_tester(jif)


def test_JDFTXInfile_expected_exceptions():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    with pytest.raises(KeyError):
        jif["barbie"] = "ken"
    with pytest.raises(ValueError, match="The initial-state tag cannot be repeated and thus cannot be appended"):
        jif.append_tag("initial-state", "$VAR")


@pytest.mark.parametrize(("infile_fname", "knowns"), [(ex_infile1_fname, ex_infile1_knowns)])
def test_JDFTXInfile_knowns_simple(infile_fname: PathLike, knowns: dict):
    jif = JDFTXInfile.from_file(infile_fname)
    for key in knowns:
        assert is_identical_jif_val(jif[key], knowns[key])


@pytest.mark.parametrize("infile_fname", [ex_infile1_fname])
def test_JDFTXInfile_self_consistency(infile_fname: PathLike):
    jif = JDFTXInfile.from_file(infile_fname)
    JDFTXInfile_self_consistency_tester(jif)


def JDFTXInfile_self_consistency_tester(jif: JDFTXInfile):
    dict_jif = jif.as_dict()
    # # Commenting out tests with jif2 due to the list representation asserted
    jif2 = JDFTXInfile.get_dict_representation(JDFTXInfile.from_dict(dict_jif))
    str_list_jif = jif.get_text_list()
    str_jif = "\n".join(str_list_jif)
    jif3 = JDFTXInfile.from_str(str_jif)
    tmp_fname = dump_files_dir / "tmp.in"
    jif.write_file(tmp_fname)
    jif4 = JDFTXInfile.from_file(tmp_fname)
    jifs = [jif, jif2, jif3, jif4]
    for i in range(len(jifs)):
        for j in range(i + 1, len(jifs)):
            print(f"{i}, {j}")
            assert is_identical_jif(jifs[i], jifs[j])
    os.remove(tmp_fname)


def is_identical_jif(jif1: JDFTXInfile | dict, jif2: JDFTXInfile | dict):
    for key in jif1:
        if key not in jif2:
            return False
        v1 = jif1[key]
        v2 = jif2[key]
        assert is_identical_jif_val(v1, v2)
    return True


def is_identical_jif_val(v1, v2):
    if not isinstance(v1, type(v2)):
        # if type(v1) != type(v2):
        return False
    if isinstance(v1, float):
        return v1 == approx(v2)
    if True in [isinstance(v1, str), isinstance(v1, int)]:
        return v1 == v2
    if True in [isinstance(v1, list)]:
        if len(v1) != len(v2):
            return False
        return all(is_identical_jif_val(v, v2[i]) for i, v in enumerate(v1))
    if True in [isinstance(v1, dict)]:
        return is_identical_jif(v1, v2)
    return None
