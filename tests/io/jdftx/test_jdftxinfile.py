from __future__ import annotations

import os
import re
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from pymatgen.core.structure import Structure
from pymatgen.io.jdftx.inputs import JDFTXInfile, JDFTXStructure
from pymatgen.io.jdftx.jdftxinfile_master_format import get_tag_object

from .conftest import assert_same_value, dump_files_dir
from .conftest import ex_in_files_dir as ex_files_dir

if TYPE_CHECKING:
    from collections.abc import Callable

    from pymatgen.util.typing import PathLike

ex_infile1_fname = ex_files_dir / "CO.in"
ex_infile1_knowns = {
    "dump-name": "$VAR",
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

ex_infile2_fname = ex_files_dir / "example_sp.in"
ex_infile3_fname = ex_files_dir / "ct_slab_001.in"


def test_jdftxinfile_structuregen():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jdftxstruct = jif.to_jdftxstructure(jif)
    assert isinstance(jdftxstruct, JDFTXStructure)


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
    with pytest.raises(ValueError, match="Phonon functionality has not been added!"):
        jif._preprocess_line("phonon idk")
    with pytest.raises(ValueError, match="Wannier functionality has not been added!"):
        jif._preprocess_line("wannier idk")
    err_str = f"The barbie tag in {['barbie', 'ken allan']} is not in MASTER_TAG_LIST and is not a comment, "
    err_str += "something is wrong with this input data!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        jif._preprocess_line("barbie ken allan")

    _filename = "barbie"
    err_str = f"The include file {_filename} ({_filename}) does not exist!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        JDFTXInfile.from_str(f"include {_filename}\n")
    filename = ex_files_dir / "barbie"
    err_str = f"The include file {_filename} ({filename}) does not exist!"
    str(err_str)
    with pytest.raises(ValueError, match=re.escape(err_str)):
        JDFTXInfile.from_str(f"include {_filename}\n", path_parent=ex_files_dir)
    with pytest.raises(ValueError, match="This input file is missing required structure tags"):
        JDFTXInfile.from_str("dump End DOS\n")
    with pytest.raises(ValueError, match="Conversion type barbie is not 'list-to-dict' or 'dict-to-list'"):
        jif._needs_conversion("barbie", ["ken"])
    tag = "exchange-params"
    value = {"blockSize": 1, "nOuterVxx": "barbie"}
    err_str = str(f"The {tag} tag with value:\n{value}\ncould not be fixed!")
    with pytest.raises(ValueError, match=re.escape(err_str)):
        # Implicitly tests validate_tags
        jif[tag] = value
    jif2 = jif.copy()
    jif2.update({tag: value})
    with pytest.raises(ValueError, match=re.escape(err_str)):
        jif2.validate_tags(try_auto_type_fix=True)
    with pytest.warns(UserWarning):
        jif2.validate_tags(try_auto_type_fix=False)
    err_str = str(f"{1.2} is not a string!")
    with pytest.raises(TypeError, match=err_str):
        jif[1.2] = 3.4


def test_JDFTXInfile_niche_cases():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    tag_object, tag, value = jif._preprocess_line("dump-only")
    assert value == ("")
    tag = "elec-ex-corr"
    tag_object = get_tag_object(tag)
    value = "gga"
    params = jif.as_dict()
    err_str = f"The '{tag}' tag appears multiple times in this input when it should not!"
    with pytest.raises(ValueError, match=err_str):
        jif._store_value(params, tag_object, tag, value)
    struct = jif.to_pmg_structure(jif)
    assert isinstance(struct, Structure)
    noneout = jif.validate_tags(return_list_rep=True)
    assert noneout is None
    jif["fluid-solvent"] = {"name": "H2O", "concentration": 0.5}
    assert len(jif["fluid-solvent"]) == 1
    jif.append_tag("fluid-solvent", {"name": "H2O", "concentration": 0.5})
    assert len(jif["fluid-solvent"]) == 2


def test_JDFTXInfile_add_method():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif2 = jif.copy()
    jif3 = jif + jif2
    assert_idential_jif(jif, jif3)
    key = "elec-ex-corr"
    val_old = jif[key]
    val_new = "lda"
    assert val_old != val_new
    jif2[key] = val_new
    err_str = f"JDFTXInfiles have conflicting values for {key}: {val_old} != {val_new}"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        jif3 = jif + jif2
    key_add = "target-mu"
    val_add = 0.5
    assert key_add not in jif
    jif2 = jif.copy()
    jif2[key_add] = val_add
    jif3 = jif + jif2
    assert jif3[key_add]["mu"] == pytest.approx(val_add)


@pytest.mark.parametrize(("infile_fname", "knowns"), [(ex_infile1_fname, ex_infile1_knowns)])
def test_JDFTXInfile_knowns_simple(infile_fname: PathLike, knowns: dict):
    jif = JDFTXInfile.from_file(infile_fname)
    for key, val in knowns.items():
        assert_same_value(jif[key], val)


@pytest.mark.parametrize("infile_fname", [ex_infile3_fname, ex_infile1_fname, ex_infile2_fname])
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
            assert_idential_jif(jifs[i], jifs[j])
    os.remove(tmp_fname)


def assert_idential_jif(jif1: JDFTXInfile | dict, jif2: JDFTXInfile | dict):
    djif1 = jif1.as_dict() if isinstance(jif1, JDFTXInfile) else jif1
    djif2 = jif2.as_dict() if isinstance(jif2, JDFTXInfile) else jif2
    assert_same_value(djif1, djif2)


def test_jdftxstructure():
    jif = JDFTXInfile.from_file(ex_infile2_fname)
    struct = jif.to_jdftxstructure(jif)
    assert isinstance(struct, JDFTXStructure)
    struc_str = str(struct)
    assert isinstance(struc_str, str)
    assert struct.natoms == 16
    with open(ex_infile2_fname) as f:
        lines = list.copy(list(f))
    data = "\n".join(lines)
    struct2 = JDFTXStructure.from_str(data)
    assert_equiv_jdftxstructure(struct, struct2)
    struct3 = JDFTXStructure.from_dict(struct.as_dict())
    assert_equiv_jdftxstructure(struct, struct3)


def test_pmg_struc():
    jif = JDFTXInfile.from_file(ex_infile2_fname)
    struc1 = jif.to_pmg_structure(jif)
    struc2 = jif.structure
    for s in [struc1, struc2]:
        assert isinstance(s, Structure)
    assert_idential_jif(struc1.as_dict(), struc2.as_dict())


def test_jdftxtructure_naming():
    struct = Structure.from_file(ex_files_dir / "Si.cif")
    jstruct = JDFTXStructure(structure=struct)
    JDFTXInfile.from_jdftxstructure(jstruct)
    JDFTXInfile.from_structure(struct)


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
