from __future__ import annotations

import re
from copy import deepcopy
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from pymatgen.core.structure import Structure
from pymatgen.io.jdftx.inputs import JDFTXInfile, JDFTXStructure
from pymatgen.io.jdftx.jdftxinfile_default_inputs import default_inputs
from pymatgen.io.jdftx.jdftxinfile_master_format import get_tag_object

from .inputs_test_utils import (
    assert_equiv_jdftxstructure,
    assert_idential_jif,
    ex_in_files_dir,
    ex_infile1_fname,
    ex_infile1_knowns,
    ex_infile2_fname,
    ex_infile3_fname,
)
from .shared_test_utils import assert_same_value

if TYPE_CHECKING:
    from collections.abc import Callable

    from pymatgen.io.jdftx.generic_tags import MultiformatTag, TagContainer
    from pymatgen.util.typing import PathLike


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


def JDFTXInfile_self_consistency_tester(jif: JDFTXInfile, tmp_path: PathLike):
    """Create an assortment of JDFTXinfile created from the same data but through different methods, and test that
    they are all equivalent through "assert_idential_jif" """
    dict_jif = jif.as_dict()
    # # Commenting out tests with jif2 due to the list representation asserted
    jif2 = JDFTXInfile.get_dict_representation(JDFTXInfile._from_dict(dict_jif))
    str_list_jif = jif.get_text_list()
    str_jif = "\n".join(str_list_jif)
    jif3 = JDFTXInfile.from_str(str_jif)
    tmp_fname = tmp_path / "tmp.in"
    jif.write_file(tmp_fname)
    jif4 = JDFTXInfile.from_file(tmp_fname)
    jifs = [jif, jif2, jif3, jif4]
    for i in range(len(jifs)):
        for j in range(i + 1, len(jifs)):
            print(f"{i}, {j}")
            assert_idential_jif(jifs[i], jifs[j])


def test_JDFTXInfile_from_dict(tmp_path) -> None:
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif_dict = jif.as_dict()
    # Test that dictionary can be modified and that _from_dict will fix set values
    jif_dict["elec-cutoff"] = 20
    jif2 = JDFTXInfile.from_dict(jif_dict)
    JDFTXInfile_self_consistency_tester(jif2, tmp_path)


@pytest.mark.parametrize("infile_fname", [ex_infile3_fname, ex_infile1_fname, ex_infile2_fname])
def test_JDFTXInfile_self_consistency_fromfile(infile_fname: PathLike, tmp_path) -> None:
    """Test that JDFTXInfile objects with different assortments of tags survive inter-conversion done within
    "JDFTXInfile_self_consistency_tester"""
    jif = JDFTXInfile.from_file(infile_fname)
    JDFTXInfile_self_consistency_tester(jif, tmp_path)


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
        ("elec-cutoff", {"Ecut": 20.0, "EcutRho": 100.0}),
        ("elec-cutoff", "20 100"),
        ("elec-cutoff", [20, 100]),
        ("elec-cutoff", 20),
    ],
)
def test_JDFTXInfile_set_values(val_key: str, val: Any, tmp_path) -> None:
    """Test value setting for various tags"""
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif[val_key] = val
    # Test that the JDFTXInfile object is still consistent
    JDFTXInfile_self_consistency_tester(jif, tmp_path)


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
def test_JDFTXInfile_append_values(val_key: str, val: Any, tmp_path) -> None:
    """Test the append_tag method"""
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    val_old = None if val_key not in jif else deepcopy(jif[val_key])
    jif.append_tag(val_key, val)
    val_new = jif[val_key]
    assert val_old != val_new
    # Test that the append_tag does not break the JDFTXInfile object
    JDFTXInfile_self_consistency_tester(jif, tmp_path)


def test_JDFTXInfile_expected_exceptions():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    with pytest.raises(KeyError):
        jif["barbie"] = "ken"
    # non-repeating tags raise value-errors when appended
    tag = "initial-state"
    with pytest.raises(ValueError, match=re.escape(f"The tag '{tag}' cannot be repeated and thus cannot be appended")):
        jif.append_tag(tag, "$VAR")
    # Phonon and Wannier tags raise value-errors at _preprocess_line
    with pytest.raises(ValueError, match="Phonon functionality has not been added!"):
        jif._preprocess_line("phonon idk")
    with pytest.raises(ValueError, match="Wannier functionality has not been added!"):
        jif._preprocess_line("wannier idk")
    # Tags not in MASTER_TAG_LIST raise value-errors at _preprocess_line
    err_str = f"The barbie tag in {['barbie', 'ken allan']} is not in MASTER_TAG_LIST and is not a comment, "
    "something is wrong with this input data!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        jif._preprocess_line("barbie ken allan")
    # include tags raise value-errors if the file cannot be found
    _filename = "barbie"
    err_str = f"The include file {_filename} ({_filename}) does not exist!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        JDFTXInfile.from_str(f"include {_filename}\n")
    # If it does exist, no error should be raised
    filename = ex_in_files_dir / "barbie"
    err_str = f"The include file {_filename} ({filename}) does not exist!"
    str(err_str)
    # If the wrong parent_path is given for a file that does exist, error
    with pytest.raises(ValueError, match=re.escape(err_str)):
        JDFTXInfile.from_str(f"include {_filename}\n", path_parent=ex_in_files_dir)
    # JDFTXInfile cannot be constructed without lattice and ion tags
    with pytest.raises(ValueError, match="This input file is missing required structure tags"):
        JDFTXInfile.from_str("dump End DOS\n")
    # "barbie" here is supposed to be "list-to-dict" or "dict-to-list"
    with pytest.raises(ValueError, match="Conversion type barbie is not 'list-to-dict' or 'dict-to-list'"):
        jif._needs_conversion("barbie", ["ken"])
    # Setting tags with unfixable values immediately raises an error
    tag = "exchange-params"
    value = {"blockSize": 1, "nOuterVxx": "barbie"}
    err_str = str(f"The {tag} tag with value:\n{value}\ncould not be fixed!")
    with pytest.raises(ValueError, match=re.escape(err_str)):
        # Implicitly tests validate_tags
        jif[tag] = value
    # Setting tags with unfixable values through "update" side-steps the error, but will raise it once
    # "validate_tags" is inevitably called
    jif2 = jif.copy()
    jif2.update({tag: value})
    with pytest.raises(ValueError, match=re.escape(err_str)):
        jif2.validate_tags(try_auto_type_fix=True)
    # The inevitable error can be reduced to a warning if you tell it not to try to fix the values
    with pytest.warns(UserWarning, match="The exchange-params tag with value"):
        jif2.validate_tags(try_auto_type_fix=False)
    # Setting a non-string tag raises an error within the JDFTXInfile object
    err_str = str(f"{1.2} is not a string!")
    with pytest.raises(TypeError, match=err_str):
        jif[1.2] = 3.4


def test_JDFTXInfile_strip_structure():
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    structural_tags = ["lattice", "ion", "coords-type"]
    assert all(tag in jif for tag in structural_tags)
    jif.strip_structure_tags()
    assert all(tag not in jif for tag in structural_tags)


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
    """Test the __add__ method"""
    # No new values are being assigned in jif2, so jif + jif2 should be the same as jif
    # Since the convenience of this method would be lost if the user has to pay special attention to duplicating
    # repeatable values, repeatable tags are not append to each other
    jif = JDFTXInfile.from_file(ex_infile1_fname)
    jif2 = jif.copy()
    assert jif2 is not jif  # Testing robustness of copy method while we are at it
    jif3 = jif + jif2
    assert_idential_jif(jif, jif3)
    # If a tag is repeated, the values must be the same since choice of value is ambiguous
    key = "elec-ex-corr"
    val_old = deepcopy(jif[key])
    val_new = "lda"
    assert val_old != val_new
    jif2[key] = val_new
    jif4 = jif + jif2
    assert_same_value(jif4[key], val_new)  # Make sure addition chooses second value for non-repeatable tags
    del jif4
    jif2.append_tag("dump", "Fluid State")
    jif.append_tag("dump", "Fluid Berry")
    jif4 = jif + jif2
    assert len(jif4["dump"]) == len(jif2["dump"]) + 1
    assert {"Fluid": {"State": True}} in jif4["dump"]
    assert {"Fluid": {"State": True}} not in jif["dump"]
    assert {"Fluid": {"Berry": True}} in jif4["dump"]
    assert {"Fluid": {"Berry": True}} not in jif2["dump"]
    # Normal expected behavior
    key_add = "target-mu"
    val_add = 0.5
    assert key_add not in jif
    jif2 = jif.copy()
    jif2[key_add] = val_add
    jif3 = jif + jif2
    assert jif3[key_add]["mu"] == pytest.approx(val_add)


@pytest.mark.parametrize(("infile_fname", "knowns"), [(ex_infile1_fname, ex_infile1_knowns)])
def test_JDFTXInfile_knowns_simple(infile_fname: PathLike, knowns: dict):
    """Test that known values that can be tested with assert_same_value are correct"""
    jif = JDFTXInfile.from_file(infile_fname)
    for key, val in knowns.items():
        assert_same_value(jif[key], val)


def test_jdftxstructure():
    """Test the JDFTXStructure object associated with the JDFTXInfile object"""
    jif = JDFTXInfile.from_file(ex_infile2_fname)
    struct = jif.to_jdftxstructure(jif)
    assert isinstance(struct, JDFTXStructure)
    struc_str = str(struct)
    assert isinstance(struc_str, str)
    newstruct = JDFTXStructure.from_str(struc_str)
    assert isinstance(newstruct, JDFTXStructure)
    # Double checking I got the column/row order right
    assert_same_value(struct.structure.lattice, newstruct.structure.lattice)
    assert struct.natoms == 16
    with open(ex_infile2_fname) as f:
        lines = list.copy(list(f))
    # Test different ways of creating a JDFTXStructure object create the same object if data is the same
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
    """Test the naming of the JDFTXStructure object.

    Test to make sure reading from a Structure with labels not exactly matching the element names
    (ie Si0, Si1, or Si+2) will still be read correctly.
    """
    struct = Structure.from_file(ex_in_files_dir / "Si.cif")
    jstruct = JDFTXStructure(structure=struct)
    JDFTXInfile.from_jdftxstructure(jstruct)
    JDFTXInfile.from_structure(struct)


@pytest.mark.parametrize(
    ("value_str", "expected_dict"),
    [
        (
            "H 1.0 1.0 1.0 0",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 0,
            },
        ),
        (
            "H 1.0 1.0 1.0 v 1.0 1.0 1.0 0",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "v": {
                    "vx0": 1.0,
                    "vx1": 1.0,
                    "vx2": 1.0,
                },
                "moveScale": 0,
            },
        ),
        (
            "H 1.0 1.0 1.0 1 Linear 1.0 1.0 1.0",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 1,
                "constraint type": "Linear",
                "d0": 1.0,
                "d1": 1.0,
                "d2": 1.0,
            },
        ),
        (
            "H 1.0 1.0 1.0 1 HyperPlane 1.0 1.0 1.0 g1 HyperPlane 1.0 1.0 1.0 g1",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 1,
                "HyperPlane": [
                    {
                        "d0": 1.0,
                        "d1": 1.0,
                        "d2": 1.0,
                        "group": "g1",
                    },
                    {
                        "d0": 1.0,
                        "d1": 1.0,
                        "d2": 1.0,
                        "group": "g1",
                    },
                ],
            },
        ),
    ],
)
def test_ion_reading(value_str: str, expected_dict: dict):
    ion_tag: MultiformatTag = get_tag_object("ion")
    i = ion_tag.get_format_index_for_str_value("ion", value_str)
    tag_object: TagContainer = ion_tag.format_options[i]
    parsed_tag = tag_object.read("ion", value_str)
    assert_same_value(parsed_tag, expected_dict)


@pytest.mark.parametrize(
    ("expected_out", "stored_dict"),
    [
        (
            "ion H 1.000000000000 1.000000000000 1.000000000000 0",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 0,
            },
        ),
        (
            "ion H 1.000000000000 1.000000000000 1.000000000000 v 1.000000000000 1.000000000000 1.000000000000 0",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "v": {
                    "vx0": 1.0,
                    "vx1": 1.0,
                    "vx2": 1.0,
                },
                "moveScale": 0,
            },
        ),
        (
            "ion H 1.000000000000 1.000000000000 1.000000000000 1 Linear 1.000000000000 1.000000000000 1.000000000000",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 1,
                "constraint type": "Linear",
                "d0": 1.0,
                "d1": 1.0,
                "d2": 1.0,
            },
        ),
        (
            "ion H 1.000000000000 1.000000000000 1.000000000000 1 HyperPlane 1.000000000000 1.000000000000 "
            "1.000000000000 g1 HyperPlane 1.000000000000 1.000000000000 1.000000000000 g1",
            {
                "species-id": "H",
                "x0": 1.0,
                "x1": 1.0,
                "x2": 1.0,
                "moveScale": 1,
                "HyperPlane": [
                    {
                        "d0": 1.0,
                        "d1": 1.0,
                        "d2": 1.0,
                        "group": "g1",
                    },
                    {
                        "d0": 1.0,
                        "d1": 1.0,
                        "d2": 1.0,
                        "group": "g1",
                    },
                ],
            },
        ),
    ],
)
def test_ion_writing(expected_out: str, stored_dict: dict):
    ion_tag: MultiformatTag = get_tag_object("ion")
    i, _ = ion_tag._determine_format_option("ion", stored_dict)
    tag_object: TagContainer = ion_tag.format_options[i]
    output = tag_object.write("ion", stored_dict)
    assert_same_value(output.strip().split(), expected_out.strip().split())


@pytest.mark.parametrize(
    ("lattice_type", "value_str"),
    [
        ("Rhombohedral", "Rhombohedral 1.0 1.0"),
        ("Triclinic", "Triclinic 1.0 1.0 1.0 1.0 1.0 1.0"),
        ("Hexagonal", "Hexagonal 1.0 1.0"),
        ("Cubic", "Body-Centered Cubic 1.0"),
        ("Cubic", "Cubic 1.0"),
    ],
)
def test_lattice_reading(lattice_type: str, value_str: str):
    mft_lattice_tag = get_tag_object("lattice")
    assert mft_lattice_tag is not None
    i = mft_lattice_tag.get_format_index_for_str_value("lattice", value_str)
    tag_object = mft_lattice_tag.format_options[i]
    parsed_tag = tag_object.read("lattice", value_str)
    assert lattice_type in parsed_tag


@pytest.mark.parametrize(
    ("value_str"),
    [
        ("Rhombohedral 1.0 1.0"),
        ("Triclinic 1.0 1.0 1.0 1.0 1.0 1.0"),
        ("Hexagonal 1.0 1.0"),
        ("Body-Centered Cubic 1.0"),
        ("Cubic 1.0"),
        ("Orthorhombic 1.0 1.0 1.0"),
        ("Base-Centered Orthorhombic 1.0 1.0 1.0"),
        ("Monoclinic 1.0 1.0 1.0 1.0"),
        ("Base-Centered Monoclinic 1.0 1.0 1.0 1.0"),
        ("Tetragonal 1.0 1.0"),
        ("Body-Centered Tetragonal 1.0 1.0"),
    ],
)
def test_lattice_writing(value_str: str):
    mft_lattice_tag = get_tag_object("lattice")
    assert mft_lattice_tag is not None
    i = mft_lattice_tag.get_format_index_for_str_value("lattice", value_str)
    tag_object = mft_lattice_tag.format_options[i]
    parsed_tag = tag_object.read("lattice", value_str)
    output = tag_object.write("lattice", parsed_tag)
    assert_same_value(
        ("lattice " + value_str).strip().split(),
        output.strip().split(),
    )


@pytest.mark.parametrize(
    ("value_str"),
    [
        ("Rhombohedral 1.0 1.0"),
        ("Triclinic 1.0 1.0 1.0 1.0 1.0 1.0"),
        ("Hexagonal 1.0 1.0"),
        ("Body-Centered Cubic 1.0"),
        ("Cubic 1.0"),
        ("Orthorhombic 1.0 1.0 1.0"),
        ("Base-Centered Orthorhombic 1.0 1.0 1.0"),
        ("Monoclinic 1.0 1.0 1.0 1.0"),
        ("Base-Centered Monoclinic 1.0 1.0 1.0 1.0"),
        ("Tetragonal 1.0 1.0"),
        ("Body-Centered Tetragonal 1.0 1.0"),
    ],
)
def test_jdftxstructure_lattice_conversion(value_str: str):
    test_vars = ["a", "b", "c", "alpha", "beta", "gamma"]
    mft_lattice_tag = get_tag_object("lattice")
    assert mft_lattice_tag is not None
    i = mft_lattice_tag.get_format_index_for_str_value("lattice", value_str)
    tag_object = mft_lattice_tag.format_options[i]
    parsed_tag = tag_object.read("lattice", value_str)
    infile = JDFTXInfile.from_str("lattice " + value_str + "\n ion H 0.0 0.0 0.0 0", dont_require_structure=True)
    if "modification" in parsed_tag:
        with pytest.raises(NotImplementedError):
            _ = infile.to_pmg_structure(infile)
    else:
        structure = infile.to_pmg_structure(infile)
        for var in test_vars:
            if var in parsed_tag:
                assert_same_value(float(getattr(structure.lattice, var)), float(parsed_tag[var]))


def test_jdftxinfile_comparison():
    jif1 = JDFTXInfile.from_file(ex_infile1_fname)
    jif2 = JDFTXInfile.from_file(ex_infile2_fname)
    assert len(jif1.get_differing_tags(jif2))  # At least one tag should be different
    jif1copy = jif1.copy()
    assert not len(jif1.get_differing_tags(jif1copy))  # No tags should be different
    default_test_tag = "davidson-band-ratio"
    default_test_val = default_inputs[default_test_tag]
    jif1[default_test_tag] = default_test_val
    assert not len(jif1.get_differing_tags(jif1copy))  # Even though jif1copy doesn't have the tag,
    # it won't be recognized as a difference since it is in the default_inputs and matches the default value
    jif1copy["elec-n-bands"] = 20001
    assert len(jif1.get_filtered_differing_tags(jif1copy))  # Change in elec-n-bands should be recognized
    assert not len(
        jif1.get_filtered_differing_tags(jif1copy, exclude_tags=["elec-n-bands"])
    )  # Specific tags can be filtered out
    assert not len(
        jif1.get_filtered_differing_tags(jif1copy, exclude_tag_categories=["electronic"])
    )  # Tag categories can be filtered out
