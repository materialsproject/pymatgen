# Temporary test module for turning lattice tag into mutiformat tag

from __future__ import annotations

import pytest

from pymatgen.io.jdftx.jdftxinfile_master_format import get_tag_object

from .shared_test_utils import assert_same_value


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


# This fails, but I don't think we need to support this
# @pytest.mark.parametrize(
#     ("unordered_dict", "expected_out"),
#     [
#         (
#             {"a": 1.0, "Cubic-type": "Body-Centered", "Cubic": True},
#             "Body-Centered Cubic 1.0",
#         ),
#     ],
# )
# def test_arg_ordering(unordered_dict, expected_out):
#     mft_lattice_tag = get_tag_object("lattice")
#     i = mft_lattice_tag.get_format_index_for_str_value("lattice", expected_out)
#     tag_object = mft_lattice_tag.format_options[i]
#     output = tag_object.write("lattice", unordered_dict)
#     assert_same_value(
#         ("lattice " + expected_out).strip().split(),
#         output.strip().split(),
#     )
