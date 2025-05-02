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
