# Temporary test module for turning lattice tag into mutiformat tag

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from pymatgen.io.jdftx.jdftxinfile_master_format import get_tag_object

from .shared_test_utils import assert_same_value

if TYPE_CHECKING:
    from pymatgen.io.jdftx.generic_tags import MultiformatTag, TagContainer


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
