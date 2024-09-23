# import pytest
from __future__ import annotations

from pymatgen.io.jdftx.generic_tags import IntTag, TagContainer

dummy_tagcontainer = TagContainer(
    allow_list_representation=True,
    can_repeat=True,
    subtags={
        "s0": IntTag(write_tagname=False, optional=True),
        "s1": IntTag(write_tagname=False, optional=True),
        "s2": IntTag(write_tagname=False, optional=True),
    },
)
dummy_tagcontainer.validate_value_type("s0", [[1]])
