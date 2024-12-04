from __future__ import annotations

import re

import pytest

from pymatgen.io.jdftx.generic_tags import (
    AbstractTag,
    BoolTag,
    BoolTagContainer,
    FloatTag,
    InitMagMomTag,
    IntTag,
    StrTag,
    TagContainer,
)
from pymatgen.io.jdftx.jdftxinfile_master_format import get_dump_tag_container, get_tag_object


def test_abstract_tag():
    with pytest.raises(TypeError):
        # AbstractTag cannot be instantiated directly
        AbstractTag()


def test_stringify():
    str_tag = StrTag(options=["ken"])
    out_str = str(str_tag)
    assert isinstance(out_str, str)
    assert len(out_str)


def test_bool_tag():
    bool_tag = BoolTag(write_value=False)
    with pytest.raises(ValueError, match="Tag object has no get_list_representation method: barbie"):
        bool_tag.get_list_representation("barbie", "ken")
    with pytest.raises(ValueError, match="Tag object has no get_dict_representation method: barbie"):
        bool_tag.get_dict_representation("barbie", "ken")
    with pytest.raises(ValueError, match="'non-empty-value sdfgsd' for barbie should not have a space in it!"):
        bool_tag.read("barbie", "non-empty-value sdfgsd")
    with pytest.raises(ValueError, match="Could not set 'non-bool-like' as True/False for barbie!"):
        bool_tag.read("barbie", "non-bool-like")
    bool_tag = BoolTag(write_value=True)
    with pytest.raises(ValueError, match="Could not set 'not-appearing-in-read-TF-options' as True/False for barbie!"):
        bool_tag.read("barbie", "not-appearing-in-read-TF-options")
    with pytest.raises(TypeError):
        bool_tag._validate_repeat("barbie", "ken")


class Unstringable:
    def __str__(self):
        raise ValueError("Cannot convert to string")


class NonIterable:
    def __iter__(self):
        raise ValueError("Cannot iterate through this object")


def test_str_tag():
    str_tag = StrTag(options=["ken"])
    with pytest.raises(ValueError, match="'ken allan' for barbie should not have a space in it!"):
        str_tag.read("barbie", "ken allan")
    with pytest.raises(ValueError, match=re.escape("Could not set (unstringable) to a str for barbie!")):
        str_tag.read("barbie", Unstringable())
    with pytest.raises(ValueError, match=re.escape(f"The 'allan' string must be one of {['ken']} for barbie")):
        str_tag.read("barbie", "allan")  # Allan is not an option
    assert str_tag.get_token_len() == 2
    str_tag = StrTag(write_tagname=False, write_value=False)
    assert str_tag.get_token_len() == 0
    str_tag = StrTag(write_tagname=True, write_value=False)
    assert str_tag.get_token_len() == 1
    str_tag = StrTag(write_tagname=False, write_value=True)
    assert str_tag.get_token_len() == 1


def test_int_tag():
    int_tag = IntTag()
    with pytest.raises(ValueError, match="'ken, allan' for barbie should not have a space in it!"):
        int_tag.read("barbie", "ken, allan")
    with pytest.raises(TypeError):
        int_tag.read("barbie", {})
    with pytest.raises(ValueError, match="Could not set 'ken' to a int for barbie!"):
        int_tag.read("barbie", "ken")  # (ken is not an integer)


def test_float_tag():
    float_tag = FloatTag()
    with pytest.raises(ValueError, match="'ken, allan' for barbie should not have a space in it!"):
        float_tag.read("barbie", "ken, allan")
    with pytest.raises(TypeError):
        float_tag.read("barbie", {})
    with pytest.raises(ValueError, match="Could not set 'ken' to a float for barbie!"):
        float_tag.read("barbie", "ken")  # (ken is not an integer)
    with pytest.raises(ValueError, match="Could not set '1.2.3' to a float for barbie!"):
        float_tag.read("barbie", "1.2.3")  # (1.2.3 cannot be a float)


def test_initmagmomtag():
    initmagmomtag = InitMagMomTag(write_tagname=True)
    with pytest.raises(ValueError, match=re.escape("Could not set (unstringable) to a str for tag!")):
        initmagmomtag.read("tag", Unstringable())
    assert initmagmomtag.read("tag", "42") == "42"
    assert initmagmomtag.write("magtag", 42) == "magtag 42 "
    initmagmomtag = InitMagMomTag(write_tagname=False)
    assert initmagmomtag.write("magtag", 42) == "42 "
    initmagmomtag = InitMagMomTag(write_tagname=True, write_value=False)
    assert initmagmomtag.write("magtag", 42) == "magtag "
    initmagmomtag = InitMagMomTag(write_tagname=False, write_value=False)
    assert initmagmomtag.write("magtag", 42) == ""
    assert initmagmomtag.get_token_len() == 0
    initmagmomtag = InitMagMomTag(write_tagname=True, write_value=True)
    assert initmagmomtag.get_token_len() == 2
    with pytest.warns(Warning):
        initmagmomtag.validate_value_type("tag", Unstringable(), try_auto_type_fix=True)


def test_tagcontainer():
    tagcontainer = TagContainer(
        can_repeat=True,
        subtags={
            "ken": StrTag(),
            "allan": IntTag(can_repeat=False),
        },
    )
    with pytest.raises(TypeError):
        tagcontainer._validate_single_entry("barbie")  # Not a dict
    val = [{"ken": 1}, {"ken": 2, "allan": 3}]
    with pytest.raises(
        ValueError,
        match=re.escape(f"The values for barbie {val} provided in a list of lists have different lengths"),
    ):
        tagcontainer.validate_value_type("barbie", val)
    with pytest.warns(Warning):
        tagcontainer.validate_value_type(
            "barbie",
            [
                {"ken": "1"},
                {"allan": "barbie"},  # Raises a warning since barbie cannot be converted to int
            ],
        )
    with pytest.raises(
        ValueError, match="Subtag allan is not allowed to repeat repeats in barbie's value allan 1 allan 2"
    ):
        tagcontainer.read("barbie", "allan 1 allan 2")
    ###
    tagcontainer = TagContainer(
        can_repeat=False,
        subtags={
            "ken": StrTag(),
            "allan": IntTag(),
        },
    )
    with pytest.warns(Warning):
        tagcontainer.validate_value_type(
            "barbie",
            {"ken": "1", "allan": "barbie"},  # Raises a warning since barbie cannot be converted to int
        )
    ###
    tagcontainer = TagContainer(
        can_repeat=True,
        subtags={
            "ken": StrTag(can_repeat=True),
            "allan": IntTag(can_repeat=False),
        },
    )
    assert isinstance(tagcontainer.read("barbie", "ken 1 ken 2 allan 3"), dict)
    assert isinstance(tagcontainer.read("barbie", "ken 1 ken 2 allan 3")["ken"], list)
    assert not isinstance(tagcontainer.read("barbie", "ken 1 ken 2 allan 3")["allan"], list)
    with pytest.raises(TypeError):
        tagcontainer.write("barbie", [{"ken": 1}])
    ###
    v1 = tagcontainer.write("barbie", {"ken": [1, 2]}).strip().split()
    v1 = [v.strip() for v in v1]
    v2 = "barbie ken 1 ken 2".split()
    assert len(v1) == len(v2)
    for i in range(len(v1)):
        assert v1[i] == v2[i]
    ###
    tagcontainer = TagContainer(
        can_repeat=True,
        write_tagname=True,
        subtags={
            "ken": StrTag(optional=False, write_tagname=True, write_value=True),
            "allan": IntTag(optional=True, write_tagname=True, write_value=True),
        },
    )
    with pytest.raises(
        ValueError, match=re.escape("The ken tag is not optional but was not populated during the read!")
    ):
        tagcontainer.read("barbie", "allan 1")
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Something is wrong in the JDFTXInfile formatting, some values were not processed: ken barbie fgfgfgf"
        ),
    ):
        tagcontainer.read("barbie", "ken barbie fgfgfgf")
    assert tagcontainer.get_token_len() == 3
    with pytest.raises(TypeError):
        tagcontainer.write("barbie", ["ken barbie"])
    ###
    tagcontainer = get_tag_object("ion")
    with pytest.warns(Warning):
        tagcontainer.read("ion", "Fe 1 1 1 1 HyperPlane")
    with pytest.raises(ValueError, match="Values for repeatable tags must be a list here"):
        tagcontainer.get_dict_representation("ion", "Fe 1 1 1 1")
    ###
    tagcontainer = TagContainer(
        can_repeat=True,
        allow_list_representation=False,
        subtags={
            "ken": StrTag(),
            "allan": IntTag(),
        },
    )
    strmatch = str([{"ken": "b"}, [["allan", 1]]])
    with pytest.raises(
        ValueError,
        match=re.escape(f"barbie with {strmatch} cannot have nested lists/dicts mixed with bool/str/int/floats!"),
    ):
        tagcontainer._check_for_mixed_nesting("barbie", [{"ken": "b"}, [["allan", 1]]])
    strmatch = str([{"ken": "b"}, {"allan": 1}])
    with pytest.raises(
        ValueError,
        match=re.escape(f"barbie with {strmatch} cannot have nested dicts mixed with bool/str/int/floats!"),
    ):
        tagcontainer._check_for_mixed_nesting("barbie", [{"ken": "b"}, {"allan": 1}])
    strmatch = str([["ken", "b"], ["allan", 1]])
    with pytest.raises(
        ValueError,
        match=re.escape(f"barbie with {strmatch} cannot have nested lists mixed with bool/str/int/floats!"),
    ):
        tagcontainer._check_for_mixed_nesting("barbie", [["ken", "b"], ["allan", 1]])
    ###
    tagcontainer = TagContainer(
        can_repeat=True,
        allow_list_representation=False,
        subtags={
            "universe": TagContainer(
                allow_list_representation=True,
                write_tagname=True,
                subtags={
                    "sun": BoolTag(
                        allow_list_representation=False,
                    ),
                    "moon": BoolTag(
                        allow_list_representation=False,
                    ),
                },
            )
        },
    )
    subtag = "universe"
    value = {"universe": "True False"}
    err_str = f"The subtag {subtag} is not a dict: '{value[subtag]}', so could not be converted"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer._make_list(value)
    value = {"universe": {"sun": "True", "moon": "False"}}
    out = tagcontainer._make_list(value)
    assert isinstance(out, list)
    # This is not actually what I would expect, but keeping for sake of coverage for now
    out_expected = ["universe", "True", "False"]
    assert len(out) == len(out_expected)
    for i in range(len(out)):
        assert out[i] == out_expected[i]
    value = [{"universe": {"sun": "True", "moon": "False"}}]
    err_str = f"The value {value} is not a dict, so could not be converted"
    with pytest.raises(TypeError, match=re.escape(err_str)):
        tagcontainer._make_list(value)
    value = {"universe": {"sun": {"True": True}, "moon": "False"}}
    with pytest.warns(Warning):
        out = tagcontainer._make_list(value)
    value = [["universe", "True", "False"]]
    out = tagcontainer.get_list_representation("barbie", value)
    for i in range(len(out[0])):
        assert out[0][i] == value[0][i]
    assert len(out) == len(value)
    tag = "barbie"
    value = [["universe", "True", "False"], {"universe": {"sun": True, "moon": False}}]
    err_str = f"The {tag} tag set to {value} must be a list of dict"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_list_representation(tag, value)
    value = {"universe": {"sun": {"True": True}, "moon": "False"}}
    err_str = "Values for repeatable tags must be a list here"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_dict_representation("barbie", value)
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_list_representation("barbie", value)


def test_dumptagcontainer():
    dtc = get_dump_tag_container()
    with pytest.raises(
        ValueError,
        match=re.escape("Something is wrong in the JDFTXInfile formatting, some values were not processed: ['barbie']"),
    ):
        dtc.read("dump", "barbie End DOS")


def test_booltagcontainer():
    btc = BoolTagContainer(
        subtags={
            "ken": BoolTag(optional=False, write_tagname=True, write_value=False),
            "allan": BoolTag(optional=True, write_tagname=True, write_value=False),
        },
    )
    with pytest.raises(ValueError, match="The ken tag is not optional but was not populated during the read!"):
        btc.read("barbie", "allan")
    with pytest.raises(
        ValueError,
        match=re.escape("Something is wrong in the JDFTXInfile formatting, some values were not processed: ['aliah']"),
    ):
        btc.read("barbie", "ken aliah")


def test_multiformattagcontainer():
    tag = "dump-fermi-density"
    value = "notabool"
    mftg = get_tag_object(tag)
    with pytest.raises(RuntimeError):
        mftg.read(tag, value)
    with pytest.raises(RuntimeError):
        mftg.write(tag, value)
    errormsg = f"No valid read format for '{tag} {value}' tag\n"
    "Add option to format_options or double-check the value string and retry!\n\n"
    with pytest.raises(ValueError, match=re.escape(errormsg)):
        mftg.get_format_index_for_str_value(tag, value)
    err_str = f"The format for {tag} for:\n{value}\ncould not be determined from the available options!"
    "Check your inputs and/or MASTER_TAG_LIST!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        mftg._determine_format_option(tag, value)
