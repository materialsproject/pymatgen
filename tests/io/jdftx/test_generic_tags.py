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
    MultiformatTag,
    StrTag,
    TagContainer,
)
from pymatgen.io.jdftx.jdftxinfile_master_format import get_dump_tag_container, get_tag_object

from .shared_test_utils import assert_same_value


class Unstringable:
    """Dummy class that cannot be converted to a string"""

    def __str__(self):
        raise ValueError("Cannot convert to string")


class NonIterable:
    """Dummy class that cannot be iterated through"""

    def __iter__(self):
        raise ValueError("Cannot iterate through this object")


def test_abstract_tag():
    with pytest.raises(TypeError):
        # AbstractTag cannot be instantiated directly
        AbstractTag()


def test_stringify():
    str_tag = StrTag(options=["ken"])
    out_str = str(str_tag)
    assert isinstance(out_str, str)
    assert out_str


def test_bool_tag():
    """Expected behavior of BoolTag is tested here"""
    bool_tag = BoolTag(write_value=False)
    # Values with spaces are impossible to interpret as bools
    tag = "barbie"
    value = "this string has spaces"
    with pytest.raises(ValueError, match=f"'{value}' for '{tag}' should not have a space in it!"):
        bool_tag.read(tag, value)
    # Value errors should be raised if value cannot be conveniently converted to a bool
    value = "non-bool-like"
    with pytest.raises(ValueError, match=f"Could not set '{value}' as True/False for tag '{tag}'!"):
        bool_tag.read(tag, value)
    # bool_tag = BoolTag(write_value=True)
    # Only values appearing in the "_TF_options" map can be converted to bool (allows for yes/no to be interpreted
    # as bools)
    value = "not-appearing-in-read-TF-options"
    with pytest.raises(ValueError, match=f"Could not set '{value}' as True/False for tag '{tag}'!"):
        bool_tag.read(tag, value)


def test_abstract_tag_inheritor():
    """Expected behavior of methods inherited from AbstractTag are tested here
    (AbstractTag cannot be directly initiated)"""
    bool_tag = BoolTag(write_value=False)
    tag = "barbie"
    value = "ken"
    # Abstract-tag inheritors should raise the following error for calling get_list_representation unless
    # they have implemented the method (bool_tag does not and should not)
    with pytest.raises(ValueError, match=f"Tag object with tag '{tag}' has no get_list_representation method"):
        bool_tag.get_list_representation(tag, value)
    # Abstract-tag inheritors should raise the following error for calling get_dict_representation unless
    # they have implemented the method (bool_tag does not and should not)
    with pytest.raises(ValueError, match=f"Tag object with tag '{tag}' has no get_dict_representation method"):
        bool_tag.get_dict_representation(tag, value)
    # "_validate_repeat" is only called if "can_repeat" is True, but must raise an error if value is not a list
    with pytest.raises(TypeError):
        bool_tag._validate_repeat(tag, value)


def test_str_tag():
    str_tag = StrTag(options=["ken"])
    # Values with spaces are rejected as spaces are occasionally used as delimiters in input files
    tag = "barbie"
    value = "ken, allan"
    with pytest.raises(ValueError, match=f"'{value}' for '{tag}' should not have a space in it!"):
        str_tag.read(tag, value)
    # Values that cannot be converted to strings have a safety net error message
    value = Unstringable()
    print_value = "(unstringable)"
    with pytest.raises(TypeError, match=re.escape(f"Value '{print_value}' for '{tag}' should be a string!")):
        str_tag.read(tag, value)
    # "str_tag" here was initiated with only "ken" as an option, so "allan" should raise an error
    # (barbie is the tagname, not the value)
    value = "allan"
    with pytest.raises(
        ValueError, match=re.escape(f"The string value '{value}' must be one of {str_tag.options} for tag '{tag}'")
    ):
        str_tag.read(tag, value)
    # If both the tagname and value are written, two tokens are to be written
    str_tag = StrTag(write_tagname=True, write_value=True)
    assert str_tag.get_token_len() == 2
    # If only the tagname is written, one token is to be written
    str_tag = StrTag(write_tagname=True, write_value=False)
    assert str_tag.get_token_len() == 1
    # If only the value is written, one token is to be written
    str_tag = StrTag(write_tagname=False, write_value=True)
    assert str_tag.get_token_len() == 1
    # If neither the tagname nor the value are written, no tokens are to be written
    # (this is useless, but it is a valid option)
    str_tag = StrTag(write_tagname=False, write_value=False)
    assert str_tag.get_token_len() == 0


def test_int_tag():
    int_tag = IntTag()
    # Values with spaces are rejected as spaces are occasionally used as delimiters in input files
    value = "ken, allan"
    tag = "barbie"
    with pytest.raises(ValueError, match=f"'{value}' for '{tag}' should not have a space in it!"):
        int_tag.read(tag, value)
    # Values passed to "read" must be strings
    value = {}
    with pytest.raises(TypeError, match=f"Value '{value}' for '{tag}' should be a string!"):
        int_tag.read(tag, value)
    # Values must be able to be type-cast to an integer
    value = "ken"  # (ken is not an integer)
    with pytest.raises(ValueError, match=f"Could not set value '{value}' to an int for tag '{tag}'!"):
        int_tag.read(tag, value)


def test_float_tag():
    float_tag = FloatTag()
    tag = "barbie"
    value = "ken, allan"
    with pytest.raises(ValueError, match=f"'{value}' for '{tag}' should not have a space in it!"):
        float_tag.read(tag, value)
    with pytest.raises(TypeError):
        float_tag.read(tag, {})
    value = "ken"  # (ken is not an number)
    with pytest.raises(ValueError, match=f"Could not set value '{value}' to a float for tag '{tag}'!"):
        float_tag.read(tag, value)
    value = "1.2.3"  # (1.2.3 cannot be a float)
    with pytest.raises(ValueError, match=f"Could not set value '{value}' to a float for tag '{tag}'!"):
        float_tag.read("barbie", "1.2.3")


def test_initmagmomtag():
    initmagmomtag = InitMagMomTag(write_tagname=True)
    tag = "tag"
    value = Unstringable()
    print_value = "(unstringable)"
    with pytest.raises(TypeError, match=re.escape(f"Value '{print_value}' for '{tag}' should be a string!")):
        initmagmomtag.read(tag, value)
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
    with pytest.warns(Warning, match="warning"):
        initmagmomtag.validate_value_type("tag", Unstringable(), try_auto_type_fix=True)


def test_tagcontainer_validation():
    tag = "barbie"
    repeatable_str_subtag = "ken"
    non_repeatable_int_subtag = "allan"
    non_intable_value = "nonintable"
    tagcontainer = TagContainer(
        can_repeat=False,
        subtags={
            f"{repeatable_str_subtag}": StrTag(can_repeat=True),
            f"{non_repeatable_int_subtag}": IntTag(can_repeat=False),
        },
    )
    # issues with converting values to the correct type should only raise a warning within validate_value_type
    with pytest.warns(Warning, match="Invalid value"):
        tagcontainer.validate_value_type(
            f"{tag}", {f"{repeatable_str_subtag}": ["1"], f"{non_repeatable_int_subtag}": f"{non_intable_value}"}
        )
    # _validate_single_entry for TagContainer should raise an error if the value is not a dict
    with pytest.raises(TypeError):
        tagcontainer._validate_single_entry("not a dict")  # Not a dict
    # inhomogeneous values for repeated TagContainers should raise an error
    # until filling with default values is implemented
    tagcontainer.can_repeat = True
    value = [{"ken": 1}, {"ken": 2, "allan": 3}]
    with pytest.raises(
        ValueError,
        match=re.escape(f"The values '{value}' for tag '{tag}' provided in a list of lists have different lengths"),
    ):
        tagcontainer.validate_value_type(tag, value)


def test_tagcontainer_mixed_nesting():
    tag = "barbie"
    str_subtag = "ken"
    int_subtag = "allan"
    tagcontainer = TagContainer(
        can_repeat=False,
        subtags={
            f"{str_subtag}": StrTag(),
            f"{int_subtag}": IntTag(),
        },
    )
    tagcontainer_mixed_nesting_tester(tagcontainer, tag, str_subtag, int_subtag)


def tagcontainer_mixed_nesting_tester(tagcontainer, tag, str_subtag, int_subtag):
    list_of_dicts_and_lists = [{str_subtag: "b"}, [[int_subtag, 1]]]
    list_of_dicts_and_lists_err = f"tag '{tag}' with value '{list_of_dicts_and_lists}' cannot have"
    " nested lists/dicts mixed with bool/str/int/floats!"
    list_of_dicts_of_strs_and_ints = [{str_subtag: "b"}, {int_subtag: 1}]
    list_of_dicts_of_strs_and_ints_err = f"tag '{tag}' with value '{list_of_dicts_of_strs_and_ints}' "
    "cannot have nested dicts mixed with bool/str/int/floats!"
    list_of_lists_of_strs_and_ints = [[str_subtag, "b"], [int_subtag, 1]]
    list_of_lists_of_strs_and_ints_err = f"tag '{tag}' with value '{list_of_lists_of_strs_and_ints}' "
    "cannot have nested lists mixed with bool/str/int/floats!"
    for value, err_str in zip(
        [list_of_dicts_and_lists, list_of_dicts_of_strs_and_ints, list_of_lists_of_strs_and_ints],
        [list_of_dicts_and_lists_err, list_of_dicts_of_strs_and_ints_err, list_of_lists_of_strs_and_ints_err],
        strict=False,
    ):
        with pytest.raises(
            ValueError,
            match=re.escape(err_str),
        ):
            tagcontainer._check_for_mixed_nesting(tag, value)


def test_tagcontainer_read():
    tag = "barbie"
    repeatable_str_subtag = "ken"
    non_repeatable_int_subtag = "allan"
    tagcontainer = TagContainer(
        can_repeat=False,
        subtags={
            f"{repeatable_str_subtag}": StrTag(can_repeat=True),
            f"{non_repeatable_int_subtag}": IntTag(can_repeat=False),
        },
    )
    # non-repeatable subtags that repeat in the value for "read" should raise an error
    value = f"{non_repeatable_int_subtag} 1 {non_repeatable_int_subtag} 2"
    with pytest.raises(
        ValueError,
        match=f"Subtag '{non_repeatable_int_subtag}' for tag '{tag}' is not allowed to repeat "
        f"but repeats value {value}",
    ):
        tagcontainer.read(tag, value)
    # output of "read" should be a dict of subtags, with list values for repeatable subtags and single values for
    # non-repeatable subtags
    assert_same_value(
        tagcontainer.read(f"{tag}", "ken a ken b allan 3"),
        {"ken": ["a", "b"], "allan": 3},
    )
    required_subtag = "ken"
    optional_subtag = "allan"
    tagcontainer = TagContainer(
        can_repeat=True,
        write_tagname=True,
        subtags={
            required_subtag: StrTag(optional=False, write_tagname=True, write_value=True),
            optional_subtag: IntTag(optional=True, write_tagname=True, write_value=True),
        },
    )
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"The subtag '{required_subtag}' for tag '{tag}' is not optional but was not populated during the read!"
        ),
    ):
        tagcontainer.read("barbie", "allan 1")
    unread_values = "fgfgfgf"
    value = f"ken a {unread_values}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Something is wrong in the JDFTXInfile formatting, the following values for tag '{tag}' "
            f"were not processed: {[unread_values]}"
        ),
    ):
        tagcontainer.read(tag, value)


def test_tagcontainer_write():
    tag = "barbie"
    repeatable_str_subtag = "ken"
    non_repeatable_int_subtag = "allan"
    tagcontainer = TagContainer(
        can_repeat=False,
        subtags={
            f"{repeatable_str_subtag}": StrTag(can_repeat=True),
            f"{non_repeatable_int_subtag}": IntTag(can_repeat=False),
        },
    )
    assert_same_value(
        tagcontainer.write(tag, {repeatable_str_subtag: ["a", "b"]}),
        f"{tag} {repeatable_str_subtag} a {repeatable_str_subtag} b ",
    )
    tagcontainer.subtags[repeatable_str_subtag].write_tagname = False
    assert_same_value(tagcontainer.write(tag, {repeatable_str_subtag: ["a", "b"]}), f"{tag} a b ")
    # Lists are reserved for repeatable tagcontainers
    value = [{"ken": 1}]
    with pytest.raises(
        TypeError,
        match=re.escape(
            f"The value '{value}' (of type {type(value)}) for tag '{tag}' must be a dict for this TagContainer!"
        ),
    ):
        tagcontainer.write(tag, value)


def test_tagcontainer_list_dict_conversion():
    top_subtag = "universe"
    bottom_subtag1 = "sun"
    bottom_subtag2 = "moon"
    tagcontainer = TagContainer(
        can_repeat=True,
        allow_list_representation=False,
        subtags={
            top_subtag: TagContainer(
                allow_list_representation=True,
                write_tagname=True,
                subtags={
                    bottom_subtag1: BoolTag(
                        write_tagname=False,
                        allow_list_representation=False,
                    ),
                    bottom_subtag2: BoolTag(
                        write_tagname=True,
                        allow_list_representation=False,
                    ),
                },
            )
        },
    )
    notadict = f"True {bottom_subtag2} False"
    value = {top_subtag: notadict}
    with pytest.raises(
        ValueError, match=re.escape(f"The subtag {top_subtag} is not a dict: '{notadict}', so could not be converted")
    ):
        tagcontainer._make_list(value)
    # Despite bottom_subtag2 have "write_tagname" set to True, it is not written in the list for _make_list
    # (obviously this is dangerous and this method should be avoided)
    assert_same_value(
        tagcontainer._make_list({top_subtag: {bottom_subtag1: "True", bottom_subtag2: "False"}}),
        [top_subtag, "True", "False"],
    )
    value = {top_subtag: {bottom_subtag1: {"True": True}, bottom_subtag2: "False"}}
    with pytest.warns(Warning, match="sun subtag does not allow list"):
        tagcontainer._make_list(value)
    value = [[top_subtag, "True", "False"]]
    assert_same_value(tagcontainer.get_list_representation("barbie", value), value)
    tag = "barbie"
    value = [["universe", "True", "False"], {"universe": {"sun": True, "moon": False}}]
    err_str = f"The tag '{tag}' set to value '{value}' must be a list of dicts when passed to "
    "'get_list_representation' since the tag is repeatable."
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_list_representation(tag, value)
    value = {"universe": {"sun": {"True": True}, "moon": "False"}}
    err_str = f"Value '{value}' must be a list when passed to 'get_dict_representation' since "
    f"tag '{tag}' is repeatable."
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_dict_representation("barbie", value)
    err_str = f"Value '{value}' must be a list when passed to 'get_list_representation' since "
    f"tag '{tag}' is repeatable."
    with pytest.raises(ValueError, match=re.escape(err_str)):
        tagcontainer.get_list_representation("barbie", value)


def test_dumptagcontainer():
    dtc = get_dump_tag_container()
    tag = "dump"
    unread_value = "barbie"
    value = f"{unread_value} End DOS"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Something is wrong in the JDFTXInfile formatting, the following values for tag '{tag}' "
            f"were not processed: {[unread_value]}"
        ),
    ):
        dtc.read(tag, value)


def test_booltagcontainer():
    tag = "barbie"
    required_subtag = "ken"
    optional_subtag = "allan"
    not_a_subtag = "alan"
    btc = BoolTagContainer(
        subtags={
            required_subtag: BoolTag(optional=False, write_tagname=True, write_value=False),
            optional_subtag: BoolTag(optional=True, write_tagname=True, write_value=False),
        },
    )
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"The subtag '{required_subtag}' for tag '{tag}' is not optional but was not populated during the read!"
        ),
    ):
        btc.read(tag, optional_subtag)
    value = f"{required_subtag} {not_a_subtag}"
    with pytest.raises(
        ValueError,
        match=re.escape(
            f"Something is wrong in the JDFTXInfile formatting, the following values for tag '{tag}' "
            f"were not processed: {[not_a_subtag]}"
        ),
    ):
        btc.read(tag, value)


def test_multiformattagcontainer():
    tag = "dump-fermi-density"
    value = "notabool"
    mftg = get_tag_object(tag)
    with pytest.raises(RuntimeError):
        mftg.read(tag, value)
    with pytest.raises(RuntimeError):
        mftg.write(tag, value)
    errormsg = f"No valid read format for tag '{tag}' with value '{value}'\n"
    "Add option to format_options or double-check the value string and retry!\n\n"
    with pytest.raises(ValueError, match=re.escape(errormsg)):
        mftg.get_format_index_for_str_value(tag, value)
    err_str = f"The format for tag '{tag}' with value '{value}' could not be determined from the available options! "
    "Check your inputs and/or MASTER_TAG_LIST!"
    with pytest.raises(ValueError, match=re.escape(err_str)):
        mftg._determine_format_option(tag, value)


def test_boundary_checking():
    # Check that non-numeric tag returns valid always
    tag = "barbie"
    value = "notanumber"
    valtag = StrTag()
    assert valtag.validate_value_bounds(tag, value)[0] is True
    # Check that numeric tags can return False
    value = 0.0
    valtag = FloatTag(lb=1.0)
    assert valtag.validate_value_bounds(tag, value)[0] is False
    valtag = FloatTag(lb=0.0, lb_incl=False)
    assert valtag.validate_value_bounds(tag, value)[0] is False
    valtag = FloatTag(ub=-1.0)
    assert valtag.validate_value_bounds(tag, value)[0] is False
    valtag = FloatTag(ub=0.0, ub_incl=False)
    assert valtag.validate_value_bounds(tag, value)[0] is False
    # Check that numeric tags can return True
    valtag = FloatTag(lb=0.0, lb_incl=True)
    assert valtag.validate_value_bounds(tag, value)[0] is True
    valtag = FloatTag(ub=0.0, ub_incl=True)
    assert valtag.validate_value_bounds(tag, value)[0] is True
    valtag = FloatTag(lb=-1.0, ub=1.0)
    assert valtag.validate_value_bounds(tag, value)[0] is True
    # Check functionality for tagcontainers
    tagcontainer = TagContainer(
        subtags={
            "ken": FloatTag(lb=0.0, lb_incl=True),
            "allan": StrTag(),
            "skipper": FloatTag(ub=1.0, ub_incl=True, lb=-1.0, lb_incl=False),
        },
    )
    valid, errors = tagcontainer.validate_value_bounds(tag, {"ken": -1.0, "allan": "notanumber", "skipper": 2.0})
    assert valid is False
    assert "allan" not in errors
    assert "ken" in errors
    assert "x >= 0.0" in errors
    assert "skipper" in errors
    assert "1.0 >= x > -1.0" in errors
    # Make sure tags will never write a value that is out of bounds
    valtag = FloatTag(lb=-1.0, ub=1.0)
    assert len(valtag.write(tag, 0.0))
    assert not len(valtag.write(tag, 2.0))


def test_tag_is_equal_to():
    strtag1 = StrTag()
    strtag2 = StrTag()
    floattag1 = FloatTag()
    floattag2 = FloatTag()
    inttag1 = IntTag()
    inttag2 = IntTag()
    booltag1 = BoolTag()
    booltag2 = BoolTag()
    tagcont1 = TagContainer(subtags={"strtag": strtag1, "floattag": floattag1, "inttag": inttag1, "booltag": booltag1})
    tagcont2 = TagContainer(subtags={"strtag": strtag2, "floattag": floattag2, "inttag": inttag2, "booltag": booltag2})
    # Make sure tags of different types are not equal
    assert not strtag1.is_equal_to("", floattag1, 1.0)
    assert not strtag1.is_equal_to("", inttag1, 1)
    assert not strtag1.is_equal_to("", booltag1, True)  # noqa: FBT003
    assert not strtag1.is_equal_to("", tagcont1, {"strtag": "a", "floattag": 1.0, "inttag": 1, "booltag": True})
    # Test some str tag equalities
    assert strtag1.is_equal_to(" a", strtag2, "a")
    assert not strtag1.is_equal_to(" a", strtag2, "b")
    # Test some float tag equalities
    assert floattag1.is_equal_to(1.0, floattag2, 1.0)
    assert not floattag1.is_equal_to(1.0, floattag2, 2.0)
    floattag3 = FloatTag(eq_atol=0.1)
    assert floattag3.is_equal_to(1.0, floattag2, 1.01)
    # Test some int tag equalities
    assert inttag1.is_equal_to(1, inttag2, 1)
    assert not inttag1.is_equal_to(1, inttag2, 2)
    # Test some bool tag equalities
    assert booltag1.is_equal_to(True, booltag2, True)  # noqa: FBT003
    assert not booltag1.is_equal_to(True, booltag2, False)  # noqa: FBT003
    # Test some tagcontainer equalities
    assert tagcont1.is_equal_to(
        {"strtag": "a", "floattag": 1.0, "inttag": 1, "booltag": True},
        tagcont2,
        {"strtag": "a", "floattag": 1.0, "inttag": 1, "booltag": True},
    )
    assert not tagcont1.is_equal_to(
        {"strtag": "a", "floattag": 1.0, "inttag": 1, "booltag": True},
        tagcont2,
        {"strtag": "b", "floattag": 1.0, "inttag": 1, "booltag": True},
    )
    assert not tagcont1.is_equal_to(
        {"strtag": "a", "inttag": 1, "booltag": True},
        tagcont2,
        {"strtag": "a", "floattag": 2.0, "inttag": 1, "booltag": True},
    )
    # Test expected error raising
    with pytest.raises(ValueError, match="Values must be in dictionary format for TagContainer comparison"):
        tagcont1.is_equal_to(
            {"strtag": "a", "floattag": 1.0, "inttag": 1, "booltag": True},
            tagcont2,
            ["strtag", "a", "floattag", 1.0, "inttag", 1, "booltag", True],
        )
    with pytest.raises(ValueError, match="Both values must be strings for StrTag comparison"):
        strtag1.is_equal_to(" a", strtag2, 1.0)
    mfgtag = MultiformatTag(
        format_options=[FloatTag(), IntTag()],
    )
    with pytest.raises(NotImplementedError):
        mfgtag.is_equal_to(" a", strtag2, 1.0)
