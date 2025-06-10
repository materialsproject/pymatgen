"""Module for class objects for containing JDFTx tags.

This module contains class objects for containing JDFTx tags. These class objects
are used to validate the type of the value for the tag, read the value string for
the tag, write the tag and its value as a string, and get the token length of the tag.
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any

import numpy as np

__author__ = "Jacob Clary, Ben Rich"


# This inheritable class is kept for use by AbstractTag instead of using pprint as
# the type conversions for tags are incredibly delicate and require strings to be
# printed in a very exact way.
class ClassPrintFormatter:
    """Generic class for printing to command line in readable format.

    Generic class for printing to command line in readable format.
    """

    def __str__(self) -> str:
        """Print the class to the command line in a readable format.

        Returns:
            str: The class in a readable format.
        """
        return f"{self.__class__}\n" + "\n".join(f"{item} = {self.__dict__[item]}" for item in sorted(self.__dict__))


@dataclass
class AbstractTag(ClassPrintFormatter, ABC):
    """Abstract base class for all tags."""

    multiline_tag: bool = (
        False  # set to True if what to print tags across multiple lines, typically like electronic-minimize
    )
    can_repeat: bool = False  # set to True for tags that can appear on multiple lines, like ion
    write_tagname: bool = True  # set to False to not print the tagname, like for subtags of elec-cutoff
    write_value: bool = True  # set to False to not print any value, like for dump-interval
    optional: bool = (
        True  # set to False if tag (usually a subtag of a TagContainer) must be set for the JDFTXInfile to be valid.
    )
    # The lattice, ion, and ion-species are the main tags that are not optional
    defer_until_struc: bool = False
    is_tag_container: bool = False
    allow_list_representation: bool = False  # if True, allow this tag to exist as a list or list of lists

    @abstractmethod
    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value.
            Defaults to False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """

    def is_equal_to(self, val1: Any | list[Any], obj2: AbstractTag, val2: Any | list[Any]) -> bool:
        """Check if the two values are equal.

        Args:
            val1 (Any): The value of this tag object.
            obj2 (AbstractTag): The other tag object.
            val2 (Any): The value of the other tag object.

        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """
        if self.can_repeat:
            if not obj2.can_repeat:
                return False
            val1 = val1 if isinstance(val1, list) else [val1]
            val2 = val2 if isinstance(val2, list) else [val2]
            if len(val1) != len(val2):
                return False
            return all(True in [self._is_equal_to(v1, obj2, v2) for v2 in val2] for v1 in val1)
        return self._is_equal_to(val1, obj2, val2)

    @abstractmethod
    def _is_equal_to(self, val1: Any, obj2: AbstractTag, val2: Any) -> bool:
        """Check if the two values are equal.

        Used to check if the two values are equal. Assumes val1 and val2 are single elements.

        Args:
            val1 (Any): The value of this tag object.
            obj2 (AbstractTag): The other tag object.
            val2 (Any): The value of the other tag object.

        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """

    def _is_same_tagtype(
        self,
        obj2: AbstractTag,
    ) -> bool:
        """Check if the two values are equal.

        Args:
            obj2 (AbstractTag): The other tag object.

        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """
        return isinstance(self, type(obj2))

    def _validate_value_type(
        self, type_check: type, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        This method is used to validate the type of the value for this tag. It is used to check if the value is of the
        correct type and if it can be fixed automatically. If the value is not of the correct type and cannot be fixed
        automatically, a warning is raised.

        Args:
            type_check (type): The type to check the value against.
            tag (str): The tag to check the value against.
            value (Any): The value to check the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value.
            Defaults to False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        if self.can_repeat:
            self._validate_repeat(tag, value)
            is_valid = all(isinstance(x, type_check) for x in value)
        else:
            is_valid = isinstance(value, type_check)

        if not is_valid and try_auto_type_fix:
            try:
                value = [self.read(tag, str(x)) for x in value] if self.can_repeat else self.read(tag, str(value))
                tag, is_valid, value = self._validate_value_type(type_check, tag, value)
            except (TypeError, ValueError):
                warning = f"Could not fix the typing for tag '{tag}'"
                try:
                    warning += f"{value}!"
                except (ValueError, TypeError):
                    warning += "(unstringable value)!"
                warnings.warn("warning", stacklevel=2)
        return tag, is_valid, value

    def _validate_repeat(self, tag: str, value: Any) -> None:
        if not isinstance(value, list):
            raise TypeError(f"The '{tag}' tag can repeat but is not a list: '{value}'")

    def validate_value_bounds(
        self,
        tag: str,
        value: Any,
    ) -> tuple[bool, str]:
        return True, ""

    @abstractmethod
    def read(self, tag: str, value_str: str) -> Any:
        """Read and parse the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value_str (str): The value string to read.

        Returns:
            Any: The parsed value.
        """

    def _general_read_validate(self, tag: str, value_str: Any) -> None:
        """General validation for values to be passed to a read method."""
        try:
            value = str(value_str)
        except (ValueError, TypeError):
            value = "(unstringable)"
        if not isinstance(value_str, str):
            raise TypeError(f"Value '{value}' for '{tag}' should be a string!")

    def _single_value_read_validate(self, tag: str, value: str) -> None:
        """Validation for values to be passed to a read method for AbstractTag inheritors that only
        read a single value."""
        self._general_read_validate(tag, value)
        if len(value.split()) > 1:
            raise ValueError(f"'{value}' for '{tag}' should not have a space in it!")

    def _check_unread_values(self, tag: str, unread_values: list[str]) -> None:
        """Check for unread values and raise an error if any are found. Used in the read method of TagContainers."""
        if len(unread_values) > 0:
            raise ValueError(
                f"Something is wrong in the JDFTXInfile formatting, the following values for tag '{tag}' "
                f"were not processed: {unread_values}"
            )

    def _check_nonoptional_subtags(self, tag: str, subdict: dict[str, Any], subtags: dict[str, AbstractTag]) -> None:
        """Check for non-optional subtags and raise an error if any are missing.
        Used in the read method of TagContainers."""
        for subtag, subtag_type in subtags.items():
            if not subtag_type.optional and subtag not in subdict:
                raise ValueError(
                    f"The subtag '{subtag}' for tag '{tag}' is not optional but was not populated during the read!"
                )

    @abstractmethod
    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """

    @abstractmethod
    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """

    def _write(self, tag: str, value: Any, multiline_override: bool = False, strip_override: bool = False) -> str:
        tag_str = f"{tag} " if self.write_tagname else ""
        if multiline_override:
            tag_str += "\\\n"
        if self.write_value:
            if not strip_override:
                tag_str += f"{value}".strip() + " "
            else:
                tag_str += f"{value}"
        return tag_str

    def _get_token_len(self) -> int:
        return int(self.write_tagname) + int(self.write_value)

    def get_list_representation(self, tag: str, value: Any) -> list | list[list]:
        """Convert the value to a list representation.

        Args:
            tag (str): The tag to convert the value to a list representation for.
            value (Any): The value to convert to a list representation.

        Returns:
            list | list[list]: The value converted to a list representation.
        """
        raise ValueError(f"Tag object with tag '{tag}' has no get_list_representation method")

    def get_dict_representation(self, tag: str, value: Any) -> dict | list[dict]:
        """Convert the value to a dict representation.

        Args:
            tag (str): The tag to convert the value to a dict representation for.
            value (Any): The value to convert to a dict representation.

        Returns:
            dict | list[dict]: The value converted to a dict representation.
        """
        raise ValueError(f"Tag object with tag '{tag}' has no get_dict_representation method")


@dataclass
class BoolTag(AbstractTag):
    """Tag for boolean values in JDFTx input files.

    Tag for boolean values in JDFTx input files.
    """

    _TF_read_options: dict[str, bool] = field(default_factory=lambda: {"yes": True, "no": False})
    _TF_write_options: dict[bool, str] = field(default_factory=lambda: {True: "yes", False: "no"})
    _TF_options: dict[str, dict] = field(init=False)

    def __post_init__(self) -> None:
        """Initialize the _TF_options attribute.

        Initialize the _TF_options attribute.
        """
        self._TF_options = {
            "read": self._TF_read_options,
            "write": self._TF_write_options,
        }

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value.
            Defaults to False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        return self._validate_value_type(bool, tag, value, try_auto_type_fix=try_auto_type_fix)

    def _is_equal_to(self, val1: Any, obj2: AbstractTag, val2: Any) -> bool:
        """Check if the two values are equal.

        Args:
            val1 (Any): The value of this tag object.
            obj2 (AbstractTag): The other tag object.
            val2 (Any): The value of the other tag object.

        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """
        return self._is_same_tagtype(obj2) and val1 == val2

    def raise_value_error(self, tag: str, value: str) -> None:
        """Raise a ValueError for the value string.

        Args:
            tag (str): The tag to raise the ValueError for.
            value (str): The value string to raise the ValueError for.
        """
        raise ValueError(f"The value '{value}' was provided to {tag}, it is not acting like a boolean")

    def read(self, tag: str, value: str) -> bool:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            bool: The parsed boolean value.
        """
        self._single_value_read_validate(tag, value)
        try:
            if not self.write_value:
                # accounts for exceptions where only the tagname is used, e.g.
                # dump-only or dump-fermi-density (sometimes) tags
                if not value:  # then the string '' was passed in because no value was provided but the tag was present
                    value = "yes"
                else:
                    self.raise_value_error(tag, value)
            return self._TF_options["read"][value]
        except (ValueError, TypeError, KeyError) as err:
            raise ValueError(f"Could not set '{value}' as True/False for tag '{tag}'!") from err

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        value2 = self._TF_options["write"][value]
        return self._write(tag, value2)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        return self._get_token_len()


@dataclass
class StrTag(AbstractTag):
    """Tag for string values in JDFTx input files.

    Tag for string values in JDFTx input files.
    """

    options: list | None = None

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value. Defaults to
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        return self._validate_value_type(str, tag, value, try_auto_type_fix=try_auto_type_fix)

    def _is_equal_to(self, val1: Any, obj2: AbstractTag, val2: Any) -> bool:
        """Check if the two values are equal.

        Args:
            val1 (Any): The value of this tag object.
            obj2 (AbstractTag): The other tag object.
            val2 (Any): The value of the other tag object.

        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """
        if self._is_same_tagtype(obj2):
            if not all(isinstance(x, str) for x in (val1, val2)):
                raise ValueError("Both values must be strings for StrTag comparison")
            return val1.strip() == val2.strip()
        return False

    def read(self, tag: str, value: str) -> str:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            str: The parsed string value.
        """
        self._single_value_read_validate(tag, value)
        if self.options is None or value in self.options:
            return value
        raise ValueError(f"The string value '{value}' must be one of {self.options} for tag '{tag}'")

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        return self._get_token_len()


@dataclass
class AbstractNumericTag(AbstractTag):
    """Abstract base class for numeric tags."""

    lb: float | None = None  # lower bound
    ub: float | None = None  # upper bound
    lb_incl: bool = True  # lower bound inclusive
    ub_incl: bool = True  # upper bound inclusive
    eq_atol: float = 1.0e-8  # absolute tolerance for equality check
    eq_rtol: float = 1.0e-5  # relative tolerance for equality check

    def val_is_within_bounds(self, value: float) -> bool:
        """Check if the value is within the bounds.

        Args:
            value (float | int): The value to check.

        Returns:
            bool: True if the value is within the bounds, False otherwise.
        """
        good = True
        if self.lb is not None:
            good = good and value >= self.lb if self.lb_incl else good and value > self.lb
        if self.ub is not None:
            good = good and value <= self.ub if self.ub_incl else good and value < self.ub
        return good

    def get_invalid_value_error_str(self, tag: str, value: float) -> str:
        """Raise a ValueError for the invalid value.

        Args:
            tag (str): The tag to raise the ValueError for.
            value (float | int): The value to raise the ValueError for.
        """
        err_str = f"Value '{value}' for tag '{tag}' is not within bounds"
        if self.ub is not None:
            err_str += f" {self.ub} >"
            if self.ub_incl:
                err_str += "="
        err_str += " x "
        if self.lb is not None:
            err_str += ">"
            if self.lb_incl:
                err_str += "="
        err_str += f" {self.lb}"
        return err_str

    def validate_value_bounds(
        self,
        tag: str,
        value: Any,
    ) -> tuple[bool, str]:
        if not self.val_is_within_bounds(value):
            return False, self.get_invalid_value_error_str(tag, value)
        return True, ""

    def _is_equal_to(self, val1, obj2, val2):
        """Check if the two values are equal.

        Used to check if the two values are equal. Doesn't need to be redefined for IntTag and FloatTag.

        Args:
            val1 (Any): The value of this tag object.
            obj2 (AbstractTag): The other tag object.
            val2 (Any): The value of the other tag object.
            rtol (float, optional): Relative tolerance. Defaults to 1.e-5.
            atol (float, optional): Absolute tolerance. Defaults to 1.e-8.
        Returns:
            bool: True if the two tag object/value pairs are equal, False otherwise.
        """
        return self._is_same_tagtype(obj2) and np.isclose(val1, val2, rtol=self.eq_rtol, atol=self.eq_atol)


@dataclass
class IntTag(AbstractNumericTag):
    """Tag for integer values in JDFTx input files.

    Tag for integer values in JDFTx input files.
    """

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value, by default
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        return self._validate_value_type(int, tag, value, try_auto_type_fix=try_auto_type_fix)

    def read(self, tag: str, value: str) -> int:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            int: The parsed integer value.
        """
        self._single_value_read_validate(tag, value)
        try:
            return int(float(value))
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set value '{value}' to an int for tag '{tag}'!") from err

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        if not self.val_is_within_bounds(value):
            return ""
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        return self._get_token_len()


@dataclass
class FloatTag(AbstractNumericTag):
    """Tag for float values in JDFTx input files.

    Tag for float values in JDFTx input files.
    """

    prec: int | None = None

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value, by default
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        return self._validate_value_type(float, tag, value, try_auto_type_fix=try_auto_type_fix)

    def read(self, tag: str, value: str) -> float:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            float: The parsed float value.
        """
        self._single_value_read_validate(tag, value)
        try:
            value_float = float(value)
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set value '{value}' to a float for tag '{tag}'!") from err
        return value_float

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        if not self.val_is_within_bounds(value):
            return ""
        # pre-convert to string: self.prec+3 is minimum room for:
        # - sign, 1 integer left of decimal, decimal, and precision.
        # larger numbers auto add places to left of decimal
        if self.prec is not None:
            value = f"{value:{self.prec + 3}.{self.prec}f}"
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        return self._get_token_len()


@dataclass
class InitMagMomTag(AbstractTag):
    """Tag for initial-magnetic-moments tag in JDFTx input files.

    Tag for initial-magnetic-moments tag in JDFTx input files.
    """

    # temporary fix to allow use of initial-magnetic-moments tag
    # requires the user to set magnetic moments as a string with no extra
    # validation. Processes input files as simply a string variable with
    # no extra type conversion

    # the formatting of this tag's value depends on the species labels in the
    # ion tags. These species labels are not necessarily element symbols.
    # There are also multiple types of formatting options of the spins
    # most robustly, this should be a MultiFormatTag, with each option
    #     being a StructureDeferredTag, because this tag needs to know the
    #     results of reading in the structure before being able to robustly
    #     parse the value of this tag
    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value, by default
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        return self._validate_value_type(str, tag, value, try_auto_type_fix=try_auto_type_fix)

    def read(self, tag: str, value: str) -> str:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            str: The parsed string value.
        """
        self._general_read_validate(tag, value)
        return str(value)

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        return self._get_token_len()

    def _is_equal_to(self, val1, obj2, val2):
        return True  # TODO: We still need to actually implement initmagmom as a multi-format tag
        # raise NotImplementedError("equality not yet implemented for InitMagMomTag")


@dataclass
class TagContainer(AbstractTag):
    """TagContainer class for handling tags that contain other tags.

    This class is used to handle tags that contain other tags. It is used to validate the type of the value for the tag,
    read the value string for the tag, write the tag and its value as a string, and get the token length of the tag.

    Note: When constructing a TagContainer, all subtags must be able to return the correct token length without any
    information about the value.
    # TODO: Remove this assumption by changing the signature of get_token_len to take the value as an argument.
    """

    linebreak_nth_entry: int | None = None  # handles special formatting for matrix tags, e.g. lattice tag
    is_tag_container: bool = (
        True  # used to ensure only TagContainers are converted between list and dict representations
    )
    subtags: dict[str, AbstractTag] = field(default_factory=dict)

    def _validate_single_entry(
        self, value: dict | list[dict], try_auto_type_fix: bool = False
    ) -> tuple[list[str], list[bool], Any]:
        if not isinstance(value, dict):
            raise TypeError(f"The value '{value}' (of type {type(value)}) must be a dict for this TagContainer!")
        tags_checked: list[str] = []
        types_checks: list[bool] = []
        updated_value = deepcopy(value)
        for subtag, subtag_value in value.items():
            subtag_object = self.subtags[subtag]
            tag, check, subtag_value2 = subtag_object.validate_value_type(
                subtag, subtag_value, try_auto_type_fix=try_auto_type_fix
            )
            if try_auto_type_fix:
                updated_value[subtag] = subtag_value2
            tags_checked.append(tag)
            types_checks.append(check)
        return tags_checked, types_checks, updated_value

    def _validate_bounds_single_entry(self, value: dict | list[dict]) -> tuple[list[str], list[bool], list[str]]:
        if not isinstance(value, dict):
            raise TypeError(f"The value '{value}' (of type {type(value)}) must be a dict for this TagContainer!")
        tags_checked: list[str] = []
        types_checks: list[bool] = []
        reported_errors: list[str] = []
        for subtag, subtag_value in value.items():
            subtag_object = self.subtags[subtag]
            check, err_str = subtag_object.validate_value_bounds(subtag, subtag_value)
            tags_checked.append(subtag)
            types_checks.append(check)
            reported_errors.append(err_str)
        return tags_checked, types_checks, reported_errors

    def validate_value_bounds(self, tag: str, value: Any) -> tuple[bool, str]:
        value_dict = value
        if self.can_repeat:
            self._validate_repeat(tag, value_dict)
            results = [self._validate_bounds_single_entry(x) for x in value_dict]
            tags_list_list: list[list[str]] = [result[0] for result in results]
            is_valids_list_list: list[list[bool]] = [result[1] for result in results]
            reported_errors_list: list[list[str]] = [result[2] for result in results]
            is_valid_out = all(all(x) for x in is_valids_list_list)
            errors_out = ",".join([",".join(x) for x in reported_errors_list])
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for i, x in enumerate(is_valids_list_list):
                    if not all(x):
                        for j, y in enumerate(x):
                            if not y:
                                warnmsg += f"{tags_list_list[i][j]} ({reported_errors_list[i][j]}) "
                warnings.warn(warnmsg, stacklevel=2)
        else:
            tags, is_valids, reported_errors = self._validate_bounds_single_entry(value_dict)
            is_valid_out = all(is_valids)
            errors_out = ",".join(reported_errors)
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for ii, xx in enumerate(is_valids):
                    if not xx:
                        warnmsg += f"{tags[ii]} ({reported_errors[ii]}) "
                warnings.warn(warnmsg, stacklevel=2)
        return is_valid_out, f"{tag}: {errors_out}"

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the type of the value for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value, by default
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        value_dict = self.get_dict_representation(tag, value)
        if self.can_repeat:
            self._validate_repeat(tag, value_dict)
            results = [self._validate_single_entry(x, try_auto_type_fix=try_auto_type_fix) for x in value_dict]
            tags_list_list: list[list[str]] = [result[0] for result in results]
            is_valids_list_list: list[list[bool]] = [result[1] for result in results]
            updated_value: Any = [result[2] for result in results]
            tag_out = ",".join([",".join(x) for x in tags_list_list])
            is_valid_out = all(all(x) for x in is_valids_list_list)
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for i, x in enumerate(is_valids_list_list):
                    if not all(x):
                        for j, y in enumerate(x):
                            if not y:
                                warnmsg += f"{tags_list_list[i][j]} "
                warnings.warn(warnmsg, stacklevel=2)
        else:
            tags, is_valids, updated_value = self._validate_single_entry(
                value_dict, try_auto_type_fix=try_auto_type_fix
            )
            tag_out = ",".join(tags)
            is_valid_out = all(is_valids)
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for ii, xx in enumerate(is_valids):
                    if not xx:
                        warnmsg += f"{tags[ii]} "
                warnings.warn(warnmsg, stacklevel=2)
        return tag_out, is_valid_out, updated_value

    def read(self, tag: str, value: str) -> dict:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            dict: The parsed value.
        """
        self._general_read_validate(tag, value)
        value_list = value.split()

        tempdict = {}  # temporarily store read tags out of order they are processed

        for subtag, subtag_type in (
            (subtag, subtag_type) for subtag, subtag_type in self.subtags.items() if subtag_type.write_tagname
        ):
            # every subtag with write_tagname=True in a TagContainer has a fixed length and can be immediately read in
            # this loop if it is present
            if subtag in value_list:  # this subtag is present in the value string
                subtag_count = value_list.count(subtag)  # Get number of times subtag appears in line
                if not subtag_type.can_repeat:
                    if subtag_count > 1:
                        raise ValueError(
                            f"Subtag '{subtag}' for tag '{tag}' is not allowed to repeat but repeats value {value}"
                        )
                    idx_start = value_list.index(subtag)
                    token_len = subtag_type.get_token_len()
                    idx_end = idx_start + token_len
                    subtag_value = " ".join(
                        value_list[(idx_start + 1) : idx_end]
                    )  # add 1 so the subtag value string excludes the subtagname
                    tempdict[subtag] = subtag_type.read(subtag, subtag_value)
                    del value_list[idx_start:idx_end]
                else:
                    tempdict[subtag] = []
                    for _ in range(subtag_count):
                        idx_start = value_list.index(subtag)
                        idx_end = idx_start + subtag_type.get_token_len()
                        subtag_value = " ".join(
                            value_list[(idx_start + 1) : idx_end]
                        )  # add 1 so the subtag value string excludes the subtagname
                        tempdict[subtag].append(subtag_type.read(subtag, subtag_value))
                        del value_list[idx_start:idx_end]

        # TODO: This breaks for TagContainers with optional StrTags that are not at the end of the line (ie
        # modification subtag for lattice tag). We need to figure out a fix.
        for subtag, subtag_type in (
            (subtag, subtag_type) for subtag, subtag_type in self.subtags.items() if not subtag_type.write_tagname
        ):
            # now try to populate remaining subtags that do not use a keyword in order of appearance.
            if len(value_list) == 0:
                break
            if subtag in tempdict or subtag_type.write_tagname:  # this tag has already been read or requires a tagname
                # keyword to be present
                continue
            tempdict[subtag] = subtag_type.read(subtag, value_list[0])
            del value_list[0]

        # reorder all tags to match order of __MASTER_TAG_LIST__ and do coarse-grained validation of read.

        subdict = {x: tempdict[x] for x in self.subtags if x in tempdict}
        self._check_nonoptional_subtags(tag, subdict, self.subtags)
        self._check_unread_values(tag, value_list)
        return subdict

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.
        """
        if not isinstance(value, dict):
            raise TypeError(
                f"The value '{value}' (of type {type(value)}) for tag '{tag}' must be a dict for this TagContainer!"
            )

        final_value = ""
        indent = "    "
        # Gather all subtag print values into a list
        print_str_list = []
        for subtag, subvalue in value.items():
            if self.subtags[subtag].can_repeat and isinstance(subvalue, list):
                print_str_list += [self.subtags[subtag].write(subtag, entry).strip() + " " for entry in subvalue]
            else:
                print_str_list.append(self.subtags[subtag].write(subtag.strip() + " ", subvalue))
        # Concatenate all subtag print values into a single string, with line breaks at appropriate
        # locations if needed
        for count, print_str in enumerate(print_str_list):
            if self.linebreak_nth_entry is not None:
                # handles special formatting with extra linebreak, e.g. for lattice tag
                i_column = count % self.linebreak_nth_entry
                if i_column == 0:
                    final_value += f"{indent}{print_str}"
                else:
                    final_value += f"{print_str}"
                if i_column == self.linebreak_nth_entry - 1:
                    final_value += "\\\n"
            else:
                final_value += f"{print_str}"
        if self.linebreak_nth_entry is not None:  # handles special formatting for lattice tag
            final_value = final_value[:-2]  # exclude final \\n from final print call

        return self._write(
            tag,
            final_value,
            multiline_override=self.linebreak_nth_entry is not None,
            strip_override=(self.linebreak_nth_entry is not None),
        )

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.
        """
        min_token_len = int(self.write_tagname)  # length of value subtags added next
        for subtag_type in self.subtags.values():
            subtag_token_len = subtag_type.get_token_len()  # recursive for nested TagContainers
            if not subtag_type.optional:  # TagContainers could be longer with optional subtags included
                min_token_len += subtag_token_len
        return min_token_len

    def _make_list(self, value: dict) -> list:
        if not isinstance(value, dict):
            raise TypeError(f"The value '{value}' is not a dict, so could not be converted")
        value_list = []
        for subtag, subtag_value in value.items():
            subtag_type = self.subtags[subtag]
            if subtag_type.allow_list_representation:
                # this block deals with making list representations of any nested TagContainers
                if not isinstance(subtag_value, dict):
                    raise ValueError(f"The subtag {subtag} is not a dict: '{subtag_value}', so could not be converted")
                subtag_value2 = subtag_type.get_list_representation(subtag, subtag_value)  # recursive list generation

                if subtag_type.write_tagname:  # needed to write 'v' subtag in 'ion' tag
                    value_list.append(subtag)
                value_list.extend(subtag_value2)
            elif isinstance(subtag_value, dict):
                # this triggers if someone sets this tag using mixed dict/list representations
                warnings.warn(
                    f"The {subtag} subtag does not allow list representation with a value {subtag_value}.\n "
                    "I added the dict to the list. Is this correct? You will not be able to convert back!",
                    stacklevel=2,
                )
                value_list.append(subtag_value)
            else:
                # the subtag is simply of form {'subtag': subtag_value} and now adds concrete values to the list
                value_list.append(subtag_value)

        # return list of lists for tags in matrix format, e.g. lattice tag
        if (ncol := self.linebreak_nth_entry) is not None:
            nrow = int(len(value_list) / ncol)
            value_list = [[value_list[row * ncol + col] for col in range(ncol)] for row in range(nrow)]
        return value_list

    def get_list_representation(self, tag: str, value: Any) -> list:
        """Convert the value to a list representation.

        Args:
            tag (str): The tag to convert the value to a list representation for.
            value (Any): The value to convert to a list representation.

        Returns:
            list: The value converted to a list representation.
        """
        # convert dict representation into list representation by writing (nested) dicts into list or list of lists.
        # there are 4 types of TagContainers in the list representation:
        # can_repeat: list of bool/str/int/float (ion-species)
        # can_repeat: list of lists (ion)
        # cannot repeat: list of bool/str/int/float (elec-cutoff)
        # cannot repeat: list of lists (lattice)
        if self.can_repeat and not isinstance(value, list):
            raise ValueError(
                f"Value '{value}' must be a list when passed to 'get_list_representation' since "
                f"tag '{tag}' is repeatable."
            )
        if self.can_repeat:
            if all(isinstance(entry, list) for entry in value):
                return value  # no conversion needed
            if any(not isinstance(entry, dict) for entry in value):
                raise ValueError(
                    f"The tag '{tag}' set to value '{value}' must be a list of dicts when passed to "
                    "'get_list_representation' since the tag is repeatable."
                )
            tag_as_list = [self._make_list(entry) for entry in value]
        else:
            tag_as_list = self._make_list(value)
        return tag_as_list

    @staticmethod
    def _check_for_mixed_nesting(tag: str, value: Any) -> None:
        has_nested_dict = any(isinstance(x, dict) for x in value)
        has_nested_list = any(isinstance(x, list) for x in value)
        if has_nested_dict and has_nested_list:
            raise ValueError(
                f"tag '{tag}' with value '{value}' cannot have nested lists/dicts mixed with bool/str/int/floats!"
            )
        if has_nested_dict:
            raise ValueError(
                f"tag '{tag}' with value '{value}' cannot have nested dicts mixed with bool/str/int/floats!"
            )
        if has_nested_list:
            raise ValueError(
                f"tag '{tag}' with value '{value}' cannot have nested lists mixed with bool/str/int/floats!"
            )

    def _make_str_for_dict(self, tag: str, value_list: list) -> str:
        """Convert the value to a string representation.

        Args:
            tag (str): The tag to convert the value to a string representation for.
            value_list (list): The value to convert to a string representation.

        Returns:
            str: The value converted to a string representation for dict representation.
        """
        value = _flatten_list(tag, value_list)
        self._check_for_mixed_nesting(tag, value)
        return " ".join([str(x) for x in value])

    def get_dict_representation(self, tag: str, value: list) -> dict | list[dict]:
        """Convert the value to a dict representation.

        Args:
            tag (str): The tag to convert the value to a dict representation for.
            value (list): The value to convert to a dict representation.

        Returns:
            dict | list[dict]: The value converted to a dict representation.
        """
        # convert list or list of lists representation into string the TagContainer can process back into (nested) dict

        if self.can_repeat and not isinstance(value, list):
            raise ValueError(
                f"Value '{value}' must be a list when passed to 'get_dict_representation' since "
                f"tag '{tag}' is repeatable."
            )
        if (
            self.can_repeat and len({len(x) for x in value}) > 1
        ):  # Creates a list of every unique length of the subdicts
            # TODO: Populate subdicts with fewer entries with JDFTx defaults to make compatible
            raise ValueError(f"The values '{value}' for tag '{tag}' provided in a list of lists have different lengths")
        value = value.tolist() if isinstance(value, np.ndarray) else value

        # there are 4 types of TagContainers in the list representation:
        # can_repeat: list of bool/str/int/float (ion-species)
        # can_repeat: list of lists (ion)
        # cannot repeat: list of bool/str/int/float (elec-cutoff)
        # cannot repeat: list of lists (lattice)

        # the .read() method automatically handles regenerating any nesting because is just like reading a file
        if self.can_repeat:
            if all(isinstance(entry, dict) for entry in value):
                return value  # no conversion needed
            string_value = [self._make_str_for_dict(tag, entry) for entry in value]
            return [self.read(tag, entry) for entry in string_value]

        if isinstance(value, dict):
            return value  # no conversion needed
        list_value = self._make_str_for_dict(tag, value)
        return self.read(tag, list_value)

    def _is_equal_to(self, val1, obj2, val2):
        """Check if the two values are equal.

        Return False if (checked in following order)
        - obj2 is not a TagContainer
        - all of val1's subtags are not in val2
        - val1 and val2 are not the same length (different number of subtags)
        - at least one subtag in val1 is not equal to the corresponding subtag in val2
        """
        if self._is_same_tagtype(obj2):
            if isinstance(val1, dict) and isinstance(val2, dict):
                if all(subtag in val2 for subtag in val1) and (len(list(val1.keys())) == len(list(val2.keys()))):
                    for subtag, subtag_type in self.subtags.items():
                        if (subtag in val1) and (
                            not subtag_type.is_equal_to(val1[subtag], obj2.subtags[subtag], val2[subtag])
                        ):
                            return False
                    return True
                return False
            raise ValueError("Values must be in dictionary format for TagContainer comparison")
        return False


# TODO: Write StructureDefferedTagContainer back in (commented out code block removed
# on 11/4/24) and make usable for tags like initial-magnetic-moments


@dataclass
class MultiformatTag(AbstractTag):
    """Class for tags with multiple format options.

    Class for tags that could have different types of input values given to them or tags where different subtag options
    directly impact how many expected arguments are provided e.g. the coulomb-truncation or van-der-waals tags.

    This class should not be used for tags with simply some combination of mandatory and optional args because the
    TagContainer class can handle those cases by itself.
    """

    format_options: list[AbstractTag] = field(default_factory=list)

    def validate_value_type(self, tag: str, value: Any, try_auto_type_fix: bool = False) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Args:
            tag (str): The tag to validate the value type for.
            value (Any): The value to validate the type of.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value. Defaults to
                False.

        Returns:
            tuple[str, bool, Any]: The tag, whether the value is of the correct type, and the possibly fixed value.
        """
        format_index, value = self._determine_format_option(tag, value, try_auto_type_fix=try_auto_type_fix)
        is_valid = format_index is not None
        return tag, is_valid, value

    def read(self, tag: str, value: str) -> None:
        """Read the value string for this tag.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Raises:
            RuntimeError: If the method is called directly on MultiformatTag.
        """
        err_str = "The read method is not supposed to be called directly on MultiformatTag."
        "Get the proper format option first."
        raise RuntimeError(err_str)

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Args:
            tag (str): The tag to write.
            value (Any): The value to write.

        Returns:
            str: The tag and its value as a string.

        Raises:
            RuntimeError: If the method is called directly on MultiformatTag.
        """
        err_str = "The write method is not supposed to be called directly on MultiformatTag."
        "Get the proper format option first."
        raise RuntimeError(err_str)

    def get_format_index_for_str_value(self, tag: str, value: str) -> int:
        """Get the format index from string representation of value.

        Args:
            tag (str): The tag to read the value string for.
            value (str): The value string to read.

        Returns:
            int: The index of the format option for the value of this tag.

        Raises:
            ValueError: If no valid read format is found for the tag.
        """
        problem_log = []
        for i, trial_format in enumerate(self.format_options):
            try:
                _ = trial_format.read(tag, value)
                return i
            except (ValueError, TypeError) as e:
                problem_log.append(f"Format {i}: {e}")
        raise ValueError(
            f"No valid read format for tag '{tag}' with value '{value}'\n"
            "Add option to format_options or double-check the value string and retry!\n\n"
        )

    def raise_invalid_format_option_error(self, tag: str, i: int) -> None:
        """Raise an error for an invalid format option.

        Args:
            tag (str): The tag to raise the error for.
            i (int): The index of the format option to raise the error for.

        Raises:
            ValueError: If the format option is invalid.
        """
        raise ValueError(f"tag '{tag}' failed to validate for option {i}")

    def _determine_format_option(self, tag: str, value_any: Any, try_auto_type_fix: bool = False) -> tuple[int, Any]:
        """Determine the format option for the value of this tag.

        Args:
            tag (str): The tag to determine the format option for.
            value_any (Any): The value to determine the format option for.
            try_auto_type_fix (bool, optional): Whether to try to automatically fix the type of the value. Defaults to
                False.

        Returns:
            tuple[int, Any]: The index of the format option for the value of this tag and the value of this tag.

        Raises:
            ValueError: If the format for the tag could not be determined from the available options.
        """
        exceptions = []
        for i, format_option in enumerate(self.format_options):
            if format_option.can_repeat:
                value = [value_any] if not isinstance(value_any, list) else value_any
            else:
                value = value_any
            try:
                _, is_tag_valid, value = format_option.validate_value_type(
                    tag, value, try_auto_type_fix=try_auto_type_fix
                )
                if not is_tag_valid:
                    self.raise_invalid_format_option_error(tag, i)
                else:
                    return i, value
            except (ValueError, TypeError, KeyError) as e:
                exceptions.append(e)
        raise ValueError(
            f"The format for tag '{tag}' with value '{value_any}' could not be determined from the available options! "
            "Check your inputs and/or MASTER_TAG_LIST!"
        )

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Returns:
            int: The token length of the tag.

        Raises:
            NotImplementedError: If the method is called directly on MultiformatTag objects.
        """
        raise NotImplementedError("This method is not supposed to be called directly on MultiformatTag objects!")

    def _is_equal_to(self, val1, obj2, val2):
        raise NotImplementedError("This method is not supposed to be called directly on MultiformatTag objects!")


@dataclass
class BoolTagContainer(TagContainer):
    """BoolTagContainer class for handling the subtags to the "dump" tag.

    This class is used to handle the subtags to the "dump" tag. All subtags are freqs for dump, and all values for these
    tags are boolean values that are read given the existence of their "var" name.
    """

    def read(self, tag: str, value_str: str) -> dict:
        """Read the value string for this tag.

        This method reads the value string for this tag. It is used to parse the value string for the tag and return the
        parsed value.

        Args:
            tag (str): The tag to read the value string for.
            value_str (str): The value string to read.

        Returns:
            dict: The parsed value.
        """
        self._general_read_validate(tag, value_str)
        value = value_str.split()
        tempdict = {}
        for subtag, subtag_type in self.subtags.items():
            if subtag in value:
                idx_start = value.index(subtag)
                idx_end = idx_start + subtag_type.get_token_len()
                subtag_value = " ".join(value[(idx_start + 1) : idx_end])
                tempdict[subtag] = subtag_type.read(subtag, subtag_value)
                del value[idx_start:idx_end]
        subdict = {x: tempdict[x] for x in self.subtags if x in tempdict}
        self._check_nonoptional_subtags(tag, subdict, self.subtags)
        self._check_unread_values(tag, value)
        return subdict


@dataclass
class DumpTagContainer(TagContainer):
    """DumpTagContainer class for handling the "dump" tag.

    This class is used to handle the "dump" tag.
    """

    def read(self, tag: str, value_str: str) -> dict:
        """Read the value string for this tag.

        This method reads the value string for this tag. It is used to parse the value string for the tag and return the
        parsed value.

        Args:
            tag (str): The tag to read the value string for.
            value_str (str): The value string to read.

        Returns:
            dict: The parsed value.
        """
        self._general_read_validate(tag, value_str)
        value = value_str.split()
        tempdict = {}
        # Each subtag is a freq, which will be a BoolTagContainer
        for subtag, subtag_type in self.subtags.items():
            if subtag in value:
                idx_start = value.index(subtag)
                subtag_value = " ".join(value[(idx_start + 1) :])
                tempdict[subtag] = subtag_type.read(subtag, subtag_value)
                del value[idx_start:]
        # reorder all tags to match order of __MASTER_TAG_LIST__ and do coarse-grained validation of read
        subdict = {x: tempdict[x] for x in self.subtags if x in tempdict}
        # There are no forced subtags for dump
        self._check_unread_values(tag, value)
        return subdict


def _flatten_list(tag: str, list_of_lists: list[Any]) -> list[Any]:
    """Flatten list of lists into a single list, then stop.

    Flatten list of lists into a single list, then stop.

    Parameters
    ----------
    tag : str
        The tag to flatten the list of lists for.
    list_of_lists : list[Any]
        The list of lists to flatten.

    Returns
    -------
    list[Any]
        The flattened list.
    """
    if not isinstance(list_of_lists, list):
        raise TypeError(f"{tag}: You must provide a list to flatten_list()!")
    flist = []
    for v in list_of_lists:
        if isinstance(v, list):
            flist.extend(_flatten_list(tag, v))
        else:
            flist.append(v)
    return flist
