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
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from pymatgen.core import Structure

__author__ = "Jacob Clary"


def flatten_list(tag: str, list_of_lists: list[Any]) -> list[Any]:
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
            flist.extend(flatten_list(tag, v))
        else:
            flist.append(v)
    return flist


class ClassPrintFormatter:
    """Generic class for printing to command line in readable format.

    Generic class for printing to command line in readable format.
    """

    def __str__(self) -> str:
        """Print the class to the command line in a readable format.

        Print the class to the command line in a readable format.

        Returns
        -------
        str
            The class in a readable format.
        """
        return f"{self.__class__}\n" + "\n".join(
            f"{item} = {self.__dict__[item]}" for item in sorted(self.__dict__)
        )


@dataclass(kw_only=True)
class AbstractTag(ClassPrintFormatter, ABC):
    """Abstract base class for all tags."""

    multiline_tag: bool = False  # set to True if what to print tags across
    # multiple lines, typically like electronic-minimize
    can_repeat: bool = (
        False  # set to True for tags that can appear on multiple lines, like ion
    )
    write_tagname: bool = (
        True  # set to False to not print the tagname, like for subtags of
        # elec-cutoff
    )
    write_value: bool = (
        True  # set to False to not print any value, like for dump-interval
    )
    optional: bool = True  # set to False if tag (usually a subtag of a
    # TagContainer) must be set for the JDFTXInfile to be valid.
    # The lattice, ion, and ion-species are the main tags that are not optional
    defer_until_struc: bool = False
    is_tag_container: bool = False
    allow_list_representation: bool = (
        False  # if True, allow this tag to exist as a list or list of lists
    )

    @abstractmethod
    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag."""

    def _validate_value_type(
        self, type_check: type, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        This method is used to validate the type of the value for this tag. It is
        used to check if the value is of the correct type and if it can be fixed
        automatically. If the value is not of the correct type and cannot be fixed
        automatically, a warning is raised.

        Parameters
        ----------
        type_check : type
            The type to check the value against.
        tag : str
            The tag to check the value against.
        value : Any
            The value to check the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to check the value against.
        is_valid : bool
            Whether the value is of the correct type.
        value : Any
            The value checked against the correct type, possibly fixed.
        """
        if self.can_repeat:
            self._validate_repeat(tag, value)
            is_valid = all(isinstance(x, type_check) for x in value)
        else:
            is_valid = isinstance(value, type_check)

        if not is_valid and try_auto_type_fix:
            try:
                if self.can_repeat:
                    value = [self.read(tag, str(x)) for x in value]
                else:
                    value = self.read(tag, str(value))
                tag, is_valid, value = self._validate_value_type(type_check, tag, value)
            except (TypeError, ValueError):
                warnings.warn(
                    f"Could not fix the typing for {tag} {value}!", stacklevel=2
                )
        return tag, is_valid, value

    def _validate_repeat(self, tag: str, value: Any) -> None:
        if not isinstance(value, list):
            raise TypeError(f"The {tag} tag can repeat but is not a list: {value}")

    @abstractmethod
    def read(self, tag: str, value_str: str) -> Any:
        """Read and parse the value string for this tag."""

    @abstractmethod
    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string."""

    @abstractmethod
    def get_token_len(self) -> int:
        """Get the token length of the tag."""

    def _write(self, tag: str, value: Any, multiline_override: bool = False) -> str:
        tag_str = f"{tag} " if self.write_tagname else ""
        if self.multiline_tag or multiline_override:
            tag_str += "\\\n"
        if self.write_value:
            tag_str += f"{value} "
        return tag_str

    def _get_token_len(self) -> int:
        return int(self.write_tagname) + int(self.write_value)

    def get_list_representation(self, tag: str, value: Any) -> list | list[list]:
        """Convert the value to a list representation."""
        raise ValueError(f"Tag object has no get_list_representation method: \
                         {tag}")

    def get_dict_representation(self, tag: str, value: Any) -> dict | list[dict]:
        """Convert the value to a dict representation."""
        raise ValueError(f"Tag object has no get_dict_representation method: \
                         {tag}")


"""
TODO:
fix dump-name and density-of-states tags

check that all ions either have or lack velocities
add validation of which tags require/forbid presence of other tags according to
JDFTx docs?

choose how DeferredTags inherit from TagContainer?? same functionality once
process the values "for real"

#possible TODO: add defaults like JDFTx does

MISC TODO:
    note which tags I've enforced a mandatory formatting,
        1. dump-name allows only 1 format
        2. debug requires 1 entry per line
        3. dump requires 1 entry per line

"""


@dataclass(kw_only=True)
class BoolTag(AbstractTag):
    """Tag for boolean values in JDFTx input files.

    Tag for boolean values in JDFTx input files.
    """

    _TF_read_options: dict[str, bool] = field(
        default_factory=lambda: {"yes": True, "no": False}
    )
    _TF_write_options: dict[bool, str] = field(
        default_factory=lambda: {True: "yes", False: "no"}
    )
    _TF_options: dict[str, dict] = field(init=False)

    def __post_init__(self) -> None:
        """Initialize the _TF_options attribute.

        Initialize the _TF_options attribute.
        """
        self._TF_options = {
            "read": self._TF_read_options,
            "write": self._TF_write_options,
        }

    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        return self._validate_value_type(
            bool, tag, value, try_auto_type_fix=try_auto_type_fix
        )

    def raise_value_error(self, tag: str, value: str) -> None:
        """Raise a ValueError for the value string.

        Raise a ValueError for the value string.

        Parameters
        ----------
        tag : str
            The tag to raise the ValueError for.
        value : str
            The value string to raise the ValueError for.
        """
        raise ValueError(
            f"The value '{value}' was provided to {tag}, it is not \
                acting like a boolean"
        )

    def read(self, tag: str, value: str) -> bool:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        bool
        """
        if len(value.split()) > 1:
            raise ValueError(f"'{value}' for {tag} should not have a space in it!")
        try:
            if not self.write_value:
                # accounts for exceptions where only the tagname is used, e.g.
                # dump-only or dump-fermi-density (sometimes) tags
                if not value:  # then the string '' was passed in because no
                    # value was provided but the tag was present
                    value = "yes"
                else:
                    self.raise_value_error(tag, value)
            return self._TF_options["read"][value]
        except (ValueError, TypeError) as err:
            raise ValueError(
                f"Could not set '{value}' as True/False for {tag}!"
            ) from err

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
        """
        value2 = self._TF_options["write"][value]
        return self._write(tag, value2)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        return self._get_token_len()


@dataclass(kw_only=True)
class StrTag(AbstractTag):
    """Tag for string values in JDFTx input files.

    Tag for string values in JDFTx input files.
    """

    options: list = None

    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        return self._validate_value_type(
            str, tag, value, try_auto_type_fix=try_auto_type_fix
        )

    def read(self, tag: str, value: str) -> str:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        str
        """
        if len(value.split()) > 1:
            raise ValueError(f"'{value}' for {tag} should not have a space in it!")
        try:
            value = str(value)
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set '{value}' to a str for {tag}!") from err
        if self.options is None or value in self.options:
            return value
        raise ValueError(
            f"The '{value}' string must be one of {self.options} for {tag}"
        )

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
            The tag and its value as a string.
        """
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        return self._get_token_len()


@dataclass(kw_only=True)
class IntTag(AbstractTag):
    """Tag for integer values in JDFTx input files.

    Tag for integer values in JDFTx input files.
    """

    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        return self._validate_value_type(
            int, tag, value, try_auto_type_fix=try_auto_type_fix
        )

    def read(self, tag: str, value: str) -> int:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        int
        """
        if len(value.split()) > 1:
            raise ValueError(f"'{value}' for {tag} should not have a space in it!")
        try:
            return int(float(value))
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set '{value}' to a int for {tag}!") from err

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
            The tag and its value as a string.
        """
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        return self._get_token_len()


@dataclass(kw_only=True)
class FloatTag(AbstractTag):
    """Tag for float values in JDFTx input files.

    Tag for float values in JDFTx input files.
    """

    prec: int = field(default=None)

    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        return self._validate_value_type(
            float, tag, value, try_auto_type_fix=try_auto_type_fix
        )

    def read(self, tag: str, value: str) -> float:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        float
        """
        if len(value.split()) > 1:
            raise ValueError(f"'{value}' for {tag} should not have a space in it!")
        try:
            value_float = float(value)
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set '{value}' to a float for {tag}!") from err
        return value_float

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
            The tag and its value as a string.
        """
        # pre-convert to string: self.prec+3 is minimum room for:
        # - sign, 1 integer left of decimal, decimal, and precision.
        # larger numbers auto add places to left of decimal
        if self.prec is not None:
            value = f"{value:{self.prec+3}.{self.prec}f}"
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        return self._get_token_len()


@dataclass(kw_only=True)
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
    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        return self._validate_value_type(
            str, tag, value, try_auto_type_fix=try_auto_type_fix
        )

    def read(self, tag: str, value: str) -> str:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        str
        """
        try:
            value = str(value)
        except (ValueError, TypeError) as err:
            raise ValueError(f"Could not set '{value}' to a str for {tag}!") from err
        return value

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
            The tag and its value as a string.
        """
        return self._write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        return self._get_token_len()


@dataclass(kw_only=True)
class TagContainer(AbstractTag):
    """TagContainer class for handling tags that contain other tags.

    This class is used to handle tags that contain other tags. It is used to
    validate the type of the value for the tag, read the value string for the tag,
    write the tag and its value as a string, and get the token length of the tag.

    """

    is_tag_container: bool = True  # used to ensure only TagContainers are
    # converted between list and dict representations
    subtags: dict[str, AbstractTag] = None
    linebreak_nth_entry: int = (
        None  # handles special formatting for matrix tags, e.g. lattice tag
    )

    def _validate_single_entry(
        self, value: dict | list[dict], try_auto_type_fix: bool = False
    ) -> tuple[list[str], list[bool], Any]:
        if not isinstance(value, dict):
            raise TypeError(f"This tag should be a dict: {value}, which is \
                             of the type {type(value)}")
        tags_checked = []
        types_checks = []
        updated_value = deepcopy(value)
        for subtag, subtag_value in value.items():
            subtag_object = self.subtags[subtag]
            tags, checks, subtag_value2 = subtag_object.validate_value_type(
                subtag, subtag_value, try_auto_type_fix=try_auto_type_fix
            )
            if try_auto_type_fix:
                updated_value[subtag] = subtag_value2
            if isinstance(checks, list):
                tags_checked.extend(tags)
                types_checks.extend(checks)
            else:
                tags_checked.append(tags)
                types_checks.append(checks)
        return tags_checked, types_checks, updated_value

    # TODO: This method violates the return signature of the AbstractTag
    # class's validate_value_type (tuple[str, bool Any]). There are enough type
    # checks in all the functions that call validate_value_type to make sure
    # this doesn't break anything, but regardless pre-commit does not allow it.
    # (TODO action is to fix the return signature and make sure it doesn't break)
    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the type of the value for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        tag : str
            The tag to validate the type of the value for.
        is_valid_out : bool
            Whether the value is of the correct type.
        updated_value : Any
            The value checked against the correct type, possibly fixed.
        """
        value_dict = self.get_dict_representation(tag, value)
        if self.can_repeat:
            self._validate_repeat(tag, value_dict)
            results = [
                self._validate_single_entry(x, try_auto_type_fix=try_auto_type_fix)
                for x in value_dict
            ]
            tags, is_valids, updated_value = [list(x) for x in list(zip(*results))]
            tag_out = ",".join([",".join(x) for x in tags])
            is_valid_out = all(all(x) for x in is_valids)
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for i, x in enumerate(is_valids):
                    if not all(x):
                        for j, y in enumerate(x):
                            if not y:
                                warnmsg += f"{tags[i][j]} "

        else:
            tags, is_valids, updated_value = self._validate_single_entry(
                value_dict, try_auto_type_fix=try_auto_type_fix
            )
            tag_out = ",".join(tags)
            is_valid_out = all(is_valids)
            if not is_valid_out:
                warnmsg = "Invalid value(s) found for: "
                for i, x in enumerate(is_valids):
                    if not x:
                        warnmsg += f"{tags[i]} "
        return tag_out, is_valid_out, updated_value

    # def read(self, tag: str, value: str | list | dict) -> dict:
    def read(self, tag: str, value: str) -> dict:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.

        Returns
        -------
        dict
        """
        value_list = value.split()
        if tag == "ion":
            special_constraints = [
                x in ["HyperPlane", "Linear", "None", "Planar"] for x in value_list
            ]
            if any(special_constraints):
                value_list = value_list[: special_constraints.index(True)]
                warnings.warn(
                    "Found special constraints reading an 'ion' tag, \
                        these were dropped; reading them has not been \
                            implemented!",
                    stacklevel=2,
                )

        tempdict = {}  # temporarily store read tags out of order they are
        # processed
        # Could this for loop only loop over subtags that have
        # write_tagname=True? This doesn't throw errors but seems dangerous.
        for subtag, subtag_type in self.subtags.items():
            # every subtag with write_tagname=True in a TagContainer has a
            # fixed length and can be immediately read in this loop if it is
            # present
            if subtag in value_list:  # this subtag is present in the value string
                subtag_count = value_list.count(subtag)
                if not subtag_type.can_repeat:
                    if subtag_count > 1:
                        raise ValueError(
                            f"Subtag {subtag} is not allowed to repeat but \
                                appears more than once in {tag}'s value {value}"
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
                        idx_start = value.index(subtag)
                        idx_end = idx_start + subtag_type.get_token_len()
                        subtag_value = " ".join(
                            value_list[(idx_start + 1) : idx_end]
                        )  # add 1 so the subtag value string excludes the
                        # subtagname
                        tempdict[subtag].append(subtag_type.read(subtag, subtag_value))
                        del value_list[idx_start:idx_end]

        for subtag, subtag_type in self.subtags.items():
            # now try to populate remaining subtags that do not use a keyword
            # in order of appearance. Since all subtags in JDFTx that are
            # TagContainers use a keyword to start their field, we know that
            # any subtags processed here are only populated with a single token.
            if len(value_list) == 0:
                break
            if (
                subtag in tempdict or subtag_type.write_tagname
            ):  # this tag has already been read or requires a tagname keyword
                # to be present
                continue
            # note that this next line breaks if the JDFTx dump-name formatting
            # is allowing dump-name would have nested repeating TagContainers,
            # which each need 2 values. You could check for which nonoptional
            # args the TagContainers need and provide those but that's not
            # general. You really need to be passing the entire value string
            # for parsing, but that changes the return args.
            tempdict[subtag] = subtag_type.read(subtag, value_list[0])
            del value_list[0]

        # reorder all tags to match order of __MASTER_TAG_LIST__ and do
        # coarse-grained validation of read.
        subdict = {x: tempdict[x] for x in self.subtags if x in tempdict}
        for subtag, subtag_type in self.subtags.items():
            if not subtag_type.optional and subtag not in subdict:
                raise ValueError(
                    f"The {subtag} tag is not optional but was not populated \
                        during the read!"
                )
        if len(value_list) > 0:
            raise ValueError(
                f"Something is wrong in the JDFTXInfile formatting, some \
                    values were not processed: {value}"
            )
        return subdict

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        Write the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
            The tag and its value as a string.
        """
        if not isinstance(value, dict):
            raise TypeError(
                f"value = {value}\nThe value to the {tag} write method must be \
                    a dict since it is a TagContainer!"
            )

        final_value = ""
        indent = "    "
        for count, (subtag, subvalue) in enumerate(value.items()):
            if self.subtags[subtag].can_repeat and isinstance(subvalue, list):
                # if a subtag.can_repeat, it is assumed that subvalue is a list
                #    the 2nd condition ensures this
                # if it is not a list, then the tag will still be printed by the else
                #    this could be relevant if someone manually sets the tag's
                #    can_repeat value to a non-list.
                print_str_list = [
                    self.subtags[subtag].write(subtag, entry) for entry in subvalue
                ]
                print_str = " ".join(print_str_list)
            else:
                print_str = self.subtags[subtag].write(subtag, subvalue)

            if self.multiline_tag:
                final_value += f"{indent}{print_str}\\\n"
            elif self.linebreak_nth_entry is not None:
                # handles special formatting with extra linebreak, e.g. for lattice tag
                i_column = count % self.linebreak_nth_entry
                if i_column == 1:
                    final_value += f"{indent}{print_str}"
                elif i_column == 0:
                    final_value += f"{print_str}\\\n"
                else:
                    final_value += f"{print_str}"
            else:
                final_value += f"{print_str}"
        if (
            self.multiline_tag or self.linebreak_nth_entry is not None
        ):  # handles special formatting for lattice tag
            final_value = final_value[:-2]  # exclude final \\n from final
            # print call

        return self._write(tag, final_value, self.linebreak_nth_entry is not None)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        min_token_len = int(self.write_tagname)  # length of value subtags
        # added next
        for subtag_type in self.subtags.values():
            subtag_token_len = (
                subtag_type.get_token_len()
            )  # recursive for nested TagContainers
            if (
                not subtag_type.optional
            ):  # TagContainers could be longer with optional subtags included
                min_token_len += subtag_token_len
        return min_token_len

    def check_representation(self, tag: str, value: Any) -> str:
        """Check the representation of the value.

        Check the representation of the value.

        Parameters
        ----------
        tag : str
            The tag to check the representation of the value for.
        value : Any
            The value to check the representation of.

        Returns
        -------
        str
            The representation of the value.
        """
        if not self.allow_list_representation:
            return "dict"
        value_list = self.get_list_representation(tag, value)
        value_dict = self.get_dict_representation(tag, value)
        if value == value_list:
            return "list"
        if value == value_dict:
            return "dict"
        raise ValueError(
            "Could not determine TagContainer representation, something is wrong"
        )

    def _make_list(self, value: dict) -> list:
        value_list = []
        for subtag in value:
            subtag_type = self.subtags[subtag]
            if subtag_type.allow_list_representation:
                # this block deals with making list representations of any
                # nested TagContainers
                if not isinstance(value[subtag], dict):
                    raise ValueError(
                        f"The subtag {subtag} is not a dict: '{value[subtag]}',\
                             so could not be converted"
                    )
                subtag_value2 = subtag_type.get_list_representation(
                    subtag, value[subtag]
                )  # recursive list generation

                if subtag_type.write_tagname:  # needed to write 'v' subtag in
                    # 'ion' tag
                    value_list.append(subtag)
                value_list.extend(subtag_value2)
            elif not subtag_type.allow_list_representation and isinstance(
                value[subtag], dict
            ):
                # this triggers if someone sets this tag using mixed dict/list
                # representations
                warnings.warn(
                    f"The {subtag} subtag does not allow list \
                              representation with a value {value[subtag]}.\n \
                              I added the dict to the list. Is this correct? \
                                You will not be able to convert back!",
                    stacklevel=2,
                )
                value_list.append(value[subtag])
            else:
                # the subtag is simply of form {'subtag': subtag_value} and now
                # adds concrete values to the list
                value_list.append(value[subtag])

        # return list of lists for tags in matrix format, e.g. lattice tag
        if (ncol := self.linebreak_nth_entry) is not None:
            nrow = int(len(value_list) / ncol)
            value_list = [
                [value_list[row * ncol + col] for col in range(ncol)]
                for row in range(nrow)
            ]
        return value_list

    def get_list_representation(self, tag: str, value: Any) -> list:
        """Convert the value to a list representation.

        Convert the value to a list representation.

        Parameters
        ----------
        tag : str
            The tag to convert the value to a list representation for.
        value : Any
            The value to convert to a list representation.

        Returns
        -------
        list
            The value converted to a list representation.
        """
        # convert dict representation into list representation by writing
        # (nested) dicts into list or list of lists.
        # there are 4 types of TagContainers in the list representation:
        # can_repeat: list of bool/str/int/float (ion-species)
        # can_repeat: list of lists (ion)
        # cannot repeat: list of bool/str/int/float (elec-cutoff)
        # cannot repeat: list of lists (lattice)
        if self.can_repeat:
            if all(isinstance(entry, list) for entry in value):
                return value  # no conversion needed
            if any(not isinstance(entry, dict) for entry in value):
                raise ValueError(f"The {tag} tag set to {value} must be a list \
                                 of dict")
            tag_as_list = [self._make_list(entry) for entry in value]
        else:
            tag_as_list = self._make_list(value)
        return tag_as_list

    @staticmethod
    def _check_for_mixed_nesting(tag: str, value: Any) -> None:
        if any(isinstance(x, (dict, list)) for x in value):
            raise ValueError(
                f"{tag} with {value} cannot have nested lists/dicts mixed with \
                    bool/str/int/floats!"
            )

    def _make_str_for_dict(self, tag: str, value_list: list) -> str:
        """Convert the value to a string representation.

        Convert the value to a string representation for dict representation.

        Parameters
        ----------
        tag : str
            The tag to convert the value to a string representation for.
        value_list : list
            The value to convert to a string representation.

        Returns
        -------
        str
            The value converted to a string representation for dict
            representation.
        """
        value = flatten_list(tag, value_list)
        self._check_for_mixed_nesting(tag, value)
        return " ".join([str(x) for x in value])

    def get_dict_representation(self, tag: str, value: list) -> dict | list[dict]:
        """Convert the value to a dict representation.

        Convert the value to a dict representation.

        Parameters
        ----------
        tag : str
            The tag to convert the value to a dict representation for.
        value : list
            The value to convert to a dict representation.

        Returns
        -------
        dict | list[dict]
            The value converted to a dict representation.
        """
        # convert list or list of lists representation into string the
        # TagContainer can process back into (nested) dict
        if self.can_repeat and len({len(x) for x in value}) > 1:  # repeated
            # tags must be in same format
            raise ValueError(
                f"The values for {tag} {value} provided in a list of lists \
                    have different lengths"
            )
        value = value.tolist() if isinstance(value, np.ndarray) else value

        # there are 4 types of TagContainers in the list representation:
        # can_repeat: list of bool/str/int/float (ion-species)
        # can_repeat: list of lists (ion)
        # cannot repeat: list of bool/str/int/float (elec-cutoff)
        # cannot repeat: list of lists (lattice)

        # the .read() method automatically handles regenerating any nesting
        # because is just like reading a file
        if self.can_repeat:
            if all(isinstance(entry, dict) for entry in value):
                return value  # no conversion needed
            string_value = [self._make_str_for_dict(tag, entry) for entry in value]
            return [self.read(tag, entry) for entry in string_value]

        if isinstance(value, dict):
            return value  # no conversion needed
        list_value = self._make_str_for_dict(tag, value)
        return self.read(tag, list_value)


@dataclass(kw_only=True)
class StructureDeferredTagContainer(TagContainer):
    """Class for tags that require a Pymatgen structure to process the value.

    This tag class accommodates tags that can have complicated values that
    depend on the number and species of atoms present. The species labels do
    not necessarily have to be elements, but just match the species given in
    the ion/ion-species tag(s). We will use the set of labels provided by the
    ion tag(s) because that is a well-defined token, while it may not be
    explicitly defined in ion-species.

    Relevant tags: add-U, initial-magnetic-moments, initial-oxidation-states,
    set-atomic-radius, setVDW
    """

    defer_until_struc: bool = True

    def read(self, tag: str, value: str, structure: Structure = None) -> None:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.
        structure : Structure, optional
            The Pymatgen structure to use for reading the value string,
            by default None.
        """
        raise NotImplementedError

        # """This method is similar to StrTag.read(), but with less validation
        # because usually will
        # get a string like 'Fe 2.0 2.5 Ni 1.0 1.1' as the value to process later

        # If this method is called separately from the JDFTXInfile processing
        # methods, a Pymatgen
        # structure may be provided directly
        # """
        # try:
        #     value = str(value)
        # except:
        #     raise ValueError(f"Could not set '{value}' to a str for {tag}!")

        # if structure is not None:
        #     value = self.read_with_structure(tag, value, structure)
        # return value

    def read_with_structure(self, tag: str, value: str, structure: Structure) -> None:
        """Read tag/value pair with a Pymatgen structure provided.

        Read tag/value pair with a Pymatgen structure provided.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.
        structure : Structure
            The Pymatgen structure to use for reading the value string.
        """
        raise NotImplementedError

        # """Fully process the value string using data from the Pymatgen
        # structure"""
        # return self._TC_read(tag, value, structure)


@dataclass(kw_only=True)
class MultiformatTag(AbstractTag):
    """Class for tags with multiple format options.

    Class for tags that could have different types of
    input values given to them or tags where different subtag options directly
    impact how many expected arguments are provided e.g. the coulomb-truncation
    or van-der-waals tags.

    This class should not be used for tags with simply some combination of
    mandatory and optional args because the TagContainer class can handle those
    cases by itself.
    """

    format_options: list[AbstractTag] = None

    def validate_value_type(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[str, bool, Any]:
        """Validate the type of the value for this tag.

        Validate the type of the value for this tag.

        Parameters
        ----------
        tag : str
            The tag to validate the value type for.
        value : Any
            The value to validate the type of.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.
        """
        format_index, value = self._determine_format_option(
            tag, value, try_auto_type_fix=try_auto_type_fix
        )
        is_valid = format_index is not None
        return tag, is_valid, value

    def read(self, tag: str, value: str) -> None:
        """Read the value string for this tag.

        Read the value string for this tag.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value : str
            The value string to read.
        """
        problem_log = []
        for i, trial_format in enumerate(
            self.format_options
        ):  # format_options is a list of AbstractTag-inheriting objects
            try:
                return trial_format.read(tag, value)
            except (ValueError, TypeError) as e:
                problem_log.append(f"Format {i}: {e}")
        errormsg = f"No valid read format for '{tag} {value}' tag\nAdd option \
            to format_options or double-check the value string and retry!\n\n"
        errormsg += "Here is the log of errors for each known \
            formatting option:\n"
        errormsg += "\n".join(
            [f"Format {x}: {problem_log[x]}" for x in range(len(problem_log))]
        )
        raise ValueError(errormsg)

    def raise_invalid_format_option_error(self, tag: str, i: int) -> None:
        """Raise an error for an invalid format option.

        Raise an error for an invalid format option.

        Parameters
        ----------
        tag : str
            The tag to raise the error for.
        i : int
            The index of the format option to raise the error for.
        """
        raise ValueError(
            f"{tag} option {i} is not it: validation \
                                     failed"
        )

    def _determine_format_option(
        self, tag: str, value: Any, try_auto_type_fix: bool = False
    ) -> tuple[int, Any]:
        """Determine the format option for the value of this tag.

        This method determines the format option for the value of this tag.

        Parameters
        ----------
        tag : str
            The tag to determine the format option for.
        value : Any
            The value to determine the format option for.
        try_auto_type_fix : bool, optional
            Whether to try to automatically fix the type of the value,
            by default False.

        Returns
        -------
        int
            The index of the format option for the value of this tag.
        Any
            The value of this tag.
        """
        exceptions = []
        for i, format_option in enumerate(self.format_options):
            try:
                _, is_tag_valid, value = format_option.validate_value_type(
                    tag, value, try_auto_type_fix=try_auto_type_fix
                )
                if not is_tag_valid:
                    self.raise_invalid_format_option_error(tag, i)
                else:
                    return i, value
            except (ValueError, TypeError) as e:  # TODO: Make sure these are all
                # the possible exceptions
                exceptions.append(e)
        raise ValueError(
            f"The format for {tag} for:\n{value}\ncould not be determined \
                from the available options! Check your inputs and/or \
                    MASTER_TAG_LIST! (exceptions: {exceptions})"
        )

    def write(self, tag: str, value: Any) -> str:
        """Write the tag and its value as a string.

        This method writes the tag and its value as a string.

        Parameters
        ----------
        tag : str
            The tag to write.
        value : Any
            The value to write.

        Returns
        -------
        str
        """
        format_index, _ = self._determine_format_option(tag, value)
        return self.format_options[format_index].write(tag, value)

    def get_token_len(self) -> int:
        """Get the token length of the tag.

        Get the token length of the tag.

        Returns
        -------
        int
            The token length of the tag.
        """
        raise NotImplementedError(
            "This method is not supposed to be called\
                                  directonly on MultiformatTag objects!"
        )


@dataclass
class BoolTagContainer(TagContainer):
    """BoolTagContainer class for handling the subtags to the "dump" tag.

    This class is used to handle the subtags to the "dump" tag. All subtags
    are freqs for dump, and all values for these tags are boolean values that
    are read given the existence of their "var" name.
    """

    # Leaving this as a warning, it confuses the hell out of pre-commit
    # even if it is correct.
    # subtags: dict[str, BoolTag] = None

    def read(self, tag: str, value_str: str) -> dict:
        """Read the value string for this tag.

        This method reads the value string for this tag. It is used to parse the
        value string for the tag and return the parsed value.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value_str : str
            The value string to read.

        Returns
        -------
        dict
        """
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
        for subtag, subtag_type in self.subtags.items():
            if not subtag_type.optional and subtag not in subdict:
                raise ValueError(
                    f"The {subtag} tag is not optional but was not populated \
                        during the read!"
                )
        if len(value) > 0:
            raise ValueError(
                f"Something is wrong in the JDFTXInfile formatting, some \
                    values were not processed: {value}"
            )
        return subdict


@dataclass
class DumpTagContainer(TagContainer):
    """DumpTagContainer class for handling the "dump" tag.

    This class is used to handle the "dump" tag.
    """

    # subtags: dict[str, BoolTagContainer] = None

    def read(self, tag: str, value_str: str) -> dict:
        """Read the value string for this tag.

        This method reads the value string for this tag. It is used to parse the
        value string for the tag and return the parsed value.

        Parameters
        ----------
        tag : str
            The tag to read the value string for.
        value_str : str
            The value string to read.

        Returns
        -------
        dict
        """
        value = value_str.split()
        tempdict = {}
        # Each subtag is a freq, which will be a BoolTagContainer
        for subtag, subtag_type in self.subtags.items():
            if subtag in value:
                idx_start = value.index(subtag)
                subtag_value = " ".join(value[(idx_start + 1) :])
                tempdict[subtag] = subtag_type.read(subtag, subtag_value)
                del value[idx_start:]
        # reorder all tags to match order of __MASTER_TAG_LIST__ and do
        # coarse-grained validation of read
        subdict = {x: tempdict[x] for x in self.subtags if x in tempdict}
        for subtag, subtag_type in self.subtags.items():
            if not subtag_type.optional and subtag not in subdict:
                raise ValueError(
                    f"The {subtag} tag is not optional but was not populated \
                        during the read!"
                )
        if len(value) > 0:
            raise ValueError(
                f"Something is wrong in the JDFTXInfile formatting, some \
                    values were not processed: {value}"
            )
        return subdict
