"""Classes for reading/manipulating/writing JDFTx input files.

Classes for reading/manipulating/writing JDFTx input files.

Note: JDFTXInfile will be moved back to its own module once a more broad inputs
class is written.
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np
import scipy.constants as const
from monty.json import MSONable

from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import Element
from pymatgen.core.units import bohr_to_ang
from pymatgen.io.jdftx.generic_tags import AbstractTag, BoolTagContainer, DumpTagContainer, MultiformatTag, TagContainer
from pymatgen.io.jdftx.jdftxinfile_default_inputs import default_inputs
from pymatgen.io.jdftx.jdftxinfile_master_format import (
    __PHONON_TAGS__,
    __TAG_LIST__,
    __WANNIER_TAGS__,
    MASTER_TAG_LIST,
    get_tag_object,
    get_tag_object_on_val,
)
from pymatgen.util.io_utils import clean_lines
from pymatgen.util.typing import SpeciesLike

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.util.typing import PathLike

__author__ = "Jacob Clary, Ben Rich"

# TODO: Add check for whether all ions have or lack velocities.
# TODO: Add default value filling like JDFTx does.
# TODO: Add more robust checking for if two repeatable tag values represent the
# same information. This is likely fixed by implementing filling of default values.
# TODO: Incorporate something to collapse repeated dump tags of the same frequency
# into a single value.


class JDFTXInfile(dict, MSONable):
    """Class for reading/writing JDFtx input files.

    JDFTxInfile object for reading and writing JDFTx input files.
    Essentially a dictionary with some helper functions.
    """

    path_parent: str | None = None  # Only gets a value if JDFTXInfile is initialized with from_file

    def __init__(self, params: dict[str, Any] | None = None) -> None:
        """
        Create a JDFTXInfile object.

        Args:
            params (dict): Input parameters as a dictionary.
        """
        super().__init__()
        if params is not None:
            self.update(params)

    def __str__(self) -> str:
        """Return str representation of JDFTXInfile.

        Returns:
            str: String representation of JDFTXInfile.
        """
        return "".join([line + "\n" for line in self.get_text_list()])

    def __add__(self, other: JDFTXInfile) -> JDFTXInfile:
        """Add existing JDFTXInfile object to method caller JDFTXInfile object.

        Add all the values of another JDFTXInfile object to this object. Facilitate the use of "standard" JDFTXInfiles.
        Repeatable tags are appended together. Non-repeatable tags are replaced with their value from the other object.

        Args:
            other (JDFTXInfile): JDFTXInfile object to add to the method caller object.

        Returns:
            JDFTXInfile: The combined JDFTXInfile object.
        """
        # Deepcopy needed here, or else in `jif1 = jif2 + jif3`, `jif1` will become a reference to `jif2`
        params: dict[str, Any] = deepcopy(dict(self.items()))
        for key, val in other.items():
            if key in self:
                if val is params[key]:
                    # Unlinking the two objects fully cannot be done by deepcopy for some reason
                    continue
                tag_object = get_tag_object(key)
                if tag_object.can_repeat:
                    if isinstance(val, list):
                        for subval in list(val):
                            # Minimum effort to avoid duplicates
                            if subval not in params[key]:
                                params[key].append(subval)
                    else:
                        params[key].append(val)
                else:
                    params[key] = val
            else:
                params[key] = val
        return type(self)(params)

    def __setitem__(self, key: str, value: Any) -> None:
        """Set an item in the JDFTXInfile.

        Set an item in the JDFTXInfile. This magic method is set explicitly to immediately validate when a user sets a
        tag's value, and to perform any conversion necessary.

        Args:
            key (str): Key to set.
            value (Any): Value to set.
        """
        if not isinstance(key, str):
            raise TypeError(f"{key} is not a string!")
        try:
            # tag_object = get_tag_object_on_val(key, value)
            tag_object = get_tag_object(key)
        except KeyError:
            raise KeyError(f"The {key} tag is not in MASTER_TAG_LIST")
        _remaining_values = None
        if tag_object.can_repeat and isinstance(value, list):
            if len(value) > 1:
                _remaining_values = value[1:]
            value = value[0]
        tag_object = get_tag_object_on_val(key, value)

        # if isinstance(tag_object, MultiformatTag):
        #     if isinstance(value, str):
        #         i = tag_object.get_format_index_for_str_value(key, value)
        #     else:
        #         i, _ = tag_object._determine_format_option(key, value)
        #     tag_object = tag_object.format_options[i]
        if tag_object.can_repeat and key in self:
            del self[key]
        if tag_object.can_repeat and not isinstance(value, list):
            value = [value]
        params: dict[str, Any] = {}
        if self._is_numeric(value):
            value = str(value)
        if not tag_object.can_repeat:
            value = [value]  # Shortcut to avoid writing a separate block for non-repeatable tags
        for v in value:
            processed_value = tag_object.read(key, v) if isinstance(v, str) else v
            params = self._store_value(params, tag_object, key, processed_value)
            self.update(params)
            self.validate_tags(try_auto_type_fix=True, error_on_failed_fix=True)
        if _remaining_values is not None:
            for v in _remaining_values:
                self.append_tag(key, v)

    @classmethod
    def from_dict(cls, d: dict[str, Any], validate_value_boundaries=True) -> JDFTXInfile:
        """Create JDFTXInfile from a dictionary.

        Args:
            d (dict): Dictionary to create JDFTXInfile from.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        instance = cls()
        for k, v in d.items():
            if k not in ("@module", "@class"):
                instance[k] = v
        if validate_value_boundaries:
            instance.validate_boundaries()
        return instance

    @classmethod
    def from_file(
        cls,
        filename: PathLike,
        dont_require_structure: bool = False,
        sort_tags: bool = True,
        assign_path_parent: bool = True,
        validate_value_boundaries: bool = True,
    ) -> JDFTXInfile:
        """Read a JDFTXInfile object from a file.

        Args:
            filename (PathLike): Filename to read from.
            dont_require_structure (bool, optional): Whether to require structure tags. Defaults to False.
            sort_tags (bool, optional): Whether to sort the tags. Defaults to True.
            assign_path_parent (bool, optional): Whether to assign the parent directory of the input file for include
                tags. Defaults to True.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        path_parent = None
        if assign_path_parent:
            path_parent = Path(filename).parents[0]
        with open(filename) as file:
            return cls.from_str(
                file.read(),
                dont_require_structure=dont_require_structure,
                sort_tags=sort_tags,
                path_parent=path_parent,
                validate_value_boundaries=validate_value_boundaries,
            )

    @classmethod
    def from_structure(
        cls,
        structure: Structure,
        selective_dynamics: ArrayLike | None = None,
        write_cart_coords: bool = False,
    ) -> JDFTXInfile:
        """Create a JDFTXInfile object from a pymatgen Structure.

        Args:
            structure (Structure): Structure to convert.
            selective_dynamics (ArrayLike, optional): Selective dynamics attribute for each site if available.
                Shape Nx1, by default None.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        jdftxstructure = JDFTXStructure(structure, selective_dynamics)
        jdftxstructure.write_cart_coords = write_cart_coords
        return cls.from_jdftxstructure(jdftxstructure)

    @classmethod
    def from_jdftxstructure(
        cls,
        jdftxstructure: JDFTXStructure,
    ) -> JDFTXInfile:
        """Create a JDFTXInfile object from a JDFTXStructure object.

        Args:
            jdftxstructure (JDFTXStructure): JDFTXStructure object to convert.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        jstr = jdftxstructure.get_str()
        return cls.from_str(jstr)

    @classmethod
    def from_str(
        cls,
        string: str,
        dont_require_structure: bool = False,
        sort_tags: bool = True,
        path_parent: Path | None = None,
        validate_value_boundaries: bool = True,
    ) -> JDFTXInfile:
        """Read a JDFTXInfile object from a string.

        Args:
            string (str): String to read from.
            dont_require_structure (bool, optional): Whether to require structure tags. Defaults to False.
            sort_tags (bool, optional): Whether to sort the tags. Defaults to True.
            path_parent (Path, optional): Path to the parent directory of the input file for include tags.
                                          Defaults to None.
            validate_value_boundaries (bool, optional): Whether to validate the value boundaries. Defaults to True.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        lines: list[str] = list(clean_lines(string.splitlines()))
        lines = cls._gather_tags(lines)

        params: dict[str, Any] = {}
        # process all tag value lines using specified tag formats in MASTER_TAG_LIST
        for line in lines:
            tag_object, tag, value = cls._preprocess_line(line)
            processed_value = tag_object.read(tag, value)
            params = cls._store_value(params, tag_object, tag, processed_value)  # this will change with tag categories

        if "include" in params:
            for filename in params["include"]:
                _filename = filename
                if not Path(_filename).exists():
                    if path_parent is not None:
                        _filename = path_parent / filename
                    if not Path(_filename).exists():
                        raise ValueError(f"The include file {filename} ({_filename}) does not exist!")
                params.update(cls.from_file(_filename, dont_require_structure=True, assign_path_parent=False))
            del params["include"]

        if (
            not dont_require_structure
            and "lattice" not in params
            and "ion" not in params
            and "ion-species" not in params
        ):
            raise ValueError("This input file is missing required structure tags")
        if sort_tags:
            params = {tag: params[tag] for tag in __TAG_LIST__ if tag in params}
        instance = cls(params)
        if validate_value_boundaries:
            instance.validate_boundaries()
        return instance

    @classmethod
    def get_list_representation(cls, jdftxinfile: JDFTXInfile) -> JDFTXInfile:
        """Convert JDFTXInfile object properties into list representation.

        Args:
            jdftxinfile (JDFTXInfile): JDFTXInfile object to convert.
        """
        reformatted_params = deepcopy(jdftxinfile.as_dict(skip_module_keys=True))
        # rest of code assumes lists are lists and not np.arrays
        reformatted_params = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in reformatted_params.items()}
        for tag, value in reformatted_params.items():
            tag_object = get_tag_object(tag)
            if all(
                [
                    tag_object.allow_list_representation,
                    tag_object.is_tag_container,
                    cls._needs_conversion("dict-to-list", value),
                ]
            ):
                reformatted_params.update({tag: tag_object.get_list_representation(tag, value)})
        return cls(reformatted_params)

    @classmethod
    def get_dict_representation(cls, jdftxinfile: JDFTXInfile) -> JDFTXInfile:
        """Convert JDFTXInfile object properties into dict representation.

        Args:
            jdftxinfile (JDFTXInfile): JDFTXInfile object to convert.
        """
        reformatted_params = deepcopy(jdftxinfile.as_dict(skip_module_keys=True))
        # Just to make sure only passing lists and no more numpy arrays
        reformatted_params = {
            k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in reformatted_params.items()
        }  # rest of code assumes lists are lists and not np.arrays
        for tag, value in reformatted_params.items():
            tag_object = get_tag_object(tag)
            if all(
                [
                    tag_object.allow_list_representation,
                    tag_object.is_tag_container,
                    cls._needs_conversion("list-to-dict", value),
                ]
            ):
                reformatted_params.update({tag: tag_object.get_dict_representation(tag, value)})
        return cls(reformatted_params)

    def as_dict(self, sort_tags: bool = True, skip_module_keys: bool = False) -> dict:
        """Return JDFTXInfile as MSONable dict.

        Args:
            sort_tags (bool, optional): Whether to sort the tags. Defaults to True.
            skip_module_keys (bool, optional): Whether to skip the module keys. Defaults to False.

        Returns:
            dict: JDFTXInfile as MSONable dict.
        """
        params = dict(self)
        if sort_tags:
            params = {tag: params[tag] for tag in __TAG_LIST__ if tag in params}
        if not skip_module_keys:
            params["@module"] = type(self).__module__
            params["@class"] = type(self).__name__
        return params

    def copy(self) -> JDFTXInfile:
        """Return a copy of the JDFTXInfile object.

        Returns:
            JDFTXInfile: Copy of the JDFTXInfile object.
        """
        return type(self)(self)

    def get_text_list(self) -> list[str]:
        """Get a list of strings representation of the JDFTXInfile.

        Returns:
            list[str]: List of strings representation of the JDFTXInfile.
        """
        self_as_dict = self.get_dict_representation(self)

        text: list[str] = []
        for tag_group in MASTER_TAG_LIST:
            added_tag_in_group = False
            for tag in MASTER_TAG_LIST[tag_group]:
                if tag not in self:
                    continue
                added_tag_in_group = True
                tag_object: AbstractTag = MASTER_TAG_LIST[tag_group][tag]
                if isinstance(tag_object, MultiformatTag):
                    i, _ = tag_object._determine_format_option(tag, self_as_dict[tag])
                    tag_object = tag_object.format_options[i]
                if tag_object.can_repeat and isinstance(self_as_dict[tag], list):
                    text += [tag_object.write(tag, entry) for entry in self_as_dict[tag]]
                else:
                    text.append(tag_object.write(tag, self_as_dict[tag]))
            if added_tag_in_group:
                text.append("")
        return text

    def write_file(self, filename: PathLike) -> None:
        """Write JDFTXInfile to an in file.

        Args:
            filename (PathLike): Filename to write to.
        """
        with open(filename, mode="w") as file:
            file.write(str(self))

    @classmethod
    def to_jdftxstructure(cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False) -> JDFTXStructure:
        """Convert JDFTXInfile to JDFTXStructure object.

        Converts JDFTx lattice, lattice-scale, ion tags into JDFTXStructure, with Pymatgen structure as attribute.

        Args:
            jdftxinfile (JDFTXInfile): JDFTXInfile object to convert.
            sort_structure (bool, optional): Whether to sort the structure. Useful if species are not grouped properly
                together. Defaults to False.
        """
        # use dict representation so it's easy to get the right column for moveScale,
        # rather than checking for velocities
        jdftxinfile_dict = cls.get_dict_representation(jdftxinfile)
        return JDFTXStructure.from_jdftxinfile(jdftxinfile_dict, sort_structure=sort_structure)

    @classmethod
    def to_pmg_structure(cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False) -> Structure | None:
        """Convert JDFTXInfile to pymatgen Structure object.

        Converts JDFTx lattice, lattice-scale, ion tags into pymatgen Structure.

        Args:
            jdftxinfile (JDFTXInfile): JDFTXInfile object to convert.
            sort_structure (bool, optional): Whether to sort the structure. Useful if species are not grouped properly
                together. Defaults to False.

        Returns:
            Structure: The created pymatgen Structure object.
        """
        # use dict representation so it's easy to get the right column for
        # moveScale, rather than checking for velocities
        print(jdftxinfile.get_dict_representation(jdftxinfile), "INPUT DICT REP")
        jdftxstructure = JDFTXStructure.from_jdftxinfile(
            jdftxinfile.get_dict_representation(jdftxinfile),
            sort_structure=sort_structure,
        )
        return jdftxstructure.structure

    def strip_structure_tags(self) -> None:
        """Strip all structural tags from the JDFTXInfile.

        Strip all structural tags from the JDFTXInfile. This is useful for preparing a JDFTXInfile object
        from a pre-existing calculation for a new structure.
        """
        strucural_tags = ["lattice", "ion", "lattice-scale", "coords-type"]
        for tag in strucural_tags:
            if tag in self:
                del self[tag]

    def append_tag(self, tag: str, value: Any) -> None:
        """Append a value to a tag.

        Append a value to a tag. Use this method instead of directly appending the list contained in the tag, such that
        the value is properly processed.

        Args:
            tag (str): Tag to append to.
            value (Any): Value to append.
        """
        tag_object = get_tag_object(tag)
        if isinstance(tag_object, MultiformatTag):
            if isinstance(value, str):
                i = tag_object.get_format_index_for_str_value(tag, value)
            else:
                i, _ = tag_object._determine_format_option(tag, value)
            tag_object = tag_object.format_options[i]
        if not tag_object.can_repeat:
            raise ValueError(f"The tag '{tag}' cannot be repeated and thus cannot be appended")
        params: dict[str, Any] = self.as_dict(skip_module_keys=True)
        processed_value = tag_object.read(tag, value) if isinstance(value, str) else value
        params = self._store_value(params, tag_object, tag, processed_value)
        self.update(params)

    # This method is called by setitem, but setitem is circumvented by update,
    # so this method's parameters is still necessary
    def validate_tags(
        self,
        try_auto_type_fix: bool = False,
        error_on_failed_fix: bool = True,
        return_list_rep: bool = False,
    ) -> None:
        """Validate the tags in the JDFTXInfile.

        Validate the tags in the JDFTXInfile. If try_auto_type_fix is True, will attempt to fix the tags. If
        error_on_failed_fix is True, will raise an error if the tags cannot be fixed. If return_list_rep is True, will
        return the tags in list representation.

        Args:
            try_auto_type_fix (bool, optional): Whether to attempt to fix the tags. Defaults to False.
            error_on_failed_fix (bool, optional): Whether to raise an error if the tags cannot be fixed.
                                                  Defaults to True.
            return_list_rep (bool, optional): Whether to return the tags in list representation. Defaults to False.
        """
        for tag in self:
            tag_object = get_tag_object(tag)
            checked_tag, is_tag_valid, value = tag_object.validate_value_type(
                tag, self[tag], try_auto_type_fix=try_auto_type_fix
            )
            should_warn = not is_tag_valid
            if return_list_rep and tag_object.allow_list_representation:
                value = tag_object.get_list_representation(tag, value)
            if error_on_failed_fix and should_warn and try_auto_type_fix:
                raise ValueError(f"The {tag} tag with value:\n{self[tag]}\ncould not be fixed!")
            if try_auto_type_fix and is_tag_valid:
                self.update({tag: value})
            if should_warn:
                warnmsg = f"The {tag} tag with value:\n{self[tag]}\nhas incorrect typing!"
                if any(isinstance(tag_object, tc) for tc in [TagContainer, DumpTagContainer, BoolTagContainer]):
                    warnmsg += "(Check earlier warnings for more details)\n"
                warnings.warn(warnmsg, stacklevel=2)

    def validate_boundaries(self) -> None:
        """Validate the boundaries of the JDFTXInfile.

        Validate the boundaries of the JDFTXInfile. This is a placeholder for future functionality.
        """
        error_strs: list[str] = []
        for tag in self:
            tag_object = get_tag_object(tag)
            is_valid, error_str = tag_object.validate_value_bounds(tag, self[tag])
            if not is_valid:
                error_strs.append(error_str)
        if len(error_strs) > 0:
            err_cat = "\n".join(error_strs)
            raise ValueError(
                f"The following boundary errors were found in the JDFTXInfile:\n{err_cat}\n"
                "\n Hint - if you are reading from a JDFTX out file, you need to set validate_value_boundaries "
                "to False, as JDFTx will dump values at non-inclusive boundaries (ie 0.0 for values strictly > 0.0)."
            )

    # TODO: Add examples in the docustring on how to use this function
    def is_comparable_to(
        self,
        other: JDFTXInfile,
        exclude_tags: list[str] | None = None,
        exclude_tag_categories: list[str] | None = None,
        ensure_include_tags: list[str] | None = None,
    ) -> bool:
        """Check if two JDFTXInfile objects are comparable.

        Args:
            other (JDFTXInfile): Other JDFTXInfile object to compare to.
            exclude_tags (list[str], optional): Tags to exclude from comparison. Defaults to None.
                Set to contain tags you are already accounting for being different to exclude them.
            exclude_tag_categories (list[str], optional): Tag categories to exclude from comparison. Defaults to None.
                If left as None, will exclude all tags in the "export", "restart", and "structure" categories.
                Add all tags of specified tag categories to be excluded in the comparison.
            ensure_include_tags (list[str], optional): Tags to ensure are included in the comparison. Defaults to None.
                If left as None, will make sure 'ion-species' is included in comparison.
                Helpful for tags in categories you want to exclude but want to ensure are included in the comparison.

        Returns:
            bool: Whether the two JDFTXInfile objects are comparable.
        """
        if exclude_tag_categories is None:
            exclude_tag_categories = ["export", "restart", "structure"]
        if ensure_include_tags is None:
            ensure_include_tags = []
        differing_tags = self.get_filtered_differing_tags(
            other,
            exclude_tags=exclude_tags,
            exclude_tag_categories=exclude_tag_categories,
            ensure_include_tags=ensure_include_tags,
        )
        return len(differing_tags) == 0

    def get_filtered_differing_tags(
        self,
        other: JDFTXInfile,
        exclude_tags: list[str] | None = None,
        exclude_tag_categories: list[str] | None = None,
        ensure_include_tags: list[str] | None = None,
    ) -> list[str]:
        """Return which tags differ between self and other JDFTXInfile.

        Returns a list of tags that differ between the two JDFTXInfile objects, excluding any tags specified in
        exclude_tags or exclude_tag_categories. This is useful for comparing two JDFTXInfile objects while ignoring
        certain tags that are expected to differ.

        Args:
            other (JDFTXInfile): Other JDFTXInfile object to compare to.
            exclude_tags (list[str], optional): Tags to exclude from comparison. Defaults to None.
                Set to contain tags you are already accounting for being different to exclude them.
            exclude_tag_categories (list[str], optional): Tag categories to exclude from comparison. Defaults to None.
                Add all tags of specified tag categories to be excluded in the comparison.
            ensure_include_tags (list[str], optional): Tags to ensure are included in the comparison. Defaults to None.
                Helpful for tags in categories you want to exclude but want to ensure are included in the comparison.
        Returns:
            list[str]: Tags that differ between the two JDFTXInfile objects.
        """
        if exclude_tags is None:
            exclude_tags = []
        if exclude_tag_categories is None:
            exclude_tag_categories = []
        if ensure_include_tags is None:
            ensure_include_tags = []
        for category in exclude_tag_categories:
            exclude_tags += list(MASTER_TAG_LIST[category].keys())
        for tag in ensure_include_tags:
            if tag in exclude_tags:
                exclude_tags.remove(tag)
        differing_tags = self.get_differing_tags(other)
        for exclude_tag in exclude_tags:
            if exclude_tag in differing_tags:
                differing_tags.remove(exclude_tag)
        return differing_tags

    # get_differing_tags for a full list of tags that differ between two JDFTXInfile objects
    def get_differing_tags(self, other: JDFTXInfile) -> list[str]:
        """Return which tags differ between self and other JDFTXInfile.

        Args:
            other (JDFTXInfile): Other JDFTXInfile object to compare to.

        Returns:
            list[str]: Tags that differ between the two JDFTXInfile objects.
        """
        list1 = self.get_differing_tags_from(other)
        list2 = other.get_differing_tags_from(self)
        return list(set(list1 + list2))

    # get_differing_tags_from for only tags contained in self that differ from other JDFTXInfile
    def get_differing_tags_from(self, other: JDFTXInfile) -> list[str]:
        """Return which self-contained tags differ from other JDFTXInfile.

        Args:
            other (JDFTXInfile): Other JDFTXInfile object to compare to.

        Returns:
            list[str]: Tags from this JDFTXInfile that differ from the other JDFTXInfile.
        """
        # TODO: `TagContainer.is_equal_to(val1, tc2, val2)`, where val2 is the fetched default
        # value for the tag, will return a false negative if val1 does not specify all the
        # possible values for a tag. To fix this, missing subtags in both vals for tagcontainer
        # comparisons should be populated with the default values for the tag.
        # (ie if val1 = {a: 1}, val2 = {b: 2}, and the default value for the tag is {a: 1, b: 2},
        # then val1 -> {a: 1, b: 2}, val2 -> {a: 1, b: 2} and the comparison will return True, but
        # if val2 = {b: 3}, then val2 -> {a: 1, b: 3} and the comparison will return False)
        differing_tags = []
        for tag in self:
            val1 = self[tag]
            obj1 = get_tag_object_on_val(tag, val1)
            val2 = None
            obj2 = None
            if tag not in other:
                if tag in ref_infile:
                    val2 = ref_infile[tag]
                    obj2 = get_tag_object_on_val(tag, val2)
                else:
                    differing_tags.append(tag)
            else:
                val2 = other[tag]
                obj2 = get_tag_object_on_val(tag, val2)
            if not (val2 is None or obj2 is None) and not obj1.is_equal_to(val1, obj2, val2):
                differing_tags.append(tag)
        return differing_tags

    @property
    def structure(self) -> Structure | None:
        """Return a pymatgen Structure object.

        Returns:
            Structure: Pymatgen structure object.
        """
        return self.to_pmg_structure(self)

    @classmethod
    def _from_dict(cls, dct: dict[str, Any]) -> JDFTXInfile:
        """Parse a dictionary to create a JDFTXInfile object.

        Args:
            dct (dict): Dictionary to parse.

        Returns:
            JDFTXInfile: The created JDFTXInfile object.
        """
        temp = cls({k: v for k, v in dct.items() if k not in ("@module", "@class")})
        temp = cls.get_dict_representation(cls.get_list_representation(temp))
        return cls.get_list_representation(temp)

    @staticmethod
    def _preprocess_line(line: str) -> tuple[AbstractTag, str, str]:
        """Preprocess a line from a JDFTXInfile, splitting it into a tag object, tag name, and value.

        Args:
            line (str): Line from the input file.

        Returns:
            tuple[AbstractTag, str, str]: Tag object, tag name, and value.
        """
        line_list = line.strip().split(maxsplit=1)
        tag: str = line_list[0].strip()
        if tag in __PHONON_TAGS__:
            raise ValueError("Phonon functionality has not been added!")
        if tag in __WANNIER_TAGS__:
            raise ValueError("Wannier functionality has not been added!")
        if tag not in __TAG_LIST__:
            err_str = f"The {tag} tag in {line_list} is not in MASTER_TAG_LIST and is not a comment, "
            err_str += "something is wrong with this input data!"
            raise ValueError(err_str)
        tag_object = get_tag_object(tag)
        value: str = ""
        if len(line_list) == 2:
            value = line_list[1].strip()
        elif len(line_list) == 1:
            value = ""  # exception for tags where only tagname is used, e.g. dump-only tag
        if isinstance(tag_object, MultiformatTag):
            i = tag_object.get_format_index_for_str_value(tag, value)
            tag_object = tag_object.format_options[i]
        return tag_object, tag, value

    @staticmethod
    def _store_value(
        params: dict[str, list | list[dict[str, dict]] | Any],
        tag_object: AbstractTag,
        tag: str,
        value: Any,
    ) -> dict:
        """Store the value in the params dictionary.

        Args:
            params (dict): Dictionary to store the value in.
            tag_object (AbstractTag): Tag object for holding tag/value pair. This tag object should be the tag object
                which the "tag" string maps to in MASTER_TAG_LIST.
            tag (str): Tag name.
            value (Any): Value to store.

        Returns:
            dict: Updated dictionary with the stored value.
        """
        if tag_object.can_repeat:  # store tags that can repeat in a list
            if tag not in params:
                params[tag] = []
            params[tag].append(value)
        else:
            if tag in params:
                raise ValueError(f"The '{tag}' tag appears multiple times in this input when it should not!")
            params[tag] = value
        return params

    @staticmethod
    def _gather_tags(lines: list[str]) -> list[str]:
        """Gather broken lines into single string for processing later.

        Args:
            lines (list[str]): List of lines from the input file.

        Returns:
            list[str]: List of strings with tags broken across lines combined into single string.
        """
        total_tag = ""
        gathered_strings = []
        for line in lines:
            if line[-1] == "\\":  # "\" indicates line continuation (the first "\" needed to be escaped)
                total_tag += line[:-1].strip() + " "  # remove \ and any extra whitespace
            elif total_tag:  # then finished with line continuations
                total_tag += line
                gathered_strings.append(total_tag)
                total_tag = ""
            else:  # then append line like normal
                gathered_strings.append(line)
        return gathered_strings

    @staticmethod
    def _needs_conversion(conversion: str, value: dict | list[dict] | list | list[list]) -> bool:
        """Determine if a value needs to be converted.

        This method is only ever called by cls.get_list/dict_representation.

        Args:
            conversion (str): Conversion type. ('dict-to-list' (value : dict | list[dict]) or
                                'list-to-dict' (value : list | list[list]))
            value (dict | list[dict] | list | list[list]): Value to check.

        Returns:
            bool: Whether the value needs to be converted.
        """
        # Check if value is not iterable
        try:
            iter(value)
        except TypeError:
            # This is triggered when JDFTXInfile is attempting to convert a non-tagcontainer to list/dict representation
            # The return boolean is meaningless in this case, so just returning False to avoid the value hitting the
            # "for x in value" loop below.
            return False
        if conversion == "list-to-dict":
            flag = False
        elif conversion == "dict-to-list":
            flag = True
        else:
            raise ValueError(f"Conversion type {conversion} is not 'list-to-dict' or 'dict-to-list'")
        if isinstance(value, dict) or all(isinstance(x, dict) for x in value):
            return flag
        return not flag

    def _is_numeric(self, value: Any) -> bool:
        """Check if a value is numeric.

        Args:
            value (Any): Value to check.

        Returns:
            bool: Whether the value is numeric.
        """
        # data-types that might accidentally be identified as numeric
        if type(value) in [bool]:
            return False
        try:
            float(value)
            is_numeric = True
        except (ValueError, TypeError):
            is_numeric = False
        return is_numeric


ref_infile = JDFTXInfile.from_dict(default_inputs, validate_value_boundaries=False)


def _allnone(val) -> bool:
    """Check if all elements in a list are None.

    Args:
        val (list | None): List to check.

    Returns:
        bool: Whether all elements in the list are None.
    """
    if val is None:
        return True
    if isinstance(val, list):
        return all(_allnone(x) for x in val)
    return False


def selective_dynamics_site_prop_to_jdftx_interpretable(selective_dynamics: Sequence[Any]) -> ArrayLike:
    """Convert selective dynamics site property to JDFTX interpretable format.

    Args:
        selective_dynamics (list[list[bool]]): Selective dynamics site property.

    Returns:
        ArrayLike: JDFTX interpretable format of selective dynamics.
    """
    _selective_dynamics = []
    for i in range(len(selective_dynamics)):
        # Not worrying about partial free axes for now, just setting movescale to 0 if any axis is fixed
        if any(not x for x in selective_dynamics[i]):
            _selective_dynamics.append(0)
        else:
            _selective_dynamics.append(1)
    return np.array(_selective_dynamics)


@dataclass
class JDFTXStructure(MSONable):
    """Object for representing the data in JDFTXStructure tags.

    Attributes:
        structure (Structure): Associated Structure.
        selective_dynamics (ArrayLike): Selective dynamics attribute for each site if available. Shape Nx1.
        sort_structure (bool): Whether to sort the structure. Useful if species are not grouped properly together.
            Defaults to False.
        write_cart_coords (bool): Whether to write cartesian coordinates. Defaults to False.
        velocities (list[np.ndarray | list[float] | None]): Velocities for each site if available. Kept separate from
            structure.site_properties["velocities"] to allow inhomogeneous shaping.
        constraint_types (list[str | None]): Constraint types for each site if available. Kept separate from
            structure.site_properties["constraint_types"] to allow inhomogeneous shaping.
        constraint_vectors (list[np.ndarray | list[np.ndarray] | None]): Constraint vectors for each site if
            available. Kept separate from structure.site_properties["constraint_vectors"] to allow inhomogeneous
            shaping. i'th entry must be a list of vectors if the i'th constraint_type is "HyperPlane". Otherwise,
            a single vector is expected.
        hyperplane_group_names (list[list[str] | None]): Group names for each site if available. Kept separate from
            structure.site_properties["group_names"] to allow inhomogeneous shaping. i'th entry must be a list of group
            names if the i'th constraint_type is "HyperPlane". Otherwise, None is expected.
    """

    structure: Structure | None = None
    selective_dynamics: ArrayLike | None = None
    sort_structure: bool = False
    write_cart_coords: bool = False
    velocities: list | None = None
    constraint_types: list | None = None
    constraint_vectors: list | None = None
    hyperplane_group_names: list | None = None

    def __post_init__(self) -> None:
        """Post init function for JDFTXStructure.

        Asserts self.structure is ordered, and adds selective dynamics if needed.
        """
        if ((self.structure is not None) and ("selective_dynamics" in self.structure.site_properties)) and (
            self.selective_dynamics is None
        ):
            self.selective_dynamics = selective_dynamics_site_prop_to_jdftx_interpretable(
                self.structure.site_properties["selective_dynamics"]
            )
        for i in range(len(self.structure.species)):
            name = ""
            if isinstance(self.structure.species[i], str):
                _name_str = self.structure.species[i]
            elif isinstance(self.structure.species[i], SpeciesLike):
                _name_str = self.structure.species[i].symbol
                if not isinstance(_name_str, str):
                    _name_str = _name_str.symbol
            else:
                raise TypeError("Species must be a string or SpeciesLike object")
            name_str = str(_name_str)
            for j in range(len(name_str)):
                if not name_str[j].isdigit():
                    name += name_str[j]
            self.structure.species[i] = name
        if self.structure.is_ordered:
            site_properties = {}
            if self.selective_dynamics is not None:
                selective_dynamics = np.array(self.selective_dynamics)
                if not selective_dynamics.all():
                    site_properties["selective_dynamics"] = selective_dynamics

            # create new copy of structure so can add selective dynamics and sort atoms if needed
            structure = Structure.from_sites(self.structure.sites)
            self.structure = structure.copy(site_properties=site_properties)
            if self.sort_structure:
                self.structure = self.structure.get_sorted_structure()
        else:
            raise ValueError("Disordered structure with partial occupancies cannot be converted into JDFTXStructure!")
        if _allnone(self.selective_dynamics):
            self.selective_dynamics = None
        if _allnone(self.velocities):
            self.velocities = None
        if _allnone(self.constraint_types):
            self.constraint_types = None
        if _allnone(self.constraint_vectors):
            self.constraint_vectors = None
        if _allnone(self.hyperplane_group_names):
            self.hyperplane_group_names = None

    def __repr__(self) -> str:
        """Return representation of JDFTXStructure file.

        Returns:
            str: Representation of JDFTXStructure file.
        """
        return f"JDFTXStructure({self.get_str()})"

    def __str__(self) -> str:
        """Return string representation of JDFTXStructure file.

        Returns:
            str: String representation of JDFTXStructure file.
        """
        return self.get_str()

    @property
    def natoms(self) -> int:
        """Return number of atoms.

        Returns:
            int: Number of sites.
        """
        return len(self.structure.species)

    @classmethod
    def from_str(cls, data: str) -> JDFTXStructure:
        """Read JDFTXStructure from string.

        Args:
            data (str): String to read from.

        Returns:
            JDFTXStructure: The created JDFTXStructure object.
        """
        return cls.from_jdftxinfile(JDFTXInfile.from_str(data))

    @classmethod
    def from_file(cls, filename: str) -> JDFTXStructure:
        """Read JDFTXStructure from file.

        Args:
            filename (str): Filename to read from.

        Returns:
            JDFTXStructure: The created JDFTXStructure object.
        """
        return cls.from_jdftxinfile(JDFTXInfile.from_file(filename))

    @classmethod
    def from_jdftxinfile(cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False) -> JDFTXStructure:
        """Get JDFTXStructure from JDFTXInfile.

        Args:
            jdftxinfile (JDFTXInfile): JDFTXInfile object.
            sort_structure (bool, optional): Whether to sort the structure. Useful if species are not grouped properly
                together as JDFTx output will have species sorted.

        Returns:
            JDFTXStructure: The created JDFTXStructure object.
        """

        lattice = _infile_to_pmg_lattice(jdftxinfile)
        atomic_symbols = [x["species-id"] for x in jdftxinfile["ion"]]
        coords = np.array([[x["x0"], x["x1"], x["x2"]] for x in jdftxinfile["ion"]])
        coords *= const.value("Bohr radius") * 10**10  # Bohr radius in Ang; convert to Ang
        selective_dynamics = np.array([x["moveScale"] for x in jdftxinfile["ion"]])
        # Gather optional tags
        velocities = [
            np.array([float(v) for v in [x["v"]["vx0"], x["v"]["vx1"], x["v"]["vx2"]]]) if "v" in x else None
            for x in jdftxinfile["ion"]
        ]
        constraint_types = [x.get("constraint-type", None) for x in jdftxinfile["ion"]]
        constraint_vectors = [
            [x["d0"], x["d1"], x["d2"]] if "constraint-type" in x else None for x in jdftxinfile["ion"]
        ]
        hyperplane_vectors = [
            [[[y["d0"], y["d1"], y["d2"]] for y in x["HyperPlane"]]] if "HyperPlane" in x else None
            for x in jdftxinfile["ion"]
        ]
        hyperplane_groups = [None for _ in range(len(hyperplane_vectors))]
        if not all(x is None for x in hyperplane_vectors):
            for i in range(len(hyperplane_vectors)):
                if hyperplane_vectors[i] is not None:
                    constraint_types[i] = "HyperPlane"
                    constraint_vectors[i] = hyperplane_vectors[i]
            hyperplane_groups = [
                [[y["group"] for y in x["HyperPlane"]] if "HyperPlane" in x else None] for x in jdftxinfile["ion"]
            ]

        coords_are_cartesian = False  # is default for JDFTx
        if "coords-type" in jdftxinfile:
            coords_are_cartesian = jdftxinfile["coords-type"] == "Cartesian"

        struct: Structure = Structure(
            lattice,
            atomic_symbols,
            coords,
            to_unit_cell=False,
            validate_proximity=False,
            coords_are_cartesian=coords_are_cartesian,
        )
        return cls(
            structure=struct,
            selective_dynamics=selective_dynamics,
            sort_structure=sort_structure,
            velocities=velocities,
            constraint_types=constraint_types,
            constraint_vectors=constraint_vectors,
            hyperplane_group_names=hyperplane_groups,
        )

    def get_str(self, in_cart_coords: bool | None = None) -> str:
        """Return a string to be written as JDFTXInfile tags.

        Allows extra options as compared to calling str(JDFTXStructure) directly.

        Args:
            in_cart_coords (bool, optional): Whether coordinates are output in direct or Cartesian.

        Returns:
            str: Representation of JDFTXInfile structure tags.
        """
        jdftx_tag_dict: dict[str, Any] = {}

        if in_cart_coords is None:
            in_cart_coords = self.write_cart_coords

        lattice = np.copy(self.structure.lattice.matrix)
        lattice = lattice.T  # transpose to get into column-vector format
        lattice /= const.value("Bohr radius") * 10**10  # Bohr radius in Ang; convert to Bohr

        jdftx_tag_dict["lattice"] = {f"R{i}{j}": lattice[i][j] for i in range(3) for j in range(3)}
        jdftx_tag_dict["ion"] = []
        jdftx_tag_dict["coords-type"] = "Cartesian" if in_cart_coords else "Lattice"
        valid_labels = [
            value.symbol for key, value in Element.__dict__.items() if not key.startswith("_") and not callable(value)
        ]
        if (self.constraint_types is not None) and (self.constraint_vectors is None):
            raise ValueError("Constraint types are present but constraint vectors are not!")
        if (self.constraint_vectors is not None) and (self.constraint_types is None):
            raise ValueError("Constraint vectors are present but constraint types are not!")
        if self.structure is None:
            return ""
        for i, site in enumerate(self.structure):
            label = site.label
            coords = site.coords * (1 / bohr_to_ang) if in_cart_coords else site.frac_coords
            sd = int(bool(cast("np.ndarray", self.selective_dynamics)[i])) if self.selective_dynamics is not None else 1
            # TODO: This is needlessly complicated, simplify this
            if label not in valid_labels:
                for varname in ["species_string", "specie.name"]:
                    if _multi_hasattr(site, varname) and _multi_getattr(site, varname) in valid_labels:
                        label = _multi_getattr(site, varname)
                        break
                if label not in valid_labels:
                    raise ValueError(f"Could not correct site label {label} for site (index {i})")
            ion_dict_rep = {
                "species-id": label,
                "x0": coords[0],
                "x1": coords[1],
                "x2": coords[2],
                "moveScale": sd,
            }
            # Gather velocities if present
            if (self.velocities is not None) and (self.velocities[i] is not None):
                ion_dict_rep["v"] = {
                    "vx0": self.velocities[i][0] * (1 / bohr_to_ang),
                    "vx1": self.velocities[i][1] * (1 / bohr_to_ang),
                    "vx2": self.velocities[i][2] * (1 / bohr_to_ang),
                }
            # Gather constraints if present
            if ((self.constraint_types is not None) and (self.constraint_vectors is not None)) and (
                self.constraint_types[i] is not None
            ):
                ctype = self.constraint_types[i]
                if ctype == "HyperPlane":
                    ion_dict_rep["HyperPlane"] = []
                    if self.constraint_vectors[i] is not None:
                        for j, hyperplane in enumerate(self.constraint_vectors[i]):
                            _hyperplane = {}
                            for k in range(3):
                                _hyperplane["d" + str(k)] = hyperplane[k]
                            _hyperplane["group"] = self.hyperplane_groups[i][j]
                            ion_dict_rep["HyperPlane"].append(_hyperplane)
                else:
                    ion_dict_rep["constraint-type"] = ctype
                    for j in range(3):
                        ion_dict_rep["d" + str(j)] = self.constraint_vectors[i][j]
            jdftx_tag_dict["ion"].append(ion_dict_rep)
        return str(JDFTXInfile.from_dict(jdftx_tag_dict))

    def write_file(self, filename: PathLike, **kwargs) -> None:
        """Write JDFTXStructure to a file.

        The supported kwargs are the same as those for the JDFTXStructure.get_str method and are passed through
        directly.

        Args:
            filename (PathLike): Filename to write to.
            **kwargs: Kwargs to pass to JDFTXStructure.get_str.
        """
        with open(filename, mode="w") as file:
            file.write(self.get_str(**kwargs))

    def as_dict(self) -> dict:
        """MSONable dict.

        Returns:
            dict: MSONable dictionary representation of the JDFTXStructure.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "selective_dynamics": np.array(self.selective_dynamics).tolist(),
            "velocities": self.velocities,
            "constraint_types": self.constraint_types,
            "constraint_vectors": self.constraint_vectors,
            "hyperplane_group_names": self.hyperplane_group_names,
        }

    @classmethod
    def from_dict(cls, params: dict) -> Self:
        """Get JDFTXStructure from dict.

        Args:
            params (dict): Serialized JDFTXStructure.

        Returns:
            JDFTXStructure: The created JDFTXStructure object.
        """
        return cls(
            Structure.from_dict(params["structure"]),
            selective_dynamics=params["selective_dynamics"],
        )


def _multi_hasattr(varbase: Any, varname: str):
    """Check if object has an attribute (capable of nesting with . splits).

    Check if object has an attribute (capable of nesting with . splits).

    Parameters
    ----------
    varbase
        Object to check.
    varname
        Attribute to check for.

    Returns
    -------
    bool
        Whether the object has the attribute.
    """
    varlist = varname.split(".")
    for i, var in enumerate(varlist):
        if i == len(varlist) - 1:
            return hasattr(varbase, var)
        if hasattr(varbase, var):
            varbase = getattr(varbase, var)
        else:
            return False
    return None


def _multi_getattr(varbase: Any, varname: str):
    """Check if object has an attribute (capable of nesting with . splits).

    Check if object has an attribute (capable of nesting with . splits).

    Parameters
    ----------
    varbase
        Object to check.
    varname
        Attribute to check for.

    Returns
    -------
    Any
        Attribute of the object.
    """
    if not _multi_hasattr(varbase, varname):
        raise AttributeError(f"{varbase} does not have attribute {varname}")
    varlist = varname.split(".")
    for var in varlist:
        varbase = getattr(varbase, var)
    return varbase


def _infile_to_pmg_lattice(jdftxinfile: JDFTXInfile) -> Lattice:
    jl = jdftxinfile["lattice"]
    if "R00" not in jl:
        return _infile_special_format_to_pmg_lattice(jdftxinfile)
    return _infile_matrix_format_to_pmg_lattice(jdftxinfile)


def _infile_special_format_to_pmg_lattice(jdftxinfile: JDFTXInfile) -> Lattice:
    jl = jdftxinfile["lattice"]
    if "modification" in jl:
        raise NotImplementedError("Special case lattices with modification not implemented yet")
    # Second boolean to check if latt-scale was passed but does nothing
    if ("latt-scale" in jl) and (not all(np.isclose(float(x), 1) in jl for x in ["s0", "s1", "s2"])):
        raise NotImplementedError("latt-scale for special case lattices not implemented yet")
    if "Monoclinic" in jl:
        lattice = Lattice.monoclinic(jl["a"], jl["b"], jl["c"], jl["beta"])
    elif "Orthorhombic" in jl:
        lattice = Lattice.orthorhombic(jl["a"], jl["b"], jl["c"])
    elif "Cubic" in jl:
        lattice = Lattice.cubic(jl["a"])
    elif "Hexagonal" in jl:
        lattice = Lattice.hexagonal(jl["a"], jl["c"])
    elif "Rhombohedral" in jl:
        lattice = Lattice.rhombohedral(jl["a"], jl["alpha"])
    elif "Tetragonal" in jl:
        lattice = Lattice.tetragonal(jl["a"], jl["c"])
    elif "Triclinic" in jl:
        lattice = Lattice.from_parameters(jl["a"], jl["b"], jl["c"], jl["alpha"], jl["beta"], jl["gamma"])
    else:
        raise ValueError(f"Unable to convert JDFTX lattice tag {jl} to pymatgen Lattice object")
    return lattice


def _infile_matrix_format_to_pmg_lattice(jdftxinfile: JDFTXInfile) -> Lattice:
    """Convert JDFTX lattice tag to pymatgen Lattice object.

    Args:
        jl (dict[str, Any]): JDFTX lattice tag.

    Returns:
        Lattice: Pymatgen Lattice object.
    """
    _lattice = np.zeros([3, 3])
    jl = jdftxinfile["lattice"]
    for i in range(3):
        for j in range(3):
            _lattice[i][j] += float(jl[f"R{i}{j}"])
    if "latt-scale" in jdftxinfile:
        latt_scale = np.array([jdftxinfile["latt-scale"][x] for x in ["s0", "s1", "s2"]])
        _lattice *= latt_scale
    _lattice = _lattice.T  # convert to row vector format
    _lattice *= const.value("Bohr radius") * 10**10  # Bohr radius in Ang; convert to Ang
    return Lattice(_lattice)
