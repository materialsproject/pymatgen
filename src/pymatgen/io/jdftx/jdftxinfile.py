"""Classes for reading/manipulating/writing JDFTx input files.

Classes for reading/manipulating/writing JDFTx input files.
"""

from __future__ import annotations

import itertools
import warnings
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import scipy.constants as const
from monty.io import zopen
from monty.json import MSONable
from pymatgen.core import Structure
from pymatgen.util.io_utils import clean_lines

from atomate2.jdftx.io.generic_tags import (
    AbstractTag,
    BoolTagContainer,
    DumpTagContainer,
    TagContainer,
)
from atomate2.jdftx.io.jdftxinfile_master_format import (
    __PHONON_TAGS__,
    __TAG_LIST__,
    __WANNIER_TAGS__,
    MASTER_TAG_LIST,
    get_tag_object,
)

if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import ArrayLike
    from pymatgen.util.typing import PathLike
    from typing_extensions import Self

__author__ = "Jacob Clary"


class JDFTXInfile(dict, MSONable):
    """Class for reading/writing JDFtx input files.

    JDFTxInfile object for reading and writing JDFTx input files.
    Essentially a dictionary with some helper functions.
    """

    path_parent: str = None  # Only gets initialized if from_file

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

        Return str representation of JDFTXInfile.

        Returns
        -------
        str
            String representation of JDFTXInfile.
        """
        return "".join([line + "\n" for line in self.get_text_list()])

    # def __add__(self, other: Self) -> Self:
    def __add__(self, other: JDFTXInfile) -> JDFTXInfile:
        """Add existing JDFTXInfile object to method caller JDFTXInfile object.

        Add all the values of another JDFTXInfile object to this object.
        Facilitate the use of "standard" JDFTXInfiles.

        Parameters
        ----------
        other : JDFTXInfile
            JDFTXInfile object to add to the method caller object.

        Returns
        -------
        JDFTXInfile
        """
        params: dict[str, Any] = dict(self.items())
        for key, val in other.items():
            if key in self and val != self[key]:
                raise ValueError(
                    f"JDFTXInfiles have conflicting values for {key}: \
                        {self[key]} != {val}"
                )
            params[key] = val
        return type(self)(params)

    def as_dict(self, sort_tags: bool = True, skip_module_keys: bool = False) -> dict:
        """Return JDFTXInfile as MSONable dict.

        Return JDFTXInfile as MSONable dict.

        Parameters
        ----------
        sort_tags : bool, optional
            Whether to sort the tags, by default True
        skip_module_keys : bool, optional
            Whether to skip the module keys, by default False

        Returns
        -------
        dict
            JDFTXInfile as MSONable dict
        """
        params = dict(self)
        if sort_tags:
            params = {tag: params[tag] for tag in __TAG_LIST__ if tag in params}
        if not skip_module_keys:
            params["@module"] = type(self).__module__
            params["@class"] = type(self).__name__
        return params

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """Parse a dictionary to create a JDFTXInfile object.

        Parse a dictionary to create a JDFTXInfile object.

        Parameters
        ----------
        dct : dict
            Dictionary to parse.

        Returns
        -------
            JDFTXInfile
        """
        temp = cls({k: v for k, v in dct.items() if k not in ("@module", "@class")})
        # since users can provide arbitrary tags and values, need to do some
        #   validation (could do more later)
        # passing through the list -> dict representation ensures that tags
        # pass through a conversion to string and then through all .read()
        # methods (happens during list->dict conversion) to help ensure correct
        # formatting the list representation is easier to look at so convert
        # back at the end
        temp = cls.get_dict_representation(cls.get_list_representation(temp))
        return cls.get_list_representation(temp)

    def copy(self) -> JDFTXInfile:
        """Return a copy of the JDFTXInfile object.

        Return a copy of the JDFTXInfile object.

        Returns
        -------
        JDFTXInfile
            Copy of the JDFTXInfile object.
        """
        return type(self)(self)

    def get_text_list(self) -> list[str]:
        """Get a list of strings representation of the JDFTXInfile.

        Get a list of strings representation of the JDFTXInfile.

        Returns
        -------
        list[str]
            List of strings representation of the JDFTXInfile.
        """
        self_as_dict = self.get_dict_representation(self)

        text: list[str] = []
        for tag_group in MASTER_TAG_LIST:
            added_tag_in_group = False
            for tag in MASTER_TAG_LIST[tag_group]:
                if tag not in self:
                    continue
                if tag in __WANNIER_TAGS__:
                    raise ValueError("Wannier functionality has not been added!")

                added_tag_in_group = True
                tag_object: AbstractTag = MASTER_TAG_LIST[tag_group][tag]
                if tag_object.can_repeat and isinstance(self_as_dict[tag], list):
                    # if a tag_object.can_repeat, it is assumed that self[tag]
                    #    is a list the 2nd condition ensures this
                    # if it is not a list, then the tag will still be printed by
                    #    the else this could be relevant if someone manually
                    #    sets the tag the can repeat's value to a non-list
                    text += [
                        tag_object.write(tag, entry) for entry in self_as_dict[tag]
                    ]
                    # for entry in self_as_dict[tag]:
                    #     text.append(tag_object.write(tag, entry))
                else:
                    text.append(tag_object.write(tag, self_as_dict[tag]))

            if added_tag_in_group:
                text.append("")
        return text

    def write_file(self, filename: PathLike) -> None:
        """Write JDFTXInfile to a file.

        Write JDFTXInfile to a file.

        Parameters
        ----------
        filename : PathLike
            Filename to write to.
        """
        with zopen(filename, mode="wt") as file:
            file.write(str(self))

    @classmethod
    def from_file(
        cls,
        filename: PathLike,
        dont_require_structure: bool = False,
        sort_tags: bool = True,
        assign_path_parent: bool = True,
    ) -> Self:
        """Read an JDFTXInfile object from a file.

        Read an JDFTXInfile object from a file.

        Parameters
        ----------
        filename : PathLike
            Filename to read from.
        dont_require_structure : bool, optional
            Whether to require structure tags, by default False
        sort_tags : bool, optional
            Whether to sort the tags, by default True
        assign_path_parent : bool, optional
            Whether to assign the parent directory of the input file for include
            tags, by default True

        Returns
        -------
            JDFTXInfile object
        """
        path_parent = None
        if assign_path_parent:
            path_parent = Path(filename).parents[0]
        with zopen(filename, mode="rt") as file:
            return cls.from_str(
                file.read(),
                dont_require_structure=dont_require_structure,
                sort_tags=sort_tags,
                path_parent=path_parent,
            )

    @staticmethod
    def _preprocess_line(line: str) -> tuple[AbstractTag, str, str]:
        """Preprocess a line from a JDFTXInfile.

        Preprocess a line from a JDFTXInfile, splitting it into a tag object,
        tag name, and value.

        Parameters
        ----------
        line : str
            Line from the input file.

        Returns
        -------
        tuple[AbstractTag, str, str]
            Tag object, tag name, and value.
        """
        line_list = line.strip().split(maxsplit=1)
        tag: str = line_list[0].strip()
        if tag in __PHONON_TAGS__:
            raise ValueError("Phonon functionality has not been added!")
        if tag in __WANNIER_TAGS__:
            raise ValueError("Wannier functionality has not been added!")
        if tag not in __TAG_LIST__:
            raise ValueError(
                f"The {tag} tag in {line_list} is not in MASTER_TAG_LIST and is\
                    not a comment, something is wrong with this input data!"
            )
        tag_object = get_tag_object(tag)
        value: str = ""
        if len(line_list) == 2:
            value = line_list[1].strip()
        elif len(line_list) == 1:
            value = (
                ""  # exception for tags where only tagname is used,
                # e.g. dump-only tag
            )
        else:
            raise ValueError(
                f"The len(line.split(maxsplit=1)) of {line_list} should never \
                    not be 1 or 2"
            )

        return tag_object, tag, value

    @staticmethod
    def _store_value(
        params: dict[str, list | list[dict[str, dict]] | Any],
        tag_object: AbstractTag,
        tag: str,
        value: Any,
    ) -> dict:
        """Store the value in the params dictionary.

        Store the value in the params dictionary.

        Parameters
        ----------
        params : dict
            Dictionary to store the value in.
        tag_object : AbstractTag
            Tag object.
        tag : str
            Tag name.
        value : Any
            Value to store.

        Returns
        -------
        dict
        """
        if tag_object.can_repeat:  # store tags that can repeat in a list
            if tag not in params:
                params[tag] = []
            if type(tag_object) not in [DumpTagContainer]:
                params[tag].append(value)
            else:  # The previous if statement will need to adapted to reference
                # a tag object flag to be stored in this manner. This manner
                # is to store all subtags as standalone dictionaries within
                # a list, but to combine alike subtags (ie the same dump freq)
                # as they appear.
                if not isinstance(value, dict):
                    raise ValueError(
                        f"The value for the {tag} tag should be a dictionary!"
                    )
                for freq in value:
                    inserted = False
                    for i, preex in enumerate(params[tag]):
                        if freq in preex:
                            params[tag][i][freq].update(value[freq])
                            inserted = True
                            break
                    if not inserted:
                        params[tag].append(value)
        else:
            if tag in params:
                raise ValueError(
                    f"The '{tag}' tag appears multiple times in this input when\
                        it should not!"
                )
            params[tag] = value
        return params

    @staticmethod
    def _gather_tags(lines: list[str]) -> list[str]:
        """Gather broken lines into single string for processing later.

        Gather all tags broken across lines into single string for processing
        later.

        Parameters
        ----------
        lines : list[str]
            List of lines from the input file.

        Returns
        -------
        gathered_strings : list[str]
            List of strings with tags broken across lines combined into single
            string.
        """
        # gather all tags broken across lines into single string for processing
        # later
        total_tag = ""
        gathered_strings = []
        for line in lines:
            if line[-1] == "\\":  # then tag is continued on next line
                total_tag += (
                    line[:-1].strip() + " "
                )  # remove \ and any extra whitespace
            elif total_tag:  # then finished with line continuations
                total_tag += line
                gathered_strings.append(total_tag)
                total_tag = ""
            else:  # then append line like normal
                gathered_strings.append(line)
        return gathered_strings

    @property
    def structure(self) -> Structure:
        """Return a pymatgen Structure object.

        Return a pymatgen Structure object.

        Returns
        -------
        structure : pymatgen.Structure
            Pymatgen structure object.
        """
        return self.to_pmg_structure(self)

    @classmethod
    def from_str(
        cls,
        string: str,
        dont_require_structure: bool = False,
        sort_tags: bool = True,
        path_parent: Path = None,
    ) -> Self:
        """Read a JDFTXInfile object from a string.

        Read a JDFTXInfile object from a string.

        Parameters
        ----------
        string : str
            String to read from.
        dont_require_structure : bool, optional
            Whether to require structure tags, by default False
        sort_tags : bool, optional
            Whether to sort the tags, by default True
        path_parent : Path, optional
            Path to the parent directory of the input file for include tags,
            by default None

        Returns
        -------
        JDFTXInfile
        """
        lines: list[str] = list(clean_lines(string.splitlines()))
        lines = cls._gather_tags(lines)

        params: dict[str, Any] = {}
        # process all tag value lines using specified tag formats in
        # MASTER_TAG_LIST
        for line in lines:
            tag_object, tag, value = cls._preprocess_line(line)
            processed_value = tag_object.read(tag, value)
            params = cls._store_value(
                params, tag_object, tag, processed_value
            )  # this will change with tag categories

        if "include" in params:
            for filename in params["include"]:
                _filename = filename
                if not Path(_filename).exists():
                    if path_parent is not None:
                        _filename = path_parent / filename
                    if not Path(_filename).exists():
                        raise ValueError(
                            f"The include file {filename} ({_filename}) does not exist!"
                        )
                params.update(
                    cls.from_file(
                        _filename, dont_require_structure=True, assign_path_parent=False
                    )
                )
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
        return cls(params)

    @classmethod
    def to_jdftxstructure(
        cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False
    ) -> JDFTXStructure:
        """Convert JDFTXInfile to JDFTXStructure object.

        Converts JDFTx lattice, lattice-scale, ion tags into JDFTXStructure,
        with Pymatgen structure as attribute.

        Parameters
        ----------
        jdftxinfile : JDFTXInfile
            JDFTXInfile object to convert.
        sort_structure : bool, optional
            Whether to sort the structure. Useful if species are not grouped
            properly together. Defaults to False.
        """
        # use dict representation so it's easy to get the right column for
        # moveScale, rather than checking for velocities
        jdftxinfile_dict = cls.get_dict_representation(jdftxinfile)
        return JDFTXStructure.from_jdftxinfile(
            jdftxinfile_dict, sort_structure=sort_structure
        )

    @classmethod
    def to_pmg_structure(
        cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False
    ) -> Structure:
        """Convert JDFTXInfile to pymatgen Structure object.

        Converts JDFTx lattice, lattice-scale, ion tags into pymatgen Structure.

        Parameters
        ----------
        jdftxinfile : JDFTXInfile
            JDFTXInfile object to convert.
        sort_structure : bool, optional
            Whether to sort the structure. Useful if species are not grouped
            properly together. Defaults to False.

        Returns
        -------
        Structure
        """
        # use dict representation so it's easy to get the right column for
        # moveScale, rather than checking for velocities
        jdftxstructure = JDFTXStructure.from_jdftxinfile(
            jdftxinfile.get_dict_representation(jdftxinfile),
            sort_structure=sort_structure,
        )
        return jdftxstructure.structure
        # JDFTXInfile = cls.get_dict_representation(JDFTXInfile)
        # return JDFTXStructure._from_jdftxinfile(
        #     JDFTXInfile, sort_structure=sort_structure
        # ).structure

    @staticmethod
    def _needs_conversion(
        conversion: str, value: dict | list[dict] | list | list[list]
    ) -> bool:
        """Determine if a value needs to be converted.

        Determine if a value needs to be converted.

        Parameters
        ----------
        conversion : str
            Conversion type.
        value : dict | list[dict] | list | list[list]
            Value to check.

        Returns
        -------
        bool
            Whether the value needs to be converted.
        """
        # value will be in one of these formats:
        #  dict-to-list:
        #    dict
        #    list[dicts] (repeat tags in dict representation)
        #  list-to-dict:
        #    list
        #    list[lists] (repeat tags in list representation or lattice in list
        #                 representation)

        if conversion == "list-to-dict":
            flag = False
        elif conversion == "dict-to-list":
            flag = True

        if isinstance(value, dict) or all(isinstance(x, dict) for x in value):
            return flag
        return not flag

    @classmethod
    def get_list_representation(cls, jdftxinfile: JDFTXInfile) -> JDFTXInfile:
        """Convert JDFTXInfile object properties into list representation.

        Convert JDFTXInfile object properties into list representation.

        Parameters
        ----------
        jdftxinfile : JDFTXInfile
            JDFTXInfile object to convert.
        """
        reformatted_params = deepcopy(jdftxinfile.as_dict(skip_module_keys=True))
        # rest of code assumes lists are lists and not np.arrays
        reformatted_params = {
            k: v.tolist() if isinstance(v, np.ndarray) else v
            for k, v in reformatted_params.items()
        }
        for tag, value in reformatted_params.items():
            tag_object = get_tag_object(tag)
            if all(
                [
                    tag_object.allow_list_representation,
                    tag_object.is_tag_container,
                    cls._needs_conversion("dict-to-list", value),
                ]
            ):
                reformatted_params.update(
                    {tag: tag_object.get_list_representation(tag, value)}
                )
        return cls(reformatted_params)

    @classmethod
    def get_dict_representation(cls, jdftxinfile: JDFTXInfile) -> JDFTXInfile:
        """Convert JDFTXInfile object properties into dict representation.

        Convert JDFTXInfile object properties into dict representation.

        Parameters
        ----------
        jdftxinfile : JDFTXInfile
            JDFTXInfile object to convert.
        """
        reformatted_params = deepcopy(jdftxinfile.as_dict(skip_module_keys=True))
        # Just to make sure only passing lists and no more numpy arrays
        reformatted_params = {
            k: v.tolist() if isinstance(v, np.ndarray) else v
            for k, v in reformatted_params.items()
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
                reformatted_params.update(
                    {tag: tag_object.get_dict_representation(tag, value)}
                )
        return cls(reformatted_params)

    def validate_tags(
        self,
        try_auto_type_fix: bool = False,
        error_on_failed_fix: bool = True,
        return_list_rep: bool = False,
    ) -> None:
        """Validate the tags in the JDFTXInfile.

        Validate the tags in the JDFTXInfile. If try_auto_type_fix is True, will
        attempt to fix the tags. If error_on_failed_fix is True, will raise an
        error if the tags cannot be fixed. If return_list_rep is True, will
        return the tags in list representation.

        Parameters
        ----------
        try_auto_type_fix : bool, optional
            Whether to attempt to fix the tags, by default False
        error_on_failed_fix : bool, optional
            Whether to raise an error if the tags cannot be fixed, by default True
        return_list_rep : bool, optional
            Whether to return the tags in list representation, by default False
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
                raise ValueError(
                    f"The {tag} tag with value:\n{self[tag]}\ncould not be \
                        fixed!"
                )
            if try_auto_type_fix and is_tag_valid:
                self.update({tag: value})
            if should_warn:
                warnmsg = f"The {tag} tag with value:\n{self[tag]}\nhas \
                    incorrect typing!\n    Subtag IsValid?\n"
                if any(
                    isinstance(tag_object, tc)
                    for tc in [TagContainer, DumpTagContainer, BoolTagContainer]
                ):
                    warnmsg += "(Check earlier warnings for more details)\n"
                warnings.warn(warnmsg, stacklevel=2)


@dataclass
class JDFTXStructure(MSONable):
    """Object for representing the data in JDFTXStructure tags.

    Object for representing the data in JDFTXStructure tags.

    Attributes
    ----------
    structure: Structure
        Associated Structure.
    selective_dynamics: ArrayLike
        Selective dynamics attribute for each site if available. Shape Nx1
    sort_structure: bool
        Whether to sort the structure. Useful if species are not grouped
        properly together. Defaults to False.
    """

    structure: Structure = None
    selective_dynamics: ArrayLike | None = None
    sort_structure: bool = False

    def __post_init__(self) -> None:
        """Post init function for JDFTXStructure.

        Post init function for JDFTXStructure. Asserts self.structure is
        ordered, and adds selective dynamics if needed.
        """
        if self.structure.is_ordered:
            site_properties = {}
            if self.selective_dynamics is not None:
                selective_dynamics = np.array(self.selective_dynamics)
                if not selective_dynamics.all():
                    site_properties["selective_dynamics"] = selective_dynamics

            # create new copy of structure so can add selective dynamics and
            # sort atoms if needed
            structure = Structure.from_sites(self.structure)
            self.structure = structure.copy(site_properties=site_properties)
            if self.sort_structure:
                self.structure = self.structure.get_sorted_structure()
        else:
            raise ValueError(
                "Disordered structure with partial occupancies cannot be \
                    converted into JDFTXStructure!"
            )

    def __repr__(self) -> str:
        """Return representation of JDFTXStructure file.

        Return representation of JDFTXStructure file.

        Returns
        -------
        str
            Representation of JDFTXStructure file.
        """
        return f"JDFTXStructure({self.get_str()})"

    def __str__(self) -> str:
        """Return string representation of JDFTXStructure file.

        Return string representation of JDFTXStructure file.

        Returns
        -------
        str
            String representation of JDFTXStructure file.
        """
        return self.get_str()

    @property
    def natoms(self) -> list[int]:
        """Return count for each atom type.

        Return sequence of number of sites of each type associated with
        JDFTXStructure

        Returns
        -------
        list[int]
            Sequence of number of sites of each type associated with
            JDFTXStructure
        """
        syms: list[str] = [site.species.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    @classmethod
    def from_str(cls, data: str) -> JDFTXStructure:
        """Read JDFTXStructure from string.

        Read JDFTXStructure from string.

        Parameters
        ----------
        data : str
            String to read from.

        Returns
        -------
        JDFTXStructure
        """
        return cls.from_jdftxinfile(JDFTXInfile.from_str(data))

    @classmethod
    def from_file(cls, filename: str) -> JDFTXStructure:
        """Read JDFTXStructure from file.

        Read JDFTXStructure from file.

        Parameters
        ----------
        filename : str
            Filename to read from.

        Returns
        -------
        JDFTXStructure
        """
        return cls.from_jdftxinfile(JDFTXInfile.from_file(filename))

    @classmethod
    def from_jdftxinfile(
        cls, jdftxinfile: JDFTXInfile, sort_structure: bool = False
    ) -> JDFTXStructure:
        """Get JDFTXStructure from JDFTXInfile.

        Get JDFTXStructure from JDFTXInfile.

        Parameters
        ----------
        jdftxinfile : JDFTXInfile
            JDFTXInfile object
        sort_structure : bool, optional
            Whether to sort the structure. Useful if species are not grouped
            properly together as JDFTx output will have species sorted.
        """
        lattice = np.array([jdftxinfile["lattice"][x] for x in jdftxinfile["lattice"]])
        if "latt-scale" in jdftxinfile:
            latt_scale = np.array(
                [[jdftxinfile["latt-scale"][x] for x in ["s0", "s1", "s2"]]]
            )
            lattice *= latt_scale
        lattice = lattice.T  # convert to row vector format
        lattice *= (
            const.value("Bohr radius") * 10**10
        )  # Bohr radius in Ang; convert to Ang

        atomic_symbols = [x["species-id"] for x in jdftxinfile["ion"]]
        coords = np.array([[x["x0"], x["x1"], x["x2"]] for x in jdftxinfile["ion"]])
        selective_dynamics = np.array([x["moveScale"] for x in jdftxinfile["ion"]])

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
        return cls(struct, selective_dynamics, sort_structure=sort_structure)

    def get_str(self, in_cart_coords: bool = False) -> str:
        """Return a string to be written as JDFTXInfile tags.

        Return a string to be written as JDFTXInfile tags. Allows extra options
        as compared to calling str(JDFTXStructure) directly.

        Parameters
        ----------
        in_cart_coords: bool
            Whether coordinates are output in direct or Cartesian

        Returns
        -------
        str
            representation of JDFTXInfile structure tags
        """
        jdftx_tag_dict = {}

        lattice = np.copy(self.structure.lattice.matrix)
        lattice = lattice.T  # transpose to get into column-vector format
        lattice /= (
            const.value("Bohr radius") * 10**10
        )  # Bohr radius in Ang; convert to Bohr

        jdftx_tag_dict["lattice"] = lattice
        jdftx_tag_dict["ion"] = []
        for i, site in enumerate(self.structure):
            coords = site.coords if in_cart_coords else site.frac_coords
            if self.selective_dynamics is not None:
                sd = self.selective_dynamics[i]
            else:
                sd = 1
            jdftx_tag_dict["ion"].append([site.label, *coords, sd])

        return str(JDFTXInfile.from_dict(jdftx_tag_dict))

    def write_file(self, filename: PathLike, **kwargs) -> None:
        """Write JDFTXStructure to a file.

        Write JDFTXStructure to file. The supported kwargs are the same as
        those for the JDFTXStructure.get_str method and are passed through
        directly.

        Parameters
        ----------
        filename : str
            Filename to write to.
        **kwargs
            Kwargs to pass to JDFTXStructure.get_str.
        """
        with zopen(filename, mode="wt") as file:
            file.write(self.get_str(**kwargs))

    def as_dict(self) -> dict:
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "structure": self.structure.as_dict(),
            "selective_dynamics": np.array(self.selective_dynamics).tolist(),
        }

    @classmethod
    def from_dict(cls, params: dict) -> Self:
        """Get JDFTXStructure from dict.

        Get JDFTXStructure from dict.

        Parameters
        ----------
        params : dict
            Serialized JDFTXStructure

        Returns
        -------
        JDFTXStructure
        """
        return cls(
            Structure.from_dict(params["structure"]),
            selective_dynamics=params["selective_dynamics"],
        )
