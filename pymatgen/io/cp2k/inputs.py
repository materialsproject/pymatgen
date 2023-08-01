"""
This module defines the building blocks of a CP2K input file. The cp2k input structure is
essentially a collection of "sections" which are similar to dictionary objects that activate
modules of the cp2k executable, and then "keywords" which adjust variables inside of those
modules. For example, FORCE_EVAL section will activate CP2K's ability to calculate forces,
and inside FORCE_EVAL, the Keyword "METHOD can be set to "QS" to set the method of force
evaluation to be the quickstep (DFT) module.

A quick overview of the module:

-- Section class defines the basis of Cp2k input and contains methods for manipulating these
   objects similarly to Dicts.
-- Keyword class defines the keywords used inside of Section objects that changes variables in
   Cp2k programs.
-- SectionList and KeywordList classes are lists of Section and Keyword objects that have
   the same dictionary key. This deals with repeated sections and keywords.
-- Cp2kInput class is special instantiation of Section that is used to represent the full cp2k
   calculation input.
-- The rest of the classes are children of Section intended to make initialization of common
   sections easier.

"""

from __future__ import annotations

import copy
import hashlib
import itertools
import os
import re
import textwrap
import typing
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Literal, Sequence

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.periodic_table import Element
from pymatgen.io.cp2k.utils import chunk, postprocessor, preprocessor
from pymatgen.io.vasp.inputs import Kpoints as VaspKpoints
from pymatgen.io.vasp.inputs import KpointsSupportedModes
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from pymatgen.core.lattice import Lattice
    from pymatgen.core.structure import Molecule, Structure

__author__ = "Nicholas Winner"
__version__ = "2.0"
__email__ = "nwinner@berkeley.edu"
__date__ = "September 2022"

MODULE_DIR = Path(__file__).resolve().parent


class Keyword(MSONable):
    """
    Class representing a keyword argument in CP2K. Within CP2K Sections, which activate features
    of the CP2K code, the keywords are arguments that control the functionality of that feature.
    For example, the section "FORCE_EVAL" activates the evaluation of forces/energies, but within
    "FORCE_EVAL" the keyword "METHOD" controls whether or not this will be done with, say,
    "Quickstep" (DFT) or "EIP" (empirical interatomic potential).
    """

    def __init__(
        self,
        name: str,
        *values,
        description: str | None = None,
        units: str | None = None,
        verbose: bool | None = True,
        repeats: bool | None = False,
    ):
        """
        Initializes a keyword. These Keywords and the value passed to them are sometimes as simple
        as KEYWORD VALUE, but can also be more elaborate such as KEYWORD [UNITS] VALUE1 VALUE2,
        which is why this class exists: to handle many values and control easy printing to an
        input file.

        Args:
            name: The name of this keyword. Must match an acceptable keyword from CP2K
            values: All non-keyword arguments after 'name' are interpreted as the values to set for
                this keyword. i.e: KEYWORD ARG1 ARG2 would provide two values to the keyword.
            description: The description for this keyword. This can make readability of
                input files easier for some. Default=None.
            units: The units for this keyword. If not specified, CP2K default units will be
                used. Consult manual for default units. Default=None.
            verbose: Whether the description should be printed with the string of this keyword
            repeats: Whether or not this keyword may be repeated. Default=False.
        """
        self.name = name
        self.values = values  # noqa: PD011
        self.description = description
        self.repeats = repeats
        self.units = units
        self.verbose = verbose

    def __str__(self):
        return (
            f"{self.name} {f'[{self.units}] ' if self.units else ''}"
            + " ".join(map(str, self.values))
            + (" ! " + self.description if (self.description and self.verbose) else "")
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Keyword):
            return NotImplemented
        if self.name.upper() == other.name.upper():
            v1 = [_.upper() if isinstance(_, str) else _ for _ in self.values]  # noqa: PD011
            v2 = [_.upper() if isinstance(_, str) else _ for _ in other.values]  # noqa: PD011
            if v1 == v2 and self.units == self.units:
                return True
        return False

    def __add__(self, other):
        return KeywordList(keywords=[self, other])

    def __getitem__(self, item):
        return self.values[item]  # noqa: PD011

    def as_dict(self):
        """Get a dictionary representation of the Keyword."""
        dct = {}
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["name"] = self.name
        dct["values"] = self.values  # noqa: PD011
        dct["description"] = self.description
        dct["repeats"] = self.repeats
        dct["units"] = self.units
        dct["verbose"] = self.verbose
        return dct

    def get_string(self):
        """String representation of Keyword."""
        return str(self)

    @classmethod
    def from_dict(cls, d):
        """Initialize from dictionary."""
        return Keyword(
            d["name"],
            *d["values"],
            description=d["description"],
            repeats=d["repeats"],
            units=d["units"],
            verbose=d["verbose"],
        )

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @staticmethod
    def from_str(s):
        """
        Initialize from a string.

        Keywords must be labeled with strings. If the postprocessor finds
        that the keywords is a number, then None is return (used by
        the file reader).

        Returns:
            Keyword or None
        """
        s = s.strip()
        if "!" in s or "#" in s:
            s, description = re.split("(?:!|#)", s)
            description = description.strip()
        else:
            description = None
        units = re.findall(r"\[(.*)\]", s) or [None]
        s = re.sub(r"\[(.*)\]", "", s)
        args = s.split()
        args = list(map(postprocessor if args[0].upper() != "ELEMENT" else str, args))
        args[0] = str(args[0])
        return Keyword(*args, units=units[0], description=description)

    def verbosity(self, v):
        """Change the printing of this keyword's description."""
        self.verbose = v


class KeywordList(MSONable):
    """
    Some keywords can be repeated, which makes accessing them via the normal dictionary
    methods a little unnatural. This class deals with this by defining a collection
    of same-named keywords that are accessed by one name.
    """

    def __init__(self, keywords: Sequence[Keyword]):
        """
        Initializes a keyword list given a sequence of keywords.

        Args:
            keywords: A list of keywords. Must all have the same name (case-insensitive)
        """
        assert all(k.name.upper() == keywords[0].name.upper() for k in keywords) if keywords else True
        self.name = keywords[0].name if keywords else None
        self.keywords = list(keywords)

    def __str__(self):
        return self.get_string()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        return all(k == o for k, o in zip(self.keywords, other.keywords))

    def __add__(self, other):
        return self.extend(other)

    def __len__(self):
        return len(self.keywords)

    def __getitem__(self, item):
        return self.keywords[item]

    def append(self, item):
        """Append the keyword list."""
        self.keywords.append(item)

    def extend(self, lst: Sequence[Keyword]) -> None:
        """Extend the keyword list."""
        self.keywords.extend(lst)

    def get_string(self, indent=0):
        """String representation of Keyword."""
        return " \n".join("\t" * indent + str(k) for k in self.keywords)

    def verbosity(self, verbosity):
        """Silence all keywords in keyword list."""
        for k in self.keywords:
            k.verbosity(verbosity)


class Section(MSONable):
    """
    Basic input representation of input to Cp2k. Activates functionality inside of the
    Cp2k executable.
    """

    def __init__(
        self,
        name: str,
        subsections: dict | None = None,
        repeats: bool = False,
        description: str | None = None,
        keywords: dict | None = None,
        section_parameters: list | tuple | None = None,
        location: str | None = None,
        verbose: bool | None = True,
        alias: str | None = None,
        **kwargs,
    ):
        """
        Basic object representing a CP2K Section. Sections activate different parts of the
        calculation. For example, FORCE_EVAL section will activate CP2K's ability to calculate
        forces.

        Args:
            name: The name of the section (must match name in CP2K)
            subsections: A dictionary of subsections that are nested in this section.
                Format is {'NAME': Section(*args, **kwargs). The name you chose for 'NAME'
                to index that subsection does not *have* to be the same as the section's true name,
                but we recommend matching them. You can specify a blank dictionary if there are
                no subsections, or if you want to insert the subsections later.
            repeats: Whether or not this section can be repeated. Most sections cannot.
                Default=False.
            description: Description of this section for easier readability
            keywords: the keywords to be set for this section. Each element should be a
                Keyword object. This can be more cumbersome than simply using kwargs for building
                a class in a script, but is more convenient for the class instantiations of CP2K
                sections (see below).
            section_parameters: the section parameters for this section. Section parameters
                are specialized keywords that modify the behavior of the section overall. Most
                sections do not have section parameters, but some do. Unlike normal Keywords,
                these are specified as strings and not as Keyword objects.
            location: the path to the section in the form 'SECTION/SUBSECTION1/SUBSECTION3',
                example for QS module: 'FORCE_EVAL/DFT/QS'. This location is used to automatically
                determine if a subsection requires a supersection to be activated.
            verbose: Controls how much is printed to Cp2k input files (Also see Keyword).
                If True, then a description of the section will be printed with it as a comment
                (if description is set). Default=True.
            alias: An alias for this class to use in place of the name.

            kwargs are interpreted as keyword, value pairs and added to the keywords array as
            Keyword objects
        """
        self.name = name
        self.subsections = subsections if subsections else {}
        self.repeats = repeats
        self.description = description
        keywords = keywords if keywords else {}
        self.keywords = keywords
        self.section_parameters = section_parameters if section_parameters else []
        self.location = location
        self.verbose = verbose
        self.alias = alias
        self.kwargs = kwargs
        for k, v in self.kwargs.items():
            self.keywords[k] = Keyword(k, v)

    def __str__(self):
        return self.get_string()

    def __eq__(self, d):
        d2 = copy.deepcopy(d)
        s2 = copy.deepcopy(self)
        d2.silence()
        s2.silence()
        return d2.as_dict() == s2.as_dict()

    def __deepcopy__(self, memodict=None):
        c = copy.deepcopy(self.as_dict())
        return getattr(__import__(c["@module"], globals(), locals(), c["@class"], 0), c["@class"]).from_dict(
            copy.deepcopy(self.as_dict())
        )

    def __getitem__(self, d):
        r = self.get_keyword(d)
        if not r:
            r = self.get_section(d)
        if r:
            return r
        raise KeyError

    def __add__(self, other):
        if isinstance(other, (Keyword, KeywordList)):
            if other.name in self.keywords:
                self.keywords[other.name] += other
            else:
                self.keywords[other.name] = other
        elif isinstance(other, (Section, SectionList)):
            self.insert(other)
        else:
            TypeError("Can only add sections or keywords.")

    def __setitem__(self, key, value):
        self.setitem(key, value)

    def setitem(self, key, value, strict=False):
        """
        Helper function for setting items. Kept separate from the double-underscore function so that
        "strict" option can be made possible.

        strict will only set values for items that already have a key entry (no insertion).
        """
        if isinstance(value, (Section, SectionList)):
            if key in self.subsections:
                self.subsections[key] = copy.deepcopy(value)
            elif not strict:
                self.insert(value)
        else:
            if not isinstance(value, (Keyword, KeywordList)):
                value = Keyword(key, value)
            match = [k for k in self.keywords if key.upper() == k.upper()]
            if match:
                del self.keywords[match[0]]
                self.keywords[key] = value
            elif not strict:
                self.keywords[key] = value

    def __delitem__(self, key):
        """
        Delete section with name matching key OR delete all keywords
        with names matching this key.
        """
        lst = [sub_sec for sub_sec in self.subsections if sub_sec.upper() == key.upper()]
        if lst:
            del self.subsections[lst[0]]
            return
        lst = [kw for kw in self.keywords if kw.upper() == key.upper()]
        if lst:
            del self.keywords[lst[0]]
            return
        raise KeyError("No section or keyword matching the given key.")

    def __sub__(self, other):
        return self.__delitem__(other)

    def add(self, other):
        """Add another keyword to the current section."""
        assert isinstance(other, (Keyword, KeywordList))
        self + other

    def get(self, d, default=None):
        """
        Similar to get for dictionaries. This will attempt to retrieve the
        section or keyword matching d. Will not raise an error if d does not
        exist.

        Args:
             d: the key to retrieve, if present
             default: what to return if d is not found
        """
        kw = self.get_keyword(d)
        if kw:
            return kw
        sec = self.get_section(d)
        if sec:
            return sec
        return default

    def get_section(self, d, default=None):
        """
        Get function, only for subsections.

        Args:
            d: Name of section to get
            default: return if d is not found in subsections
        """
        for k, v in self.subsections.items():
            if str(k).upper() == str(d).upper():
                return v
        return default

    def get_keyword(self, d, default=None):
        """
        Get function, only for subsections.

        Args:
            d: Name of keyword to get
            default: return if d is not found in keyword list
        """
        for k, v in self.keywords.items():
            if str(k).upper() == str(d).upper():
                return v
        return default

    def update(self, d: dict, strict=False):
        """
        Update the Section according to a dictionary argument. This is most useful
        for providing user-override settings to default parameters. As you pass a
        dictionary the class variables like "description", "location", or "repeats"
        are not included. Therefore, it is recommended that this be used to modify
        existing Section objects to a user's needs, but not used for the creation
        of new Section child-classes.

        Args:
            d (dict): A dictionary containing the update information. Should use nested dictionaries
                to specify the full path of the update. If a section or keyword does not exist, it
                will be created, but only with the values that are provided in "d", not using
                default values from a Section object.
                Example: {
                    'SUBSECTION1': {
                        'SUBSEC2': {'NEW_KEYWORD': 'NEW_VAL'},
                        'NEW_SUBSEC': {'NEW_KWD': 'NEW_VAL'}
                        }
                    }

            strict (bool): If true, only update existing sections and keywords. If false, allow
                new sections and keywords. Default: False
        """
        Section._update(self, d, strict=strict)

    @staticmethod
    def _update(d1, d2, strict=False):
        """Helper method for self.update(d) method (see above)."""
        for k, v in d2.items():
            if isinstance(v, (str, float, int, bool)):
                d1.setitem(k, Keyword(k, v), strict=strict)
            elif isinstance(v, (Keyword, KeywordList)):
                d1.setitem(k, v, strict=strict)
            elif isinstance(v, dict):
                tmp = [_ for _ in d1.subsections if k.upper() == _.upper()]
                if not tmp:
                    if strict:
                        continue
                    d1.insert(Section(k, subsections={}))
                    Section._update(d1.subsections[k], v, strict=strict)
                else:
                    Section._update(d1.subsections[tmp[0]], v, strict=strict)
            elif isinstance(v, Section):
                if not strict:
                    d1.insert(v)
            else:
                raise TypeError(f"Unrecognized type: {type(v)}")

    def set(self, d: dict):
        """Alias for update. Used by custodian."""
        self.update(d)

    def safeset(self, d: dict):
        """Alias for update with strict (no insertions). Used by custodian."""
        self.update(d, strict=True)

    def unset(self, d: dict):
        """Dict based deletion. Used by custodian."""
        for k, v in d.items():
            if isinstance(v, (str, float, int, bool)):
                del self[k][v]
            elif isinstance(v, (Keyword, Section, KeywordList, SectionList)):
                del self[k][v.name]
            elif isinstance(v, dict):
                self[k].unset(v)
            else:
                TypeError("Can only add sections or keywords.")

    def inc(self, d: dict):
        """Mongo style dict modification. Include."""
        for k, v in d.items():
            if isinstance(v, (str, float, bool, int, list)):
                v = Keyword(k, v)
            if isinstance(v, (Keyword, Section, KeywordList, SectionList)):
                self.add(v)
            elif isinstance(v, dict):
                self[k].inc(v)
            else:
                TypeError("Can only add sections or keywords.")

    def insert(self, d):
        """Insert a new section as a subsection of the current one."""
        assert isinstance(d, (Section, SectionList))
        self.subsections[d.alias or d.name] = copy.deepcopy(d)

    def check(self, path: str):
        """
        Check if section exists within the current using a path. Can be useful for cross-checking
        whether or not required dependencies have been satisfied, which CP2K does not enforce.

        Args:
            path (str): Path to section of form 'SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST'
        """
        _path = path.split("/")
        s = self.subsections
        for p in _path:
            tmp = [_ for _ in s if p.upper() == _.upper()]
            if tmp:
                s = s[tmp[0]].subsections
            else:
                return False
        return True

    def by_path(self, path: str):
        """
        Access a sub-section using a path. Used by the file parser.

        Args:
            path (str): Path to section of form 'SUBSECTION1/SUBSECTION2/SUBSECTION_OF_INTEREST'

        """
        _path = path.split("/")
        if _path[0].upper() == self.name.upper():
            _path = _path[1:]
        s = self
        for p in _path:
            s = s.get_section(p)
        return s

    def get_string(self):
        """Get string representation of Section."""
        return Section._get_string(self)

    @staticmethod
    def _get_string(d, indent=0):
        """
        Helper function to return a pretty string of the section. Includes indentation and
        descriptions (if present).
        """
        string = ""
        if d.description and d.verbose:
            filled = textwrap.fill(
                d.description,
                initial_indent="\t" * indent + "! ",
                subsequent_indent="\t" * indent + "! ",
                width=50,
            )
            string += f"\n{filled}\n"
        string += "\t" * indent + "&" + d.name
        string += f" {' '.join(map(str, d.section_parameters))}\n"

        for v in d.keywords.values():
            if isinstance(v, KeywordList):
                string += v.get_string(indent=indent + 1) + "\n"
            else:
                string += "\t" * (indent + 1) + v.get_string() + "\n"
        for v in d.subsections.values():
            string += v._get_string(v, indent + 1)
        string += "\t" * indent + "&END " + d.name + "\n"

        return string

    def verbosity(self, verbosity):
        """
        Change the section verbosity recursively by turning on/off the printing of descriptions.
        Turning off descriptions may reduce the appealing documentation of input files, but also
        helps de-clutter them.
        """
        self.verbose = verbosity
        for v in self.keywords.values():
            v.verbosity(verbosity)
        for v in self.subsections.values():
            v.verbosity(verbosity)

    def silence(self):
        """Recursively delete all print sections so that only defaults are printed out."""
        if self.subsections:
            if self.subsections.get("PRINT"):
                del self.subsections["PRINT"]
            for v in self.subsections.values():
                v.silence()


class SectionList(MSONable):
    """Section list."""

    def __init__(self, sections: Sequence[Section]):
        """
        Initializes a SectionList object using a sequence of sections.

        Args:
            sections: A list of keywords. Must all have the same name (case-insensitive)
        """
        assert all(k.name.upper() == sections[0].name.upper() for k in sections) if sections else True
        self.name = sections[0].name if sections else None
        self.alias = sections[0].alias if sections else None
        self.sections = list(sections)

    def __str__(self):
        return self.get_string()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SectionList):
            return NotImplemented
        return all(k == o for k, o in zip(self.sections, other.sections))

    def __add__(self, other):
        self.append(other)
        return self

    def __len__(self):
        return len(self.sections)

    def __getitem__(self, item):
        return self.sections[item]

    def __deepcopy__(self, memodict=None):
        return SectionList(sections=[d.__deepcopy__() for d in self.sections])

    @staticmethod
    def _get_string(d, indent=0):
        return " \n".join(s._get_string(s, indent) for s in d)

    def get_string(self):
        """Return string representation of section list."""
        return SectionList._get_string(self.sections)

    def get(self, d, index=-1):
        """
        Get for section list. If index is specified, return the section at that index.
        Otherwise, return a get on the last section.
        """
        return self.sections[index].get(d)

    def append(self, item) -> None:
        """Append the section list."""
        self.sections.append(item)

    def extend(self, lst: list) -> None:
        """Extend the section list."""
        self.sections.extend(lst)

    def verbosity(self, verbosity) -> None:
        """Silence all sections in section list."""
        for k in self.sections:
            k.verbosity(verbosity)


class Cp2kInput(Section):
    """
    Special instance of 'Section' class that is meant to represent the overall cp2k input.
    Distinguishes itself from Section by overriding get_string() to not print this section's
    title and by implementing the file i/o.
    """

    def __init__(self, name: str = "CP2K_INPUT", subsections: dict | None = None, **kwargs):
        """Initialize Cp2kInput by calling the super."""
        self.name = name
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = "CP2K Input"
        super().__init__(
            name,
            repeats=False,
            description=description,
            section_parameters=[],
            subsections=subsections,
            **kwargs,
        )

    def get_string(self):
        """Get string representation of the Cp2kInput."""
        s = ""
        for v in self.subsections.values():
            s += v.get_string()
        return s

    @classmethod
    def _from_dict(cls, d):
        """Initialize from a dictionary."""
        return Cp2kInput(
            "CP2K_INPUT",
            subsections=getattr(
                __import__(d["@module"], globals(), locals(), d["@class"], 0),
                d["@class"],
            )
            .from_dict(d)
            .subsections,
        )

    @staticmethod
    def from_file(file: str):
        """Initialize from a file."""
        with zopen(file, "rt") as f:
            txt = preprocessor(f.read(), os.path.dirname(f.name))
            return Cp2kInput.from_str(txt)

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @staticmethod
    def from_str(s: str):
        """Initialize from a string."""
        lines = s.splitlines()
        lines = [line.replace("\t", "") for line in lines]
        lines = [line.strip() for line in lines]
        lines = [line for line in lines if line]
        return Cp2kInput.from_lines(lines)

    @classmethod
    def from_lines(cls, lines: list | tuple):
        """Helper method to read lines of file."""
        cp2k_input = Cp2kInput("CP2K_INPUT", subsections={})
        Cp2kInput._from_lines(cp2k_input, lines)
        return cp2k_input

    def _from_lines(self, lines):
        """Helper method, reads lines of text to get a Cp2kInput."""
        current = self.name
        description = ""
        for line in lines:
            if line.startswith(("!", "#")):
                description += line[1:].strip()
            elif line.upper().startswith("&END"):
                current = "/".join(current.split("/")[:-1])
            elif line.startswith("&"):
                name, subsection_params = line.split()[0][1:], line.split()[1:]
                subsection_params = (
                    []
                    if len(subsection_params) == 1 and subsection_params[0].upper() in ("T", "TRUE", "F", "FALSE", "ON")
                    else subsection_params
                )
                alias = name + " " + " ".join(subsection_params) if subsection_params else None
                s = Section(
                    name,
                    section_parameters=subsection_params,
                    alias=alias,
                    subsections={},
                    description=description,
                )
                description = ""
                tmp = self.by_path(current).get_section(s.alias or s.name)
                if tmp:
                    if isinstance(tmp, SectionList):
                        self.by_path(current)[s.alias or s.name].append(s)
                    else:
                        self.by_path(current)[s.alias or s.name] = SectionList(sections=[tmp, s])
                else:
                    self.by_path(current).insert(s)
                current = current + "/" + alias if alias else current + "/" + name
            else:
                kwd = Keyword.from_str(line)
                tmp = self.by_path(current).get(kwd.name)
                if tmp:
                    if isinstance(tmp, KeywordList):
                        self.by_path(current).get(kwd.name).append(kwd)
                    elif isinstance(self.by_path(current), SectionList):
                        self.by_path(current)[-1][kwd.name] = KeywordList(keywords=[tmp, kwd])
                    else:
                        self.by_path(current)[kwd.name] = KeywordList(keywords=[kwd, tmp])
                elif isinstance(self.by_path(current), SectionList):
                    self.by_path(current)[-1].keywords[kwd.name] = kwd
                else:
                    self.by_path(current).keywords[kwd.name] = kwd

    def write_file(
        self,
        input_filename: str = "cp2k.inp",
        output_dir: str = ".",
        make_dir_if_not_present: bool = True,
    ):
        """Write input to a file.

        Args:
            input_filename (str, optional): Defaults to "cp2k.inp".
            output_dir (str, optional): Defaults to ".".
            make_dir_if_not_present (bool, optional): Defaults to True.
        """
        if not os.path.isdir(output_dir) and make_dir_if_not_present:
            os.mkdir(output_dir)
        filepath = os.path.join(output_dir, input_filename)
        with open(filepath, "w") as f:
            f.write(self.get_string())


class Global(Section):
    """Controls 'global' settings for cp2k execution such as RUN_TYPE and PROJECT_NAME."""

    def __init__(
        self,
        project_name: str = "CP2K",
        run_type: str = "ENERGY_FORCE",
        keywords: dict | None = None,
        **kwargs,
    ):
        """Initialize the global section.

        Args:
            project_name: Defaults to "CP2K".
            run_type: what type of calculation to run
            keywords: Additional keywords to add
        """
        self.project_name = project_name
        self.run_type = run_type
        keywords = keywords if keywords else {}

        description = (
            "Section with general information regarding which kind of simulation to perform an general settings"
        )

        _keywords = {
            "PROJECT_NAME": Keyword("PROJECT_NAME", project_name),
            "RUN_TYPE": Keyword("RUN_TYPE", run_type),
            "EXTENDED_FFT_LENGTHS": Keyword("EXTENDED_FFT_LENGTHS", True),
        }
        keywords.update(_keywords)
        super().__init__(
            "GLOBAL",
            description=description,
            keywords=keywords,
            subsections={},
            **kwargs,
        )


class ForceEval(Section):
    """Controls the calculation of energy and forces in Cp2k."""

    def __init__(self, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        """Initialize the ForceEval section."""
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}

        description = "Parameters needed to calculate energy and forces and describe the system you want to analyze."

        _keywords = {
            "METHOD": Keyword("METHOD", kwargs.get("METHOD", "QS")),
            "STRESS_TENSOR": Keyword("STRESS_TENSOR", kwargs.get("STRESS_TENSOR", "ANALYTICAL")),
        }
        keywords.update(_keywords)
        super().__init__(
            "FORCE_EVAL",
            repeats=True,
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Dft(Section):
    """Controls the DFT parameters in Cp2k."""

    def __init__(
        self,
        basis_set_filenames: Iterable = ("BASIS_MOLOPT",),
        potential_filename="GTH_POTENTIALS",
        uks: bool = True,
        wfn_restart_file_name: str | None = None,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """Initialize the DFT section.

        Args:
            basis_set_filenames: Name of the file that contains the basis set
                information. Defaults to "BASIS_MOLOPT".
            potential_filename: Name of the file that contains the pseudopotential
                information. Defaults to "GTH_POTENTIALS".
            uks: Whether to run unrestricted Kohn Sham (spin polarized).
                Defaults to True.
            wfn_restart_file_name: Defaults to None.
            keywords: additional keywords to add.
            subsections: Any subsections to initialize with. Defaults to None.
        """
        self.basis_set_filenames = basis_set_filenames
        self.potential_filename = potential_filename
        self.uks = uks
        self.wfn_restart_file_name = wfn_restart_file_name
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}

        description = "Parameter needed by dft programs"

        _keywords = {
            "BASIS_SET_FILE_NAME": KeywordList([Keyword("BASIS_SET_FILE_NAME", k) for k in basis_set_filenames]),
            "POTENTIAL_FILE_NAME": Keyword("POTENTIAL_FILE_NAME", potential_filename),
            "UKS": Keyword(
                "UKS",
                uks,
                description="Whether to run unrestricted Kohn Sham (i.e. spin polarized)",
            ),
        }

        if wfn_restart_file_name:
            _keywords["WFN_RESTART_FILE_NAME"] = Keyword("WFN_RESTART_FILE_NAME", wfn_restart_file_name)

        keywords.update(_keywords)
        super().__init__(
            "DFT",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Subsys(Section):
    """Controls the definition of the system to be simulated."""

    def __init__(self, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        """Initialize the subsys section."""
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "A subsystem: coordinates, topology, molecules and cell"
        super().__init__("SUBSYS", keywords=keywords, description=description, subsections=subsections, **kwargs)


class QS(Section):
    """Controls the quickstep settings (DFT driver)."""

    def __init__(
        self,
        method: str = "GPW",
        eps_default: float = 1e-10,
        eps_pgf_orb: float | None = None,
        extrapolation: str = "ASPC",
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the QS Section.

        Args:
            method ("GPW" | "GAPW"): What DFT methodology to use. GPW (Gaussian Plane Waves) for
                DFT with pseudopotentials or GAPW (Gaussian Augmented Plane Waves) for all
                electron calculations.
            eps_default (float): The default level of convergence accuracy. NOTE: This is a
                global value for all the numerical value of all EPS_* values in QS module.
                It is not the same as EPS_SCF, which sets convergence accuracy of the SCF cycle
                alone.
            eps_pgf_orb: Precision for the overlap matrix. Default is to use sqrt(eps_default)
            extrapolation ("PS" | "ASPC"): Method use for extrapolation. If using
                gamma-point-only calculation, then one should either PS
                or ASPC (ASPC especially for MD runs). See the manual for other options.
            keywords: Additional keywords to add
            subsections: Subsections to initialize with.
        """
        self.method = method
        self.eps_default = eps_default
        self.eps_pgf_orb = eps_pgf_orb
        self.extrapolation = extrapolation
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Parameters needed to set up the Quickstep framework"

        _keywords = {
            "METHOD": Keyword("METHOD", self.method),
            "EPS_DEFAULT": Keyword("EPS_DEFAULT", self.eps_default, description="Base precision level (in Ha)"),
            "EXTRAPOLATION": Keyword(
                "EXTRAPOLATION", self.extrapolation, description="WFN extrapolation between steps"
            ),
        }
        if eps_pgf_orb:
            _keywords["EPS_PGF_ORB"] = Keyword("EPS_PGF_ORB", self.eps_pgf_orb, description="Overlap matrix precision")
        keywords.update(_keywords)
        super().__init__(
            "QS",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Scf(Section):
    """Controls the self consistent field loop."""

    def __init__(
        self,
        max_scf: int = 50,
        eps_scf: float = 1e-6,
        scf_guess: str = "RESTART",
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the Scf section.

        Args:
            max_scf (int): Maximum number of SCF loops before terminating. Defaults to 50.
            eps_scf (float): Convergence criteria for SCF loop. Defaults to 1e-6.
            scf_guess: Initial guess for SCF loop.
                "ATOMIC": Generate an atomic density using the atomic code
                "CORE": Diagonalize the core Hamiltonian for an initial guess.
                "HISTORY_RESTART": Extrapolated from previous RESTART files.
                "MOPAC": Use same guess as MOPAC for semi-empirical methods or a simple
                    diagonal density matrix for other methods.
                "NONE": Skip initial guess (only for NON-SCC DFTB).
                "RANDOM": Use random wavefunction coefficients.
                "RESTART": Use the RESTART file as an initial guess (and ATOMIC if not present).
                "SPARSE": Generate a sparse wavefunction using the atomic code (for OT based
                    methods).
            keywords: Additional keywords
            subsections: Additional subsections
        """
        self.max_scf = max_scf
        self.eps_scf = eps_scf
        self.scf_guess = scf_guess

        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}

        description = "Parameters needed to perform an SCF run."

        _keywords = {
            "MAX_SCF": Keyword("MAX_SCF", max_scf, description="Max number of steps for an inner SCF loop"),
            "EPS_SCF": Keyword("EPS_SCF", eps_scf, description="Convergence threshold for SCF"),
            "SCF_GUESS": Keyword("SCF_GUESS", scf_guess, description="How to initialize the density matrix"),
            "MAX_ITER_LUMO": Keyword(
                "MAX_ITER_LUMO",
                kwargs.get("max_iter_lumo", 400),
                description="Iterations for solving for unoccupied levels when running OT",
            ),
        }
        keywords.update(_keywords)
        super().__init__(
            "SCF",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Mgrid(Section):
    """Controls the multigrid for numerical integration."""

    def __init__(
        self,
        cutoff: int | float = 1200,
        rel_cutoff: int | float = 80,
        ngrids: int = 5,
        progression_factor: int = 3,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the MGRID section.

        Args:
            cutoff: Cutoff energy (in Rydbergs for historical reasons) defining how find of
                Gaussians will be used
            rel_cutoff: The relative cutoff energy, which defines how to map the Gaussians onto
                the multigrid. If the the value is too low then, even if you have a high cutoff
                with sharp Gaussians, they will be mapped to the course part of the multigrid
            ngrids: number of grids to use
            progression_factor: divisor that decides how to map Gaussians the multigrid after
                the highest mapping is decided by rel_cutoff
            keywords: additional keywords
            subsections: additional subsections
        """
        self.cutoff = cutoff
        self.rel_cutoff = rel_cutoff
        self.ngrids = ngrids
        self.progression_factor = progression_factor

        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = (
            "Multigrid information. Multigrid allows for sharp gaussians and diffuse "
            "gaussians to be treated on different grids, where the spacing of FFT integration "
            "points can be tailored to the degree of sharpness/diffusiveness"
        )

        _keywords = {
            "CUTOFF": Keyword("CUTOFF", cutoff, description="Cutoff in [Ry] for finest level of the MG."),
            "REL_CUTOFF": Keyword(
                "REL_CUTOFF",
                rel_cutoff,
                description="Controls which gaussians are mapped to which level of the MG",
            ),
            "NGRIDS": Keyword("NGRIDS", ngrids, description="Number of grid levels in the MG"),
            "PROGRESSION_FACTOR": Keyword("PROGRESSION_FACTOR", progression_factor),
        }
        keywords.update(_keywords)
        super().__init__(
            "MGRID",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Diagonalization(Section):
    """Controls diagonalization settings (if using traditional diagonalization)."""

    def __init__(
        self,
        eps_adapt: float = 0,
        eps_iter: float = 1e-8,
        eps_jacobi: float = 0,
        jacobi_threshold: float = 1e-7,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """Initialize the diagronalization section."""
        self.eps_adapt = eps_adapt
        self.eps_iter = eps_iter
        self.eps_jacobi = eps_jacobi
        self.jacobi_threshold = jacobi_threshold
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        location = "CP2K_INPUT/FORCE_EVAL/DFT/SCF/DIAGONALIZATION"
        description = "Settings for the SCF's diagonalization routines"

        _keywords = {
            "EPS_ADAPT": Keyword("EPS_ADAPT", eps_adapt),
            "EPS_ITER": Keyword("EPS_ITER", eps_iter),
            "EPS_JACOBI": Keyword("EPS_JACOBI", eps_jacobi),
            "JACOBI_THRESHOLD": Keyword("JACOBI_THRESHOLD", jacobi_threshold),
        }
        keywords.update(_keywords)
        super().__init__(
            "DIAGONALIZATION",
            keywords=keywords,
            repeats=False,
            location=location,
            description=description,
            subsections=subsections,
            **kwargs,
        )


class Davidson(Section):
    """Parameters for davidson diagonalization."""

    def __init__(
        self,
        new_prec_each: int = 20,
        preconditioner: str = "FULL_SINGLE_INVERSE",
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Args:
            new_prec_each (int): How often to recalculate the preconditioner.
            preconditioner (str): Preconditioner to use.
                "FULL_ALL": Most effective state selective preconditioner based on diagonalization,
                    requires the ENERGY_GAP parameter to be an underestimate of the HOMO-LUMO gap.
                    This preconditioner is recommended for almost all systems, except very large
                    systems where make_preconditioner would dominate the total computational cost.
                "FULL_KINETIC": Cholesky inversion of S and T, fast construction, robust, use for
                    very large systems.
                "FULL_SINGLE": Based on H-eS diagonalisation, not as good as FULL_ALL, but
                    somewhat cheaper to apply.
                "FULL_SINGLE_INVERSE": Based on H-eS cholesky inversion, similar to FULL_SINGLE
                    in preconditioning efficiency but cheaper to construct, might be somewhat
                    less robust. Recommended for large systems.
                "FULL_S_INVERSE": Cholesky inversion of S, not as good as FULL_KINETIC,
                    yet equally expensive.
                "NONE": skip preconditioning
            keywords: additional keywords
            subsections: additional subsections.
        """
        self.new_prec_each = new_prec_each
        self.preconditioner = preconditioner
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        _keywords = {
            "NEW_PREC_EACH": Keyword("NEW_PREC_EACH", new_prec_each),
            "PRECONDITIONER": Keyword("PRECONDITIONER", preconditioner),
        }
        keywords.update(_keywords)
        super().__init__(
            "DAVIDSON",
            keywords=keywords,
            repeats=False,
            location=None,
            subsections=subsections,
            **kwargs,
        )


class OrbitalTransformation(Section):
    """
    Turns on the Orbital Transformation scheme for diagonalizing the Hamiltonian. Often faster
    and with guaranteed convergence compared to normal diagonalization, but requires the system
    to have a band gap.

    NOTE: OT has poor convergence for metallic systems and cannot use SCF mixing or smearing.
    Therefore, you should not use it for metals or systems with 'small' band gaps. In that
    case, use normal diagonalization
    """

    def __init__(
        self,
        minimizer: str = "CG",
        preconditioner: str = "FULL_ALL",
        algorithm: str = "STRICT",
        rotation: bool = False,
        occupation_preconditioner: bool = False,
        energy_gap: float = -1,
        linesearch: str = "2PNT",
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the OT section.

        Args:
            minimizer: The minimizer to use with the OT method. Default is conjugate gradient
                method, which is more robust, but more well-behaved systems should use DIIS, which
                can be as much as 50% faster.
            preconditioner: Preconditioner to use for OT, FULL_ALL tends to be most robust,
                but is not always most efficient. For difficult systems, FULL_SINGLE_INVERSE can be
                more robust, and is reasonably efficient with large systems. For huge, but well
                behaved, systems, where construction of the preconditioner can take a very long
                time, FULL_KINETIC can be a good choice.
            algorithm: What algorithm to use for OT. 'Strict': Taylor or diagonalization
                based algorithm. IRAC: Orbital Transformation based Iterative Refinement of the
                Approximative Congruence transformation (OT/IR).
            rotation: Introduce additional variables to allow subspace rotations (i.e fractional
                occupations)
            occupation_preconditioner: include the fractional occupation in the preconditioning
            energy_gap: Guess for the band gap. For FULL_ALL, should be smaller than the
                actual band gap, so simply using 0.01 is a robust value. Choosing a larger value
                will help if you start with a bad initial guess though. For FULL_SINGLE_INVERSE,
                energy_gap is treated as a lower bound. Values lower than 0.05 in this case can
                lead to stability issues.
            linesearch (str): From the manual: 1D line search algorithm to be used with the OT
                minimizer, in increasing order of robustness and cost. MINIMIZER CG combined with
                LINESEARCH GOLD should always find an electronic minimum. Whereas the 2PNT
                minimizer is almost always OK, 3PNT might be needed for systems in which successive
                OT CG steps do not decrease the total energy.
            keywords: additional keywords
            subsections: additional subsections
        """
        self.minimizer = minimizer
        self.preconditioner = preconditioner
        self.algorithm = algorithm
        self.rotation = rotation
        self.occupation_preconditioner = occupation_preconditioner
        self.energy_gap = energy_gap
        self.linesearch = linesearch
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}

        description = (
            "Sets the various options for the orbital transformation (OT) method. "
            "Default settings already provide an efficient, yet robust method. Most "
            "systems benefit from using the FULL_ALL preconditioner combined with a small "
            "value (0.001) of ENERGY_GAP. Well-behaved systems might benefit from using "
            "a DIIS minimizer. Advantages: It's fast, because no expensive diagonalization"
            "is performed. If preconditioned correctly, method guaranteed to find minimum. "
            "Disadvantages: Sensitive to preconditioning. A good preconditioner can be "
            "expensive. No smearing, or advanced SCF mixing possible: POOR convergence for "
            "metallic systems."
        )

        _keywords = {
            "MINIMIZER": Keyword("MINIMIZER", minimizer),
            "PRECONDITIONER": Keyword("PRECONDITIONER", preconditioner),
            "ENERGY_GAP": Keyword("ENERGY_GAP", energy_gap),
            "ALGORITHM": Keyword("ALGORITHM", algorithm),
            "LINESEARCH": Keyword("LINESEARCH", linesearch),
            "ROTATION": Keyword("ROTATION", rotation),
            "OCCUPATION_PRECONDITIONER": Keyword("OCCUPATION_PRECONDITIONER", occupation_preconditioner),
        }
        keywords.update(_keywords)
        super().__init__(
            "OT",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Cell(Section):
    """Defines the simulation cell (lattice)."""

    def __init__(self, lattice: Lattice, keywords: dict | None = None, **kwargs):
        """
        Initialize the cell section.

        Args:
            lattice: pymatgen lattice object
            keywords: additional keywords
        """
        self.lattice = lattice
        keywords = keywords if keywords else {}
        description = "Lattice parameters and optional settings for creating a the CELL"

        _keywords = {
            "A": Keyword("A", *lattice.matrix[0]),
            "B": Keyword("B", *lattice.matrix[1]),
            "C": Keyword("C", *lattice.matrix[2]),
        }
        keywords.update(_keywords)
        super().__init__("CELL", description=description, keywords=keywords, subsections={}, **kwargs)


class Kind(Section):
    """Specifies the information for the different atom types being simulated."""

    def __init__(
        self,
        specie: str,
        alias: str | None = None,
        magnetization: float = 0.0,
        basis_set: GaussianTypeOrbitalBasisSet | str | None = "GTH_BASIS",
        potential: GthPotential | str | None = "GTH_POTENTIALS",
        ghost: bool = False,
        aux_basis: GaussianTypeOrbitalBasisSet | str | None = None,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize a KIND section.

        Args:
            specie: Object representing the atom.
            alias: Alias for the atom, can be used for specifying modifications
                to certain atoms but not all, e.g. Mg_1 and Mg_2 to force difference
                oxidation states on the two atoms.
            magnetization: From the CP2K Manual: The magnetization used
                in the atomic initial guess. Adds magnetization/2 spin-alpha
                electrons and removes magnetization/2 spin-beta electrons.
            basis_set: Basis set for this atom, accessible from the
                basis set file specified
            potential: Pseudopotential for this atom, accessible from the
                potential file
            ghost: Turn this into ghost atom (disaple the potential)
            aux_basis: Auxiliary basis to use with ADMM
            keywords: additional keywords
            subsections: additional subsections
            kwargs: Additional kwargs to pass to Section()
        """
        self.name = "KIND"
        self.specie = specie
        self.alias = alias
        self.magnetization = magnetization or 0  # if None, set 0
        self.basis_set = basis_set
        self.potential = potential
        self.ghost = ghost or False  # if None, set False
        self.aux_basis = aux_basis
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "The description of this kind of atom including basis sets, element, etc."

        # Special case for closed-shell elements. Cannot impose magnetization in cp2k.
        if Element(self.specie).Z in {
            2,
            4,
            10,
            12,
            18,
            20,
            30,
            36,
            38,
            48,
            54,
            56,
            70,
            80,
            86,
            88,
            102,
            112,
            118,
        }:
            self.magnetization = 0

        _keywords = {
            "ELEMENT": Keyword("ELEMENT", specie.__str__()),
            "MAGNETIZATION": Keyword("MAGNETIZATION", magnetization),
            "GHOST": Keyword("GHOST", ghost),
        }
        if basis_set:
            _keywords["BASIS_SET"] = (
                Keyword("BASIS_SET", basis_set) if isinstance(basis_set, str) else basis_set.get_keyword()
            )
        if potential:
            _keywords["POTENTIAL"] = (
                Keyword("POTENTIAL", potential) if isinstance(potential, str) else potential.get_keyword()
            )
        if aux_basis:
            _keywords["BASIS_SET"] += (
                Keyword("BASIS_SET", f"BASIS_SET AUX_FIT {aux_basis}")
                if isinstance(aux_basis, str)
                else aux_basis.get_keyword()
            )

        kind_name = alias if alias else specie.__str__()
        alias = kind_name

        section_parameters = [kind_name]
        location = "FORCE_EVAL/SUBSYS/KIND"
        keywords.update(_keywords)
        super().__init__(
            name=self.name,
            subsections=subsections,
            description=description,
            keywords=keywords,
            section_parameters=section_parameters,
            alias=alias,
            location=location,
            verbose=True,
            **kwargs,
        )


class DftPlusU(Section):
    """Controls DFT+U for an atom kind."""

    def __init__(
        self,
        eps_u_ramping=1e-5,
        init_u_ramping_each_scf=False,
        l=-1,  # noqa: E741
        u_minus_j=0,
        u_ramping=0,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the DftPlusU section.

        Args:
            eps_u_ramping: (float) SCF convergence threshold at which to start ramping the U value
            init_u_ramping_each_scf: (bool) Whether or not to do u_ramping each scf cycle
            l: (int) angular moment of the orbital to apply the +U correction
            u_minus_j: (float) the effective U parameter, Ueff = U-J
            u_ramping: (float) stepwise amount to increase during ramping until u_minus_j is reached
            keywords: additional keywords
            subsections: additional subsections
        """
        name = "DFT_PLUS_U"
        self.eps_u_ramping = 1e-5
        self.init_u_ramping_each_scf = False
        self.l = l
        self.u_minus_j = u_minus_j
        self.u_ramping = u_ramping
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Settings for on-site Hubbard +U correction for this atom kind."

        _keywords = {
            "EPS_U_RAMPING": Keyword("EPS_U_RAMPING", eps_u_ramping),
            "INIT_U_RAMPING_EACH_SCF": Keyword("INIT_U_RAMPING_EACH_SCF", init_u_ramping_each_scf),
            "L": Keyword("L", l),
            "U_MINUS_J": Keyword("U_MINUS_J", u_minus_j),
            "U_RAMPING": Keyword("U_RAMPING", u_ramping),
        }
        keywords.update(_keywords)
        super().__init__(name=name, subsections=None, description=description, keywords=keywords, **kwargs)


class Coord(Section):
    """Specifies the coordinates of the atoms using a pymatgen structure object."""

    def __init__(
        self,
        structure: Structure | Molecule,
        aliases: dict | None = None,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure object
            alias (bool): whether or not to identify the sites by Element + number so you can do
                things like assign unique magnetization do different elements.
            keywords: additional keywords
            subsections: additional subsections.
        """
        self.structure = structure
        self.aliases = aliases
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = (
            "The coordinates for simple systems (like small QM cells) are specified "
            "here by default using explicit XYZ coordinates. More complex systems "
            "should be given via an external coordinate file in the SUBSYS%TOPOLOGY section."
        )

        if aliases:
            aliases = {index: kind for kind, indices in aliases.items() for index in indices}

        keywords = {
            i: Keyword(aliases[i] if aliases else structure[i].species_string, *structure[i].coords)
            for i in range(len(structure))
        }
        super().__init__(
            name="COORD",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class DOS(Section):
    """Controls printing of the density of states."""

    def __init__(self, ndigits: int = 6, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        """
        Initialize the DOS section.

        Args:
            ndigits: how many digits of precision to print. As of 2022.1,
                this is necessary to not lose information.
            keywords: additional keywords
            subsections: additional subsections
        """
        self.ndigits = ndigits
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Controls printing of the overall density of states"
        _keywords = {"NDIGITS": Keyword("NDIGITS", ndigits)}
        keywords.update(_keywords)
        super().__init__("DOS", description=description, keywords=keywords, subsections=subsections, **kwargs)


class PDOS(Section):
    """
    Controls printing of projected density of states onto the different atom KINDS
    (elemental decomposed DOS).
    """

    def __init__(self, nlumo: int = -1, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        """
        Initialize the PDOS section.

        Args:
            nlumo: how many unoccupied orbitals to include (-1==ALL)
            keywords: additional keywords
            subsections: additional subsections
        """
        self.nlumo = nlumo
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Controls printing of the projected density of states"

        _keywords = {
            "NLUMO": Keyword("NLUMO", nlumo),
            "COMPONENTS": Keyword("COMPONENTS"),
        }
        keywords.update(_keywords)
        super().__init__("PDOS", description=description, keywords=keywords, subsections=subsections, **kwargs)


class LDOS(Section):
    """Controls printing of the LDOS (List-Density of states). i.e. projects onto specific atoms."""

    def __init__(
        self,
        index: int = 1,
        alias: str | None = None,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """
        Initialize the LDOS section.

        Args:
            index: Index of the atom to project onto
            alias: section alias
            keywords: additional keywords
            subsections: additional subsections
        """
        self.index = index
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Controls printing of the projected density of states decomposed by atom type"
        _keywords = {"COMPONENTS": Keyword("COMPONENTS"), "LIST": Keyword("LIST", index)}
        keywords.update(_keywords)
        super().__init__(
            "LDOS",
            subsections=subsections,
            alias=alias,
            description=description,
            keywords=keywords,
            **kwargs,
        )


class V_Hartree_Cube(Section):
    """Controls printing of the hartree potential as a cube file."""

    def __init__(self, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            "Note: by convention the potential has opposite sign than the expected physical one."
        )
        super().__init__(
            "V_HARTREE_CUBE",
            subsections=subsections,
            description=description,
            keywords=keywords,
            **kwargs,
        )


class MO_Cubes(Section):
    """Controls printing of the molecular orbital eigenvalues."""

    def __init__(
        self,
        write_cube: bool = False,
        nhomo: int = 1,
        nlumo: int = 1,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        """Initialize the MO_CUBES section."""
        self.write_cube = write_cube
        self.nhomo = nhomo
        self.nlumo = nlumo
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            "Note: by convention the potential has opposite sign than the expected physical one."
        )

        _keywords = {
            "WRITE_CUBES": Keyword("WRITE_CUBE", write_cube),
            "NHOMO": Keyword("NHOMO", nhomo),
            "NLUMO": Keyword("NLUMO", nlumo),
        }
        keywords.update(_keywords)
        super().__init__(
            "MO_CUBES",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs,
        )


class E_Density_Cube(Section):
    """Controls printing of the electron density cube file."""

    def __init__(self, keywords: dict | None = None, subsections: dict | None = None, **kwargs):
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = (
            "Controls the printing of cube files with the electronic density and, for LSD "
            "calculations, the spin density."
        )

        super().__init__(
            "E_DENSITY_CUBE",
            subsections=subsections,
            description=description,
            keywords=keywords,
            **kwargs,
        )


class Smear(Section):
    """Control electron smearing."""

    def __init__(
        self,
        elec_temp: int | float = 300,
        method: str = "FERMI_DIRAC",
        fixed_magnetic_moment: float = -1e2,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        self.elec_temp = elec_temp
        self.method = method
        self.fixed_magnetic_moment = fixed_magnetic_moment
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        description = "Activates smearing of electron occupations"

        _keywords = {
            "ELEC_TEMP": Keyword("ELEC_TEMP", elec_temp),
            "METHOD": Keyword("METHOD", method),
            "FIXED_MAGNETIC_MOMENT": Keyword("FIXED_MAGNETIC_MOMENT", fixed_magnetic_moment),
        }
        keywords.update(_keywords)
        super().__init__(
            "SMEAR",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class BrokenSymmetry(Section):
    """
    Define the required atomic orbital occupation assigned in initialization
    of the density matrix, by adding or subtracting electrons from specific
    angular momentum channels. It works only with GUESS ATOMIC.
    """

    def __init__(
        self,
        l_alpha: Sequence = (-1,),
        n_alpha: Sequence = (0,),
        nel_alpha: Sequence = (-1,),
        l_beta: Sequence = (-1,),
        n_beta: Sequence = (0,),
        nel_beta: Sequence = (-1,),
    ):
        """
        Initialize the broken symmetry section.

        Args:
            l_alpha: Angular momentum quantum number of the orbitals whose occupation is changed
            n_alpha: Principal quantum number of the orbitals whose occupation is changed.
                Default is the first not occupied
            nel_alpha: Orbital occupation change per angular momentum quantum number. In
                unrestricted calculations applied to spin alpha
            l_beta: Same as L_alpha for beta channel
            n_beta: Same as N_alpha for beta channel
            nel_beta: Same as NEL_alpha for beta channel
        """
        self.l_alpha = l_alpha
        self.n_alpha = n_alpha
        self.nel_alpha = nel_alpha
        self.l_beta = l_beta
        self.n_beta = n_beta
        self.nel_beta = nel_beta
        description = (
            "Define the required atomic orbital occupation assigned in initialization "
            "of the density matrix, by adding or subtracting electrons from specific "
            "angular momentum channels. It works only with GUESS ATOMIC"
        )

        keywords_alpha = {
            "L": Keyword("L", *map(int, l_alpha)),
            "N": Keyword("N", *map(int, n_alpha)),
            "NEL": Keyword("NEL", *map(int, nel_alpha)),
        }
        alpha = Section("ALPHA", keywords=keywords_alpha, subsections={}, repeats=False)

        keywords_beta = {
            "L": Keyword("L", *map(int, l_beta)),
            "N": Keyword("N", *map(int, n_beta)),
            "NEL": Keyword("NEL", *map(int, nel_beta)),
        }
        beta = Section("BETA", keywords=keywords_beta, subsections={}, repeats=False)

        super().__init__(
            "BS",
            description=description,
            subsections={"ALPHA": alpha, "BETA": beta},
            keywords={},
            repeats=False,
        )

    @classmethod
    def from_el(cls, el, oxi_state=0, spin=0):
        """Create section from element, oxidation state, and spin."""
        el = el if isinstance(el, Element) else Element(el)

        def f(x):
            return {"s": 0, "p": 1, "d": 2, "f": 4}.get(x)

        def f2(x):
            return {0: 2, 1: 6, 2: 10, 3: 14}.get(x)

        def f3(x):
            return {0: 2, 1: 6, 2: 10, 3: 14}.get(x)

        es = el.electronic_structure
        esv = [(int(_[0]), f(_[1]), int(_[2:])) for _ in es.split(".") if "[" not in _]
        esv.sort(key=lambda x: (x[0], x[1]), reverse=True)

        tmp = oxi_state
        l_alpha = []
        l_beta = []
        nel_alpha = []
        nel_beta = []
        n_alpha = []
        n_beta = []
        unpaired_orbital = None
        while tmp:
            tmp2 = -min((esv[0][2], tmp)) if tmp > 0 else min((f2(esv[0][1]) - esv[0][2], -tmp))
            l_alpha.append(esv[0][1])
            l_beta.append(esv[0][1])
            nel_alpha.append(tmp2)
            nel_beta.append(tmp2)
            n_alpha.append(esv[0][0])
            n_beta.append(esv[0][0])
            tmp += tmp2
            unpaired_orbital = esv[0][0], esv[0][1], esv[0][2] + tmp2
            esv.pop(0)

        if spin == "low-up":
            spin = unpaired_orbital[2] % 2
        elif spin == "low-down":
            spin = -(unpaired_orbital[2] % 2)
        elif spin == "high-up":
            spin = unpaired_orbital[2] % (f2(unpaired_orbital[1]) // 2)
        elif spin == "high-down":
            spin = -(unpaired_orbital[2] % (f2(unpaired_orbital[1]) // 2))

        if spin:
            for i in reversed(range(len(nel_alpha))):
                nel_alpha[i] += min((spin, f3(l_alpha[i]) - oxi_state))
                nel_beta[i] -= min((spin, f3(l_beta[i]) - oxi_state))
                if spin > 0:
                    spin -= min((spin, f3(l_alpha[i]) - oxi_state))
                else:
                    spin += min((spin, f3(l_beta[i]) - oxi_state))

        return BrokenSymmetry(
            l_alpha=l_alpha,
            l_beta=l_beta,
            nel_alpha=nel_alpha,
            nel_beta=nel_beta,
            n_beta=n_beta,
            n_alpha=n_alpha,
        )


class Xc_Functional(Section):
    """Defines the XC functional(s) to use."""

    def __init__(
        self,
        functionals: Iterable | None = None,
        keywords: dict | None = None,
        subsections: dict | None = None,
        **kwargs,
    ):
        self.functionals = functionals if functionals else []
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}
        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL"

        for functional in self.functionals:
            subsections[functional] = Section(functional, subsections={}, repeats=False)

        super().__init__(
            "XC_FUNCTIONAL",
            subsections=subsections,
            keywords=keywords,
            location=location,
            repeats=False,
            **kwargs,
        )


class PBE(Section):
    """Info about the PBE functional."""

    def __init__(
        self,
        parameterization: str = "ORIG",
        scale_c: float | int = 1,
        scale_x: float | int = 1,
        keywords: dict | None = None,
        subsections: dict | None = None,
    ):
        """
        Args:
            parameterization (str):
                ORIG: original PBE
                PBESOL: PBE for solids/surfaces
                REVPBE: revised PBE
            scale_c (float): scales the correlation part of the functional.
            scale_x (float): scales the exchange part of the functional.
            keywords: additional keywords
            subsections: additional subsections.
        """
        self.parameterization = parameterization
        self.scale_c = scale_c
        self.scale_x = scale_x
        keywords = keywords if keywords else {}
        subsections = subsections if subsections else {}

        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL/PBE"

        _keywords = {
            "PARAMETRIZATION": Keyword("PARAMETRIZATION", parameterization),
            "SCALE_C": Keyword("SCALE_C", scale_c),
            "SCALE_X": Keyword("SCALE_X", scale_x),
        }
        keywords.update(_keywords)
        super().__init__(
            "PBE",
            subsections=subsections,
            repeats=False,
            location=location,
            section_parameters=[],
            keywords=keywords,
        )


class Kpoints(Section):
    """Description of the k-points to use for the calculation."""

    def __init__(
        self,
        kpts: Sequence | Sequence[Sequence[int]],
        weights: Sequence | None = None,
        eps_geo: float = 1e-6,
        full_grid: bool = False,
        parallel_group_size: int = -1,
        scheme: str = "MONKHORST-PACK",
        symmetry: bool = False,
        units: str = "B_VECTOR",
        verbose: bool = False,
        wavefunctions: str = "COMPLEX",
    ):
        """
        Args:
            kpts (list, tuple): a 2D array for the kpoints of the form
                [(1,1,1),]. If len(kpts) == 1. Then it is taken as subdivisions
                for automatic kpoint scheme. If it has more entries, it is
                taken as manual entries for kpoints.
            weights (list, tuple): a weight for each kpoint. Default is to
                weigh each by 1
            eps_geo (float): tolerance for symmetry. Default=1e-6
            full_grid (bool): use full (not reduced) kpoint grid. Default=False.
            parallel_group_size (int): from cp2k manual: Number of processors
                to be used for a single kpoint. This number must divide the
                total number of processes. The number of groups must divide
                the total number of kpoints. Value=-1 (smallest possible
                number of processes per group, satisfying the constraints).
                Value=0 (all processes). Value=n (exactly n processes).
                Default=-1.
            scheme (str): kpoint generation scheme. Default='Monkhorst-Pack'
            symmetry (bool): Use symmetry to reduce the number of kpoints.
                Default=False.
            units (str): Units for the kpoint coordinates (reciprocal coordinates
                or cartesian). Default='B_VECTOR' (reciprocal)
            verbose (bool): verbose output for kpoints. Default=False
            wavefunctions (str): Whether to use complex or real valued wavefunctions
                (if available). Default='complex'.
        """
        description = "Sets up the kpoints"
        keywords = {}

        self.kpts = kpts
        self.weights = weights if weights else [1] * len(kpts)
        assert len(self.kpts) == len(self.weights)
        self.eps_geo = eps_geo
        self.full_grid = full_grid
        self.parallel_group_size = parallel_group_size
        self.scheme = scheme
        self.symmetry = symmetry
        self.units = units
        self.verbose = verbose
        self.wavefunctions = wavefunctions

        if len(kpts) == 1:
            keywords["SCHEME"] = Keyword("SCHEME", scheme, *kpts[0])
        elif len(kpts) > 1:
            keywords["KPOINT"] = KeywordList([Keyword("KPOINT", *k, w) for k, w in zip(self.kpts, self.weights)])
        else:
            raise ValueError("No k-points provided!")

        keywords.update(
            {
                "SCHEME": Keyword("SCHEME", scheme),
                "EPS_GEO": Keyword("EPS_GEO", eps_geo),
                "FULL_GRID": Keyword("FULL_GRID", full_grid),
                "PARALLEL_GROUP_SIZE": Keyword("PARALLEL_GROUP_SIZE", parallel_group_size),
                "SYMMETRY": Keyword("SYMMETRY", symmetry),
                "UNITS": Keyword("UNITS", units),
                "VERBOSE": Keyword("VERBOSE", verbose),
                "WAVEFUNCTIONS": Keyword("WAVEFUNCTIONS", wavefunctions),
            }
        )

        super().__init__(
            name="KPOINTS",
            subsections=None,
            repeats=False,
            description=description,
            keywords=keywords,
        )

    @classmethod
    def from_kpoints(cls, kpoints: VaspKpoints, structure=None):
        """
        Initialize the section from a Kpoints object (pymatgen.io.vasp.inputs). CP2K
        does not have an automatic gamma-point constructor, so this is generally used
        to get the number of divisions from a kpoint static constructor and then
        build a Monkhorst-Pack grid, which is sufficient for gamma-recommended systems
        so long as the grid is fine enough.

        Args:
            kpoints: A pymatgen kpoints object.
            structure: Pymatgen structure object. Required for automatically performing
                symmetry analysis and reducing the kpoint grid.
            reduce: whether or not to reduce the grid using symmetry. CP2K itself cannot
                do this automatically without spglib present at execution time.
        """
        kpts = kpoints.kpts
        weights = kpoints.kpts_weights

        if kpoints.style == KpointsSupportedModes.Monkhorst:
            k = kpts[0]
            if isinstance(k, (int, float)):
                x, y, z = k, k, k
            else:
                x, y, z = k
            scheme = f"MONKHORST-PACK {x} {y} {z}"
            units = "B_VECTOR"
        elif kpoints.style == KpointsSupportedModes.Reciprocal:
            units = "B_VECTOR"
            scheme = "GENERAL"
        elif kpoints.style == KpointsSupportedModes.Cartesian:
            units = "CART_ANGSTROM"
            scheme = "GENERAL"
        elif kpoints.style == KpointsSupportedModes.Gamma:
            if (isinstance(kpts[0], Iterable) and tuple(kpts[0]) == (1, 1, 1)) or (
                isinstance(kpts[0], (float, int)) and int(kpts[0]) == 1
            ):
                scheme = "GAMMA"
                units = "B_VECTOR"
            elif not structure:
                raise ValueError(
                    "No cp2k automatic gamma constructor. A structure is required to construct from spglib"
                )
            else:
                sga = SpacegroupAnalyzer(structure)
                _kpts, weights = zip(*sga.get_ir_reciprocal_mesh(mesh=kpts))
                kpts = list(itertools.chain.from_iterable(_kpts))
                scheme = "GENERAL"
                units = "B_VECTOR"
        elif kpoints.style == KpointsSupportedModes.Line_mode:
            scheme = "GENERAL"
            units = "B_VECTOR"
        else:
            raise ValueError("Unrecognized k-point style")
        return Kpoints(kpts=kpts, weights=weights, scheme=scheme, units=units)


class Kpoint_Set(Section):
    """Specifies a kpoint line to be calculated between special points."""

    def __init__(self, npoints: int, kpoints: Iterable, units: str = "B_VECTOR") -> None:
        """
        Args:
            npoints (int): Number of kpoints along the line.
            kpoints: A dictionary of {label: kpoint} kpoints defining the path
            units (str): Units for the kpoint coordinates.
                Options: "B_VECTOR" (reciprocal coordinates)
                         "CART_ANGSTROM" (units of 2*Pi/Angstrom)
                         "CART_BOHR" (units of 2*Pi/Bohr).
        """
        self.npoints = npoints
        self.kpoints = kpoints
        self.units = units

        keywords = {
            "NPOINTS": Keyword("NPOINTS", npoints),
            "UNITS": Keyword("UNITS", units),
            "SPECIAL_POINT": KeywordList(
                [
                    Keyword("SPECIAL_POINT", "Gamma" if label.upper() == "\\GAMMA" else label, *kpt)
                    for label, kpt in kpoints
                ]
            ),
        }

        super().__init__(
            name="KPOINT_SET",
            subsections=None,
            repeats=True,
            description="Specifies a single k-point line for band structure calculations",
            keywords=keywords,
        )


class Band_Structure(Section):
    """Specifies high symmetry paths for outputting the band structure in CP2K."""

    def __init__(
        self,
        kpoint_sets: Sequence[Kpoint_Set],
        filename: str = "BAND.bs",
        added_mos: int = -1,
        keywords: dict | None = None,
        subsections: dict | None = None,
    ):
        """
        Args:
            kpoint_sets: Sequence of Kpoint_Set objects for the band structure calculation.
            filename: Filename for the band structure output
            added_mos: Added (unoccupied) molecular orbitals for the calculation.
            keywords: additional keywords
            subsections: additional subsections.
        """
        self.kpoint_sets = SectionList(kpoint_sets)
        self.filename = filename
        self.added_mos = added_mos
        keywords = keywords if keywords else {}
        _keywords = {
            "FILE_NAME": Keyword("FILE_NAME", filename),
            "ADDED_MOS": Keyword("ADDED_MOS", added_mos),
        }
        keywords.update(_keywords)
        super().__init__(
            name="BAND_STRUCTURE",
            subsections={"KPOINT_SET": self.kpoint_sets},
            repeats=False,
            description="Controls printing of band structure calculation",
            keywords=keywords,
        )

    # TODO kpoints objects are defined in the vasp module instead of a code agnostic module
    # if this changes in the future as other codes are added, then this will need to change
    @staticmethod
    def from_kpoints(kpoints: VaspKpoints, kpoints_line_density=20):
        """
        Initialize band structure section from a line-mode Kpoint object.

        Args:
            kpoints: a kpoint object from the vasp module, which was constructed in line mode
            kpoints_line_density: Number of kpoints along each path
        """
        if kpoints.style == KpointsSupportedModes.Line_mode:

            def pairwise(iterable):
                a = iter(iterable)
                return zip(a, a)

            kpoint_sets = [
                Kpoint_Set(
                    npoints=kpoints_line_density,
                    kpoints=[(lbls[0], kpts[0]), (lbls[1], kpts[1])],
                    units="B_VECTOR",
                )
                for lbls, kpts in zip(pairwise(kpoints.labels), pairwise(kpoints.kpts))
            ]
        elif kpoints.style in (
            KpointsSupportedModes.Reciprocal,
            KpointsSupportedModes.Cartesian,
        ):
            kpoint_sets = [
                Kpoint_Set(
                    npoints=1,
                    kpoints=[("None", kpts) for kpts in kpoints.kpts],
                    units="B_VECTOR" if kpoints.coord_type == "Reciprocal" else "CART_ANGSTROM",
                )
            ]
        else:
            raise ValueError(
                "Unsupported k-point style. Must be line-mode or explicit k-points (reciprocal/cartesian)."
            )
        return Band_Structure(kpoint_sets=kpoint_sets, filename="BAND.bs")


@dataclass
class BasisInfo(MSONable):
    """
    Summary info about a basis set.

    Attributes:
        electrons: Number of electrons
        core: Number of basis functions per core electron
        valence: Number of basis functions per valence electron OR number of exp if it
            is a FIT formatted admm basis
        polarization: Number of polarization functions
        diffuse: Number of added, diffuse/augmentation functions
        cc: Correlation consistent
        pc: Polarization consistent
        sr: Short-range optimized
        molopt: Optimized for molecules/solids
        admm: Whether this is an auxiliary basis set for ADMM
        lri: Whether this is a local resolution of identity auxiliary basis
        contracted: Whether this basis set is contracted
        xc: Exchange correlation functional used for creating this potential
    """

    electrons: int | None = None
    core: int | None = None
    valence: int | None = None
    polarization: int | None = None
    diffuse: int | None = None
    cc: bool | None = False
    pc: bool | None = False
    sr: bool | None = False
    molopt: bool | None = False
    admm: bool | None = False
    lri: bool | None = False
    contracted: bool | None = None
    xc: str | None = None

    def softmatch(self, other):
        """
        Soft matching to see if two basis sets match.

        Will only match those attributes which *are* defined for this basis info object (one way checking)
        """
        if not isinstance(other, BasisInfo):
            return False
        d1 = self.as_dict()
        d2 = other.as_dict()
        return all(not (v is not None and v != d2[k]) for k, v in d1.items())

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string: str) -> BasisInfo:
        """Get summary info from a string."""
        string = string.upper()
        data: dict[str, Any] = {}
        data["cc"] = "CC" in string
        string = string.replace("CC", "")
        data["pc"] = "PC" in string
        string = string.replace("PC", "")
        data["sr"] = "SR" in string
        string = string.replace("SR", "")
        data["molopt"] = "MOLOPT" in string
        string = string.replace("MOLOPT", "")
        for x in ("LDA", "PADE", "MGGA", "GGA", "HF", "PBE0", "PBE", "BP", "BLYP", "B3LYP", "SCAN"):
            if x in string:
                data["xc"] = x
                string = string.replace(x, "")
                break

        tmp = {"S": 1, "D": 2, "T": 3, "Q": 4}
        if "ADMM" in string or "FIT" in string:
            data["admm"] = True
            bool_core = False
            data["contracted"] = "C" in string
            nums = "".join(s for s in string if s.isnumeric())
            data["valence"] = int(nums) if nums else None
        else:
            data["admm"] = False
            if "LRI" in string:
                data["lri"] = True
            bool_core = "V" not in string or "ALL" in string

        data["polarization"] = string.count("P")
        data["diffuse"] = string.count("X")
        string = "#" + string
        for i, s in enumerate(string):
            if s == "Z":
                z = int(tmp.get(string[i - 1], string[i - 1]))
                data["core"] = z if bool_core else None
                data["valence"] = z
            elif s == "P" and string[i - 1].isnumeric():
                data["polarization"] = int(string[i - 1])
            elif s == "X" and string[i - 1].isnumeric():
                data["diffuse"] = int(string[i - 1])
            elif s == "Q" and string[i + 1].isnumeric():
                data["electrons"] = int("".join(_ for _ in string[i + 1 :] if _.isnumeric()))

        if not data["diffuse"]:
            data["diffuse"] = string.count("AUG")

        return cls(**data)


@dataclass
class AtomicMetadata(MSONable):
    """
    Metadata for basis sets and potentials in cp2k.

    Attributes:
        info: Info about this object
        element: Element for this object
        potential: The potential for this object
        name: Name of the object
        alias_names: Optional aliases
        filename: Name of the file containing this object
        version: Version
    """

    info: BasisInfo | PotentialInfo | None = None
    element: Element | None = None
    potential: Literal["All Electron", "Pseudopotential"] | None = None
    name: str | None = None
    alias_names: list = field(default_factory=list)
    filename: str | None = None
    version: str | None = None

    def softmatch(self, other):
        """
        Soft matching to see if a desired basis/potential matches requirements.

        Does soft matching on the "info" attribute first. Then soft matches against the
        element and name/aliases.
        """
        if not isinstance(other, type(self)):
            return False
        if self.info and not self.info.softmatch(other.info):
            return False
        if self.element is not None and self.element != other.element:
            return False
        if self.potential is not None and self.potential != other.potential:
            return False
        this_names = [self.name]
        if self.alias_names:
            this_names.extend(self.alias_names)
        other_names = [other.name]
        if other.alias_names:
            other_names.extend(other.alias_names)
        return all(not (nm is not None and nm not in other_names) for nm in this_names)

    def get_hash(self) -> str:
        """Get a hash of this object."""
        # usedforsecurity=False needed in FIPS mode (Federal Information Processing Standards)
        # https://github.com/materialsproject/pymatgen/issues/2804
        md5 = hashlib.new("md5", usedforsecurity=False)  # hashlib.md5(usedforsecurity=False) is py39+
        md5.update(self.get_string().lower().encode("utf-8"))
        return md5.hexdigest()

    def get_string(self):
        """Get string representation."""
        return str(self)


@dataclass
class GaussianTypeOrbitalBasisSet(AtomicMetadata):
    """
    Model definition of a GTO basis set.

    Attributes:
        info: Cardinality of this basis
        nset: Number of exponent sets
        n: Principle quantum number for each set
        lmax: Maximum angular momentum quantum number for each set
        lmin: Minimum angular momentum quantum number for each set
        nshell: Number of shells for angular momentum l for each set
        exponents: Exponents for each set
        coefficients: Contraction coefficients for each set. Dict[exp->l->shell]
    """

    info: BasisInfo | None = None
    nset: int | None = None
    n: list[int] | None = None
    lmax: list[int] | None = None
    lmin: list[int] | None = None
    nshell: list[dict[int, int]] | None = None
    exponents: list[list[float]] | None = None
    coefficients: list[dict[int, dict[int, dict[int, float]]]] | None = None

    def __post_init__(self) -> None:
        if self.info and self.potential == "All Electron" and self.element:
            self.info.electrons = self.element.Z
        if self.name == "ALLELECTRON":
            self.name = "ALL"  # cp2k won't parse ALLELECTRON for some reason

        def cast(d):
            new = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    v = cast(v)
                new[int(k)] = v
            return new

        if self.nshell:
            self.nshell = [cast(n) for n in self.nshell]
        if self.coefficients:
            self.coefficients = [cast(c) for c in self.coefficients]

    def get_keyword(self) -> Keyword:
        """Convert basis to keyword object."""
        if not self.name:
            raise ValueError("No name attribute. Cannot create keyword")
        vals: Any = []
        if self.info and self.info.admm:
            vals.append("AUX_FIT")
        vals.append(self.name)
        return Keyword("BASIS_SET", *vals)

    @property
    def nexp(self):
        """Number of exponents."""
        return [len(e) for e in self.exponents]

    @typing.no_type_check
    def get_string(self) -> str:
        """Get standard cp2k GTO formatted string."""
        if any(
            getattr(self, x, None) is None
            for x in ("info", "nset", "n", "lmax", "lmin", "nshell", "exponents", "coefficients")
        ):
            raise ValueError("Must have all attributes defined to get string representation")

        out = f"{self.element} {self.name} {' '.join(self.alias_names)}\n"
        out += f"{self.nset}\n"
        for set_index in range(self.nset):
            out += (
                f"{self.n[set_index]} {self.lmin[set_index]} {self.lmax[set_index]} {self.nexp[set_index]} "
                f"{' '.join(map(str, self.nshell[set_index].values()))}\n"
            )
            for exp in self.coefficients[set_index]:
                out += f"\t {self.exponents[set_index][exp]: .14f} "
                for ll in self.coefficients[set_index][exp]:
                    for shell in self.coefficients[set_index][exp][ll]:
                        out += f"{self.coefficients[set_index][exp][ll][shell]: .14f} "
                out += "\n"
        return out

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string: str) -> GaussianTypeOrbitalBasisSet:
        """Read from standard cp2k GTO formatted string."""
        lines = [line for line in string.split("\n") if line]
        firstline = lines[0].split()
        element = Element(firstline[0])
        names = firstline[1:]
        name, aliases = names[0], names[1:]
        _info = BasisInfo.from_str(name).as_dict()
        for alias in aliases:
            for k, v in BasisInfo.from_str(alias).as_dict().items():
                if _info[k] is None:
                    _info[k] = v
        info = BasisInfo.from_dict(_info)
        potential: Literal["All Electron", "Pseudopotential"]
        if any("ALL" in x for x in [name, *aliases]):
            info.electrons = element.Z
            potential = "All Electron"
        else:
            potential = "Pseudopotential"
        nset = int(lines[1].split()[0])
        n = []
        lmin = []
        lmax = []
        nshell = []
        exponents: list[list[float]] = []
        coefficients: list[dict[int, dict[int, dict[int, float]]]] = []

        line_index = 2
        for set_index in range(nset):
            setinfo = lines[line_index].split()
            _n, _lmin, _lmax, _nexp = map(int, setinfo[0:4])
            n.append(_n)
            lmin.append(_lmin)
            lmax.append(_lmax)

            _nshell = map(int, setinfo[4:])
            nshell.append({ll: int(next(_nshell, 0)) for ll in range(_lmin, _lmax + 1)})
            exponents.append([])
            coefficients.append({i: {ll: {} for ll in range(_lmin, _lmax + 1)} for i in range(_nexp)})
            line_index += 1

            for ii in range(_nexp):
                line = lines[line_index].split()
                exponents[set_index].append(float(line[0]))
                coeffs = list(map(float, line[1:]))

                jj = 0
                for ll in range(_lmin, _lmax + 1):
                    for shell in range(nshell[set_index][ll]):
                        coefficients[set_index][ii][ll][shell] = coeffs[jj]
                        jj += 1

                line_index += 1

        return cls(
            element=element,
            name=name,
            alias_names=aliases,
            info=info,
            potential=potential,
            nset=nset,
            n=n,
            lmin=lmin,
            lmax=lmax,
            nshell=nshell,
            exponents=exponents,
            coefficients=coefficients,
        )


@dataclass
class PotentialInfo(MSONable):
    """
    Metadata for this potential.

    Attributes:
        electrons: Total number of electrons
        potential_type: Potential type (e.g. GTH)
        nlcc: Nonlinear core corrected potential
        xc: Exchange correlation functional used for creating this potential
    """

    electrons: int | None = None
    potential_type: str | None = None
    nlcc: bool | None = None
    xc: str | None = None

    def softmatch(self, other):
        """
        Soft matching to see if two potentials match.

        Will only match those attributes which *are* defined for this basis info object (one way checking)
        """
        if not isinstance(other, PotentialInfo):
            return False
        d1 = self.as_dict()
        d2 = other.as_dict()
        return all(not (v is not None and v != d2[k]) for k, v in d1.items())

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string):
        """Get a cp2k formatted string representation."""
        string = string.upper()
        data = {}
        if "NLCC" in string:
            data["nlcc"] = True
        if "GTH" in string:
            data["potential_type"] = "GTH"
        for i, s in enumerate(string):
            if s == "Q" and string[i + 1].isnumeric():
                data["electrons"] = int("".join(_ for _ in string[i + 1 :] if _.isnumeric()))

        for x in ("LDA", "PADA", "MGGA", "GGA", "HF", "PBE0", "PBE", "BP", "BLYP", "B3LYP", "SCAN"):
            if x in string:
                data["xc"] = x
                break

        return cls(**data)


@dataclass
class GthPotential(AtomicMetadata):
    """
    Representation of GTH-type (pseudo)potential.

    Attributes:
        info: Info about this potential
        n_elecs: Number of electrons for each quantum number
        r_loc: Radius of local projectors
        nexp_ppl: Number of the local pseudopotential functions
        c_exp_ppl: Sequence = field(None, description="Coefficients of the local pseudopotential functions
        radii: Radius of the nonlocal part for angular momentum quantum number l defined by the Gaussian
            function exponents alpha_prj_ppnl
        nprj: Number of projectors
        nprj_ppnl: Number of the non-local projectors for the angular momentum quantum number
        hprj_ppnl: Coefficients of the non-local projector functions. Coeff ij for ang momentum l
        )
    """

    info: PotentialInfo
    n_elecs: dict[int, int] | None = None
    r_loc: float | None = None
    nexp_ppl: int | None = None
    c_exp_ppl: Sequence | None = None
    radii: dict[int, float] | None = None
    nprj: int | None = None
    nprj_ppnl: dict[int, int] | None = None
    hprj_ppnl: dict[int, dict[int, dict[int, float]]] | None = None

    def __post_init__(self) -> None:
        if self.potential == "All Electron" and self.element:
            self.info.electrons = self.element.Z
        if self.name == "ALLELECTRON":
            self.name = "ALL"  # cp2k won't parse ALLELECTRON for some reason

        def cast(d):
            new = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    v = cast(v)
                new[int(k)] = v
            return new

        if self.n_elecs:
            self.n_elecs = cast(self.n_elecs)
        if self.radii:
            self.radii = cast(self.radii)
        if self.nprj_ppnl:
            self.nprj_ppnl = cast(self.nprj_ppnl)
        if self.hprj_ppnl:
            self.hprj_ppnl = cast(self.hprj_ppnl)

    def get_keyword(self) -> Keyword:
        """Get keyword object for the potential."""
        if self.name is None:
            raise ValueError("Cannot get keyword without name attribute")

        return Keyword("POTENTIAL", self.name)

    def get_section(self) -> Section:
        """Convert model to a GTH-formatted section object for input files."""
        if self.name is None:
            raise ValueError("Cannot get section without name attribute")

        keywords = {"POTENTIAL": Keyword("", self.get_string())}
        return Section(
            name=self.name,
            section_parameters=None,
            subsections=None,
            description="Manual definition of GTH Potential",
            keywords=keywords,
        )

    @classmethod
    def from_section(cls, section: Section) -> GthPotential:
        """Extract GTH-formatted string from a section and convert it to model."""
        sec = copy.deepcopy(section)
        sec.verbosity(False)
        s = sec.get_string()
        s = [_ for _ in s.split("\n") if not _.startswith("&")]
        s = "\n".join(s)
        return cls.from_str(s)

    def get_string(self):
        """Convert model to a GTH-formatted string."""
        if None in (
            self.info,
            self.n_elecs,
            self.r_loc,
            self.nexp_ppl,
            self.c_exp_ppl,
            self.radii,
            self.nprj,
            self.nprj_ppnl,
            self.hprj_ppnl,
        ):
            raise ValueError("Must initialize all attributes in order to get string")

        out = f"{self.element} {self.name} {' '.join(self.alias_names)}\n"
        out += f"{' '.join(str(self.n_elecs[i]) for i in range(len(self.n_elecs)))}\n"
        out += f"{self.r_loc: .14f} {self.nexp_ppl} "
        for i in range(self.nexp_ppl):
            out += f"{self.c_exp_ppl[i]: .14f} "
        out += "\n"
        out += f"{self.nprj} \n"
        for idx in range(self.nprj):
            total_fill = self.nprj_ppnl[idx] * 20 + 24
            tmp = f"{self.radii[idx]: .14f} {self.nprj_ppnl[idx]: d}"
            out += f"{tmp:>{''}{24}}"
            for i in range(self.nprj_ppnl[idx]):
                k = total_fill - 24 if i == 0 else total_fill
                tmp = " ".join(f"{v: .14f}" for v in self.hprj_ppnl[idx][i].values())
                out += f"{tmp:>{''}{k}}"
                out += "\n"
        return out

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls, string):
        """Initialize model from a GTH formatted string."""
        lines = [line for line in string.split("\n") if line]
        firstline = lines[0].split()
        element, name, aliases = firstline[0], firstline[1], firstline[2:]
        info = PotentialInfo.from_str(name).as_dict()
        for alias in aliases:
            for k, v in PotentialInfo.from_str(alias).as_dict().items():
                if info[k] is None:
                    info[k] = v
        info = PotentialInfo.from_dict(info)
        potential: Literal["All Electron", "Pseudopotential"]
        if any("ALL" in x for x in [name, *aliases]):
            potential = "All Electron"
            info.electrons = Element(element).Z
        else:
            potential = "Pseudopotential"
        nelecs = {i: int(n) for i, n in enumerate(lines[1].split())}
        info.electrons = sum(nelecs.values())  # override, more reliable than name
        thirdline = lines[2].split()
        r_loc, nexp_ppl, c_exp_ppl = (
            float(thirdline[0]),
            int(thirdline[1]),
            list(map(float, thirdline[2:])),
        )
        nprj = int(lines[3].split()[0]) if len(lines) > 3 else 0

        radii = {}
        nprj_ppnl = {}
        hprj_ppnl = {}
        lines = lines[4:]
        i = 0
        ll = 0
        L = 0

        while ll < nprj:
            line = lines[i].split()
            radii[ll] = float(line[0])
            nprj_ppnl[ll] = int(line[1])
            hprj_ppnl[ll] = {x: {} for x in range(nprj_ppnl[ll])}
            line = list(map(float, line[2:]))
            hprj_ppnl[ll][0] = {j: float(ln) for j, ln in enumerate(line)}

            L = 1
            i += 1
            while nprj_ppnl[ll] > L:
                line2 = list(map(float, lines[i].split()))
                hprj_ppnl[ll][L] = {j: float(ln) for j, ln in enumerate(line2)}
                i += 1
                L += 1
            ll += 1

        return cls(
            element=Element(element),
            name=name,
            alias_names=aliases,
            potential=potential,
            n_elecs=nelecs,
            r_loc=r_loc,
            nexp_ppl=nexp_ppl,
            c_exp_ppl=c_exp_ppl,
            info=info,
            radii=radii,
            nprj=nprj,
            nprj_ppnl=nprj_ppnl,
            hprj_ppnl=hprj_ppnl,
        )


@dataclass
class DataFile(MSONable):
    """A data file for a cp2k calc."""

    objects: Sequence | None = None

    @classmethod
    def from_file(cls, fn):
        """Load from a file."""
        with open(fn) as f:
            data = cls.from_str(f.read())
            for obj in data.objects:
                obj.filename = fn
            return data

    @classmethod
    @np.deprecate(message="Use from_str instead")
    def from_string(cls, *args, **kwargs):
        return cls.from_str(*args, **kwargs)

    @classmethod
    def from_str(cls):
        """Initialize from a string."""
        raise NotImplementedError

    def write_file(self, fn):
        """Write to a file."""
        with open(fn, "w") as f:
            f.write(self.get_string())

    def get_string(self):
        """Get string representation."""
        return "\n".join(b.get_string() for b in self.objects)

    def __str__(self):
        return self.get_string()


@dataclass
class BasisFile(DataFile):
    """Data file for basis sets only."""

    @classmethod
    def from_str(cls, string):
        """Initialize from a string representation."""
        basis_sets = [GaussianTypeOrbitalBasisSet.from_str(c) for c in chunk(string)]
        return cls(objects=basis_sets)


@dataclass
class PotentialFile(DataFile):
    """Data file for potentials only."""

    @classmethod
    def from_str(cls, string):
        """Initialize from a string representation."""
        basis_sets = [GthPotential.from_str(c) for c in chunk(string)]
        return cls(objects=basis_sets)
