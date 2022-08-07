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
import itertools
import os
import re
import textwrap
import warnings
from typing import Iterable, Sequence

from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Molecule, Structure
from pymatgen.io.cp2k.utils import _postprocessor, _preprocessor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Nicholas Winner"
__version__ = "1.0"
__email__ = "nwinner@berkeley.edu"
__date__ = "March 2022"


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
        description: str = None,
        units: str = None,
        verbose: bool = True,
        repeats: bool = False,
    ):
        """
        Initializes a keyword. These Keywords and the value passed to them are sometimes as simple
        as KEYWORD VALUE, but can also be more elaborate such as KEYWORD [UNITS] VALUE1 VALUE2,
        which is why this class exists: to handle many values and control easy printing to an
        input file.

        Args:
            name (str): The name of this keyword. Must match an acceptable keyword from CP2K
            args: All non-keyword arguments after 'name' are interpreted as the values to set for
                this keyword. i.e: KEYWORD ARG1 ARG2 would provide two values to the keyword.
            description (str): The description for this keyword. This can make readability of
                input files easier for some. Default=None.
            units (str): The units for this keyword. If not specified, CP2K default units will be
                used. Consult manual for default units. Default=None.
            repeats (bool): Whether or not this keyword may be repeated. Default=False.
        """
        self.name = name
        self.values = values
        self.description = description
        self.repeats = repeats
        self.units = units
        self.verbose = verbose

    def __str__(self):
        return (
            str(self.name)
            + " "
            + (f"[{self.units}] " if self.units else "")
            + " ".join(map(str, self.values))
            + (" ! " + self.description if (self.description and self.verbose) else "")
        )

    def __eq__(self, other):
        if self.name.upper() == other.name.upper():
            v1 = [_.upper() if isinstance(_, str) else _ for _ in self.values]
            v2 = [_.upper() if isinstance(_, str) else _ for _ in other.values]
            if v1 == v2:
                if self.units == self.units:
                    return True
        return False

    def __add__(self, other):
        return KeywordList(keywords=[self, other])

    def __getitem__(self, item):
        return self.values[item]

    def as_dict(self):
        """
        Get a dictionary representation of the Keyword
        """
        d = {}
        d["@module"] = type(self).__module__
        d["@class"] = type(self).__name__
        d["name"] = self.name
        d["values"] = self.values
        d["description"] = self.description
        d["repeats"] = self.repeats
        d["units"] = self.units
        d["verbose"] = self.verbose
        return d

    def get_string(self):
        """
        String representation of Keyword
        """
        return str(self)

    @classmethod
    def from_dict(cls, d):
        """
        Initialise from dictionary
        """
        return Keyword(
            d["name"],
            *d["values"],
            description=d["description"],
            repeats=d["repeats"],
            units=d["units"],
            verbose=d["verbose"],
        )

    @staticmethod
    def from_string(s):
        """
        Initialise from a string.

        Keywords must be labeled with strings. If the postprocessor finds
        that the keywords is a number, then None is return (used by
        the file reader).

        returns:
            Keyword or None
        """
        s = s.strip()
        if ("!" in s) or ("#" in s):
            s, description = re.split("(?:!|#)", s)
            description = description.strip()
        else:
            description = None
        units = re.findall(r"\[(.*)\]", s) or [None]
        s = re.sub(r"\[(.*)\]", "", s)
        args = list(map(_postprocessor, s.split()))
        args[0] = str(args[0])
        return Keyword(*args, units=units[0], description=description)

    def verbosity(self, v):
        """
        Change the printing of this keyword's description.
        """
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
        self.keywords = keywords

    def __str__(self):
        return self.get_string()

    def __eq__(self, other):
        return all(k == o for k, o in zip(self.keywords, other.keywords))

    def __add__(self, other):
        return self.extend(other)

    def __len__(self):
        return len(self.keywords)

    def __getitem__(self, item):
        return self.keywords[item]

    def append(self, item):
        """
        append the keyword list
        """
        self.keywords.append(item)

    def extend(self, l):
        """
        extend the keyword list
        """
        self.keywords.extend(l)

    def get_string(self, indent=0):
        """
        String representation of Keyword
        """
        return " \n".join(["\t" * indent + str(k) for k in self.keywords])

    def verbosity(self, verbosity):
        """
        Silence all keywords in keyword list
        """
        for k in self.keywords:
            k.verbosity(verbosity)


class Section(MSONable):

    """
    Basic input representation of input to Cp2k. Activates functionality inside of the
    Cp2k executable.
    """

    required_sections: tuple = ()
    required_keywords: tuple = ()

    def __init__(
        self,
        name: str,
        subsections: dict = None,
        repeats: bool = False,
        description: str | None = None,
        keywords: dict = None,
        section_parameters: Sequence[str] = None,
        location: str = None,
        verbose: bool = True,
        alias: str = None,
        **kwargs,
    ):
        """
        Basic object representing a CP2K Section. Sections activate different parts of the
        calculation. For example, FORCE_EVAL section will activate CP2K's ability to calculate
        forces.

        Args:
            name (str): The name of the section (must match name in CP2K)
            subsections (dict): A dictionary of subsections that are nested in this section.
                Format is {'NAME': Section(*args, **kwargs). The name you chose for 'NAME'
                to index that subsection does not *have* to be the same as the section's true name,
                but we recommend matching them. You can specify a blank dictionary if there are
                no subsections, or if you want to insert the subsections later.
            repeats (bool): Whether or not this section can be repeated. Most sections cannot.
                Default=False.
            description (str): Description of this section for easier readability
            keywords (list): the keywords to be set for this section. Each element should be a
                Keyword object. This can be more cumbersome than simply using kwargs for building
                a class in a script, but is more convenient for the class instantiations of CP2K
                sections (see below).
            section_parameters (list): the section parameters for this section. Section parameters
                are specialized keywords that modify the behavior of the section overall. Most
                sections do not have section parameters, but some do. Unlike normal Keywords,
                these are specified as strings and not as Keyword objects.
            location (str): the path to the section in the form 'SECTION/SUBSECTION1/SUBSECTION3',
                example for QS module: 'FORCE_EVAL/DFT/QS'. This location is used to automatically
                determine if a subsection requires a supersection to be activated.
            verbose (str): Controls how much is printed to Cp2k input files (Also see Keyword).
                If True, then a description of the section will be printed with it as a comment
                (if description is set). Default=True.

            kwargs are interpreted as keyword, value pairs and added to the keywords array as
            Keyword objects
        """

        self.name = name
        self.subsections = subsections if subsections is not None else {}
        self.repeats = repeats
        self.description = description
        self.keywords = keywords or {}
        self.section_parameters = section_parameters or []
        self.location = location
        self.verbose = verbose
        self.alias = alias
        self.kwargs = kwargs
        for k, v in self.kwargs.items():
            self.keywords[k] = Keyword(k, v)

        for k in self.required_sections:
            if not self.check(k):
                raise UserWarning(f"WARNING: REQUIRED SECTION {k} HAS NOT BEEN INITIALIZED")
        for k in self.required_keywords:
            if k not in self.keywords:
                raise UserWarning(f"WARNING: REQUIRED KEYWORD {k} HAS NOT BEEN PROVIDED")

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

    def __setitem__(self, key, value, strict=False):
        if isinstance(value, (Section, SectionList)):
            if key in self.subsections:
                self.subsections[key] = value.__deepcopy__()
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
        with names matching this key
        """
        l = [s for s in self.subsections if s.upper() == key.upper()]
        if l:
            del self.subsections[l[0]]
            return
        l = [k for k in self.keywords if k.upper() == key.upper()]
        if l:
            del self.keywords[l[0]]
            return
        raise KeyError("No section or keyword matching the given key.")

    def __sub__(self, other):
        return self.__delitem__(other)

    def add(self, other):
        """
        Add another keyword to the current section
        """
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
        r = self.get_keyword(d)
        if r:
            return r
        r = self.get_section(d)
        if r:
            return r
        return default

    def get_section(self, d, default=None):
        """
        Get function, only for subsections

        Args:
            d: Name of section to get
            default: return if d is not found in subsections
        """
        for k in self.subsections:
            if str(k).upper() == str(d).upper():
                return self.subsections[k]
        return default

    def get_keyword(self, d, default=None):
        """
        Get function, only for subsections

        Args:
            d: Name of keyword to get
            default: return if d is not found in keyword list
        """
        for k in self.keywords:
            if str(k).upper() == str(d).upper():
                return self.keywords[k]
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
        """
        Helper method for self.update(d) method (see above).
        """
        for k, v in d2.items():
            if isinstance(v, (str, float, int, bool)):
                d1.__setitem__(k, Keyword(k, v), strict=strict)  # pylint: disable=C2801
            elif isinstance(v, (Keyword, KeywordList)):
                d1.__setitem__(k, v, strict=strict)  # pylint: disable=C2801
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
        """
        Alias for update. Used by custodian.
        """
        self.update(d)

    def safeset(self, d: dict):
        """
        Alias for update with strict (no insertions). Used by custodian.
        """
        self.update(d, strict=True)

    def unset(self, d: dict):
        """
        Dict based deletion. Used by custodian.
        """
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
        """
        Mongo style dict modification. Include.
        """
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
        """
        Insert a new section as a subsection of the current one
        """
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
        """
        Get string representation of Section
        """
        return Section._get_string(self)

    @staticmethod
    def _get_string(d, indent=0):
        """
        Helper function to return a pretty string of the section. Includes indentation and
        descriptions (if present).
        """
        string = ""
        if d.description and d.verbose:
            string += (
                "\n"
                + textwrap.fill(
                    d.description,
                    initial_indent="\t" * indent + "! ",
                    subsequent_indent="\t" * indent + "! ",
                    width=50,
                )
                + "\n"
            )
        string += "\t" * indent + "&" + d.name
        string += " " + " ".join(map(str, d.section_parameters)) + "\n"

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
        Change the section verbossity recursively by turning on/off the printing of descriptions.
        Turning off descriptions may reduce the appealing documentation of input files, but also
        helps de-clutter them.
        """
        self.verbose = verbosity
        for v in self.keywords.values():
            v.verbosity(verbosity)
        for v in self.subsections.values():
            v.verbosity(verbosity)

    def silence(self):
        """
        Recursively delete all print sections so that only defaults are printed out.
        """
        if self.subsections:
            if self.subsections.get("PRINT"):
                del self.subsections["PRINT"]
            for _s in self.subsections:
                self.subsections[_s].silence()


class SectionList(MSONable):

    """
    Section list
    """

    def __init__(self, sections: Sequence[Section]):
        """
        Initializes a SectionList object using a sequence of sections.

        Args:
            sections: A list of keywords. Must all have the same name (case-insensitive)
        """
        assert all(k.name.upper() == sections[0].name.upper() for k in sections) if sections else True
        self.name = sections[0].name if sections else None
        self.alias = sections[0].alias if sections else None
        self.sections = sections

    def __str__(self):
        return self.get_string()

    def __eq__(self, other):
        return all(k == o for k, o in zip(self.sections, other.sections))

    def __add__(self, other):
        return self.extend(other)

    def __len__(self):
        return len(self.sections)

    def __getitem__(self, item):
        return self.sections[item]

    def __deepcopy__(self, memodict=None):
        return SectionList(sections=[d.__deepcopy__() for d in self.sections])

    @staticmethod
    def _get_string(d, indent=0):
        return " \n".join([s._get_string(s, indent) for s in d])

    def get_string(self):
        """
        Return string representation of section list
        """
        return SectionList._get_string(self.sections)

    def get(self, d, index=-1):
        """
        Get for section list. If index is specified, return the section at that index.
        Otherwise, return a get on the last section.
        """
        return self.sections[index].get(d)

    def append(self, item):
        """
        append the keyword list
        """
        self.sections.append(item)

    def extend(self, l):
        """
        extend the keyword list
        """
        self.sections.extend(l)

    def verbosity(self, verbosity):
        """
        Silence all keywords in keyword list
        """
        for k in self.sections:
            k.verbosity(verbosity)


class Cp2kInput(Section):

    """
    Special instance of 'Section' class that is meant to represent the overall cp2k input.
    Distinguishes itself from Section by overriding get_string() to not print this section's
    title and by implementing the file i/o.
    """

    def __init__(self, name: str = "CP2K_INPUT", subsections: dict = None, **kwargs):
        """
        Initialize Cp2kInput by calling the super
        """
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
        """
        Get string representation of the Cp2kInput
        """
        s = ""
        for k in self.subsections:
            s += self.subsections[k].get_string()
        return s

    @classmethod
    def _from_dict(cls, d):
        """
        Initialize from a dictionary
        """

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
        """
        Initialize from a file
        """
        with zopen(file, "rt") as f:
            txt = _preprocessor(f.read(), os.path.dirname(f.name))
            return Cp2kInput.from_string(txt)

    @staticmethod
    def from_string(s: str):
        """
        Initialize from a string
        """
        lines = s.splitlines()
        lines = [line.replace("\t", "") for line in lines]
        lines = [line.strip() for line in lines]
        lines = [line for line in lines if line]
        return Cp2kInput.from_lines(lines)

    @classmethod
    def from_lines(cls, lines: list | tuple):
        """
        Helper method to read lines of file
        """
        cp2k_input = Cp2kInput("CP2K_INPUT", subsections={})
        Cp2kInput._from_lines(cp2k_input, lines)
        return cp2k_input

    def _from_lines(self, lines):
        """
        Helper method, reads lines of text to get a Cp2kInput
        """
        current = self.name
        description = ""
        for line in lines:
            if line.startswith("!") or line.startswith("#"):
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
                kwd = Keyword.from_string(line)
                tmp = self.by_path(current).get(kwd.name)
                if tmp:
                    if isinstance(tmp, KeywordList):
                        self.by_path(current).get(kwd.name).append(kwd)
                    else:
                        if isinstance(self.by_path(current), SectionList):
                            self.by_path(current)[-1][kwd.name] = KeywordList(keywords=[tmp, kwd])
                        else:
                            self.by_path(current)[kwd.name] = KeywordList(keywords=[kwd, tmp])
                else:
                    if isinstance(self.by_path(current), SectionList):
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

    """
    Controls 'global' settings for cp2k execution such as RUN_TYPE and PROJECT_NAME
    """

    def __init__(
        self,
        project_name: str = "CP2K",
        run_type: str = "ENERGY_FORCE",
        **kwargs,
    ):
        """Initialize the global section

        Args:
            project_name (str, optional): Defaults to "CP2K".
            run_type (str, optional) what type of calculation to run
        """

        self.project_name = project_name
        self.run_type = run_type
        self.kwargs = kwargs

        description = (
            "Section with general information regarding which kind of simulation" + " to perform an general settings"
        )

        keywords = {
            "PROJECT_NAME": Keyword("PROJECT_NAME", project_name),
            "RUN_TYPE": Keyword("RUN_TYPE", run_type),
            "EXTENDED_FFT_LENGTHS": Keyword("EXTENDED_FFT_LENGTHS", True),
        }

        super().__init__(
            "GLOBAL",
            description=description,
            keywords=keywords,
            subsections={},
            **kwargs,
        )


class ForceEval(Section):

    """
    Controls the calculation of energy and forces in Cp2k
    """

    def __init__(self, subsections: dict = None, **kwargs):
        """Initialize the ForceEval section

        Args:
            subsections (dict, optional): Defaults to None.
        """

        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = (
            "Parameters needed to calculate energy and forces" + " and describe the system you want to analyze."
        )

        keywords = {
            "METHOD": Keyword("METHOD", kwargs.get("METHOD", "QS")),
            "STRESS_TENSOR": Keyword("STRESS_TENSOR", kwargs.get("STRESS_TENSOR", "ANALYTICAL")),
        }

        super().__init__(
            "FORCE_EVAL",
            repeats=True,
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Dft(Section):

    """
    Controls the DFT parameters in Cp2k
    """

    def __init__(
        self,
        basis_set_filenames="BASIS_MOLOPT",
        potential_filename="GTH_POTENTIALS",
        uks: bool = True,
        wfn_restart_file_name: str = None,
        subsections: dict = None,
        **kwargs,
    ):
        """Initialize the DFT section.

        Args:
            basis_set_filenames (str, optional): Name of the file that contains the basis set
                information. Defaults to "BASIS_MOLOPT".
            potential_filename (str, optional): Name of the file that contains the pseudopotential
                information. Defaults to "GTH_POTENTIALS".
            uks (bool, optional): Whether to run unrestricted Kohn Sham (spin polarized).
                Defaults to True.
            wfn_restart_file_name (str, optional): Defaults to None.
            subsections (dict, optional): Any subsections to initialize with. Defaults to None.
        """

        self.basis_set_filenames = basis_set_filenames
        self.potential_filename = potential_filename
        self.uks = uks
        self.wfn_restart_file_name = wfn_restart_file_name
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = "Parameter needed by dft programs"

        keywords = {
            "BASIS_SET_FILE_NAME": KeywordList([Keyword("BASIS_SET_FILE_NAME", k) for k in basis_set_filenames]),
            "POTENTIAL_FILE_NAME": Keyword("POTENTIAL_FILE_NAME", potential_filename),
            "UKS": Keyword(
                "UKS",
                uks,
                description="Whether to run unrestricted Kohn Sham (i.e. spin polarized)",
            ),
        }

        if wfn_restart_file_name:
            keywords["WFN_RESTART_FILE_NAME"] = Keyword("WFN_RESTART_FILE_NAME", wfn_restart_file_name)

        super().__init__(
            "DFT",
            description=description,
            keywords=keywords,
            subsections=self.subsections,
            **kwargs,
        )


class Subsys(Section):

    """
    Controls the definition of the system to be simulated
    """

    def __init__(self, subsections: dict = None, **kwargs):
        """
        Initialize the subsys section
        """
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs
        description = "A subsystem: coordinates, topology, molecules and cell"
        super().__init__("SUBSYS", description=description, subsections=subsections, **kwargs)


class QS(Section):

    """
    Controls the quickstep settings (DFT driver)
    """

    def __init__(
        self,
        method: str = "GPW",
        eps_default: float = 1e-10,
        extrapolation: str = "ASPC",
        subsections: dict = None,
        **kwargs,
    ):
        """
        Initialize the QS Section

        Args:
            method ("GPW" | "GAPW"): What DFT methodology to use. GPW (Gaussian Plane Waves) for
                DFT with pseudopotentials or GAPW (Gaussian Augmented Plane Waves) for all
                electron calculations.
            eps_default (float): The default level of convergence accuracy. NOTE: This is a
                global value for all the numerical value of all EPS_* values in QS module.
                It is not the same as EPS_SCF, which sets convergence accuracy of the SCF cycle
                alone.
            extrapolation ("PS" | "ASPC"): Method use for extrapolation. If using
                gamma-point-only calculation, then one should either PS
                or ASPC (ASPC especially for MD runs). See the manual for other options.
            subsections (dict): Subsections to initialize with.
        """

        self.method = method
        self.eps_default = eps_default
        self.extrapolation = extrapolation
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = "Parameters needed to set up the Quickstep framework"

        keywords = {
            "METHOD": Keyword("METHOD", method),
            "EPS_DEFAULT": Keyword("EPS_DEFAULT", eps_default, description="Base precision level (in Ha)"),
            "EXTRAPOLATION": Keyword("EXTRAPOLATION", extrapolation, description="WFN extrapolation between steps"),
        }

        super().__init__(
            "QS",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Scf(Section):

    """
    Controls the self consistent field loop
    """

    def __init__(
        self,
        max_scf: int = 50,
        eps_scf: float = 1e-6,
        scf_guess: str = "RESTART",
        subsections: dict = None,
        **kwargs,
    ):
        """
        Initialize the Scf section

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
        """

        self.max_scf = max_scf
        self.eps_scf = eps_scf
        self.scf_guess = scf_guess
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = "Parameters needed to perform an SCF run."

        keywords = {
            "MAX_SCF": Keyword("MAX_SCF", max_scf, description="Max number of steps for an inner SCF loop"),
            "EPS_SCF": Keyword("EPS_SCF", eps_scf, description="Convergence threshold for SCF"),
            "SCF_GUESS": Keyword("SCF_GUESS", scf_guess, description="How to initialize the density matrix"),
            "MAX_ITER_LUMO": Keyword(
                "MAX_ITER_LUMO",
                kwargs.get("max_iter_lumo", 400),
                description="Iterations for solving for unoccupied levels when running OT",
            ),
        }

        super().__init__(
            "SCF",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Mgrid(Section):

    """
    Controls the multigrid for numerical integration
    """

    def __init__(
        self,
        cutoff: int | float = 1200,
        rel_cutoff: int | float = 80,
        ngrids: int = 5,
        progression_factor: int = 3,
        subsections: dict = None,
        **kwargs,
    ):
        """
        Initialize the MGRID section

        Args:
            cutoff: Cutoff energy (in Rydbergs for historical reasons) defining how find of
                Gaussians will be used
            rel_cutoff: The relative cutoff energy, which defines how to map the Gaussians onto
                the multigrid. If the the value is too low then, even if you have a high cutoff
                with sharp Gaussians, they will be mapped to the course part of the multigrid
            ngrids: number of grids to use
            progression_factor: divisor that decides how to map Gaussians the multigrid after
                the highest mapping is decided by rel_cutoff
        """

        self.cutoff = cutoff
        self.rel_cutoff = rel_cutoff
        self.ngrids = ngrids
        self.progression_factor = progression_factor
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = (
            "Multigrid information. Multigrid allows for sharp gaussians and diffuse "
            + "gaussians to be treated on different grids, where the spacing of FFT integration "
            + "points can be tailored to the degree of sharpness/diffusiveness"
        )

        keywords = {
            "CUTOFF": Keyword("CUTOFF", cutoff, description="Cutoff in [Ry] for finest level of the MG."),
            "REL_CUTOFF": Keyword(
                "REL_CUTOFF",
                rel_cutoff,
                description="Controls which gaussians are mapped to which level of the MG",
            ),
            "NGRIDS": Keyword("NGRIDS", ngrids, description="Number of grid levels in the MG"),
            "PROGRESSION_FACTOR": Keyword("PROGRESSION_FACTOR", progression_factor),
        }

        super().__init__(
            "MGRID",
            description=description,
            keywords=keywords,
            subsections=subsections,
            **kwargs,
        )


class Diagonalization(Section):

    """
    Controls diagonalization settings (if using traditional diagonalization).
    """

    def __init__(
        self,
        eps_adapt: float = 0,
        eps_iter: float = 1e-8,
        eps_jacobi: float = 0,
        jacobi_threshold: float = 1e-7,
        subsections: dict = None,
        **kwargs,
    ):
        """
        Initialize the diagronalization section
        """

        self.eps_adapt = eps_adapt
        self.eps_iter = eps_iter
        self.eps_jacobi = eps_jacobi
        self.jacobi_threshold = jacobi_threshold
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs
        self.location = "CP2K_INPUT/FORCE_EVAL/DFT/SCF/DIAGONALIZATION"
        self.description = "Settings for the SCF's diagonalization routines"

        keywords = {
            "EPS_ADAPT": Keyword("EPS_ADAPT", eps_adapt),
            "EPS_ITER": Keyword("EPS_ITER", eps_iter),
            "EPS_JACOBI": Keyword("EPS_JACOBI", eps_jacobi),
            "JACOBI_THRESHOLD": Keyword("JACOBI_THRESHOLD", jacobi_threshold),
        }

        super().__init__(
            "DIAGONALIZATION",
            keywords=keywords,
            repeats=False,
            location=self.location,
            description=self.description,
            subsections=self.subsections,
            **kwargs,
        )


class Davidson(Section):

    """
    Parameters for davidson diagonalization
    """

    def __init__(
        self,
        new_prec_each: int = 20,
        preconditioner: str = "FULL_SINGLE_INVERSE",
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
        """
        self.new_prec_each = new_prec_each
        self.preconditioner = preconditioner

        keywords = {
            "NEW_PREC_EACH": Keyword("NEW_PREC_EACH", new_prec_each),
            "PRECONDITIONER": Keyword("PRECONDITIONER", preconditioner),
        }

        super().__init__(
            "DAVIDSON",
            keywords=keywords,
            repeats=False,
            location=None,
            subsections={},
            **kwargs,
        )


class OrbitalTransformation(Section):

    """
    Turns on the Orbital Transformation scheme for diagonalizing the Hamiltonian. Much faster and
    with guaranteed convergence compared to normal diagonalization, but requires the system
    to have a band gap.

    NOTE: OT has poor convergence for metallic systems and cannot use SCF mixing or smearing.
    Therefore, you should not use it for metals or systems with 'small' band gaps. In that
    case, use normal diagonalization, which will be slower, but will converge properly.
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
        subsections: dict = None,
        **kwargs,
    ):
        """
        Initialize the OT section

        Args:
            minimizer (str): The minimizer to use with the OT method. Default is conjugate gradient
                method, which is more robust, but more well-behaved systems should use DIIS, which
                can be as much as 50% faster.
            preconditioner (str): Preconditioner to use for OT, FULL_ALL tends to be most robust, but is
                not always most efficient. For difficult systems, FULL_SINGLE_INVERSE can be more robust,
                and is reasonably efficient with large systems. For huge, but well behaved, systems,
                where construction of the preconditioner can take a very long time, FULL_KINETIC can be a good
                choice.
            energy_gap (float): Guess for the band gap. For FULL_ALL, should be smaller than the actual band gap,
                so simply using 0.01 is a robust value. Choosing a larger value will help if you start with a bad
                initial guess though. For FULL_SINGLE_INVERSE, energy_gap is treated as a lower bound. Values lower
                than 0.05 in this case can lead to stability issues.
            algorithm (str): What algorithm to use for OT. 'Strict': Taylor or diagonalization based algorithm.
                IRAC: Orbital Transformation based Iterative Refinement of the Approximative Congruence
                transformation (OT/IR).
            linesearch (str): From the manual: 1D line search algorithm to be used with the OT minimizer,
                in increasing order of robustness and cost. MINIMIZER CG combined with LINESEARCH
                GOLD should always find an electronic minimum. Whereas the 2PNT minimizer is almost always OK,
                3PNT might be needed for systems in which successive OT CG steps do not decrease the total energy.
        """

        self.minimizer = minimizer
        self.preconditioner = preconditioner
        self.algorithm = algorithm
        self.rotation = rotation
        self.occupation_preconditioner = occupation_preconditioner
        self.energy_gap = energy_gap
        self.linesearch = linesearch
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = (
            "Sets the various options for the orbital transformation (OT) method. "
            + "Default settings already provide an efficient, yet robust method. Most "
            + "systems benefit from using the FULL_ALL preconditioner combined with a small "
            + "value (0.001) of ENERGY_GAP. Well-behaved systems might benefit from using "
            + "a DIIS minimizer. Advantages: It's fast, because no expensive diagonalization"
            + "is performed. If preconditioned correctly, method guaranteed to find minimum. "
            + "Disadvantages: Sensitive to preconditioning. A good preconditioner can be "
            + "expensive. No smearing, or advanced SCF mixing possible: POOR convergence for "
            + "metallic systems."
        )

        keywords = {
            "MINIMIZER": Keyword("MINIMIZER", minimizer),
            "PRECONDITIONER": Keyword("PRECONDITIONER", preconditioner),
            "ENERGY_GAP": Keyword("ENERGY_GAP", energy_gap),
            "ALGORITHM": Keyword("ALGORITHM", algorithm),
            "LINESEARCH": Keyword("LINESEARCH", linesearch),
            "ROTATION": Keyword("ROTATION", rotation),
            "OCCUPATION_PRECONDITIONER": Keyword("OCCUPATION_PRECONDITIONER", occupation_preconditioner),
        }

        super().__init__(
            "OT",
            description=description,
            keywords=keywords,
            subsections=self.subsections,
            **kwargs,
        )


class Cell(Section):

    """
    Defines the simulation cell (lattice)
    """

    def __init__(self, lattice: Lattice, **kwargs):
        """
        Initialize the cell section.

        Args:
            lattice: pymatgen lattice object
        """

        self.lattice = lattice
        self.kwargs = kwargs

        description = "Lattice parameters and optional settings for creating a the CELL"

        keywords = {
            "A": Keyword("A", *lattice.matrix[0]),
            "B": Keyword("B", *lattice.matrix[1]),
            "C": Keyword("C", *lattice.matrix[2]),
        }

        super().__init__("CELL", description=description, keywords=keywords, subsections={}, **kwargs)


class Kind(Section):

    """
    Specifies the information for the different atom types being simulated.
    """

    def __init__(
        self,
        specie: str,
        alias: str | None = None,
        magnetization: float = 0.0,
        subsections: dict = None,
        basis_set: str = "GTH_BASIS",
        potential: str = "GTH_POTENTIALS",
        ghost: bool = False,
        aux_basis: str | None = None,
        **kwargs,
    ):
        """
        Initialize a KIND section

        Args:
            specie (Species or Element): Object representing the atom.
            alias (str): Alias for the atom, can be used for specifying modifications
                to certain atoms but not all, e.g. Mg_1 and Mg_2 to force difference
                oxidation states on the two atoms.
            magnetization (float): From the CP2K Manual: The magnetization used
                in the atomic initial guess. Adds magnetization/2 spin-alpha
                electrons and removes magnetization/2 spin-beta electrons.
            basis_set (str): Basis set for this atom, accessible from the
                basis set file specified
            potential (str): Pseudopotential for this atom, accessible from the
                potential file
            kwargs: Additional kwargs to pass to Section()
        """

        self.name = "KIND"
        self.specie = specie
        self.alias = alias
        self.magnetization = magnetization
        self.subsections = subsections if subsections else {}
        self.basis_set = basis_set
        self.potential = potential
        self.ghost = ghost
        self.aux_basis = aux_basis
        self.kwargs = kwargs

        self.description = "The description of this kind of atom including basis sets, element, etc."

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

        keywords = {
            "ELEMENT": Keyword("ELEMENT", str(specie)),
            "MAGNETIZATION": Keyword("MAGNETIZATION", magnetization),
            "BASIS_SET": Keyword("BASIS_SET", basis_set),
            "POTENTIAL": Keyword("POTENTIAL", potential),
            "GHOST": Keyword("GHOST", ghost),
        }
        if aux_basis:
            keywords["BASIS_SET"] += Keyword("BASIS_SET", "AUX_FIT", aux_basis)

        kind_name = alias if alias else str(specie)
        self.alias = kind_name

        self.section_parameters = [kind_name]
        self.location = "FORCE_EVAL/SUBSYS/KIND"
        self.verbose = True
        self.repeats = False

        super().__init__(
            name=self.name,
            subsections=self.subsections,
            description=self.description,
            keywords=keywords,
            section_parameters=self.section_parameters,
            alias=self.alias,
            location=self.location,
            verbose=self.verbose,
            **self.kwargs,
        )


class DftPlusU(Section):

    """
    Controls DFT+U for an atom kind
    """

    def __init__(
        self,
        eps_u_ramping=1e-5,
        init_u_ramping_each_scf=False,
        l=-1,
        u_minus_j=0,
        u_ramping=0,
    ):
        """
        Initialize the DftPlusU section.

        Args:
            eps_u_ramping: (float) SCF convergence threshold at which to start ramping the U value
            init_u_ramping_each_scf: (bool) Whether or not to do u_ramping each scf cycle
            l: (int) angular moment of the orbital to apply the +U correction
            u_minus_j: (float) the effective U parameter, Ueff = U-J
            u_ramping: (float) stepwise amount to increase during ramping until u_minus_j is reached
        """

        self.name = "DFT_PLUS_U"
        self.eps_u_ramping = 1e-5
        self.init_u_ramping_each_scf = False
        self.l = l
        self.u_minus_j = u_minus_j
        self.u_ramping = u_ramping
        self.description = "Settings for on-site Hubbard +U correction for this atom kind."

        keywords = {
            "EPS_U_RAMPING": Keyword("EPS_U_RAMPING", eps_u_ramping),
            "INIT_U_RAMPING_EACH_SCF": Keyword("INIT_U_RAMPING_EACH_SCF", init_u_ramping_each_scf),
            "L": Keyword("L", l),
            "U_MINUS_J": Keyword("U_MINUS_J", u_minus_j),
            "U_RAMPING": Keyword("U_RAMPING", u_ramping),
        }

        super().__init__(
            name=self.name,
            subsections=None,
            description=self.description,
            keywords=keywords,
            section_parameters=self.section_parameters,
            alias=None,
            location=None,
        )


class Coord(Section):

    """
    Specifies the coordinates of the atoms using a pymatgen structure object.
    """

    def __init__(
        self,
        structure: Structure | Molecule,
        aliases: dict | None = None,
        subsections: dict = None,
        **kwargs,
    ):
        """
        Args:
            structure: Pymatgen structure object
            alias (bool): whether or not to identify the sites by Element + number so you can do things like
                assign unique magnetization do different elements.
        """

        self.structure = structure
        self.aliases = aliases
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        description = (
            "The coordinates for simple systems (like small QM cells) are specified "
            + "here by default using explicit XYZ coordinates. More complex systems "
            + "should be given via an external coordinate file in the SUBSYS%TOPOLOGY section."
        )
        if aliases:
            keywords = {
                k[0]: KeywordList([Keyword(k[0], *structure[i].coords) for i in k[1]])
                for k in sorted(aliases.items(), key=lambda x: x[1])
            }
        else:
            keywords = {
                ss: KeywordList([Keyword(s.specie.symbol, *s.coords) for s in structure.sites if s.specie.symbol == ss])
                for ss in structure.symbol_set
            }
        super().__init__(
            name="COORD",
            description=description,
            keywords=keywords,
            alias=None,
            subsections=self.subsections,
            **kwargs,
        )


class PDOS(Section):

    """
    Controls printing of projected density of states onto the different atom KINDS
    (elemental decomposed DOS).
    """

    def __init__(self, nlumo: int = -1, **kwargs):
        """
        Initialize the PDOS section

        Args:
            nlumo: how many unoccupied orbitals to include (-1==ALL)
        """

        self.nlumo = nlumo
        self.kwargs = kwargs

        description = "Controls printing of the projected density of states"

        keywords = {
            "NLUMO": Keyword("NLUMO", nlumo),
            "COMPONENTS": Keyword("COMPONENTS"),
        }

        super().__init__("PDOS", description=description, keywords=keywords, subsections={}, **kwargs)


class LDOS(Section):

    """
    Controls printing of the LDOS (List-Density of states). i.e. projects onto specific atoms.
    """

    def __init__(self, index: int = 1, alias: str | None = None, **kwargs):
        """
        Initialize the LDOS section

        Args:
            index: Index of the atom to project onto
        """

        self.index = index
        self.alias = alias
        self.kwargs = kwargs

        description = "Controls printing of the projected density of states decomposed by atom type"

        keywords = {"COMPONENTS": Keyword("COMPONENTS"), "LIST": Keyword("LIST", index)}

        super().__init__(
            "LDOS",
            subsections={},
            alias=alias,
            description=description,
            keywords=keywords,
            **kwargs,
        )


class V_Hartree_Cube(Section):

    """
    Controls printing of the hartree potential as a cube file.
    """

    def __init__(self, keywords=None, **kwargs):
        """
        Initialize the V_HARTREE_CUBE section
        """

        self.keywords = keywords if keywords else {}
        self.kwargs = kwargs

        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            + "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            + "Note that by convention the potential has opposite sign than the expected physical one."
        )

        super().__init__(
            "V_HARTREE_CUBE",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs,
        )


class MO_Cubes(Section):

    """
    Controls printing of the molecular orbital eigenvalues
    """

    def __init__(self, write_cube: bool = False, nhomo: int = 1, nlumo: int = 1, **kwargs):
        """
        Initialize the MO_CUBES section
        """

        self.write_cube = write_cube
        self.nhomo = nhomo
        self.nlumo = nlumo
        self.kwargs = kwargs

        description = (
            "Controls the printing of a cube file with eletrostatic potential generated by "
            + "the total density (electrons+ions). It is valid only for QS with GPW formalism. "
            + "Note that by convention the potential has opposite sign than the expected physical one."
        )

        keywords = {
            "WRITE_CUBES": Keyword("WRITE_CUBE", write_cube),
            "NHOMO": Keyword("NHOMO", nhomo),
            "NLUMO": Keyword("NLUMO", nlumo),
        }

        super().__init__(
            "MO_CUBES",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs,
        )


class E_Density_Cube(Section):

    """
    Controls printing of the electron density cube file
    """

    def __init__(self, keywords=None, **kwargs):
        """
        Initialize the E_DENSITY_CUBE Section
        """

        self.keywords = keywords if keywords else {}
        self.kwargs = kwargs

        description = (
            "Controls the printing of cube files with the electronic density and, for LSD "
            + "calculations, the spin density."
        )

        super().__init__(
            "E_DENSITY_CUBE",
            subsections={},
            description=description,
            keywords=keywords,
            **kwargs,
        )


class Smear(Section):

    """
    Control electron smearing
    """

    def __init__(
        self,
        elec_temp: int | float = 300,
        method: str = "FERMI_DIRAC",
        fixed_magnetic_moment: float = -1e2,
        **kwargs,
    ):
        """
        Initialize the SMEAR section
        """

        self.elec_temp = elec_temp
        self.method = method
        self.fixed_magnetic_moment = fixed_magnetic_moment
        self.kwargs = kwargs

        description = "Activates smearing of electron occupations"

        keywords = {
            "ELEC_TEMP": Keyword("ELEC_TEMP", elec_temp),
            "METHOD": Keyword("METHOD", method),
            "FIXED_MAGNETIC_MOMENT": Keyword("FIXED_MAGNETIC_MOMENT", fixed_magnetic_moment),
        }

        super().__init__(
            "SMEAR",
            description=description,
            keywords=keywords,
            subsections={},
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
        Initialize the broken symmetry section

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
            + "of the density matrix, by adding or subtracting electrons from specific "
            + "angular momentum channels. It works only with GUESS ATOMIC"
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
        """
        Create section from element, oxidation state, and spin.
        """
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
            if tmp > 0:
                tmp2 = -min((esv[0][2], tmp))
            else:
                tmp2 = min((f2(esv[0][1]) - esv[0][2], -tmp))
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

    """
    Defines the XC functional(s) to use.
    """

    def __init__(self, functionals: Iterable = None, subsections: dict = None, **kwargs):
        """
        Initialize the XC_FUNCTIONAL class
        """

        self.functionals = functionals or []
        self.subsections = subsections if subsections else {}
        self.kwargs = kwargs

        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL"

        for functional in self.functionals:
            self.subsections[functional] = Section(functional, subsections={}, repeats=False)

        super().__init__(
            "XC_FUNCTIONAL",
            subsections=self.subsections,
            location=location,
            repeats=False,
            **kwargs,
        )


class PBE(Section):

    """
    Info about the PBE functional.
    """

    def __init__(
        self,
        parameterization: str = "ORIG",
        scale_c: float | int = 1,
        scale_x: float | int = 1,
    ):
        """
        Args:
            parameterization (str):
                ORIG: original PBE
                PBESOL: PBE for solids/surfaces
                REVPBE: revised PBE
            scale_c (float): scales the correlation part of the functional.
            scale_x (float): scales the exchange part of the functional.
        """

        self.parameterization = parameterization
        self.scale_c = scale_c
        self.scale_x = scale_x

        location = "CP2K_INPUT/FORCE_EVAL/DFT/XC/XC_FUNCTIONAL/PBE"

        keywords = {
            "PARAMETRIZATION": Keyword("PARAMETRIZATION", parameterization),
            "SCALE_C": Keyword("SCALE_C", scale_c),
            "SCALE_X": Keyword("SCALE_X", scale_x),
        }

        super().__init__(
            "PBE",
            subsections={},
            repeats=False,
            location=location,
            section_parameters=[],
            keywords=keywords,
        )


class Kpoints(Section):

    """
    Description of the k-points to use for the calculation.
    """

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
                (if available). Default='complex'
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
    def from_kpoints(cls, kpoints, structure=None, reduce=True):
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
        k = kpoints.as_dict()
        kpoints = k["kpoints"]
        weights = k["kpts_weights"]
        scheme = k["generation_style"]

        if reduce and structure:
            sga = SpacegroupAnalyzer(structure)
            kpoints, weights = zip(*sga.get_ir_reciprocal_mesh(mesh=kpoints))
            kpoints = list(itertools.chain.from_iterable(kpoints))
            scheme = "GENERAL"
        elif scheme.lower() == "monkhorst":
            scheme = "MONKHORST-PACK"
        else:
            warnings.warn("No automatic constructor for this scheme" "defaulting to monkhorst-pack grid.")
            scheme = "MONKHORST-PACK"
        units = k["coord_type"]
        if k["coord_type"]:
            if k["coord_type"].lower() == "reciprocal":
                units = "B_VECTOR"
            elif k["coord_type"].lower() == "cartesian":
                units = "CART_ANGSTROM"
        else:
            units = "B_VECTOR"

        return Kpoints(kpts=kpoints, weights=weights, scheme=scheme, units=units)
