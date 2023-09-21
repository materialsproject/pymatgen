"""Classes and methods related to the Structure Notation Language (SNL)."""

from __future__ import annotations

import datetime
import json
import re
import sys
from collections import namedtuple
from io import StringIO

from monty.json import MontyDecoder, MontyEncoder
from monty.string import remove_non_ascii

try:
    from pybtex import errors
    from pybtex.database.input import bibtex
except ImportError:
    pybtex = bibtex = None

from typing import TYPE_CHECKING

from pymatgen.core.structure import Molecule, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence

__author__ = "Anubhav Jain, Shyue Ping Ong"
__credits__ = "Dan Gunter"


MAX_HNODE_SIZE = 64000  # maximum size (bytes) of SNL HistoryNode
MAX_DATA_SIZE = 256000  # maximum size (bytes) of SNL data field
MAX_HNODES = 100  # maximum number of HistoryNodes in SNL file
MAX_BIBTEX_CHARS = 20000  # maximum number of characters for BibTeX reference


def is_valid_bibtex(reference: str) -> bool:
    """Use pybtex to validate that a reference is in proper BibTeX format.

    Args:
        reference: A String reference in BibTeX format.

    Returns:
        Boolean indicating if reference is valid bibtex.
    """
    # str is necessary since pybtex seems to have an issue with unicode. The
    # filter expression removes all non-ASCII characters.
    sio = StringIO(remove_non_ascii(reference))
    parser = bibtex.Parser()
    errors.set_strict_mode(enable=False)
    bib_data = parser.parse_stream(sio)
    return len(bib_data.entries) > 0


class HistoryNode(namedtuple("HistoryNode", ["name", "url", "description"])):
    """A HistoryNode represents a step in the chain of events that lead to a
    Structure. HistoryNodes leave 'breadcrumbs' so that you can trace back how
    a Structure was created. For example, a HistoryNode might represent pulling
    a Structure from an external database such as the ICSD or CSD. Or, it might
    represent the application of a code (e.g. pymatgen) to the Structure, with
    a custom description of how that code was applied (e.g. a site removal
    Transformation was applied).

    A HistoryNode contains three fields:

    Attributes:
        name (str): The name of a code or resource that this Structure encountered in its history.
        url (str): The URL of that code/resource.
        description (dict): A free-form description of how the code/resource is related to the Structure.
    """

    __slots__ = ()

    def as_dict(self) -> dict[str, str]:
        """Returns: Dict."""
        return {"name": self.name, "url": self.url, "description": self.description}

    @staticmethod
    def from_dict(dct: dict[str, str]) -> HistoryNode:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            HistoryNode
        """
        return HistoryNode(dct["name"], dct["url"], dct["description"])

    @staticmethod
    def parse_history_node(h_node):
        """Parses a History Node object from either a dict or a tuple.

        Args:
            h_node: A dict with name/url/description fields or a 3-element
                tuple.

        Returns:
            History node.
        """
        if isinstance(h_node, dict):
            return HistoryNode.from_dict(h_node)

        if len(h_node) != 3:
            raise ValueError(f"Invalid History node, should be dict or (name, version, description) tuple: {h_node}")
        return HistoryNode(h_node[0], h_node[1], h_node[2])


class Author(namedtuple("Author", ["name", "email"])):
    """An Author contains two fields: name and email. It is meant to represent
    the author of a Structure or the author of a code that was applied to a Structure.
    """

    __slots__ = ()

    def __str__(self):
        """String representation of an Author."""
        return f"{self.name} <{self.email}>"

    def as_dict(self):
        """Returns: MSONable dict."""
        return {"name": self.name, "email": self.email}

    @staticmethod
    def from_dict(d):
        """
        Args:
            d (dict): Dict representation.

        Returns:
            Author
        """
        return Author(d["name"], d["email"])

    @staticmethod
    def parse_author(author):
        """Parses an Author object from either a String, dict, or tuple.

        Args:
            author: A String formatted as "NAME <email@domain.com>",
                (name, email) tuple, or a dict with name and email keys.

        Returns:
            An Author object.
        """
        if isinstance(author, str):
            # Regex looks for whitespace, (any name), whitespace, <, (email),
            # >, whitespace
            m = re.match(r"\s*(.*?)\s*<(.*?@.*?)>\s*", author)
            if not m or m.start() != 0 or m.end() != len(author):
                raise ValueError(f"Invalid author format! {author}")
            return Author(m.groups()[0], m.groups()[1])
        if isinstance(author, dict):
            return Author.from_dict(author)
        if len(author) != 2:
            raise ValueError(f"Invalid author, should be String or (name, email) tuple: {author}")
        return Author(author[0], author[1])


class StructureNL:
    """The Structure Notation Language (SNL, pronounced 'snail') is a container for a pymatgen
    Structure/Molecule object with some additional fields for enhanced provenance.

    It is meant to be imported/exported in a JSON file format with the following structure:
        - sites
        - lattice (optional)
        - about
            - created_at
            - authors
            - projects
            - references
            - remarks
            - data
            - history
    """

    def __init__(
        self,
        struct_or_mol,
        authors,
        projects=None,
        references="",
        remarks=None,
        data=None,
        history=None,
        created_at=None,
    ):
        """
        Args:
            struct_or_mol: A pymatgen Structure/Molecule object
            authors: *List* of {"name":'', "email":''} dicts,
                *list* of Strings as 'John Doe <johndoe@gmail.com>',
                or a single String with commas separating authors
            projects: List of Strings ['Project A', 'Project B']
            references: A String in BibTeX format
            remarks: List of Strings ['Remark A', 'Remark B']
            data: A free form dict. Namespaced at the root level with an
                underscore, e.g. {"_materialsproject": <custom data>}
            history: List of dicts - [{'name':'', 'url':'', 'description':{}}]
            created_at: A datetime object.
        """
        # initialize root-level structure keys
        self.structure = struct_or_mol

        # turn authors into list of Author objects
        authors = authors.split(",") if isinstance(authors, str) else authors
        self.authors = [Author.parse_author(a) for a in authors]

        # turn projects into list of Strings
        projects = projects or []
        self.projects = [projects] if isinstance(projects, str) else projects

        # check that references are valid BibTeX
        if not isinstance(references, str):
            raise ValueError("Invalid format for SNL reference! Should be empty string or BibTeX string.")
        if references and not is_valid_bibtex(references):
            raise ValueError("Invalid format for SNL reference! Should be BibTeX string.")
        if len(references) > MAX_BIBTEX_CHARS:
            raise ValueError(
                f"The BibTeX string must be fewer than {MAX_BIBTEX_CHARS} chars, you have {len(references)}"
            )

        self.references = references

        # turn remarks into list of Strings
        remarks = remarks or []
        self.remarks = [remarks] if isinstance(remarks, str) else remarks

        # check remarks limit
        for remark in self.remarks:
            if len(remark) > 140:
                raise ValueError(f"The remark exceeds the maximum size of 140 characters: {len(remark)}")

        # check data limit
        self.data = data or {}
        if not sys.getsizeof(self.data) < MAX_DATA_SIZE:
            raise ValueError(
                f"The data dict exceeds the maximum size limit of {MAX_DATA_SIZE} "
                f"bytes (you have {sys.getsizeof(data)})"
            )

        for key in self.data:
            if not key.startswith("_"):
                raise ValueError(
                    "data must contain properly namespaced data with keys starting with an underscore. "
                    f"{key=} does not start with an underscore."
                )

        # check for valid history nodes
        history = history or []  # initialize null fields
        if len(history) > MAX_HNODES:
            raise ValueError(f"A maximum of {MAX_HNODES} History nodes are supported, you have {len(history)}!")
        self.history = [HistoryNode.parse_history_node(h) for h in history]
        if not all(sys.getsizeof(h) < MAX_HNODE_SIZE for h in history):
            raise ValueError(f"One or more history nodes exceeds the maximum size limit of {MAX_HNODE_SIZE} bytes")

        self.created_at = created_at or datetime.datetime.utcnow()

    def as_dict(self):
        """Returns: MSONable dict."""
        dct = self.structure.as_dict()
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["about"] = {
            "authors": [a.as_dict() for a in self.authors],
            "projects": self.projects,
            "references": self.references,
            "remarks": self.remarks,
            "history": [h.as_dict() for h in self.history],
            "created_at": json.loads(json.dumps(self.created_at, cls=MontyEncoder)),
        }
        dct["about"].update(json.loads(json.dumps(self.data, cls=MontyEncoder)))
        return dct

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation.

        Returns:
            Class
        """
        a = d["about"]
        dec = MontyDecoder()

        created_at = dec.process_decoded(a.get("created_at"))
        data = {k: v for k, v in d["about"].items() if k.startswith("_")}
        data = dec.process_decoded(data)

        structure = Structure.from_dict(d) if "lattice" in d else Molecule.from_dict(d)
        return cls(
            structure,
            a["authors"],
            projects=a.get("projects"),
            references=a.get("references", ""),
            remarks=a.get("remarks"),
            data=data,
            history=a.get("history"),
            created_at=created_at,
        )

    @classmethod
    def from_structures(
        cls,
        structures: Sequence[Structure],
        authors: Sequence[dict[str, str]],
        projects=None,
        references="",
        remarks=None,
        data=None,
        histories=None,
        created_at=None,
    ):
        """A convenience method for getting a list of StructureNL objects by
        specifying structures and metadata separately. Some of the metadata
        is applied to all of the structures for ease of use.

        Args:
            structures: A list of Structure objects
            authors: *List* of {"name":'', "email":''} dicts,
                *list* of Strings as 'John Doe <johndoe@gmail.com>',
                or a single String with commas separating authors
            projects: List of Strings ['Project A', 'Project B']. This
                applies to all structures.
            references: A String in BibTeX format. Again, this applies to all
                structures.
            remarks: List of Strings ['Remark A', 'Remark B']
            data: A list of free form dict. Namespaced at the root level
                with an underscore, e.g. {"_materialsproject":<custom data>}
                . The length of data should be the same as the list of
                structures if not None.
            histories: List of list of dicts - [[{'name':'', 'url':'',
                'description':{}}], ...] The length of histories should be the
                same as the list of structures if not None.
            created_at: A datetime object
        """
        data = [{}] * len(structures) if data is None else data
        histories = [[]] * len(structures) if histories is None else histories

        snl_list = []
        for idx, struct in enumerate(structures):
            snl = StructureNL(
                struct,
                authors,
                projects=projects,
                references=references,
                remarks=remarks,
                data=data[idx],
                history=histories[idx],
                created_at=created_at,
            )
            snl_list.append(snl)

        return snl_list

    def __str__(self):
        return "\n".join(
            [
                f"{key}\n{getattr(self, key)}"
                for key in (
                    "structure",
                    "authors",
                    "projects",
                    "references",
                    "remarks",
                    "data",
                    "history",
                    "created_at",
                )
            ]
        )

    def __eq__(self, other: object) -> bool:
        needed_attrs = ("structure", "authors", "projects", "references", "remarks", "data", "history", "created_at")

        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        return all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)
