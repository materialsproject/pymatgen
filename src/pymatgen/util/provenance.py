"""Classes and methods related to the Structure Notation Language (SNL)."""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from typing import TYPE_CHECKING, NamedTuple

from monty.json import MontyDecoder, MontyEncoder

from pymatgen.core.structure import Molecule, Structure

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from typing_extensions import Self


__author__ = "Anubhav Jain, Shyue Ping Ong"
__credits__ = "Dan Gunter"


MAX_HNODE_SIZE: int = 64_000  # maximum size (bytes) of SNL HistoryNode
MAX_DATA_SIZE: int = 256_000  # maximum size (bytes) of SNL data field
MAX_HNODES: int = 100  # maximum number of HistoryNodes in SNL file
MAX_BIBTEX_CHARS: int = 20_000  # maximum number of characters for BibTeX reference


def is_valid_bibtex(reference: str) -> bool:
    """Validate that a reference is in proper BibTeX format.

    Args:
        reference (str): Reference in BibTeX format.

    Returns:
        bool: True if reference is valid BibTeX.
    """
    from bibtexparser.bparser import BibTexParser

    parser = BibTexParser()
    try:
        bib_database = parser.parse(reference)
        return bool(bib_database.entries)
    except Exception:
        return False


class HistoryNode(NamedTuple):
    """A HistoryNode represents a step in the chain of events that lead to a
    Structure. HistoryNodes leave 'breadcrumbs' so that you can trace back how
    a Structure was created. For example, a HistoryNode might represent pulling
    a Structure from an external database such as the ICSD or CSD. Or, it might
    represent the application of a code (e.g. pymatgen) to the Structure, with
    a custom description of how that code was applied (e.g. a site removal
    Transformation was applied).

    Attributes:
        name (str): The name of a code or resource that this Structure encountered in its history.
        url (str): The URL of that code/resource.
        description (str): A free-form description of how the code/resource is related to the Structure.
    """

    name: str
    url: str
    description: str

    def as_dict(self) -> dict[str, str]:
        """Get MSONable dict."""
        return {"name": self.name, "url": self.url, "description": self.description}

    @classmethod
    def from_dict(cls, dct: dict[str, str]) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            HistoryNode
        """
        return cls(dct["name"], dct["url"], dct["description"])

    @classmethod
    def parse_history_node(cls, h_node) -> Self:
        """Parse a History Node object from either a dict or a tuple.

        Args:
            h_node: A dict with name/url/description fields or a 3-element tuple.

        Returns:
            HistoryNode
        """
        if isinstance(h_node, dict):
            return cls.from_dict(h_node)

        if len(h_node) != 3:
            raise ValueError(f"Invalid History node, should be dict or (name, version, description) tuple: {h_node}")
        return cls(h_node[0], h_node[1], h_node[2])


class Author(NamedTuple):
    """An Author contains two fields: name and email. It is meant to represent
    the author of a Structure or the author of a code that was applied to a Structure.
    """

    name: str
    email: str

    def __str__(self) -> str:
        """String representation of an Author."""
        return f"{self.name} <{self.email}>"

    def as_dict(self) -> dict[str, str]:
        """Get MSONable dict."""
        return {"name": self.name, "email": self.email}

    @classmethod
    def from_dict(cls, dct: dict[str, str]) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Author
        """
        return cls(dct["name"], dct["email"])

    @classmethod
    def parse_author(
        cls,
        author: str | dict[str, str] | tuple[str, str],
    ) -> Self:
        """Parse an Author object from either a String, dict, or tuple.

        Args:
            author: in one of the three accepted formats:
                - A string formatted as "NAME <email@domain.com>".
                - A (name, email) tuple.
                - A dict with name and email as keys.

        Returns:
            An Author object.
        """
        if isinstance(author, str):
            # Regex looks for whitespace, (any name), whitespace, <, (email),
            # >, whitespace
            match = re.match(r"\s*(.*?)\s*<(.*?@.*?)>\s*", author)
            if not match or match.start() != 0 or match.end() != len(author):
                raise ValueError(f"Invalid author format! {author}")
            return cls(match.groups()[0], match.groups()[1])

        if isinstance(author, dict):
            return cls.from_dict(author)

        if len(author) != 2:
            raise ValueError(f"Invalid author, should be String or (name, email) tuple: {author}")

        return cls(author[0], author[1])


class StructureNL:
    """The Structure Notation Language (SNL, pronounced "snail") is a container for a pymatgen
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
        struct_or_mol: Structure | Molecule,
        authors: str | list[dict] | list[str],
        projects: list[str] | str | None = None,
        references: str = "",
        remarks: list[str] | str | None = None,
        data: dict | None = None,
        history: list[dict] | None = None,
        created_at: datetime | None = None,
    ) -> None:
        """
        Args:
            struct_or_mol: A pymatgen Structure/Molecule object
            authors (list | str):
                - List of {"name":"", "email":""} dicts.
                - List of strings as "John Doe <johndoe@gmail.com>".
                - A single String with commas separating authors.
            projects (list[str] | str): e.g. ["Project A", "Project B"]
            references (str): A string in BibTeX format
            remarks (list[str]): e.g. ["Remark A", "Remark B"]
            data (dict): A free form dict. Namespaced at the root level with an
                underscore, e.g. {"_materialsproject": <custom data>}
            history (list[dict]): e.g. [{"name":"", "url":"", "description":{}}]
            created_at (datetime): Creation date.
        """
        # Initialize root-level structure keys
        self.structure = struct_or_mol

        # Turn `authors` into list of `Author` objects
        _authors: list = authors.split(",") if isinstance(authors, str) else authors
        self.authors: list[Author] = [Author.parse_author(a) for a in _authors]

        # Turn `projects` into list of strings
        projects = projects or []
        self.projects: list[str] = [projects] if isinstance(projects, str) else projects

        # Check that references are valid BibTeX string
        if not isinstance(references, str):
            raise TypeError("Invalid format for SNL reference! Should be empty string or BibTeX string.")

        if references and not is_valid_bibtex(references):
            raise ValueError("Invalid format for SNL reference! Should be BibTeX string.")

        if len(references) > MAX_BIBTEX_CHARS:
            raise ValueError(
                f"The BibTeX string must be fewer than {MAX_BIBTEX_CHARS} chars, you have {len(references)}"
            )

        self.references = references

        # Turn `remarks` into list of strings
        remarks = remarks or []
        self.remarks: list[str] = [remarks] if isinstance(remarks, str) else remarks

        # Check `remarks` length limit
        for remark in self.remarks:
            if len(remark) > 140:
                raise ValueError(f"The remark exceeds the maximum size of 140 characters: {len(remark)}")

        # Check data length limit
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

        # Check for valid history nodes
        history = history or []  # initialize null fields
        if len(history) > MAX_HNODES:
            raise ValueError(f"A maximum of {MAX_HNODES} History nodes are supported, you have {len(history)}!")
        self.history = [HistoryNode.parse_history_node(h) for h in history]
        if not all(sys.getsizeof(h) < MAX_HNODE_SIZE for h in history):
            raise ValueError(f"One or more history nodes exceeds the maximum size limit of {MAX_HNODE_SIZE} bytes")

        self.created_at = created_at or f"{datetime.now(tz=timezone.utc):%Y-%m-%d %H:%M:%S.%f%z}"

    def __str__(self) -> str:
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
        """Check for equality between two StructureNL objects."""
        needed_attrs = (
            "structure",
            "authors",
            "projects",
            "references",
            "remarks",
            "data",
            "history",
            "created_at",
        )

        if not all(hasattr(other, attr) for attr in needed_attrs):
            return NotImplemented

        return all(getattr(self, attr) == getattr(other, attr) for attr in needed_attrs)

    def as_dict(self) -> dict[str, Any]:
        """Get MSONable dict."""
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
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Class
        """
        about = dct["about"]

        created_at = MontyDecoder().process_decoded(about.get("created_at"))
        data = {k: v for k, v in dct["about"].items() if k.startswith("_")}
        data = MontyDecoder().process_decoded(data)

        structure = Structure.from_dict(dct) if "lattice" in dct else Molecule.from_dict(dct)
        return cls(
            structure,
            about["authors"],
            projects=about.get("projects"),
            references=about.get("references", ""),
            remarks=about.get("remarks"),
            data=data,
            history=about.get("history"),
            created_at=created_at,
        )

    @classmethod
    def from_structures(
        cls,
        structures: Sequence[Structure],
        authors,
        projects: list[str] | None = None,
        references: str = "",
        remarks: list[str] | None = None,
        data: list[dict] | None = None,
        histories: list[list[dict]] | None = None,
        created_at: datetime | None = None,
    ) -> list[Self]:
        """A convenience method for getting a list of StructureNL objects by
        specifying structures and metadata separately. Some of the metadata
        is applied to all of the structures for ease of use.

        Args:
            structures (Sequence[Structure]): Structure objects
            authors: *List* of {"name":'', "email":''} dicts,
                *list* of Strings as 'John Doe <johndoe@gmail.com>',
                or a single String with commas separating authors
            projects (list[str]): e.g. ["Project A", "Project B"]. This
                applies to all structures.
            references (str): A String in BibTeX format. This applies to all
                structures.
            remarks (list[str]): e.g. ["Remark A", "Remark B"].
            data (list[dict]): Free form dicts. Namespaced at the root level
                with an underscore, e.g. {"_materialsproject":<custom data>}
                . The length of data should be the same as the list of
                structures if not None.
            histories(list[list[dict]]): e.g. [[{"name":"", "url":"",
                "description":{}}], ...] The length of histories should be the
                same as the list of structures if not None.
            created_at (datetime): Creation date.
        """
        data = [{}] * len(structures) if data is None else data
        histories = [[]] * len(structures) if histories is None else histories

        snl_list = []
        for idx, struct in enumerate(structures):
            snl = cls(
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
