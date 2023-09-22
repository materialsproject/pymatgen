"""This module provides various representations of transformed structures. A
TransformedStructure is a structure that has been modified by undergoing a
series of transformations.
"""

from __future__ import annotations

import datetime
import json
import re
from typing import TYPE_CHECKING, Any
from warnings import warn

import numpy as np
from monty.json import MSONable, jsanitize

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, VaspInputSet
from pymatgen.util.provenance import StructureNL

if TYPE_CHECKING:
    from pymatgen.alchemy.filters import AbstractStructureFilter
    from pymatgen.transformations.transformation_abc import AbstractTransformation


class TransformedStructure(MSONable):
    """Container object for new structures that include history of
    transformations.

    Each transformed structure is made up of a sequence of structures with
    associated transformation history.
    """

    def __init__(
        self,
        structure: Structure,
        transformations: list[AbstractTransformation] | None = None,
        history: list[AbstractTransformation | dict[str, Any]] | None = None,
        other_parameters: dict[str, Any] | None = None,
    ) -> None:
        """Initializes a transformed structure from a structure.

        Args:
            structure (Structure): Input structure
            transformations (list[Transformation]): List of transformations to
                apply.
            history (list[Transformation]): Previous history.
            other_parameters (dict): Additional parameters to be added.
        """
        self.final_structure = structure
        self.history = history or []
        self.other_parameters = other_parameters or {}
        self._undone: list[tuple[AbstractTransformation | dict[str, Any], Structure]] = []

        transformations = transformations or []
        for trafo in transformations:
            self.append_transformation(trafo)

    def undo_last_change(self) -> None:
        """Undo the last change in the TransformedStructure.

        Raises:
            IndexError: If already at the oldest change.
        """
        if len(self.history) == 0:
            raise IndexError("No more changes to undo")
        if "input_structure" not in self.history[-1]:
            raise IndexError("Can't undo. Latest history has no input_structure")
        h = self.history.pop()
        self._undone.append((h, self.final_structure))
        struct = h["input_structure"]
        if isinstance(struct, dict):
            struct = Structure.from_dict(struct)
        self.final_structure = struct

    def redo_next_change(self) -> None:
        """Redo the last undone change in the TransformedStructure.

        Raises:
            IndexError: If already at the latest change.
        """
        if len(self._undone) == 0:
            raise IndexError("No more changes to redo")
        h, s = self._undone.pop()
        self.history.append(h)
        self.final_structure = s

    def __getattr__(self, name) -> Any:
        struct = object.__getattribute__(self, "final_structure")
        return getattr(struct, name)

    def __len__(self) -> int:
        return len(self.history)

    def append_transformation(
        self, transformation, return_alternatives: bool = False, clear_redo: bool = True
    ) -> list[TransformedStructure] | None:
        """Appends a transformation to the TransformedStructure.

        Args:
            transformation: Transformation to append
            return_alternatives: Whether to return alternative
                TransformedStructures for one-to-many transformations.
                return_alternatives can be a number, which stipulates the
                total number of structures to return.
            clear_redo: Boolean indicating whether to clear the redo list.
                By default, this is True, meaning any appends clears the
                history of undoing. However, when using append_transformation
                to do a redo, the redo list should not be cleared to allow
                multiple redos.
        """
        if clear_redo:
            self._undone = []

        if return_alternatives and transformation.is_one_to_many:
            ranked_list = transformation.apply_transformation(
                self.final_structure, return_ranked_list=return_alternatives
            )

            input_structure = self.final_structure.as_dict()
            alts = []
            for x in ranked_list[1:]:
                struct = x.pop("structure")
                actual_transformation = x.pop("transformation", transformation)
                h_dict = actual_transformation.as_dict()
                h_dict["input_structure"] = input_structure
                h_dict["output_parameters"] = x
                self.final_structure = struct
                d = self.as_dict()
                d["history"].append(h_dict)
                d["final_structure"] = struct.as_dict()
                alts.append(TransformedStructure.from_dict(d))

            x = ranked_list[0]
            struct = x.pop("structure")
            actual_transformation = x.pop("transformation", transformation)
            h_dict = actual_transformation.as_dict()
            h_dict["input_structure"] = self.final_structure.as_dict()
            h_dict["output_parameters"] = x
            self.history.append(h_dict)
            self.final_structure = struct
            return alts

        struct = transformation.apply_transformation(self.final_structure)
        h_dict = transformation.as_dict()
        h_dict["input_structure"] = self.final_structure.as_dict()
        h_dict["output_parameters"] = {}
        self.history.append(h_dict)
        self.final_structure = struct
        return None

    def append_filter(self, structure_filter: AbstractStructureFilter) -> None:
        """Adds a filter.

        Args:
            structure_filter (StructureFilter): A filter implementing the
                AbstractStructureFilter API. Tells transmuter what structures to retain.
        """
        h_dict = structure_filter.as_dict()
        h_dict["input_structure"] = self.final_structure.as_dict()
        self.history.append(h_dict)

    def extend_transformations(
        self, transformations: list[AbstractTransformation], return_alternatives: bool = False
    ) -> None:
        """Extends a sequence of transformations to the TransformedStructure.

        Args:
            transformations: Sequence of Transformations
            return_alternatives: Whether to return alternative
                TransformedStructures for one-to-many transformations.
                return_alternatives can be a number, which stipulates the
                total number of structures to return.
        """
        for t in transformations:
            self.append_transformation(t, return_alternatives=return_alternatives)

    def get_vasp_input(self, vasp_input_set: type[VaspInputSet] = MPRelaxSet, **kwargs) -> dict[str, Any]:
        """Returns VASP input as a dict of VASP objects.

        Args:
            vasp_input_set (pymatgen.io.vasp.sets.VaspInputSet): input set
                to create VASP input files from structures
            **kwargs: All keyword args supported by the VASP input set.
        """
        dct = vasp_input_set(self.final_structure, **kwargs).get_vasp_input()
        dct["transformations.json"] = json.dumps(self.as_dict())
        return dct

    def write_vasp_input(
        self,
        vasp_input_set: type[VaspInputSet] = MPRelaxSet,
        output_dir: str = ".",
        create_directory: bool = True,
        **kwargs,
    ) -> None:
        """Writes VASP input to an output_dir.

        Args:
            vasp_input_set: pymatgen.io.vasp.sets.VaspInputSet like object that creates vasp input files from
                structures.
            output_dir: Directory to output files
            create_directory: Create the directory if not present. Defaults to
                True.
            **kwargs: All keyword args supported by the VASP input set.
        """
        vasp_input_set(self.final_structure, **kwargs).write_input(output_dir, make_dir_if_not_present=create_directory)
        with open(f"{output_dir}/transformations.json", "w") as fp:
            json.dump(self.as_dict(), fp)

    def __str__(self) -> str:
        output = ["Current structure", "------------", str(self.final_structure), "\nHistory", "------------"]
        for h in self.history:
            h.pop("input_structure", None)
            output.append(str(h))
        output.append("\nOther parameters")
        output.append("------------")
        output.append(str(self.other_parameters))
        return "\n".join(output)

    def set_parameter(self, key: str, value: Any) -> None:
        """Set a parameter.

        :param key: The string key
        :param value: The value.
        """
        self.other_parameters[key] = value

    @property
    def was_modified(self) -> bool:
        """Boolean describing whether the last transformation on the structure
        made any alterations to it one example of when this would return false
        is in the case of performing a substitution transformation on the
        structure when the specie to replace isn't in the structure.
        """
        return self.final_structure != self.structures[-2]

    @property
    def structures(self) -> list[Structure]:
        """Copy of all structures in the TransformedStructure. A
        structure is stored after every single transformation.
        """
        h_structs = [Structure.from_dict(s["input_structure"]) for s in self.history if "input_structure" in s]
        return [*h_structs, self.final_structure]

    @classmethod
    @np.deprecate(message="Use from_cif_str instead")
    def from_cif_string(cls, *args, **kwargs):  # noqa: D102
        return cls.from_cif_str(*args, **kwargs)

    @classmethod
    def from_cif_str(
        cls,
        cif_string: str,
        transformations: list[AbstractTransformation] | None = None,
        primitive: bool = True,
        occupancy_tolerance: float = 1.0,
    ) -> TransformedStructure:
        """Generates TransformedStructure from a cif string.

        Args:
            cif_string (str): Input cif string. Should contain only one
                structure. For CIFs containing multiple structures, please use
                CifTransmuter.
            transformations (list[Transformation]): Sequence of transformations
                to be applied to the input structure.
            primitive (bool): Option to set if the primitive cell should be
                extracted. Defaults to True. However, there are certain
                instances where you might want to use a non-primitive cell,
                e.g., if you are trying to generate all possible orderings of
                partial removals or order a disordered structure.
            occupancy_tolerance (float): If total occupancy of a site is
                between 1 and occupancy_tolerance, the occupancies will be
                scaled down to 1.

        Returns:
            TransformedStructure
        """
        parser = CifParser.from_str(cif_string, occupancy_tolerance=occupancy_tolerance)
        raw_string = re.sub(r"'", '"', cif_string)
        cif_dict = parser.as_dict()
        cif_keys = list(cif_dict)
        struct = parser.get_structures(primitive)[0]
        partial_cif = cif_dict[cif_keys[0]]
        if "_database_code_ICSD" in partial_cif:
            source = partial_cif["_database_code_ICSD"] + "-ICSD"
        else:
            source = "uploaded cif"
        source_info = {
            "source": source,
            "datetime": str(datetime.datetime.now()),
            "original_file": raw_string,
            "cif_data": cif_dict[cif_keys[0]],
        }
        return cls(struct, transformations, history=[source_info])

    @classmethod
    @np.deprecate(message="Use from_poscar_str instead")
    def from_poscar_string(cls, *args, **kwargs):  # noqa: D102
        return cls.from_poscar_str(*args, **kwargs)

    @classmethod
    def from_poscar_str(
        cls, poscar_string: str, transformations: list[AbstractTransformation] | None = None
    ) -> TransformedStructure:
        """Generates TransformedStructure from a poscar string.

        Args:
            poscar_string (str): Input POSCAR string.
            transformations (list[Transformation]): Sequence of transformations
                to be applied to the input structure.
        """
        poscar = Poscar.from_str(poscar_string)
        if not poscar.true_names:
            raise ValueError(
                "Transformation can be created only from POSCAR strings with proper VASP5 element symbols."
            )
        raw_string = re.sub(r"'", '"', poscar_string)
        struct = poscar.structure
        source_info = {
            "source": "POSCAR",
            "datetime": str(datetime.datetime.now()),
            "original_file": raw_string,
        }
        return cls(struct, transformations, history=[source_info])

    def as_dict(self) -> dict[str, Any]:
        """Dict representation of the TransformedStructure."""
        dct = self.final_structure.as_dict()
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["history"] = jsanitize(self.history)
        dct["last_modified"] = str(datetime.datetime.utcnow())
        dct["other_parameters"] = jsanitize(self.other_parameters)
        return dct

    @classmethod
    def from_dict(cls, d) -> TransformedStructure:
        """Creates a TransformedStructure from a dict."""
        struct = Structure.from_dict(d)
        return cls(struct, history=d["history"], other_parameters=d.get("other_parameters"))

    def to_snl(self, authors, **kwargs) -> StructureNL:
        """Generate SNL from TransformedStructure.

        :param authors: List of authors
        :param **kwargs: All kwargs supported by StructureNL.

        Returns:
            StructureNL
        """
        if self.other_parameters:
            warn("Data in TransformedStructure.other_parameters discarded during type conversion to SNL")
        hist = []
        for h in self.history:
            snl_metadata = h.pop("_snl", {})
            hist.append(
                {
                    "name": snl_metadata.pop("name", "pymatgen"),
                    "url": snl_metadata.pop("url", "http://pypi.python.org/pypi/pymatgen"),
                    "description": h,
                }
            )

        return StructureNL(self.final_structure, authors, history=hist, **kwargs)

    @classmethod
    def from_snl(cls, snl: StructureNL) -> TransformedStructure:
        """Create TransformedStructure from SNL.

        Args:
            snl (StructureNL): Starting snl

        Returns:
            TransformedStructure
        """
        hist = []
        for h in snl.history:
            d = h.description
            d["_snl"] = {"url": h.url, "name": h.name}
            hist.append(d)
        return cls(snl.structure, history=hist)
