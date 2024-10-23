"""This module provides various representations of transformed structures. A
TransformedStructure is a structure that has been modified by undergoing a
series of transformations.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from typing import TYPE_CHECKING
from warnings import warn

from monty.json import MSONable, jsanitize

from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, VaspInputSet
from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.util.provenance import StructureNL

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from typing_extensions import Self

    from pymatgen.alchemy.filters import AbstractStructureFilter


class TransformedStructure(MSONable):
    """Container for new structures that include history of transformations.

    Each transformed structure is made up of a sequence of structures with
    associated transformation history.
    """

    def __init__(
        self,
        structure: Structure,
        transformations: (AbstractTransformation | Sequence[AbstractTransformation] | None) = None,
        history: list[AbstractTransformation | dict[str, Any]] | None = None,
        other_parameters: dict[str, Any] | None = None,
    ) -> None:
        """Initialize a transformed structure from a structure.

        Args:
            structure (Structure): Input structure
            transformations (list[Transformation]): Transformations to apply.
            history (list[Transformation]): Previous history.
            other_parameters (dict): Additional parameters to be added.
        """
        self.final_structure = structure
        self.history = history or []
        self.other_parameters = other_parameters or {}
        self._undone: list[tuple[AbstractTransformation | dict[str, Any], Structure]] = []

        if isinstance(transformations, AbstractTransformation):
            transformations = [transformations]
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
        hist, struct = self._undone.pop()
        self.history.append(hist)
        self.final_structure = struct

    def __getattr__(self, name: str) -> Any:
        # Don't use getattr(self.final_structure, name) here to avoid infinite recursion if name = "final_structure"
        struct = self.__getattribute__("final_structure")
        return getattr(struct, name)

    def __len__(self) -> int:
        return len(self.history)

    def append_transformation(
        self, transformation, return_alternatives: bool = False, clear_redo: bool = True
    ) -> list[TransformedStructure] | None:
        """Append a transformation to the TransformedStructure.

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
                dct = self.as_dict()
                dct["history"].append(h_dict)
                dct["final_structure"] = struct.as_dict()
                alts.append(TransformedStructure.from_dict(dct))

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
        """Add a filter.

        Args:
            structure_filter (StructureFilter): A filter implementing the
                AbstractStructureFilter API. Tells transmuter what structures to retain.
        """
        h_dict = structure_filter.as_dict()
        h_dict["input_structure"] = self.final_structure.as_dict()
        self.history.append(h_dict)

    def extend_transformations(
        self,
        transformations: list[AbstractTransformation],
        return_alternatives: bool = False,
    ) -> None:
        """Extend a sequence of transformations to the TransformedStructure.

        Args:
            transformations: Sequence of Transformations
            return_alternatives: Whether to return alternative
                TransformedStructures for one-to-many transformations.
                return_alternatives can be a number, which stipulates the
                total number of structures to return.
        """
        for trafo in transformations:
            self.append_transformation(trafo, return_alternatives=return_alternatives)

    def get_vasp_input(self, vasp_input_set: type[VaspInputSet] = MPRelaxSet, **kwargs) -> dict[str, Any]:
        """Get VASP input as a dict of VASP objects.

        Args:
            vasp_input_set (VaspInputSet): input set
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
        """Write VASP input to an output_dir.

        Args:
            vasp_input_set: pymatgen.io.vasp.sets.VaspInputSet like object that creates vasp input files from
                structures.
            output_dir: Directory to output files
            create_directory: Create the directory if not present. Defaults to
                True.
            **kwargs: All keyword args supported by the VASP input set.
        """
        vasp_input_set(self.final_structure, **kwargs).write_input(output_dir, make_dir_if_not_present=create_directory)
        with open(f"{output_dir}/transformations.json", mode="w", encoding="utf-8") as file:
            json.dump(self.as_dict(), file)

    def __str__(self) -> str:
        output = [
            "Current structure",
            "------------",
            str(self.final_structure),
            "\nHistory",
            "------------",
        ]
        for hist in self.history:
            hist.pop("input_structure", None)
            output.append(str(hist))
        output += ("\nOther parameters", "------------", str(self.other_parameters))
        return "\n".join(output)

    def set_parameter(self, key: str, value: Any) -> TransformedStructure:
        """Set a parameter.

        Args:
            key (str): The string key.
            value (Any): The value.

        Returns:
            TransformedStructure
        """
        self.other_parameters[key] = value

        return self

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
    def from_cif_str(
        cls,
        cif_string: str,
        transformations: list[AbstractTransformation] | None = None,
        primitive: bool = True,
        occupancy_tolerance: float = 1.0,
    ) -> Self:
        """Generate TransformedStructure from a CIF string.

        Args:
            cif_string (str): Input CIF string. Should contain only one
                structure. For CIFs containing multiple structures, please use
                CifTransmuter.
            transformations (list[Transformation]): Sequence of transformations
                to be applied to the input structure.
            primitive (bool): Option to set if the primitive cell should be
                extracted. Defaults to True. However, there are certain
                instances where you might want to use a non-primitive cell,
                e.g. if you are trying to generate all possible orderings of
                partial removals or order a disordered structure. Defaults to True.
            occupancy_tolerance (float): If total occupancy of a site is
                between 1 and occupancy_tolerance, the occupancies will be
                scaled down to 1.

        Returns:
            TransformedStructure
        """
        parser = CifParser.from_str(cif_string, occupancy_tolerance=occupancy_tolerance)
        raw_str = re.sub(r"'", '"', cif_string)
        cif_dict = parser.as_dict()
        cif_keys = list(cif_dict)
        struct = parser.parse_structures(primitive=primitive)[0]
        partial_cif = cif_dict[cif_keys[0]]
        if "_database_code_ICSD" in partial_cif:
            source = partial_cif["_database_code_ICSD"] + "-ICSD"
        else:
            source = "uploaded cif"
        source_info = {
            "source": source,
            "datetime": str(datetime.now(tz=timezone.utc)),
            "original_file": raw_str,
            "cif_data": cif_dict[cif_keys[0]],
        }
        return cls(struct, transformations, history=[source_info])

    @classmethod
    def from_poscar_str(
        cls,
        poscar_string: str,
        transformations: list[AbstractTransformation] | None = None,
    ) -> Self:
        """Generate TransformedStructure from a poscar string.

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
        raw_str = re.sub(r"'", '"', poscar_string)
        struct = poscar.structure
        source_info = {
            "source": "POSCAR",
            "datetime": str(datetime.now(tz=timezone.utc)),
            "original_file": raw_str,
        }
        return cls(struct, transformations, history=[source_info])

    def as_dict(self) -> dict[str, Any]:
        """Dict representation of the TransformedStructure."""
        dct = self.final_structure.as_dict()
        dct["@module"] = type(self).__module__
        dct["@class"] = type(self).__name__
        dct["history"] = jsanitize(self.history)
        dct["last_modified"] = str(datetime.now(timezone.utc))
        dct["other_parameters"] = jsanitize(self.other_parameters)
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create a TransformedStructure from a dict."""
        struct = Structure.from_dict(dct)
        return cls(struct, history=dct["history"], other_parameters=dct.get("other_parameters"))

    def to_snl(self, authors: list[str], **kwargs) -> StructureNL:
        """Generate a StructureNL from TransformedStructure.

        Args:
            authors (List[str]): Authors contributing to the generated StructureNL.
            **kwargs (Any): All kwargs supported by StructureNL.

        Returns:
            StructureNL: The generated StructureNL object.
        """
        if self.other_parameters:
            warn("Data in TransformedStructure.other_parameters discarded during type conversion to SNL")
        history = []
        for hist in self.history:
            snl_metadata = hist.pop("_snl", {})
            history += [
                {
                    "name": snl_metadata.pop("name", "pymatgen"),
                    "url": snl_metadata.pop("url", "http://pypi.python.org/pypi/pymatgen"),
                    "description": hist,
                }
            ]

        return StructureNL(self.final_structure, authors, history=history, **kwargs)

    @classmethod
    def from_snl(cls, snl: StructureNL) -> Self:
        """Create TransformedStructure from SNL.

        Args:
            snl (StructureNL): Starting snl

        Returns:
            TransformedStructure
        """
        history: list[dict] = []
        for hist in snl.history:
            dct = hist.description
            dct["_snl"] = {"url": hist.url, "name": hist.name}
            history.append(dct)
        return cls(snl.structure, history=history)
