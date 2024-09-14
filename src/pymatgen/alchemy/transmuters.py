"""This module implements various transmuter classes.
Transmuters are essentially classes that generate TransformedStructures from
various data sources. They enable the high-throughput generation of new
structures and input files.

It also includes the helper function, batch_write_vasp_input to generate an
entire directory of vasp input files for running.
"""

from __future__ import annotations

import os
import re
from multiprocessing import Pool
from typing import TYPE_CHECKING

from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.io.vasp.sets import MPRelaxSet, VaspInputSet

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from typing_extensions import Self

    from pymatgen.alchemy.filters import AbstractStructureFilter

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 4, 2012"


class StandardTransmuter:
    """An example of a Transmuter object, which performs a sequence of
    transformations on many structures to generate TransformedStructures.

    Attributes:
        transformed_structures (list[Structure]): All transformed structures.
    """

    def __init__(
        self,
        transformed_structures: list[TransformedStructure],
        transformations=None,
        extend_collection: int = 0,
        ncores: int | None = None,
    ) -> None:
        """Initialize a transmuter from an initial list of
        pymatgen.alchemy.materials.TransformedStructure.

        Args:
            transformed_structures ([TransformedStructure]): Input transformed
                structures
            transformations ([Transformations]): New transformations to be
                applied to all structures.
            extend_collection (int): Whether to use more than one output
                structure from one-to-many transformations. extend_collection
                can be an int, which determines the maximum branching for each
                transformation.
            ncores (int): Number of cores to use for applying transformations.
                Uses multiprocessing.Pool. Default is None, which implies
                serial.
        """
        self.transformed_structures = transformed_structures
        self.ncores = ncores
        if transformations is not None:
            for trans in transformations:
                self.append_transformation(trans, extend_collection=extend_collection)

    def __getitem__(self, index):
        return self.transformed_structures[index]

    def __getattr__(self, name):
        return [getattr(ts, name) for ts in self.transformed_structures]

    def __len__(self):
        return len(self.transformed_structures)

    def __str__(self):
        output = ["Current structures", "------------"]
        for ts in self.transformed_structures:
            output.append(str(ts.final_structure))
        return "\n".join(output)

    def undo_last_change(self) -> None:
        """Undo the last transformation in the TransformedStructure.

        Raises:
            IndexError if already at the oldest change.
        """
        for ts in self.transformed_structures:
            ts.undo_last_change()

    def redo_next_change(self) -> None:
        """Redo the last undone transformation in the TransformedStructure.

        Raises:
            IndexError if already at the latest change.
        """
        for ts in self.transformed_structures:
            ts.redo_next_change()

    def append_transformation(self, transformation, extend_collection=False, clear_redo=True) -> list[bool]:
        """Append a transformation to all TransformedStructures.

        Args:
            transformation: Transformation to append
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations. extend_collection can be a
                number, which determines the maximum branching for each transformation.
            clear_redo (bool): Whether to clear the redo list. By default,
                this is True, meaning any appends clears the history of
                undoing. However, when using append_transformation to do a
                redo, the redo list should not be cleared to allow multiple redos.

        Returns:
            list[bool]: Each list item is True if the transformation altered the structure
                with the corresponding index.
        """
        if self.ncores and transformation.use_multiprocessing:
            with Pool(self.ncores) as pool:
                # need to condense arguments into single tuple to use map
                z = ((ts, transformation, extend_collection, clear_redo) for ts in self.transformed_structures)
                trafo_new_structs = pool.map(_apply_transformation, z, 1)
                self.transformed_structures = []
                for ts in trafo_new_structs:
                    self.transformed_structures.extend(ts)
        else:
            new_structures = []
            for ts in self.transformed_structures:
                new = ts.append_transformation(transformation, extend_collection, clear_redo=clear_redo)
                if new is not None:
                    new_structures += new
            self.transformed_structures += new_structures

        # len(ts) > 1 checks if the structure has history
        return [len(ts) > 1 for ts in self.transformed_structures]

    def extend_transformations(self, transformations):
        """Extend a sequence of transformations to the TransformedStructure.

        Args:
            transformations: Sequence of Transformations
        """
        for trafo in transformations:
            self.append_transformation(trafo)

    def apply_filter(self, structure_filter: AbstractStructureFilter):
        """Apply a structure_filter to the list of TransformedStructures
        in the transmuter.

        Args:
            structure_filter: StructureFilter to apply.
        """
        self.transformed_structures = list(
            filter(lambda ts: structure_filter.test(ts.final_structure), self.transformed_structures)
        )
        for ts in self.transformed_structures:
            ts.append_filter(structure_filter)

    def write_vasp_input(self, **kwargs):
        """Batch write vasp input for a sequence of transformed structures to
        output_dir, following the format output_dir/{formula}_{number}.

        Args:
            kwargs: All kwargs supported by batch_write_vasp_input.
        """
        batch_write_vasp_input(self.transformed_structures, **kwargs)

    def set_parameter(self, key, value):
        """Add parameters to the transmuter. Additional parameters are stored in
        the as_dict() output.

        Args:
            key: The key for the parameter.
            value: The value for the parameter.
        """
        for struct in self.transformed_structures:
            struct.other_parameters[key] = value

    def add_tags(self, tags):
        """Add tags for the structures generated by the transmuter.

        Args:
            tags: A sequence of tags. Note that this should be a sequence of
                strings, e.g. ["My awesome structures", "Project X"].
        """
        self.set_parameter("tags", tags)

    def append_transformed_structures(self, trafo_structs_or_transmuter):
        """Overloaded to accept either a list of transformed structures
        or transmuter, it which case it appends the second transmuter's structures.

        Args:
            trafo_structs_or_transmuter: A list of transformed structures or a transmuter.
        """
        if not isinstance(trafo_structs_or_transmuter, self.__class__) and not all(
            isinstance(ts, TransformedStructure) for ts in trafo_structs_or_transmuter
        ):
            raise TypeError("Some transformed structure has incorrect type.")

        self.transformed_structures += trafo_structs_or_transmuter

    @classmethod
    def from_structures(cls, structures, transformations=None, extend_collection=0) -> Self:
        """Alternative constructor from structures rather than
        TransformedStructures.

        Args:
            structures: Sequence of structures
            transformations: New transformations to be applied to all
                structures
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations. extend_collection can be a
                number, which determines the maximum branching for each
                transformation.

        Returns:
            StandardTransmuter
        """
        t_struct = [TransformedStructure(s, []) for s in structures]
        return cls(t_struct, transformations, extend_collection)


class CifTransmuter(StandardTransmuter):
    """Generate a Transmuter from a CIF string, possibly containing multiple structures."""

    def __init__(self, cif_string, transformations=None, primitive=True, extend_collection=False):
        """Generate a Transmuter from a CIF string, possibly
        containing multiple structures.

        Args:
            cif_string: A string containing a CIF or a series of CIFs
            transformations: New transformations to be applied to all
                structures
            primitive: Whether to generate the primitive cell from the CIF.
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations. extend_collection can be a
                number, which determines the maximum branching for each
                transformation.
        """
        transformed_structures = []
        lines = cif_string.split("\n")
        structure_data: list = []
        read_data = False
        for line in lines:
            if re.match(r"^\s*data", line):
                structure_data.append([])
                read_data = True
            if read_data:
                structure_data[-1].append(line)
        for data in structure_data:
            t_struct = TransformedStructure.from_cif_str("\n".join(data), [], primitive)
            transformed_structures.append(t_struct)
        super().__init__(transformed_structures, transformations, extend_collection)

    @classmethod
    def from_filenames(cls, filenames, transformations=None, primitive=True, extend_collection=False) -> Self:
        """Generate a TransformedStructureCollection from a cif, possibly
        containing multiple structures.

        Args:
            filenames (list[str]): The CIF file paths.
            transformations: New transformations to be applied to all
                structures
            primitive: Same meaning as in __init__.
            extend_collection: Same meaning as in __init__.
        """
        cif_files = []
        for filename in filenames:
            with open(filename, encoding="utf-8") as file:
                cif_files.append(file.read())
        return cls(
            "\n".join(cif_files),
            transformations,
            primitive=primitive,
            extend_collection=extend_collection,
        )


class PoscarTransmuter(StandardTransmuter):
    """Generate a transmuter from a sequence of POSCARs."""

    def __init__(self, poscar_string, transformations=None, extend_collection=False):
        """
        Args:
            poscar_string (list[str]): POSCAR strings.
            transformations: New transformations to be applied to all
                structures.
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations.
        """
        t_struct = TransformedStructure.from_poscar_str(poscar_string, [])
        super().__init__([t_struct], transformations, extend_collection=extend_collection)

    @classmethod
    def from_filenames(cls, poscar_filenames, transformations=None, extend_collection=False) -> StandardTransmuter:
        """Convenient constructor to generates a POSCAR transmuter from a list of
        POSCAR filenames.

        Args:
            poscar_filenames (list[str]): The POSCAR file paths.
            transformations: New transformations to be applied to all
                structures.
            extend_collection:
                Same meaning as in __init__.
        """
        trafo_structs = []
        for filename in poscar_filenames:
            with open(filename, encoding="utf-8") as file:
                trafo_structs.append(TransformedStructure.from_poscar_str(file.read(), []))
        return StandardTransmuter(trafo_structs, transformations, extend_collection=extend_collection)


def batch_write_vasp_input(
    transformed_structures: Sequence[TransformedStructure],
    vasp_input_set: type[VaspInputSet] = MPRelaxSet,
    output_dir: str = ".",
    create_directory: bool = True,
    subfolder: Callable[[TransformedStructure], str] | None = None,
    include_cif: bool = False,
    **kwargs,
):
    """Batch write vasp input for a sequence of transformed structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        transformed_structures: Sequence of TransformedStructures.
        vasp_input_set: pymatgen.io.vasp.sets.VaspInputSet to creates
            vasp input files from structures.
        output_dir: Directory to output files
        create_directory (bool): Create the directory if not present.
            Defaults to True.
        subfolder: Function to create subdirectory name from transformed_structure.
            E.g. lambda x: x.other_parameters["tags"][0] to use the first tag.
        include_cif (bool): Pass True to output a CIF as well. CIF files are generally
            better supported in visualization programs.
        **kwargs: Any kwargs supported by vasp_input_set.
    """
    for idx, struct in enumerate(transformed_structures):
        formula = re.sub(r"\s+", "", struct.final_structure.formula)
        if subfolder is not None:
            subdir = subfolder(struct)
            dirname = f"{output_dir}/{subdir}/{formula}_{idx}"
        else:
            dirname = f"{output_dir}/{formula}_{idx}"
        struct.write_vasp_input(vasp_input_set, dirname, create_directory=create_directory, **kwargs)
        if include_cif:
            from pymatgen.io.cif import CifWriter

            writer = CifWriter(struct.final_structure)
            writer.write_file(os.path.join(dirname, f"{formula}.cif"))


def _apply_transformation(inputs):
    """Helper method for multiprocessing of apply_transformation. Must not be
    in the class so that it can be pickled.

    Args:
        inputs: Tuple containing the transformed structure, the transformation
            to be applied, a boolean indicating whether to extend the
            collection, and a boolean indicating whether to clear the redo

    Returns:
        list[Structure]: the modified initial structure, plus
            any new structures created by a one-to-many transformation
    """
    ts, transformation, extend_collection, clear_redo = inputs
    new = ts.append_transformation(transformation, extend_collection, clear_redo=clear_redo)
    out = [ts]
    if new:
        out += new
    return out
