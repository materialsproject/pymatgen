# coding: utf-8

from __future__ import unicode_literals

"""
This module implements various transmuter classes.
Transmuters are essentially classes that generate TransformedStructures from
various data sources. They enable the high-throughput generation of new
structures and input files.

It also includes the helper function, batch_write_vasp_input to generate an
entire directory of vasp input files for running.
"""

from six.moves import filter, map

__author__ = "Shyue Ping Ong, Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 4, 2012"

import os
import re
import warnings

from multiprocessing import Pool
from pymatgen.alchemy.materials import TransformedStructure


class StandardTransmuter(object):
    """
    An example of a Transmuter object, which performs a sequence of
    transformations on many structures to generate TransformedStructures.

    .. attribute: transformed_structures

        List of all transformed structures.
    """

    def __init__(self, transformed_structures, transformations=None,
                 extend_collection=0, ncores=None):
        """
        Initializes a transmuter from an initial list of
        :class:`pymatgen.alchemy.materials.TransformedStructure`.

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
                self.append_transformation(trans,
                                           extend_collection=extend_collection)

    def get_transformed_structures(self):
        """
        Returns all TransformedStructures.

        .. deprecated:: v2.1.0

            Use transformed_structures attribute instead. Will be removed in
            next version.
        """
        warnings.warn("Use transformed_structures attribute instead.",
                      DeprecationWarning)
        return self.transformed_structures

    def __getitem__(self, index):
        return self.transformed_structures[index]

    def __getattr__(self, name):
        return [getattr(x, name) for x in self.transformed_structures]

    def undo_last_change(self):
        """
        Undo the last transformation in the TransformedStructure.

        Raises:
            IndexError if already at the oldest change.
        """
        for x in self.transformed_structures:
            x.undo_last_change()

    def redo_next_change(self):
        """
        Redo the last undone transformation in the TransformedStructure.

        Raises:
            IndexError if already at the latest change.
        """
        for x in self.transformed_structures:
            x.redo_next_change()

    def __len__(self):
        return len(self.transformed_structures)

    def append_transformation(self, transformation, extend_collection=False,
                              clear_redo=True):
        """
        Appends a transformation to all TransformedStructures.

        Args:
            transformation: Transformation to append
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations. extend_collection can be a
                number, which determines the maximum branching for each
                transformation.
            clear_redo (bool): Whether to clear the redo list. By default,
                this is True, meaning any appends clears the history of
                undoing. However, when using append_transformation to do a
                redo, the redo list should not be cleared to allow multiple
                redos.

        Returns:
            List of booleans corresponding to initial transformed structures
            each boolean describes whether the transformation altered the
            structure
        """
        if self.ncores and transformation.use_multiprocessing:
            p = Pool(self.ncores)
            #need to condense arguments into single tuple to use map
            z = map(
                lambda x: (x, transformation, extend_collection, clear_redo),
                self.transformed_structures)
            new_tstructs = p.map(_apply_transformation, z, 1)
            self.transformed_structures = []
            for ts in new_tstructs:
                self.transformed_structures.extend(ts)
        else:
            new_structures = []
            for x in self.transformed_structures:
                new = x.append_transformation(transformation,
                                              extend_collection,
                                              clear_redo=clear_redo)
                if new is not None:
                    new_structures.extend(new)
            self.transformed_structures.extend(new_structures)

    def extend_transformations(self, transformations):
        """
        Extends a sequence of transformations to the TransformedStructure.

        Args:
            transformations: Sequence of Transformations
        """
        for t in transformations:
            self.append_transformation(t)

    def apply_filter(self, structure_filter):
        """
        Applies a structure_filter to the list of TransformedStructures
        in the transmuter.

        Args:
            structure_filter: StructureFilter to apply.
        """

        def test_transformed_structure(ts):
            return structure_filter.test(ts.final_structure)

        self.transformed_structures = list(filter(test_transformed_structure,
                                                  self.transformed_structures))
        for ts in self.transformed_structures:
            ts.append_filter(structure_filter)

    def write_vasp_input(self, vasp_input_set, output_dir,
                         create_directory=True, subfolder=None,
                         include_cif=False):
        """
        Batch write vasp input for a sequence of transformed structures to
        output_dir, following the format output_dir/{formula}_{number}.

        Args:
            vasp_input_set: pymatgen.io.vaspio_set.VaspInputSet to create
                vasp input files from structures
            output_dir: Directory to output files
            create_directory (bool): Create the directory if not present.
                Defaults to True.
            subfolder: Callable to create subdirectory name from
                transformed_structure. e.g.,
                lambda x: x.other_parameters["tags"][0] to use the first tag.
            include_cif (bool): Whether to output a CIF as well. CIF files
                are generally better supported in visualization programs.
        """
        batch_write_vasp_input(self.transformed_structures, vasp_input_set,
                               output_dir, create_directory, subfolder,
                               include_cif)

    def set_parameter(self, key, value):
        """
        Add parameters to the transmuter. Additional parameters are stored in
        the as_dict() output.

        Args:
            key: The key for the parameter.
            value: The value for the parameter.
        """
        for x in self.transformed_structures:
            x.other_parameters[key] = value

    def add_tags(self, tags):
        """
        Add tags for the structures generated by the transmuter.

        Args:
            tags: A sequence of tags. Note that this should be a sequence of
                strings, e.g., ["My awesome structures", "Project X"].
        """
        self.set_parameter("tags", tags)

    def __str__(self):
        output = ["Current structures", "------------"]
        for x in self.transformed_structures:
            output.append(str(x.final_structure))
        return "\n".join(output)

    def append_transformed_structures(self, tstructs_or_transmuter):
        """
        Method is overloaded to accept either a list of transformed structures
        or transmuter, it which case it appends the second transmuter"s
        structures.

        Args:
            tstructs_or_transmuter: A list of transformed structures or a
                transmuter.
        """
        if isinstance(tstructs_or_transmuter, self.__class__):
            self.transformed_structures.extend(tstructs_or_transmuter
                                               .transformed_structures)
        else:
            for ts in tstructs_or_transmuter:
                assert isinstance(ts, TransformedStructure)
            self.transformed_structures.extend(tstructs_or_transmuter)

    @staticmethod
    def from_structures(structures, transformations=None, extend_collection=0):
        """
        Alternative constructor from structures rather than
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
        tstruct = [TransformedStructure(s, []) for s in structures]
        return StandardTransmuter(tstruct, transformations, extend_collection)


class CifTransmuter(StandardTransmuter):
    """
    Generates a Transmuter from a cif string, possibly containing multiple
    structures.
    """

    def __init__(self, cif_string, transformations=None, primitive=True,
                 extend_collection=False):
        """
        Generates a Transmuter from a cif string, possibly
        containing multiple structures.

        Args:
            cif_string: A string containing a cif or a series of cifs
            transformations: New transformations to be applied to all
                structures
            primitive: Whether to generate the primitive cell from the cif.
            extend_collection: Whether to use more than one output structure
                from one-to-many transformations. extend_collection can be a
                number, which determines the maximum branching for each
                transformation.
        """
        transformed_structures = []
        lines = cif_string.split("\n")
        structure_data = []
        read_data = False
        for line in lines:
            if re.match("^\s*data", line):
                structure_data.append([])
                read_data = True
            if read_data:
                structure_data[-1].append(line)
        for data in structure_data:
            tstruct = TransformedStructure.from_cif_string("\n".join(data), [],
                                                           primitive)
            transformed_structures.append(tstruct)
        StandardTransmuter.__init__(self, transformed_structures,
                                    transformations, extend_collection)

    @staticmethod
    def from_filenames(filenames, transformations=None, primitive=True,
                       extend_collection=False):
        """
        Generates a TransformedStructureCollection from a cif, possibly
        containing multiple structures.

        Args:
            filenames: List of strings of the cif files
            transformations: New transformations to be applied to all
                structures
            primitive: Same meaning as in __init__.
            extend_collection: Same meaning as in __init__.
        """

        allcifs = []
        for fname in filenames:
            with open(fname, "r") as f:
                allcifs.append(f.read())
        return CifTransmuter("\n".join(allcifs), transformations,
                             primitive=primitive,
                             extend_collection=extend_collection)


class PoscarTransmuter(StandardTransmuter):
    """
    Generates a transmuter from a sequence of POSCARs.

    Args:
        poscar_string: List of POSCAR strings
        transformations: New transformations to be applied to all
            structures.
        extend_collection: Whether to use more than one output structure
            from one-to-many transformations.
    """

    def __init__(self, poscar_string, transformations=None,
                 extend_collection=False):
        tstruct = TransformedStructure.from_poscar_string(poscar_string, [])
        StandardTransmuter.__init__(self, [tstruct], transformations,
                                    extend_collection=extend_collection)

    @staticmethod
    def from_filenames(poscar_filenames, transformations=None,
                       extend_collection=False):
        """
        Convenient constructor to generates a POSCAR transmuter from a list of
        POSCAR filenames.

        Args:
            poscar_filenames: List of POSCAR filenames
            transformations: New transformations to be applied to all
                structures.
            extend_collection:
                Same meaning as in __init__.
        """
        tstructs = []
        for filename in poscar_filenames:
            with open(filename, "r") as f:
                tstructs.append(TransformedStructure
                                .from_poscar_string(f.read(), []))
        return StandardTransmuter(tstructs, transformations,
                                  extend_collection=extend_collection)


def batch_write_vasp_input(transformed_structures, vasp_input_set, output_dir,
                           create_directory=True, subfolder=None,
                           include_cif=False):
    """
    Batch write vasp input for a sequence of transformed structures to
    output_dir, following the format output_dir/{group}/{formula}_{number}.

    Args:
        transformed_structures: Sequence of TransformedStructures.
        vasp_input_set: pymatgen.io.vaspio_set.VaspInputSet to creates
            vasp input files from structures.
        output_dir: Directory to output files
        create_directory (bool): Create the directory if not present.
            Defaults to True.
        subfolder: Function to create subdirectory name from
            transformed_structure.
            e.g., lambda x: x.other_parameters["tags"][0] to use the first
            tag.
        include_cif (bool): Boolean indication whether to output a CIF as
            well. CIF files are generally better supported in visualization
            programs.
    """
    for i, s in enumerate(transformed_structures):
        formula = re.sub("\s+", "", s.final_structure.formula)
        if subfolder is not None:
            subdir = subfolder(s)
            dirname = os.path.join(output_dir, subdir,
                                   "{}_{}".format(formula, i))
        else:
            dirname = os.path.join(output_dir, "{}_{}".format(formula, i))
        s.write_vasp_input(vasp_input_set, dirname,
                           create_directory=create_directory)
        if include_cif:
            from pymatgen.io.cifio import CifWriter

            writer = CifWriter(s.final_structure)
            writer.write_file(os.path.join(dirname, "{}.cif".format(formula)))


def _apply_transformation(inputs):
    """
    Helper method for multiprocessing of apply_transformation. Must not be
    in the class so that it can be pickled.

    Args:
        inputs: Tuple containing the transformed structure, the transformation
            to be applied, a boolean indicating whether to extend the
            collection, and a boolean indicating whether to clear the redo

    Returns:
        List of output structures (the modified initial structure, plus
        any new structures created by a one-to-many transformation)
    """
    ts, transformation, extend_collection, clear_redo = inputs
    new = ts.append_transformation(transformation, extend_collection,
                                   clear_redo=clear_redo)
    o = [ts]
    if new:
        o.extend(new)
    return o
