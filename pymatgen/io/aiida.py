# coding: utf-8
# Copyright (c) AiiDA Development Team.
# Distributed under the terms of the MIT License.

from __future__ import absolute_import

"""
This module provides conversion between AiiDA StructureData object and
pymatgen Molecule/Structure objects.
"""

__author__ = "Andrius Merkys"
__copyright__ = "Copyright 2015, AiiDA Development Team"
__version__ = "1.0"
__maintainer__ = "Andrius Merkys"
__email__ = "andrius.merkys@gmail.com"
__date__ = "Oct 9, 2015"

from monty.dev import requires

try:
    from aiida.orm import DataFactory
    from aiida.common.exceptions import MissingPluginError
    try:
        StructureData = DataFactory('structure')
    except MissingPluginError:
        raise ImportError
    aiida_loaded = True
except ImportError:
    aiida_loaded = False


@requires(aiida_loaded, "To use the AiidaStructureAdaptor, you need to have "
                        "aiida installed.")
class AiidaStructureAdaptor(object):
    """
    Adaptor serves as a bridge between AiiDA StructureData and pymatgen
    Molecule/Structure objects.
    """

    @staticmethod
    def get_structuredata(structure):
        """
        Returns AiiDA StructureData object from pymatgen structure or
        molecule.

        Args:
            structure: pymatgen.core.structure.Structure or
                       pymatgen.core.structure.Molecule

        Returns:
            AiiDA StructureData object
        """
        return StructureData(pymatgen=structure)

    @staticmethod
    def get_molecule(structuredata):
        """
        Returns pymatgen molecule from AiiDA StructureData.

        Args:
            structuredata: AiiDA StructureData object

        Returns:
            pymatgen.core.structure.Molecule
        """
        return structuredata.get_pymatgen_molecule()

    @staticmethod
    def get_structure(structuredata):
        """
        Returns pymatgen structure from AiiDA StructureData.

        Args:
            structuredata: AiiDA StructureData object

        Returns:
            pymatgen.core.structure.Structure
        """
        return structuredata.get_pymatgen_structure()
