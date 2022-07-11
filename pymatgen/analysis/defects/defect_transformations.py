# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines classes for point defect transformations on structures
"""

from pymatgen.transformations.transformation_abc import AbstractTransformation

__author__ = "Danny Broberg, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"


class DefectTransformation(AbstractTransformation):
    """
    Generates Defect structures based on pymatgen Defect Core classes
    """

    def __init__(self, scaling_matrix, defect):
        """
        :param scaling_matrix: Supercell scaling matrix
        :param defect: Defect pymatgen object
            NOTE: defect.bulk_structure should be same as provided structure in the apply_transformation step
        """
        self.scaling_matrix = scaling_matrix
        self.defect = defect

    def apply_transformation(self, structure):
        """
        :param structure: (bulk structure to be scaled up - typically conventional unit cell)
        :return: defect_structure, with charge applied
        """
        if structure != self.defect.bulk_structure:
            raise ValueError("Defect bulk_structure is not the same as input structure.")

        def_structure = self.defect.generate_defect_structure(self.scaling_matrix)

        return def_structure

    def __str__(self):
        inp_args = [
            f"Supercell scaling matrix = {self.scaling_matrix}",
            f"Defect = {self.defect}",
        ]
        return "Defect Transformation : " + ", ".join(inp_args)

    def __repr__(self):
        return self.__str__()

    @property
    def inverse(self):
        """
        Not implemented
        """
        raise NotImplementedError()

    @property
    def is_one_to_many(self):
        """
        Returns: False
        """
        return False
