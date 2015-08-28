# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

"""
Defines an abstract base class contract for Transformation object.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Sep 23, 2011"

import abc

from pymatgen.serializers.json_coders import PMGSONable
import six


class AbstractTransformation(six.with_metaclass(abc.ABCMeta, PMGSONable)):
    """
    Abstract transformation class.
    """

    @abc.abstractmethod
    def apply_transformation(self, structure):
        """
        Applies the transformation to a structure. Depending on whether a
        transformation is one-to-many, there may be an option to return a
        ranked list of structures.

        Args:
            structure:
                input structure
            return_ranked_list:
                Boolean stating whether or not multiple structures are
                returned. If return_ranked_list is a number, that number of
                structures is returned.

        Returns:
            depending on returned_ranked list, either a transformed structure
            or
            a list of dictionaries, where each dictionary is of the form
            {'structure' = .... , 'other_arguments'}
            the key 'transformation' is reserved for the transformation that
            was actually applied to the structure.
            This transformation is parsed by the alchemy classes for generating
            a more specific transformation history. Any other information will
            be stored in the transformation_parameters dictionary in the
            transmuted structure class.
        """
        return

    @abc.abstractproperty
    def inverse(self):
        """
        Returns the inverse transformation if available.
        Otherwise, should return None.
        """
        return

    @abc.abstractproperty
    def is_one_to_many(self):
        """
        Determines if a Transformation is a one-to-many transformation. If a
        Transformation is a one-to-many transformation, the
        apply_transformation method should have a keyword arg
        "return_ranked_list" which allows for the transformed structures to be
        returned as a ranked list.
        """
        return False

    @property
    def use_multiprocessing(self):
        """
        Indicates whether the transformation can be applied by a
        subprocessing pool. This should be overridden to return True for
        transformations that the transmuter can parallelize.
        """
        return False
