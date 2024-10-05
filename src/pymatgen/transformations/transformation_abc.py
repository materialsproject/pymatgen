"""Abstract base class for structure transformations."""

from __future__ import annotations

import abc
from typing import TYPE_CHECKING, Any

from monty.json import MSONable

if TYPE_CHECKING:
    from pymatgen.core import Structure

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Sep 23, 2011"


class AbstractTransformation(MSONable, abc.ABC):
    """Abstract transformation class."""

    @abc.abstractmethod
    def apply_transformation(self, structure: Structure) -> Structure | list[dict[str, Any]]:
        """Apply the transformation to a structure. Depending on whether a
        transformation is one-to-many, there may be an option to return a
        ranked list of structures.

        Args:
            structure:
                input structure
            return_ranked_list (bool | int, optional): If return_ranked_list is int, that number of structures

                is returned. If False, only the single lowest energy structure is returned. Defaults to False.

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

    @property
    def inverse(self) -> AbstractTransformation | None:
        """The inverse transformation if available.
        Otherwise, should return None. Defaults to None, so only need to
        override if applicable.
        """

    @property
    def is_one_to_many(self) -> bool:
        """Determine if a Transformation is a one-to-many transformation. In that case, the
        apply_transformation method should have a keyword arg "return_ranked_list" which
        allows for the transformed structures to be returned as a ranked list.
        Defaults to False, so only need to override if True.
        """
        return False

    @property
    def use_multiprocessing(self) -> bool:
        """Indicates whether the transformation can be applied by a
        subprocessing pool. This should be overridden to return True for
        transformations that the transmuter can parallelize.
        """
        return False
