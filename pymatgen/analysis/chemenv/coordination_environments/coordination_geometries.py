# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module contains the class describing the coordination geometries that can exist in a given structure. These
"model" coordination geometries are described in the following articles :
 - Pure Appl. Chem., Vol. 79, No. 10, pp. 1779--1799, 2007.
 - Acta Cryst. A, Vol. 46, No. 1, pp. 1--11, 1990.
The module also contains descriptors of part of these geometries (plane of separation, ...) that are used in the
identification algorithms.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"

import abc
import itertools
import json
import os

import numpy as np
from monty.json import MontyDecoder, MSONable
from scipy.special import factorial

module_dir = os.path.dirname(os.path.abspath(__file__))

UNKNOWN_ENVIRONMENT_SYMBOL = "UNKNOWN"
UNCLEAR_ENVIRONMENT_SYMBOL = "UNCLEAR"
EXPLICIT_PERMUTATIONS = "EXPLICIT_PERMUTATIONS"
SEPARATION_PLANE = "SEPARATION_PLANE"


class AbstractChemenvAlgorithm(MSONable, metaclass=abc.ABCMeta):
    """
    Base class used to define a Chemenv algorithm used to identify the correct permutation for the computation
    of the Continuous Symmetry Measure.
    """

    def __init__(self, algorithm_type):
        """
        Base constructor for ChemenvAlgorithm.

        Args:
            algorithm_type (str): Type of algorithm.
        """
        self._algorithm_type = algorithm_type

    @abc.abstractmethod
    def as_dict(self):
        """
        A JSON serializable dict representation of the algorithm
        """
        pass

    @property
    def algorithm_type(self):
        """
        Return the type of algorithm.

        Returns: Type of the algorithm
        """
        return self._algorithm_type

    @abc.abstractmethod
    def __str__(self):
        return ""


class ExplicitPermutationsAlgorithm(AbstractChemenvAlgorithm):
    """
    Class representing the algorithm doing the explicit permutations for the calculation of
    the Continuous Symmetry Measure.
    """

    def __init__(self, permutations):
        """
            Initializes a separation plane for a given perfect coordination geometry.

        Args:
            permutations: Permutations used for this algorithm.
        """
        super().__init__(algorithm_type=EXPLICIT_PERMUTATIONS)
        self._permutations = permutations

    def __str__(self):
        return self.algorithm_type

    @property
    def permutations(self):
        """
        Return the permutations to be performed for this algorithm.

        Returns: Permutations to be performed.
        """
        return self._permutations

    @property
    def as_dict(self):
        """
        Return the JSON serializable dict representation of this ExplicitPermutationsAlgorithm algorithm.

        Returns: a JSON serializable dict representation of this ExplicitPermutationsAlgorithm algorithm.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "permutations": self._permutations,
        }

    @classmethod
    def from_dict(cls, dd):
        """
        Reconstructs the ExplicitPermutationsAlgorithm algorithm from its JSON serializable dict representation.

        Args:
            dd: a JSON serializable dict representation of an ExplicitPermutationsAlgorithm algorithm.

        Returns: an ExplicitPermutationsAlgorithm algorithm.
        """
        return cls(dd["permutations"])


class SeparationPlane(AbstractChemenvAlgorithm):
    """
    Class representing the algorithm using separation planes for the calculation of
    the Continuous Symmetry Measure.
    """

    def __init__(
        self,
        plane_points,
        mirror_plane=False,
        ordered_plane=False,
        point_groups=None,
        ordered_point_groups=None,  # include_inverted_plane=False,
        # do_inverse_pt_gp_permutations=False, plane_type='MIRROR',
        explicit_permutations=None,
        minimum_number_of_points=None,
        explicit_optimized_permutations=None,
        multiplicity=None,
        other_plane_points=None,
    ):  # , plane_safe_permutations=False):
        """
            Initializes a separation plane for a given perfect coordination geometry

        Args:
            plane_points: Indices of the points that are in the plane in the perfect structure (and should be
                found in the defective one as well).
            mirror_plane: True if the separation plane is a mirror plane, in which case there is a correspondence
                of the points in each point_group (can reduce the number of permutations).
            ordered_plane: True if the order of the points in the plane can be taken into account to reduce the
                number of permutations.
            point_groups: Indices of the points in the two groups of points separated by the plane.
            ordered_point_groups: Whether the order of the points in each group of points can be taken into account to
                reduce the number of permutations.
            explicit_permutations: Explicit permutations to be performed in this separation plane algorithm.
            minimum_number_of_points: Minimum number of points needed to initialize a separation plane
                for this algorithm.
            explicit_optimized_permutations: Optimized set of explicit permutations to be performed in this
                separation plane algorithm.
            multiplicity: Number of such planes in the model geometry.
            other_plane_points: Indices of the points that are in the plane in the perfect structure for the other
                planes. The multiplicity should be equal to the length of this list + 1 ("main" separation plane +
                the other ones).
        """
        super().__init__(algorithm_type=SEPARATION_PLANE)
        self.mirror_plane = mirror_plane
        self.plane_points = plane_points
        self.point_groups = point_groups
        if len(point_groups[0]) > len(point_groups[1]):
            raise RuntimeError(
                "The number of points in the first group should be\n"
                "less than or equal to the number of points in the second group"
            )
        self._hash = 10000 * len(plane_points) + 100 * len(point_groups[0]) + len(point_groups[1])
        self.ordered_plane = ordered_plane
        self.ordered_point_groups = [False, False] if ordered_point_groups is None else ordered_point_groups
        # self._ordered_indices = list(point_groups[0])
        # self._ordered_indices.extend(plane_points)
        # self._ordered_indices.extend(point_groups[1])
        # self._inv_ordered_indices = np.argsort(self._ordered_indices)
        self.explicit_permutations = explicit_permutations
        self.explicit_optimized_permutations = explicit_optimized_permutations
        self._safe_permutations = None
        if self.explicit_optimized_permutations is not None:
            self._permutations = self.explicit_optimized_permutations
        elif self.explicit_permutations is not None:
            self._permutations = self.explicit_permutations
        self.multiplicity = multiplicity
        self.other_plane_points = other_plane_points
        self.minimum_number_of_points = minimum_number_of_points
        self.maximum_number_of_points = len(self.plane_points)
        self._ref_separation_perm = list(self.point_groups[0])
        self._ref_separation_perm.extend(list(self.plane_points))
        self._ref_separation_perm.extend(list(self.point_groups[1]))
        self._argsorted_ref_separation_perm = list(np.argsort(self._ref_separation_perm))
        self.separation = (
            len(point_groups[0]),
            len(plane_points),
            len(point_groups[1]),
        )

    # @property
    # def ordered_indices(self):
    #     """
    #     Ordered indices of the separation plane.
    #
    #     Examples:
    #          For a separation plane of type 2|4|3, with plane_points indices [0, 3, 5, 8] and
    #          point_groups indices [1, 4] and [2, 7, 6], the list of ordered indices is :
    #          [0, 3, 5, 8, 1, 4, 2, 7, 6].
    #
    #     Returns: list of ordered indices of this separation plane.
    #     """
    #     return self._ordered_indices
    #
    # @property
    # def inv_ordered_indices(self):
    #     return self._inv_ordered_indices

    @property
    def permutations(self):
        """
        Permutations used for this separation plane algorithm.

        Returns: List of permutations to be performed.
        """
        return self._permutations

    @property
    def ref_separation_perm(self):
        """
        Ordered indices of the separation plane.

        Examples:
             For a separation plane of type 2|4|3, with plane_points indices [0, 3, 5, 8] and
             point_groups indices [1, 4] and [2, 7, 6], the list of ordered indices is :
             [0, 3, 5, 8, 1, 4, 2, 7, 6].

        Returns: list of ordered indices of this separation plane.
        """
        return self._ref_separation_perm

    @property
    def argsorted_ref_separation_perm(self):
        """
        "Arg sorted" ordered indices of the separation plane.

        This is used in the identification of the final permutation to be used.

        Returns: list of the "arg sorted" ordered indices of the separation plane.
        """
        return self._argsorted_ref_separation_perm

    def safe_separation_permutations(self, ordered_plane=False, ordered_point_groups=None, add_opposite=False):
        """
        Simple and safe permutations for this separation plane.

        This is not meant to be used in production. Default configuration for ChemEnv does not use this method.
        Args:
            ordered_plane: Whether the order of the points in the plane can be used to reduce the
                number of permutations.
            ordered_point_groups: Whether the order of the points in each point group can be used to reduce the
                number of permutations.
            add_opposite: Whether to add the permutations from the second group before the first group as well.

        Returns: List of safe permutations.
        """
        s0 = list(range(len(self.point_groups[0])))
        plane = list(
            range(
                len(self.point_groups[0]),
                len(self.point_groups[0]) + len(self.plane_points),
            )
        )
        s1 = list(
            range(
                len(self.point_groups[0]) + len(self.plane_points),
                len(self.point_groups[0]) + len(self.plane_points) + len(self.point_groups[1]),
            )
        )
        ordered_point_groups = [False, False] if ordered_point_groups is None else ordered_point_groups

        def rotate(s, n):
            return s[-n:] + s[:-n]

        if ordered_plane and self.ordered_plane:
            plane_perms = [rotate(plane, ii) for ii in range(len(plane))]
            inv_plane = plane[::-1]
            plane_perms.extend([rotate(inv_plane, ii) for ii in range(len(inv_plane))])
        else:
            plane_perms = list(itertools.permutations(plane))
        if ordered_point_groups[0] and self.ordered_point_groups[0]:
            s0_perms = [rotate(s0, ii) for ii in range(len(s0))]
            inv_s0 = s0[::-1]
            s0_perms.extend([rotate(inv_s0, ii) for ii in range(len(inv_s0))])
        else:
            s0_perms = list(itertools.permutations(s0))
        if ordered_point_groups[1] and self.ordered_point_groups[1]:
            s1_perms = [rotate(s1, ii) for ii in range(len(s1))]
            inv_s1 = s1[::-1]
            s1_perms.extend([rotate(inv_s1, ii) for ii in range(len(inv_s1))])
        else:
            s1_perms = list(itertools.permutations(s1))
        if self._safe_permutations is None:
            self._safe_permutations = []
            for perm_side1 in s0_perms:
                for perm_sep_plane in plane_perms:
                    for perm_side2 in s1_perms:
                        perm = list(perm_side1)
                        perm.extend(list(perm_sep_plane))
                        perm.extend(list(perm_side2))
                        self._safe_permutations.append(perm)
                        if add_opposite:
                            perm = list(perm_side2)
                            perm.extend(list(perm_sep_plane))
                            perm.extend(list(perm_side1))
                            self._safe_permutations.append(perm)
        return self._safe_permutations

    @property
    def as_dict(self):
        """
        Return the JSON serializable dict representation of this SeparationPlane algorithm.

        Returns: a JSON serializable dict representation of this SeparationPlane algorithm.
        """
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "plane_points": self.plane_points,
            "mirror_plane": self.mirror_plane,
            "ordered_plane": self.ordered_plane,
            "point_groups": self.point_groups,
            "ordered_point_groups": self.ordered_point_groups,
            "explicit_permutations": [eperm.tolist() for eperm in self.explicit_permutations]
            if self.explicit_permutations is not None
            else None,
            "explicit_optimized_permutations": [eoperm.tolist() for eoperm in self.explicit_optimized_permutations]
            if self.explicit_optimized_permutations is not None
            else None,
            "multiplicity": self.multiplicity,
            "other_plane_points": self.other_plane_points,
            "minimum_number_of_points": self.minimum_number_of_points,
        }

    @classmethod
    def from_dict(cls, dd):
        """
        Reconstructs the SeparationPlane algorithm from its JSON serializable dict representation.

        Args:
            dd: a JSON serializable dict representation of an SeparationPlane algorithm.

        Returns: a SeparationPlane algorithm.
        """
        eop = (
            [np.array(eoperm) for eoperm in dd["explicit_optimized_permutations"]]
            if ("explicit_optimized_permutations" in dd and dd["explicit_optimized_permutations"] is not None)
            else None
        )
        return cls(
            plane_points=dd["plane_points"],
            mirror_plane=dd["mirror_plane"],
            ordered_plane=dd["ordered_plane"],
            point_groups=dd["point_groups"],
            ordered_point_groups=dd["ordered_point_groups"],
            explicit_permutations=[np.array(eperm) for eperm in dd["explicit_permutations"]],
            explicit_optimized_permutations=eop,
            multiplicity=dd["multiplicity"] if "multiplicity" in dd else None,
            other_plane_points=dd["other_plane_points"] if "other_plane_points" in dd else None,
            minimum_number_of_points=dd["minimum_number_of_points"],
        )

    def __str__(self):
        out = "Separation plane algorithm with the following reference separation :\n"
        out += "[{}] | [{}] | [{}]".format(
            "-".join(str(pp) for pp in [self.point_groups[0]]),
            "-".join(str(pp) for pp in [self.plane_points]),
            "-".join(str(pp) for pp in [self.point_groups[1]]),
        )
        return out


class CoordinationGeometry:
    """
    Class used to store the ideal representation of a chemical environment or "coordination geometry".
    """

    # Default value of continuous symmetry measure beyond which no further
    # search is performed for the separation plane algorithms
    CSM_SKIP_SEPARATION_PLANE_ALGO = 10.0

    class NeighborsSetsHints:
        """
        Class used to describe neighbors sets hints.

        This allows to possibly get a lower coordination from a capped-like model polyhedron.
        """

        ALLOWED_HINTS_TYPES = ["single_cap", "double_cap", "triple_cap"]

        def __init__(self, hints_type, options):
            """
            Constructor for this NeighborsSetsHints.

            Args:
                hints_type: type of hint (single, double or triple cap)
                options: options for the "hinting", e.g. the maximum csm value beyond which no additional
                    neighbors set could be found from a "cap hint".
            """
            if hints_type not in self.ALLOWED_HINTS_TYPES:
                raise ValueError('Type "{}" for NeighborsSetsHints is not allowed'.format(type))
            self.hints_type = hints_type
            self.options = options

        def hints(self, hints_info):
            """
            Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
            neighbors set.

            Args:
                hints_info: Info needed to build new "hinted" neighbors set.

            Returns: Voronoi indices of the new "hinted" neighbors set.
            """
            if hints_info["csm"] > self.options["csm_max"]:
                return []
            return object.__getattribute__(self, "{}_hints".format(self.hints_type))(hints_info)

        def single_cap_hints(self, hints_info):
            """
            Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
            neighbors set, in case of a "Single cap" hint.

            Args:
                hints_info: Info needed to build new "hinted" neighbors set.

            Returns: Voronoi indices of the new "hinted" neighbors set.
            """
            cap_index_perfect = self.options["cap_index"]
            nb_set = hints_info["nb_set"]
            permutation = hints_info["permutation"]
            nb_set_voronoi_indices_perfect_aligned = nb_set.get_neighb_voronoi_indices(permutation=permutation)
            cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[cap_index_perfect]
            new_site_voronoi_indices = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices.remove(cap_voronoi_index)
            return [new_site_voronoi_indices]

        def double_cap_hints(self, hints_info):
            """
            Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
            neighbors set, in case of a "Double cap" hint.

            Args:
                hints_info: Info needed to build new "hinted" neighbors set.

            Returns: Voronoi indices of the new "hinted" neighbors set.
            """
            first_cap_index_perfect = self.options["first_cap_index"]
            second_cap_index_perfect = self.options["second_cap_index"]
            nb_set = hints_info["nb_set"]
            permutation = hints_info["permutation"]
            nb_set_voronoi_indices_perfect_aligned = nb_set.get_neighb_voronoi_indices(permutation=permutation)
            first_cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[first_cap_index_perfect]
            second_cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[second_cap_index_perfect]
            new_site_voronoi_indices1 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices2 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices3 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices1.remove(first_cap_voronoi_index)
            new_site_voronoi_indices2.remove(second_cap_voronoi_index)
            new_site_voronoi_indices3.remove(first_cap_voronoi_index)
            new_site_voronoi_indices3.remove(second_cap_voronoi_index)
            return [
                new_site_voronoi_indices1,
                new_site_voronoi_indices2,
                new_site_voronoi_indices3,
            ]

        def triple_cap_hints(self, hints_info):
            """
            Return hints for an additional neighbors set, i.e. the voronoi indices that constitute this new
            neighbors set, in case of a "Triple cap" hint.

            Args:
                hints_info: Info needed to build new "hinted" neighbors set.

            Returns: Voronoi indices of the new "hinted" neighbors set.
            """
            first_cap_index_perfect = self.options["first_cap_index"]
            second_cap_index_perfect = self.options["second_cap_index"]
            third_cap_index_perfect = self.options["third_cap_index"]
            nb_set = hints_info["nb_set"]
            permutation = hints_info["permutation"]
            nb_set_voronoi_indices_perfect_aligned = nb_set.get_neighb_voronoi_indices(permutation=permutation)
            first_cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[first_cap_index_perfect]
            second_cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[second_cap_index_perfect]
            third_cap_voronoi_index = nb_set_voronoi_indices_perfect_aligned[third_cap_index_perfect]
            new_site_voronoi_indices1 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices2 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices3 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices4 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices5 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices6 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices7 = list(nb_set.site_voronoi_indices)
            new_site_voronoi_indices1.remove(first_cap_voronoi_index)
            new_site_voronoi_indices2.remove(second_cap_voronoi_index)
            new_site_voronoi_indices3.remove(third_cap_voronoi_index)
            new_site_voronoi_indices4.remove(second_cap_voronoi_index)
            new_site_voronoi_indices4.remove(third_cap_voronoi_index)
            new_site_voronoi_indices5.remove(first_cap_voronoi_index)
            new_site_voronoi_indices5.remove(third_cap_voronoi_index)
            new_site_voronoi_indices6.remove(first_cap_voronoi_index)
            new_site_voronoi_indices6.remove(second_cap_voronoi_index)
            new_site_voronoi_indices7.remove(first_cap_voronoi_index)
            new_site_voronoi_indices7.remove(second_cap_voronoi_index)
            new_site_voronoi_indices7.remove(third_cap_voronoi_index)
            return [
                new_site_voronoi_indices1,
                new_site_voronoi_indices2,
                new_site_voronoi_indices3,
                new_site_voronoi_indices4,
                new_site_voronoi_indices5,
                new_site_voronoi_indices6,
                new_site_voronoi_indices7,
            ]

        def as_dict(self):
            """
            A JSON serializable dict representation of this NeighborsSetsHints
            """
            return {"hints_type": self.hints_type, "options": self.options}

        @classmethod
        def from_dict(cls, dd):
            """
            Reconstructs the NeighborsSetsHints from its JSON serializable dict representation.

            Args:
                dd: a JSON serializable dict representation of a NeighborsSetsHints.

            Returns: a NeighborsSetsHints.
            """
            return cls(hints_type=dd["hints_type"], options=dd["options"])

    def __init__(
        self,
        mp_symbol,
        name,
        alternative_names=None,
        IUPAC_symbol=None,
        IUCr_symbol=None,
        coordination=None,
        central_site=np.zeros(3),
        points=None,
        solid_angles=None,
        permutations_safe_override=False,
        deactivate=False,
        faces=None,
        edges=None,
        algorithms=None,
        equivalent_indices=None,
        neighbors_sets_hints=None,
    ):
        """
        Initializes one "coordination geometry" according to [Pure Appl. Chem., Vol. 79, No. 10, pp. 1779--1799, 2007]
        and [Acta Cryst. A, Vol. 46, No. 1, pp. 1--11, 1990].

        Args:
            mp_symbol: Symbol used internally for the coordination geometry.
            name: Name of the coordination geometry.
            alternative_names: Alternative names for this coordination geometry.
            IUPAC_symbol: The IUPAC symbol of this coordination geometry.
            IUCr_symbol: The IUCr symbol of this coordination geometry.
            coordination: The coordination number of this coordination geometry (number of neighboring atoms).
            central_site: The coordinates of the central site of this coordination geometry.
            points: The list of the coordinates of all the points of this coordination geometry.
            solid_angles: The list of solid angles for each neighbor in this coordination geometry.
            permutations_safe_override: Computes all the permutations if set to True (overrides the plane separation
                algorithms or any other algorithm, for testing purposes)
            deactivate: Whether to deactivate this coordination geometry
            faces: List of the faces with their vertices given in a clockwise or anticlockwise order, for drawing
                purposes.
            edges: List of edges, for drawing purposes.
            algorithms: Algorithms used to identify this coordination geometry.
            equivalent_indices: The equivalent sets of indices in this coordination geometry (can be used to skip
                equivalent permutations that have already been performed).
            neighbors_sets_hints: Neighors sets hints for this coordination geometry.
        """
        self._mp_symbol = mp_symbol
        self.name = name
        self.alternative_names = alternative_names if alternative_names is not None else []
        self.IUPACsymbol = IUPAC_symbol
        self.IUCrsymbol = IUCr_symbol
        self.coordination = coordination
        self.central_site = np.array(central_site)
        self.points = points
        self._solid_angles = solid_angles
        self.permutations_safe_override = permutations_safe_override
        # self.plane_safe_permutations = plane_safe_permutations
        # self.setup_permutations(permutations)
        self.deactivate = deactivate
        self._faces = faces
        self._edges = edges
        self._algorithms = algorithms
        if points is not None:
            self.centroid = np.mean(np.array(points), axis=0)
        else:
            self.centroid = None
        self.equivalent_indices = equivalent_indices
        self.neighbors_sets_hints = neighbors_sets_hints
        self._pauling_stability_ratio = None

    def as_dict(self):
        """
        A JSON serializable dict representation of this CoordinationGeometry.
        """
        return {
            "mp_symbol": self._mp_symbol,
            "name": self.name,
            "alternative_names": self.alternative_names,
            "IUPAC_symbol": self.IUPACsymbol,
            "IUCr_symbol": self.IUCrsymbol,
            "coordination": self.coordination,
            "central_site": [float(xx) for xx in self.central_site],
            "points": [[float(xx) for xx in pp] for pp in self.points] if self.points is not None else None,
            "solid_angles": [float(ang) for ang in self._solid_angles] if self._solid_angles is not None else None,
            "deactivate": self.deactivate,
            "_faces": self._faces,
            "_edges": self._edges,
            "_algorithms": [algo.as_dict for algo in self._algorithms] if self._algorithms is not None else None,
            "equivalent_indices": self.equivalent_indices,
            "neighbors_sets_hints": [nbsh.as_dict() for nbsh in self.neighbors_sets_hints]
            if self.neighbors_sets_hints is not None
            else None,
        }

    @classmethod
    def from_dict(cls, dd):
        """
        Reconstructs the CoordinationGeometry from its JSON serializable dict representation.

        Args:
            dd: a JSON serializable dict representation of a CoordinationGeometry.

        Returns: a CoordinationGeometry.
        """
        dec = MontyDecoder()
        return cls(
            mp_symbol=dd["mp_symbol"],
            name=dd["name"],
            alternative_names=dd["alternative_names"],
            IUPAC_symbol=dd["IUPAC_symbol"],
            IUCr_symbol=dd["IUCr_symbol"],
            coordination=dd["coordination"],
            central_site=dd["central_site"],
            points=dd["points"],
            solid_angles=(
                dd["solid_angles"] if "solid_angles" in dd else [4.0 * np.pi / dd["coordination"]] * dd["coordination"]
            ),
            deactivate=dd["deactivate"],
            faces=dd["_faces"],
            edges=dd["_edges"],
            algorithms=[dec.process_decoded(algo_d) for algo_d in dd["_algorithms"]]
            if dd["_algorithms"] is not None
            else None,
            equivalent_indices=dd["equivalent_indices"] if "equivalent_indices" in dd else None,
            neighbors_sets_hints=[cls.NeighborsSetsHints.from_dict(nbshd) for nbshd in dd["neighbors_sets_hints"]]
            if ("neighbors_sets_hints" in dd and dd["neighbors_sets_hints"] is not None)
            else None,
        )

    def __str__(self):
        symbol = ""
        if self.IUPAC_symbol is not None:
            symbol += " (IUPAC: {s}".format(s=self.IUPAC_symbol)
            if self.IUCr_symbol is not None:
                symbol += " || IUCr: {s})".format(s=self.IUCr_symbol)
            else:
                symbol += ")"
        elif self.IUCr_symbol is not None:
            symbol += " (IUCr: {s})".format(s=self.IUCr_symbol)
        outs = [
            "Coordination geometry type : {n}{s}\n".format(n=self.name, s=symbol),
            "  - coordination number : {c}".format(c=self.coordination),
        ]
        if self.points is None:
            outs.append("... not yet implemented")
        else:
            outs.append("  - list of points :")
            for pp in self.points:
                outs.append("    - {p}".format(p=pp))
        outs.append("------------------------------------------------------------")
        outs.append("")

        return "\n".join(outs)

    def __repr__(self):
        symbol = ""
        if self.IUPAC_symbol is not None:
            symbol += " (IUPAC: {s}".format(s=self.IUPAC_symbol)
            if self.IUCr_symbol is not None:
                symbol += " || IUCr: {s})".format(s=self.IUCr_symbol)
            else:
                symbol += ")"
        elif self.IUCr_symbol is not None:
            symbol += " (IUCr: {s})".format(s=self.IUCr_symbol)
        outs = [
            "Coordination geometry type : {n}{s}\n".format(n=self.name, s=symbol),
            "  - coordination number : {c}".format(c=self.coordination),
        ]
        outs.append("------------------------------------------------------------")
        outs.append("")
        return "\n".join(outs)

    def __len__(self):
        return self.coordination

    def set_permutations_safe_override(self, permutations_safe_override):
        """
        Setup ChemEnv so that a safe set of permutations are used.

        Args:
            permutations_safe_override: Whether to use safe permutations.
        """
        self.permutations_safe_override = permutations_safe_override
        # self.setup_permutations()

    # @property
    # def csm_skip_algo(self):
    #     return self.CSM_SKIP_SEPARATION_PLANE_ALGO

    @property
    def distfactor_max(self):
        """
        The maximum distfactor for the perfect CoordinationGeometry.

        Returns: Maximum distfactor for the perfect CoordinationGeometry (usually 1.0 for symmetric polyhedrons).
        """
        dists = [np.linalg.norm(pp - self.central_site) for pp in self.points]
        return np.max(dists) / np.min(dists)

    @property
    def coordination_number(self):
        """
        Returns the coordination number of this coordination geometry.
        """
        return self.coordination

    @property
    def pauling_stability_ratio(self):
        """
        Returns the theoretical Pauling stability ratio (rC/rA) for this environment.
        """
        if self._pauling_stability_ratio is None:
            if self.ce_symbol in ["S:1", "L:2"]:
                self._pauling_stability_ratio = 0.0
            else:
                mindist_anions = 1000000.0
                mindist_cation_anion = 1000000.0
                for ipt1 in range(len(self.points)):
                    pt1 = np.array(self.points[ipt1])
                    mindist_cation_anion = min(mindist_cation_anion, np.linalg.norm(pt1 - self.central_site))
                    for ipt2 in range(ipt1 + 1, len(self.points)):
                        pt2 = np.array(self.points[ipt2])
                        mindist_anions = min(mindist_anions, np.linalg.norm(pt1 - pt2))
                anion_radius = mindist_anions / 2.0
                cation_radius = mindist_cation_anion - anion_radius
                self._pauling_stability_ratio = cation_radius / anion_radius
        return self._pauling_stability_ratio

    @property
    def mp_symbol(self):
        """
        Returns the MP symbol of this coordination geometry.
        """
        return self._mp_symbol

    @property
    def ce_symbol(self):
        """
        Returns the symbol of this coordination geometry.
        """
        return self._mp_symbol

    def get_coordination_number(self):
        """
        Returns the coordination number of this coordination geometry.
        """
        return self.coordination

    def is_implemented(self):
        """
        Returns True if this coordination geometry is implemented.
        """
        return bool(self.points)

    def get_name(self):
        """
        Returns the name of this coordination geometry.
        """
        return self.name

    @property
    def IUPAC_symbol(self):
        """
        Returns the IUPAC symbol of this coordination geometry.
        """
        return self.IUPACsymbol

    @property
    def IUPAC_symbol_str(self):
        """
        Returns a string representation of the IUPAC symbol of this coordination geometry.
        """
        return str(self.IUPACsymbol)

    @property
    def IUCr_symbol(self):
        """
        Returns the IUCr symbol of this coordination geometry.
        """
        return self.IUCrsymbol

    @property
    def IUCr_symbol_str(self):
        """
        Returns a string representation of the IUCr symbol of this coordination geometry.
        """
        return str(self.IUCrsymbol)

    @property
    def number_of_permutations(self):
        """
        Returns the number of permutations of this coordination geometry.
        """
        if self.permutations_safe_override:
            return factorial(self.coordination)
        if self.permutations is None:  # pylint: disable=E1101
            return factorial(self.coordination)
        return len(self.permutations)  # pylint: disable=E1101

    def ref_permutation(self, permutation):
        """
        Returns the reference permutation for a set of equivalent permutations.

        Can be useful to skip permutations that have already been performed.

        Args:
            permutation: Current permutation

        Returns: Reference permutation of the perfect CoordinationGeometry.
        """
        perms = []
        for eqv_indices in self.equivalent_indices:
            perms.append(tuple(permutation[ii] for ii in eqv_indices))
        perms.sort()
        return perms[0]

    @property
    def algorithms(self):
        """
        Returns the list of algorithms that are used to identify this coordination geometry.
        """
        return self._algorithms

    def get_central_site(self):
        """
        Returns the central site of this coordination geometry.
        """
        return self.central_site

    def faces(self, sites, permutation=None):
        """
        Returns the list of faces of this coordination geometry. Each face is given as a
        list of its vertices coordinates.
        """
        if permutation is None:
            coords = [site.coords for site in sites]
        else:
            coords = [sites[ii].coords for ii in permutation]
        return [[coords[ii] for ii in f] for f in self._faces]

    def edges(self, sites, permutation=None, input="sites"):
        """
        Returns the list of edges of this coordination geometry. Each edge is given as a
        list of its end vertices coordinates.
        """
        if input == "sites":
            coords = [site.coords for site in sites]
        elif input == "coords":
            coords = sites
        # if permutation is None:
        #     coords = [site.coords for site in sites]
        # else:
        #     coords = [sites[ii].coords for ii in permutation]
        if permutation is not None:
            coords = [coords[ii] for ii in permutation]
        return [[coords[ii] for ii in e] for e in self._edges]

    def solid_angles(self, permutation=None):
        """
        Returns the list of "perfect" solid angles Each edge is given as a
        list of its end vertices coordinates.
        """
        if permutation is None:
            return self._solid_angles
        return [self._solid_angles[ii] for ii in permutation]

    def get_pmeshes(self, sites, permutation=None):
        """
        Returns the pmesh strings used for jmol to show this geometry.
        """
        pmeshes = []
        # _vertices = [site.coords for site in sites]
        if permutation is None:
            _vertices = [site.coords for site in sites]
        else:
            _vertices = [sites[ii].coords for ii in permutation]
        _face_centers = []
        number_of_faces = 0
        for face in self._faces:
            if len(face) in [3, 4]:
                number_of_faces += 1
            else:
                number_of_faces += len(face)

            _face_centers.append(
                np.array([np.mean([_vertices[face_vertex][ii] for face_vertex in face]) for ii in range(3)])
            )

        out = "{}\n".format(len(_vertices) + len(_face_centers))
        for vv in _vertices:
            out += "{:15.8f} {:15.8f} {:15.8f}\n".format(vv[0], vv[1], vv[2])
        for fc in _face_centers:
            out += "{:15.8f} {:15.8f} {:15.8f}\n".format(fc[0], fc[1], fc[2])
        out += "{:d}\n".format(number_of_faces)
        for iface, face in enumerate(self._faces):
            if len(face) == 3:
                out += "4\n"
            elif len(face) == 4:
                out += "5\n"
            else:
                for ii, f in enumerate(face):
                    out += "4\n"
                    out += "{:d}\n".format(len(_vertices) + iface)
                    out += "{:d}\n".format(f)
                    out += "{:d}\n".format(face[np.mod(ii + 1, len(face))])
                    out += "{:d}\n".format(len(_vertices) + iface)
            if len(face) in [3, 4]:
                for face_vertex in face:
                    out += "{:d}\n".format(face_vertex)
                out += "{:d}\n".format(face[0])
        pmeshes.append({"pmesh_string": out})
        return pmeshes


class AllCoordinationGeometries(dict):
    """
    Class used to store all the reference "coordination geometries" (list with instances of the CoordinationGeometry
    classes)
    """

    def __init__(self, permutations_safe_override=False, only_symbols=None):
        """
        Initializes the list of Coordination Geometries.

        Args:
            permutations_safe_override: Whether to use safe permutations.
            only_symbols: Whether to restrict the list of environments to be identified.
        """
        dict.__init__(self)
        self.cg_list = list()
        if only_symbols is None:
            f = open("{}/coordination_geometries_files/allcg.txt".format(module_dir), "r")
            data = f.readlines()
            f.close()
            for line in data:
                cg_file = "{}/{}".format(module_dir, line.strip())
                f = open(cg_file, "r")
                dd = json.load(f)
                f.close()
                self.cg_list.append(CoordinationGeometry.from_dict(dd))
        else:
            for symbol in only_symbols:
                fsymbol = symbol.replace(":", "#")
                cg_file = "{}/coordination_geometries_files/{}.json".format(module_dir, fsymbol)
                f = open(cg_file, "r")
                dd = json.load(f)
                f.close()
                self.cg_list.append(CoordinationGeometry.from_dict(dd))

        self.cg_list.append(CoordinationGeometry(UNKNOWN_ENVIRONMENT_SYMBOL, "Unknown environment", deactivate=True))
        self.cg_list.append(CoordinationGeometry(UNCLEAR_ENVIRONMENT_SYMBOL, "Unclear environment", deactivate=True))
        if permutations_safe_override:
            for cg in self.cg_list:
                cg.set_permutations_safe_override(True)

        self.minpoints = {}
        self.maxpoints = {}
        self.separations_cg = {}
        for cn in range(6, 21):
            for cg in self.get_implemented_geometries(coordination=cn):
                if only_symbols is not None and cg.ce_symbol not in only_symbols:
                    continue
                if cn not in self.separations_cg:
                    self.minpoints[cn] = 1000
                    self.maxpoints[cn] = 0
                    self.separations_cg[cn] = {}
                for algo in cg.algorithms:
                    sep = (
                        len(algo.point_groups[0]),
                        len(algo.plane_points),
                        len(algo.point_groups[1]),
                    )
                    if sep not in self.separations_cg[cn]:
                        self.separations_cg[cn][sep] = []
                    self.separations_cg[cn][sep].append(cg.mp_symbol)
                    self.minpoints[cn] = min(self.minpoints[cn], algo.minimum_number_of_points)
                    self.maxpoints[cn] = max(self.maxpoints[cn], algo.maximum_number_of_points)
        self.maxpoints_inplane = {cn: max([sep[1] for sep in seps.keys()]) for cn, seps in self.separations_cg.items()}

    def __getitem__(self, key):
        return self.get_geometry_from_mp_symbol(key)

    def __contains__(self, item):
        try:
            self[item]
            return True
        except LookupError:
            return False

    def __repr__(self):
        """
        Returns a string with the list of coordination geometries.
        """
        outs = [
            "",
            "#=================================#",
            "# List of coordination geometries #",
            "#=================================#",
            "",
        ]
        for cg in self.cg_list:
            outs.append(repr(cg))

        return "\n".join(outs)

    def __str__(self):
        """
        Returns a string with the list of coordination geometries that are implemented.
        """
        outs = [
            "",
            "#=======================================================#",
            "# List of coordination geometries currently implemented #",
            "#=======================================================#",
            "",
        ]
        for cg in self.cg_list:
            if cg.is_implemented():
                outs.append(str(cg))

        return "\n".join(outs)

    def get_geometries(self, coordination=None, returned="cg"):
        """
        Returns a list of coordination geometries with the given coordination number.

        Args:
            coordination: The coordination number of which the list of coordination geometries are returned.
            returned: Type of objects in the list.
        """
        geom = list()
        if coordination is None:
            for gg in self.cg_list:
                if returned == "cg":
                    geom.append(gg)
                elif returned == "mp_symbol":
                    geom.append(gg.mp_symbol)
        else:
            for gg in self.cg_list:
                if gg.get_coordination_number() == coordination:
                    if returned == "cg":
                        geom.append(gg)
                    elif returned == "mp_symbol":
                        geom.append(gg.mp_symbol)
        return geom

    def get_symbol_name_mapping(self, coordination=None):
        """
        Return a dictionary mapping the symbol of a CoordinationGeometry to its name.

        Args:
            coordination: Whether to restrict the dictionary to a given coordination.

        Returns: Dictionary mapping the symbol of a CoordinationGeometry to its name.
        """
        geom = {}
        if coordination is None:
            for gg in self.cg_list:
                geom[gg.mp_symbol] = gg.name
        else:
            for gg in self.cg_list:
                if gg.get_coordination_number() == coordination:
                    geom[gg.mp_symbol] = gg.name
        return geom

    def get_symbol_cn_mapping(self, coordination=None):
        """
        Return a dictionary mapping the symbol of a CoordinationGeometry to its coordination.

        Args:
            coordination: Whether to restrict the dictionary to a given coordination.

        Returns: Dictionary mapping the symbol of a CoordinationGeometry to its coordination.
        """
        geom = {}
        if coordination is None:
            for gg in self.cg_list:
                geom[gg.mp_symbol] = gg.coordination_number
        else:
            for gg in self.cg_list:
                if gg.get_coordination_number() == coordination:
                    geom[gg.mp_symbol] = gg.coordination_number
        return geom

    def get_implemented_geometries(self, coordination=None, returned="cg", include_deactivated=False):
        """
        Returns a list of the implemented coordination geometries with the given coordination number.

        Args:
            coordination: The coordination number of which the list of implemented coordination geometries
                are returned.
            returned: Type of objects in the list.
            include_deactivated: Whether to include CoordinationGeometry that are deactivated.
        """
        geom = list()
        if coordination is None:
            for gg in self.cg_list:
                if gg.points is not None and ((not gg.deactivate) or include_deactivated):
                    if returned == "cg":
                        geom.append(gg)
                    elif returned == "mp_symbol":
                        geom.append(gg.mp_symbol)
        else:
            for gg in self.cg_list:
                if (
                    gg.get_coordination_number() == coordination
                    and gg.points is not None
                    and ((not gg.deactivate) or include_deactivated)
                ):
                    if returned == "cg":
                        geom.append(gg)
                    elif returned == "mp_symbol":
                        geom.append(gg.mp_symbol)
        return geom

    def get_not_implemented_geometries(self, coordination=None, returned="mp_symbol"):
        """
        Returns a list of the implemented coordination geometries with the given coordination number.

        Args:
            coordination: The coordination number of which the list of implemented coordination geometries
                are returned.
            returned: Type of objects in the list.
        """
        geom = list()
        if coordination is None:
            for gg in self.cg_list:
                if gg.points is None:
                    if returned == "cg":
                        geom.append(gg)
                    elif returned == "mp_symbol":
                        geom.append(gg.mp_symbol)
        else:
            for gg in self.cg_list:
                if gg.get_coordination_number() == coordination and gg.points is None:
                    if returned == "cg":
                        geom.append(gg)
                    elif returned == "mp_symbol":
                        geom.append(gg.mp_symbol)
        return geom

    def get_geometry_from_name(self, name):
        """
        Returns the coordination geometry of the given name.

        Args:
            name: The name of the coordination geometry.
        """
        for gg in self.cg_list:
            if gg.name == name or name in gg.alternative_names:
                return gg
        raise LookupError('No coordination geometry found with name "{name}"'.format(name=name))

    def get_geometry_from_IUPAC_symbol(self, IUPAC_symbol):
        """
        Returns the coordination geometry of the given IUPAC symbol.

        Args:
            IUPAC_symbol: The IUPAC symbol of the coordination geometry.
        """
        for gg in self.cg_list:
            if gg.IUPAC_symbol == IUPAC_symbol:
                return gg
        raise LookupError('No coordination geometry found with IUPAC symbol "{symbol}"'.format(symbol=IUPAC_symbol))

    def get_geometry_from_IUCr_symbol(self, IUCr_symbol):
        """
        Returns the coordination geometry of the given IUCr symbol.

        Args:
            IUCr_symbol: The IUCr symbol of the coordination geometry.
        """
        for gg in self.cg_list:
            if gg.IUCr_symbol == IUCr_symbol:
                return gg
        raise LookupError('No coordination geometry found with IUCr symbol "{symbol}"'.format(symbol=IUCr_symbol))

    def get_geometry_from_mp_symbol(self, mp_symbol):
        """
        Returns the coordination geometry of the given mp_symbol.

        Args:
            mp_symbol: The mp_symbol of the coordination geometry.
        """
        for gg in self.cg_list:
            if gg.mp_symbol == mp_symbol:
                return gg
        raise LookupError('No coordination geometry found with mp_symbol "{symbol}"'.format(symbol=mp_symbol))

    def is_a_valid_coordination_geometry(self, mp_symbol=None, IUPAC_symbol=None, IUCr_symbol=None, name=None, cn=None):
        """
        Checks whether a given coordination geometry is valid (exists) and whether the parameters are coherent with
        each other.

        Args:
            mp_symbol: The mp_symbol of the coordination geometry.
            IUPAC_symbol: The IUPAC_symbol of the coordination geometry.
            IUCr_symbol: The IUCr_symbol of the coordination geometry.
            name: The name of the coordination geometry.
            cn: The coordination of the coordination geometry.
        """
        if name is not None:
            raise NotImplementedError("is_a_valid_coordination_geometry not implemented for the name")
        if mp_symbol is None and IUPAC_symbol is None and IUCr_symbol is None:
            raise SyntaxError(
                "missing argument for is_a_valid_coordination_geometry : at least one of mp_symbol, "
                "IUPAC_symbol and IUCr_symbol must be passed to the function"
            )
        if mp_symbol is not None:
            try:
                cg = self.get_geometry_from_mp_symbol(mp_symbol)
                if IUPAC_symbol is not None:
                    if IUPAC_symbol != cg.IUPAC_symbol:
                        return False
                if IUCr_symbol is not None:
                    if IUCr_symbol != cg.IUCr_symbol:
                        return False
                if cn is not None:
                    if int(cn) != int(cg.coordination_number):
                        return False
                return True
            except LookupError:
                return False
        elif IUPAC_symbol is not None:
            try:
                cg = self.get_geometry_from_IUPAC_symbol(IUPAC_symbol)
                if IUCr_symbol is not None:
                    if IUCr_symbol != cg.IUCr_symbol:
                        return False
                if cn is not None:
                    if cn != cg.coordination_number:
                        return False
                return True
            except LookupError:
                return False
        elif IUCr_symbol is not None:
            try:
                cg = self.get_geometry_from_IUCr_symbol(IUCr_symbol)
                if cn is not None:
                    if cn != cg.coordination_number:
                        return False
                return True
            except LookupError:
                return True
        raise Exception("Should not be here !")

    def pretty_print(self, type="implemented_geometries", maxcn=8, additional_info=None):
        """
        Return a string with a list of the Coordination Geometries.

        Args:
            type: Type of string to be returned (all_geometries, all_geometries_latex_images, all_geometries_latex,
                implemented_geometries).
            maxcn: Maximum coordination.
            additional_info: Whether to add some additional info for each coordination geometry.

        Returns: String describing the list of coordination geometries.
        """
        if type == "all_geometries_latex_images":
            mystring = ""
            for cn in range(1, maxcn + 1):
                mystring += "\\section*{{Coordination {cn}}}\n\n".format(cn=cn)
                for cg in self.get_implemented_geometries(coordination=cn, returned="cg"):
                    mystring += "\\subsubsection*{{{mp} : {name}}}\n\n".format(mp=cg.mp_symbol, name=cg.get_name())
                    mystring += "IUPAC : {iupac}\n\nIUCr : {iucr}\n\n".format(
                        iupac=cg.IUPAC_symbol, iucr=cg.IUCr_symbol
                    )
                    mystring += "\\begin{center}\n"
                    mystring += "\\includegraphics[scale=0.15]{{images/{let}_{cif}.png}}\n".format(
                        let=cg.mp_symbol.split(":")[0], cif=cg.mp_symbol.split(":")[1]
                    )
                    mystring += "\\end{center}\n\n"
                for cg in self.get_not_implemented_geometries(cn, returned="cg"):
                    mystring += "\\subsubsection*{{{mp} : {name}}}\n\n".format(mp=cg.mp_symbol, name=cg.get_name())
                    mystring += "IUPAC : {iupac}\n\nIUCr : {iucr}\n\n".format(
                        iupac=cg.IUPAC_symbol, iucr=cg.IUCr_symbol
                    )
        elif type == "all_geometries_latex":
            mystring = ""
            for cn in range(1, maxcn + 1):
                mystring += "\\subsection*{{Coordination {cn}}}\n\n".format(cn=cn)
                mystring += "\\begin{itemize}\n"
                for cg in self.get_implemented_geometries(coordination=cn, returned="cg"):
                    mystring += "\\item {mp} $\\rightarrow$ {name} ".format(
                        mp=cg.mp_symbol.replace("_", "\\_"), name=cg.get_name()
                    )
                    mystring += "(IUPAC : {iupac} - IUCr : {iucr})\n".format(
                        iupac=cg.IUPAC_symbol_str,
                        iucr=cg.IUCr_symbol_str.replace("[", "$[$").replace("]", "$]$"),
                    )
                for cg in self.get_not_implemented_geometries(cn, returned="cg"):
                    mystring += "\\item {mp} $\\rightarrow$ {name} ".format(
                        mp=cg.mp_symbol.replace("_", "\\_"), name=cg.get_name()
                    )
                    mystring += "(IUPAC : {iupac} - IUCr : {iucr})\n".format(
                        iupac=cg.IUPAC_symbol_str,
                        iucr=cg.IUCr_symbol_str.replace("[", "$[$").replace("]", "$]$"),
                    )
                mystring += "\\end{itemize}\n\n"
        else:
            mystring = "+-------------------------+\n| Coordination geometries |\n+-------------------------+\n\n"
            for cn in range(1, maxcn + 1):
                mystring += "==>> CN = {cn} <<==\n".format(cn=cn)
                if type == "implemented_geometries":
                    for cg in self.get_implemented_geometries(coordination=cn):
                        if additional_info is not None:
                            if "nb_hints" in additional_info:
                                if cg.neighbors_sets_hints is not None:
                                    addinfo = " *"
                                else:
                                    addinfo = ""
                            else:
                                addinfo = ""
                        else:
                            addinfo = ""
                        mystring += " - {mp} : {name}{addinfo}\n".format(
                            mp=cg.mp_symbol, name=cg.get_name(), addinfo=addinfo
                        )
                elif type == "all_geometries":
                    for cg in self.get_geometries(coordination=cn):
                        mystring += " - {mp} : {name}\n".format(mp=cg.mp_symbol, name=cg.get_name())
                mystring += "\n"
        return mystring
