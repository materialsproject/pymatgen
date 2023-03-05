# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module defines filters for Transmuter object.
"""

from __future__ import annotations

import abc
from collections import defaultdict

from monty.json import MSONable

from pymatgen.analysis.structure_matcher import ElementComparator, StructureMatcher
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class AbstractStructureFilter(MSONable, metaclass=abc.ABCMeta):
    """
    AbstractStructureFilter that defines an API to perform testing of
    Structures. Structures that return True to a test are retained during
    transmutation while those that return False are removed.
    """

    @abc.abstractmethod
    def test(self, structure):
        """
        Method to execute the test.

        Args:
            structure (Structure): Input structure to test

        Returns:
            (bool) Structures that return true are kept in the Transmuter
            object during filtering.
        """
        return


class ContainsSpecieFilter(AbstractStructureFilter):
    """
    Filter for structures containing certain elements or species.
    By default compares by atomic number.
    """

    def __init__(self, species, strict_compare=False, AND=True, exclude=False):
        """
        Args:
            species ([Species/Element]): list of species to look for
            AND: whether all species must be present to pass (or fail) filter.
            strict_compare: if true, compares objects by specie or element
                object if false, compares atomic number
            exclude: If true, returns false for any structures with the specie
                (excludes them from the Transmuter)
        """
        self._species = list(map(get_el_sp, species))
        self._strict = strict_compare
        self._AND = AND
        self._exclude = exclude

    def test(self, structure):
        """
        Method to execute the test.

        Returns: True if structure do not contain specified species.
        """
        # set up lists to compare
        if not self._strict:
            # compare by atomic number
            filter_set = {sp.Z for sp in self._species}
            structure_set = {sp.Z for sp in structure.composition.elements}
        else:
            # compare by specie or element object
            filter_set = set(self._species)
            structure_set = set(structure.composition.elements)

        if self._AND and filter_set <= structure_set:
            # return true if we aren't excluding since all are in structure
            return not self._exclude
        if (not self._AND) and filter_set & structure_set:
            # return true if we aren't excluding since one is in structure
            return not self._exclude
        # return false if we aren't excluding otherwise
        return self._exclude

    def __repr__(self):
        return "\n".join(
            [
                "ContainsSpecieFilter with parameters:",
                f"species = {self._species}",
                f"strict_compare = {self._strict}",
                f"AND = {self._AND}",
                f"exclude = {self._exclude}",
            ]
        )

    def as_dict(self):
        """
        Returns: MSONAble dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "init_args": {
                "species": [str(sp) for sp in self._species],
                "strict_compare": self._strict,
                "AND": self._AND,
                "exclude": self._exclude,
            },
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            Filter
        """
        return cls(**d["init_args"])


class SpecieProximityFilter(AbstractStructureFilter):
    """
    This filter removes structures that have certain species that are too close
    together.
    """

    def __init__(self, specie_and_min_dist_dict):
        """
        Args:
            specie_and_min_dist_dict (dict): A species string to float mapping. For
                example, {"Na+": 1} means that all Na+ ions must be at least 1
                Angstrom away from each other. Multiple species criteria can be
                applied. Note that the testing is done based on the actual object
                . If you have a structure with Element, you must use {"Na":1}
                instead to filter based on Element and not Species.
        """
        self.specie_and_min_dist = {get_el_sp(k): v for k, v in specie_and_min_dist_dict.items()}

    def test(self, structure):
        """
        Method to execute the test.

        Args:
            structure (Structure): Input structure to test

        Returns: True if structure does not contain species within specified
            distances.
        """
        all_species = set(self.specie_and_min_dist)
        for site in structure:
            species = set(site.species)
            sp_to_test = species.intersection(all_species)
            if sp_to_test:
                max_r = max(self.specie_and_min_dist[sp] for sp in sp_to_test)
                nn = structure.get_neighbors(site, max_r)
                for sp in sp_to_test:
                    for nn_site, dist, *_ in nn:
                        if sp in nn_site.species and dist < self.specie_and_min_dist[sp]:
                            return False
        return True

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "init_args": {"specie_and_min_dist_dict": {str(sp): v for sp, v in self.specie_and_min_dist.items()}},
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            Filter
        """
        return cls(**d["init_args"])


class RemoveDuplicatesFilter(AbstractStructureFilter):
    """
    This filter removes exact duplicate structures from the transmuter.
    """

    def __init__(self, structure_matcher: dict | StructureMatcher | None = None, symprec: float = None) -> None:
        """
        Remove duplicate structures based on the structure matcher
        and symmetry (if symprec is given).

        Args:
            structure_matcher (dict | StructureMatcher, optional): Provides a structure matcher to be used for
                structure comparison.
            symprec (float, optional): The precision in the symmetry finder algorithm if None (
                default value), no symmetry check is performed and only the
                structure matcher is used. A recommended value is 1e-5.
        """
        self.symprec = symprec
        self.structure_list: dict[str, list[Structure]] = defaultdict(list)
        if not isinstance(structure_matcher, (dict, StructureMatcher, type(None))):
            raise ValueError(f"structure_matcher must be a dict, StructureMatcher or None, got {structure_matcher}")
        if isinstance(structure_matcher, dict):
            self.structure_matcher = StructureMatcher.from_dict(structure_matcher)
        else:
            self.structure_matcher = structure_matcher or StructureMatcher(comparator=ElementComparator())

    def test(self, structure):
        """
        Args:
            structure (Structure): Input structure to test

        Returns: True if structure is not in list.
        """
        hash = self.structure_matcher._comparator.get_hash(structure.composition)
        if not self.structure_list[hash]:
            self.structure_list[hash].append(structure)
            return True

        def get_spg_num(struct: structure) -> int:
            finder = SpacegroupAnalyzer(struct, symprec=self.symprec)
            return finder.get_space_group_number()

        for s in self.structure_list[hash]:
            if (self.symprec is None or get_spg_num(s) == get_spg_num(structure)) and self.structure_matcher.fit(
                s, structure
            ):
                return False

        self.structure_list[hash].append(structure)
        return True


class RemoveExistingFilter(AbstractStructureFilter):
    """
    This filter removes structures existing in a given list from the transmuter.
    """

    def __init__(self, existing_structures, structure_matcher=None, symprec=None):
        """
        Remove existing structures based on the structure matcher
        and symmetry (if symprec is given).

        Args:
            existing_structures: List of existing structures to compare with
            structure_matcher: Provides a structure matcher to be used for
                structure comparison.
            symprec: The precision in the symmetry finder algorithm if None (
                default value), no symmetry check is performed and only the
                structure matcher is used. A recommended value is 1e-5.
        """
        self.symprec = symprec
        self.structure_list = []
        self.existing_structures = existing_structures
        if isinstance(structure_matcher, dict):
            self.structure_matcher = StructureMatcher.from_dict(structure_matcher)
        else:
            self.structure_matcher = structure_matcher or StructureMatcher(comparator=ElementComparator())

    def test(self, structure):
        """
        Method to execute the test.

        Args:
            structure (Structure): Input structure to test

        Returns: True if structure is not in existing list.
        """

        def get_sg(s):
            finder = SpacegroupAnalyzer(s, symprec=self.symprec)
            return finder.get_space_group_number()

        for s in self.existing_structures:
            if (
                self.structure_matcher._comparator.get_hash(structure.composition)
                == self.structure_matcher._comparator.get_hash(s.composition)
                and self.symprec is None
                or get_sg(s) == get_sg(structure)
            ) and self.structure_matcher.fit(s, structure):
                return False

        self.structure_list.append(structure)
        return True

    def as_dict(self):
        """
        Returns: MSONable dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "init_args": {"structure_matcher": self.structure_matcher.as_dict()},
        }


class ChargeBalanceFilter(AbstractStructureFilter):
    """
    This filter removes structures that are not charge balanced from the
    transmuter. This only works if the structure is oxidation state
    decorated, as structures with only elemental sites are automatically
    assumed to have net charge of 0.
    """

    def __init__(self):
        """
        No args required.
        """

    def test(self, structure):
        """
        Method to execute the test.

        Args:
            structure (Structure): Input structure to test

        Returns: True if structure is neutral.
        """
        return structure.charge == 0.0


class SpeciesMaxDistFilter(AbstractStructureFilter):
    """
    This filter removes structures that do have two particular species that are
    not nearest neighbors by a predefined max_dist. For instance, if you are
    analyzing Li battery materials, you would expect that each Li+ would be
    nearest neighbor to lower oxidation state transition metal for
    electrostatic reasons. This only works if the structure is oxidation state
    decorated, as structures with only elemental sites are automatically
    assumed to have net charge of 0.
    """

    def __init__(self, sp1, sp2, max_dist):
        """
        Args:
            sp1 (Species): First specie
            sp2 (Species): Second specie
            max_dist (float): Maximum distance between species.
        """
        self.sp1 = get_el_sp(sp1)
        self.sp2 = get_el_sp(sp2)
        self.max_dist = max_dist

    def test(self, structure):
        """
        Method to execute the test.

        Args:
            structure (Structure): Input structure to test

        Returns: True if structure does not contain the two species are distances
            greater than max_dist.
        """
        sp1_indices = [i for i, site in enumerate(structure) if site.specie == self.sp1]
        sp2_indices = [i for i, site in enumerate(structure) if site.specie == self.sp2]
        fcoords = structure.frac_coords
        fcoords1 = fcoords[sp1_indices, :]
        fcoords2 = fcoords[sp2_indices, :]
        lattice = structure.lattice
        dists = lattice.get_all_distances(fcoords1, fcoords2)
        return all(any(row) for row in dists < self.max_dist)
