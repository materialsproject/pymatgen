# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import abc

import six
from six.moves import map

from pymatgen.core.periodic_table import get_el_sp
from monty.json import MSONable
from pymatgen.analysis.structure_matcher import StructureMatcher,\
    ElementComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from collections import defaultdict

"""
This module defines filters for Transmuter object.
"""


__author__ = "Will Richards, Shyue Ping Ong, Stephen Dacek"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Will Richards"
__email__ = "wrichards@mit.edu"
__date__ = "Sep 25, 2012"


class AbstractStructureFilter(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    AbstractStructureFilter that defines an API to perform testing of
    Structures. Structures that return True to a test are retained during
    transmutation while those that return False are removed.
    """

    @abc.abstractmethod
    def test(self, structure):
        """
        Method to execute the test.

        Returns:
            (bool) Structures that return true are kept in the Transmuter
            object during filtering.
        """
        return


class ContainsSpecieFilter(AbstractStructureFilter):

    def __init__(self, species, strict_compare=False, AND=True, exclude=False):
        """
        Filter for structures containing certain elements or species.
        By default compares by atomic number

        Args:
            species ([Specie/Element]): list of species to look for
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
        #set up lists to compare
        if not self._strict:
            #compare by atomic number
            atomic_number = lambda x: x.Z
            filter_set = set(map(atomic_number, self._species))
            structure_set = set(map(atomic_number,
                                structure.composition.elements))
        else:
            #compare by specie or element object
            filter_set = set(self._species)
            structure_set = set(structure.composition.elements)

        if self._AND and filter_set <= structure_set:
            #return true if we aren't excluding since all are in structure
            return not self._exclude
        elif (not self._AND) and filter_set & structure_set:
            #return true if we aren't excluding since one is in structure
            return not self._exclude
        else:
            #return false if we aren't excluding otherwise
            return self._exclude

    def __repr__(self):
        return "\n".join(["ContainsSpecieFilter with parameters:",
                          "species = {}".format(self._species),
                          "strict_compare = {}".format(self._strict),
                          "AND = {}".format(self._AND),
                          "exclude = {}".format(self._exclude)])

    def as_dict(self):
        return {"version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "init_args": {"species": [str(sp) for sp in self._species],
                              "strict_compare": self._strict,
                              "AND": self._AND,
                              "exclude": self._exclude}}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])


class SpecieProximityFilter(AbstractStructureFilter):
    """
    This filter removes structures that have certain species that are too close
    together.

    Args:
        specie_and_min_dist_dict: A species string to float mapping. For
            example, {"Na+": 1} means that all Na+ ions must be at least 1
            Angstrom away from each other. Multiple species criteria can be
            applied. Note that the testing is done based on the actual object
            . If you have a structure with Element, you must use {"Na":1}
            instead to filter based on Element and not Specie.

    """

    def __init__(self, specie_and_min_dist_dict):
        self.specie_and_min_dist = {get_el_sp(k): v
                                    for k, v
                                    in specie_and_min_dist_dict.items()}

    def test(self, structure):
        all_species = set(self.specie_and_min_dist.keys())
        for site in structure:
            species = site.species_and_occu.keys()
            sp_to_test = set(species).intersection(all_species)
            if sp_to_test:
                max_r = max([self.specie_and_min_dist[sp]
                             for sp in sp_to_test])
                nn = structure.get_neighbors(site, max_r)
                for sp in sp_to_test:
                    for (nnsite, dist) in nn:
                        if sp in nnsite.species_and_occu.keys():
                            if dist < self.specie_and_min_dist[sp]:
                                return False
        return True

    def as_dict(self):
        return {"version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "init_args": {"specie_and_min_dist_dict":
                              {str(sp): v
                               for sp, v in self.specie_and_min_dist.items()}}}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])


class RemoveDuplicatesFilter(AbstractStructureFilter):
    """
    This filter removes exact duplicate structures from the transmuter.
    """

    def __init__(self, structure_matcher=StructureMatcher(
                 comparator=ElementComparator()), symprec=None):
        """
        Remove duplicate structures based on the structure matcher
        and symmetry (if symprec is given).

        Args:
            structure_matcher: Provides a structure matcher to be used for
                structure comparison.
            symprec: The precision in the symmetry finder algorithm if None (
                default value), no symmetry check is performed and only the
                structure matcher is used. A recommended value is 1e-5.
        """
        self.symprec = symprec
        self.structure_list = defaultdict(list)
        if isinstance(structure_matcher, dict):
            self.structure_matcher = StructureMatcher.from_dict(structure_matcher)
        else:
            self.structure_matcher = structure_matcher

    def test(self, structure):
        h = self.structure_matcher._comparator.get_hash(structure.composition)
        if not self.structure_list[h]:
            self.structure_list[h].append(structure)
            return True

        def get_sg(s):
            finder = SpacegroupAnalyzer(s, symprec=self.symprec)
            return finder.get_space_group_number()

        for s in self.structure_list[h]:
            if self.symprec is None or \
                    get_sg(s) == get_sg(structure):
                if self.structure_matcher.fit(s, structure):
                    return False

        self.structure_list[h].append(structure)
        return True


class RemoveExistingFilter(AbstractStructureFilter):
    """
    This filter removes structures existing in a given list from the transmuter.
    """

    def __init__(self, existing_structures, structure_matcher=StructureMatcher(
                 comparator=ElementComparator()), symprec=None):
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
            self.structure_matcher = structure_matcher

    def test(self, structure):

        def get_sg(s):
            finder = SpacegroupAnalyzer(s, symprec=self.symprec)
            return finder.get_space_group_number()

        for s in self.existing_structures:
            if self.structure_matcher._comparator.get_hash(structure.composition) ==\
                    self.structure_matcher._comparator.get_hash(s.composition):
                if self.symprec is None or \
                        get_sg(s) == get_sg(structure):
                    if self.structure_matcher.fit(s, structure):
                        return False

        self.structure_list.append(structure)
        return True

    def as_dict(self):
        return {"version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "init_args": {"structure_matcher": self.structure_matcher.as_dict()}}


class ChargeBalanceFilter(AbstractStructureFilter):
    """
    This filter removes structures that are not charge balanced from the
    transmuter. This only works if the structure is oxidation state
    decorated, as structures with only elemental sites are automatically
    assumed to have net charge of 0.
    """
    def __init__(self):
        pass

    def test(self, structure):
        if structure.charge == 0.0:
            return True
        else:
            return False


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
        self.sp1 = get_el_sp(sp1)
        self.sp2 = get_el_sp(sp2)
        self.max_dist = max_dist

    def test(self, structure):
        sp1_indices = [i for i, site in enumerate(structure) if
                       site.specie == self.sp1]
        sp2_indices = [i for i, site in enumerate(structure) if
                       site.specie == self.sp2]
        fcoords = structure.frac_coords
        fcoords1 = fcoords[sp1_indices, :]
        fcoords2 = fcoords[sp2_indices, :]
        lattice = structure.lattice
        dists = lattice.get_all_distances(fcoords1, fcoords2)
        return all([any(row) for row in dists < self.max_dist])
