#!/usr/bin/env python

"""
This module defines filters for Transmuter object.
"""

from __future__ import division

__author__ = "Will Richards, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Will Richards"
__email__ = "wrichards@mit.edu"
__date__ = "Sep 25, 2012"


from pymatgen.core.periodic_table import smart_element_or_specie
from pymatgen.serializers.json_coders import MSONable
import abc


class AbstractStructureFilter(MSONable):
    """
    Abstract structure filter class.
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @abc.abstractmethod
    def test(self, structure):
        '''
        Returns a boolean for any structure. Structures that return true are
        kept in the Transmuter object during filtering.
        '''
        return

    @staticmethod
    def from_dict(d):
        for trans_modules in ['filters']:
            mod = __import__('pymatgen.alchemy.' + trans_modules,
                             globals(), locals(), [d['@class']], -1)
            if hasattr(mod, d['@class']):
                trans = getattr(mod, d['@class'])
                return trans(**d['init_args'])
        raise ValueError("Invalid filter dict")


class ContainsSpecieFilter(AbstractStructureFilter):
    """
    Filter for structures containing certain elements or species.
    By default compares by atomic number
    """

    def __init__(self, species, strict_compare=False, AND=True, exclude=False):
        """
        Args:
            species:
                list of species to look for
            AND:
                whether all species must be present to pass (or fail)
                filter.
            strict_compare:
                if true, compares objects by specie or element object
                if false, compares atomic number
            exclude:
                if true, returns false for any structures with the specie
                (excludes them from the Transmuter)
        """
        self._species = map(smart_element_or_specie, species)
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
        output = ["ContainsSpecieFilter with parameters:"]
        output.append("species = {}".format(self._species))
        output.append("strict_compare = {}".format(self._strict))
        output.append("AND = {}".format(self._AND))
        output.append("exclude = {}".format(self._exclude))
        return  "\n".join(output)

    @property
    def to_dict(self):
        d = {"version": __version__}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["init_args"] = {"species": [str(sp) for sp in self._species],
                          "strict_compare": self._strict,
                          "AND": self._AND,
                          "exclude": self._exclude}
        return d


class SpecieProximityFilter(AbstractStructureFilter):
    """
    This filter removes structures that have certain species that are too close
    together.
    """

    def __init__(self, specie_and_min_dist_dict):
        """
        Args:
            specie_and_min_dist_dict:
                A species string to float mapping. For example, {"Na+": 1}
                means that all Na+ ions must be at least 1 Angstrom away from
                each other. Multiple species criteria can be applied. Note that
                the testing is done based on the actual object. If you have a
                structure with Element, you must use {"Na":1} instead to filter
                based on Element and not Specie.
        """
        self.specie_and_min_dist = {smart_element_or_specie(k): v
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

    @property
    def to_dict(self):
        d = {"version": __version__}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["init_args"] = {"specie_and_min_dist_dict":
                          {str(sp): v
                           for sp, v in self.specie_and_min_dist.items()}}
        return d
