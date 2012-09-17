"""
Defines filters for Transmuter object
"""

from pymatgen.core.periodic_table import smart_element_or_specie

import abc

class AbstractStructureFilter(object):
    """
    Abstract structure filter class.
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self):
        pass
    
    @abc.abstractmethod
    def test(self, structure):
        '''
        Returns a boolean for any structure. Structures that 
        return true are kept in the Transmuter object during 
        filtering
        '''
        return


class ContainsSpecieFilter(AbstractStructureFilter):
    """
    Filter for structures containing certain elements or species.
    By default compares by atomic number
    """
    def __init__(self, species, strict_compare = False
                 , AND = True
                 , exclude = False):
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
                if true, returns false for any structures with 
                the specie (excludes them from the Transmuter)
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
            structure_set = set(map(atomic_number
                                    , structure.composition.elements))
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

    
