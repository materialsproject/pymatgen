#!/usr/bin/env python

'''
This module defines standard transformations which transforms a structure into another structure.
All transformations should inherit the AbstractTransformation ABC.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import json
import itertools
import warnings
import numpy as np
import copy

from pymatgen.transformations.transformation_abc import AbstractTransformation
from pymatgen.core.structure import Structure
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor, SupercellMaker
from pymatgen.core.periodic_table import smart_element_or_specie
from pymatgen.analysis.ewald import EwaldSummation

class IdentityTransformation(AbstractTransformation):
    """
    This is a demo transformation which does nothing, i.e. just return the same structure.
    """
    
    def __init__(self):
        pass
    
    def apply_transformation(self, structure):
        return Structure(structure.lattice, structure.species_and_occu, structure.frac_coords)
    
    def __str__(self):
        return "Identity Transformation"
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return self
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__, 'init_args': {} }
        return output
      
class RotationTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """
    
    def __init__(self, axis, angle, angle_in_radians = False):
        """
        Arguments:
            axis - Axis of rotation, e.g., [1, 0, 0]
            angle - angle to rotate
            angle_in_radians - Set to True if angle is supplied in radians. Else degrees are assumed.
        """
        self._axis = axis
        self._angle = angle
        self._angle_in_radians = angle_in_radians
        self._symmop = SymmOp.from_axis_angle_and_translation(self._axis, self._angle, self._angle_in_radians)
    
    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        editor.apply_operation(self._symmop)
        return editor.modified_structure
    
    def __str__(self):
        return "Rotation Transformation about axis %s with angle = %.4f %s" % (str(self._axis), self._angle, "radians" if self._angle_in_radians else "degrees")
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return RotationTransformation(self._axis, - self._angle, self._angle_in_radians)
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'axis': self._axis, 'angle':self._angle, 'angle_in_radians':self._angle_in_radians}
        return output
        
class SupercellTransformation(AbstractTransformation):
    """
    The RotationTransformation applies a rotation to a structure.
    """
    
    def __init__(self, scaling_matrix = [[1,0,0],[0,1,0],[0,0,1]]):
        """
        Arguments:
            scaling_matrix - Set to True if angle is supplied in radians. Else degrees are assumed.
        """
        self._matrix = scaling_matrix
        
    def apply_transformation(self, structure):
        maker = SupercellMaker(structure, self._matrix)
        return maker.modified_structure
    
    def __str__(self):
        return "Supercell Transformation with scaling matrix %s" % (str(self._matrix))
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        raise NotImplementedError()
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'scaling_matrix': self._matrix}
        return output
        
    
class SubstitutionTransformation(AbstractTransformation):
    """
    This transformation substitutes species for one another.
    """
    def __init__(self, species_map):
        """
        Arguments:
            species_map - a dict containing the species mapping in string-string pairs. E.g., { "Li":"Na"} or {"Fe2+","Mn2+"}. Multiple substitutions can be done.
        """
        self._species_map = species_map
    
    def apply_transformation(self, structure):
        species_map = {smart_element_or_specie(k): smart_element_or_specie(v) for k, v in self._species_map.items()}
        editor = StructureEditor(structure)
        editor.replace_species(species_map)
        return editor.modified_structure
    
    def __str__(self):
        return "Substitution Transformation :" + ", ".join([k+"->"+v for k, v in self._species_map.items()])
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return SubstitutionTransformation({v:k for k, v in self._species_map.items()})
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'species_map': self._species_map}
        return output
        
     
class RemoveSpeciesTransformation(AbstractTransformation):
    """
    Remove all occurrences of some species from a structure.
    """
    def __init__(self, species_to_remove):
        """
        Arguments:
            species_to_remove - List of species to remove. E.g., ["Li", "Mn"] 
        """
        self._species = species_to_remove
    
    def apply_transformation(self, structure):
        editor = StructureEditor(structure)
        map(editor.remove_species, [[smart_element_or_specie(sp)] for sp in self._species])
        return editor.modified_structure
    
    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._species)
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'species_to_remove': self._species}
        return output
        
class PartialRemoveSpecieTransformation(AbstractTransformation):
    """
    Remove fraction of specie from a structure. Requires an oxidation state decorated structure for ewald sum to be computed.
    """
    def __init__(self, specie_to_remove, fraction_to_remove):
        """
        Arguments:
            specie_to_remove - Specie to remove. Must have oxidation state E.g., "Li1+"
            fraction_to_remove - Fraction of specie to remove. E.g., 0.5
        """
        self._specie = specie_to_remove
        self._frac = fraction_to_remove
    
    def apply_transformation(self, structure):
        sp = smart_element_or_specie(self._specie)
        num_to_remove = structure.composition[sp] * self._frac
        if abs(num_to_remove - int(num_to_remove)) > 1e-8:
            raise ValueError("Fraction to remove must be consistent with integer amounts in structure.")
        else:
            num_to_remove = int(round(num_to_remove))
        
        specie_indices = [i for i in xrange(len(structure)) if structure[i].specie == sp]
                
        lowestewald = float('inf')
        lowestenergy_s = None
        all_structures = list()
        for indices in itertools.combinations(specie_indices, num_to_remove):
            mod = StructureEditor(structure)
            mod.delete_sites(indices)
            s_new = mod.modified_structure.get_sorted_structure()
            all_structures.append(s_new)
            ewaldsum = EwaldSummation(s_new)
            if ewaldsum.total_energy < lowestewald:
                lowestewald = ewaldsum.total_energy
                lowestenergy_s = s_new
        
        self.all_structures = all_structures        
        
        return lowestenergy_s
    
    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._specie)
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'specie_to_remove': self._specie, 'fraction_to_remove': self._frac}
        return output
    
class PartialRemoveSpecieTransformation2(AbstractTransformation):
    """
    Remove fraction of specie from a structure. Requires an oxidation state decorated structure for ewald sum to be computed.
    This method uses the matrix form of ewaldsum to calculate the ewald sums of the potential structures.
    This is on the order of 4 orders of magnitude faster when there are large numbers of permutations to consider.
    There are further optimizations possible (doing a smarter search of permutations for example), but this wont make a difference
    until the number of permutations is on the order of 30,000.
    """
    def __init__(self, specie_to_remove, fraction_to_remove):
        """
        Arguments:
            specie_to_remove - Specie to remove. Must have oxidation state E.g., "Li1+"
            fraction_to_remove - Fraction of specie to remove. E.g., 0.5
        """
        self._specie = specie_to_remove
        self._frac = fraction_to_remove
        
    def _minimize_matrix(self, matrix, indices, num_to_remove, CURRENT_MINIMUM): #for whatever reason a default value wont work for CURRENT_MINIMUM
                                                                            #must pass [float('inf')] VERY IMPORTANT
                                                                            #current minimum works like a global variable since it is a list
        '''minimize a matrix by removing a specific number of rows and columns (if row 4 is removed, column 4 must also be removed)
        This method looks for short circuits to a brute force search by looking at best and worse case scenarios (which can be computed quickly)
        It is about 1000 times faster than brute force'''
        
        if num_to_remove>len(indices):  #if we've kept too many rows
            return[float('inf'), []]    #abandon the branch
        
        if num_to_remove==0:                        #if we don't have to remove any more rows
            matrix_sum = sum(sum(matrix))
            if matrix_sum<CURRENT_MINIMUM[0]:
                CURRENT_MINIMUM[0] = matrix_sum
            return [matrix_sum, []]                 #return the sum of the matrix
        
        indices = list(indices)         #make a copy of the indices so recursion doesn't alter them
        matrix2=copy.deepcopy(matrix)
        
        max_index = None
        
        #compute the best case sum for removing rows
        if num_to_remove>1: #lets not do this if we're finding the minimum in the next step anyway
            index_sum = [sum(matrix[i,:])+sum(matrix[:,i])-matrix[i,i] for i in indices]    #compute the value associated with an index assuming no other indices are removed
            max_index = indices[index_sum.index(max(index_sum))]                            #get the index of the maximum value (we'll try removing this row first)
            index_sum_sorted = list(index_sum)
            index_sum_sorted.sort(reverse = True)
            all_interactions = list(matrix[indices,:][:,indices].flatten())                 #get the cells that we could be double counting when removing multiple index
            all_interactions.sort()
            
            #sum the matrix - the rows with maximum 
            best_case = sum(sum(matrix))-sum(index_sum_sorted[:num_to_remove])+sum(all_interactions[:(num_to_remove*(num_to_remove-1))])
            
            if best_case > CURRENT_MINIMUM[0]:
                return [float('inf'), []]   #if the best case doesn't beat the minimum abandon the branch
            
            #try to find rows that should definitely be removed or definitely ignored based on best and worse case performances
            most_positive = []
            most_negative = []
            for i in range(len(indices)):
                index = indices[i]
                interactions = [matrix[index,x]+matrix[x,index] for x in indices if not x == index]
                interactions.sort()
                most_positive.append(index_sum[i]-sum(interactions[:(num_to_remove-1)]))
                most_negative.append(index_sum[i]-sum(interactions[-(num_to_remove-1):]))
                
            most_positive_sorted = sorted(most_positive, reverse = True)
            most_negative_sorted = sorted(most_negative, reverse = True)    
            
            deletion_indices = []
            ignore_indices = []
            for i in range(len(indices)):
                if most_negative[i] > most_positive_sorted[num_to_remove-1]:
                    deletion_indices.append(indices[i])
                    pass
                if most_positive[i] < most_negative_sorted[num_to_remove-1]:
                    ignore_indices.append(indices[i])
                    pass
                    
            if deletion_indices + ignore_indices:
                for r_index in deletion_indices:
                    matrix2[:,r_index] = 0
                    matrix2[r_index,:] = 0
                    num_to_remove -= 1
                    indices.remove(r_index)
                for x in ignore_indices:
                    indices.remove(x)
                output = self._minimize_matrix(matrix2, indices, num_to_remove, CURRENT_MINIMUM)
                output[1] = output[1]+deletion_indices
                output[1].sort()
                return output
        
        #if no shortcuts could be found, recurse down one level by both removing and ignoring one row
        if max_index:
            r_index = max_index
            indices.remove(max_index)
        else:    
            r_index = indices.pop()
        matrix2[:,r_index] = 0
        matrix2[r_index,:] = 0
        sum2 = self._minimize_matrix(matrix2, indices, num_to_remove-1, CURRENT_MINIMUM)
        sum1 = self._minimize_matrix(matrix, indices, num_to_remove, CURRENT_MINIMUM)
        if sum1[0]<sum2[0]:
            return sum1
        else:
            sum2[1].append(r_index)
            sum2[1].sort()
            return sum2
    
    def apply_transformation(self, structure):
        sp = smart_element_or_specie(self._specie)
        num_to_remove = structure.composition[sp] * self._frac
        if abs(num_to_remove - int(num_to_remove)) > 1e-8:
            raise ValueError("Fraction to remove must be consistent with integer amounts in structure.")
        else:
            num_to_remove = int(round(num_to_remove))
        
        specie_indices = [i for i in xrange(len(structure)) if structure[i].specie == sp]
        
        ewaldmatrix = EwaldSummation(structure).total_energy_matrix
        
        lowestenergy_indices = self._minimize_matrix(ewaldmatrix, specie_indices, num_to_remove, [float('inf')])[1]
        
        mod = StructureEditor(structure)
        mod.delete_sites(lowestenergy_indices)
        return mod.modified_structure.get_sorted_structure()
    
    def __str__(self):
        return "Remove Species Transformation :" + ", ".join(self._specie)
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {'specie_to_remove': self._specie, 'fraction_to_remove': self._frac}
        return output

class OrderDisorderedStructureTransformation(AbstractTransformation):
    """
    Order a disordered structure. The disordered structure must be oxidation state decorated for ewald sum to be computed.
    Please note that the current form uses a "dumb" algorithm of completely enumerating all possible combinations of
    partially occupied sites.  No attempt is made to perform symmetry determination to reduce the number of combinations.
    Hence, attempting to performing ordering on a large number of disordered sites may be extremely expensive.  Also, simple
    rounding of the occupancies are performed, with no attempt made to achieve a target composition.  This is usually not a problem
    for most ordering problems, but there can be times where rounding errors may result in structures that do not have the desired composition.
    This second step will be implemented in the next iteration of the code.
    
    USE WITH CARE.
    """
    def __init__(self):
        pass
    
    def apply_transformation(self, structure, max_iterations = 100):
        """
        For this transformation, the apply_transformation method will return only the ordered
        structure with the lowest Ewald energy, to be consistent with the method signature of the other transformations.  
        However, all structures are stored in the  all_structures attribute in the transformation object for easy access.
        Arguments:
            structure:
                Oxidation state decorated disordered structure to order
            max_iterations:
                Maximum number of structures to consider.  Defaults to 100. This is useful if there are a large number of sites 
                and there are too many orderings to enumerate.
        """
        ordered_sites = []
        
        sites_to_order = {}
        for site in structure:
            species_and_occu = site.species_and_occu
            if sum(species_and_occu.values()) == 1 and len(species_and_occu) == 1:
                ordered_sites.append(site)
            else:
                spec = tuple([(sp, occu) for sp, occu in species_and_occu.items()])
                if spec not in sites_to_order:
                    sites_to_order[spec] = [site]
                else:
                    sites_to_order[spec].append(site)

        allselections = []
        species = []
        for spec, sites in sites_to_order.items():
            total_sites = len(sites)
            for (sp, fraction) in spec:
                num_to_select = int(round(fraction * total_sites))
                if num_to_select == 0:
                    raise ValueError("Fraction not consistent with selection of at least a single site.  Make a supercell before proceeding further.")
                allselections.append(itertools.combinations(sites, num_to_select))
                species.append(sp)
        
        all_ordered_s = {}
        count = 0
        
        def in_coords(allcoords, coord):
            for test_coord in allcoords:
                if all(coord == test_coord):
                    return True
            return False
        
        for selection in itertools.product(*allselections):
            all_species = [site.species_and_occu for site in ordered_sites]
            all_coords = [site.frac_coords for site in ordered_sites]
            
            contains_dupes = False
            for i in xrange(len(selection)):
                subsel = selection[i]
                sp = species[i]
                for site in subsel:
                    if not in_coords(all_coords, site.frac_coords):
                        all_species.append(sp)
                        all_coords.append(site.frac_coords)
                    else:
                        contains_dupes = True
                        break
                if contains_dupes:
                    break
            
            if not contains_dupes:
                s = Structure(structure.lattice, all_species, all_coords, False).get_sorted_structure()
                ewaldsum = EwaldSummation(s)
                ewald_energy = ewaldsum.total_energy
                all_ordered_s[s] = ewald_energy
                count += 1
                if count == max_iterations:
                    warnings.warn("Maximum number of iterations reached.  Structures will be ordered based on " + str(max_iterations) + " structures.")
                    break

        self.all_structures = all_ordered_s
        sorted_structures = sorted(all_ordered_s.keys(), key=lambda a: all_ordered_s[a])
                
        return sorted_structures[0]
    
    def __str__(self):
        return "Order disordered structure transformation"
     
    def __repr__(self):
        return self.__str__()
    
    @property
    def inverse(self):
        return None
    
    @property
    def to_dict(self):
        output = {'name' : self.__class__.__name__}
        output['init_args'] = {}
        return output


def transformation_from_json(json_string):
    """
    A helper function that can simply get a transformation from a json representation.
    
    Arguments:
        json_string:
            A json string representation of a transformation with init args.
    
    Returns:
        A properly initialized Transformation object
    """
    jsonobj = json.loads(json_string)
    trans = globals()[jsonobj['name']]
    return trans(**jsonobj['init_args'])


    