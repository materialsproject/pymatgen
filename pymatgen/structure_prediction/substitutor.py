from __future__ import division

__author__ = "Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__date__ = "Aug 31, 2012"

from pymatgen.serializers.json_coders import MSONable
from pymatgen.structure_prediction.substitution_probability \
    import SubstitutionProbability
    
from operator import mul

class Substitutor(MSONable):
    def __init__(self, substitution_probability = None
                 , threshold = 1e-3):
        """
        This substitutor uses the substitution probability class to 
        find good substitutions for a given chemistry or structure. 
        Args:
            substitution_probability:
                a substitutionprobability object for finding the 
                probabilities
            threshold:
                probability threshold for predictions
        """
        if substitution_probability:
            self._sp = substitution_probability
        else:
            self._sp = SubstitutionProbability()
        self._threshold = threshold
        
    def pred_from_list(self, species_list):
        """
        Args:
            species_list:
                list of species in the starting structure
        Returns:
            list of species dictionaries (predictions)

        There are an exceptionally large number of substitutions to 
        look at (260^n), where n is the number of species in the 
        list. We need a more efficient than brute force way of going 
        through these possibilities. The brute force method would be:
        
        output = []
        for p in itertools.product(self._sp.species_list
                                   , repeat = len(species_list)):
            if self._sp.conditional_probability_list(p, species_list) 
                                   > self._threshold:
                output.append(dict(zip(species_list,p)))
        return output
        
        Instead of that we do a branch and bound
        """
        #calculate the highest probabilities to help us stop the 
        #recursion
        max_probabilities = []
        for s2 in species_list:
            max_p = 0
            for s1 in self._sp.species_list:
                max_p = max([self._sp.cond_prob(s1, s2),max_p])
            max_probabilities.append(max_p)
        output = []
        
        def _recurse(output_prob, output_species):
            best_case_prob = list(max_probabilities)
            best_case_prob[:len(output_prob)] = output_prob
            if reduce(mul, best_case_prob) > self._threshold:
                if len(output_species) == len(species_list):
                    output.append(dict(zip(species_list, output_species)))
                    return
                for sp in self._sp.species_list:
                    i = len(output_prob)
                    prob = self._sp.cond_prob(sp, species_list[i])
                    _recurse(output_prob+[prob], output_species + [sp])
                    
        _recurse([],[])
        return output
    
    def pred_from_comp(self, composition):
        """
        Similar to pred_from_list except this method returns a list 
        after checking that compositions are charge balanced
        """
        output = []
        predictions = self.pred_from_list(composition.elements)
        for p in predictions:
            charge = 0
            for i_el in composition.elements:
                f_el = p[i_el]
                charge += f_el.oxi_state * composition[i_el]
            if charge == 0:
                output.append(p)
        return output
    
    @property
    def to_dict(self):
        d = {"name": self.__class__.__name__, "version": __version__}
        d["substitution_probability"] = self._sp.to_dict
        d["threshold"] = self._threshold
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
    
    @staticmethod
    def from_dict(d):
        t = d['threshold']
        spd = d['substitution_probability']
        sp = SubstitutionProbability.from_dict(spd) 
        return Substitutor(sp, t)
        
    