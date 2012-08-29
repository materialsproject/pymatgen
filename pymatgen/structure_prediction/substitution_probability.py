from pymatgen.core.periodic_table import Specie

import json
import math
import os
import pymatgen

class SubstitutionProbability():
    """
    This class finds substitution probabilities given lists of
    atoms to substitute. needs A LOT of work (though it replicates burp fairly well)
    """
    def __init__(self, lambda_dict, alpha, Z):
        self._lambda = lambda_dict
        self._alpha = alpha
        self._Z = Z
        
    def probability(self, s1, s2):
        """
        get the probability of 2 species substitution
        Returns:
            probability of s1 and s2 substitution
        """
        try:
            l = self._lambda[frozenset([s1, s2])]
        except:
            l = self._alpha
        return math.exp(l)/self._Z
    
    def list_probability(self, l1, l2):
        """
        Find the probabilities of 2 lists. These should include ALL species
        though I have NO idea why. Really any unchanged species shouldn't affect
        the probability!!!! (eg any substitution on a structure that contains Zr3+
        will have a lower probability than one that contains Zr4+, even if the 
        substitution is the same!!!!!!) WTF?
        eg: try the substitution Na1+ Cl1- to Na1+ Br1-
            and the substitution Ag2+ Cl1- to Ag2+ Cl1-
            They get VERY different results. 
        Args:
            l1, l2:
                lists of species
        Returns:
            the probability
            
        """
        assert len(l1) == len(l2)
        p = 1.
        for i, s1 in enumerate(l1):
            s2 = l2[i]
            p *= self.probability(s1, s2)
        return p
    
    @staticmethod
    def from_defaults(alpha = 1e-4):
        module_dir = os.path.dirname(pymatgen.__file__)
        json_file = os.path.join(module_dir, 'structure_prediction'
                                 , 'data', 'lambda.json')
        
        with open(json_file) as f:
            table = json.load(f)
            
        l = {}
        sp_set = set()
        Z = 0
        for row in table:
            if not row[0] == 'D1+' and not row[1] == 'D1+':
                s1 = Specie.from_string(row[0])
                s2 = Specie.from_string(row[1])
                sp_set.add(s1)
                sp_set.add(s2)
                l[frozenset([s1, s2])] = float(row[2])
                if s1 == s2:
                    Z += math.exp(float(row[2]))
                else:
                    Z += 2*math.exp(float(row[2]))
        
        nbZeros = len(sp_set)**2+len(sp_set)-2*len(l)

        Z += nbZeros*math.exp(alpha)
        
        #Something that seems to be weird about the SubsProbaMRF file:
        #alpha is VERY high. alpha essentially replaces the lambda value for never-observed
        #compounds, but is positive even though many (most?) are negative numbers
        
        return SubstitutionProbability(l, alpha, Z)
    
    
    