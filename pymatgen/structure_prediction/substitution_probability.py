from pymatgen.core.periodic_table import Specie
from datetime import datetime

import itertools
import json
import math
import os
import pymatgen

class SubstitutionProbability():
    """
    This class finds substitution probabilities given lists of atoms
    to substitute.
    The inputs make more sense if you look through the from_defaults 
    static method
    Args:
        lambda_dict:
            dictionary of the weight functions lambda (theta in BURP)
        alpha:
            weight function for never observed substitutions
        Z:
            partition function (constant)
        px:
            dictionary of partition function for individual species
    """
    def __init__(self, lambda_dict, alpha, Z, px):
        self._lambda = lambda_dict
        self._alpha = alpha
        self._Z = Z
        self._px = px
        self.species_list = list(px.keys())
        
    def prob(self, s1, s2):
        """
        get the probability of 2 species substitution
        Returns:
            probability of s1 and s2 substitution
        Not used by the structure predictor. Replicates some of BURP
        object functionality
        """
        l = self._lambda.get(frozenset([s1, s2]), self._alpha)
        return math.exp(l)/self._Z
    
    def cond_prob(self, s1, s2):
        """
        Args:
            s1: The VARIABLE specie
            s2: The FIXED specie
        this is the conditional probability used by burp structure 
        predictor (which is why its in a wierd order)
        """
        l = self._lambda.get(frozenset([s1, s2]), self._alpha)
        return math.exp(l)/self._px[s2]
    
    def pair_corr(self, s1, s2):
        """
        returns the pair correlation of 2 species
        """
        l = self._lambda.get(frozenset([s1, s2]), self._alpha)
        return math.exp(l)*self._Z/(self._px[s1] * self._px[s2])
        
    def cond_prob_list(self, l1, l2):
        """
        Find the probabilities of 2 lists. These should include ALL 
        species. This is the probability conditional on l2
        Args:
            l1, l2:
                lists of species
        Returns:
            the conditional probability (assuming these species are in
            l2)
            
        """
        assert len(l1) == len(l2)
        p = 1.
        for i, s1 in enumerate(l1):
            s2 = l2[i]
            p *= self.cond_prob(s1, s2)
        return p
    
    @staticmethod
    def from_defaults(alpha = 1e-4):
        #Something seems to be weird about the SubsProbaMRF file:
        #alpha is VERY high. alpha essentially replaces the lambda 
        #value for never-observed compounds, but is positive even 
        #though many (most?) are negative numbers
        module_dir = os.path.dirname(pymatgen.__file__)
        json_file = os.path.join(module_dir, 'structure_prediction'
                                 , 'data', 'lambda.json')
        
        with open(json_file) as f:
            table = json.load(f)
        
        #build map of specie pairs to lambdas
        #build set of species
        l = {}
        sp_set = set()
        for row in table:
            if not row[0] == 'D1+' and not row[1] == 'D1+':
                s1 = Specie.from_string(row[0])
                s2 = Specie.from_string(row[1])
                l[frozenset([s1, s2])] = float(row[2])
                sp_set.add(s1)
                sp_set.add(s2) 

        #calculate Z and the individual species partition function 
        #(px)
        px = dict.fromkeys(sp_set, 0.)
        Z=0
        for s1, s2 in itertools.product(sp_set,repeat = 2):
            value = math.exp(l.get(frozenset([s1, s2]), alpha))
            #not sure why the factor of 2 is here but it matches up 
            #with BURP. BURP may actually be missing a factor of 2, 
            #but it doesn't have a huge effect
            px[s1] += value/2 
            px[s2] += value/2 
            Z += value

        return SubstitutionProbability(l, alpha, Z, px)
    
    
    