from __future__ import division

__author__ = "Will Richards"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.2"
__maintainer__ = "Will Richards"
__email__ = "wrichard@mit.edu"
__date__ = "Aug 31, 2012"

from pymatgen.core.periodic_table import Specie

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
        lambda_table:
            json table of the weight functions lambda (theta 
            in BURP)
            if None, will use the default lambda.json table
        alpha:
            weight function for never observed substitutions
    """
    def __init__(self, lambda_table = None, alpha = 1e-4):
        #Something seems to be weird about the SubsProbaMRF file:
        #alpha is VERY high. alpha essentially replaces the lambda 
        #value for never-observed compounds, but is positive even 
        #though many (most?) are negative numbers
        
        #store the input table for the to_dict method
        self._lambda_table = lambda_table
        
        if not lambda_table:
            module_dir = os.path.dirname(pymatgen.__file__)
            json_file = os.path.join(module_dir, 'structure_prediction'
                                     , 'data', 'lambda.json')
            with open(json_file) as f:
                lambda_table = json.load(f)
            
        #build map of specie pairs to lambdas
        l = {}
        for row in lambda_table:
            if not row[0] == 'D1+' and not row[1] == 'D1+':
                s1 = Specie.from_string(row[0])
                s2 = Specie.from_string(row[1])
                l[frozenset([s1, s2])] = float(row[2])
        
        self._lambda = l 
        self._alpha = alpha
        self.update_normalizations()
         
    def update_normalizations(self):
        """
        updates the partition functions Z and px,
        and the species list
        """
        sp_set = set()
        for key in self._lambda.keys():
            sp_set.update(key)
        px = dict.fromkeys(sp_set, 0.)
        Z=0
        for s1, s2 in itertools.product(sp_set,repeat = 2):
            value = math.exp(self._lambda.get(frozenset([s1, s2])
                                              , self._alpha))
            #not sure why the factor of 2 is here but it matches up 
            #with BURP. BURP may actually be missing a factor of 2, 
            #but it doesn't have a huge effect
            px[s1] += value/2 
            px[s2] += value/2 
            Z += value
        
        self._Z = Z
        self._px = px
        self.species_list = list(sp_set)
        
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
    
    @property
    def to_dict(self):
        d = {"name": self.__class__.__name__, "version": __version__}
        d["init_args"] = {"lambda_table": self._lambda_table,
                          "alpha": self._alpha}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d
    
    @staticmethod
    def from_dict(d):
        return SubstitutionProbability(**d['init_args'])
    
    @staticmethod
    def _with_test_table(alpha = 1e-4):
        """
        loads a lightweight lambda table for use in unit tests to 
        reduce initialization time, and make unit tests insensitive
        to changes in the default lambda table
        """
        module_dir = os.path.dirname(pymatgen.__file__)
        json_file = os.path.join(module_dir, 'structure_prediction'
                                 , 'tests', 'test_data'
                                 , 'test_lambda.json')
        with open(json_file) as f:
            lambda_table = json.load(f)
            
        return SubstitutionProbability(lambda_table, alpha)
    