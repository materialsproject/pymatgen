#!/usr/bin/env python

"""
This module provides classes that define a chemical reaction.
"""

from __future__ import division

__author__="Shyue Ping Ong, Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import logging
import itertools
import numpy as np

from pymatgen.core.structure import Composition

logger = logging.getLogger(__name__)

class Reaction(object):
    """
    A class representing a Reaction.
    """
    
    TOLERANCE = 1e-6 #: Tolerance for determining if a particular component fraction is > 0.
    
    def __init__(self, reactants, products):
        """
        Reactants and products to be specified as list of pymatgen.core.structure.Composition.  
        e.g., [comp1, comp2]
        
        Args:
            reactants : List of reactants.
            products : List of products.
        """
        all_comp = reactants[:]
        all_comp.extend(products[:])
        els = set()
        for c in all_comp:
            els.update(c.elements)
        els = tuple(els)
        
        num_constraints = len(all_comp)
        num_els = len(els)
        dim = max(num_els, num_constraints)
        logger.debug('num_els = {}'.format(num_els))
        logger.debug('num_constraints = {}'.format(num_constraints))
        logger.debug('dim = {}'.format(dim))
        
        if num_constraints < 2:
            raise ReactionError("A reaction cannot be formed with just one composition.")
        elif num_constraints == 2:
            if all_comp[0].reduced_formula != all_comp[1].reduced_formula:
                raise ReactionError("%s and %s cannot be matched.  Reaction cannot be balanced. " % (all_comp[0].formula, all_comp[1].formula) )
            else:
                coeffs = [-all_comp[1][els[0]]/all_comp[0][els[0]], 1]
        else:
            comp_matrix = np.zeros((dim, dim))
            count = 0
            if num_constraints < num_els:
                for i in range(num_constraints, num_els):
                    all_comp.append(Composition({els[i]:1}))
            for c in all_comp:
                for i in range(num_els):
                    comp_matrix[i][count] = c[els[i]]
                count += 1
            
            if num_constraints > num_els:
                
                #Try two schemes for making the comp matrix non-singular.
                for i in range(num_els, num_constraints):
                    for j in range(num_els):
                        comp_matrix[i][j] = 0
                    comp_matrix[i][i] = 1
                count = 0
                if abs(np.linalg.det(comp_matrix)) < self.TOLERANCE:
                    for i in range(num_els, num_constraints):
                        for j in range(num_els):
                            comp_matrix[i][j] = count
                            count += 1
                        comp_matrix[i][i] = count
                
                ans_matrix = np.zeros(num_constraints)
                ans_matrix[num_els:num_constraints] = 1
                coeffs = np.linalg.solve(comp_matrix, ans_matrix)
            elif num_constraints <= num_els:
                if abs(np.linalg.det(comp_matrix)) < self.TOLERANCE:
                    logger.debug('Linear solution possible. Trying various permutations.')
                    comp_matrix = comp_matrix[0:num_els][:,0:num_constraints]
                    logger.debug('comp_matrix = {}'.format(comp_matrix))
                    ans_found = False
                    for perm_matrix in itertools.permutations(comp_matrix):
                        logger.debug('Testing permuted matrix = {}'.format(perm_matrix))
                        for m in xrange(num_constraints):
                            submatrix = [[perm_matrix[i][j] for j in xrange(num_constraints) if j != m] for i in xrange(num_constraints) if i != m]
                            logger.debug('Testing submatrix = {}'.format(submatrix))
                            if abs(np.linalg.det(submatrix)) > self.TOLERANCE:
                                logger.debug('Possible sol')
                                subansmatrix = [perm_matrix[i][m] for i in xrange(num_constraints) if i != m]
                                coeffs = - np.linalg.solve(submatrix, subansmatrix)
                                coeffs = [c for c in coeffs]
                                coeffs.insert(m,1)
                                #Check if final coeffs are valid
                                overall_mat = np.dot(perm_matrix, coeffs)
                                if (abs(overall_mat) < 1e-8).all():
                                    ans_found = True
                                    break
                    if not ans_found:
                        raise ReactionError("Reaction is ill-formed and cannot be balanced.")
                else:
                    raise ReactionError("Reaction is ill-formed and cannot be balanced.")
            else:
                raise ReactionError("Reaction is ill-formed and cannot be balanced.")
        for i in xrange(len(coeffs)-1,-1,-1):
            if coeffs[i] != 0:
                normfactor = coeffs[i]
                break
        #Invert negative solutions and scale to final product
        coeffs = [c/normfactor for c in coeffs]
        self._els = els
        self._all_comp = all_comp
        self._coeffs = coeffs
        self._num_comp = num_constraints
    
    def copy(self):
        """
        Returns a copy of the Reaction object.
        """
        return Reaction(self.reactants, self.products)
    
    def calculate_energy(self, energies):
        """
        Calculates the energy of the reaction.
         
        Args:
            energies - dict of {comp:energy}.  e.g., {comp1: energy1, comp2: energy2}.
        
        Returns:
            reaction energy as a float.
        """
        return sum([self._coeffs[i] * energies[self._all_comp[i]] for i in range(self._num_comp)])
    
    def normalize_to(self, comp, factor = 1):
        """
        Normalizes the reaction to one of the compositions.
        By default, normalizes such that the composition given has a coefficient of 1.
        Another factor can be specified.
        """
        scale_factor = abs(1/self._coeffs[self._all_comp.index(comp)] * factor)
        self._coeffs = [c * scale_factor for c in self._coeffs]
    
    def normalize_to_element(self, element, target_amount = 1):
        """
        Normalizes the reaction to one of the elements.
        By default, normalizes such that the amount of the element is 1.
        Another factor can be specified.
        """
        current_element_amount = sum([self._all_comp[i][element] * abs(self._coeffs[i]) for i in xrange(len(self._all_comp))]) / 2
        scale_factor = target_amount / current_element_amount
        self._coeffs = [c * scale_factor for c in self._coeffs]
    
    @property
    def elements(self):
        """
        List of elements in the reaction
        """
        return self._els[:]
    
    @property
    def coeffs(self):
        """
        Final coefficients of the calculated reaction
        """
        return self._coeffs[:]
    
    @property
    def all_comp(self):
        """
        List of all compositions in the reaction.
        """
        return self._all_comp
    
    @property
    def reactants(self):
        """List of reactants"""
        
        return [self._all_comp[i] for i in xrange(len(self._all_comp)) if self._coeffs[i] < 0]
    
    @property
    def products(self):
        """List of products"""
        return [self._all_comp[i] for i in xrange(len(self._all_comp)) if self._coeffs[i] > 0]
    
    def get_coeff(self, comp):
        """Returns coefficient for a particular composition"""
        return self._coeffs[self._all_comp.index(comp)]
    
    def normalized_repr_and_factor(self):
        """
        Normalized representation for a reaction
        For example, ``4 Li + 2 O -> 2Li2O`` becomes ``2 Li + O -> Li2O``
        """
        reactant_str= []
        product_str = []
        scaled_coeffs = []
        reduced_formulas = []
        for i in range(self._num_comp):
            comp = self._all_comp[i]
            coeff = self._coeffs[i]
            (reduced_formula, scale_factor) = comp.get_reduced_formula_and_factor()
            scaled_coeffs.append(coeff * scale_factor)
            reduced_formulas.append(reduced_formula)
        
        count = 0
        while sum([abs(coeff) % 1 for coeff in scaled_coeffs]) > 1e-8:
            norm_factor = 1/smart_float_gcd(scaled_coeffs)
            scaled_coeffs = [c / norm_factor for c in scaled_coeffs]
            count += 1
            if count > 10: #Prevent an infinite loop
                break
            
        for i in range(self._num_comp):  
            if scaled_coeffs[i] == -1:
                reactant_str.append(reduced_formulas[i])
            elif scaled_coeffs[i] == 1:
                product_str.append(reduced_formulas[i])
            elif scaled_coeffs[i] < 0:
                reactant_str.append("%d %s" % (-scaled_coeffs[i], reduced_formulas[i]))
            elif scaled_coeffs[i] > 0:
                product_str.append("%d %s" % (scaled_coeffs[i], reduced_formulas[i]))
        factor = scaled_coeffs[0] / self._coeffs[0]
        
        return (" + ".join(reactant_str) + " -> " + " + ".join(product_str), factor)
    
    @property
    def normalized_repr(self):
        return self.normalized_repr_and_factor()[0]
    
    def __repr__(self):
        return self.__str__()
    
    def __str__(self):
        reactant_str= []
        product_str = []
        for i in range(self._num_comp):
            comp = self._all_comp[i]
            coeff = self._coeffs[i]
            red_comp = Composition.from_formula(comp.reduced_formula)
            scale_factor = comp.num_atoms / red_comp.num_atoms
            scaled_coeff = coeff * scale_factor
            if scaled_coeff < 0:
                reactant_str.append("%.3f %s" % (-scaled_coeff, comp.reduced_formula))
            elif scaled_coeff > 0:
                product_str.append("%.3f %s" % (scaled_coeff, comp.reduced_formula))
        
        return " + ".join(reactant_str) + " -> " + " + ".join(product_str)

        
def smart_float_gcd(list_of_floats):
    """
    Determines the great common denominator (gcd).  Works on floats as well as integers.
    
    Args:
        list_of_floats: List of floats to determine gcd.
    """
    mult_factor = 1.0
    all_remainders = sorted([abs(f - int(f)) for f in list_of_floats])
    for i in range(len(all_remainders)):
        if all_remainders[i] > 1e-5:
            mult_factor *= all_remainders[i]
            all_remainders = [ f2 / all_remainders[i] - int(f2 / all_remainders[i]) for f2 in all_remainders]
    return 1/mult_factor
    

class ReactionError(Exception):
    '''
    Exception class for Reactions. Allows more information exception messages to cover situations not
    covered by standard exception classes.
    '''

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "Query Error : " + self.msg
