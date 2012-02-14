#!/usr/bin/env python

"""
This module provides classes for calculating the ewald sum of a structure.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = "Christopher Fischer"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import numpy as np
from math import pi, sqrt, log, exp, cos, sin, erfc
import scipy.constants as sc
from scipy.misc import comb
from pymatgen.core.structure import Structure
from copy import deepcopy, copy
import bisect


def compute_average_oxidation_state(site):
    """
    Calculates the average oxidation state of a site
    
    Args:
        Site to compute average oxidation state
        
    Returns:
        Average oxidation state of site.
    """
    return sum([sp.oxi_state * occu for sp, occu in site.species_and_occu.items() if sp != None])

def minimize_matrix( matrix, indices, num_to_remove, CURRENT_MINIMUM = None): #current minimum works like a global variable since it is a list
    '''minimize a matrix by removing a specific number of rows and columns (if row 4 is removed, column 4 must also be removed)
    This method looks for short circuits to a brute force search by looking at best and worse case scenarios (which can be computed quickly)
    It is about 1000 times faster than brute force'''
    
    if not CURRENT_MINIMUM:
        CURRENT_MINIMUM = [float('inf')]
    
    if num_to_remove>len(indices):  #if we've kept too many rows
        return[float('inf'), []]    #abandon the branch
    
    if num_to_remove==0:                        #if we don't have to remove any more rows
        matrix_sum = sum(sum(matrix))
        if matrix_sum<CURRENT_MINIMUM[0]:
            CURRENT_MINIMUM[0] = matrix_sum
        return [matrix_sum, []]                 #return the sum of the matrix
    
    indices = list(indices)         #make a copy of the indices so recursion doesn't alter them
    matrix2=deepcopy(matrix)
    
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
            output = minimize_matrix(matrix2, indices, num_to_remove, CURRENT_MINIMUM)
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
    sum2 = minimize_matrix(matrix2, indices, num_to_remove-1, CURRENT_MINIMUM)
    sum1 = minimize_matrix(matrix, indices, num_to_remove, CURRENT_MINIMUM)
    if sum1[0]<sum2[0]:
        return sum1
    else:
        sum2[1].append(r_index)
        sum2[1].sort()
        return sum2

class EwaldSummation:
    """
    Calculates the electrostatic energy of a periodic array of charges using the Ewald technique. 
    References : http://www.ee.duke.edu/~ayt/ewaldpaper/ewaldpaper.html
    
    This matrix can be used to do fast calculations of ewald sums after species removal
    
    E = E_recip + E_real + E_point
    
    Atomic units used in the code, then converted to eV.
    """

    # taken from convasp. converts unit of q*q/r into eV 
    CONV_FACT = 1e10 * sc.e / (4 * pi * sc.epsilon_0)
    
    def __init__(self, structure, real_space_cut = -1.0, recip_space_cut = -1.0, eta = -1.0, acc_factor = 8.0):
        """
        Initializes and calculates the Ewald sum. Default convergence parameters have been
        specified, but you can override them if you wish.
        
        Args:
            structure: input structure that must have proper Specie on all sites, i.e.
                       Element with oxidation state. Use OxidationStateDecorator in 
                       pymatgen.core.structure_modifier for example.
            real_space_cut: Real space cutoff radius dictating how many terms are used in the 
                            real space sum. negative means determine automagically using the
                            formula given in gulp 3.1 documentation.
            recip_space_cut: Reciprocal space cutoff radius. negative means determine automagically using
                             formula given in gulp 3.1 documentation.
            eta: The screening parameter. negative means determine automatically.
                 calculate_forces: Set to true if forces are desired
            acc_factor: No. of significant figures each sum is converged to. See the gulp manual. 
        """
        self._s = structure
        self._vol = structure.volume
        
        self._acc_factor = acc_factor

        # set screening length
        self._eta = eta if eta > 0 else (len(structure) * 0.01 / self._vol) ** (1/3) * pi        
        self._sqrt_eta = sqrt(self._eta)
        
        # acc factor used to automatically determine the optimal real and 
        # reciprocal space cutoff radii
        self._accf = sqrt(log(10 ** acc_factor))  
        
        self._rmax = real_space_cut if real_space_cut  > 0 else self._accf/self._sqrt_eta 
        self._gmax = recip_space_cut if recip_space_cut > 0 else 2* self._sqrt_eta * self._accf
        self._oxi_states = [compute_average_oxidation_state(site) for site in structure]
        self._total_q = sum(self._oxi_states)        
        self._forces = np.zeros((len(structure),3))
        self._calc_recip()
        self._calc_real_and_point()
        
    @property
    def reciprocal_space_energy(self):
        return sum(sum(self._recip))
    
    @property
    def reciprocal_space_energy_matrix(self):
        return self._recip
    
    @property
    def real_space_energy(self):
        return sum(sum(self._real))
    
    @property
    def real_space_energy_matrix(self):
        return self._real
    
    @property
    def point_energy(self):
        return sum(self._point)
    
    @property
    def point_energy_matrix(self):
        return self._point
    
    @property
    def total_energy(self):
        return sum(sum(self._recip)) + sum(sum(self._real)) + sum(self._point)
    
    @property
    def total_energy_matrix(self):
        totalenergy = self._recip + self._real
        for i in range(len(self._point)):
            totalenergy[i,i] += self._point[i]
        return totalenergy
    
    @property
    def forces(self):
        return self._forces
    
    def _calc_recip(self):
        """
        Perform the reciprocal space summation. Calculates the quantity
        E_recip = 1/(2PiV) sum_{G < Gmax} exp(-(G.G/4/eta))/(G.G) S(G)S(-G) where
        S(G) = sum_{k=1,N} q_k exp(-i G.r_k)
        S(G)S(-G) = |S(G)|**2
        """
        
        prefactor = 2 * pi / self._vol
        erecip = np.zeros( (self._s.num_sites,self._s.num_sites) )
                
        recip = Structure(self._s.lattice.reciprocal_lattice, ["H"], [np.array([0,0,0])])
        nn = recip.get_neighbors(recip[0], self._gmax)
        for (n,dist) in nn:
            gvect = n.coords
            gsquare = np.linalg.norm(gvect) ** 2
            expval = exp(-1.0 * gsquare / (4.0 * self._eta))
                        
            #calculate the structure factor
            sfactor = np.zeros( (self._s.num_sites,self._s.num_sites) )
            for i in range(self._s.num_sites):
                sitei = self._s[i]
                qi = self._oxi_states[i]
                for j in range(self._s.num_sites):
                    sitej = self._s[j]
                    qj = self._oxi_states[j]
                    exparg = np.dot(gvect, sitei.coords-sitej.coords) 
                    sfactor[i,j] = qi * qj * (cos(exparg) + sin(exparg))
                    
            erecip += expval/gsquare*sfactor
            
            #do forces if necessary
            sreal = 0.0
            simag = 0.0
            for i in range(self._s.num_sites):
                site = self._s[i]
                exparg = np.dot(gvect, site.coords)
                qj = self._oxi_states[i]  
                sreal += qj * cos(exparg)
                simag += qj * sin(exparg)
                
            for i in range(self._s.num_sites):    
                site = self._s[i]
                exparg = np.dot(gvect, site.coords) 
                qj = self._oxi_states[i]  
                pref = 2 * expval / gsquare * qj
                self._forces[i] += prefactor * pref * gvect * (sreal * sin(exparg) - simag * cos(exparg)) * EwaldSummation.CONV_FACT
            
        self._recip = erecip * prefactor * EwaldSummation.CONV_FACT 
    

    def _calc_real_and_point(self):
        """
         * determines the self energy
         * 
         * -(eta/pi)**(1/2) * sum_{i=1}^{N} q_i**2
         * 
         * if cell is charged a compensating background is added (i.e. a G=0 term)
         *  
         * @return
         */
        """
        
        all_nn = self._s.get_all_neighbors(self._rmax, True)
        forcepf = 2.0 * self._sqrt_eta/sqrt(pi)
        ereal = np.zeros( (self._s.num_sites,self._s.num_sites) )
        epoint = np.zeros( (self._s.num_sites) )
        
        for i in range(self._s.num_sites):
            site = self._s[i]
            nn = all_nn[i] #self._s.get_neighbors(site, self._rmax)
            qi = self._oxi_states[i]
            epoint[i] = qi*qi
            epoint[i] *= -1.0 * sqrt(self._eta/pi)
            epoint[i] += qi * pi / (2.0 * self._vol * self._eta)   #add jellium term
            for j in range(len(nn)):  #for (nsite, rij)  in nn:
                nsite = nn[j][0]
                rij = nn[j][1]
                qj = compute_average_oxidation_state(nn[j][0])  #nsite 
                erfcval = erfc(self._sqrt_eta*rij)
                ereal[nn[j][2],i] += erfcval * qi * qj / rij
                fijpf = qj / rij / rij / rij* (erfcval + forcepf * rij * exp(-self._eta * rij * rij)) 
                self._forces[i] += fijpf * (site.coords - nsite.coords) * qi * EwaldSummation.CONV_FACT
                    
        self._real = ereal * 0.5 * EwaldSummation.CONV_FACT
        self._point = epoint * EwaldSummation.CONV_FACT 
    
    @property
    def eta(self):
        return self._eta
    
    def __str__(self):
        output = ["Real = " + str(self.real_space_energy)]
        output.append("Reciprocal = " + str(self.reciprocal_space_energy))
        output.append("Point = " + str(self.point_energy))
        output.append("Total = " + str(self.total_energy))
        output.append("Forces:\n" + str(self.forces))
        return "\n".join(output)
    
class EwaldMinimizer:
    '''This class determines the manipulations that will minimize an ewald matrix, given a list of possible manipulations
    This class does not perform the manipulations on a structure, but will return the list of manipulations that should be
    done on one to produce the minimal structure. It returns the manipulations for the n lowest energy orderings. This
    class should be used to perform fractional species substitution or fractional species removal to produce a new structure.
    These manipulations create large numbers of candidate structures, and this class can be used to pick out those with the lowest ewald sum.
    An alternative (possibly more intuitive) interface to this class is the order disordered structure transformation.
    
    Author: Will Richards'''
    
    def __init__(self, matrix, m_list, num_to_return = 1, fast = True):
        '''arguments:
        matrix      a matrix of the ewald sum interaction energies. This is stored in the class as a diagonally symmetric array
                    and so self._matrix will not be the same as the input matrix
        m_list      list of manipulations
                    each item is of the form (multiplication fraction, number_of_indices, indices, species)
                    These are sorted such that the first manipulation contains the most permutations. this is actually
                    evaluated last in the recursion since I'm using pop
        num_to_return the minimizer will find the number_returned lowest energy structures. This is likely to return a number of
                    duplicate structures so it may be necessary to overestimate and then remove the duplicates later. (duplicate checking
                    in this process is extremely expensive)
        '''
        self._matrix = copy(matrix)
        for i in range(len(self._matrix)): #make the matrix diagonally symmetric (so matrix[i,:] == matrix[:,j])
            for j in range(i,len(self._matrix)):
                value = (self._matrix[i,j] + self._matrix[j,i])/2
                self._matrix[i,j] = value
                self._matrix[j,i] = value
        self._m_list = deepcopy(sorted(m_list, key = lambda x: comb(len(x[2]), x[1]), reverse = True)) #sort the m_list based on number of permutations
        self._current_minimum = float('inf')
        
        self._output_lists = []
        self._num_to_return = num_to_return
        
        n_combinations = 1  #calculate number of permutations to calculate
        for m in self._m_list:
            n_combinations *=comb(len(m[2]), m[1])
            
        if n_combinations <= num_to_return:
            fast = False                    #If we're returning all sums, don't bother with short circuit algorithm
        
        self.minimize_matrix(fast)
        
        self._best_m_list = self._output_lists[0][1]
        self._minimized_sum = self._output_lists[0][0]
        
    def add_m_list(self, matrix_sum, m_list):
        '''
        this adds an m_list to the output_lists and updates the current minimum if the list is full.
        '''
        if self._output_lists is None:
            self._output_lists = [[matrix_sum, m_list]]
        else:
            bisect.insort(self._output_lists, [matrix_sum, m_list])
        if len(self._output_lists)>self._num_to_return:
            self._current_minimum.pop()
        if len(self._output_lists) == self._num_to_return:
            self._current_minimum = self._output_lists[-1][0]
        
        
    def best_case(self, matrix, manipulation, indices_left):
        '''
        computes a best case given a matrix and manipulation list. This is only used for when there is only one manipulation left
        calculating a best case when there are multiple fractions remaining is much more complex (as sorting and dot products have to be done
        on each row and it just generally scales pretty badly.
        
        arguments:
        matrix: the current matrix (with some permutations already performed
        manipulation: (multiplication fraction, number_of_indices, indices, species) describing the manipulation
        indices: set of indices which haven't had a permutation performed on them 
        '''
        indices = list(indices_left.intersection(manipulation[2])) #only look at the indices that are in the manipulation and haven't been used
        fraction = manipulation[0]
        num = manipulation[1]
        
        sums = [2*sum(matrix[i,:])-(1-fraction)*matrix[i,i] for i in indices]   #compute the sum of each row. the weird prefactor on the self 
                                                                                #interaction term is to make the math work out because it actually 
                                                                                #is multiplied by f**2vrather than f like all the other terms
        
        if fraction < 1:    #we want do the thing that we guess will be most minimizing (to get to a lower current minimum faster). Whether 
                            #we want to use the most positive or negative row depends on whether we're increasing or decreasing that row
            next_index = indices[sums.index(max(sums))]
        else:
            next_index = indices[sums.index(min(sums))]

        interactions = list(matrix[indices,:][:,indices].flatten()) #compute a list of things that we may have double counted
        
        if fraction <=1:    #we use the first part of sorted sums. Whether we want highest or lowest depends on the fraction
            sums.sort(reverse = True)
        else:    
            sums.sort()
        if 2*fraction-fraction**2 <= 1:
            interactions.sort(reverse = True)
        else:
            interactions.sort()
        
        #compute a best case using the most minimizing row sums and the most minimizing interactions
        best_case = sum(sum(matrix))+(fraction-1)*sum(sums[:num])+(2*fraction-fraction**2-1)*sum(interactions[:num*(num-1)])
        
        return best_case, next_index
    
    def get_next_index(self, matrix, manipulation, indices_left): #returns an index that should have the most negative effect on the matrix sum
        f = manipulation[0]
        indices = list(indices_left.intersection(manipulation[2]))
        sums = [2*sum(matrix[i,:])-(1-f)*matrix[i,i] for i in indices]
        if f < 1:
            next_index = indices[sums.index(max(sums))]
        else:
            next_index = indices[sums.index(min(sums))]
        
        return next_index

    def _recurse(self, matrix, m_list, indices, fast, output_m_list = []):
        '''
        this method recursively finds the minimal permutations using a binary tree search strategy
        arguments:
        matrix: the current matrix (with some permutations already performed
        m_list: the list of permutations still to be performed
        indices: set of indices which haven't had a permutation performed on them
        fast: boolean indicating whether shortcuts by looking at the best case are to be used, the default is true
        

        returns:
        
        [minimal value, [list of replacements]]
        each replacement is a list [index of replaced specie, specie inserted at that index]
        '''
        
        if m_list[-1][1] == 0:
            m_list = copy(m_list)
            m_list.pop()
            if not m_list:
                matrix_sum = sum(sum(matrix))
                if matrix_sum < self._current_minimum:
                    self.add_m_list(matrix_sum, output_m_list)
                return 
            
        if m_list[-1][1]>len(indices.intersection(m_list[-1][2])):
            return
        
        index = None
        
        if fast and len(m_list)==1 and m_list[-1][1]>1:
            best_case = self.best_case(copy(matrix), m_list[-1], indices)
            index = best_case[1]
            if best_case[0] > self._current_minimum:
                return

        
        if index is None and fast:
            index = self.get_next_index(matrix, m_list[-1], indices)
        elif index is None:
            index = m_list[-1][2][-1]    
        
        m_list[-1][2].remove(index)
            
        matrix2 = copy(matrix)
        m_list2 = deepcopy(m_list)
        output_m_list2 = deepcopy(output_m_list)
        
        matrix2[index,:] *= m_list[-1][0]
        matrix2[:,index] *= m_list[-1][0]
        output_m_list2.append([index, m_list[-1][3]])
        indices2= copy(indices)
        m_list2[-1][1] -= 1
        
        if index in indices2:
            indices2.remove(index)
            output2 = self._recurse(matrix2, m_list2, indices2, fast, output_m_list2)
            
        output1 = self._recurse(matrix, m_list, indices, fast, output_m_list)
        
        
    def minimize_matrix(self, fast = True):
        '''
        this method finds and returns the permutations that produce the lowest ewald sum
        calls recursive function to iterate through permutations
        '''
        
        return self._recurse( self._matrix, self._m_list, set(range(len(self._matrix))), fast)
    
    @property
    def best_m_list(self):
        return self._best_m_list
    
    @property
    def minimized_sum(self):
        return self._minimized_sum
    
    @property
    def output_lists(self):
        return self._output_lists
        
    
    
    
   