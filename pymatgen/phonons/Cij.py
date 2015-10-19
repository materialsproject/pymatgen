import warnings, sys, operator, os, unittest
import pymatgen
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp import Vasprun
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import *
import copy
import numpy as np

__author__="Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Anubhav Jain, Mark Asta"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ ="March 20, 2012"


class CijFitter(object): 
    """
    Class that takes a dict in the form {IndependentStrain:stress matrix} and fits the elastic tensor from this.

    Args:
        - strain_stress_dict: dict containing stress matrices
    """
    
    def __init__(self, strain_stress_dict):
        self._sig_eps = strain_stress_dict

    def _chain_stresses(self, stress, ind1, ind2):
        s = []
        for el in stress:
            s.append(el[ind1,ind2])
        return s

    def fitCij(self, tol=0.1):
        stress_dict = self._sig_eps
        # Initialize Cij array
        Cij = np.zeros((6, 6))
        # Upper triangular indices
        # inds = zip(*np.triu_indices(3))
        inds = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]

        for n1 in range(0,6):
            strain = []
            stress = []
            count1 = 0

            for c in stress_dict:
                if c.i == inds[n1][0] and c.j== inds[n1][1]:
                    strain.append(c[c.i, c.j])
                    stress.append(stress_dict[c])

            for k in inds:
                true_data = self._chain_stresses(stress, k[0], k[1])
                p1 = np.polyfit(strain, true_data, 1)
                #p2 = np.polyfit(strain, true_data, 2)
                #p3 = np.polyfit(strain, true_data, 3)

                f1 = np.polyval(p1, strain)
                #f2 = np.polyval(p2, strain)
                #f3 = np.polyval(p3, strain)

                Cij[count1, n1] = -0.10*p1[0]
                count1 += 1                

        for n1 in range(0,6):
            for n2 in range(0,6):
                if np.abs(Cij[n1, n2]) < tol:
                    Cij[n1, n2] = 0

                if n2 > 2:
                    Cij[n1, n2] = Cij[n1, n2]*0.50

        return Cij


    def fitCij2(self, tol=0.1):
        # an extension of fitCij, taking into account the possibility of not all stress-calculations having finished yet
        # this method will attempt fitting the elastic constants in those cases of incomplete calculations
        stress_dict = self._sig_eps
        Cij = np.zeros((6, 6))

        inds = [[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]]

        for n1 in range(0, 6):
            strain = []
            stress = []
            count1 = 0

            for c in stress_dict:
                if c.i == inds[n1][0] and c.j== inds[n1][1]:
                    strain.append(c[c.i, c.j])
                    stress.append(stress_dict[c])

            for k in inds:
                true_data = self._chain_stresses(stress, k[0], k[1])
                straink = copy.copy(strain)

                #print true_data, strain


                naninds = []

                for kk in range(0, len(true_data)):
                    if np.isnan(true_data[kk]) == True:
                        naninds.append(kk)


                for kk in naninds[::-1]:
                    #print kk
                    true_data.pop(kk)
                    straink.pop(kk)

                #print straink, true_data
                p1 = np.polyfit(straink, true_data, 1)
                f1 = np.polyval(p1, straink)

                Cij[count1, n1] = -0.10*p1[0]
                count1 += 1


        for n1 in range(0,6):
            for n2 in range(0,6):
                if np.abs(Cij[n1, n2]) < tol:
                    Cij[n1, n2] = 0

                if n2 > 2:
                    Cij[n1, n2] = Cij[n1, n2]*0.50

        return Cij
    
    def fitCij3(self, tol=0.1):
        """
        Quick check to make sure JHM understands the algorithm
        and tests to see if it can be shortened.
        """
        stress_dict = self._sig_eps
        # Initialize Cij array
        inds = [(0,0), (1,1), (2,2), (1,2), (0,2), (0,1)]
        Cij = np.array([[np.polyfit([strain[ind1] for strain in stress_dict.keys() 
                                        if (strain.i,strain.j)==ind1],
                                    [stress_dict[strain][ind2] for strain in stress_dict.keys()
                                        if (strain.i,strain.j)==ind1],1)[0] 
                           for ind1 in inds] for ind2 in inds])
        Cij = -0.1*Cij
        Cij[3:,3:] = 0.5*Cij[3:,3:]
        return Cij



