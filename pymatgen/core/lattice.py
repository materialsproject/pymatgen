#!/usr/bin/env python

"""
This module defines the classes relating to 3D lattices.
"""

from __future__ import division

__author__="Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import numpy as np
import numpy.linalg as npl
from numpy import pi, cos, sin

class Lattice(object):
    '''
    A lattice object.  Essentially a matrix with conversion matrices.
    '''
    
    def __init__(self, matrix):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence is assumed to be read
        one row at a time.  Each row represents one lattice vector.
        
        Args:
            matrix - sequence of numbers in any form. Examples of acceptable input.
                i) An actual numpy array.
                ii) [[1, 0, 0],[0, 1, 0], [0, 0, 1]]
                iii) [1,0,0,0,1,0,0,0,1]
                iv) (1,0,0,0,1,0,0,0,1)
        """
        
        self._matrix = np.array(matrix).reshape((3,3))
        #Store these matrices for faster access
        self._md2c = np.transpose(self._matrix)
        self._mc2d = npl.inv(np.transpose(self._matrix))
    
    @property
    def md2c(self):
        '''Matrix for converting direct to cartesian coordinates'''
        return np.copy(self._md2c)
    
    @property
    def mc2d(self):
        '''Matrix for converting cartesian to direct coordinates'''
        return np.copy(self._mc2d)
    
    @property
    def matrix(self):
        '''Copy of matrix representing the Lattice'''
        return np.copy(self._matrix)
      
    def get_cartesian_coords(self, fractional_coords):
        """
        Returns the cartesian coordinates given fractional coordinates.
        
        Args:
            fractional_coords : Fractional coords.
            
        Returns:
            cartesian coordinates
        """
        return np.transpose(np.dot(self._md2c, np.transpose(fractional_coords)))
    
    def get_fractional_coords(self,cart_coords):
        """
        Returns the fractional coordinates given cartesian coordinates
        
        Args:
            cartesian_coords : Cartesian coords.
            
        Returns:
            fractional coordinates
        """
        return np.transpose(np.dot(self._mc2d, np.transpose(cart_coords)))
    
    @staticmethod
    def cubic(a):
        """
        Returns a new cubic lattice of dimensions a x a x a.
        
        Args:
            a - *a* lattice parameter
            
        Returns:
            Cubic lattice of lattice parameter a.
        """
        return Lattice(np.array([[a, 0.0, 0.0],[0.0, a, 0.0],[0.0, 0.0, a]]))
    
    @staticmethod
    def tetragonal(a,c):
        """
        Returns a new tetragonal lattice of dimensions a x a x c
        """
        return Lattice.from_parameters(a,a,c,90,90,90)

    @staticmethod
    def orthorhombic(a,b,c):
        """
        Returns a new orthorhombic lattice of dimensions a x b x c
        """
        return Lattice.from_parameters(a,b,c,90,90,90)

    @staticmethod
    def monoclinic(a,b,c,alpha):
        """
        Returns a new monoclinic lattice of dimensions a x b x c with angle alpha between lattice vectors b and c
        """
        return Lattice.from_parameters(a,b,c,alpha,90,90)

    @staticmethod
    def hexagonal(a,c):
        """
        Returns a new hexagonal lattice of dimensions a x a x c.
        """
        return Lattice.from_parameters(a,a,c,90.0,90.0,120.0)

    @staticmethod
    def rhombohedral(a):
        """
        Returns a new rhombohedral lattice of dimensions a x a x a.
        """
        return Lattice.from_parameters(a,a,a,60.0,60.0,60.0)    

    @staticmethod
    def from_lengths_and_angles(abc,ang):
        '''
        Create a Lattice using unit cell lengths and angles (in degrees).
        
        Args:
            abc: lattice parameters, e.g. (4,4,5)
            ang: lattice angles in degrees, e.g., (90,90,120)
        '''
        return Lattice.from_parameters(abc[0],abc[1],abc[2],ang[0],ang[1],ang[2])
    
    @staticmethod
    def from_parameters(a,b,c,alpha,beta,gamma):
        '''
        Create a Lattice using unit cell lengths and angles (in degrees).
        a, b, c: lattice parameters
        alpha, beta, gamma : lattice angles
        '''
        to_r = lambda degrees: np.radians(degrees)
        
        alpha_r = to_r(alpha)
        beta_r = to_r(beta)
        gamma_r = to_r(gamma)
        
        gamma_star = np.arccos((np.cos(alpha_r) * np.cos(beta_r)
                               - np.cos(gamma_r)) / (np.sin(alpha_r)
                               * np.sin(beta_r)))
        vector_a = [a * np.sin(to_r(beta)), 0.0, a * np.cos(to_r(beta))]
        vector_b = [-b * np.sin(to_r(alpha)) * np.cos(gamma_star), b
                    * np.sin(to_r(alpha)) * np.sin(gamma_star), b
                    * np.cos(to_r(alpha))]
        vector_c = np.array([0.0, 0.0, float(c)])
        return Lattice([vector_a, vector_b, vector_c])
    
    @staticmethod
    def from_dict(args):
        """
        Create a Lattice from a dictionary containing the a, b, c, alpha, beta, and gamma parameters.
        """
        a = args["a"]
        b = args["b"]
        c = args["c"]
        alpha = args["alpha"]
        beta = args["beta"]
        gamma = args["gamma"]
        return Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    
    @property
    def angles(self):
        """
        returns the angles (alpha, beta, gamma) of the lattice
        """
        return self.lengths_and_angles[1]
        
    @property
    def a(self):
        """
        a, i.e., [1,0,0] lattice parameter
        """
        return self.abc[0]
    
    @property
    def b(self):
        """
        b, i.e., [0,1,0] lattice parameter
        """
        return self.abc[1]
    
    @property
    def c(self):
        """
        c, i.e., [0,0,1] lattice parameter
        """
        return self.abc[2]
    
    @property
    def abc(self):
        """
        Lengths of the lattice vectors
        """
        return self.lengths_and_angles[0]
    
    def _angle_between(self, x, y):
        """ 
        internal method to calculate the angle
        between two vectors
        """
        angle_between = lambda x, y: 180.0 / np.pi * np.arccos(np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y)))
        return angle_between(x, y)
    
    @property
    def alpha(self):
        """
        Angle alpha of lattice
        """
        return self._angle_between(self._matrix[1], self._matrix[2])
    
    @property
    def beta(self):
        """
        Angle beta of lattice
        """
        return self._angle_between(self._matrix[0], self._matrix[2])
    
    @property
    def gamma(self):
        """
        Angle gamma of lattice
        """
        return self._angle_between(self._matrix[0], self._matrix[1])
       
    @property
    def volume(self):
        """
        Volume of the unit cell
        """
        return npl.det(self._matrix)
    
    @property
    def lengths_and_angles(self):
        '''
        Returns (lattice lengths, lattice angles)
        '''
        prim = self._matrix
        lengths= np.sum(prim**2,axis=1)**0.5
        angles = np.zeros((3),float)
        angles[0] = np.arccos(np.dot(prim[1],prim[2])/(lengths[1]*lengths[2]))*180./pi
        angles[1] = np.arccos(np.dot(prim[2],prim[0])/(lengths[2]*lengths[0]))*180./pi
        angles[2] = np.arccos(np.dot(prim[0],prim[1])/(lengths[0]*lengths[1]))*180./pi
        return lengths, np.around(angles,9)
    
    @property
    def reciprocal_lattice(self):
        """
        return the reciprocal lattice
        """
        v = 2 * np.pi / self.volume
        k1 = np.cross(self._matrix[1], self._matrix[2]) * v
        k2 = np.cross(self._matrix[2], self._matrix[0]) * v
        k3 = np.cross(self._matrix[0], self._matrix[1]) * v
        return Lattice([k1, k2, k3])

    def __repr__(self):
        f = lambda x: '%0.6f' % x
        outs = []
        outs.append('Lattice')
        outs.append('    abc : ' + ' '.join(map(f, self.abc)))
        outs.append(' angles : ' + ' '.join(map(f, self.angles)))
        outs.append(' volume : %0.4f' % self.volume)
        outs.append('      A : ' + ' '.join(map(f, self._matrix[0])))
        outs.append('      B : ' + ' '.join(map(f, self._matrix[1])))
        outs.append('      C : ' + ' '.join(map(f, self._matrix[2])))
        return "\n".join(outs)
    
    def __eq__(self,other):
        if other == None:
            return False
        return np.allclose(self._matrix, other._matrix)
    
    def __hash__(self):
        return 7
        
    def __str__(self):
        return '\n'.join([' '.join(["%.6f" % i for i in row]) for row in self._matrix])
    
    @property
    def to_dict(self):
        '''dict representation of the Lattice'''
        d = {
        'matrix': self._matrix.tolist(),
        'a': self.a,
        'b': self.b,
        'c': self.c,
        'alpha': self.alpha,
        'beta': self.beta,
        'gamma': self.gamma,
        'volume': self.volume,
        }
        return d
    
    def get_primitive_lattice(self, lattice_type):
        if lattice_type == 'P':
            return Lattice(self._matrix)
        conv_to_prim = {
            'R': np.array([[2/3, 1/3, 1/3],[-1/3, 1/3, 1/3],[-1/3, -2/3, 1/3]]),
            'A': np.array([[1, 0, 0],[0, 1/2, 1/2],[0, -1/2, 1/2]]),
            'B': np.array([[1/2, 0, 1/2],[0, 1, 0],[-1/2, 0, 1/2]]),
            'C': np.array([[1/2, 1/2, 0],[-1/2, 1/2, 0],[0, 0, 1]]),
            'I': np.array([[-1/2, 1/2, 1/2],[1/2, -1/2, 1/2],[1/2, 1/2, -1/2]]),
            'F': np.array([[1/2, 1/2, 0],[0, 1/2, 1/2],[1/2, 0, 1/2]])
            }
        return Lattice(np.dot(conv_to_prim[lattice_type], self._matrix))
    
    def get_most_compact_basis_on_lattice(self):
        """
        this method get the an alternative basis corresponding to the shortest 3
        linearly independent translational operations permitted.
        This tends to create larger angles for every elongated cells and is 
        beneficial for viewing crystal structure (especially when they are Niggli cells)
        Geoffroy Hautier adapted this code from a java method in jcmc coded by Charles Moore
        """
        matrix=self.matrix
        a=matrix[0]
        b=matrix[1]
        c=matrix[2]
        while True:
            anychange=False
            """
            take care of c
            """
            if(np.dot(a,b)>0):
                diffvector=np.subtract(a,b)
            else:
                diffvector=np.add(a,b)
            if(np.linalg.norm(diffvector)<np.linalg.norm(a)
               or np.linalg.norm(diffvector)<np.linalg.norm(b)):
                if(np.linalg.norm(a)<np.linalg.norm(b)):
                    b=diffvector
                else:
                    a=diffvector
            anychange=True
            """
            take care of b
            """
            if(np.dot(a,c)>0):
                diffvector=np.subtract(a,c)
            else:
                diffvector=np.add(a,c)
            if(np.linalg.norm(diffvector)<np.linalg.norm(a)
               or np.linalg.norm(diffvector)<np.linalg.norm(c)):
                if(np.linalg.norm(a)<np.linalg.norm(c)):
                    c=diffvector
                else:
                    a=diffvector
            anychange=True       
            """
            take care of a
            """
            if(np.dot(c,b)>0):
                diffvector=np.subtract(c,b)
            else:
                diffvector=np.add(c,b)
            if(np.linalg.norm(diffvector)<np.linalg.norm(c)
               or np.linalg.norm(diffvector)<np.linalg.norm(b)):
                if(np.linalg.norm(c)<np.linalg.norm(b)):
                    b=diffvector
                else:
                    c=diffvector
            anychange=True
            if anychange==True:
                break
        return Lattice([a,b,c])
