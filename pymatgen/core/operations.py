#!/usr/bin/env python

"""
This module provides classes that operate on points or vectors in 3D space.
"""

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import numpy as np
from math import sin, cos, pi

class SymmOp (object):
    """
    A symmetry operation in cartesian space.
    Consists of a rotation plus a translation.
    Implementation is as an affine transformation matrix of rank 4 for efficiency.
    Read: http://en.wikipedia.org/wiki/Affine_transformation
    """
    
    def __init__(self, affine_transformation_matrix, tol = 0.01):
        """
        Initializes the SymmOp from a 4x4 affine transformation matrix.
        In general, this constructor should not be used unless you are 
        transferring rotations.  Use the static constructors instead to
        generate a SymmOp from proper rotations and translation.
        
        Args:
            affine_transformation_matrix: A 4x4 numpy.array representing an affine transformation.
            tol: Tolerance for determining if matrices are equal.
        """
        if affine_transformation_matrix.shape != (4,4):
            raise ValueError("Affine Matrix must be a 4x4 numpy array!")
        self._matrix = affine_transformation_matrix
        self._tol = tol
    
    @staticmethod
    def from_rotation_matrix_and_translation_vector(rotation_matrix, translation_vec, tol = 0.1):
        """
        Creates a symmetry operation from a rotatino matrix and a translation vector.
        
        Args:
            rotation_matrix: A 3x3 numpy.array specifying a rotation matrix
            translation_vec: A rank 1 numpy.array specifying a translation vector
            tol: tolerance to determine if rotation matrix is valid
        """
        if rotation_matrix.shape != (3,3):
            raise ValueError("Rotation Matrix must be a 3x3 numpy array.")
        if translation_vec.shape != (3,):
            raise ValueError("Translation vector must be a rank 1 numpy array with 3 elements.")
        affine_matrix = np.eye(4)
        affine_matrix[0:3][:,0:3] = rotation_matrix
        affine_matrix[0:3][:,3] = translation_vec
        return SymmOp(affine_matrix, tol)
               
    def __eq__(self, other):
        return (abs(self._matrix - other._matrix) < self._tol).all()

    def __hash__(self):
        return 7

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        output = ["Rot:"]
        output.append(str(self._matrix[0:3][:,0:3]))
        output.append("tau")
        output.append(str(self._matrix[0:3][:,3]))
        return "\n".join(output)
    
    def operate(self, point):
        """
        Apply the operation on a point.
        
        Args:
            point - a cartesian coordinate represented as a rank 1 numpy array of 3 elements.
        """
        affine_point = np.array([point[0],point[1],point[2],1])
        affine_point[0:3] = point
        return np.dot(self._matrix, affine_point)[0:3]

    def apply_rotation_only(self, vector):
        """
        Vectors should only be operated by the rotation matrix and not the translation vector
        
        Args:
            vector - a rank 1 numpy array of 3 elements representing a vector.
        """
        return np.dot(self.rotation_matrix, vector)
    
    def are_symmetrically_related(self, point_a, point_b, tol = 0.001):
        """
        Checks if two points are symmetrically related.
        
        Args:
            point_a - Point a
            point_b - Point b
            tol - tolerance for checking.
        Returns:
            True if self.operate(point_a) == point_b or vice versa.
        """
        if (abs(self.operate(point_a) - point_b) < tol).all():
            return True
        if (abs(self.operate(point_b) - point_a) < tol).all():
            return True
        return False
    
    @property
    def affine_matrix(self):
        """
        A 4x4 numpy.array representing the symmetry operation.
        """
        return self._matrix
    
    @property
    def rotation_matrix(self):
        """
        A 3x3 numpy.array representing the rotation matrix
        """
        return self._matrix[0:3][:,0:3]
    
    @property
    def translation_vector(self):
        """
        A rank 1 numpy.array of dim 3 representing the translation vector.
        """
        return self._matrix[0:3][:,3]
    
    def __mul__(self, other):
        """
        Returns a new SymmOp which is equivalent to apply the "other" SymmOp followed by this one.
        """
        new_matrix = np.dot(self._matrix, other._matrix)
        return SymmOp(new_matrix)

    @property
    def inverse(self):
        """
        Returns inverse of transformation.
        """
        invr = np.linalg.inv(self._matrix)
        return SymmOp(invr)

    @staticmethod
    def from_axis_angle_and_translation(axis, angle, angle_in_radians = False, translation_vec = np.zeros(3)):
        """
        Generates a SymmOp for a rotation about a given axis plus a translation.
        
        Args:
            axis - The axis of rotation in cartesian space. For example, [1,0,0] indicates rotation about x-axis.
            angle - The angle of rotation.
            angle_in_radians - Set to True if angles are given in radians.  Else, units of degrees is assumed.
            translation_vec - A translation vector.  Defaults to zero.
            
        """
        if isinstance(axis, (tuple, list)):
            axis = np.array(axis)
        a = angle if angle_in_radians else angle * pi / 180
        cosa = cos(a)
        sina = sin(a) 
        u = axis / np.linalg.norm(axis)
        r = np.zeros((3,3)) 
        r[0,0] = cosa + u[0] ** 2 * (1-cosa)
        r[0,1] = u[0] * u[1] * (1-cosa) - u[2] * sina
        r[0,2] = u[0] * u[2] * (1-cosa) + u[1] * sina
        r[1,0] = u[0] * u[1] * (1-cosa) + u[2] * sina     
        r[1,1] = cosa + u[1]**2 * (1-cosa)     
        r[1,2] = u[1] * u[2] * (1-cosa) - u[0] * sina     
        r[2,0] = u[0] * u[2] * (1-cosa) - u[1] * sina     
        r[2,1] = u[1] * u[2] * (1-cosa) + u[0] * sina     
        r[2,2] = cosa + u[2]**2 * (1-cosa)    
        
        return SymmOp.from_rotation_matrix_and_translation_vector(r, translation_vec)
