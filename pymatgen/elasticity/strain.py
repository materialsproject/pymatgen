from __future__ import absolute_import
from pymatgen.core.lattice import Lattice
from pymatgen.elasticity.tensors import SQTensor
import warnings
import numpy as np
from six.moves import zip

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Mark Asta, Anubhav Jain"
__version__ = "1.0"
__maintainer__ = "Maarten de Jong"
__email__ = "maartendft@gmail.com"
__status__ = "Development"
__date__ = "March 13, 2012"


class Deformation(SQTensor):
    """
    Subclass of SQTensor that describes the deformation gradient tensor
    """

    def __new__(cls, deformation_gradient):
        """
        Create a Deformation object.  Note that the constructor uses __new__ 
        rather than __init__ according to the standard method of subclassing 
        numpy ndarrays.

        Args:
            deformation_gradient (3x3 array-like): the 3x3 array-like
                representing the deformation gradient
        """

        obj = SQTensor(deformation_gradient).view(cls)
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

    def __repr__(self):
        return "Deformation({})".format(self.__str__())

    def check_independent(self):
        """
        checks to determine whether the deformation matrix represents an
        independent deformation, raises a ValueError if not.  If so, returns 
        the indices of the deformation gradient entry representing the 
        independent component
        """
        indices = list(zip(*np.asarray(self - np.eye(3)).nonzero()))
        if len(indices) != 1:
            raise ValueError("One and only one independent deformation"
                             "must be applied.")
        return indices[0]

    @property
    def green_lagrange_strain(self):
        """
        calculates the euler-lagrange strain from
        the deformation gradient
        """
        return Strain.from_deformation(self)

    def apply_to_structure(self, structure):
        """
        Apply the deformation gradient to a structure.
        
        Args:
            structure (Structure object): the structure object to
                be modified by the deformation
        """
        def_struct = structure.copy()
        def_struct.modify_lattice(Lattice(np.dot(def_struct.lattice.matrix,
                                                 self)))
        return def_struct

    @classmethod
    def from_index_amount(cls, matrixpos, amt):
        """
        Factory method for constructing a Deformation object
        from a matrix position and amount

        Args:
            matrixpos (tuple): tuple corresponding the matrix position to
                have a perturbation added
            amt (float): amount to add to the identity matrix at position 
                matrixpos
        """
        f = np.identity(3)
        f[matrixpos] += amt
        return cls(f)


class DeformedStructureSet(object):
    """
    class that generates a set of independently deformed structures that
    can be used to calculate linear stress-strain response
    """

    def __init__(self, rlxd_str, nd=0.01, ns=0.08,
                 num_norm=4, num_shear=4, symmetry=False):
        """
        constructs the deformed geometries of a structure.  Generates
        m + n deformed structures according to the supplied parameters.

        Args:
            rlxd_str (structure): structure to undergo deformation, if 
                fitting elastic tensor is desired, should be a geometry 
                optimized structure
            nd (float): maximum perturbation applied to normal deformation
            ns (float): maximum perturbation applied to shear deformation
            m (int): number of deformation structures to generate for 
                normal deformation, must be even
            n (int): number of deformation structures to generate for 
                shear deformation, must be even
        """

        if num_norm % 2 != 0:
            raise ValueError("Number of normal deformations (num_norm)"
                             " must be even.")
        if num_shear % 2 != 0:
            raise ValueError("Number of shear deformations (num_shear)"
                             " must be even.")

        norm_deformations = np.linspace(-nd, nd, num=num_norm + 1)
        norm_deformations = norm_deformations[norm_deformations.nonzero()]
        shear_deformations = np.linspace(-ns, ns, num=num_shear + 1)
        shear_deformations = shear_deformations[shear_deformations.nonzero()]

        self.undeformed_structure = rlxd_str
        self.deformations = []
        self.def_structs = []
        if symmetry:
            raise NotImplementedError("Symmetry reduction of deformed "
                                      "structure set not yet implemented")
        else:
            self.symmetry = None
            # Determine normal deformation gradients
            # Apply normal deformations
            for ind in [(0, 0), (1, 1), (2, 2)]:
                for amount in norm_deformations:
                    defo = Deformation.from_index_amount(ind, amount)
                    self.deformations.append(defo)
                    self.def_structs.append(defo.apply_to_structure(rlxd_str))

            # Apply shear deformations 
            for ind in [(0, 1), (0, 2), (1, 2)]:
                for amount in shear_deformations:
                    defo = Deformation.from_index_amount(ind, amount)
                    self.deformations.append(defo)
                    self.def_structs.append(defo.apply_to_structure(rlxd_str))

    def __iter__(self):
        return iter(self.def_structs)

    def as_strain_dict(self):
        """
        Returns dictionary of deformed structures indexed by independent
        strain objects in accordance with legacy behavior of elasticity
        package
        """
        strains = [IndependentStrain(defo) for defo in self.deformations]
        return dict(zip(strains, self.def_structs))


class Strain(SQTensor):
    """
    Subclass of SQTensor that describes the Green-Lagrange strain tensor.
    """

    def __new__(cls, strain_matrix, dfm=None):
        """
        Create a Strain object.  Note that the constructor uses __new__ 
        rather than __init__ according to the standard method of 
        subclassing numpy ndarrays.  Note also that the default constructor
        does not include the deformation gradient

        Args:
            strain_matrix (3x3 array-like): the 3x3 array-like
                representing the Green-Lagrange strain
        """

        obj = SQTensor(strain_matrix).view(cls)
        obj._dfm = dfm
        if not obj.is_symmetric():
            raise ValueError("Strain objects must be initialized "
                             "with a symmetric array-like.")

        if dfm is None:
            warnings.warn("Constructing a strain object without a deformation "
                          "matrix makes many methods unusable.  Use "
                          "Strain.from_deformation to construct a Strain object"
                          " from a deformation gradient.")
        elif (np.array(dfm) - obj < 1e-5).all():
            warnings.warn("Warning: deformation matrix does not correspond "
                          "to input strain_matrix value")
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._dfm = getattr(obj, "_dfm", None)

    def __repr__(self):
        return "Strain({})".format(self.__str__())

    @classmethod
    def from_deformation(cls, deformation):
        """
        Factory method that returns a Strain object from a deformation
        gradient

        Args:
            deformation (3x3 array-like):
        """
        dfm = Deformation(deformation)
        return cls(0.5 * (np.dot(dfm.T,dfm) - np.eye(3)), dfm)

    @property
    def deformation_matrix(self):
        """
        returns the deformation matrix
        """
        return self._dfm

    @property
    def independent_deformation(self):
        """
        determines whether the deformation matrix represents an
        independent deformation, raises a value error if not.
        Returns the index of the deformation gradient corresponding
        to the independent deformation
        """
        if self._dfm is None:
            raise ValueError("No deformation matrix supplied "
                             "for this strain tensor.")
        return self._dfm.check_independent()

    @property
    def voigt(self):
        """
        translates a strain tensor into a voigt notation vector
        """
        return [self[0, 0], self[1, 1], self[2, 2],
                2. * self[1, 2], 2. * self[0, 2], 2. * self[0, 1]]


class IndependentStrain(Strain):
    """
    Class for independent strains intended for use with old Materials Project
    elasticity workflow.  Note that the default constructor constructs from 
    a deformation matrix, rather than an array representing the strain, to 
    emulate the legacy behavior.
    """

    def __new__(cls, deformation_gradient):
        """
        Create an Independent Strain object.  Note that the constructor uses 
        __new__ rather than __init__ according to the standard method of 
        subclassing numpy ndarrays.  Note also that, unlike the Strain class,
        the default constructor of IndependentStrain takes the deformation 
        gradient as input, rather than an array representing the Green-Lagrange 
        strain.

        Args:
            deformation_gradient (3x3 array-like): the 3x3 array-like
                representing the deformation gradient
        """

        obj = Strain.from_deformation(deformation_gradient).view(cls)
        (obj._i, obj._j) = obj.independent_deformation
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._dfm = getattr(obj, "_dfm", None)
        self._i = getattr(obj, "_i", None)
        self._j = getattr(obj, "_j", None)

    def __repr__(self):
        return "IndependentStrain({})".format(self.__str__())
 
    @property
    def i(self):
        return self._i

    @property
    def j(self):
        return self._j
