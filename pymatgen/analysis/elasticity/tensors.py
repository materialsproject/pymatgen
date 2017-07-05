# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from scipy.linalg import polar
import numpy as np
import sympy as sp
import itertools
import warnings
import collections
import string
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.core.lattice import Lattice

"""
This module provides a base class for tensor-like objects and methods for
basic tensor manipulation.  It also provides a class, SquareTensor,
that provides basic methods for creating and manipulating rank 2 tensors
"""

__author__ = "Maarten de Jong"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = ("Joseph Montoya, Shyam Dwaraknath, Wei Chen, "
               "Mark Asta, Anubhav Jain, Terence Lew")
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "March 22, 2012"

voigt_map = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
reverse_voigt_map = np.array([[0, 5, 4],
                              [5, 1, 3],
                              [4, 3, 2]])


class Tensor(np.ndarray):
    """
    Base class for doing useful general operations on Nth order tensors,
    without restrictions on the type (stress, elastic, strain, piezo, etc.)
    """

    def __new__(cls, input_array, vscale=None, check_rank=None):
        """
        Create a Tensor object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            input_array: (array-like with shape 3^N): array-like representing
                a tensor quantity in standard (i. e. non-voigt) notation
            vscale: (N x M array-like): a matrix corresponding
                to the coefficients of the voigt-notation tensor
        """
        obj = np.asarray(input_array).view(cls)
        obj.rank = len(obj.shape)

        if check_rank and check_rank != obj.rank:
            raise ValueError("{} input must be rank {}".format(
                obj.__class__.__name__, check_rank))

        vshape = tuple([3] * (obj.rank % 2) + [6] * (obj.rank // 2))
        obj._vscale = np.ones(vshape)
        if vscale is not None:
            obj._vscale = vscale
        if obj._vscale.shape != vshape:
            raise ValueError("Voigt scaling matrix must be the shape of the "
                             "voigt notation matrix or vector.")
        if not all([i == 3 for i in obj.shape]):
            raise ValueError("Pymatgen only supports 3-dimensional tensors, "
                             "and default tensor constructor uses standard "
                             "notation.  To construct from voigt notation, use"
                             " {}.from_voigt".format(obj.__class__.__name__))
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.rank = getattr(obj, 'rank', None)
        self._vscale = getattr(obj, '_vscale', None)
        self._vdict = getattr(obj, '_vdict', None)

    def __array_wrap__(self, obj):
        """
        Overrides __array_wrap__ methods in ndarray superclass to avoid errors
        associated with functions that return scalar values
        """

        if len(obj.shape) == 0:
            return obj[()]
        else:
            return np.ndarray.__array_wrap__(self, obj)

    def __hash__(self):
        """
        define a hash function, since numpy arrays
        have their own __eq__ method
        """
        return hash(self.tostring())

    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,
                               self.__str__())

    def zeroed(self, tol=1e-3):
        """
        returns the matrix with all entries below a certain threshold
        (i.e. tol) set to zero
        """
        new_tensor = self.copy()
        new_tensor[abs(new_tensor) < tol] = 0
        return new_tensor

    def transform(self, symm_op):
        """
        Applies a transformation (via a symmetry operation) to a tensor. 

        Args:
            symm_op (SymmOp): a symmetry operation to apply to the tensor
        """
        return self.__class__(symm_op.transform_tensor(self))

    def rotate(self, matrix, tol=1e-3):
        """
        Applies a rotation directly, and tests input matrix to ensure a valid
        rotation.

        Args:
            matrix (3x3 array-like): rotation matrix to be applied to tensor
            tol (float): tolerance for testing rotation matrix validity
        """
        matrix = SquareTensor(matrix)
        if not matrix.is_rotation(tol):
            raise ValueError("Rotation matrix is not valid.")
        sop = SymmOp.from_rotation_and_translation(matrix,
                                                   [0., 0., 0.])
        return self.transform(sop)

    def einsum_sequence(self, other_arrays, einsum_string=None):
        """
        Calculates the result of an einstein summation expression
        """
        if not isinstance(other_arrays, list):
            raise ValueError("other tensors must be list of "
                             "tensors or tensor input")

        other_arrays = [np.array(a) for a in other_arrays]
        if not einsum_string:
            lc = string.ascii_lowercase
            einsum_string = lc[:self.rank]
            other_indices = ''
            other_ranks = [len(a.shape) for a in other_arrays]
            idx = self.rank - sum(other_ranks)
            for length in other_ranks:
                einsum_string += ',' + lc[idx:idx + length]
                idx += length

        einsum_args = [self] + list(other_arrays)
        return np.einsum(einsum_string, *einsum_args)

    @property
    def symmetrized(self):
        """
        Returns a generally symmetrized tensor, calculated by taking 
        the sum of the tensor and its transpose with respect to all 
        possible permutations of indices
        """
        perms = list(itertools.permutations(range(self.rank)))
        return sum([np.transpose(self, ind) for ind in perms]) / len(perms)

    @property
    def voigt_symmetrized(self):
        """
        Returns a "voigt"-symmetrized tensor, i. e. a voigt-notation
        tensor such that it is invariant wrt permutation of indices
        """
        if not (self.rank % 2 == 0 and self.rank > 2):
            raise ValueError("V-symmetrization requires rank even and > 2")

        v = self.voigt
        perms = list(itertools.permutations(range(len(v.shape))))
        new_v = sum([np.transpose(v, ind) for ind in perms]) / len(perms)
        return self.__class__.from_voigt(new_v)

    def is_symmetric(self, tol=1e-5):
        """
        Tests whether a tensor is symmetric or not based on the residual
        with its symmetric part, from self.symmetrized

        Args:
            tol (float): tolerance to test for symmetry
        """
        return (self - self.symmetrized < tol).all()

    def fit_to_structure(self, structure, symprec=0.1):
        """
        Returns a tensor that is invariant with respect to symmetry
        operations corresponding to a structure

        Args: 
            structure (Structure): structure from which to generate 
                symmetry operations
            symprec (float): symmetry tolerance for the Spacegroup Analyzer
                used to generate the symmetry operations
        """
        sga = SpacegroupAnalyzer(structure, symprec)
        symm_ops = sga.get_symmetry_operations(cartesian=True)
        return sum([self.transform(symm_op)
                    for symm_op in symm_ops]) / len(symm_ops)

    def is_fit_to_structure(self, structure, tol=1e-2):
        """
        Tests whether a tensor is invariant with respect to the
        symmetry operations of a particular structure by testing
        whether the residual of the symmetric portion is below a 
        tolerance
        
        Args:
            structure (Structure): structure to be fit to
            tol (float): tolerance for symmetry testing
        """
        return (self - self.fit_to_structure(structure) < tol).all()

    @property
    def voigt(self):
        """
        Returns the tensor in Voigt notation
        """
        v_matrix = np.zeros(self._vscale.shape, dtype=self.dtype)
        this_voigt_map = self.get_voigt_dict(self.rank)
        for ind in this_voigt_map:
            v_matrix[this_voigt_map[ind]] = self[ind]
        if not self.is_voigt_symmetric():
            warnings.warn("Tensor is not symmetric, information may "
                          "be lost in voigt conversion.")
        return v_matrix * self._vscale

    def is_voigt_symmetric(self, tol=1e-6):
        """
        Tests symmetry of tensor to that necessary for voigt-conversion
        by grouping indices into pairs and constructing a sequence of
        possible permutations to be used in a tensor transpose
        """
        transpose_pieces = [[[0 for i in range(self.rank % 2)]]]
        transpose_pieces += [[range(j, j + 2)] for j in
                             range(self.rank % 2, self.rank, 2)]
        for n in range(self.rank % 2, len(transpose_pieces)):
            if len(transpose_pieces[n][0]) == 2:
                transpose_pieces[n] += [transpose_pieces[n][0][::-1]]
        for trans_seq in itertools.product(*transpose_pieces):
            trans_seq = list(itertools.chain(*trans_seq))
            if (self - self.transpose(trans_seq) > tol).any():
                return False
        return True

    @staticmethod
    def get_voigt_dict(rank):
        """
        Returns a dictionary that maps indices in the tensor to those
        in a voigt representation based on input rank

        Args:
            rank (int): Tensor rank to generate the voigt map
        """
        vdict = {}
        for ind in itertools.product(*[range(3)] * rank):
            v_ind = ind[:rank % 2]
            for j in range(rank // 2):
                pos = rank % 2 + 2 * j
                v_ind += (reverse_voigt_map[ind[pos:pos + 2]],)
            vdict[ind] = v_ind
        return vdict

    @classmethod
    def from_voigt(cls, voigt_input):
        """
        Constructor based on the voigt notation vector or matrix.

        Args: 
            voigt_input (array-like): voigt input for a given tensor
        """
        voigt_input = np.array(voigt_input)
        rank = sum(voigt_input.shape) // 3
        t = cls(np.zeros([3] * rank))
        if voigt_input.shape != t._vscale.shape:
            raise ValueError("Invalid shape for voigt matrix")
        voigt_input = voigt_input / t._vscale
        this_voigt_map = t.get_voigt_dict(rank)
        for ind in this_voigt_map:
            t[ind] = voigt_input[this_voigt_map[ind]]
        return cls(t)

    def convert_to_ieee(self, structure, initial_fit=True):
        """
        Given a structure associated with a tensor, attempts a
        calculation of the tensor in IEEE format according to
        the 1987 IEEE standards.

        Args:
            structure (Structure): a structure associated with the
                tensor to be converted to the IEEE standard
            initial_fit (bool): flag to indicate whether initial
                tensor is fit to the symmetry of the structure.
                Defaults to true. Note that if false, inconsistent
                results may be obtained due to symmetrically
                equivalent, but distinct transformations
                being used in different versions of spglib.
        """

        # Check conventional setting:
        sga = SpacegroupAnalyzer(structure)
        dataset = sga.get_symmetry_dataset()
        trans_mat = dataset['transformation_matrix']
        conv_latt = Lattice(np.transpose(np.dot(np.transpose(
            structure.lattice.matrix), np.linalg.inv(trans_mat))))
        xtal_sys = sga.get_crystal_system()

        vecs = conv_latt.matrix
        lengths = np.array(conv_latt.abc)
        angles = np.array(conv_latt.angles)
        rotation = np.zeros((3, 3))

        # IEEE rules: a,b,c || x1,x2,x3
        if xtal_sys == "cubic":
            rotation = [vecs[i] / lengths[i] for i in range(3)]

        # IEEE rules: a=b in length; c,a || x3, x1
        elif xtal_sys == "tetragonal":
            rotation = np.array([vec / mag for (mag, vec) in
                                 sorted(zip(lengths, vecs),
                                        key=lambda x: x[0])])
            if abs(lengths[2] - lengths[1]) < abs(lengths[1] - lengths[0]):
                rotation[0], rotation[2] = rotation[2], rotation[0].copy()
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[0]))

        # IEEE rules: c<a<b; c,a || x3,x1
        elif xtal_sys == "orthorhombic":
            rotation = [vec / mag for (mag, vec) in sorted(zip(lengths, vecs))]
            rotation = np.roll(rotation, 2, axis=0)

        # IEEE rules: c,a || x3,x1, c is threefold axis
        # Note this also includes rhombohedral crystal systems
        elif xtal_sys in ("trigonal", "hexagonal"):
            # find threefold axis:
            tf_index = np.argmin(abs(angles - 120.))
            non_tf_mask = np.logical_not(angles == angles[tf_index])
            rotation[2] = get_uvec(vecs[tf_index])
            rotation[0] = get_uvec(vecs[non_tf_mask][0])
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[0]))

        # IEEE rules: b,c || x2,x3; alpha=beta=90, c<a
        elif xtal_sys == "monoclinic":
            # Find unique axis
            u_index = np.argmax(abs(angles - 90.))
            n_umask = np.logical_not(angles == angles[u_index])
            rotation[1] = get_uvec(vecs[u_index])
            # Shorter of remaining lattice vectors for c axis
            c = [vec / mag for (mag, vec) in
                 sorted(zip(lengths[n_umask], vecs[n_umask]))][0]
            rotation[2] = np.array(c)
            rotation[0] = np.cross(rotation[1], rotation[2])

        # IEEE rules: c || x3
        elif xtal_sys == "triclinic":
            rotation = [vec / mag for (mag, vec) in sorted(zip(lengths, vecs))]
            rotation = np.roll(rotation, 2, axis=0)
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[1]))
            rotation[0] = np.cross(rotation[1], rotation[2])

        result = self.copy()
        if initial_fit:
            result = result.fit_to_structure(structure)
        return result.rotate(rotation, tol=1e-2)

    @classmethod
    def from_values_indices(cls, values, indices, populate=False,
                            structure=None, voigt_rank=None, 
                            vsym=True, verbose=False):
        """
        Creates a tensor from values and indices, with options
        for populating the remainder of the tensor.
        
        Args:
            values (floats): numbers to place at indices
            indices (array-likes): indices to place values at
            populate (bool): whether to populate the tensor
            structure (Structure): structure to base population
                or fit_to_structure on
            voigt_rank (int): full tensor rank to indicate the
                shape of the resulting tensor.  This is necessary
                if one provides a set of indices more minimal than
                the shape of the tensor they want, e.g. 
                Tensor.from_values_indices((0, 0), 100)
            vsym (bool): whether to voigt symmetrize during the
                optimization procedure
        """
        # auto-detect voigt notation
        # TODO: refactor rank inheritance to make this easier
        indices = np.array(indices)
        if voigt_rank:
            shape = ([3]*(voigt_rank % 2) + [6]*(voigt_rank // 2))
        else:
            shape = np.ceil(np.max(indices+1, axis=0) / 3.) * 3
        base = np.zeros(shape.astype(int))
        for v, idx in zip(values, indices):
            base[tuple(idx)] = v
        if 6 in shape:
            obj = cls.from_voigt(base)
        else:
            obj = cls(base)
        if populate:
            assert structure, "Populate option must include structure input"
            obj = obj.populate(structure, vsym=vsym, verbose=verbose)
        elif structure:
            obj = obj.fit_to_structure(structure)
        return obj

    def populate(self, structure, prec=1e-5, maxiter=200, verbose=False,
                 precond=True, vsym=True):
        """
        Takes a partially populated tensor, and populates the non-zero
        entries according to the following procedure, iterated until
        the desired convergence (specified via prec) is achieved.

        1. Find non-zero entries
        2. Symmetrize the tensor with respect to crystal symmetry and
           (optionally) voigt symmetry
        3. Reset the non-zero entries of the original tensor
        
        Args:
            structure (structure object)
        """
        if precond:
            # Generate the guess from populated
            sops = SpacegroupAnalyzer(structure).get_symmetry_operations()
            guess = Tensor(np.zeros(self.shape))
            mask = abs(self) > prec
            guess[mask] = self[mask]

            def merge(old, new):
                gmask = np.abs(old) > prec
                nmask = np.abs(new) > prec
                new_mask = np.logical_not(gmask)*nmask
                avg_mask = gmask*nmask
                old[avg_mask] = (old[avg_mask] + new[avg_mask]) / 2.
                old[new_mask] = new[new_mask]
            if verbose:
                print("Preconditioning for {} symmops".format(len(sops)))
            for sop in sops:
                rot = guess.transform(sop)
                # Store non-zero entries of new that weren't previously
                # in the guess in the guess
                merge(guess, rot)
            if verbose:
                print("Preconditioning for voigt symmetry".format(len(sops)))
            if vsym:
                v = guess.voigt
                perms = list(itertools.permutations(range(len(v.shape))))
                for perm in perms:
                    vtrans = np.transpose(v, perm)
                    merge(v, vtrans)
                guess = Tensor.from_voigt(v)
        else:
            guess = np.zeros(self.shape)

        assert guess.shape == self.shape, "Guess must have same shape"
        converged = False
        test_new, test_old = [guess.copy()]*2
        for i in range(maxiter):
            test_new = test_old.fit_to_structure(structure)
            if vsym:
                test_new = test_new.voigt_symmetrized
            diff = np.abs(test_old - test_new) 
            converged = (diff < prec).all()
            if converged:
                break
            test_new[mask] = self[mask]
            test_old = test_new
            if verbose:
                print("Iteration {}: {}".format(i, np.max(diff)))
        if not converged:
            max_diff = np.max(np.abs(self - test_new))
            warnings.warn("Warning, populated tensor is not converged "
                          "with max diff of {}".format(max_diff))
        return self.__class__(test_new)



class TensorCollection(collections.Sequence):
    """
    A sequence of tensors that can be used for fitting data
    or for having a tensor expansion
    """
    def __init__(self, tensor_list, base_class=Tensor):
        self.tensors = [base_class(t) if not isinstance(t, base_class)
                        else t for t in tensor_list]

    def __len__(self):
        return len(self.tensors)

    def __getitem__(self, ind):
        return self.tensors[ind]

    def __iter__(self):
        return self.tensors.__iter__()

    def zeroed(self, tol=1e-3):
        return self.__class__([t.zeroed(tol) for t in self])

    def transform(self, symm_op):
        return self.__class__([t.transform(symm_op) for t in self])

    def rotate(self, matrix, tol=1e-3):
        return self.__class__([t.rotate(matrix, tol) for t in self])

    @property
    def symmetrized(self):
        return self.__class__([t.symmetrized for t in self])

    def is_symmetric(self, tol=1e-5):
        return all([t.is_symmetric(tol) for t in self])

    def fit_to_structure(self, structure, symprec=0.1):
        return self.__class__([t.fit_to_structure(structure, symprec) 
                               for t in self])
    
    def is_fit_to_structure(self, structure, tol=1e-2):
        return all([t.is_fit_to_structure(structure, tol) for t in self])

    @property
    def voigt(self):
        return [t.voigt for t in self]
    
    @property
    def ranks(self):
        return [t.rank for t in self]

    def is_voigt_symmetric(self, tol=1e-6):
        return all([t.is_voigt_symmetric(tol) for t in self])

    @classmethod
    def from_voigt(cls, voigt_input_list, base_class=Tensor):
        return cls([base_class.from_voigt(v) for v in voigt_input_list])

    def convert_to_ieee(self, structure):
        return self.__class__([t.convert_to_ieee(structure) for t in self])


class SquareTensor(Tensor):
    """
    Base class for doing useful general operations on second rank tensors
    (stress, strain etc.).
    """

    def __new__(cls, input_array, vscale=None):
        """
        Create a SquareTensor object.  Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.  Error is thrown when the class is
        initialized with non-square matrix.

        Args:
            input_array (3x3 array-like): the 3x3 array-like
                representing the content of the tensor
            vscale (6x1 array-like): 6x1 array-like scaling the
                voigt-notation vector with the tensor entries
        """

        obj = super(SquareTensor, cls).__new__(cls, input_array, vscale,
                                               check_rank=2)
        return obj.view(cls)

    @property
    def trans(self):
        """
        shorthand for transpose on SquareTensor
        """
        return SquareTensor(np.transpose(self))

    @property
    def inv(self):
        """
        shorthand for matrix inverse on SquareTensor
        """
        if self.det == 0:
            raise ValueError("SquareTensor is non-invertible")
        return SquareTensor(np.linalg.inv(self))

    @property
    def det(self):
        """
        shorthand for the determinant of the SquareTensor
        """
        return np.linalg.det(self)

    def is_rotation(self, tol=1e-3, include_improper=True):
        """
        Test to see if tensor is a valid rotation matrix, performs a
        test to check whether the inverse is equal to the transpose
        and if the determinant is equal to one within the specified
        tolerance

        Args:
            tol (float): tolerance to both tests of whether the
                the determinant is one and the inverse is equal
                to the transpose
            include_improper (bool): whether to include improper
                rotations in the determination of validity
        """
        det = np.abs(np.linalg.det(self))
        if include_improper:
            det = np.abs(det)
        return (np.abs(self.inv - self.trans) < tol).all() \
            and (np.abs(det - 1.) < tol)

    def get_scaled(self, scale_factor):
        """
        Scales the tensor by a certain multiplicative scale factor

        Args:
            scale_factor (float): scalar multiplier to be applied to the
                SquareTensor object
        """
        return SquareTensor(self * scale_factor)

    @property
    def principal_invariants(self):
        """
        Returns a list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic
        polynomial for the matrix
        """
        return np.poly(self)[1:] * np.array([-1, 1, -1])

    def polar_decomposition(self, side='right'):
        """
        calculates matrices for polar decomposition
        """
        return polar(self, side=side)


def get_uvec(vec):
    """ Gets a unit vector parallel to input vector"""
    l = np.linalg.norm(vec)
    if l < 1e-8:
        return vec
    return vec / l

def symmetry_reduce(tensors, structure, tol=1e-8, **kwargs):
    """
    Function that converts a list of tensors corresponding to a structure
    and returns a dictionary consisting of unique tensor keys with symmop
    values corresponding to transformations that will result in derivative
    tensors from the original list

    Args:
        tensors (list of tensors): list of Tensor objects to test for
            symmetrically-equivalent duplicates
        structure (Structure): structure from which to get symmetry
        tol (float): tolerance for tensor equivalence
        kwargs: keyword arguments for the SpacegroupAnalyzer

    returns:
        dictionary consisting of unique tensors with symmetry operations
        corresponding to those which will reconstruct the remaining
        tensors as values
    """
    sga = SpacegroupAnalyzer(structure, **kwargs)
    symmops = sga.get_symmetry_operations(cartesian=True)
    unique_tdict = {}
    for tensor in tensors:
        is_unique = True
        for unique_tensor, symmop in itertools.product(unique_tdict, symmops):
            if (np.abs(unique_tensor.transform(symmop) - tensor) < tol).all():
                unique_tdict[unique_tensor].append(symmop)
                is_unique = False
                break
        if is_unique:
            unique_tdict[tensor] = []
    return unique_tdict

#TODO: abstract this
def diff_fit(independent, dependent, eq_state=None, order=2, tol=1e-10):
    """
    nth order elastic constant fitting function based on 
    central-difference derivatives with respect to distinct
    strain states.  The algorithm is summarized as follows:

    1. Identify distinct strain states as sets of indices 
       for which nonzero strain values exist, typically
       [(0), (1), (2), (3), (4), (5), (0, 1) etc.]
    2. For each strain state, find and sort strains and
       stresses by strain value.
    3. Find first, second .. nth derivatives of each stress
       with respect to scalar variable corresponding to
       the smallest perturbation in the strain.
    4. Use the pseudoinverse of a matrix-vector expression 
       corresponding to the parameterized stress-strain
       relationship and multiply that matrix by the respective 
       calculated first or second derivatives from the
       previous step.
    5. Place the calculated nth-order elastic 
       constants appropriately.

    Args:
        order (int): order of the elastic tensor set to return
        strains (nx3x3 array-like): Array of 3x3 strains
            to use in fitting of ECs
        stresses (nx3x3 array-like): Array of 3x3 stresses
            to use in fitting ECs.  These should be PK2 stresses.
        eq_stress (3x3 array-like): stress corresponding to
            equilibrium strain (i. e. "0" strain state).
            If not specified, function will try to find
            the state in the list of provided stresses
            and strains.  If not found, defaults to 0.
        tol (float): value for which strains below
            are ignored in identifying strain states.

    Returns:
        Set of tensors corresponding to nth order expansion of
        the stress/strain relation
    """
    ustate_dict = get_ind_state_dict(
        independents, dependents, eq_state=eq_state, tol=tol,
        add_eq=True, sort=True)

    # Collect derivative data
    c_list = []
    dyi_dsi = np.zeros((order - 1, 6, len(ustate_dict)))
    for n, (strain_state, data) in enumerate(ustate_dict.items()):
        hvec = data["indeps"][:, strain_state.index(1)]
        for i in range(1, order):
            coef = get_diff_coeff(hvec, i)
            dei_dsi[i-1, :, n] = np.dot(coef, data["deps"])

    m, absent = generate_pseudo(list(strain_state_dict.keys()), order)
    for i in range(1, order):
        cvec, carr = get_symbol_list(i+1)
        svec = np.ravel(dei_dsi[i-1].T)
        cmap = dict(zip(cvec, np.dot(m[i-1], svec)))
        c_list.append(v_subs(carr, cmap))
    return [Tensor.from_voigt(c) for c in c_list]


def find_eq_stress(strains, stresses, tol=1e-10):
    """
    Finds stress corresponding to zero strain state in stress-strain list

    Args:
        strains (Nx3x3 array-like): array corresponding to strains
        stresses (Nx3x3 array-like): array corresponding to stresses
        tol (float): tolerance to find zero strain state
    """
    stress_array = np.array(stresses)
    strain_array = np.array(strains)
    eq_stress = stress_array[np.all(abs(strain_array)<tol, axis=(1,2))]

    if eq_stress.size != 0:
        all_same = (abs(eq_stress - eq_stress[0]) < 1e-8).all()
        if len(eq_stress) > 1 and not all_same:
            raise ValueError("Multiple stresses found for equilibrium strain"
                             " state, please specify equilibrium stress or  "
                             " remove extraneous stresses.")
        eq_stress = eq_stress[0]
    else:
        warnings.warn("No eq state found, returning zero voigt stress")
        eq_stress = Stress(np.zeros((3, 3)))
    return eq_stress


def get_ind_state_dict(independents, dependents, eq_state=None,
                       tol=1e-10, add_eq=True, sort=True):
    """
    Creates a dictionary of voigt-notation stress-strain sets 
    keyed by "strain state", i. e. a tuple corresponding to 
    the non-zero entries in ratios to the lowest nonzero value, 
    e.g. [0, 0.1, 0, 0.2, 0, 0] -> (0,1,0,2,0,0)
    This allows strains to be collected in stencils as to
    evaluate parameterized finite difference derivatives

    Args:
        strains (Nx3x3 array-like): strain matrices
        stresses (Nx3x3 array-like): stress matrices
        eq_stress (Nx3x3 array-like): equilibrium stress
        tol (float): tolerance for sorting strain states
        add_eq (bool): flag for whether to add eq_strain
            to stress-strain sets for each strain state
        sort (bool): flag for whether to sort strain states

    Returns:
        OrderedDict with strain state keys and dictionaries
        with stress-strain data corresponding to strain state 
    """
    # Recast inputs to voigt if not already
    if not all([isinstance(s, Tensor) for s in independents + dependents]):
        raise ValueError("Inputs must be lists of Tensor objects.  Cast "
                         "inputs to relevant tensor objects")
    vindeps = np.array([i.zeroed(tol).voigt for i in independents])
    vdeps = np.array([d.zeroed(tol).voigt for d in dependents])

    # Collect unique strain states:
    unique_states = set([tuple(np.nonzero(vindep)[0].tolist())
                         for vindep in vindeps])
    ustate_dict = OrderedDict()
    if add_eq:
        if eq_state is not None:
            veq_state = indeps[0].__class__(eq_state).voigt
        else:
            veq_state = find_eq_state(independents, dependents).voigt

    for n, us in enumerate(unique_states):
        # match strains with templates
        template = np.zeros(vindeps.shape, dtype=bool)
        np.put(template, us, True)
        template = np.tile(template, [vdeps.shape[0], 1])
        mode = (template == (np.abs(vindeps) > 1e-10)).all(axis=1)
        mindeps = vindeps[mode]
        mdeps = vdeps[mode]

        if add_eq:
            # add zero strain state
            mindeps = np.vstack([mindeps, np.zeros(mindeps.shape[1:])])
            mdeps = np.vstack([mdeps, veq_state])
        # sort strains/stresses by strain values
        if sort:
            mdeps = mdeps[mindeps[:, ind[0]].argsort()]
            mindeps = mindeps[mindeps[:, ind[0]].argsort()]
        # Get "strain state", i.e. ratio of each value to minimum strain
        ustate = mindeps[-1] / np.min(np.take(mindeps[-1], ind))
        ustate = tuple(state)
        ustate_dict[ustate] = {"indeps": mindeps,
                               "deps": mdeps}
    return ustate_dict


def generate_pseudo(ustates, shape_deps=6, order=2):
    """
    Generates the pseudoinverse for a given set of strains.
    
    Args:
        ustates (6xN array like): a list of voigt-notation 
            "strain-states", i. e. perturbed indices of the strain
            as a function of the smallest strain e. g. (0, 1, 0, 0, 1, 0)
        order (int): order of pseudoinverse to calculate
    
    Returns:
        mis: pseudo inverses for each order tensor, these can 
            be multiplied by the central difference derivative 
            of the stress with respect to the strain state
        absent_syms: symbols of the tensor absent from the PI
            expression
    """
    s = sp.Symbol('s')
    nstates = len(ustates)
    ni = np.array(ustates)*s
    mis, absent_syms = [], []
    for degree in range(2, order + 1):
        tvec, tarr = get_symbol_list(degree)
        sarr = np.zeros((nstates, len(ustates[0])), dtype=object)
        for n, us_vec in enumerate(ni):
            # Get expressions
            exps = carr.copy()
            for i in range(degree - 1):
                exps = np.dot(exps, us_vec)
            exps /= np.math.factorial(degree - 1)
            sarr[n] = [sp.diff(exp, s, degree - 1) for exp in exps]
        svec = sarr.ravel()
        present_syms = set.union(*[exp.atoms(sp.Symbol) for exp in svec])
        absent_syms += [set(cvec) - present_syms]
        m = np.zeros((6*nstates, len(cvec)))
        for n, c in enumerate(cvec):
            m[:, n] = v_diff(svec, c)
        mis.append(np.linalg.pinv(m))
    return mis, absent_syms


def get_symbol_list(shape=(6,6)):
    """
    Returns a symbolic representation of the voigt-notation
    tensor that places identical symbols for entries related
    by index transposition, i. e. C_1121 = C_1211 etc.

    Args:
        dim (int): dimension of matrix/tensor, e. g. 6 for
            voigt notation and 3 for standard
        rank (int): rank of tensor, e. g. 3 for third-order ECs

    Returns:
        t_vec (array): array representing distinct indices
        t_arr (array): array representing tensor with equivalent
            indices assigned as above
    """
    # Enforce symmetry if even rank
    if all([s==6 for s in shape]):
        indices = list(
            itertools.combinations_with_replacement(6, r=len(shape)))
    else:
        indices = itertools.product(*[range(s) for s in shape])
    t_vec = np.zeros(len(indices), dtype=object)
    t_arr = np.zeros(shape, dtype=object)
    for n, idx in enumerate(indices):
        t_vec[n] = sp.Symbol('c_'+''.join([str(i) for i in idx]))
        if symm:
            for perm in itertools.permutations(idx):
                t_arr[perm] = t_vec[n]
        else:
            t_arr[idx] = t_vec[n]
    return t_vec, t_arr


def subs(entry, tmap):
    """
    Sympy substitution function, primarily for the purposes
    of numpy vectorization

    Args:
        entry (symbol or exp): sympy expr to undergo subs
        tmap (dict): map for symbols to values to use in subs

    Returns:
        Evaluated expression with substitution
    """
    return entry.subs(tmap)

# Vectorized functions
v_subs = np.vectorize(subs)
v_diff = np.vectorize(sp.diff)


def get_diff_coeff(hvec, n=1):
    """
    Helper function to find difference coefficients of an
    derivative on an arbitrary mesh.

    Args:
        hvec (1D array-like): sampling stencil
        n (int): degree of derivative to find
    """
    hvec = np.array(hvec, dtype=np.float)
    acc = len(hvec)
    exp = np.column_stack([np.arange(acc)]*acc)
    a = np.vstack([hvec] * acc) ** exp
    b = np.zeros(acc)
    b[n] = factorial(n)
    return np.linalg.solve(a, b)
