"""This module provides a base class Tensor for tensor-like objects and
methods for basic tensor manipulation. It also provides SquareTensor,
which provides basic methods for creating and manipulating rank 2 tensors.
"""

from __future__ import annotations

import collections
import itertools
import os
import string
import warnings
from typing import TYPE_CHECKING

import numpy as np
from monty.json import MSONable
from monty.serialization import loadfn
from scipy.linalg import polar

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from numpy.typing import NDArray
    from typing_extensions import Self

    from pymatgen.core import Structure

__author__ = "Joseph Montoya"
__credits__ = "Maarten de Jong, Shyam Dwaraknath, Wei Chen, Mark Asta, Anubhav Jain, Terence Lew"


DEFAULT_QUAD = loadfn(os.path.join(os.path.dirname(__file__), "quad_data.json"))


class Tensor(np.ndarray, MSONable):
    """Base class for doing useful general operations on Nth order tensors,
    without restrictions on the type (stress, elastic, strain, piezo, etc.).
    """

    symbol = "T"

    def __new__(
        cls,
        input_array: NDArray,
        vscale: NDArray | None = None,
        check_rank: int | None = None,
    ) -> Self:
        """Create a Tensor object. Note that the constructor uses __new__
        rather than __init__ according to the standard method of
        subclassing numpy ndarrays.

        Args:
            input_array: (array-like with shape 3^N): array-like representing
                a tensor quantity in standard (i.e. non-Voigt) notation
            vscale: (N x M array-like): a matrix corresponding
                to the coefficients of the Voigt-notation tensor
            check_rank: (int): If not None, checks that input_array's rank == check_rank.
                Defaults to None.
        """
        obj = np.asarray(input_array).view(cls)
        obj.rank = len(obj.shape)

        if check_rank and check_rank != obj.rank:
            raise ValueError(f"{type(obj).__name__} input must be rank {check_rank}")

        vshape = tuple([3] * (obj.rank % 2) + [6] * (obj.rank // 2))
        obj._vscale = np.ones(vshape)
        if vscale is not None:
            obj._vscale = vscale
        if obj._vscale.shape != vshape:
            raise ValueError("Voigt scaling matrix must be the shape of the Voigt notation matrix or vector.")
        if any(dim != 3 for dim in obj.shape):
            raise ValueError(
                "Pymatgen only supports 3-dimensional tensors, and default tensor constructor uses standard "
                f"notation. To construct from Voigt notation, use {type(obj).__name__}.from_voigt"
            )
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.rank = getattr(obj, "rank", None)
        self._vscale = getattr(obj, "_vscale", None)
        self._vdict = getattr(obj, "_vdict", None)

    def __array_wrap__(self, obj):
        """Overrides __array_wrap__ methods in ndarray superclass to avoid errors
        associated with functions that return scalar values.
        """
        if len(obj.shape) == 0:
            return obj[()]
        return np.ndarray.__array_wrap__(self, obj)

    def __hash__(self) -> int:
        """Define a hash function, since numpy arrays have their own __eq__ method."""
        return hash(self.tostring())

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self})"

    def zeroed(self, tol: float = 1e-3) -> Self:
        """Get the matrix with all entries below a certain threshold (i.e. tol) set to zero."""
        new_tensor = self.copy()
        new_tensor[abs(new_tensor) < tol] = 0
        return new_tensor

    def transform(self, symm_op: SymmOp) -> Self:
        """Apply a transformation (via a symmetry operation) to a tensor.

        Args:
            symm_op (SymmOp): a symmetry operation to apply to the tensor
        """
        return type(self)(symm_op.transform_tensor(self))

    def rotate(self, matrix: NDArray, tol: float = 1e-3) -> Self:
        """Apply a rotation directly, and tests input matrix to ensure a valid
        rotation.

        Args:
            matrix (3x3 array-like): rotation matrix to be applied to tensor
            tol (float): tolerance for testing rotation matrix validity
        """
        matrix = SquareTensor(matrix)
        if not matrix.is_rotation(tol):
            raise ValueError("Rotation matrix is not valid.")
        symm_op = SymmOp.from_rotation_and_translation(matrix, [0.0, 0.0, 0.0])
        return self.transform(symm_op)

    def einsum_sequence(
        self,
        other_arrays: NDArray,
        einsum_string: str | None = None,
    ) -> NDArray:
        """Calculate the result of an einstein summation expression."""
        if not isinstance(other_arrays, list):
            raise TypeError("other tensors must be list of tensors or tensor input")

        other_arrays = [np.array(a) for a in other_arrays]
        if not einsum_string:
            lc = string.ascii_lowercase
            einsum_string = lc[: self.rank]
            other_ranks = [len(a.shape) for a in other_arrays]
            idx = self.rank - sum(other_ranks)
            for length in other_ranks:
                einsum_string += f",{lc[idx : idx + length]}"
                idx += length

        einsum_args = [self, *other_arrays]
        return np.einsum(einsum_string, *einsum_args)

    def project(self, n: NDArray) -> Self:
        """Project a tensor into a vector. Returns the tensor
        dotted into a unit vector along the input n.

        Args:
            n (3x1 array-like): direction to project onto

        Returns:
            float: scalar value corresponding to the projection of
                the tensor into the vector
        """
        unit_vec = get_uvec(n)
        return self.einsum_sequence([unit_vec] * self.rank)

    def average_over_unit_sphere(self, quad: dict | None = None) -> Self:
        """Average the tensor projection over the unit with option for custom quadrature.

        Args:
            quad (dict): quadrature for integration, should be
                dictionary with "points" and "weights" keys defaults
                to quadpy.sphere.Lebedev(19) as read from file

        Returns:
            Average of tensor projected into vectors on the unit sphere
        """
        quad = quad or DEFAULT_QUAD
        weights, points = quad["weights"], quad["points"]
        return sum(w * self.project(n) for w, n in zip(weights, points, strict=True))

    def get_grouped_indices(self, voigt: bool = False, **kwargs) -> list[list]:
        """Get index sets for equivalent tensor values.

        Args:
            voigt (bool): whether to get grouped indices
                of voigt or full notation tensor, defaults
                to false
            **kwargs: keyword args for np.isclose. Can take atol
                and rtol for absolute and relative tolerance, e.g.

                >>> tensor.group_array_indices(atol=1e-8)

                or

                >>> tensor.group_array_indices(rtol=1e-5)

        Returns:
            list of index groups where tensor values are equivalent to
            within tolerances
        """
        array = self.voigt if voigt else self

        indices = list(itertools.product(*(range(n) for n in array.shape)))
        remaining = indices.copy()
        # Start with everything near zero
        grouped = [list(zip(*np.where(np.isclose(array, 0, **kwargs)), strict=True))]
        remaining = [i for i in remaining if i not in grouped[0]]
        # Iteratively run through remaining indices
        while remaining:
            new = list(zip(*np.where(np.isclose(array, array[remaining[0]], **kwargs)), strict=True))
            grouped.append(new)
            remaining = [i for i in remaining if i not in new]
        # Don't return any empty lists
        return [g for g in grouped if g]

    def get_symbol_dict(
        self,
        voigt: bool = True,
        zero_index: bool = False,
        **kwargs,
    ) -> dict[str, NDArray]:
        """Create a summary dict for tensor with associated symbol.

        Args:
            voigt (bool): whether to get symbol dict for voigt
                notation tensor, as opposed to full notation,
                defaults to true
            zero_index (bool): whether to set initial index to zero,
                defaults to false, since tensor notations tend to use
                one-indexing, rather than zero indexing like python
            **kwargs: keyword args for np.isclose. Can take atol
                and rtol for absolute and relative tolerance, e.g.

                >>> tensor.get_symbol_dict(atol=1e-8)

                or

                >>> tensor.get_symbol_dict(rtol=1e-5)

        Returns:
            list of index groups where tensor values are equivalent to
            within tolerances
        """
        dct = {}
        array = self.voigt if voigt else self
        grouped = self.get_grouped_indices(voigt=voigt, **kwargs)
        p = 0 if zero_index else 1
        for indices in grouped:
            sym_string = f"{self.symbol}_"
            sym_string += "".join(str(i + p) for i in indices[0])
            value = array[indices[0]]
            if not np.isclose(value, 0):
                dct[sym_string] = array[indices[0]]
        return dct

    def round(self, decimals: int = 0) -> Self:
        """Wrapper around numpy.round to ensure object
        of same type is returned.

        Args:
            decimals: Number of decimal places to round to (default: 0).
                If decimals is negative, it specifies the number of
                positions to the left of the decimal point.

        Returns:
            Tensor: rounded tensor of same type
        """
        return type(self)(np.round(self, decimals=decimals))

    @property
    def symmetrized(self) -> Self:
        """A generally symmetrized tensor, calculated by taking
        the sum of the tensor and its transpose with respect to all
        possible permutations of indices.
        """
        perms = list(itertools.permutations(range(self.rank)))
        return sum(np.transpose(self, ind) for ind in perms) / len(perms)

    @property
    def voigt_symmetrized(self) -> Self:
        """A "voigt"-symmetrized tensor, i.e. a Voigt-notation
        tensor such that it is invariant w.r.t. permutation of indices.
        """
        if self.rank % 2 != 0 or self.rank < 2:
            raise ValueError("V-symmetrization requires rank even and >= 2")

        v = self.voigt
        perms = list(itertools.permutations(range(len(v.shape))))
        new_v = sum(np.transpose(v, ind) for ind in perms) / len(perms)
        return type(self).from_voigt(new_v)

    def is_symmetric(self, tol: float = 1e-5) -> bool:
        """Test whether a tensor is symmetric or not based on the residual
        with its symmetric part, from self.symmetrized.

        Args:
            tol (float): tolerance to test for symmetry
        """
        return (self - self.symmetrized < tol).all()

    def fit_to_structure(
        self,
        structure: Structure,
        symprec: float = 0.1,
    ):
        """Get a tensor that is invariant with respect to symmetry
        operations corresponding to a structure.

        Args:
            structure (Structure): structure from which to generate
                symmetry operations
            symprec (float): symmetry tolerance for the Spacegroup Analyzer
                used to generate the symmetry operations
        """
        sga = SpacegroupAnalyzer(structure, symprec)
        symm_ops = sga.get_symmetry_operations(cartesian=True)
        return sum(self.transform(symm_op) for symm_op in symm_ops) / len(symm_ops)

    def is_fit_to_structure(self, structure: Structure, tol: float = 1e-2) -> bool:
        """Test whether a tensor is invariant with respect to the
        symmetry operations of a particular structure by testing
        whether the residual of the symmetric portion is below a
        tolerance.

        Args:
            structure (Structure): structure to be fit to
            tol (float): tolerance for symmetry testing
        """
        return (self - self.fit_to_structure(structure) < tol).all()

    @property
    def voigt(self) -> NDArray:
        """The tensor in Voigt notation."""
        v_matrix = np.zeros(self._vscale.shape, dtype=self.dtype)
        this_voigt_map = self.get_voigt_dict(self.rank)
        for ind, v in this_voigt_map.items():
            v_matrix[v] = self[ind]
        if not self.is_voigt_symmetric():
            warnings.warn("Tensor is not symmetric, information may be lost in Voigt conversion.")
        return v_matrix * self._vscale

    def is_voigt_symmetric(self, tol: float = 1e-6) -> bool:
        """Test symmetry of tensor to that necessary for voigt-conversion
        by grouping indices into pairs and constructing a sequence of
        possible permutations to be used in a tensor transpose.
        """
        transpose_pieces = [[[0 for _ in range(self.rank % 2)]]]
        transpose_pieces += [[list(range(j, j + 2))] for j in range(self.rank % 2, self.rank, 2)]
        for n in range(self.rank % 2, len(transpose_pieces)):
            if len(transpose_pieces[n][0]) == 2:
                transpose_pieces[n] += [transpose_pieces[n][0][::-1]]
        for trans_seq in itertools.product(*transpose_pieces):
            transpose_seq = list(itertools.chain(*trans_seq))
            if (self - self.transpose(transpose_seq) > tol).any():
                return False
        return True

    @staticmethod
    def get_voigt_dict(rank: int) -> dict[tuple[int, ...], tuple[int, ...]]:
        """Get a dictionary that maps indices in the tensor to those
        in a voigt representation based on input rank.

        Args:
            rank (int): Tensor rank to generate the voigt map
        """
        reverse_voigt_map = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])

        voigt_dict = {}
        for ind in itertools.product(*[range(3)] * rank):
            v_ind = ind[: rank % 2]
            for j in range(rank // 2):
                pos = rank % 2 + 2 * j
                v_ind += (reverse_voigt_map[ind[pos : pos + 2]],)
            voigt_dict[ind] = v_ind
        return voigt_dict

    @classmethod
    def from_voigt(cls, voigt_input: NDArray) -> Self:
        """Constructor based on the voigt notation vector or matrix.

        Args:
            voigt_input (array-like): voigt input for a given tensor
        """
        voigt_input = np.array(voigt_input)
        rank = sum(voigt_input.shape) // 3
        t = cls(np.zeros([3] * rank))
        if voigt_input.shape != t._vscale.shape:
            raise ValueError("Invalid shape for Voigt matrix")
        voigt_input = voigt_input / t._vscale  # (ruff-preview) noqa: PLR6104
        this_voigt_map = t.get_voigt_dict(rank)
        for ind, v in this_voigt_map.items():
            t[ind] = voigt_input[v]
        return cls(t)

    @staticmethod
    def get_ieee_rotation(
        structure: Structure,
        refine_rotation: bool = True,
    ) -> SquareTensor:
        """Given a structure associated with a tensor, determines
        the rotation matrix for IEEE conversion according to
        the 1987 IEEE standards.

        Args:
            structure (Structure): a structure associated with the
                tensor to be converted to the IEEE standard
            refine_rotation (bool): whether to refine the rotation
                using SquareTensor.refine_rotation
        """
        # Check conventional setting:
        sga = SpacegroupAnalyzer(structure)
        dataset = sga.get_symmetry_dataset()
        trans_mat = dataset.transformation_matrix
        conv_latt = Lattice(np.transpose(np.dot(np.transpose(structure.lattice.matrix), np.linalg.inv(trans_mat))))
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
            rotation = np.array(
                [vec / mag for (mag, vec) in sorted(zip(lengths, vecs, strict=True), key=lambda x: x[0])]
            )
            if abs(lengths[2] - lengths[1]) < abs(lengths[1] - lengths[0]):
                rotation[0], rotation[2] = rotation[2], rotation[0].copy()
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[0]))

        # IEEE rules: c<a<b; c,a || x3,x1
        elif xtal_sys == "orthorhombic":
            rotation = [vec / mag for (mag, vec) in sorted(zip(lengths, vecs, strict=True))]
            rotation = np.roll(rotation, 2, axis=0)

        # IEEE rules: c,a || x3,x1, c is threefold axis
        # Note this also includes rhombohedral crystal systems
        elif xtal_sys in ("trigonal", "hexagonal"):
            # find threefold axis:
            tf_index = np.argmin(abs(angles - 120.0))
            non_tf_mask = np.logical_not(angles == angles[tf_index])
            rotation[2] = get_uvec(vecs[tf_index])
            rotation[0] = get_uvec(vecs[non_tf_mask][0])
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[0]))

        # IEEE rules: b,c || x2,x3; alpha=beta=90, c<a
        elif xtal_sys == "monoclinic":
            # Find unique axis
            u_index = np.argmax(abs(angles - 90.0))
            n_umask = np.logical_not(angles == angles[u_index])
            rotation[1] = get_uvec(vecs[u_index])
            # Shorter of remaining lattice vectors for c axis
            c = next(vec / mag for (mag, vec) in sorted(zip(lengths[n_umask], vecs[n_umask], strict=True)))
            rotation[2] = np.array(c)
            rotation[0] = np.cross(rotation[1], rotation[2])

        # IEEE rules: c || x3, x2 normal to ac plane
        elif xtal_sys == "triclinic":
            rotation = [vec / mag for (mag, vec) in sorted(zip(lengths, vecs, strict=True))]
            rotation[1] = get_uvec(np.cross(rotation[2], rotation[0]))
            rotation[0] = np.cross(rotation[1], rotation[2])

        rotation = SquareTensor(rotation)
        if refine_rotation:
            rotation = rotation.refine_rotation()

        return rotation

    def convert_to_ieee(
        self,
        structure: Structure,
        initial_fit: bool = True,
        refine_rotation: bool = True,
    ) -> Self:
        """Given a structure associated with a tensor, attempts a
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
            refine_rotation (bool): whether to refine the rotation
                produced by the ieee transform generator, default True
        """
        rotation = self.get_ieee_rotation(structure, refine_rotation)
        result = self.copy()
        if initial_fit:
            result = result.fit_to_structure(structure)
        return result.rotate(rotation, tol=1e-2)

    def structure_transform(
        self,
        original_structure: Structure,
        new_structure: Structure,
        refine_rotation: bool = True,
    ) -> Self:
        """Transforms a tensor from one basis for an original structure
        into a new basis defined by a new structure.

        Args:
            original_structure (Structure): structure corresponding
                to the basis of the current tensor
            new_structure (Structure): structure corresponding to the
                desired basis
            refine_rotation (bool): whether to refine the rotations
                generated in get_ieee_rotation

        Returns:
            Tensor that has been transformed such that its basis
            corresponds to the new_structure's basis
        """
        sm = StructureMatcher()
        if not sm.fit(original_structure, new_structure):
            warnings.warn("original and new structures do not match!")
        trans_1 = self.get_ieee_rotation(original_structure, refine_rotation)
        trans_2 = self.get_ieee_rotation(new_structure, refine_rotation)
        # Get the ieee format tensor
        new = self.rotate(trans_1)
        # Reverse the ieee format rotation for the second structure
        return new.rotate(np.transpose(trans_2))

    @classmethod
    def from_values_indices(
        cls,
        values: list[float],
        indices: NDArray,
        populate: bool = False,
        structure: Structure | None = None,
        voigt_rank: int | None = None,
        vsym: bool = True,
        verbose: bool = False,
    ) -> Self:
        """Create a tensor from values and indices, with options
        for populating the remainder of the tensor.

        Args:
            values (floats): numbers to place at indices
            indices (array-likes): indices to place values at
            populate (bool): whether to populate the tensor
            structure (Structure): structure to base population
                or fit_to_structure on
            voigt_rank (int): full tensor rank to indicate the
                shape of the resulting tensor. This is necessary
                if one provides a set of indices more minimal than
                the shape of the tensor they want, e.g.
                Tensor.from_values_indices((0, 0), 100)
            vsym (bool): whether to voigt symmetrize during the
                optimization procedure
            verbose (bool): whether to populate verbosely
        """
        # Auto-detect voigt notation
        # TODO: refactor rank inheritance to make this easier
        indices = np.array(indices)
        if voigt_rank:
            shape = np.array([3] * (voigt_rank % 2) + [6] * (voigt_rank // 2))
        else:
            shape = np.ceil(np.max(indices + 1, axis=0) / 3.0) * 3

        base = np.zeros(shape.astype(int))
        for v, idx in zip(values, indices, strict=True):
            base[tuple(idx)] = v
        obj = cls.from_voigt(base) if 6 in shape else cls(base)

        if populate:
            if not structure:
                raise ValueError("Populate option must include structure input")
            obj = obj.populate(structure, vsym=vsym, verbose=verbose)
        elif structure:
            obj = obj.fit_to_structure(structure)

        return obj

    def populate(
        self,
        structure: Structure,
        prec: float = 1e-5,
        maxiter: int = 200,
        verbose: bool = False,
        precond: bool = True,
        vsym: bool = True,
    ) -> Self:
        """Takes a partially populated tensor, and populates the non-zero
        entries according to the following procedure, iterated until
        the desired convergence (specified via prec) is achieved.

        1. Find non-zero entries
        2. Symmetrize the tensor with respect to crystal symmetry and
           (optionally) voigt symmetry
        3. Reset the non-zero entries of the original tensor

        Args:
            structure (Structure): structure to base population on
            prec (float): precision for determining a non-zero value. Defaults to 1e-5.
            maxiter (int): maximum iterations for populating the tensor
            verbose (bool): whether to populate verbosely
            precond (bool): whether to precondition by cycling through
                all symmops and storing new nonzero values, default True
            vsym (bool): whether to enforce voigt symmetry, defaults
                to True

        Returns:
            Tensor: Populated tensor
        """
        guess = type(self)(np.zeros(self.shape))
        mask = None
        if precond:
            # Generate the guess from populated
            sops = SpacegroupAnalyzer(structure).get_symmetry_operations()
            mask = abs(self) > prec
            guess[mask] = self[mask]

            def merge(old, new) -> None:
                gmask = np.abs(old) > prec
                nmask = np.abs(new) > prec
                new_mask = np.logical_not(gmask) * nmask
                avg_mask = gmask * nmask
                old[avg_mask] = (old[avg_mask] + new[avg_mask]) / 2.0
                old[new_mask] = new[new_mask]

            if verbose:
                print(f"Preconditioning for {len(sops)} symmops")
            for sop in sops:
                rot = guess.transform(sop)
                # Store non-zero entries of new that weren't previously
                # in the guess in the guess
                merge(guess, rot)
            if verbose:
                print("Preconditioning for voigt symmetry")
            if vsym:
                v = guess.voigt
                perms = list(itertools.permutations(range(len(v.shape))))
                for perm in perms:
                    vtrans = np.transpose(v, perm)
                    merge(v, vtrans)
                guess = type(self).from_voigt(v)

        if guess.shape != self.shape:
            raise ValueError("Guess must have same shape")
        converged = False
        test_new, test_old = [guess.copy()] * 2
        for idx in range(maxiter):
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
                print(f"Iteration {idx}: {np.max(diff)}")
        if not converged:
            max_diff = np.max(np.abs(self - test_new))
            warnings.warn(f"Warning, populated tensor is not converged with max diff of {max_diff}")
        return type(self)(test_new)

    def as_dict(self, voigt: bool = False) -> dict:
        """Serializes the tensor object.

        Args:
            voigt (bool): flag for whether to store entries in Voigt notation.
                Defaults to false, as information may be lost in conversion.

        Returns:
            dict: serialized format tensor object
        """
        input_array = self.voigt if voigt else self
        dct = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "input_array": input_array.tolist(),
        }
        if voigt:
            dct["voigt"] = voigt
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Instantiate Tensors from dicts (using MSONable API).

        Returns:
            Tensor: hydrated tensor object
        """
        if dct.get("voigt"):
            return cls.from_voigt(dct["input_array"])
        return cls(dct["input_array"])


class TensorCollection(collections.abc.Sequence, MSONable):
    """A sequence of tensors that can be used for fitting data
    or for having a tensor expansion.
    """

    def __init__(self, tensor_list: Sequence, base_class=Tensor) -> None:
        """
        Args:
            tensor_list: List of tensors.
            base_class: Class to be used.
        """
        self.tensors = [tensor if isinstance(tensor, base_class) else base_class(tensor) for tensor in tensor_list]

    def __len__(self) -> int:
        return len(self.tensors)

    def __getitem__(self, ind):
        return self.tensors[ind]

    def __iter__(self):
        return iter(self.tensors)

    def zeroed(self, tol: float = 1e-3) -> Self:
        """
        Args:
            tol: Tolerance.

        Returns:
            TensorCollection where small values are set to 0.
        """
        return type(self)([tensor.zeroed(tol) for tensor in self])

    def transform(self, symm_op: SymmOp) -> Self:
        """Transforms TensorCollection with a symmetry operation.

        Args:
            symm_op: SymmetryOperation.

        Returns:
            TensorCollection.
        """
        return type(self)([tensor.transform(symm_op) for tensor in self])

    def rotate(self, matrix, tol: float = 1e-3) -> Self:
        """Rotates TensorCollection.

        Args:
            matrix: Rotation matrix.
            tol: tolerance.

        Returns:
            TensorCollection.
        """
        return type(self)([tensor.rotate(matrix, tol) for tensor in self])

    @property
    def symmetrized(self) -> Self:
        """TensorCollection where all tensors are symmetrized."""
        return type(self)([tensor.symmetrized for tensor in self])

    def is_symmetric(self, tol: float = 1e-5) -> bool:
        """
        Args:
            tol: tolerance.

        Returns:
            Whether all tensors are symmetric.
        """
        return all(tensor.is_symmetric(tol) for tensor in self)

    def fit_to_structure(
        self,
        structure: Structure,
        symprec: float = 0.1,
    ) -> Self:
        """Fit all tensors to a Structure.

        Args:
            structure: Structure
            symprec: symmetry precision.

        Returns:
            TensorCollection.
        """
        return type(self)([tensor.fit_to_structure(structure, symprec) for tensor in self])

    def is_fit_to_structure(
        self,
        structure: Structure,
        tol: float = 1e-2,
    ) -> bool:
        """
        Args:
            structure: Structure
            tol: tolerance.

        Returns:
            Whether all tensors are fitted to Structure.
        """
        return all(tensor.is_fit_to_structure(structure, tol) for tensor in self)

    @property
    def voigt(self) -> list[NDArray]:
        """TensorCollection where all tensors are in Voigt form."""
        return [tensor.voigt for tensor in self]

    @property
    def ranks(self) -> list:
        """Ranks for all tensors."""
        return [tensor.rank for tensor in self]

    def is_voigt_symmetric(self, tol: float = 1e-6) -> bool:
        """
        Args:
            tol: tolerance.

        Returns:
            Whether all tensors are voigt symmetric.
        """
        return all(tensor.is_voigt_symmetric(tol) for tensor in self)

    @classmethod
    def from_voigt(
        cls,
        voigt_input_list: list[Tensor],
        base_class=Tensor,
    ) -> Self:
        """Create TensorCollection from voigt form.

        Args:
            voigt_input_list: List of voigt tensors
            base_class: Class for tensor.

        Returns:
            TensorCollection.
        """
        return cls([base_class.from_voigt(v) for v in voigt_input_list])

    def convert_to_ieee(
        self,
        structure: Structure,
        initial_fit: bool = True,
        refine_rotation: bool = True,
    ) -> Self:
        """Convert all tensors to IEEE.

        Args:
            structure: Structure
            initial_fit: Whether to perform an initial fit.
            refine_rotation: Whether to refine the rotation.

        Returns:
            TensorCollection.
        """
        return type(self)([tensor.convert_to_ieee(structure, initial_fit, refine_rotation) for tensor in self])

    def round(self, *args, **kwargs) -> Self:
        """Round all tensors.

        Args:
            args: Passthrough to Tensor.round
            kwargs: Passthrough to Tensor.round

        Returns:
            TensorCollection.
        """
        return type(self)([tensor.round(*args, **kwargs) for tensor in self])

    @property
    def voigt_symmetrized(self) -> Self:
        """TensorCollection where all tensors are voigt symmetrized."""
        return type(self)([tensor.voigt_symmetrized for tensor in self])

    def as_dict(self, voigt: bool = False) -> dict:
        """
        Args:
            voigt: Whether to use Voigt form.

        Returns:
            Dict representation of TensorCollection.
        """
        tensor_list = self.voigt if voigt else self
        dct: dict[str, Any] = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "tensor_list": [tensor.tolist() for tensor in tensor_list],
        }
        if voigt:
            dct["voigt"] = voigt
        return dct

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """Create TensorCollection from dict.

        Args:
            dct: dict

        Returns:
            TensorCollection
        """
        if dct.get("voigt"):
            return cls.from_voigt(dct["tensor_list"])
        return cls(dct["tensor_list"])


class SquareTensor(Tensor):
    """Base class for doing useful general operations on second rank tensors
    (stress, strain etc.).
    """

    def __new__(
        cls,
        input_array: NDArray,
        vscale: NDArray | None = None,
    ) -> Self:
        """Create a SquareTensor object. Note that the constructor uses __new__ rather than
        __init__ according to the standard method of subclassing numpy ndarrays. Error
        is thrown when the class is initialized with non-square matrix.

        Args:
            input_array (3x3 array-like): the 3x3 array-like
                representing the content of the tensor
            vscale (6x1 array-like): 6x1 array-like scaling the
                Voigt-notation vector with the tensor entries
        """
        obj = super().__new__(cls, input_array, vscale, check_rank=2)
        return obj.view(cls)

    @property
    def trans(self) -> Self:
        """Shorthand for transpose on SquareTensor."""
        return type(self)(np.transpose(self))

    @property
    def inv(self) -> Self:
        """Shorthand for matrix inverse on SquareTensor."""
        if self.det == 0:
            raise ValueError("SquareTensor is non-invertible")
        return type(self)(np.linalg.inv(self))

    @property
    def det(self) -> Self:
        """Shorthand for the determinant of the SquareTensor."""
        return np.linalg.det(self)

    def is_rotation(
        self,
        tol: float = 1e-3,
        include_improper: bool = True,
    ) -> bool:
        """Test to see if tensor is a valid rotation matrix, performs a
        test to check whether the inverse is equal to the transpose
        and if the determinant is equal to one within the specified
        tolerance.

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
        return (np.abs(self.inv - self.trans) < tol).all() and (np.abs(det - 1.0) < tol)

    def refine_rotation(self) -> Self:
        """Helper method for refining rotation matrix by ensuring
        that second and third rows are perpendicular to the first.
        Gets new y vector from an orthogonal projection of x onto y
        and the new z vector from a cross product of the new x and y.

        Args:
            tol to test for rotation

        Returns:
            new rotation matrix
        """
        new_x, y = get_uvec(self[0]), get_uvec(self[1])
        # Get a projection on y
        new_y = y - np.dot(new_x, y) * new_x
        new_z = np.cross(new_x, new_y)
        return type(self)([new_x, new_y, new_z])

    def get_scaled(self, scale_factor: float) -> Self:
        """Scales the tensor by a certain multiplicative scale factor.

        Args:
            scale_factor (float): scalar multiplier to be applied to the
                SquareTensor object
        """
        return type(self)(self * scale_factor)

    @property
    def principal_invariants(self) -> NDArray:
        """A list of principal invariants for the tensor,
        which are the values of the coefficients of the characteristic
        polynomial for the matrix.
        """
        return np.poly(self)[1:] * np.array([-1, 1, -1])

    def polar_decomposition(self, side: str = "right") -> tuple:
        """Calculate matrices for polar decomposition."""
        return polar(self, side=side)


def get_uvec(vec: NDArray) -> NDArray:
    """Get a unit vector parallel to input vector."""
    norm = np.linalg.norm(vec)
    return vec if norm < 1e-8 else vec / norm


def symmetry_reduce(
    tensors,
    structure: Structure,
    tol: float = 1e-8,
    **kwargs,
) -> TensorMapping:
    """Convert a list of tensors corresponding to a structure
    and returns a dictionary consisting of unique tensor keys with SymmOp
    values corresponding to transformations that will result in derivative
    tensors from the original list.

    Args:
        tensors (list of tensors): list of Tensor objects to test for
            symmetrically-equivalent duplicates
        structure (Structure): structure from which to get symmetry
        tol (float): tolerance for tensor equivalence
        kwargs: keyword arguments for the SpacegroupAnalyzer

    Returns:
        dictionary consisting of unique tensors with symmetry operations
        corresponding to those which will reconstruct the remaining
        tensors as values
    """
    sga = SpacegroupAnalyzer(structure, **kwargs)
    symm_ops = sga.get_symmetry_operations(cartesian=True)
    unique_mapping = TensorMapping([tensors[0]], [[]], tol=tol)
    for tensor in tensors[1:]:
        is_unique = True
        for unique_tensor, symm_op in itertools.product(unique_mapping, symm_ops):
            if np.allclose(unique_tensor.transform(symm_op), tensor, atol=tol):
                unique_mapping[unique_tensor].append(symm_op)
                is_unique = False
                break
        if is_unique:
            unique_mapping[tensor] = []
    return unique_mapping


class TensorMapping(collections.abc.MutableMapping):
    """Base class for tensor mappings, which function much like
    a dictionary, but use numpy routines to determine approximate
    equality to keys for getting and setting items.

    This is intended primarily for convenience with things like
    stress-strain pairs and fitting data manipulation. In general,
    it is significantly less robust than a typical hashing
    and should be used with care.
    """

    def __init__(
        self,
        tensors: Sequence[Tensor] = (),
        values: Sequence = (),
        tol: float = 1e-5,
    ) -> None:
        """Initialize a TensorMapping.

        Args:
            tensors (Sequence[Tensor], optional): Defaults to (,).
            values (Sequence, optional): Values to be associated with tensors. Defaults to (,).
            tol (float, optional): an absolute tolerance for getting and setting items in the mapping.
                Defaults to 1e-5.

        Raises:
            ValueError: if tensors and values are not the same length
        """
        if len(values) != len(tensors):
            raise ValueError("TensorMapping must be initialized with tensors and values of equivalent length")
        self._tensor_list = list(tensors)  # needs to be a list
        self._value_list = list(values)  # needs to be a list
        self.tol = tol

    def __getitem__(self, item):
        index = self._get_item_index(item)
        if index is None:
            raise KeyError(f"{item} not found in mapping.")
        return self._value_list[index]

    def __setitem__(self, key, value) -> None:
        index = self._get_item_index(key)
        if index is None:
            self._tensor_list.append(key)
            self._value_list.append(value)
        else:
            self._value_list[index] = value

    def __delitem__(self, key) -> None:
        index = self._get_item_index(key)
        self._tensor_list.pop(index)
        self._value_list.pop(index)

    def __len__(self) -> int:
        return len(self._tensor_list)

    def __iter__(self):
        yield from self._tensor_list

    def __contains__(self, item) -> bool:
        return self._get_item_index(item) is not None

    def values(self):
        """Values in mapping."""
        return self._value_list

    def items(self):
        """Items in mapping."""
        return zip(self._tensor_list, self._value_list, strict=True)

    def _get_item_index(self, item):
        if len(self._tensor_list) == 0:
            return None
        item = np.array(item)
        axis = tuple(range(1, len(item.shape) + 1))
        mask = np.all(np.abs(np.array(self._tensor_list) - item) < self.tol, axis=axis)
        indices = np.where(mask)[0]
        if len(indices) > 1:
            raise ValueError("Tensor key collision.")

        return None if len(indices) == 0 else indices[0]
