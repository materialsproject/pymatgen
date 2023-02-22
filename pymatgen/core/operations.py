# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes that operate on points or vectors in 3D space.
"""

from __future__ import annotations

import re
import string
import warnings
from math import cos, pi, sin, sqrt
from typing import Any

import numpy as np
from monty.json import MSONable

from pymatgen.electronic_structure.core import Magmom
from pymatgen.util.string import transformation_to_string
from pymatgen.util.typing import ArrayLike

__author__ = "Shyue Ping Ong, Shyam Dwaraknath, Matthew Horton"


class SymmOp(MSONable):
    """
    A symmetry operation in Cartesian space. Consists of a rotation plus a
    translation. Implementation is as an affine transformation matrix of rank 4
    for efficiency. Read: http://en.wikipedia.org/wiki/Affine_transformation.

    .. attribute:: affine_matrix

        A 4x4 numpy.array representing the symmetry operation.
    """

    def __init__(self, affine_transformation_matrix: ArrayLike, tol: float = 0.01) -> None:
        """
        Initializes the SymmOp from a 4x4 affine transformation matrix.
        In general, this constructor should not be used unless you are
        transferring rotations. Use the static constructors instead to
        generate a SymmOp from proper rotations and translation.

        Args:
            affine_transformation_matrix (4x4 array): Representing an
                affine transformation.
            tol (float): Tolerance for determining if matrices are equal.
        """
        affine_transformation_matrix = np.array(affine_transformation_matrix)
        if affine_transformation_matrix.shape != (4, 4):
            raise ValueError("Affine Matrix must be a 4x4 numpy array!")
        self.affine_matrix = affine_transformation_matrix
        self.tol = tol

    @staticmethod
    def from_rotation_and_translation(
        rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        translation_vec: ArrayLike = (0, 0, 0),
        tol: float = 0.1,
    ) -> SymmOp:
        """
        Creates a symmetry operation from a rotation matrix and a translation
        vector.

        Args:
            rotation_matrix (3x3 array): Rotation matrix.
            translation_vec (3x1 array): Translation vector.
            tol (float): Tolerance to determine if rotation matrix is valid.

        Returns:
            SymmOp object
        """
        rotation_matrix = np.array(rotation_matrix)
        translation_vec = np.array(translation_vec)
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Rotation Matrix must be a 3x3 numpy array.")
        if translation_vec.shape != (3,):
            raise ValueError("Translation vector must be a rank 1 numpy array with 3 elements.")
        affine_matrix = np.eye(4)
        affine_matrix[0:3][:, 0:3] = rotation_matrix
        affine_matrix[0:3][:, 3] = translation_vec
        return SymmOp(affine_matrix, tol)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SymmOp):
            return NotImplemented
        return np.allclose(self.affine_matrix, other.affine_matrix, atol=self.tol)

    def __hash__(self) -> int:
        return 7

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        output = [
            "Rot:",
            str(self.affine_matrix[0:3][:, 0:3]),
            "tau",
            str(self.affine_matrix[0:3][:, 3]),
        ]
        return "\n".join(output)

    def operate(self, point: ArrayLike) -> np.ndarray:
        """
        Apply the operation on a point.

        Args:
            point: Cartesian coordinate.

        Returns:
            Coordinates of point after operation.
        """
        affine_point = np.array([*point, 1])  # type: ignore
        return np.dot(self.affine_matrix, affine_point)[0:3]

    def operate_multi(self, points: ArrayLike) -> np.ndarray:
        """
        Apply the operation on a list of points.

        Args:
            points: List of Cartesian coordinates

        Returns:
            Numpy array of coordinates after operation
        """
        points = np.array(points)
        affine_points = np.concatenate([points, np.ones(points.shape[:-1] + (1,))], axis=-1)
        return np.inner(affine_points, self.affine_matrix)[..., :-1]

    def apply_rotation_only(self, vector: ArrayLike) -> np.ndarray:
        """
        Vectors should only be operated by the rotation matrix and not the
        translation vector.

        Args:
            vector (3x1 array): A vector.
        """
        return np.dot(self.rotation_matrix, vector)

    def transform_tensor(self, tensor: np.ndarray) -> np.ndarray:
        """
        Applies rotation portion to a tensor. Note that tensor has to be in
        full form, not the Voigt form.

        Args:
            tensor (numpy array): a rank n tensor

        Returns:
            Transformed tensor.
        """
        dim = tensor.shape
        rank = len(dim)
        assert all(i == 3 for i in dim)
        # Build einstein sum string
        lc = string.ascii_lowercase
        indices = lc[:rank], lc[rank : 2 * rank]
        einsum_string = ",".join(a + i for a, i in zip(*indices))
        einsum_string += f",{indices[::-1][0]}->{indices[::-1][1]}"
        einsum_args = [self.rotation_matrix] * rank + [tensor]

        return np.einsum(einsum_string, *einsum_args)

    def are_symmetrically_related(self, point_a: ArrayLike, point_b: ArrayLike, tol: float = 0.001) -> bool:
        """
        Checks if two points are symmetrically related.

        Args:
            point_a (3x1 array): First point.
            point_b (3x1 array): Second point.
            tol (float): Absolute tolerance for checking distance.

        Returns:
            True if self.operate(point_a) == point_b or vice versa.
        """
        if np.allclose(self.operate(point_a), point_b, atol=tol):
            return True
        if np.allclose(self.operate(point_b), point_a, atol=tol):
            return True
        return False

    def are_symmetrically_related_vectors(
        self,
        from_a: ArrayLike,
        to_a: ArrayLike,
        r_a: ArrayLike,
        from_b: ArrayLike,
        to_b: ArrayLike,
        r_b: ArrayLike,
        tol: float = 0.001,
    ) -> tuple[bool, bool]:
        """
        Checks if two vectors, or rather two vectors that connect two points
        each are symmetrically related. r_a and r_b give the change of unit
        cells. Two vectors are also considered symmetrically equivalent if starting
        and end point are exchanged.

        Args:
            from_a (3x1 array): Starting point of the first vector.
            to_a (3x1 array): Ending point of the first vector.
            from_b (3x1 array): Starting point of the second vector.
            to_b (3x1 array): Ending point of the second vector.
            r_a (3x1 array): Change of unit cell of the first vector.
            r_b (3x1 array): Change of unit cell of the second vector.
            tol (float): Absolute tolerance for checking distance.

        Returns:
            (are_related, is_reversed)
        """
        from_c = self.operate(from_a)
        to_c = self.operate(to_a)

        floored = np.floor([from_c, to_c])
        is_too_close = np.abs([from_c, to_c] - floored) > 1 - tol
        floored[is_too_close] += 1

        r_c = self.apply_rotation_only(r_a) - floored[0] + floored[1]
        from_c = from_c % 1
        to_c = to_c % 1

        if np.allclose(from_b, from_c, atol=tol) and np.allclose(to_b, to_c) and np.allclose(r_b, r_c, atol=tol):
            return (True, False)
        if np.allclose(to_b, from_c, atol=tol) and np.allclose(from_b, to_c) and np.allclose(r_b, -r_c, atol=tol):
            return (True, True)
        return (False, False)

    @property
    def rotation_matrix(self) -> np.ndarray:
        """
        A 3x3 numpy.array representing the rotation matrix.
        """
        return self.affine_matrix[0:3][:, 0:3]

    @property
    def translation_vector(self) -> np.ndarray:
        """
        A rank 1 numpy.array of dim 3 representing the translation vector.
        """
        return self.affine_matrix[0:3][:, 3]

    def __mul__(self, other):
        """
        Returns a new SymmOp which is equivalent to apply the "other" SymmOp
        followed by this one.
        """
        new_matrix = np.dot(self.affine_matrix, other.affine_matrix)
        return SymmOp(new_matrix)

    @property
    def inverse(self) -> SymmOp:
        """
        Returns inverse of transformation.
        """
        invr = np.linalg.inv(self.affine_matrix)
        return SymmOp(invr)

    @staticmethod
    def from_axis_angle_and_translation(
        axis: ArrayLike, angle: float, angle_in_radians: bool = False, translation_vec: ArrayLike = (0, 0, 0)
    ) -> SymmOp:
        """
        Generates a SymmOp for a rotation about a given axis plus translation.

        Args:
            axis: The axis of rotation in Cartesian space. For example,
                [1, 0, 0]indicates rotation about x-axis.
            angle (float): Angle of rotation.
            angle_in_radians (bool): Set to True if angles are given in
                radians. Or else, units of degrees are assumed.
            translation_vec: A translation vector. Defaults to zero.

        Returns:
            SymmOp for a rotation about given axis and translation.
        """
        if isinstance(axis, (tuple, list)):
            axis = np.array(axis)

        vec = np.array(translation_vec)

        a = angle if angle_in_radians else angle * pi / 180
        cosa = cos(a)
        sina = sin(a)
        u = axis / np.linalg.norm(axis)  # type: ignore
        r = np.zeros((3, 3))
        r[0, 0] = cosa + u[0] ** 2 * (1 - cosa)  # type: ignore
        r[0, 1] = u[0] * u[1] * (1 - cosa) - u[2] * sina  # type: ignore
        r[0, 2] = u[0] * u[2] * (1 - cosa) + u[1] * sina  # type: ignore
        r[1, 0] = u[0] * u[1] * (1 - cosa) + u[2] * sina  # type: ignore
        r[1, 1] = cosa + u[1] ** 2 * (1 - cosa)  # type: ignore
        r[1, 2] = u[1] * u[2] * (1 - cosa) - u[0] * sina  # type: ignore
        r[2, 0] = u[0] * u[2] * (1 - cosa) - u[1] * sina  # type: ignore
        r[2, 1] = u[1] * u[2] * (1 - cosa) + u[0] * sina  # type: ignore
        r[2, 2] = cosa + u[2] ** 2 * (1 - cosa)  # type: ignore

        return SymmOp.from_rotation_and_translation(r, vec)

    @staticmethod
    def from_origin_axis_angle(
        origin: ArrayLike, axis: ArrayLike, angle: float, angle_in_radians: bool = False
    ) -> SymmOp:
        """
        Generates a SymmOp for a rotation about a given axis through an
        origin.

        Args:
            origin (3x1 array): The origin which the axis passes through.
            axis (3x1 array): The axis of rotation in Cartesian space. For
                example, [1, 0, 0]indicates rotation about x-axis.
            angle (float): Angle of rotation.
            angle_in_radians (bool): Set to True if angles are given in
                radians. Or else, units of degrees are assumed.

        Returns:
            SymmOp.
        """
        theta = angle * pi / 180 if not angle_in_radians else angle
        a = origin[0]  # type: ignore
        b = origin[1]  # type: ignore
        c = origin[2]  # type: ignore
        u = axis[0]  # type: ignore
        v = axis[1]  # type: ignore
        w = axis[2]  # type: ignore
        # Set some intermediate values.
        u2 = u * u  # type: ignore
        v2 = v * v  # type: ignore
        w2 = w * w  # type: ignore
        cos_t = cos(theta)
        sin_t = sin(theta)
        l2 = u2 + v2 + w2  # type: ignore
        l = sqrt(l2)  # type: ignore

        # Build the matrix entries element by element.
        m11 = (u2 + (v2 + w2) * cos_t) / l2  # type: ignore
        m12 = (u * v * (1 - cos_t) - w * l * sin_t) / l2  # type: ignore
        m13 = (u * w * (1 - cos_t) + v * l * sin_t) / l2  # type: ignore
        m14 = (  # type: ignore
            a * (v2 + w2)  # type: ignore
            - u * (b * v + c * w)  # type: ignore
            + (u * (b * v + c * w) - a * (v2 + w2)) * cos_t  # type: ignore
            + (b * w - c * v) * l * sin_t  # type: ignore
        ) / l2

        m21 = (u * v * (1 - cos_t) + w * l * sin_t) / l2  # type: ignore
        m22 = (v2 + (u2 + w2) * cos_t) / l2  # type: ignore
        m23 = (v * w * (1 - cos_t) - u * l * sin_t) / l2  # type: ignore
        m24 = (  # type: ignore
            b * (u2 + w2)  # type: ignore
            - v * (a * u + c * w)  # type: ignore
            + (v * (a * u + c * w) - b * (u2 + w2)) * cos_t  # type: ignore
            + (c * u - a * w) * l * sin_t  # type: ignore
        ) / l2

        m31 = (u * w * (1 - cos_t) - v * l * sin_t) / l2  # type: ignore
        m32 = (v * w * (1 - cos_t) + u * l * sin_t) / l2  # type: ignore
        m33 = (w2 + (u2 + v2) * cos_t) / l2  # type: ignore
        m34 = (  # type: ignore
            c * (u2 + v2)  # type: ignore
            - w * (a * u + b * v)  # type: ignore
            + (w * (a * u + b * v) - c * (u2 + v2)) * cos_t  # type: ignore
            + (a * v - b * u) * l * sin_t  # type: ignore
        ) / l2

        return SymmOp(
            [  # type: ignore
                [m11, m12, m13, m14],
                [m21, m22, m23, m24],
                [m31, m32, m33, m34],
                [0, 0, 0, 1],
            ]
        )

    @staticmethod
    def reflection(normal: ArrayLike, origin: ArrayLike = (0, 0, 0)) -> SymmOp:
        """
        Returns reflection symmetry operation.

        Args:
            normal (3x1 array): Vector of the normal to the plane of
                reflection.
            origin (3x1 array): A point in which the mirror plane passes
                through.

        Returns:
            SymmOp for the reflection about the plane
        """
        # Normalize the normal vector first.
        n = np.array(normal, dtype=float) / np.linalg.norm(normal)

        u, v, w = n

        translation = np.eye(4)
        translation[0:3, 3] = -np.array(origin)

        xx = 1 - 2 * u**2
        yy = 1 - 2 * v**2
        zz = 1 - 2 * w**2
        xy = -2 * u * v
        xz = -2 * u * w
        yz = -2 * v * w
        mirror_mat = [[xx, xy, xz, 0], [xy, yy, yz, 0], [xz, yz, zz, 0], [0, 0, 0, 1]]

        if np.linalg.norm(origin) > 1e-6:
            mirror_mat = np.dot(np.linalg.inv(translation), np.dot(mirror_mat, translation))
        return SymmOp(mirror_mat)

    @staticmethod
    def inversion(origin: ArrayLike = (0, 0, 0)) -> SymmOp:
        """
        Inversion symmetry operation about axis.

        Args:
            origin (3x1 array): Origin of the inversion operation. Defaults
                to [0, 0, 0].

        Returns:
            SymmOp representing an inversion operation about the origin.
        """
        mat = -np.eye(4)
        mat[3, 3] = 1
        mat[0:3, 3] = 2 * np.array(origin)
        return SymmOp(mat)

    @staticmethod
    def rotoreflection(axis: ArrayLike, angle: float, origin: ArrayLike = (0, 0, 0)) -> SymmOp:
        """
        Returns a roto-reflection symmetry operation

        Args:
            axis (3x1 array): Axis of rotation / mirror normal
            angle (float): Angle in degrees
            origin (3x1 array): Point left invariant by roto-reflection.
                Defaults to (0, 0, 0).

        Return:
            Roto-reflection operation
        """
        rot = SymmOp.from_origin_axis_angle(origin, axis, angle)
        refl = SymmOp.reflection(axis, origin)
        m = np.dot(rot.affine_matrix, refl.affine_matrix)
        return SymmOp(m)

    def as_dict(self) -> dict[str, Any]:
        """
        :return: MSONAble dict.
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "matrix": self.affine_matrix.tolist(),
            "tolerance": self.tol,
        }

    def as_xyz_string(self) -> str:
        """
        Returns a string of the form 'x, y, z', '-x, -y, z',
        '-y+1/2, x+1/2, z+1/2', etc. Only works for integer rotation matrices
        """
        # test for invalid rotation matrix
        if not np.all(np.isclose(self.rotation_matrix, np.round(self.rotation_matrix))):
            warnings.warn("Rotation matrix should be integer")

        return transformation_to_string(self.rotation_matrix, translation_vec=self.translation_vector, delim=", ")

    @staticmethod
    def from_xyz_string(xyz_string: str) -> SymmOp:
        """
        Args:
            xyz_string: string of the form 'x, y, z', '-x, -y, z',
                '-2y+1/2, 3x+1/2, z-y+1/2', etc.

        Returns:
            SymmOp
        """
        rot_matrix = np.zeros((3, 3))
        trans = np.zeros(3)
        toks = xyz_string.strip().replace(" ", "").lower().split(",")
        re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
        re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
        for i, tok in enumerate(toks):
            # build the rotation matrix
            for m in re_rot.finditer(tok):
                factor = -1.0 if m.group(1) == "-" else 1.0
                if m.group(2) != "":
                    factor *= float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                j = ord(m.group(4)) - 120
                rot_matrix[i, j] = factor
            # build the translation vector
            for m in re_trans.finditer(tok):
                factor = -1 if m.group(1) == "-" else 1
                num = float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                trans[i] = num * factor
        return SymmOp.from_rotation_and_translation(rot_matrix, trans)

    @classmethod
    def from_dict(cls, d) -> SymmOp:
        """
        :param d: dict
        :return: SymmOp from dict representation.
        """
        return cls(d["matrix"], d["tolerance"])


class MagSymmOp(SymmOp):
    """
    Thin wrapper around SymmOp to extend it to support magnetic symmetry by including a time
    reversal operator. Magnetic symmetry is similar to conventional crystal symmetry, except
    symmetry is reduced by the addition of a time reversal operator which acts on an atom's magnetic
    moment.
    """

    def __init__(self, affine_transformation_matrix: ArrayLike, time_reversal: int, tol: float = 0.01):
        """
        Initializes the MagSymmOp from a 4x4 affine transformation matrix and time reversal
        operator. In general, this constructor should not be used unless you are transferring
        rotations. Use the static constructors instead to generate a SymmOp from proper rotations
        and translation.

        Args:
            affine_transformation_matrix (4x4 array): Representing an
                affine transformation.
            time_reversal (int): 1 or -1
            tol (float): Tolerance for determining if matrices are equal.
        """
        SymmOp.__init__(self, affine_transformation_matrix, tol=tol)
        if time_reversal not in (-1, 1):
            raise Exception(f"Time reversal operator not well defined: {time_reversal}, {type(time_reversal)}")
        self.time_reversal = time_reversal

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SymmOp):
            return NotImplemented
        return np.allclose(self.affine_matrix, other.affine_matrix, atol=self.tol) and (
            self.time_reversal == other.time_reversal
        )

    def __str__(self):
        return self.as_xyzt_string()

    def __repr__(self):
        output = [
            "Rot:",
            str(self.affine_matrix[0:3][:, 0:3]),
            "tau",
            str(self.affine_matrix[0:3][:, 3]),
            "Time reversal:",
            str(self.time_reversal),
        ]
        return "\n".join(output)

    def __hash__(self):
        # useful for obtaining a set of unique MagSymmOps
        hashable_value = (*tuple(self.affine_matrix.flatten()), self.time_reversal)
        return hash(hashable_value)

    def operate_magmom(self, magmom):
        """
        Apply time reversal operator on the magnetic moment. Note that
        magnetic moments transform as axial vectors, not polar vectors.

        See 'Symmetry and magnetic structures', Rodríguez-Carvajal and
        Bourée for a good discussion. DOI: 10.1051/epjconf/20122200010

        Args:
            magmom: Magnetic moment as electronic_structure.core.Magmom
            class or as list or np array-like

        Returns:
            Magnetic moment after operator applied as Magmom class
        """
        magmom = Magmom(magmom)  # type casting to handle lists as input

        transformed_moment = (
            self.apply_rotation_only(magmom.global_moment) * np.linalg.det(self.rotation_matrix) * self.time_reversal
        )

        # retains input spin axis if different from default
        return Magmom.from_global_moment_and_saxis(transformed_moment, magmom.saxis)

    @classmethod
    def from_symmop(cls, symmop: SymmOp, time_reversal) -> MagSymmOp:
        """
        Initialize a MagSymmOp from a SymmOp and time reversal operator.

        Args:
            symmop (SymmOp): SymmOp
            time_reversal (int): Time reversal operator, +1 or -1.

        Returns:
            MagSymmOp object
        """
        mag_symmop = cls(symmop.affine_matrix, time_reversal, symmop.tol)
        return mag_symmop

    @staticmethod
    def from_rotation_and_translation_and_time_reversal(
        rotation_matrix: ArrayLike = ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        translation_vec: ArrayLike = (0, 0, 0),
        time_reversal: int = 1,
        tol: float = 0.1,
    ) -> MagSymmOp:
        """
        Creates a symmetry operation from a rotation matrix, translation
        vector and time reversal operator.

        Args:
            rotation_matrix (3x3 array): Rotation matrix.
            translation_vec (3x1 array): Translation vector.
            time_reversal (int): Time reversal operator, +1 or -1.
            tol (float): Tolerance to determine if rotation matrix is valid.

        Returns:
            MagSymmOp object
        """
        symm_op = SymmOp.from_rotation_and_translation(
            rotation_matrix=rotation_matrix, translation_vec=translation_vec, tol=tol
        )
        return MagSymmOp.from_symmop(symm_op, time_reversal)

    @staticmethod
    def from_xyzt_string(xyzt_string: str) -> MagSymmOp:
        """
        Args:
            xyz_string: string of the form 'x, y, z, +1', '-x, -y, z, -1',
                '-2y+1/2, 3x+1/2, z-y+1/2, +1', etc.

        Returns:
            MagSymmOp object
        """
        symm_op = SymmOp.from_xyz_string(xyzt_string.rsplit(",", 1)[0])
        try:
            time_reversal = int(xyzt_string.rsplit(",", 1)[1])
        except Exception:
            raise Exception("Time reversal operator could not be parsed.")
        return MagSymmOp.from_symmop(symm_op, time_reversal)

    def as_xyzt_string(self) -> str:
        """
        Returns a string of the form 'x, y, z, +1', '-x, -y, z, -1',
        '-y+1/2, x+1/2, z+1/2, +1', etc. Only works for integer rotation matrices
        """
        xyzt_string = SymmOp.as_xyz_string(self)
        return f"{xyzt_string}, {self.time_reversal:+}"

    def as_dict(self) -> dict[str, Any]:
        """
        :return: MSONABle dict
        """
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "matrix": self.affine_matrix.tolist(),
            "tolerance": self.tol,
            "time_reversal": self.time_reversal,
        }

    @classmethod
    def from_dict(cls, d: dict) -> MagSymmOp:
        """
        :param d: dict
        :return: MagneticSymmOp from dict representation.
        """
        return cls(d["matrix"], tol=d["tolerance"], time_reversal=d["time_reversal"])
