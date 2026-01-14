"""This module provides core classes to define electronic structure,
including Spin, Orbital and Magmom.
"""

from __future__ import annotations

from enum import Enum, unique
from typing import TYPE_CHECKING, Literal, cast

import numpy as np
from monty.json import MSONable

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import ArrayLike, NDArray
    from typing_extensions import Self

    from pymatgen.core import Lattice
    from pymatgen.util.typing import MagMomentLike


@unique
class Spin(Enum):
    """Enum type for Spin. Only up and down. Usage: Spin.up, Spin.down."""

    up, down = 1, -1

    def __int__(self) -> Literal[-1, 1]:
        return cast("Literal[-1, 1]", self.value)

    def __float__(self) -> float:
        return float(self.value)

    def __str__(self) -> Literal["-1", "1"]:
        return cast("Literal['-1', '1']", str(self.value))


@unique
class OrbitalType(Enum):
    """Enum type for orbital type. Indices are the azimuthal quantum number l."""

    s = 0
    p = 1
    d = 2
    f = 3

    def __str__(self) -> Literal["s", "p", "d", "f"]:
        return cast("Literal['s', 'p', 'd', 'f']", str(self.name))


@unique
class Orbital(Enum):
    """Enum type for specific orbitals. The indices are the order in
    which the orbitals are reported in VASP and has no special meaning.
    """

    s = 0
    py = 1
    pz = 2
    px = 3
    dxy = 4
    dyz = 5
    dz2 = 6
    dxz = 7
    dx2 = 8
    f_3 = 9
    f_2 = 10
    f_1 = 11
    f0 = 12
    f1 = 13
    f2 = 14
    f3 = 15

    def __int__(self) -> int:
        return self.value

    def __str__(self) -> str:
        return str(self.name)

    @property
    def orbital_type(self) -> OrbitalType:
        """OrbitalType of an orbital."""
        return OrbitalType[self.name[0]]


class Magmom(MSONable):
    """In active development. Use with caution, feedback is appreciated.

    Class to handle magnetic moments. Define the magnetic moment of a
    site or species relative to a spin quantization axis. Designed for
    use in electronic structure calculations.

        * For the general case, Magmom can be specified by a 3D vector,
        e.g. m = Magmom([1.0, 1.0, 2.0]), and indexing will work as
        expected, e.g. m[0] gives 1.0.

        * For collinear calculations, Magmom can assumed to be float-like,
        e.g. m = Magmom(5.0) will work as expected, e.g. float(m) gives 5.0.

    Both cases should be safe and shouldn't give any surprise,
    and more advanced functionality is available if required.

    There are also static methods for sequences of magmoms:

        * Magmom.are_collinear(magmoms) - If True, a collinear electronic
        structure calculation can be safely initialized, with float(Magmom)
        giving the expected scalar magnetic moment value.

        * Magmom.get_consistent_set_and_saxis(magmoms) - For non-collinear
        electronic structure calculations, a global and consistent spin axis
        has to be used. This method returns a list of Magmoms which all
        share a common spin axis, along with the global spin axis.

    All methods that take sequence of magmoms will accept either Magmom
    objects, or as scalars/lists and will automatically convert to Magmom
    representations internally.

    The following methods are also useful for VASP calculations:
        - Magmom.get_xyz_magmom_with_001_saxis()
        - Magmom.get_00t_magmom_with_xyz_saxis()

        See VASP documentation for more information:
            https://cms.mpi.univie.ac.at/wiki/index.php/SAXIS
    """

    def __init__(
        self,
        moment: MagMomentLike,
        saxis: ArrayLike = (0, 0, 1),
    ) -> None:
        """
        Args:
            moment (float | Sequence[float] | NDArray, Magmom): Magnetic moment.
            saxis (tuple[float, float, float]): Spin axis, and will be converted to unit
                vector (default is (0, 0, 1)).
        """
        # Init from another Magmom instance
        if isinstance(moment, type(self)):
            saxis = moment.saxis  # type:ignore[has-type]
            moment = moment.moment  # type:ignore[has-type]

        magmom = np.array(moment, dtype="d")
        if magmom.ndim == 0:
            magmom = magmom * (0, 0, 1)  # (ruff-preview) noqa: PLR6104

        self.moment = magmom

        saxis = np.array(saxis, dtype="d")

        self.saxis = saxis / np.linalg.norm(saxis)

    def __getitem__(self, key):
        return self.moment[key]

    def __iter__(self):
        return iter(self.moment)

    def __abs__(self) -> float:
        return np.linalg.norm(self.moment)  # type: ignore[return-value]

    def __eq__(self, other: object) -> bool:
        """Whether global magnetic moments are the same, saxis can differ."""
        try:
            other_magmom = type(self)(other)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return NotImplemented

        return np.allclose(self.global_moment, other_magmom.global_moment)

    def __lt__(self, other: Self) -> bool:
        return abs(self) < abs(other)

    def __neg__(self) -> Self:
        return type(self)(-self.moment, saxis=self.saxis)

    def __hash__(self) -> int:
        return hash(tuple(self.moment) + tuple(self.saxis))

    def __float__(self) -> float:
        """Get magnitude of magnetic moment with a sign with respect to
        an arbitrary direction.

        Should give unsurprising output if Magmom is treated like a
        float or if a set of Magmoms describes a collinear structure.

        Implemented this way rather than simpler abs(self) so that
        moments will have a consistent sign in case of e.g.
        antiferromagnetic collinear structures without additional
        user intervention.

        However, should be used with caution for non-collinear
        structures and might give nonsensical results except in the case
        of only slightly non-collinear structures (e.g. small canting).

        This method is also used to obtain "diff" VolumetricDensity
        in pymatgen.io.vasp.outputs.VolumetricDensity when processing
        CHGCARs from SOC calculations.
        """
        return float(self.get_00t_magmom_with_xyz_saxis()[2])

    def __str__(self) -> str:
        return str(float(self))

    def __repr__(self) -> str:
        if np.allclose(self.saxis, (0, 0, 1)):
            return f"Magnetic moment {self.moment}"
        return f"Magnetic moment {self.moment} (spin axis = {self.saxis})"

    @classmethod
    def from_global_moment_and_saxis(
        cls,
        global_moment: MagMomentLike,
        saxis: tuple[float, float, float],
    ) -> Self:
        """Initialize Magmom from a given global magnetic moment,
        i.e. magnetic moment with saxis=(0, 0, 1), and provided saxis.

        Method is useful if you do not know the components of your
        magnetic moment in frame of your desired spin axis.

        Args:
            global_moment (MagMomentLike): Global magnetic moment.
            saxis (tuple[float, float, float]): Spin axis.
        """
        magmom = cls(global_moment)
        return cls(magmom.get_moment(saxis=saxis), saxis=saxis)

    def get_moment(self, saxis: tuple[float, float, float] = (0, 0, 1)) -> NDArray:
        """Get magnetic moment relative to a given spin quantization axis.
        If no axis is provided, moment will be given relative to the
        Magmom's internal spin quantization axis, i.e. equivalent to
        Magmom.moment.

        Args:
            saxis (tuple[float, float, float]): Spin quantization axis.

        Returns:
            NDArray of length 3.
        """

        def get_transformation_matrix(
            saxis: tuple[float, float, float],
        ) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
            """Get the matrix to transform spin axis to z-axis."""
            saxis /= np.linalg.norm(saxis)

            alpha = np.arctan2(saxis[1], saxis[0])
            beta = np.arctan2(np.sqrt(saxis[0] ** 2 + saxis[1] ** 2), saxis[2])

            cos_a = np.cos(alpha)
            cos_b = np.cos(beta)
            sin_a = np.sin(alpha)
            sin_b = np.sin(beta)

            return (
                (cos_b * cos_a, -sin_a, sin_b * cos_a),
                (cos_b * sin_a, cos_a, sin_b * sin_a),
                (-sin_b, 0, cos_b),
            )

        def get_transformation_matrix_inv(
            saxis: tuple[float, float, float],
        ) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
            """Get the inverse of matrix to transform spin axis to z-axis."""
            saxis /= np.linalg.norm(saxis)

            alpha = np.arctan2(saxis[1], saxis[0])
            beta = np.arctan2(np.sqrt(saxis[0] ** 2 + saxis[1] ** 2), saxis[2])

            cos_a = np.cos(alpha)
            cos_b = np.cos(beta)
            sin_a = np.sin(alpha)
            sin_b = np.sin(beta)

            return (
                (cos_b * cos_a, cos_b * sin_a, -sin_b),
                (-sin_a, cos_a, 0),
                (sin_b * cos_a, sin_b * sin_a, cos_b),
            )

        # Transform to moment with spin axis (0, 0, 1)
        trafo_mat_inv = get_transformation_matrix_inv(self.saxis)
        moment = np.matmul(self.moment, trafo_mat_inv)

        # Transform to new saxis
        trafo_mat = get_transformation_matrix(saxis)
        moment = np.matmul(moment, trafo_mat)

        # Round small values to zero
        moment[np.abs(moment) < 1e-8] = 0

        return moment

    @property
    def global_moment(self) -> NDArray:
        """The magnetic moment defined in an arbitrary global reference frame,
        as a np.array of length 3.
        """
        return self.get_moment()

    @property
    def projection(self) -> float:
        """Project moment along spin quantization axis.

        Useful for obtaining collinear approximation for slightly non-collinear magmoms.

        Returns:
            float: The projected moment.
        """
        return np.dot(self.moment, self.saxis)

    def get_xyz_magmom_with_001_saxis(self) -> Self:
        """Get a Magmom in the default setting of saxis = (0, 0, 1) and
        the magnetic moment rotated as required.

        Returns:
            Magmom
        """
        return type(self)(self.get_moment())

    def get_00t_magmom_with_xyz_saxis(self) -> Self:
        """For internal implementation reasons, the non-collinear calculations
        in VASP prefer the following.

            MAGMOM = 0 0 total_magnetic_moment
            SAXIS = x y z

        to an equivalent:

            MAGMOM = x y z
            SAXIS = 0 0 1

        Returns:
            Magmom: With magnetic moment (0, 0, t), where t is the total magnetic
                moment, and saxis rotated as required.

                A consistent direction of saxis is applied such that t might be
                positive or negative depending on the direction of the initial moment.
                This is useful in the case of collinear structures, rather than
                assuming t is always positive.
        """
        # Reference direction gives sign of moment arbitrarily,
        # there can be a pathological case where a consistent sign
        # is not possible if the magnetic moments are aligned along the
        # reference direction, but in practice this is unlikely to happen.
        ref_direction = np.array([1.01, 1.02, 1.03])
        total_magmom = abs(self)
        if total_magmom != 0:
            new_saxis = self.moment / np.linalg.norm(self.moment)
            if np.dot(ref_direction, new_saxis) < 0:
                total_magmom = -total_magmom
                new_saxis = -new_saxis
            return type(self)([0, 0, total_magmom], saxis=new_saxis)
        return type(self)(self)

    @staticmethod
    def have_consistent_saxis(magmoms: Sequence[MagMomentLike]) -> bool:
        """Check whether all Magmoms have a consistent spin quantization axis.

        To write MAGMOM tags to a VASP INCAR, a consistent global SAXIS value for
        all magmoms has to be used.

        If spin axes are inconsistent, can create a consistent set with:
            Magmom.get_consistent_set(magmoms).

        Args:
            magmoms (Sequence[MagMomentLike]): Magmoms.

        Returns:
            bool
        """
        _magmoms: list[Magmom] = [Magmom(magmom) for magmom in magmoms]
        ref_saxis = _magmoms[0].saxis
        match_ref = [magmom.saxis == ref_saxis for magmom in _magmoms]
        return bool(np.all(match_ref))

    @staticmethod
    def get_consistent_set_and_saxis(
        magmoms: Sequence[MagMomentLike],
        saxis: tuple[float, float, float] | None = None,
    ) -> tuple[list[NDArray], NDArray]:
        """Ensure magmoms use the same spin axis.

        Args:
            magmoms (Sequence[MagMomentLike]): Magmoms, floats or vectors.
            saxis (tuple[float, float, float]): An optional global spin axis.

        Returns:
            tuple[list[Magmom], NDArray]: Magmoms and their global spin axes.
        """
        _magmoms: list[Magmom] = [Magmom(magmom) for magmom in magmoms]
        _saxis: NDArray = Magmom.get_suggested_saxis(_magmoms) if saxis is None else saxis / np.linalg.norm(saxis)
        moments: list[NDArray] = [magmom.get_moment(saxis=_saxis) for magmom in _magmoms]  # type: ignore[arg-type]
        return moments, _saxis

    @staticmethod
    def get_suggested_saxis(magmoms: Sequence[MagMomentLike]) -> NDArray:
        """Get a suggested spin axis for magmoms, taking the largest magnetic
        moment as the reference. For calculations with collinear spins,
        this would give a sensible saxis for a NCL calculation.

        Args:
            magmoms (Sequence[MagMomentLike]): Magmoms, floats or vectors.

        Returns:
            NDArray of length 3
        """
        # Heuristic, will pick largest magmom as the reference.
        # Useful for creating collinear approximations of
        # e.g. slightly canted magnetic structures.
        # For fully collinear structures, will return expected result.

        _magmoms: list[Magmom] = [Magmom(magmom) for magmom in magmoms]
        # Keep non-zero magmoms only
        _magmoms = [magmom for magmom in _magmoms if abs(magmom)]  # type: ignore[arg-type]
        _magmoms.sort(reverse=True)

        if _magmoms:
            return _magmoms[0].get_00t_magmom_with_xyz_saxis().saxis
        return np.array([0, 0, 1], dtype="d")

    @staticmethod
    def are_collinear(magmoms: Sequence[MagMomentLike]) -> bool:
        """Check if a list of magnetic moments are collinear with each other.

        Args:
            magmoms (Sequence[MagMomentLike]): Magmoms, floats or vectors.

        Returns:
            bool.
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        if not Magmom.have_consistent_saxis(magmoms):
            magmoms = Magmom.get_consistent_set_and_saxis(magmoms)[0]

        # Convert to numpy array for convenience
        magmoms = np.array([list(cast("Magmom", magmom)) for magmom in magmoms])  # type: ignore[assignment]
        magmoms = magmoms[np.any(magmoms, axis=1)]  # type:ignore[assignment,arg-type]
        if len(magmoms) == 0:
            return True

        # Use first moment as reference to compare against
        ref_magmom = magmoms[0]
        # Magnitude of cross products != 0 if non-collinear with reference
        num_ncl = np.count_nonzero(np.linalg.norm(np.cross(ref_magmom, magmoms), axis=1))  # type: ignore[arg-type]
        return num_ncl == 0

    @classmethod
    def from_moment_relative_to_crystal_axes(
        cls,
        moment: tuple[float, float, float],
        lattice: Lattice,
    ) -> Self:
        """Obtain a Magmom object from a magnetic moment provided
        relative to crystal axes.

        Used for obtaining moments from magCIF file.

        Args:
            moment (tuple[float, float, float]): Magnetic moment.
            lattice (Lattice): Lattice.

        Returns:
            Magmom
        """
        # Get matrix representing unit lattice vectors
        unit_m = lattice.matrix / np.linalg.norm(lattice.matrix, axis=1)[:, None]
        _moment: NDArray = np.matmul(list(moment), unit_m)
        # Round small values to zero
        _moment[np.abs(_moment) < 1e-8] = 0
        return cls(_moment)

    def get_moment_relative_to_crystal_axes(self, lattice: Lattice) -> tuple[float, float, float]:
        """If scalar magmoms, moments will be given arbitrarily along z.

        Used for writing moments to magCIF file.

        Args:
            lattice (Lattice): The lattice.

        Returns:
            tuple[float, float, float]
        """
        # Get matrix representing unit lattice vectors
        unit_m = lattice.matrix / np.linalg.norm(lattice.matrix, axis=1)[:, None]
        moment = np.matmul(self.global_moment, np.linalg.inv(unit_m))
        # Round small values to zero
        moment[np.abs(moment) < 1e-8] = 0
        return moment
