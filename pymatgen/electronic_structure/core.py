"""
This module provides core classes needed by all define electronic structure,
such as the Spin, Orbital, etc.
"""

from __future__ import annotations

from enum import Enum, unique

import numpy as np
from monty.json import MSONable


@unique
class Spin(Enum):
    """
    Enum type for Spin. Only up and down.
    Usage: Spin.up, Spin.down.
    """

    up, down = 1, -1

    def __int__(self):
        return self.value

    def __float__(self):
        return float(self.value)

    def __str__(self):
        return str(self.value)


@unique
class OrbitalType(Enum):
    """
    Enum type for orbital type. Indices are basically the azimuthal quantum
    number, l.
    """

    s = 0
    p = 1
    d = 2
    f = 3

    def __str__(self):
        return str(self.name)


@unique
class Orbital(Enum):
    """
    Enum type for specific orbitals. The indices are basically the order in
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

    def __int__(self):
        return self.value

    def __str__(self):
        return str(self.name)

    @property
    def orbital_type(self):
        """Returns OrbitalType of an orbital."""
        # pylint: disable=E1136
        return OrbitalType[self.name[0]]


class Magmom(MSONable):
    """
    New class in active development. Use with caution, feedback is
    appreciated.

    Class to handle magnetic moments. Defines the magnetic moment of a
    site or species relative to a spin quantization axis. Designed for
    use in electronic structure calculations.

    * For the general case, Magmom can be specified by a vector,
      e.g. m = Magmom([1.0, 1.0, 2.0]), and subscripts will work as
      expected, e.g. m[0] gives 1.0

    * For collinear calculations, Magmom can assumed to be scalar-like,
      e.g. m = Magmom(5.0) will work as expected, e.g. float(m) gives 5.0

    Both of these cases should be safe and shouldn't give any surprises,
    but more advanced functionality is available if required.

    There also exist useful static methods for lists of magmoms:

    * Magmom.are_collinear(magmoms) - if true, a collinear electronic
      structure calculation can be safely initialized, with float(Magmom)
      giving the expected scalar magnetic moment value

    * Magmom.get_consistent_set_and_saxis(magmoms) - for non-collinear
      electronic structure calculations, a global, consistent spin axis
      has to be used. This method returns a list of Magmoms which all
      share a common spin axis, along with the global spin axis.

    All methods that take lists of magmoms will accept magmoms either as
    Magmom objects or as scalars/lists and will automatically convert to
    a Magmom representation internally.

    The following methods are also particularly useful in the context of
    VASP calculations:

    * Magmom.get_xyz_magmom_with_001_saxis()
    * Magmom.get_00t_magmom_with_xyz_saxis()

    See VASP documentation for more information:

    https://cms.mpi.univie.ac.at/wiki/index.php/SAXIS
    """

    def __init__(self, moment, saxis=(0, 0, 1)):
        """
        :param moment: magnetic moment, supplied as float or list/np.ndarray
        :param saxis: spin axis, supplied as list/np.ndarray, parameter will
            be converted to unit vector (default is [0, 0, 1])
        :return: Magmom object
        """
        # to init from another Magmom instance
        if isinstance(moment, Magmom):
            saxis = moment.saxis
            moment = moment.moment

        moment = np.array(moment, dtype="d")
        if moment.ndim == 0:
            moment = moment * [0, 0, 1]

        self.moment = moment

        saxis = np.array(saxis, dtype="d")

        self.saxis = saxis / np.linalg.norm(saxis)

    @classmethod
    def from_global_moment_and_saxis(cls, global_moment, saxis):
        """
        Convenience method to initialize Magmom from a given global
        magnetic moment, i.e. magnetic moment with saxis=(0,0,1), and
        provided saxis.

        Method is useful if you do not know the components of your
        magnetic moment in frame of your desired saxis.

        :param global_moment:
        :param saxis: desired saxis
        :return:
        """
        magmom = Magmom(global_moment)
        return cls(magmom.get_moment(saxis=saxis), saxis=saxis)

    @classmethod
    def _get_transformation_matrix(cls, saxis):
        saxis = saxis / np.linalg.norm(saxis)

        alpha = np.arctan2(saxis[1], saxis[0])
        beta = np.arctan2(np.sqrt(saxis[0] ** 2 + saxis[1] ** 2), saxis[2])

        cos_a = np.cos(alpha)
        cos_b = np.cos(beta)
        sin_a = np.sin(alpha)
        sin_b = np.sin(beta)

        return [
            [cos_b * cos_a, -sin_a, sin_b * cos_a],
            [cos_b * sin_a, cos_a, sin_b * sin_a],
            [-sin_b, 0, cos_b],
        ]

    @classmethod
    def _get_transformation_matrix_inv(cls, saxis):
        saxis = saxis / np.linalg.norm(saxis)

        alpha = np.arctan2(saxis[1], saxis[0])
        beta = np.arctan2(np.sqrt(saxis[0] ** 2 + saxis[1] ** 2), saxis[2])

        cos_a = np.cos(alpha)
        cos_b = np.cos(beta)
        sin_a = np.sin(alpha)
        sin_b = np.sin(beta)

        return [
            [cos_b * cos_a, cos_b * sin_a, -sin_b],
            [-sin_a, cos_a, 0],
            [sin_b * cos_a, sin_b * sin_a, cos_b],
        ]

    def get_moment(self, saxis=(0, 0, 1)):
        """
        Get magnetic moment relative to a given spin quantization axis.
        If no axis is provided, moment will be given relative to the
        Magmom's internal spin quantization axis, i.e. equivalent to
        Magmom.moment.

        :param saxis: (list/numpy array) spin quantization axis
        :return: np.ndarray of length 3
        """
        # transform back to moment with spin axis [0, 0, 1]
        m_inv = self._get_transformation_matrix_inv(self.saxis)
        moment = np.matmul(self.moment, m_inv)

        # transform to new saxis
        m = self._get_transformation_matrix(saxis)
        moment = np.matmul(moment, m)

        # round small values to zero
        moment[np.abs(moment) < 1e-8] = 0

        return moment

    @property
    def global_moment(self):
        """
        Get the magnetic moment defined in an arbitrary global reference frame.

        :return: np.ndarray of length 3
        """
        return self.get_moment()

    @property
    def projection(self):
        """
        Projects moment along spin quantisation axis. Useful for obtaining
        collinear approximation for slightly non-collinear magmoms.

        :return: float
        """
        return np.dot(self.moment, self.saxis)

    def get_xyz_magmom_with_001_saxis(self):
        """
        Returns a Magmom in the default setting of saxis = [0, 0, 1] and
        the magnetic moment rotated as required.

        :return: Magmom
        """
        return Magmom(self.get_moment())

    def get_00t_magmom_with_xyz_saxis(self):
        """
        For internal implementation reasons, in non-collinear calculations VASP prefers the following.

            MAGMOM = 0 0 total_magnetic_moment
            SAXIS = x y z

        to an equivalent:

            MAGMOM = x y z
            SAXIS = 0 0 1

        This method returns a Magmom object with magnetic moment [0, 0, t],
        where t is the total magnetic moment, and saxis rotated as required.

        A consistent direction of saxis is applied such that t might be positive
        or negative depending on the direction of the initial moment. This is useful
        in the case of collinear structures, rather than constraining assuming
        t is always positive.

        :return: Magmom
        """
        # reference direction gives sign of moment
        # entirely arbitrary, there will always be a pathological case
        # where a consistent sign is not possible if the magnetic moments
        # are aligned along the reference direction, but in practice this
        # is unlikely to happen
        ref_direction = np.array([1.01, 1.02, 1.03])
        t = abs(self)
        if t != 0:
            new_saxis = self.moment / np.linalg.norm(self.moment)
            if np.dot(ref_direction, new_saxis) < 0:
                t = -t
                new_saxis = -new_saxis
            return Magmom([0, 0, t], saxis=new_saxis)
        return Magmom(self)

    @staticmethod
    def have_consistent_saxis(magmoms):
        """
        This method checks that all Magmom objects in a list have a
        consistent spin quantization axis. To write MAGMOM tags to a
        VASP INCAR, a global SAXIS value for all magmoms has to be used.
        If saxis are inconsistent, can create consistent set with:
        Magmom.get_consistent_set(magmoms).

        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :return: bool
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        ref_saxis = magmoms[0].saxis
        match_ref = [magmom.saxis == ref_saxis for magmom in magmoms]
        if np.all(match_ref):
            return True
        return False

    @staticmethod
    def get_consistent_set_and_saxis(magmoms, saxis=None):
        """
        Method to ensure a list of magmoms use the same spin axis.
        Returns a tuple of a list of Magmoms and their global spin axis.

        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :param saxis: can provide a specific global spin axis
        :return: (list of Magmoms, global spin axis) tuple
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        saxis = Magmom.get_suggested_saxis(magmoms) if saxis is None else saxis / np.linalg.norm(saxis)
        magmoms = [magmom.get_moment(saxis=saxis) for magmom in magmoms]
        return magmoms, saxis

    @staticmethod
    def get_suggested_saxis(magmoms):
        """
        This method returns a suggested spin axis for a set of magmoms,
        taking the largest magnetic moment as the reference. For calculations
        with collinear spins, this would give a sensible saxis for a ncl
        calculation.

        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :return: np.ndarray of length 3
        """
        # heuristic, will pick largest magmom as reference
        # useful for creating collinear approximations of
        # e.g. slightly canted magnetic structures
        # for fully collinear structures, will return expected
        # result

        magmoms = [Magmom(magmom) for magmom in magmoms]
        # filter only non-zero magmoms
        magmoms = [magmom for magmom in magmoms if abs(magmom)]
        magmoms.sort(reverse=True)
        if len(magmoms) > 0:
            return magmoms[0].get_00t_magmom_with_xyz_saxis().saxis
        return np.array([0, 0, 1], dtype="d")

    @staticmethod
    def are_collinear(magmoms) -> bool:
        """
        Method checks to see if a set of magnetic moments are collinear
        with each other.
        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :return: bool.
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        if not Magmom.have_consistent_saxis(magmoms):
            magmoms = Magmom.get_consistent_set_and_saxis(magmoms)[0]

        # convert to numpy array for convenience
        magmoms = np.array([list(magmom) for magmom in magmoms])
        magmoms = magmoms[np.any(magmoms, axis=1)]  # remove zero magmoms
        if len(magmoms) == 0:
            return True

        # use first moment as reference to compare against
        ref_magmom = magmoms[0]
        # magnitude of cross products != 0 if non-collinear with reference
        num_ncl = np.count_nonzero(np.linalg.norm(np.cross(ref_magmom, magmoms), axis=1))
        if num_ncl > 0:
            return False
        return True

    @classmethod
    def from_moment_relative_to_crystal_axes(cls, moment, lattice):
        """
        Obtaining a Magmom object from a magnetic moment provided
        relative to crystal axes.

        Used for obtaining moments from magCIF file.
        :param moment: list of floats specifying vector magmom
        :param lattice: Lattice
        :return: Magmom
        """
        # get matrix representing unit lattice vectors
        unit_m = lattice.matrix / np.linalg.norm(lattice.matrix, axis=1)[:, None]
        moment = np.matmul(list(moment), unit_m)
        # round small values to zero
        moment[np.abs(moment) < 1e-8] = 0
        return cls(moment)

    def get_moment_relative_to_crystal_axes(self, lattice):
        """
        If scalar magmoms, moments will be given arbitrarily along z.
        Used for writing moments to magCIF file.

        :param lattice: Lattice
        :return: vector as list of floats
        """
        # get matrix representing unit lattice vectors
        unit_m = lattice.matrix / np.linalg.norm(lattice.matrix, axis=1)[:, None]
        # note np.matmul() requires numpy version >= 1.10
        moment = np.matmul(self.global_moment, np.linalg.inv(unit_m))
        # round small values to zero
        moment[np.abs(moment) < 1e-8] = 0
        return moment

    def __getitem__(self, key):
        return self.moment[key]

    def __iter__(self):
        return iter(self.moment)

    def __abs__(self):
        return np.linalg.norm(self.moment)

    def __eq__(self, other: object) -> bool:
        """Equal if 'global' magnetic moments are the same, saxis can differ."""
        try:
            other_magmom = Magmom(other)
        except (TypeError, ValueError):
            return NotImplemented

        return np.allclose(self.global_moment, other_magmom.global_moment)

    def __lt__(self, other):
        return abs(self) < abs(other)

    def __neg__(self):
        return Magmom(-self.moment, saxis=self.saxis)

    def __hash__(self) -> int:
        return (tuple(self.moment) + tuple(self.saxis)).__hash__()

    def __float__(self):
        """
        Returns magnitude of magnetic moment with a sign with respect to
        an arbitrary direction.

        Should give unsurprising output if Magmom is treated like a
        scalar or if a set of Magmoms describes a collinear structure.

        Implemented this way rather than simpler abs(self) so that
        moments will have a consistent sign in case of e.g.
        antiferromagnetic collinear structures without additional
        user intervention.

        However, should be used with caution for non-collinear
        structures and might give nonsensical results except in the case
        of only slightly non-collinear structures (e.g. small canting).

        This approach is also used to obtain "diff" VolumetricDensity
        in pymatgen.io.vasp.outputs.VolumetricDensity when processing
        Chgcars from SOC calculations.
        """
        return float(self.get_00t_magmom_with_xyz_saxis()[2])

    def __str__(self):
        return str(float(self))

    def __repr__(self):
        if np.allclose(self.saxis, (0, 0, 1)):
            return f"Magnetic moment {self.moment}"
        return f"Magnetic moment {self.moment} (spin axis = {self.saxis})"
