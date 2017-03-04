# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
from monty.json import MSONable

"""
This module provides core classes needed by all define electronic structure,
such as the Spin, Orbital, etc.
"""


__author__ = "Shyue Ping Ong, Matthew Horton"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "2.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__status__ = "Production"
__date__ = "Sep 23, 2011"

from enum import Enum, unique
import numpy as np

@unique
class Spin(Enum):
    """
    Enum type for Spin.  Only up and down.
    Usage: Spin.up, Spin.down.
    """
    up, down = (1, -1)

    def __int__(self):
        return self.value

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
        return self.name


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
        return self.name

    @property
    def orbital_type(self):
        """
        Returns OrbitalType of an orbital.
        """
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

    def __init__(self, moment, saxis=None):
        """
        :param moment: magnetic moment, supplied as float or list/np.ndarray
        :param saxis: spin axis, supplied as list/np.ndarray, parameter will
        be converted to unit vector (default is [0, 0, 1])
        :return: Magmom object
        """

        # defaults with spin axis along z axis
        # TODO: VASP default has saxis[0] = 0+, could use np.nextafter(0,1) to be rigorous?
        self._default_saxis = np.array([0, 0, 1], dtype='d')
        self._default_saxis.setflags(write=False)

        # to init from another Magmom instance
        # TODO: is there a more Pythonic way of doing this?
        if isinstance(moment, type(self)) and saxis==None:
            saxis = moment.saxis
            moment = moment.moment

        if isinstance(moment, int) or isinstance(moment, float):
            self._moment = np.array([0, 0, moment], dtype='d')
        elif isinstance(moment, list) or isinstance(moment, np.ndarray):
            if len(moment) == 3:
                self._moment = np.array(moment, dtype='d')
            else:
                raise ValueError("Moment should be vector of length 3.")
        else:
            raise ValueError("Moment should be scalar or vector.")

        if saxis is None:
            self._saxis = self._default_saxis
        else:
            if len(saxis) != 3:
                raise ValueError("Spin axis should be vector of length 3.")
            else:
                self._saxis = np.array(saxis, dtype='d')
                # make unit vector
                self._saxis /= np.linalg.norm(self._saxis)

        # ideally, should treat Magmom as immutable
        self._moment.setflags(write=False)
        self._saxis.setflags(write=False)

    @property
    def saxis(self):
        """
        Get the magnetic moment's spin quantization axis.

        :return: np.ndarray of length 3
        """
        return self._saxis

    @property
    def moment(self):
        """
        Get the magnetic moment relative to its spin quantization axis.

        :return: np.ndarray of length 3
        """
        return self._moment

    def _get_transformation_matrix(self, saxis):

        saxis = saxis / np.linalg.norm(saxis)

        alpha = np.arctan2(saxis[1], saxis[0])
        beta = np.arctan2(np.sqrt(saxis[0]**2 + saxis[1]**2), saxis[2])

        cos_a = np.cos(alpha)
        cos_b = np.cos(beta)
        sin_a = np.sin(alpha)
        sin_b = np.sin(beta)

        m = [[cos_b*cos_a,   -sin_a,    sin_b*cos_a],
             [cos_b*sin_a,    cos_a,    sin_b*sin_a],
             [-sin_b,             0,          cos_b]]

        return m

    def _get_transformation_matrix_inv(self, saxis):

        saxis = saxis / np.linalg.norm(saxis)

        alpha = np.arctan2(saxis[1], saxis[0])
        beta = np.arctan2(np.sqrt(saxis[0]**2 + saxis[1]**2), saxis[2])

        cos_a = np.cos(alpha)
        cos_b = np.cos(beta)
        sin_a = np.sin(alpha)
        sin_b = np.sin(beta)

        m = [[cos_b*cos_a,  cos_b*sin_a, -sin_b],
             [-sin_a,             cos_a,      0],
             [sin_b*cos_a,  sin_b*sin_a,  cos_b]]

        return m


    def get_moment(self, saxis=None):
        """
        Get magnetic moment relative to a given spin quantization axis.
        If no axis is provided, moment will be given relative to the
        Magmom's internal spin quantization axis, i.e. equivalent to
        Magmom.moment

        :param axis: (list/numpy array) spin quantization axis
        :return: np.ndarray of length 3
        """

        if saxis is None:
            return self._moment
        elif isinstance(saxis, list) or isinstance(saxis, np.ndarray):
            if len(saxis) != 3:
                raise ValueError("Spin axis should be vector of length 3.")
        else:
            raise ValueError("Ill-defined spin axis.")

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
        # TODO: change _default_saxis to [0,0,1] everywhere
        return self.get_moment(saxis=self._default_saxis)

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
        return Magmom(self.get_moment(saxis=self._default_saxis), saxis=self._default_saxis)


    def get_00t_magmom_with_xyz_saxis(self):
        """
        For internal implementation reasons, in non-collinear calculations
        VASP prefers:

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
        ref_direction = np.array([0, 0, 1])
        t = abs(self)
        if t != 0:
            new_saxis = self.moment/np.linalg.norm(self.moment)
            if np.dot(ref_direction, new_saxis) < 0:
                t = -t
                new_saxis = -new_saxis
            return Magmom([0, 0, t], saxis=new_saxis)
        else:
            return Magmom(self)

    @staticmethod
    def have_consistent_saxis(magmoms):
        """
        This method checks that all Magmom objects in a list have a
        consistent spin quantization axis. To write MAGMOM tags to a
        VASP INCAR, a global SAXIS value for all magmoms has to be used.
        If saxis are inconsistent, can create consistent set with:
        Magmom.get_consistent_set(magmoms)

        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :return: bool
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        ref_saxis = magmoms[0].saxis
        match_ref = [magmom.saxis == ref_saxis for magmom in magmoms]
        if np.all(match_ref):
            return True
        else:
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
        if saxis is None:
            saxis = Magmom.get_suggested_saxis(magmoms)
        else:
            saxis = saxis/np.linalg.norm(saxis)
        magmoms = [magmom.get_moment(saxis=saxis) for magmom in magmoms]
        return (magmoms, saxis)


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
        else:
            return self._default_saxis

    @staticmethod
    def are_collinear(magmoms):
        """
        Method checks to see if a set of magnetic moments are collinear
        with each other.
        :param magmoms: list of magmoms (Magmoms, scalars or vectors)
        :return: bool
        """
        magmoms = [Magmom(magmom) for magmom in magmoms]
        if not Magmom.have_consistent_saxis(magmoms):
            magmoms = Magmom.get_consistent_set(magmoms)

        # convert to numpy array for convenience
        magmoms = np.array([list(magmom) for magmom in magmoms])
        magmoms = magmoms[np.any(magmoms, axis=1)] # remove zero magmoms
        if len(magmoms) == 0:
            return True

        # use first moment as reference to compare against
        ref_magmom = magmoms[0]
        # magnitude of cross products != 0 if non-collinear with reference
        num_ncl = np.count_nonzero(np.linalg.norm(np.cross(ref_magmom, magmoms), axis=1))
        if num_ncl > 0:
            return False
        else:
            return True

    @staticmethod
    def get_moments_relative_to_crystal_axes(magmoms, lattice):
        """
        Used for writing moments to magCIF file.
        If scalar magmoms, moments will be arbitrarily given along z.

        :param magmoms: list of magmoms (Magmoms, scalars, or vectors)
        :param lattice: Lattice
        :return: list of lists
        """
        return NotImplementedError

    def __getitem__(self, key):
        return self.moment[key]

    def __iter__(self):
        return iter(self.moment)

    def __abs__(self):
        return np.linalg.norm(self.moment)

    def __eq__(self, other):
        """
        Equal if 'global' magnetic moments are the same, saxis can differ.
        """
        #TODO: is this sensible?
        return self.global_moment == other.global_moment

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return abs(self) < abs(other)

    def __hash__(self):
        return (tuple(self.moment)+tuple(self.saxis)).__hash__()

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
        structures and might give non-sensical results except in the case
        of only slightly non-collinear structures (e.g. small canting).
        """
        return self.get_00t_magmom_with_xyz_saxis()[2]

    def __str__(self):
        return str(float(self))

    def __repr__(self):
        return 'Magnetic moment {0} (spin axis = {1})'.format(self.moment, self.saxis)