# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements a TEM pattern calculator.
"""

from __future__ import division, print_function, unicode_literals
import json
import os
from collections import namedtuple
from fractions import Fraction
from typing import List, Dict, Tuple
import numpy as np  # type: ignore
import scipy.constants as sc  # type: ignore
import pandas as pd  # type: ignore
from pymatgen import Structure, Element  # type: ignore
from pymatgen.analysis.diffraction.core import AbstractDiffractionPatternCalculator  # type: ignore

with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)

# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
# Credit to Dr. Shyue Ping Ong for the template of the calculator

__author__ = "Frank Wan, modified by Jason L"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Frank Wan respect for S.P.O"
__email__ = "fwan@berkeley.edu, yhljason@berkeley.edu"
__date__ = "06/19/2019, updated 10/2019"


class TEMCalculator(AbstractDiffractionPatternCalculator):
    """
    Computes the TEM pattern of a crystal structure for multiple Laue zones.
    """

    def __init__(self, wavelength_cache: Dict[float, float] = {}, symprec: float = None, voltage: float = 200,
                 beam_direction: List[int] = [0, 0, 1], camera_length: int = 160,
                 debye_waller_factors: Dict[str, float] = None, cs: float = 1) -> None:
        """
        Initializes the TEM calculator with a given radiation.
        Args:
            symprec (float): Symmetry precision for structure refinement. If
            set to 0, no refinement is done. Otherwise, refinement is
            performed using spglib with provided precision.
            voltage (float): The wavelength is a function of the TEM microscope's
            voltage. By default, set to 200 kV. Units in kV.
            beam_direction: The direction of the electron beam fired onto the sample.
            By default, set to [0,0,1], which corresponds to the normal direction
            of the sample plane.
            camera_length (int): The distance from the sample to the projected diffraction pattern.
            By default, set to 160 cm. Units in cm.
            debye_waller_factors ({element symbol: float}): Allows the
            specification of Debye-Waller factors. Note that these
            factors are temperature dependent.
            cs (float): the chromatic aberration coefficient. set by default to 1 mm.
            later on: may want "number of iterations", "magnification", "critical value of beam",
            "twin direction" for certain materials, "sample thickness", and "excitation error s"
        """
        self.wavelength_cache = wavelength_cache
        self.symprec = symprec
        self.voltage = voltage
        self.beam_direction = beam_direction
        self.camera_length = camera_length
        self.debye_waller_factors = debye_waller_factors or {}
        self.cs = cs

    def wavelength_rel(self) -> float:
        """
        Calculates the wavelength of the electron beam with relativstic kinematic effects taken
        into account (electrons are way faster than X-rays, so you can't neglect these effects).
        Args:
            none
        Returns:
            relativisticWavelength (in meters)
        """
        if self.voltage in self.wavelength_cache:
            return self.wavelength_cache[self.voltage]
        wavelength_rel = sc.h / np.sqrt(2 * sc.m_e * sc.e * 1000 * self.voltage *
                                        (1 + (sc.e * 1000 * self.voltage) / (2 * sc.m_e * sc.c ** 2)))
        self.wavelength_cache[self.voltage] = wavelength_rel
        return wavelength_rel

    def generate_points(self, coord_left: int = -10, coord_right: int = 10) -> np.ndarray:
        """
        Generates a bunch of 3D points that span a cube.
        Args:
            coord_left (int): The minimum coordinate value.
            coord_right (int): The maximum coordinate value.
        Returns:
            Numpy 2d array
        """
        points = [0, 0, 0]
        coord_values = np.arange(coord_left, coord_right + 1)
        points[0], points[1], points[2] = np.meshgrid(coord_values, coord_values, coord_values)
        points_matrix = (np.ravel(points[i]) for i in range(0, 3))
        result = np.vstack(list(points_matrix)).transpose()
        return result

    def zone_axis_filter(self, points: np.ndarray, laue_zone: int = 0) -> List[Tuple[int, int, int]]:
        """
        Filters out all points that exist within the specified Laue zone according to the zone axis rule.
        Args:
            points (np.ndarray): The list of points to be filtered.
            laue_zone (int): The desired Laue zone.
        Returns:
            list of 3-tuples
        """
        if any(isinstance(n, tuple) for n in points):
            return points
        if len(points) == 0:
            return []
        filtered = np.where(np.dot(np.array(self.beam_direction), np.transpose(points)) == laue_zone)
        result = points[filtered]
        result_tuples = [tuple(x) for x in result.tolist()]
        return result_tuples  # type: ignore

    def get_interplanar_spacings(self, structure: Structure, points: List[Tuple[int, int, int]]) \
            -> Dict[Tuple[int, int, int], float]:
        """
        Args:
            structure (Structure): the input structure.
            points (tuple): the desired hkl indices.
        Returns:
            Dict of hkl to its interplanar spacing (float).
        """
        points_filtered = self.zone_axis_filter(points)
        if (0, 0, 0) in points_filtered:
            points_filtered.remove((0, 0, 0))
        interplanar_spacings_val = np.array(list(map(lambda x: structure.lattice.d_hkl(x), points_filtered)))
        interplanar_spacings = dict(zip(points_filtered, interplanar_spacings_val))
        return interplanar_spacings

    def bragg_angles(self, interplanar_spacings: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], float]:
        """
        Gets the Bragg angles for every hkl point passed in (where n = 1).
        Args:
            interplanar_spacings (dict): dictionary of hkl to interplanar spacing
        Returns:
            dict of hkl plane (3-tuple) to Bragg angle in radians (float)
        """
        plane = list(interplanar_spacings.keys())
        interplanar_spacings_val = np.array(list(interplanar_spacings.values()))
        bragg_angles_val = np.arcsin(self.wavelength_rel() / (2 * interplanar_spacings_val))
        bragg_angles = dict(zip(plane, bragg_angles_val))
        return bragg_angles

    def get_s2(self, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], float]:
        """
        Calculates the s squared parameter (= square of sin theta over lambda) for each hkl plane.
        Args:
            bragg_angles (Dict): The bragg angles for each hkl plane.
        Returns:
            Dict of hkl plane to s2 parameter, calcs the s squared parameter (= square of sin theta over lambda).
        """
        plane = list(bragg_angles.keys())
        bragg_angles_val = np.array(list(bragg_angles.values()))
        s2_val = (np.sin(bragg_angles_val) / self.wavelength_rel()) ** 2
        s2 = dict(zip(plane, s2_val))
        return s2

    def x_ray_factors(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Element, Dict]:
        """
        Calculates x-ray factors, which are required to calculate atomic scattering factors. Method partially inspired
        by the equivalent process in the xrd module.
        Args:
            structure (Structure): The input structure.
            bragg_angles (Dict): Dictionary of hkl plane to Bragg angle.
        Returns:
            A Dict of atomic symbol to another dict of hkl plane to x-ray factor
        """
        x_ray_factors = {}
        s2 = self.get_s2(bragg_angles)
        atoms = structure.composition.elements
        scattering_factors_for_atom = {}
        for atom in atoms:
            coeffs = np.array(ATOMIC_SCATTERING_PARAMS[atom.symbol])
            for plane in bragg_angles:
                scattering_factor_curr = atom.Z - 41.78214 * s2[plane] * np.sum(coeffs[:, 0]
                                                                                * np.exp(-coeffs[:, 1] * s2[plane]),
                                                                                axis=None)
                scattering_factors_for_atom[plane] = scattering_factor_curr
            x_ray_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}
        return x_ray_factors

    def electron_scattering_factors(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Element, Dict]:
        """
        Calculates atomic scattering factors for electrons using the Mott-Bethe formula (1st order Born approximation).
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict from atomic symbol to another dict of hkl plane to factor
        """
        electron_scattering_factors = {}
        x_ray_factors = self.x_ray_factors(structure, bragg_angles)
        s2 = self.get_s2(bragg_angles)
        atoms = structure.composition.elements
        prefactor = sc.e / (16 * (np.pi ** 2) * sc.epsilon_0)
        scattering_factors_for_atom = {}
        for atom in atoms:
            for plane in bragg_angles:
                scattering_factor_curr = prefactor * (atom.Z - x_ray_factors[atom.symbol][plane]) / s2[plane]
                scattering_factors_for_atom[plane] = scattering_factor_curr
            electron_scattering_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}
        return electron_scattering_factors

    def cell_scattering_factors(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], int]:
        """
        Calculates the scattering factor for the whole cell.
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict of hkl plane (3-tuple) to scattering factor
        """
        cell_scattering_factors = {}
        electron_scattering_factors = self.electron_scattering_factors(structure, bragg_angles)
        scattering_factor_curr = 0
        for plane in bragg_angles:
            for site in structure:
                for sp, occu in site.species.items():
                    g_dot_r = np.dot(np.array(plane), np.transpose(site.frac_coords))
                    scattering_factor_curr += electron_scattering_factors[sp.symbol][plane] * np.exp(
                        2j * np.pi * g_dot_r)
            cell_scattering_factors[plane] = scattering_factor_curr
            scattering_factor_curr = 0
        return cell_scattering_factors

    def cell_intensity(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], float]:
        """
        Calculates cell intensity for each hkl plane. For simplicity's sake, take I = |F|**2.
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict of hkl plane to cell intensity
        """
        csf = self.cell_scattering_factors(structure, bragg_angles)
        plane = bragg_angles.keys()
        csf_val = np.array(list(csf.values()))
        cell_intensity_val = (csf_val * csf_val.conjugate()).real
        cell_intensity = dict(zip(plane, cell_intensity_val))
        return cell_intensity

    def get_pattern(self, structure: Structure, scaled: bool = True, two_theta_range: tuple = (0, 90)) \
            -> pd.DataFrame:
        """
            Returns all relevant TEM DP info in a pandas dataframe.
            Args:
                structure (Structure): The input structure.
            Returns:
                PandasDataFrame
        """
        points = self.generate_points(-10, 11)
        TEM_dots = self.TEM_dots(structure, points)
        field_names = ["Pos", "(hkl)", "Intnsty (norm)", "Film rad", "Interplanar Spacing"]
        rows_list = []
        for dot in TEM_dots:
            dict1 = {}
            dict1 = {'Pos': dot.position, '(hkl)': dot.hkl, 'Intnsty (norm)': dot.intensity,
                     'Film rad': dot.film_radius, 'Interplanar Spacing': dot.d_spacing}
            rows_list.append(dict1)
        df = pd.DataFrame(rows_list, columns=field_names)
        return df

    def normalized_cell_intensity(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], float]:
        """
        Normalizes the cell_intensity dict to 1, for use in plotting.
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict of hkl plane to normalized cell intensity
        """
        normalized_cell_intensity = {}
        cell_intensity = self.cell_intensity(structure, bragg_angles)
        max_intensity = max([v for v in cell_intensity.values()])
        norm_factor = 1 / max_intensity
        for plane in cell_intensity:
            normalized_cell_intensity[plane] = cell_intensity[plane] * norm_factor
        return normalized_cell_intensity

    def is_parallel(self, plane: Tuple[int, int, int], other_plane: Tuple[int, int, int]) \
            -> bool:
        """
        Checks if two hkl planes are parallel in reciprocal space.
        Args:
            plane (3-tuple): The first plane to be compared.
            other_plane (3-tuple): The other plane to be compared.
        Returns:
            boolean
        """
        return np.array_equal(np.cross(np.asarray(plane), np.asarray(other_plane)), np.array([0, 0, 0]))

    def get_first_point(self, structure: Structure, points: list) -> Dict[Tuple[int, int, int], float]:
        """
        Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            A dict of a hkl plane to max interplanar distance.
        """
        max_d = -100.0
        max_d_plane = (0, 0, 1)
        points = self.zone_axis_filter(points)
        spacings = self.get_interplanar_spacings(structure, points)
        for plane in sorted(spacings.keys()):
            if spacings[plane] > max_d:
                max_d_plane = plane
                max_d = spacings[plane]
        return {max_d_plane: max_d}

    def get_plot_coeffs(self, p1: Tuple[int, int, int], p2: Tuple[int, int, int], p3: Tuple[int, int, int],
                        denom: float, init_denom_0: bool) -> list:
        """
        Calculates coefficients of the vector addition required to generate positions for each DP point.
        Args:
            p1 (3-tuple): The first point. Fixed.
            p2 (3-tuple): The second point. Fixed.
            p3 (3-tuple): The point whose coefficients are to be calculted.
            denom (float): The denominator in the matrix calculation.
            init_denom_0 (boolean): Whether or not the first calculated denominator was 0.
        Returns:
            list of length 2 [x-coefficient, y-coefficient]
        """
        coeffs = []
        if init_denom_0:
            a_num = np.array([[p3[0], p3[2]], [p2[0], p2[2]]])
            b_num = np.array([[p1[0], p1[2]], [p3[0], p3[2]]])
        else:
            a_num = np.array([[p3[0], p3[1]], [p2[0], p2[1]]])
            b_num = np.array([[p1[0], p1[1]], [p3[0], p3[1]]])
        coeffs_0 = np.linalg.det(a_num) / denom
        coeffs_1 = np.linalg.det(b_num) / denom
        coeffs.append(coeffs_0)
        coeffs.append(coeffs_1)
        return coeffs

    def get_positions(self, structure: Structure, points: list) -> Dict[Tuple[int, int, int], list]:
        """
        Calculates all the positions of each hkl point in the 2D diffraction pattern. Distance in centimeters.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            dict of hkl plane to xy-coordinates.
        """
        positions = {}
        points = self.zone_axis_filter(points)
        # first is the max_d, min_r
        first_point_dict = self.get_first_point(structure, points)
        for point in first_point_dict:
            first_point = point
            first_d = first_point_dict[point]
        spacings = self.get_interplanar_spacings(structure, points)
        # second is the first non-parallel-to-first-point vector when sorted.
        # note 000 is "parallel" to every plane vector.
        for plane in sorted(spacings.keys()):
            second_point, second_d = plane, spacings[plane]
            if not self.is_parallel(first_point, second_point):
                break
        p1 = first_point
        p2 = second_point
        if (0, 0, 0) in points:
            points.remove((0, 0, 0))
        points.remove(first_point)
        points.remove(second_point)
        positions[(0, 0, 0)] = np.array([0, 0])
        r1 = 10 ** 10 * self.wavelength_rel() * self.camera_length / first_d
        positions[first_point] = np.array([r1, 0])
        r2 = 10 ** 10 * self.wavelength_rel() * self.camera_length / second_d
        phi = np.arccos(
            np.dot(p1, np.transpose(p2)) / (np.sqrt(np.dot(p1, np.transpose(p1)) * np.dot(p2, np.transpose(p2)))))
        positions[second_point] = np.array([r2 * np.cos(phi), r2 * np.sin(phi)])
        denom = np.linalg.det(np.array([[p1[0], p1[1]], [p2[0], p2[1]]]))
        init_denom_0 = (denom == 0)
        if init_denom_0:
            denom = np.linalg.det(np.array([[p1[0], p1[2]], [p2[0], p2[2]]]))
        for plane in points:
            coeffs = self.get_plot_coeffs(p1, p2, plane, denom, init_denom_0)
            pos = np.array([coeffs[0] * positions[first_point][0] + coeffs[1] * positions[second_point][0],
                            coeffs[0] * positions[first_point][1] + coeffs[1] * positions[second_point][1]])
            positions[plane] = pos
        points.append((0, 0, 0))
        points.append(first_point)
        points.append(second_point)

        return positions

    def TEM_dots(self, structure: Structure, points: list) -> list:
        """
        Generates all TEM_dot as named tuples that will appear on the 2D diffraction pattern.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            list of TEM_dots
        """
        dots = []
        interplanar_spacings = self.get_interplanar_spacings(structure, points)
        bragg_angles = self.bragg_angles(interplanar_spacings)
        cell_intensity = self.normalized_cell_intensity(structure, bragg_angles)
        positions = self.get_positions(structure, points)
        for plane in cell_intensity.keys():
            dot = namedtuple('TEM_dot', ['position', 'hkl', 'intensity', 'film_radius', 'd_spacing'])
            position = positions[plane]
            hkl = plane
            intensity = cell_intensity[plane]
            film_radius = 0.91 * (10 ** -3 * self.cs * self.wavelength_rel() ** 3) ** Fraction('1/4')
            d_spacing = interplanar_spacings[plane]
            TEM_dot = dot(position, hkl, intensity, film_radius, d_spacing)
            dots.append(TEM_dot)
        return dots
