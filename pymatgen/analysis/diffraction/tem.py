from __future__ import division, print_function, unicode_literals
import json
import os
from fractions import Fraction
from typing import List, Dict, Tuple
import numpy as np
import plotly.graph_objs as go
import plotly.offline as poff
import scipy as sc
import scipy.constants as sc
from IPython.display import set_matplotlib_formats
from prettytable import PrettyTable
from pymatgen import Structure, Element
from pymatgen.analysis.diffraction.core import DiffractionPattern, AbstractDiffractionPatternCalculator, \
    get_unique_families

with open(os.path.join(os.path.dirname(__file__),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)
set_matplotlib_formats('retina')
poff.init_notebook_mode(connected=True)
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
# tempattern inherits from diffractionpattern

"""
This module implements a TEM pattern calculator.
"""
# Credit to Dr. Shyue Ping Ong for the template of the calculator
__author__ = "Frank Wan, modified by JasonL"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Frank Wan respect for S.P.O"
__email__ = "fwan@berkeley.edu, yhljason@berkeley.edu"
__date__ = "06/19/2019, updated 10/2019"


class TEMDot:
    """
    Instantatiates a point on the TEM diffraction pattern.
    """

    def __init__(self, position: List[float], hkl: List[int], intensity: float, film_radius: float,
                 d_spacing: float) -> None:
        """
        Args:
              hkl (3-tuple): The hkl plane that the point corresponds/is indexed to.
              d_spacing (float): The interplanar spacing of the dot.
              film_radius (float): The radius of the dot on the film. Determined
              by microscope aberration equations (ie Cs corrections and other such
              aberrations)
              intensity (float): The intensity of the dot. Determines its brightness
              in the pattern, relative to those of the other dots.
              position (Position): The xy-coordinates of the dot on the plot.
        """
        self.position = position
        self.hkl = hkl
        self.intensity = intensity
        self.film_radius = film_radius
        self.d_spacing = d_spacing


class TEMCalculator(AbstractDiffractionPatternCalculator):
    """
    Computes the TEM pattern of a crystal structure for multiple Laue zones.
    """
    wavelength_cache = {}

    def __init__(self, symprec: float = None, voltage: float = 200, beam_direction: List[int] = [0, 0, 1],
                 camera_length: int = 160, debye_waller_factors: Dict[str, float] = None, cs: float = 1) -> None:
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
        if any(type(n) is tuple for n in points):
            return points
        if len(points) == 0:
            return []
        filtered = np.where(np.dot(np.array(self.beam_direction), np.transpose(points)) == laue_zone)
        result = points[filtered]
        result_tuples = [tuple(x) for x in result.tolist()]
        return result_tuples

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
            -> Dict[Tuple[int, int, int], Dict]:
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
            -> Dict[Tuple[int, int, int], Dict]:
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
            -> DiffractionPattern:
        """
        Calculates the diffraction pattern for a structure. As you'll find out if you try to run this method,
        xrd-relevant info is tem-irrelevant. only included to satisfy the requirements of a subclass. Also,
        runtime is a bit long for this method.
        Args:
            structure (Structure): Input structure
            scaled (bool): Whether to return scaled intensities. The maximum
            peak is set to a value of 100. Defaults to True. Use False if
            you need the absolute values to combine XRD plots.
            two_theta_range ([float of length 2]): Tuple for range of
            two_thetas to calculate in degrees. Defaults to (0, 90). Set to
            None if you want all diffracted beams within the limiting
            sphere of radius 2 / wavelength.
        Returns:
            (XRDPattern)
        """
        points = self.generate_points(-10, 11)
        points_filtered = self.zone_axis_filter(points)
        interplanar_spacings = self.get_interplanar_spacings(structure, points_filtered)
        bragg_angles = self.bragg_angles(interplanar_spacings)
        cell_intensity = self.cell_intensity(structure, bragg_angles)
        max_intensity = max([v for v in cell_intensity.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        # creates a dict of 2thetas to cell_intensities
        xy_pairs = {}
        hkls_crude = []
        for plane in cell_intensity:
            xy_pairs[2 * bragg_angles[plane]] = cell_intensity[plane]
            hkls_crude.append(plane)
        for k in sorted(xy_pairs.keys()):
            v = xy_pairs[k]
            if v / max_intensity * 100 > AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v)
        fam = get_unique_families(hkls_crude)
        hkls.append([{"hkl": hkl, "multiplicity": mult}
                     for hkl, mult in fam.items()])
        for plane in fam:
            d_hkls.append(bragg_angles[plane])
        tem = DiffractionPattern(x, y, hkls, d_hkls)
        if scaled:
            tem.normalize(mode="max", value=100)
        return tem

    def normalized_cell_intensity(self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]) \
            -> Dict[Tuple[int, int, int], Dict]:
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

    def get_first_point(self, structure: Structure, points: list) -> Dict[Tuple[int, int, int], int]:
        """
        Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            A dict of a hkl plane to max interplanar distance.
        """
        max_d = -100
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
        p1 = list(first_point)
        p2 = list(second_point)
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
        Generates all TEM_dot objects that will appear on the 2D diffraction pattern.
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
            position = positions[plane]
            hkl = plane
            intensity = cell_intensity[plane]
            film_radius = 0.91 * (10 ** -3 * self.cs * self.wavelength_rel() ** 3) ** Fraction('1/4')
            d_spacing = interplanar_spacings[plane]
            dot = TEMDot(position, hkl, intensity, film_radius, d_spacing)
            dots.append(dot)
        return dots

    def show_plot_2d(self, structure: Structure):
        """
        Generates the 2D diffraction pattern of the input structure.
        Args:
            structure (Structure): The input structure.
        Returns:
            none (shows 2D DP)
        """
        points = self.generate_points(-10, 11)
        TEM_dots = self.TEM_dots(structure, points)
        xs = []
        ys = []
        hkls = []
        intensities = []

        for dot in TEM_dots:
            position = np.array([dot.position[0], dot.position[1]])
            xs.append(dot.position[0])
            ys.append(dot.position[1])
            hkls.append(dot.hkl)
            intensities.append(dot.intensity)

        data = [
            go.Scatter(
                x=xs,
                y=ys,
                text=hkls,
                hoverinfo='text',
                mode='markers',
                marker=dict(
                    size=8,
                    cmax=1,
                    cmin=0,
                    color=intensities,
                    colorbar=dict(
                        title='Colorbar',
                        yanchor='top'
                    ),
                    colorscale=[[0, 'black'], [1.0, 'white']]
                ),
                showlegend=False
            ), go.Scatter(
                x=[0],
                y=[0],
                text="(0, 0, 0): Direct beam",
                hoverinfo='text',
                mode='markers',
                marker=dict(
                    size=14,
                    cmax=1,
                    cmin=0,
                    color='white'
                ),
            )
        ]
        layout = go.Layout(
            title='2D Diffraction Pattern<br>Beam Direction: ' + ''.join(str(e) for e in self.beam_direction),
            font=dict(
                family='Comic Sans, monospace',
                size=18,
                color='#7f7f7f'),
            hovermode='closest',
            xaxis=dict(
                autorange=True,
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks='',
                showticklabels=False
            ),
            yaxis=dict(
                autorange=True,
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks='',
                showticklabels=False,
            ),
            width=600,
            height=600,
            paper_bgcolor='rgba(100,110,110,0.5)',
            plot_bgcolor='black'
        )

        fig = go.Figure(data=data, layout=layout)
        poff.iplot(fig, filename='stuff')

    def get_pattern_2d(self, structure: Structure) -> PrettyTable:
        """
        Returns all relevant TEM DP info in a PrettyTable.
        Args:
            structure (Structure): The input structure.
        Returns:
            PrettyTable
        """
        points = self.generate_points(-10, 11)
        TEM_dots = self.TEM_dots(structure, points)
        table = PrettyTable()
        table.field_names = ["Pos", "(hkl)", "Intnsty (norm)", "Film rad", "Interplanar Spacing"]
        for dot in TEM_dots:
            position = np.array([dot.position[0], dot.position[1]])
            table.add_row([position, dot.hkl, dot.intensity, dot.film_radius, dot.d_spacing])
        return table
