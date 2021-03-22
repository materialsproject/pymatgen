# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
# Credit to Dr. Shyue Ping Ong for the template of the calculator

"""
This module implements a TEM pattern calculator.
"""

import json
import os
from collections import namedtuple
from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Tuple, cast, Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import scipy.constants as sc

from pymatgen.analysis.diffraction.core import AbstractDiffractionPatternCalculator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.string import latexify_spacegroup, unicodeify_spacegroup

with open(os.path.join(os.path.dirname(__file__), "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)

__author__ = "Frank Wan, Jason Liang"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.22"
__maintainer__ = "Jason Liang"
__email__ = "fwan@berkeley.edu, yhljason@berkeley.edu"
__date__ = "03/31/2020"


class TEMCalculator(AbstractDiffractionPatternCalculator):
    """
    Computes the TEM pattern of a crystal structure for multiple Laue zones.
    Code partially inspired from XRD calculation implementation. X-ray factor to electron factor
        conversion based on the International Table of Crystallography.
    #TODO: Could add "number of iterations", "magnification", "critical value of beam",
            "twin direction" for certain materials, "sample thickness", and "excitation error s"
    """

    def __init__(
        self,
        symprec: float = None,
        voltage: float = 200,
        beam_direction: Tuple[int, int, int] = (0, 0, 1),
        camera_length: int = 160,
        debye_waller_factors: Dict[str, float] = None,
        cs: float = 1,
    ) -> None:
        """
        Args:
            symprec (float): Symmetry precision for structure refinement. If
                set to 0, no refinement is done. Otherwise, refinement is
                performed using spglib with provided precision.
            voltage (float): The wavelength is a function of the TEM microscope's
                voltage. By default, set to 200 kV. Units in kV.
            beam_direction (tuple): The direction of the electron beam fired onto the sample.
                By default, set to [0,0,1], which corresponds to the normal direction
                of the sample plane.
            camera_length (int): The distance from the sample to the projected diffraction pattern.
                By default, set to 160 cm. Units in cm.
            debye_waller_factors ({element symbol: float}): Allows the
                specification of Debye-Waller factors. Note that these
                factors are temperature dependent.
            cs (float): the chromatic aberration coefficient. set by default to 1 mm.
        """
        self.symprec = symprec
        self.voltage = voltage
        self.beam_direction = beam_direction
        self.camera_length = camera_length
        self.debye_waller_factors = debye_waller_factors or {}
        self.cs = cs

    @lru_cache(1)
    def wavelength_rel(self) -> float:
        """
        Calculates the wavelength of the electron beam with relativistic kinematic effects taken
            into account.
        Args:
            none
        Returns:
            Relativistic Wavelength (in angstroms)
        """
        wavelength_rel = (
            sc.h
            / np.sqrt(
                2 * sc.m_e * sc.e * 1000 * self.voltage * (1 + (sc.e * 1000 * self.voltage) / (2 * sc.m_e * sc.c ** 2))
            )
            * (10 ** 10)
        )
        return wavelength_rel

    @staticmethod
    def generate_points(coord_left: int = -10, coord_right: int = 10) -> np.ndarray:
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

    def zone_axis_filter(
        self, points: Union[List[Tuple[int, int, int]], np.ndarray], laue_zone: int = 0
    ) -> Union[List[Tuple[int, int, int]]]:
        """
        Filters out all points that exist within the specified Laue zone according to the zone axis rule.
        Args:
            points (np.ndarray): The list of points to be filtered.
            laue_zone (int): The desired Laue zone.
        Returns:
            list of 3-tuples
        """
        if any(isinstance(n, tuple) for n in points):
            return list(points)
        if len(points) == 0:
            return []
        filtered = np.where(np.dot(np.array(self.beam_direction), np.transpose(points)) == laue_zone)
        result = points[filtered]
        result_tuples = cast(List[Tuple[int, int, int]], [tuple(x) for x in result.tolist()])
        return result_tuples

    def get_interplanar_spacings(
        self, structure: Structure, points: Union[List[Tuple[int, int, int]], np.ndarray]
    ) -> Dict[Tuple[int, int, int], float]:
        """
        Args:
            structure (Structure): the input structure.
            points (tuple): the desired hkl indices.
        Returns:
            Dict of hkl to its interplanar spacing, in angstroms (float).
        """
        points_filtered = self.zone_axis_filter(points)
        if (0, 0, 0) in points_filtered:
            points_filtered.remove((0, 0, 0))
        interplanar_spacings_val = np.array(list(map(lambda x: structure.lattice.d_hkl(x), points_filtered)))
        interplanar_spacings = dict(zip(points_filtered, interplanar_spacings_val))
        return interplanar_spacings

    def bragg_angles(
        self, interplanar_spacings: Dict[Tuple[int, int, int], float]
    ) -> Dict[Tuple[int, int, int], float]:
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

    def get_s2(self, bragg_angles: Dict[Tuple[int, int, int], float]) -> Dict[Tuple[int, int, int], float]:
        """
        Calculates the s squared parameter (= square of sin theta over lambda) for each hkl plane.
        Args:
            bragg_angles (Dict): The bragg angles for each hkl plane.
        Returns:
            Dict of hkl plane to s2 parameter, calculates the s squared parameter
                (= square of sin theta over lambda).
        """
        plane = list(bragg_angles.keys())
        bragg_angles_val = np.array(list(bragg_angles.values()))
        s2_val = (np.sin(bragg_angles_val) / self.wavelength_rel()) ** 2
        s2 = dict(zip(plane, s2_val))
        return s2

    def x_ray_factors(
        self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]
    ) -> Dict[str, Dict[Tuple[int, int, int], float]]:
        """
        Calculates x-ray factors, which are required to calculate atomic scattering factors. Method partially inspired
        by the equivalent process in the xrd module.
        Args:
            structure (Structure): The input structure.
            bragg_angles (Dict): Dictionary of hkl plane to Bragg angle.
        Returns:
            dict of atomic symbol to another dict of hkl plane to x-ray factor (in angstroms).
        """
        x_ray_factors = {}
        s2 = self.get_s2(bragg_angles)
        atoms = structure.composition.elements
        scattering_factors_for_atom = {}
        for atom in atoms:
            coeffs = np.array(ATOMIC_SCATTERING_PARAMS[atom.symbol])
            for plane in bragg_angles:
                scattering_factor_curr = atom.Z - 41.78214 * s2[plane] * np.sum(
                    coeffs[:, 0] * np.exp(-coeffs[:, 1] * s2[plane]), axis=None
                )
                scattering_factors_for_atom[plane] = scattering_factor_curr
            x_ray_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}
        return x_ray_factors

    def electron_scattering_factors(
        self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]
    ) -> Dict[str, Dict[Tuple[int, int, int], float]]:
        """
        Calculates atomic scattering factors for electrons using the Mott-Bethe formula (1st order Born approximation).
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict from atomic symbol to another dict of hkl plane to factor (in angstroms)
        """
        electron_scattering_factors = {}
        x_ray_factors = self.x_ray_factors(structure, bragg_angles)
        s2 = self.get_s2(bragg_angles)
        atoms = structure.composition.elements
        prefactor = 0.023934
        scattering_factors_for_atom = {}
        for atom in atoms:
            for plane in bragg_angles:
                scattering_factor_curr = prefactor * (atom.Z - x_ray_factors[atom.symbol][plane]) / s2[plane]
                scattering_factors_for_atom[plane] = scattering_factor_curr
            electron_scattering_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}
        return electron_scattering_factors

    def cell_scattering_factors(
        self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]
    ) -> Dict[Tuple[int, int, int], int]:
        """
        Calculates the scattering factor for the whole cell.
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.
        Returns:
            dict of hkl plane (3-tuple) to scattering factor (in angstroms).
        """
        cell_scattering_factors = {}
        electron_scattering_factors = self.electron_scattering_factors(structure, bragg_angles)
        scattering_factor_curr = 0
        for plane in bragg_angles:
            for site in structure:
                for sp, occu in site.species.items():
                    g_dot_r = np.dot(np.array(plane), np.transpose(site.frac_coords))
                    scattering_factor_curr += electron_scattering_factors[sp.symbol][plane] * np.exp(
                        2j * np.pi * g_dot_r
                    )
            cell_scattering_factors[plane] = scattering_factor_curr
            scattering_factor_curr = 0
        return cell_scattering_factors

    def cell_intensity(
        self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]
    ) -> Dict[Tuple[int, int, int], float]:
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

    def get_pattern(
        self,
        structure: Structure,
        scaled: bool = None,
        two_theta_range: Tuple[float, float] = None,
    ) -> pd.DataFrame:
        """
        Returns all relevant TEM DP info in a pandas dataframe.
        Args:
            structure (Structure): The input structure.
            scaled (boolean): Required value for inheritance, does nothing in TEM pattern
            two_theta_range (Tuple): Required value for inheritance, does nothing in TEM pattern
        Returns:
            PandasDataFrame
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()
        points = self.generate_points(-10, 11)
        tem_dots = self.tem_dots(structure, points)
        field_names = [
            "Position",
            "(hkl)",
            "Intensity (norm)",
            "Film radius",
            "Interplanar Spacing",
        ]
        rows_list = []
        for dot in tem_dots:
            dict1 = {
                "Pos": dot.position,
                "(hkl)": dot.hkl,
                "Intnsty (norm)": dot.intensity,
                "Film rad": dot.film_radius,
                "Interplanar Spacing": dot.d_spacing,
            }
            rows_list.append(dict1)
        df = pd.DataFrame(rows_list, columns=field_names)
        return df

    def normalized_cell_intensity(
        self, structure: Structure, bragg_angles: Dict[Tuple[int, int, int], float]
    ) -> Dict[Tuple[int, int, int], float]:
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
        max_intensity = max(cell_intensity.values())
        norm_factor = 1 / max_intensity
        for plane in cell_intensity:
            normalized_cell_intensity[plane] = cell_intensity[plane] * norm_factor
        return normalized_cell_intensity

    def is_parallel(
        self,
        structure: Structure,
        plane: Tuple[int, int, int],
        other_plane: Tuple[int, int, int],
    ) -> bool:
        """
        Checks if two hkl planes are parallel in reciprocal space.
        Args:
            structure (Structure): The input structure.
            plane (3-tuple): The first plane to be compared.
            other_plane (3-tuple): The other plane to be compared.
        Returns:
            boolean
        """
        phi = self.get_interplanar_angle(structure, plane, other_plane)
        return phi in (180, 0) or np.isnan(phi)

    def get_first_point(self, structure: Structure, points: list) -> Dict[Tuple[int, int, int], float]:
        """
        Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            dict of a hkl plane to max interplanar distance.
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

    @staticmethod
    def get_interplanar_angle(structure: Structure, p1: Tuple[int, int, int], p2: Tuple[int, int, int]) -> float:
        """
        Returns the interplanar angle (in degrees) between the normal of two crystal planes.
        Formulas from International Tables for Crystallography Volume C pp. 2-9.
        Args:
            structure (Structure): The input structure.
            p1 (3-tuple): plane 1
            p2 (3-tuple): plane 2
        Returns:
            float
        """
        a, b, c = structure.lattice.a, structure.lattice.b, structure.lattice.c
        alpha, beta, gamma = (
            np.deg2rad(structure.lattice.alpha),
            np.deg2rad(structure.lattice.beta),
            np.deg2rad(structure.lattice.gamma),
        )
        v = structure.lattice.volume
        a_star = b * c * np.sin(alpha) / v
        b_star = a * c * np.sin(beta) / v
        c_star = a * b * np.sin(gamma) / v
        cos_alpha_star = (np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / (np.sin(beta) * np.sin(gamma))
        cos_beta_star = (np.cos(alpha) * np.cos(gamma) - np.cos(beta)) / (np.sin(alpha) * np.sin(gamma))
        cos_gamma_star = (np.cos(alpha) * np.cos(beta) - np.cos(gamma)) / (np.sin(alpha) * np.sin(beta))
        r1_norm = np.sqrt(
            p1[0] ** 2 * a_star ** 2
            + p1[1] ** 2 * b_star ** 2
            + p1[2] ** 2 * c_star ** 2
            + 2 * p1[0] * p1[1] * a_star * b_star * cos_gamma_star
            + 2 * p1[0] * p1[2] * a_star * c_star * cos_beta_star
            + 2 * p1[1] * p1[2] * b_star * c_star * cos_gamma_star
        )
        r2_norm = np.sqrt(
            p2[0] ** 2 * a_star ** 2
            + p2[1] ** 2 * b_star ** 2
            + p2[2] ** 2 * c_star ** 2
            + 2 * p2[0] * p2[1] * a_star * b_star * cos_gamma_star
            + 2 * p2[0] * p2[2] * a_star * c_star * cos_beta_star
            + 2 * p2[1] * p2[2] * b_star * c_star * cos_gamma_star
        )
        r1_dot_r2 = (
            p1[0] * p2[0] * a_star ** 2
            + p1[1] * p2[1] * b_star ** 2
            + p1[2] * p2[2] * c_star ** 2
            + (p1[0] * p2[1] + p2[0] * p1[1]) * a_star * b_star * cos_gamma_star
            + (p1[0] * p2[2] + p2[0] * p1[1]) * a_star * c_star * cos_beta_star
            + (p1[1] * p2[2] + p2[1] * p1[2]) * b_star * c_star * cos_alpha_star
        )
        phi = np.arccos(r1_dot_r2 / (r1_norm * r2_norm))
        return np.rad2deg(phi)

    @staticmethod
    def get_plot_coeffs(
        p1: Tuple[int, int, int],
        p2: Tuple[int, int, int],
        p3: Tuple[int, int, int],
    ) -> np.ndarray:
        """
        Calculates coefficients of the vector addition required to generate positions for each DP point
        by the Moore-Penrose inverse method.
        Args:
            p1 (3-tuple): The first point. Fixed.
            p2 (3-tuple): The second point. Fixed.
            p3 (3-tuple): The point whose coefficients are to be calculted.
        Returns:
            Numpy array
        """
        a = np.array([[p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]]])
        b = np.array([[p3[0], p3[1], p3[2]]]).T
        a_pinv = np.linalg.pinv(a)
        x = np.dot(a_pinv, b)
        return np.ravel(x)

    def get_positions(self, structure: Structure, points: list) -> Dict[Tuple[int, int, int], np.ndarray]:
        """
        Calculates all the positions of each hkl point in the 2D diffraction pattern by vector addition.
        Distance in centimeters.
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
            if not self.is_parallel(structure, first_point, second_point):
                break
        p1 = first_point
        p2 = second_point
        if (0, 0, 0) in points:
            points.remove((0, 0, 0))
        points.remove(first_point)
        points.remove(second_point)
        positions[(0, 0, 0)] = np.array([0, 0])
        r1 = self.wavelength_rel() * self.camera_length / first_d
        positions[first_point] = np.array([r1, 0])
        r2 = self.wavelength_rel() * self.camera_length / second_d
        phi = np.deg2rad(self.get_interplanar_angle(structure, first_point, second_point))
        positions[second_point] = np.array([r2 * np.cos(phi), r2 * np.sin(phi)])
        for plane in points:
            coeffs = self.get_plot_coeffs(p1, p2, plane)
            pos = np.array(
                [
                    coeffs[0] * positions[first_point][0] + coeffs[1] * positions[second_point][0],
                    coeffs[0] * positions[first_point][1] + coeffs[1] * positions[second_point][1],
                ]
            )
            positions[plane] = pos
        points.append((0, 0, 0))
        points.append(first_point)
        points.append(second_point)

        return positions

    def tem_dots(self, structure: Structure, points) -> List:
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
            dot = namedtuple("dot", ["position", "hkl", "intensity", "film_radius", "d_spacing"])
            position = positions[plane]
            hkl = plane
            intensity = cell_intensity[plane]
            film_radius = 0.91 * (10 ** -3 * self.cs * self.wavelength_rel() ** 3) ** Fraction("1/4")
            d_spacing = interplanar_spacings[plane]
            tem_dot = dot(position, hkl, intensity, film_radius, d_spacing)
            dots.append(tem_dot)
        return dots

    def get_plot_2d(self, structure: Structure) -> go.Figure:
        """
        Generates the 2D diffraction pattern of the input structure.
        Args:
            structure (Structure): The input structure.
        Returns:
            Figure
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()
        points = self.generate_points(-10, 11)
        tem_dots = self.tem_dots(structure, points)
        xs = []
        ys = []
        hkls = []
        intensities = []
        for dot in tem_dots:
            xs.append(dot.position[0])
            ys.append(dot.position[1])
            hkls.append(str(dot.hkl))
            intensities.append(dot.intensity)
        hkls = list(map(unicodeify_spacegroup, list(map(latexify_spacegroup, hkls))))
        data = [
            go.Scatter(
                x=xs,
                y=ys,
                text=hkls,
                hoverinfo="text",
                mode="markers",
                marker=dict(
                    size=8,
                    cmax=1,
                    cmin=0,
                    color=intensities,
                    colorscale=[[0, "black"], [1.0, "white"]],
                ),
                showlegend=False,
            ),
            go.Scatter(
                x=[0],
                y=[0],
                text="(0, 0, 0): Direct beam",
                hoverinfo="text",
                mode="markers",
                marker=dict(size=14, cmax=1, cmin=0, color="white"),
                showlegend=False,
            ),
        ]
        layout = go.Layout(
            title="2D Diffraction Pattern<br>Beam Direction: " + "".join(str(e) for e in self.beam_direction),
            font=dict(size=14, color="#7f7f7f"),
            hovermode="closest",
            xaxis=dict(
                range=[-4, 4],
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks="",
                showticklabels=False,
            ),
            yaxis=dict(
                range=[-4, 4],
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks="",
                showticklabels=False,
            ),
            width=550,
            height=550,
            paper_bgcolor="rgba(100,110,110,0.5)",
            plot_bgcolor="black",
        )
        fig = go.Figure(data=data, layout=layout)
        return fig

    def get_plot_2d_concise(self, structure: Structure) -> go.Figure:
        """
        Generates the concise 2D diffraction pattern of the input structure of a smaller size and without layout.
        Does not display.
        Args:
            structure (Structure): The input structure.
        Returns:
            Figure
        """
        if self.symprec:
            finder = SpacegroupAnalyzer(structure, symprec=self.symprec)
            structure = finder.get_refined_structure()
        points = self.generate_points(-10, 11)
        tem_dots = self.tem_dots(structure, points)
        xs = []
        ys = []
        hkls = []
        intensities = []
        for dot in tem_dots:
            if dot.hkl != (0, 0, 0):
                xs.append(dot.position[0])
                ys.append(dot.position[1])
                hkls.append(dot.hkl)
                intensities.append(dot.intensity)
        data = [
            go.Scatter(
                x=xs,
                y=ys,
                text=hkls,
                mode="markers",
                hoverinfo="skip",
                marker=dict(
                    size=4,
                    cmax=1,
                    cmin=0,
                    color=intensities,
                    colorscale=[[0, "black"], [1.0, "white"]],
                ),
                showlegend=False,
            )
        ]
        layout = go.Layout(
            xaxis=dict(
                range=[-4, 4],
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks="",
                showticklabels=False,
            ),
            yaxis=dict(
                range=[-4, 4],
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks="",
                showticklabels=False,
            ),
            plot_bgcolor="black",
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
            width=121,
            height=121,
        )
        fig = go.Figure(data=data, layout=layout)
        fig.layout.update(showlegend=False)
        return fig
