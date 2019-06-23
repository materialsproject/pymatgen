from __future__ import division, unicode_literals
import pymatgen as pmg
from pymatgen import Lattice, Structure
from prettytable import PrettyTable
import matplotlib.colors as mc
import colorsys

from typing import List, Set, Dict, Tuple, Optional

import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline as poff
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor, Button
from IPython.display import set_matplotlib_formats, Image, display
set_matplotlib_formats('retina')
poff.init_notebook_mode(connected=True)

# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from math import *
from fractions import Fraction 

import os
import json
import math
import numpy as np
import random

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from pymatgen.analysis.diffraction.core import DiffractionPattern, AbstractDiffractionPatternCalculator, \
    get_unique_families

from pymatgen.analysis.local_env import site_is_of_motif_type, MinimumDistanceNN

with open(os.path.join(os.path.dirname("__file__"),
                       "atomic_scattering_params.json")) as f:
    ATOMIC_SCATTERING_PARAMS = json.load(f)

"""
This module implements a TEM pattern calculator.
"""

#Credit to Dr. Shyue Ping Ong for the template of the calculator
__author__ = "Frank Wan"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Frank Wan respect for S.P.O"
__email__ = "fwan@berkeley.edu"
__date__ = "06/19/2019"

class TEMDot():
    """
    Instantatiates a point on the TEM diffraction pattern.
    """
# TODO: make this MSONable
# TODO: add type hints for every method
# TODO: memoize methods to improve runtime
# TODO: standardize input argument convention (i.e., should things like bragg_angles be passed in, or calculated ad hoc?)
# TODO: add support for showing different Laue zones.
    def __init__(self, position: List[float], hkl: List[int], intensity: float, film_radius: float, d_spacing: float) -> None:
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

    """
    some important constants needed for calculations
    """
    planck = 6.62607004 * 10**-34
    restMassElec = 9.1093835611 * 10**-31
    chargeElec = 1.60217662 * 10**-19 
    speedLight = 299792458
    vacuumPerm = 8.8541878176 * 10**-12
    wavelength_cache = {}
    
    def __init__(self, symprec: float=None, voltage: float=300, beam_direction: List[int]=[0,0,1], 
                 camera_length: int=160, debye_waller_factors: Dict[str, float]=None, cs: float=1) -> None:
        """
        Initializes the TEM calculator with a given radiation.
        Args:
            symprec (float): Symmetry precision for structure refinement. If
                set to 0, no refinement is done. Otherwise, refinement is
                performed using spglib with provided precision.
            voltage (float): The wavelength is a function of the TEM microscope's
            	voltage. By default, set to 300 kV. Units in kV.
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
        wavelength_rel = self.planck/math.sqrt(2*self.restMassElec*self.chargeElec*1000*self.voltage*
        	(1+(self.chargeElec*1000*self.voltage)/(2*self.restMassElec*self.speedLight**2)))
        self.wavelength_cache[self.voltage] = wavelength_rel
        return wavelength_rel

    def generate_points(self, coord_left: int = -10, coord_right: int = 10) -> List[Tuple[int, int, int]]:
        """
        Generates a bunch of 3D points that span a cube.
        Args:
            coord_left (int): The minimum coordinate value.
            coord_right (int): The maximum coordinate value.
        Returns:
            list of 3-tuples
        """
        points = []
        coord_values = range(coord_left, coord_right + 1)
        for x in coord_values:
            for y in coord_values:
                for z in coord_values:
                    points.append((x,y,z))
        return points
            
    def zone_axis_filter(self, points: List[Tuple[int, int, int]], laue_zone: int = 0) -> List[Tuple[int, int, int]]:
        """
        Filters out all points that exist within the specified Laue zone according to the zone axis rule.
        Args:
            points (List[Tuple[int, int, int]]): The list of points to be filtered.
            laue_zone (int): The desired Laue zone.
        Returns:
            list of 3-tuples
        """
        # TODO: memoize this as a dict of laue zone to pointlist. prob want a system where if you input a diff range of 
        # initial points, you edit the cache entry.

        observed_points = []
        for point in points:
            if np.dot(self.beam_direction, np.transpose(point)) == laue_zone:
                observed_points.append(point)
        return observed_points

    def dist_cubic(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in a cubic crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            cubic interplanar distance (float)
        """
        return structure.lattice.a/sqrt(np.dot(point, np.transpose(point)))

    def dist_tetragonal(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in a tetragonal crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            tetragonal interplanar distance (float)
        """
        h = point[0]
        k = point[1]
        l = point[2]
        a = structure.lattice.a
        c = structure.lattice.c
        return (((h*h + k*k) / (a*a)) + (l*l)/(c*c))**(-0.5)

    def dist_hexagonal(self, structure, point):
        h = point[0]
        k = point[1]
        l = point[2]
        a = structure.lattice.a
        c = structure.lattice.c
        return ((4/3) * ((h*h + h*k + l*l) / (a*a)) + ((l*l)/(c*c))) ** (-0.5)
    
    def dist_rhombohedral(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in a rhombohedral crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            rhombohedral interplanar distance (float)
        """
        h = point[0]
        k = point[1]
        l = point[2]
        a = structure.lattice.a
        alpha = structure.lattice.alpha
        term_1 = np.dot(point, np.transpose(point)) * (np.sin(alpha))**2
        term_2 = 2 * (h*k + k*l + h*l) * ((np.cos(alpha))**2 - np.cos(alpha))
        term_3 = a*a*(1 - 3*(np.cos(alpha))**2 + 2 * (np.cos(alpha))**3)
        return ((term_1 + term_2) / term_3) ** (-0.5)
        
    def dist_orthorhombic(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in an orthorhombic crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            rhombohedral interplanar distance (float)
        """
        h_2 = point[0] * point[0]
        k_2 = point[1] * point[1]
        l_2 = point[2] * point[2]
        a_2 = structure.lattice.a * structure.lattice.a
        b_2 = structure.lattice.b * structure.lattice.b
        c_2 = structure.lattice.c * structure.lattice.c
        return ((h_2 / a_2) + (k_2 / b_2) + (l_2 / c_2))**(-0.5)
    
    def dist_monoclinic(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in a monoclinic crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            monoclinic interplanar distance (float)
        """
        h = point[0]
        k = point[1]
        l = point[2]
        a = structure.lattice.a
        b = structure.lattice.b
        c = structure.lattice.c
        alpha = structure.lattice.alpha
        beta = structure.lattice.beta
        gamma = structure.lattice.gamma
        return ((((h*h) / (a*a)) + ((k*k*(np.sin(beta))**2)/(b*b)) + ((l*l)/(c*c)) - ((2*h*l*np.cos(beta))/(a*c))) * (1/ ((np.sin(beta))**2)))**(-0.5)
        
    def dist_triclinic(self, structure, point):
        """
        Calculates the interplanar distance of an hkl plane in a triclinic crystal.
        Args:
            structure (Structure): The structure in question.
            point (3-tuple): The hkl plane in question.
        Returns:
            triclinic interplanar distance (float)
        """
        h = point[0]
        k = point[1]
        l = point[2]
        a = structure.lattice.a
        b = structure.lattice.b
        c = structure.lattice.c
        alpha = structure.lattice.alpha
        beta = structure.lattice.beta
        gamma = structure.lattice.gamma
        elm_1 = (h * np.sin(alpha))/a
        elm_2 = (k * np.sin(beta))/b
        elm_3 = (l * np.sin(gamma))/c
        vect_1 = (elm_1, elm_2, elm_3)
        term_1 = np.dot(vect_1, np.transpose(vect_1))
        
        term_2 = ((2*k*l)/(b*c)) * (np.cos(beta)*np.cos(gamma) - np.cos(alpha))
        term_3 = ((2*h*l)/(a*c)) * (np.cos(gamma)*np.cos(alpha) - np.cos(beta))
        term_4 = ((2*h*k)/(a*b)) * (np.cos(alpha)*np.cos(beta) - np.cos(gamma))
        numerator = term_1 + term_2 + term_3 + term_4
        
        denominator = 1 - np.cos(alpha)*np.cos(alpha) - np.cos(beta)*np.cos(beta) - np.cos(gamma)*np.cos(gamma) + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) 
        return (numerator / denominator)**(-0.5)
        
    def get_interplanar_spacings(self, structure: Structure, points: List[Tuple[int, int, int]]) -> Dict[Tuple[int, int, int], float]:
        """
        Gets the interplanar spacings for every hkl point passed in.
        Args:
            structure (Structure): The structure in question.
            points (3-tuple list): The hkl points in question.
        Returns:
            dict of hkl plane (3-tuple) to interplanar spacing (float)
        """        
        interplanar_spacings = {}
        analy = SpacegroupAnalyzer(structure, 0.1)
        points_filtered = self.zone_axis_filter(points)
        for point in points_filtered:
            if point != (0,0,0):
                if analy.get_lattice_type() == 'cubic':
                    interplanar_spacings[point] = self.dist_cubic(structure, point)
                if analy.get_lattice_type() == 'tetragonal':
                    interplanar_spacings[point] = self.dist_tetragonal(structure, point)
                if analy.get_lattice_type() == 'hexagonal':
                    interplanar_spacings[point] = self.dist_hexagonal(structure, point)
                if analy.get_lattice_type() == 'rhombohedral':
                    interplanar_spacings[point] = self.dist_rhombohedral(structure, point)
                if analy.get_lattice_type() == 'orthorhombic':
                    interplanar_spacings[point] = self.dist_orthorhombic(structure, point)
                if analy.get_lattice_type() == 'monoclinic':
                    interplanar_spacings[point] = self.dist_monoclinic(structure, point)
                if analy.get_lattice_type() == 'triclinic':
                    interplanar_spacings[point] = self.dist_triclinic(structure, point)
        return interplanar_spacings   
        
    def bragg_angles(self, interplanar_spacings: Dict[Tuple[int, int, int], float]) -> Dict[Tuple[int, int, int], float]:
        """
        Gets the Bragg angles for every hkl point passed in (where n = 1).
        Args:
            structure (Structure): The structure in question.
            points (3-tuple list): The hkl points in question.
        Returns:
            dict of hkl plane (3-tuple) to Bragg angle in radians (float)
        """        
        bragg_angles = {}
        for plane in interplanar_spacings:
            bragg_angles[plane] = np.arcsin(self.wavelength_rel()/(2*interplanar_spacings[plane]))
        return bragg_angles    
    
    def get_atoms(self, structure: Structure):
        """
        Gets all unique atoms present in the input structure.
        Args:
            structure (Structure): The structure in question.
        Returns:
            list of Atom objects
        """
        atoms = []
        for site in structure:
            for sp, occu in site.species.items():
                atoms.append(sp)
        return atoms
    
    def get_s2(self, structure, bragg_angles):
        """
        Calculates the s squared parameter (= square of sin theta over lambda) for each hkl plane.
        Args:
            structure (Structure): The structure in question.
            bragg_angles (Dict): The bragg angles for each hkl plane.
        Returns:
            Dict of hkl plane to s2 parameter
        """
        #calcs the s squared parameter (= square of sin theta over lambda).
        s2 = {}
        
        for plane in bragg_angles:
            s2[plane] = (np.sin(bragg_angles[plane])/self.wavelength_rel())**2
        return s2
    
    def x_ray_factors(self, structure, bragg_angles):
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
        s2 = self.get_s2(structure, bragg_angles)
        atoms = self.get_atoms(structure)
        coeffs = []
        scattering_factors_for_atom = {}
        scattering_factor_curr = 0

        for atom in atoms:
            coeffs = np.array(ATOMIC_SCATTERING_PARAMS[atom.symbol])
            for plane in bragg_angles:
                scattering_factor_curr = atom.Z - 41.78214 * s2[plane] * np.sum(coeffs[:, 0] 
                                                                                * np.exp(-coeffs[:, 1] * s2[plane]), axis=None)
                scattering_factors_for_atom[plane] = scattering_factor_curr
            x_ray_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}        
        return x_ray_factors
                
    def electron_scattering_factors(self, structure, bragg_angles):
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
        s2 = self.get_s2(structure, bragg_angles)
        atoms = self.get_atoms(structure)
        prefactor = self.chargeElec/(16*(pi**2)*self.vacuumPerm) 
        
        scattering_factors_for_atom = {}
        scattering_factor_curr = 0
        for atom in atoms:
            for plane in bragg_angles:
                scattering_factor_curr = prefactor * (atom.Z - x_ray_factors[atom.symbol][plane]) / s2[plane]
                scattering_factors_for_atom[plane] = scattering_factor_curr
            electron_scattering_factors[atom.symbol] = scattering_factors_for_atom
            scattering_factors_for_atom = {}
        return electron_scattering_factors
    
    def cell_scattering_factors(self, structure, bragg_angles):
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
        s2 = self.get_s2(structure, bragg_angles)
        atoms = self.get_atoms(structure)
        scattering_factor_curr = 0
        
        for plane in bragg_angles:
            for site in structure:
                for sp, occu in site.species.items(): #depending on how this iterates it may increase scatt by factor of 2.
                    g_dot_r = np.dot(np.array(plane), np.transpose(site.frac_coords))
                    scattering_factor_curr += electron_scattering_factors[sp.symbol][plane] * np.exp(2j * pi * g_dot_r)
            cell_scattering_factors[plane] = scattering_factor_curr
            scattering_factor_curr = 0
        return cell_scattering_factors    
    
    def cell_intensity(self, structure, bragg_angles):
        """
        Calculates cell intensity for each hkl plane. For simplicity's sake, take I = |F|**2.
        Args:
            structure (Structure): The input structure.
            bragg_angles (dict of 3-tuple to float): The Bragg angles for each hkl plane.  
        Returns:
            dict of hkl plane to cell intensity 
        """
        cell_intensity = {}
        csf = self.cell_scattering_factors(structure, bragg_angles)
        for plane in bragg_angles:
            cell_intensity[plane] = (csf[plane] * csf[plane].conjugate()).real
        return cell_intensity
    
    def get_pattern(self, structure, scaled=True, two_theta_range=(0, 90)):
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
        points = self.generate_points(-10,11)
        points_filtered = self.zone_axis_filter(points)
        interplanar_spacings = self.get_interplanar_spacings(structure, points_filtered)
        bragg_angles = self.bragg_angles(interplanar_spacings)
        cell_intensity = self.cell_intensity(structure, bragg_angles)
        
        max_intensity = max([v for v in cell_intensity.values()])
        x = []
        y = []
        hkls = []
        d_hkls = []
        
        #creates a dict of 2thetas to cell_intensities
        xy_pairs = {}
        hkls_crude = []
        
        for plane in cell_intensity:
            xy_pairs[2*bragg_angles[plane]] = cell_intensity[plane]
            hkls_crude.append(plane)
        
        for k in sorted(xy_pairs.keys()):
            v = xy_pairs[k]
            fam = get_unique_families(hkls_crude)
            if v / max_intensity * 100 > AbstractDiffractionPatternCalculator.SCALED_INTENSITY_TOL:
                x.append(k)
                y.append(v)
        hkls.append([{"hkl": hkl, "multiplicity": mult}
                     for hkl, mult in fam.items()])
        for plane in fam:
            d_hkls.append(bragg_angles[plane])
            
        tem = DiffractionPattern(x, y, hkls, d_hkls)
        if scaled:
            tem.normalize(mode="max", value=100)      
        return tem

    def normalized_cell_intensity(self, structure, bragg_angles):
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
    
    def is_parallel(self, plane, other_plane):
        """
        Checks if two hkl planes are parallel in reciprocal space.
        Args:
            plane (3-tuple): The first plane to be compared.
            other_plane (3-tuple): The other plane to be compared.
        Returns:
            boolean
        """
        return np.array_equal(np.cross(np.asarray(plane), np.asarray(other_plane)), np.array([0,0,0]))

    def get_first_point(self, structure, points):
        """
        Gets the first point to be plotted in the 2D DP, corresponding to maximum d/minimum R.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            dict of hkl plane to interplanar distance.
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
    
    def get_plot_coeffs(self, p1, p2, p3, denom, init_denom_0):
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
        if (init_denom_0):
            a_num = np.array([[p3[0], p3[2]], [p2[0],p2[2]]])
            b_num = np.array([[p1[0],p1[2]], [p3[0], p3[2]]])  
        else:    
            a_num = np.array([[p3[0], p3[1]], [p2[0],p2[1]]])
            b_num = np.array([[p1[0],p1[1]], [p3[0], p3[1]]])
        coeffs_0 = np.linalg.det(a_num) / denom
        coeffs_1 = np.linalg.det(b_num) / denom
        coeffs.append(coeffs_0)
        coeffs.append(coeffs_1)
        return coeffs
        
    def get_positions(self, structure, points):
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
        #first is the max_d, min_r
        first_point_dict = self.get_first_point(structure, points)
        for point in first_point_dict:
            first_point = point
            first_d = first_point_dict[point]
        spacings = self.get_interplanar_spacings(structure, points)
        
        #second is the first non-parallel-to-first-point vector when sorted. note 000 is "parallel" to every plane vector.
        for plane in sorted(spacings.keys()):
            second_point, second_d = plane, spacings[plane]
            if not self.is_parallel(first_point, second_point):
                break

        p1 = list(first_point)
        p2 = list(second_point)
        points.remove((0,0,0))            
        points.remove(first_point)
        points.remove(second_point)  


        positions[(0, 0, 0)] = np.array([0,0])
        
        #factor of 10**10 needed because first_d is in Angstroms (since first_d's calc is with lattice parameter which
        #in pymatgen is angstroms by default). WLoG, put first point on x-axis
        r1 = 10**10 * self.wavelength_rel() * self.camera_length / first_d
        positions[first_point] = np.array([r1,0])  

        #gets position of the second point. WLoG, assume it is located an angle phi (calculated by formula below) 
        #counterclockwise to the first point.  
        r2 = 10**10 * self.wavelength_rel() * self.camera_length / second_d
        phi = np.arccos(np.dot(p1, np.transpose(p2)) / (sqrt(np.dot(p1, np.transpose(p1)) * np.dot(p2, np.transpose(p2)))))
        positions[second_point] = np.array([r2*np.cos(phi), r2*np.sin(phi)])

        #in theory you have to check satisfaction of z3 = a*z1 + b*z2. in practice, the "physical realness" of
        #electron diffraction "ensures" that you don't. 
        #you also HAVE to make sure your denominator is nonzero. that is, that x and y for one of the two points are NOT
        #both zero. if one of them is, then then do this function again but with xz/yz coords and not xy. by the physical
        #realness stated above, this MUST work.
        denom = np.linalg.det(np.array([[p1[0],p1[1]], [p2[0],p2[1]]]))
        init_denom_0 = (denom == 0)
        if (init_denom_0):
            denom = np.linalg.det(np.array([[p1[0],p1[2]], [p2[0],p2[2]]]))
        for plane in points:
            coeffs = self.get_plot_coeffs(p1, p2, plane, denom, init_denom_0)
            pos = np.array([coeffs[0]*positions[first_point][0] + coeffs[1]*positions[second_point][0], 
                            coeffs[0]*positions[first_point][1] + coeffs[1]*positions[second_point][1]])
            positions[plane] = pos           
        points.append((0,0,0))            
        points.append(first_point)
        points.append(second_point) 
                
        return positions
    
    def TEM_dots(self, structure, points):
        """
        Generates all TEM_dot objects that will appear on the 2D diffraction pattern.
        Args:
            structure (Structure): The input structure.
            points (list): All points to be checked.
        Returns:
            list of TEM_dots
        """
        dots = []
        points_filtered = self.zone_axis_filter(points)
        interplanar_spacings = self.get_interplanar_spacings(structure, points_filtered)
        bragg_angles = self.bragg_angles(interplanar_spacings)
        cell_intensity = self.normalized_cell_intensity(structure, bragg_angles)
        positions = self.get_positions(structure, points)
        
        #just realized that the "canonical" lens aberration radius formula doesn't depend on the plane examined. weird.
        #TODO: look into lens aberration formula
        for plane in cell_intensity.keys():
            position = positions[plane]
            hkl = plane
            intensity = cell_intensity[plane]
            film_radius = 0.91 * (10**-3 * self.cs * self.wavelength_rel()**3)**Fraction('1/4')
            d_spacing = interplanar_spacings[plane]
            dot = TEMDot(position, hkl, intensity, film_radius, d_spacing)
            dots.append(dot)
        return dots    
    
    def show_plot_2d(self, structure):
        """
        Generates the 2D diffraction pattern of the input structure.
        Args:
            structure (Structure): The input structure.
        Returns:
            none (shows 2D DP)
        """
        points = self.generate_points(-10, 11)
        TEM_dots = self.TEM_dots(structure, points)
        film_radius = 0.91 * (10**-3 * self.cs * self.wavelength_rel()**3)**Fraction('1/4')

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
                x = xs,
                y = ys,
                text = hkls,
                hoverinfo = 'text',
                mode = 'markers',
                marker=dict(
                    size=8,
                    cmax=1,
                    cmin=0,
                    color=intensities,
                    colorbar=dict(
                        title='Colorbar',
                        yanchor = 'top'
                    ),
                    colorscale= [[0, 'black'], [1.0, 'white']]
                ),
                showlegend = False
            ), go.Scatter(
                x = [0],
                y = [0],
                text = "(0, 0, 0): Direct beam",
                hoverinfo = 'text',
                mode = 'markers',
                marker = dict(
                    size = 14,
                    cmax = 1,
                    cmin = 0,
                    color = 'white'
                    ),
            )
        ]
        layout = go.Layout(
            title = '2D Diffraction Pattern<br>Beam Direction: ' + ''.join(str(e) for e in self.beam_direction),
            font=dict(
                family='Comic Sans, monospace',
                size=18,
                color='#7f7f7f'),
            hovermode = 'closest',
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
            width = 600,
            height = 600,
            paper_bgcolor = 'rgba(100,110,110,0.5)',
            plot_bgcolor = 'black'
        )

        fig = go.Figure(data=data, layout=layout)
        poff.iplot(fig, filename='stuff')

    def get_pattern_2d(self, structure):
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
        