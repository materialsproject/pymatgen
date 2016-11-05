# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
from __future__ import absolute_import, print_function

"""
This module provides classes used to enumerate surface sites
and to find adsorption sites on slabs
"""

import numpy as np
from six.moves import range
from pymatgen import Structure, Lattice
import tempfile
import sys
import subprocess
import itertools
from pyhull.delaunay import DelaunayTri
from pyhull.voronoi import VoronoiTess
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.analyzer import generate_full_symmops
from pymatgen.util.coord_utils import in_coord_list, in_coord_list_pbc
from pymatgen.core.sites import PeriodicSite
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from pymatgen.core.surface import generate_all_slabs

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "December 2, 2015"

class AdsorbateSiteFinder(object):
    """
    This class finds adsorbate sites on slabs and generates
    adsorbate structures
    """

    def __init__(self, slab, selective_dynamics=False, 
                 height_threshold=None):
        """
        Create an AdsorbateSiteFinder object.

        Args:
            slab (Slab): slab object for which to find adsorbate
            sites
        """
        self.mi_string = ''.join([str(i) for i in slab.miller_index])
        # get surface normal from miller index
        self.mvec = mi_vec(slab.miller_index)
        slab = self.assign_site_properties(slab)
        if selective_dynamics:
            slab = self.assign_selective_dynamics(slab)
        self.slab = slab

    @classmethod
    def from_bulk_and_miller(cls, structure, miller_index, min_slab_size=5.0,
                             min_vacuum_size=10.0, max_normal_search=None, 
                             center_slab = True, selective_dynamics=False):
        """
        This method constructs the adsorbate site finder 
        from a bulk structure and a miller index, which
        allows the surface sites to be determined from
        the difference in bulk and slab coordination
        
        Args:
            structure (Structure): structure from which slab
                input to the ASF is constructed
            miller_index (3-tuple or list): index for slab
        """
        # TODO: for some reason this is buggy with primitive cells, 
        # might be a problem elsewhere
        vcf_bulk = VoronoiCoordFinder(structure)
        bulk_coords = [len(vcf_bulk.get_coordinated_sites(n))
                       for n in range(len(structure))]
        struct = structure.copy(site_properties = {'bulk_coordinations':bulk_coords})
        # This should probably be a slab generation method
        #   rather than selecting from generate_all_slabs
        slabs = generate_all_slabs(struct, max_index=max(miller_index), 
                                   min_slab_size=min_slab_size,
                                   min_vacuum_size=min_vacuum_size,
                                   max_normal_search = max_normal_search,
                                   center_slab = center_slab)

        slab_dict = {slab.miller_index:slab for slab in slabs}
        
        if miller_index not in slab_dict:
            raise ValueError("Miller index not in slab dict")

        this_slab = slab_dict[miller_index]

        vcf_surface = VoronoiCoordFinder(this_slab)
        surf_props = []
        this_mi_vec = mi_vec(this_slab.miller_index)
        mi_mags = [np.dot(this_mi_vec, site.coords) for site in this_slab]
        average_mi_mag = np.average(mi_mags)
        for n, site in enumerate(this_slab):
            bulk_coord = this_slab.site_properties['bulk_coordinations'][n]
            surf_coord = len(vcf_surface.get_coordinated_sites(n))
            mi_mag = np.dot(this_mi_vec, site.coords)
            if surf_coord != bulk_coord and mi_mags[n] > average_mi_mag:
                surf_props += ['surface']
            else:
                surf_props += ['subsurface']
        new_site_properties = {'surface_properties':surf_props}
        new_slab = this_slab.copy(site_properties=new_site_properties)
        return cls(new_slab, selective_dynamics)

    def find_surface_sites_by_height(self, slab, window = 0.9):
        """
        This method finds surface sites by determining which sites are within
        a threshold value in height from the topmost site in a list of sites

        Args:
            site_list (list): list of sites from which to select surface sites
            window (float): threshold in angstroms of distance from topmost
                site in slab along the slab c-vector to include in surface 
                site determination

        Returns:
            list of sites selected to be within a threshold of the highest
        """

        # Get projection of coordinates along the miller index
        m_projs = np.array([np.dot(site.coords, self.mvec)
                            for site in slab.sites])
        # Mask based on window threshold along the miller index
        mask = (m_projs - np.amax(m_projs)) >= -window
        return [slab.sites[n] for n in np.where(mask)[0]]

    def find_surface_sites_by_coordination(self, slab):
        """
        This method finds surface sites by testing the coordination against the
        bulk coordination
        """
        if 'bulk_coordination' not in slab.site_properties:
            raise ValueError("Input slabs must have bulk_coordination assigned."
                             "Use adsorption.generate_decorated_slabs to assign.")
        pass

    def assign_site_properties(self, slab, height=0.1):
        """
        Assigns site properties.
        """
        if 'surface_properties' in slab.site_properties.keys():
            return slab
        else:
            surf_sites = self.find_surface_sites_by_height(slab)
        surf_props = ['surface' if site in surf_sites
                      else 'subsurface' for site in slab.sites]
        return slab.copy(
            site_properties = {'surface_properties': surf_props})

    def get_extended_surface_mesh(self, radius = 6.0, window = 1.0):
        """
        """
        surf_str = Structure.from_sites(self.surface_sites)
        surface_mesh = []
        for site in surf_str.sites:
            surface_mesh += [site]
            surface_mesh += [s[0] for s in surf_str.get_neighbors(site,
                                                                  radius)]
        return list(set(surface_mesh))

    @property
    def surface_sites(self):
        """
        convenience method to return a list of surface sites
        """
        return [site for site in self.slab.sites 
                if site.properties['surface_properties']=='surface']

    def find_adsorption_sites(self, distance = 2.0, put_inside = True,
                              symm_reduce = True, near_reduce = True,
                              positions = ['ontop', 'bridge', 'hollow'],
                              near_reduce_threshold = 1e-2, 
                              no_right_triangles = True):
        """
        """
        # find on-top sites
        ads_sites = []
        if 'ontop' in positions:
            ads_sites += [s.coords for s in self.surface_sites]
        # Get bridge sites via DelaunayTri of extended surface mesh
        mesh = self.get_extended_surface_mesh()
        sop = get_rot(self.slab)
        dt = DelaunayTri([sop.operate(m.coords)[:2] for m in mesh])
        for v in dt.vertices:
            # TODO: bridge/hollow should be refactored at some point
            #   to properly account for multi-coordinated sites
            if -1 not in v:
                vecs = []
                for data in itertools.combinations(v, 2):
                    # Add bridge sites at midpoints of edges of D. Tri
                    if 'bridge' in positions:
                        ads_sites += [self.ensemble_center(mesh, data, 
                                                           cartesian = True)]
                    vecs += [mesh[data[1]].coords - mesh[data[0]].coords]
                # Prevent addition of hollow sites in right triangles
                if no_right_triangles:
                    dots = np.array([np.dot(vec1, vec2) for vec1, vec2 
                                     in itertools.combinations(vecs, 2)])
                    if (np.abs(dots) < 1e-10).any():
                        continue
                # Add hollow sites at centers of D. Tri faces
                if 'hollow' in positions:
                    ads_sites += [self.ensemble_center(mesh, v, 

        if near_reduce:
            ads_sites = self.near_reduce(ads_sites, 
                                         threshold=near_reduce_threshold)
        
        if symm_reduce:
            ads_sites = self.symm_reduce(ads_sites)
        ads_sites = [ads_site + distance*self.mvec 
                     for ads_site in ads_sites]
        return ads_sites

    def symm_reduce(self, coords_set, threshold = 1e-6):
        """
        Reduces the set of adsorbate sites by finding removing
        symmetrically equivalent duplicates

        Args:
            coords_set: coordinate set in cartesian coordinates
            threshold: tolerance for distance equivalence, used
                as input to in_coord_list_pbc for dupl. checking
        """
        surf_sg = SpacegroupAnalyzer(self.slab, 0.1)
        symm_ops = surf_sg.get_symmetry_operations()
        unique_coords = []
        # Convert to fractional
        coords_set = [cart_to_frac(self.slab.lattice, coords) 
                      for coords in coords_set]
        for coords in coords_set:
            incoord = False
            for op in symm_ops:
                if in_coord_list_pbc(unique_coords, op.operate(coords), 
                                     atol = threshold):
                    incoord = True
                    break
            if not incoord:
                unique_coords += [coords]
        # convert back to cartesian
        return [frac_to_cart(self.slab.lattice, coords) 
                for coords in unique_coords]

    def near_reduce(self, coords_set, threshold = 1e-4):
        """
        Prunes coordinate set for coordinates that are within 
        threshold
        
        Args:
            coords_set (Nx3 array-like): list or array of coordinates
            threshold (float): threshold value for distance
        """
        unique_coords = []
        coords_set = [cart_to_frac(self.slab.lattice, coords) 
                      for coords in coords_set]
        for coord in coords_set:
            if not in_coord_list_pbc(unique_coords, coord, threshold):
                unique_coords += [coord]
        return [frac_to_cart(self.slab.lattice, coords) 
                for coords in unique_coords]

    def ensemble_center(self, site_list, indices, cartesian = True):
        """
        Finds the center of an ensemble of sites selected from
        a list of sites.  Helper method for the find_adsorption_sites
        algorithm.
        """
        if cartesian:
            return np.average([site_list[i].coords for i in indices], 
                              axis = 0)
        else:
            return np.average([site_list[i].frac_coords for i in indices], 
                              axis = 0)

    def add_adsorbate(self, molecule, ads_coord, repeat = None):
        """
        Adds an adsorbate at a particular coordinate
        """
        struct = self.slab.copy()
        if repeat:
            struct.make_supercell(repeat)
        if 'surface_properties' in struct.site_properties.keys():
            molecule.add_site_property("surface_properties",
                                       ["adsorbate"] * molecule.num_sites)
        if 'selective_dynamics' in struct.site_properties.keys():
            molecule.add_site_property("selective_dynamics",
                                       [[True, True, True]] * molecule.num_sites)
        for site in molecule:
            struct.append(site.specie, ads_coord + site.coords, coords_are_cartesian = True,
                          properties = site.properties)
        return struct
        
    def assign_selective_dynamics(self, slab):
        sd_list = []
        sd_list = [[False, False, False]  if site.properties['surface_properties']=='subsurface' 
                   else [True, True, True] for site in slab.sites]
        new_sp = slab.site_properties
        new_sp['selective_dynamics'] = sd_list
        return slab.copy(site_properties = new_sp)

    def generate_adsorption_structures(self, molecule, repeat = None, 
                                       min_lw = 5.0):
        if repeat is None:
            xrep = np.ceil(min_lw / np.linalg.norm(self.slab.lattice.matrix[0]))
            yrep = np.ceil(min_lw / np.linalg.norm(self.slab.lattice.matrix[1]))
            repeat = [xrep, yrep, 1]
        structs = []
        for coords in self.find_adsorption_sites():
            structs += [self.add_adsorbate(molecule, coords,
                                      repeat = repeat)]
        return structs

def mi_vec(mi_index):
    """
    Convenience function which returns the unit vector aligned 
    with the miller index.
    """
    mvec = np.array([1./n if n!=0 else 0 
                     for n in mi_index])
    return mvec / np.linalg.norm(mvec)

def get_rot(structure):
    """
    Gets the transformation to rotate the z axis into the miller index
    """
    a, b, c = structure.lattice.matrix
    new_x = a / np.linalg.norm(a)
    new_y = (b - np.dot(new_x, b) * new_x) / \
            np.linalg.norm(b - np.dot(new_x, b) * new_x)
    new_z = np.cross(new_x, new_y)
    if np.dot(new_z, c) < 0.:
        new_z = -new_z
    x, y, z = np.eye(3)
    rot_matrix = np.array([np.dot(*el) for el in 
                           itertools.product([x, y, z], 
                                   [new_x, new_y, new_z])]).reshape(3,3)
    rot_matrix = np.transpose(rot_matrix)
    sop = SymmOp.from_rotation_and_translation(rot_matrix)
    return sop

def reorient_z(structure):
    """
    reorients a structure such that the z axis is concurrent with the 
    normal to the A-B plane
    """
    struct = structure.copy()
    sop = get_rot(struct)
    struct.apply_operation(sop)
    return struct

def frac_to_cart(lattice, frac_coord):
    """
    converts fractional coordinates to cartesian
    """
    return np.dot(np.transpose(lattice.matrix), frac_coord)

def cart_to_frac(lattice, cart_coord):
    """
    converts cartesian coordinates to fractional
    """
    return np.dot(np.linalg.inv(np.transpose(lattice.matrix)), cart_coord)
