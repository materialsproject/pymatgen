# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Created on Jan 25, 2012
"""

__author__ = "Jimmy Shen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"
__date__ = "April 1, 2019"

from pprint import pprint
from collections import defaultdict
import math
from copy import deepcopy
import os
import logging
import sys
from monty.serialization import loadfn
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Composition, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import NearNeighbors
from pymatgen.io.vasp import Chgcar, VolumetricData
import operator
import networkx as nx
import numpy as np
import pandas as pd
from pandas.io.json import json_normalize

logger = logging.getLogger(__name__)

grouper_sm = StructureMatcher()


def generic_groupby(list_in, comp=operator.eq, lab_num=True):
    """
    Group a list of unhasable objects

    Args:
      list_in:
      comp: (Default value = operator.eq)
      lab_num: (Default value = True)

    Returns:

    """
    list_out = ['TODO'] * len(list_in)
    cnt = 0
    for i1, ls1 in enumerate(list_out):
        if ls1 != 'TODO':
            continue

        if not lab_num:
            list_out[i1] = alpha[cnt]
        else:
            list_out[i1] = cnt
        cnt += 1
        for i2, ls2 in enumerate(list_out):
            if comp(list_in[i1], list_in[i2]):
                list_out[i2] = list_out[i1]
    return list_out


class ConnectSitesNN(NearNeighbors):
    """
    Local environment class to help look for migration paths through the lattice.
    This allows us to use the StructureGraph.with_local_env_strategy function without any modifications
    Since we are interested in many migration pathways through the material,
    we need a more connections than other local_env classes offer.

    Args:
      cutoff(float): cutoff radius in Angstrom to look for trial
    near-neighbor sites (default: 10.0).

    Returns:

    """

    def __init__(self, cutoff=5.0):
        """


        Args:
          cutoff:  (Default value = 5.0)

        Returns:

        """
        self.cutoff = cutoff

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest neighbor
        distance-based method.

        Args:
          structure(Structure): input structure.
          n(integer): index of site for which to determine near
        neighbors.

        Returns:
          list of tuples (Site: tuples, each one
          list of tuples (Site: tuples, each one
          of which represents a neighbor site, its image location,
          list of tuples (Site: tuples, each one
          of which represents a neighbor site, its image location,
          and its weight.

        """

        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        min_dist = min([dist for neigh, dist in neighs_dists])

        siw = []
        for s, dist in neighs_dists:
            w = min_dist / dist
            siw.append({
                'site': s,
                'image': self._get_image(structure, s),
                'weight': w,
                'site_index': self._get_original_site(structure, s)
            })
        return siw


class MigrationPathAnalyzer():
    """
    Methods for analyzing the migration path of a material using the following scheme:
    - Map the relaxed sites of a material back to the empty host lattice
    - Apply symmetry operations of the empty lattice to obtain the other positions of the intercollated atom
    - Get the symmetry inequivalent hops
    - Get the migration barriers for each inequivalent hop

    Args:

    Returns:

    """

    def __init__(self,
                 base_entry,
                 single_cat_entries,
                 base_aeccar=None,
                 cation='Li',
                 ltol=0.2,
                 stol=0.3,
                 angle_tol=5):
        """
        Pass in a entries for analysis

        Args:
          base_entry: the structure without a working ion for us to analyze the migration
          single_cat_entries: list of structures containing a single cation at different positions
          base_aeccar: Chgcar object that contains the AECCAR0 + AECCAR2 (Default value = None)
          cation: a String symbol or Element for the cation. (Default value = 'Li')
          ltol: parameter for StructureMatcher (Default value = 0.2)
          stol: parameter for StructureMatcher (Default value = 0.3)
          angle_tol: parameter for StructureMatcher (Default value = 5)

        Returns:

        """

        self.sm = StructureMatcher(
            comparator=ElementComparator(),
            primitive_cell=False,
            ignored_species=[cation],
            ltol=ltol,
            stol=stol,
            angle_tol=angle_tol)

        self.single_cat_entries = single_cat_entries
        self.cation = cation
        self.base_entry = base_entry
        self.base_aeccar = base_aeccar

        logger.debug('See if the structures all match')
        for ent in self.single_cat_entries:
            assert (self.sm.fit(self.base_entry.structure, ent.structure))

        self.translated_single_cat_entries = list(
            map(self.match_ent_to_base, self.single_cat_entries))

        self.full_sites = self.get_full_sites()
        self.base_structure_full_sites = self.full_sites.copy()
        self.base_structure_full_sites.sites.extend(
            self.base_entry.structure.sites)

    def match_ent_to_base(self, ent):
        """
        Transform the structure of one entry to match the base structure

        Args:
          ent:

        Returns:
          ComputedStructureEntry: entry with modified structure

        """
        new_ent = deepcopy(ent)
        new_struct = self.sm.get_s2_like_s1(self.base_entry.structure,
                                            ent.structure)
        new_ent.structure = new_struct
        return new_ent

    def get_all_sym_sites(self, ent):
        """
        Return all of the symmetry equivalent sites

        Args:
          ent: ComputedStructureEntry that contains cation

        Returns:
          Structure: Structure containing all of the symmetry equivalent sites

        """

        sa = SpacegroupAnalyzer(
            self.base_entry.structure, symprec=0.3, angle_tolerance=10)
        host_allsites = self.base_entry.structure.copy()
        host_allsites.remove_species(host_allsites.species)
        pos_Li = list(
            filter(lambda isite: isite.species_string == 'Li',
                   ent.structure.sites))

        for isite in pos_Li:
            host_allsites.insert(
                0,
                'Li',
                isite.frac_coords,
                properties={'inserted_energy': ent.energy})

        for op in sa.get_space_group_operations():
            struct_tmp = host_allsites.copy()
            struct_tmp.apply_operation(symmop=op, fractional=True)
            for isite in struct_tmp.sites:
                if isite.species_string == "Li":
                    host_allsites.insert(
                        0,
                        'Li',
                        isite.frac_coords,
                        properties={'inserted_energy': ent.energy})

        host_allsites.merge_sites(
            mode='average'
        )  # keeps only one position but average the properties

        return host_allsites

    def get_full_sites(self):
        """
        Get each group of symmetry inequivalent sites and combine them

        Args:

        Returns:
          Structure: Structure with all possible Li sites, the enregy of the structure is stored as a site property

        """
        res = []
        for itr in self.translated_single_cat_entries:
            res.extend(self.get_all_sym_sites(itr).sites)
        res = Structure.from_sites(res)
        res.merge_sites(tol=1.0, mode='average')
        return res

    def _get_graph(self):
        """Construct the self.s_graph object"""
        # Generate the graph edges between these sites
        self.s_graph = StructureGraph.with_local_env_strategy(
            self.full_sites, ConnectSitesNN())
        self.s_graph.set_node_attributes()

    def compare_edges(self, edge1, edge2):
        """
        Given two idexs in the self.edgelist dataframe, say if they are the same object

        Args:
          edge1(int): index of one row in the self.edgelist dataframe
          edge2(int): index of one row in the self.edgelist dataframe

        Returns:
          bool: if the two structures, each contains two cations, match

        """
        # Test
        temp_struct1 = self.base_entry.structure.copy()
        temp_struct2 = self.base_entry.structure.copy()

        temp_struct1.insert(0, self.cation, self.edgelist.iloc[edge1]['i_pos'])
        temp_struct1.insert(0, self.cation, self.edgelist.iloc[edge1]['f_pos'])
        temp_struct2.insert(0, self.cation, self.edgelist.iloc[edge2]['i_pos'])
        temp_struct2.insert(0, self.cation, self.edgelist.iloc[edge2]['f_pos'])
        return grouper_sm.fit(temp_struct1, temp_struct2)

    def _setup_grids(self):
        """Populate the internal varialbes used for defining the grid points in the charge density analysis"""

        def _shift_grid(vv):
            """
            Move the grid points by half a step so that they sit in the center
            """
            step = vv[1] - vv[0]
            vv += step / 2.

        # set up the grid
        aa = np.linspace(
            0, 1, len(self.base_aeccar.get_axis_grid(0)), endpoint=False)
        bb = np.linspace(
            0, 1, len(self.base_aeccar.get_axis_grid(1)), endpoint=False)
        cc = np.linspace(
            0, 1, len(self.base_aeccar.get_axis_grid(2)), endpoint=False)
        # move the grid points to the center
        _shift_grid(aa)
        _shift_grid(bb)
        _shift_grid(cc)

        # mesh grid for each unit cell
        AA, BB, CC = np.meshgrid(aa, bb, cc, indexing='ij')

        # mesh grid for 3x3x3 set of nearby cells
        IMA, IMB, IMC = np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1],
                                    indexing='ij')

        # store these
        self._uc_grid_shape = AA.shape
        self._fcoords = np.vstack([AA.flatten(), BB.flatten(), CC.flatten()]).T
        self._images = np.vstack([IMA.flatten(),
                                  IMB.flatten(),
                                  IMC.flatten()]).T

    def _get_chg_between_sites_tube(self, edge_index, mask_file_seedname=None):
        """
        Calculate the amount of charge that a cation has to move through in order to complete a hop

        Args:
          edge_index: the index value to read from self.edgelist
          mask_file_seedname: seedname for output of the migration path masks (for debugging and visualization) (Default value = None)

        Returns:
          float: The total charge density in a tube that connects two sites of a given edges of the graph

        """
        try:
            self._tube_radius
        except NameError:
            logger.error(
                "The radius of the tubes for charge analysis need to be first."
            )

        pbc_mask = np.zeros(self._uc_grid_shape).flatten()
        e0 = self.edgelist.iloc[edge_index].i_pos.astype('float64')
        e1 = self.edgelist.iloc[edge_index].f_pos.astype('float64')

        cart_e0 = np.dot(e0, self.base_aeccar.structure.lattice.matrix)
        cart_e1 = np.dot(e1, self.base_aeccar.structure.lattice.matrix)
        pbc_mask = np.zeros(self._uc_grid_shape, dtype=bool).flatten()
        for img in self._images:
            grid_pos = np.dot(self._fcoords + img,
                              self.base_aeccar.structure.lattice.matrix)
            proj_on_line = np.dot(grid_pos - cart_e0, cart_e1 - cart_e0) / (
                np.linalg.norm(cart_e1 - cart_e0))
            dist_to_line = np.linalg.norm(
                np.cross(grid_pos - cart_e0, cart_e1 - cart_e0) /
                (np.linalg.norm(cart_e1 - cart_e0)),
                axis=-1)

            mask = (proj_on_line >= 0) * (proj_on_line < np.linalg.norm(
                cart_e1 - cart_e0)) * (dist_to_line < self._tube_radius)
            pbc_mask = pbc_mask + mask
        pbc_mask = pbc_mask.reshape(self._uc_grid_shape)

        if mask_file_seedname:
            mask_out = VolumetricData(
                structure=self.base_aeccar.structure.copy(),
                data={'total': self.base_aeccar.data['total']})
            sites_idx = self.edgelist.loc[edge_index, ['isite', 'fsite'
                                                       ]].values
            mask_out.structure.insert(
                0, "X",
                self.full_sites.sites[sites_idx[0]].frac_coords)
            mask_out.structure.insert(
                0, "X",
                self.full_sites.sites[sites_idx[1]].frac_coords)
            mask_out.data['total'] = pbc_mask
            mask_out.write_file('{}_{}.vasp'.format(
                mask_file_seedname, self.edgelist.iloc[edge_index].edge_tuple))

        return self.base_aeccar.data['total'][pbc_mask].sum(
        ) / self.base_aeccar.ngridpts / self.base_aeccar.structure.volume

    def _populate_unique_edge_df(self):
        """Populate the self.unique_edges dataframe which is a subset of self.edgelist that are unique"""

        if not hasattr(self, 'gt'):
            self._get_graph()
        d = [{
            "isite": u,
            "fsite": v,
            "to_jimage": d['to_jimage'],
            'edge_tuple': (u, v)
        } for u, v, d in self.s_graph.graph.edges(data=True)]
        self.edgelist = pd.DataFrame(d)

        self.edgelist['i_pos'] = self.edgelist.apply(
            lambda u: self.full_sites.sites[u.isite].frac_coords, axis=1)
        self.edgelist['f_pos'] = self.edgelist.apply(
            lambda u: self.full_sites.sites[u.fsite].frac_coords + u.to_jimage,
            axis=1)

        edge_lab = generic_groupby(
            self.edgelist.index.values, comp=self.compare_edges)
        self.edgelist.loc[:, 'edge_label'] = edge_lab

        # write the image
        self.unique_edges = self.edgelist.drop_duplicates(
            'edge_label', keep='first').copy()

    def get_chg_values_for_unique_hops(self,
                                       tube_radius=0.5,
                                       mask_file_seedname=None):
        """
        Populate the charge density values for each unique hop in the crystal

        Args:
          tube_radius: The radius (in Angstroms) of the tube used to evaluate the total charge density for one hop (Default value = 0.5)

        Returns:

        """
        if not hasattr(self, '_uc_grid_shape'):
            self._setup_grids()
        if not hasattr(self, 'tube_radius'):
            self._tube_radius = tube_radius

        total_chg = self.unique_edges.apply(
            lambda row: self._get_chg_between_sites_tube(
                row.name, mask_file_seedname=mask_file_seedname),
            axis=1)
        self.unique_edges.loc[:, 'chg_total_tube'] = total_chg

    def _read_value_from_uniq(self, row, key):
        """
        Get the value from the self.unique_edges edges dataframe

        Args:
          row: the row from another dataframe from which we read the edge_label
          key: The key of the quantity of interest

        Returns:
          : the value in the matching row of the self.unique_edges associated with the given key

        """
        select_uniq = self.unique_edges.edge_label == row.edge_label
        return self.unique_edges.loc[select_uniq, key].values[0]

    def get_chg_values_for_all_hops(self):
        """Populate the charge density values for each hop in the full list using the unque hops as a lookup table"""

        total_chg = self.edgelist.apply(
            lambda row: self._read_value_from_uniq(row, 'chg_total_tube'),
            axis=1)
        self.edgelist['chg_total_tube'] = total_chg
