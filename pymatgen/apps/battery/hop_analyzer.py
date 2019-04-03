# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
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
from pymatgen.io.vasp import Chgcar
import operator
import networkx as nx
import numpy as np
import pandas as pd
from pandas.io.json import json_normalize
__author__ = "Jimmy Shen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"
__date__ = "April 1, 2019"

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

grouper_sm = StructureMatcher()

def generic_groupby(list_in, comp=operator.eq, lab_num=True):
    '''
    Group a list of unhasable objects
    Return a list of labels that represent which group the entry is in
    '''
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
    Local environment class to help look for migration paths through the lattice

    Since we are interested in many migration pathways through the material, we need a more connections than other local_env classes offer

    Args:
        cutoff (float): cutoff radius in Angstrom to look for trial
            near-neighbor sites (default: 10.0).
    """

    def __init__(self, cutoff=5.0):
        self.cutoff = cutoff

    def get_nn_info(self, structure, n):
        """
        Get all near-neighbor sites as well as the associated image locations
        and weights of the site with index n using the closest neighbor
        distance-based method.

        Args:
            structure (Structure): input structure.
            n (integer): index of site for which to determine near
                neighbors.

        Returns:
            siw (list of tuples (Site, array, float)): tuples, each one
                of which represents a neighbor site, its image location,
                and its weight.
        """

        site = structure[n]
        neighs_dists = structure.get_neighbors(site, self.cutoff)
        min_dist = min([dist for neigh, dist in neighs_dists])

        siw = []
        for s, dist in neighs_dists:
            w = min_dist / dist
            siw.append({'site': s,
                        'image': self._get_image(structure, s),
                        'weight': w,
                        'site_index': self._get_original_site(structure, s)})
        return siw


class MigrationPathAnalyzer():
    """
    Methods for analyzing the migration path of a material using the following scheme:
    - Map the relaxed sites of a material back to the empty host lattice
    - Apply symmetry operations of the empty lattice to obtain the other positions of the intercollated atom
    - Get the symmetry inequivalent hops
    - Get the migration barriers for each inequivalent hop

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

        :param base_entry: the structure without a working ion for us to analyze the migration
        :param base_aeccar: Chgcar object that contains the AECCAR0 + AECCAR2
        :param single_cat_entries: list of structures containing a single cation at different positions
        :param cation: a String symbol or Element for the cation. It must be positively charged, but can be 1+/2+/3+ etc.
        :param ltol: parameter for StructureMatcher
        :param stol: parameter for StructureMa
        :param angle_tol: parameter for StructureMa

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
        self.base_structure_full_sites.sites.extend(self.base_entry.structure.sites)

    def match_ent_to_base(self, ent):
        """
        Transform the structure of one entry to match the base structure

        :param ent: inserted structure with cation atoms
        :returns: entry with modified structure
        :rtype: ComputedStructureEntry

        """
        new_ent = deepcopy(ent)
        new_struct = self.sm.get_s2_like_s1(self.base_entry.structure,
                                            ent.structure)
        new_ent.structure = new_struct
        return new_ent

    def get_all_sym_sites(self, ent):
        """
        Return all of the symmetry equivalent sites

        :param ent: ComputedStructureEntry that contains cation
        :returns: Structure containing all of the symmetry equivalent sites
        :rtype: Structure

        """

        sa = SpacegroupAnalyzer(
            self.base_entry.structure, symprec=0.3, angle_tolerance=10)
        host_allsites = self.base_entry.structure.copy()
        host_allsites.remove_species(host_allsites.species)
        pos_Li = list(
            filter(lambda isite: isite.species_string == 'Li',
                   ent.structure.sites))

        for isite in pos_Li:
            host_allsites.insert(0, 'Li', isite.frac_coords, properties={'inserted_energy' : ent.energy})

        for op in sa.get_space_group_operations():
            struct_tmp = host_allsites.copy()
            struct_tmp.apply_operation(symmop=op, fractional=True)
            for isite in struct_tmp.sites:
                if isite.species_string == "Li":
                    host_allsites.insert(0, 'Li', isite.frac_coords, properties={'inserted_energy' : ent.energy})

        host_allsites.merge_sites(mode='average') # keeps only one position but average the properties

        return host_allsites

    def get_full_sites(self):
        """
        Get each group of symmetry inequivalent sites and combine them

        :returns: Structure with all possible Li sites, the enregy of the structure is stored as a site property
        :rtype: Structure

        """
        res = []
        for itr in self.translated_single_cat_entries:
            res.extend(self.get_all_sym_sites(itr).sites)
        res = Structure.from_sites(res)
        res.merge_sites(tol=1.0, mode='average')
        return res

    def get_graph(self):
        # Generate the graph edges between these sites
        self.gt = StructureGraph.with_local_env_strategy(self.full_sites, ConnectSitesNN())
        self.gt.set_node_attributes()

    def compare_edges(self, edge1, edge2):
        # Test
        temp_struct1 = self.base_entry.structure.copy()
        temp_struct2 = self.base_entry.structure.copy()

        temp_struct1.insert(0, self.cation, self._edgelist.iloc[edge1]['i_pos'])
        temp_struct1.insert(0, self.cation, self._edgelist.iloc[edge1]['f_pos'])
        temp_struct2.insert(0, self.cation, self._edgelist.iloc[edge2]['i_pos'])
        temp_struct2.insert(0, self.cation, self._edgelist.iloc[edge2]['f_pos'])
        return grouper_sm.fit(temp_struct1, temp_struct2)

    def get_edges_labels(self, mask_file=None):

        d = [{"isite" : u, "fsite" : v, "to_jimage" : d['to_jimage'], 'edge_tuple' : (u, v)} for u, v, d in self.gt.graph.edges(data=True)]
        self._edgelist = pd.DataFrame(d)

        self._edgelist['i_pos'] = self._edgelist.apply(lambda u : self.full_sites.sites[u.isite].frac_coords, axis=1)
        self._edgelist['f_pos'] = self._edgelist.apply(lambda u : self.full_sites.sites[u.fsite].frac_coords + u.to_jimage, axis=1)

        edge_lab = generic_groupby(self._edgelist.index.values, comp = self.compare_edges)
        self._edgelist.loc[:, 'edge_label'] = edge_lab

        # write the image
        self.unique_edges = self._edgelist.drop_duplicates('edge_label', keep='first').copy()

        # set up the grid
        aa = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(0)),
                         endpoint=False)
        bb = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(1)),
                         endpoint=False)
        cc = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(2)),
                         endpoint=False)
        AA, BB, CC = np.meshgrid(aa, bb, cc, indexing='ij')
        fcoords = np.vstack([AA.flatten(), BB.flatten(), CC.flatten()]).T

        IMA, IMB, IMC = np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing='ij')
        images = np.vstack([IMA.flatten(), IMB.flatten(), IMC.flatten()]).T

        # get the charge density masks for each hop (for plotting and sanity check purposes)
        idx_pbc_mask = np.zeros_like(AA)
        surf_idx=0
        total_chg=[]
        if mask_file:
            mask_out = copy(self.base_aeccar)
            mask_out.data['total'] = np.zeros_like(AA)

        for _, row in self.unique_edges.iterrows():
            pbc_mask = np.zeros_like(AA).flatten()
            e0 = row.i_pos.astype('float64')
            e1 = row.f_pos.astype('float64')

            cart_e0 = np.dot(e0, self.base_aeccar.structure.lattice.matrix)
            cart_e1 = np.dot(e1, self.base_aeccar.structure.lattice.matrix)
            pbc_mask = np.zeros_like(AA,dtype=bool).flatten()
            for img in images:
                grid_pos = np.dot(fcoords + img, self.base_aeccar.structure.lattice.matrix)
                proj_on_line = np.dot(grid_pos - cart_e0, cart_e1 - cart_e0) / (np.linalg.norm(cart_e1 - cart_e0))
                dist_to_line = np.linalg.norm(
                    np.cross(grid_pos - cart_e0, cart_e1 - cart_e0) / (np.linalg.norm(cart_e1 - cart_e0)), axis=-1)

                mask = (proj_on_line >= 0) * (proj_on_line < np.linalg.norm(cart_e1 - cart_e0)) * (dist_to_line < 0.5)
                pbc_mask = pbc_mask + mask
            pbc_mask = pbc_mask.reshape(AA.shape)
            if mask_file:
                mask_out.data['total'] = pbc_mask
                mask_out.write_file('{}_{}.vasp'.format(mask_file,row['edge_tuple']))


            total_chg.append(self.base_aeccar.data['total'][pbc_mask].sum() / self.base_aeccar.ngridpts / self.base_aeccar.structure.volume)

        self.complete_mask=idx_pbc_mask
        print('test')
        self.unique_edges.loc[self.unique_edges.index, 'chg_total'] = total_chg

    def _setup_grids(self):
        """
        Populate the internal varialbes used for defining the grid points in the charge density analysis
        """
        # set up the grid
        aa = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(0)),
                         endpoint=False)
        bb = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(1)),
                         endpoint=False)
        cc = np.linspace(0, 1, len(self.base_aeccar.get_axis_grid(2)),
                         endpoint=False)
        AA, BB, CC = np.meshgrid(aa, bb, cc, indexing='ij')

        IMA, IMB, IMC = np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing='ij')

        # store these
        self.uc_grid_shape = AA.shape
        self.fcoords = np.vstack([AA.flatten(), BB.flatten(), CC.flatten()]).T
        self.images = np.vstack([IMA.flatten(), IMB.flatten(), IMC.flatten()]).T

    def _get_chg_between_sites_tube(self, edge_index, mask_file_seedname=None):
        """Calculate the amount of charge that a cation has to move through in order to complete a hop

        :param edge_index: the index value to read from self._edgelist
        :param mask_file_seedname: seedname for output of the migration path masks (for debugging and visualization)
        :returns: total amount of chg in a tube connecting the endpoints
        :rtype: float

        # TODO use the cell centers as the positions for evaluating the tubes

        """
        try:
            self._tube_radius
        except NameError:
            logger.error("The radius of the tubes for charge analysis need to be first.")

        pbc_mask = np.zeros(self.uc_grid_shape).flatten()
        e0 = self._edgelist.iloc[edge_index].i_pos.astype('float64')
        e1 = self._edgelist.iloc[edge_index].f_pos.astype('float64')

        cart_e0 = np.dot(e0, self.base_aeccar.structure.lattice.matrix)
        cart_e1 = np.dot(e1, self.base_aeccar.structure.lattice.matrix)
        pbc_mask = np.zeros(self.uc_grid_shape, dtype=bool).flatten()
        for img in self.images:
            grid_pos = np.dot(self.fcoords + img, self.base_aeccar.structure.lattice.matrix)
            proj_on_line = np.dot(grid_pos - cart_e0, cart_e1 - cart_e0) / (np.linalg.norm(cart_e1 - cart_e0))
            dist_to_line = np.linalg.norm(
                np.cross(grid_pos - cart_e0, cart_e1 - cart_e0) / (np.linalg.norm(cart_e1 - cart_e0)), axis=-1)

            mask = (proj_on_line >= 0) * (proj_on_line < np.linalg.norm(cart_e1 - cart_e0)) * (dist_to_line < self._tube_radius)
            pbc_mask = pbc_mask + mask
        pbc_mask = pbc_mask.reshape(self.uc_grid_shape)

        if mask_file_seedname:
            mask_out.data['total'] = pbc_mask
            mask_out.write_file('{}_{}.vasp'.format(mask_file,self._edgelist.iloc[edge_index].edge_tuple))

        return self.base_aeccar.data['total'][pbc_mask].sum() / self.base_aeccar.ngridpts / self.base_aeccar.structure.volume

    def _populate_unique_edge_df(self):
        """
        Populate the self.unique_edges dataframe which is a subset of self._edgelist that are unique
        """

        if not hasattr(self, 'gt'):
            self.get_graph()
        d = [{"isite" : u, "fsite" : v, "to_jimage" : d['to_jimage'], 'edge_tuple' : (u, v)} for u, v, d in self.gt.graph.edges(data=True)]
        self._edgelist = pd.DataFrame(d)

        self._edgelist['i_pos'] = self._edgelist.apply(lambda u : self.full_sites.sites[u.isite].frac_coords, axis=1)
        self._edgelist['f_pos'] = self._edgelist.apply(lambda u : self.full_sites.sites[u.fsite].frac_coords + u.to_jimage, axis=1)

        edge_lab = generic_groupby(self._edgelist.index.values, comp = self.compare_edges)
        self._edgelist.loc[:, 'edge_label'] = edge_lab

        # write the image
        self.unique_edges = self._edgelist.drop_duplicates('edge_label', keep='first').copy()

    def _get_chg_values_for_unique_hops(self, tube_radius=0.5):
        """Populate the charge density values for each 

        :param tube_radius: The radius (in Angstroms) of the tube used to evaluate the total charge density for one hop

        """
        self._tube_radius = tube_radius
        self._setup_grids()
        total_chg = self.unique_edges.apply(lambda row : self._get_chg_between_sites_tube(row.name), axis=1)
        self.unique_edges.loc[:, 'chg_total_tube'] = total_chg

    def _get_chg_values_for_all_hops(self):
        """
        Populate the charge density values for each hop in the full list using the unque hops as a lookup table
        """
        def read_value_from_uniq(row):
            select_uniq = self.unique_edges.edge_label==row.edge_label
            return self.unique_edges.loc[select_uniq, 'chg_total_tube'].values[0]

        total_chg = self._edgelist.apply(read_value_from_uniq, axis=1)
        self._edgelist['chg_total_tube'] = total_chg

if __name__ == "__main__":
    main()
