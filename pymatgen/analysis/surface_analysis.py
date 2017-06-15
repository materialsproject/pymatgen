# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
import numpy as np

"""
This module provides classes to perform analysis slab models.
"""

__author__ = "Richard Tran"
__copyright__ = "Copyright 2017, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__status__ = "Production"
__date__ = "April 17, 2017"


EV_PER_ANG2_TO_JOULES_PER_M2 = 16.0217656


class SurfaceAnalysis(object):
    def __init__(self, slab, top=True, tol=2):

        self.ucell = slab.oriented_unit_cell
        self.tol = tol
        self.top = top
        self.slab = slab

    def bulk_coordination(self):

        cn_bulk = []
        ucell_coord_finder = VoronoiCoordFinder(self.ucell)
        for i, site in enumerate(self.ucell):
            cn = np.around(ucell_coord_finder.get_coordination_number(i), decimals=self.tol)
            if cn not in cn_bulk:
                cn_bulk.append(cn)

        return cn_bulk

    def surface_coordination(self):

        cn_bulk = self.bulk_coordination()

        cn_slab = []
        slab_coord_finder = VoronoiCoordFinder(self.slab, allow_pathological=True)
        for i, site in enumerate(self.slab):
            if self.top and site.frac_coords[2] > 0.5:
                cn = np.around(slab_coord_finder.get_coordination_number(i), decimals=self.tol)
                if cn not in cn_bulk:
                    cn_slab.append([cn, site])
            elif not self.top and site.frac_coords[2] < 0.5:
                cn = np.around(slab_coord_finder.get_coordination_number(i), decimals=self.tol)
                if cn not in cn_bulk:
                    cn_slab.append([cn, site])

        return cn_slab

    def surface_coordination_nn(self):

        cn_bulk = self.bulk_coordination()
        blength = get_bond_length(self.ucell)

        cn_slab = []
        for i, site in enumerate(self.slab):

            if self.top and site.frac_coords[2] > 0.5:
                neighbors = self.slab.get_neighbors(site, blength + 0.1)
                cn = len(neighbors)
                if cn not in cn_bulk:
                    cn_slab.append([cn, site, neighbors])

            elif not self.top and site.frac_coords[2] < 0.5:
                neighbors = self.slab.get_neighbors(site, blength + 0.1)
                cn = len(neighbors)
                if cn not in cn_bulk:
                    cn_slab.append([cn, site, neighbors])

        return cn_slab

    def roughness(self, from_nn=True):

        # for now, only handles one coordination number in a bulk
        if from_nn:
            cn_slab = self.surface_coordination_nn()
        else:
            cn_slab = self.surface_coordination()
        bulk_cn = self.bulk_coordination()[0]
        broken_bonds = [bulk_cn - cn[0] for cn in cn_slab]

        return sum(broken_bonds) / (2 * self.slab.surface_area)

    def bb_model_surface_energy(self):
        #         cn_slab = self.surface_coordination_nn()
        return self.roughness() / (self.bulk_coordination())[0]


def bb_manual(slab, bonds, bondlength=3):

    # check if slab is centered
    slab_pos = np.mean([site.frac_coords[2] for site in slab])
    if 0.55 < slab_pos or 0.45 > slab_pos:
        slab = center_slab(slab)

    broken_bonds = {}
    nbb = 0
    for site in slab:
        nn = slab.get_neighbors(site, bondlength, include_index=True)
        if len(nn) != cn and site.frac_coords[2] > 0.5:
            nbb += cn - len(nn)
    return nbb

def center_slab(slab):
    indices = [i for i, site in enumerate(slab)]
    ave_c = np.mean([site.frac_coords[2] for site in slab])
    slab.translate_sites(indices, [0,0,abs(0.5-ave_c)])
    return slab


def bb_manual(slab, bonds, top=True):
    # check if slab is centered
    slab_pos = np.mean([site.frac_coords[2] for site in slab])
    if 0.55 < slab_pos or 0.45 > slab_pos:
        slab = center_slab(slab)

    broken_bonds = {}
    for bond in bonds.keys():

        # First we check the cn of the bulk for each type of bond
        center_ion = bond[0]
        mean_cn = []
        ucell = slab.oriented_unit_cell
        for site in ucell:
            cn = 0
            if str(site.specie) == center_ion:
                nn = ucell.get_neighbors(site, bonds[bond],
                                         include_index=True)

                for n in nn:

                    if str(n[0].specie) == bond[1]:
                        cn += 1
                mean_cn.append(cn)
        cn = int(np.mean(mean_cn))

        # Next we use the cn of the bulk as
        # reference to find the number of broken bonds
        nbb = 0
        for site in slab:
            if str(site.specie) == center_ion:
                nn = slab.get_neighbors(site, bonds[bond],
                                        include_index=True)

                def count_nbb(nbb):
                    slab_cn = 0
                    for n in nn:
                        if str(n[0].specie) == bond[1]:
                            slab_cn += 1
                    nbb += cn - slab_cn
                    return nbb

                if top and site.frac_coords[2] > 0.5:
                    nbb = count_nbb(nbb)
                if not top and site.frac_coords[2] < 0.5:
                    nbb = count_nbb(nbb)

        broken_bonds[bond] = nbb / slab.surface_area

    return broken_bonds


def center_slab(slab):
    indices = [i for i, site in enumerate(slab)]
    ave_c = np.mean([site.frac_coords[2] for site in slab])
    slab.translate_sites(indices, [0, 0, abs(0.5 - ave_c)])
    return slab