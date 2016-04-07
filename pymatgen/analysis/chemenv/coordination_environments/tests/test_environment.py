#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder, AbstractGeometry
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries, UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import StructureEnvironments
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from pymatgen.analysis.chemenv.utils.scripts_utils import visualize
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.io.cif import CifParser
from pymatgen.core.sites import PeriodicSite
from math import factorial

import numpy as np
import itertools
from random import shuffle


if __name__ == '__main__':

    cg_symbol = 'T:4'

    allcg = AllCoordinationGeometries()
    cg = allcg[cg_symbol]

    lgf = LocalGeometryFinder()
    lgf.setup_parameters(centering_type='centroid', structure_refinement=lgf.STRUCTURE_REFINEMENT_NONE,
                         include_central_site_in_centroid=False)


    myindices = range(cg.coordination_number)

    test = raw_input('Enter if you want to test all possible permutations ("all" or "a") or a given number of random permutations (i.e. "25")')

    if test == 'all' or test == 'a':
        perms_iterator = itertools.permutations(myindices)
        nperms = factorial(cg.coordination_number)
    else:
        try:
            nperms = int(test)
        except:
            raise ValueError('Could not turn {} into integer ...'.format(test))
        perms_iterator = []
        for ii in range(nperms):
            shuffle(myindices)
            perms_iterator.append(list(myindices))

    iperm = 1
    for indices_perm in perms_iterator:


        lgf.setup_test_perfect_environment(cg_symbol, indices=indices_perm)

        lgf.perfect_geometry = AbstractGeometry.from_cg(cg=cg)

        se = lgf.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.2)

        iperm += 1

    exit()


    strategy = SimplestChemenvStrategy()

    se = lgf.compute_structure_environments_detailed_voronoi(maximum_distance_factor=1.2)

    strategy.set_structure_environments(se)

    site = se.equivalent_sites[0][0]
    isite = se.structure.index(site)
    ces = strategy.get_site_coordination_environment(site=site)
    neighbors_map = strategy.structure_environments.voronoi.neighbors_map(isite, strategy.distance_cutoff,
                                                                          strategy.angle_cutoff,
                                                                          strategy.additional_condition)
    this_site_this_map_neighbors_list = (strategy.structure_environments.voronoi.neighbors_lists
                                             [isite]
                                             [neighbors_map['i_distfactor']]
                                             [neighbors_map['i_angfactor']]
                                             [neighbors_map['i_additional_condition']])

    print(len(this_site_this_map_neighbors_list))
    for nlist in this_site_this_map_neighbors_list:
        print(nlist)
    print(neighbors_map)
    print(strategy.structure_environments.voronoi.voronoi_list[isite][0])
    print(strategy.structure_environments.voronoi.voronoi_list[isite][0][0])
    print(strategy.structure_environments.voronoi.voronoi_list[isite][0][1])

    cn_map = (strategy.structure_environments.voronoi.parameters_to_unique_coordinated_neighbors_map
                  [isite]
                  [neighbors_map['i_distfactor']][neighbors_map['i_angfactor']]
                  [neighbors_map['i_additional_condition']])

    print(cn_map)


    print(ces)

    if True:
        vis = StructureVis(show_polyhedron=False, show_unit_cell=True)
        vis.show_help = False
        vis.set_structure(se.structure)
        for eqslist in se.equivalent_sites:
            site = eqslist[0]
            isite = se.structure.index(site)
            try:
                if strategy.uniquely_determines_coordination_environments:
                    ces = strategy.get_site_coordination_environments(site=site)
                else:
                    ces = strategy.get_site_coordination_environments_fractions(site)
            except NeighborsNotComputedChemenvError:
                continue
            if ces is None:
                continue
            if len(ces) == 0:
                continue
            comp = site.species_and_occu
            #ce = strategy.get_site_coordination_environment(site)
            if strategy.uniquely_determines_coordination_environments:
                ce = ces[0]
                if ce is None:
                    continue
                thecg = allcg.get_geometry_from_mp_symbol(ce[0])
                mystring = 'Environment for site #{} {} ({}) : {} ({})\n'.format(str(isite),
                                                                                 comp.get_reduced_formula_and_factor()[0],
                                                                                 str(comp),
                                                                                 thecg.name,
                                                                                 ce[0])
            else:
                mystring = 'Environments for site #{} {} ({}) : \n'.format(str(isite),
                                                                           comp.get_reduced_formula_and_factor()[0],
                                                                           str(comp))
                for ce in ces:
                    cg = allcg.get_geometry_from_mp_symbol(ce[0])
                    mystring += ' - {} ({}): {:.2f} % (csm : {:2f})\n'.format(cg.name, cg.mp_symbol,
                                                                              100.0*ce[2],
                                                                              ce[1]['symmetry_measure'])
            print(mystring)
            mydeltas = []
            multiplicity = [1, 1, 1]
            for i0 in range(int(multiplicity[0])):
                for i1 in range(int(multiplicity[1])):
                    for i2 in range(int(multiplicity[2])):
                        mydeltas.append(np.array([1.0*i0, 1.0*i1, 1.0*i2], np.float))


        strategy.set_structure_environments(se)
        for isite, site in enumerate(se.structure):
            try:
                ces = strategy.get_site_coordination_environments(site)
            except NeighborsNotComputedChemenvError:
                continue
            if len(ces) == 0:
                continue
            #ce = ces[0]
            ce = strategy.get_site_coordination_environment(site)
            if ce is not None and ce[0] != UNCLEAR_ENVIRONMENT_SYMBOL:
                for mydelta in mydeltas:
                    psite = PeriodicSite(site._species, site._fcoords + mydelta, site._lattice,
                                         properties=site._properties)
                    vis.add_site(psite)
                    neighbors = strategy.get_site_neighbors(psite)
                    print(ce[1])
                    draw_cg(vis, psite, neighbors, cg=lgf.cg.get_geometry_from_mp_symbol(ce[0]),
                            # perfect2local_map=ce[1]['perfect2local_map'])
                            perm=ce[1]['permutation'])
                    # neighbors = strategy.get_site_neighbors(site)
                    # draw_cg(vis, site, neighbors, cg=lgf.cg.get_geometry_from_mp_symbol(ce[0]),
                    #        perm=ce[1]['permutation'])
        vis.show()

