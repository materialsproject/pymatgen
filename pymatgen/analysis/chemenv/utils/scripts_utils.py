# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module contains some utils for the main script of the chemenv package.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2012, The Materials Project"
__credits__ = "Geoffroy Hautier"
__version__ = "2.0"
__maintainer__ = "David Waroquiers"
__email__ = "david.waroquiers@gmail.com"
__date__ = "Feb 20, 2016"


from pymatgen import MPRester, Structure
from pymatgen.io.cif import CifParser
try:
    from pymatgen.vis.structure_vtk import StructureVis
    no_vis = False
except ImportError:
    StructureVis = None
    no_vis = True

try:
    input = raw_input
except NameError:
    pass

from pymatgen.core.sites import PeriodicSite
import pymongo
import re
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimpleAbundanceChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import TargettedPenaltiedAbundanceChemenvStrategy

from pymatgen.core.structure import Molecule
from collections import OrderedDict
import numpy as np


strategies_class_lookup = OrderedDict()
strategies_class_lookup['SimplestChemenvStrategy'] = SimplestChemenvStrategy
strategies_class_lookup['SimpleAbundanceChemenvStrategy'] = SimpleAbundanceChemenvStrategy
strategies_class_lookup['TargettedPenaltiedAbundanceChemenvStrategy'] = TargettedPenaltiedAbundanceChemenvStrategy


def draw_cg(vis, site, neighbors, cg=None, perm=None, perfect2local_map=None):
    if perm is not None and perfect2local_map is not None:
        raise ValueError('Only "perm" or "perfect2local_map" should be provided in draw_cg, not both')
    vis.add_bonds(neighbors, site)
    for n in neighbors:
        vis.add_site(n)
    if len(neighbors) < 3:
        vis.add_bonds(neighbors, site, color=[0.0, 1.0, 0.0], opacity=0.4, radius=0.175)
    else:
        if perm is not None:
            faces = cg.faces(neighbors, permutation=perm)
            edges = cg.edges(neighbors, permutation=perm)
        elif perfect2local_map is not None:
            faces = cg.faces(neighbors, perfect2local_map=perfect2local_map)
            edges = cg.edges(neighbors, perfect2local_map=perfect2local_map)
        else:
            faces = cg.faces(neighbors)
            edges = cg.edges(neighbors)
        symbol = list(site.species_and_occu.keys())[0].symbol
        mycolor = [float(i) / 255 for i in vis.el_color_mapping[symbol]]
        vis.add_faces(faces, mycolor, opacity=0.4)
        vis.add_edges(edges)




# Visualizing a coordination geometry
def visualize(cg, zoom=None, vis=None, myfactor=1.0):
    if vis is None:
        vis = StructureVis(show_polyhedron=False, show_unit_cell=False)
    myspecies = ["O"] * (cg.coordination_number+1)
    myspecies[0] = "Cu"
    coords = [np.zeros(3, np.float) + cg.central_site]

    for pp in cg.points:
        coords.append(np.array(pp) + cg.central_site)
    coords = [cc * myfactor for cc in coords]
    structure = Molecule(species=myspecies, coords=coords)
    vis.set_structure(structure=structure, reset_camera=True)
    # neighbors_list = coords[1:]
    draw_cg(vis, site=structure[0], neighbors=structure[1:], cg=cg)
    for ineighbor, neighbor in enumerate(structure[1:]):
        vis.add_text(neighbor.coords, '{}'.format(ineighbor), color=(0, 0, 0))
    if zoom is not None:
        vis.zoom(zoom)
    return vis


def welcome(chemenv_config):
    print('Chemical Environment package (ChemEnv)')
    print(chemenv_config.package_options_description())


def thankyou():
    print('Thank you for using the ChemEnv package')


def compute_environments(chemenv_configuration):
    string_sources = {'cif': {'string': 'a Cif file', 'regexp': '.*\.cif$'},
                      'mp': {'string': 'the Materials Project database', 'regexp': 'mp-[0-9]+$'}}
    questions = {'c': 'cif'}
    if chemenv_configuration.has_materials_project_access:
        questions['m'] = 'mp'
    lgf = LocalGeometryFinder()
    lgf.setup_parameters()
    allcg = AllCoordinationGeometries()
    strategy_class = strategies_class_lookup[chemenv_configuration.package_options['default_strategy']['strategy']]
    #TODO: Add the possibility to change the parameters and save them in the chemenv_configuration
    default_strategy = strategy_class()
    default_strategy.setup_options(chemenv_configuration.package_options['default_strategy']['strategy_options'])
    firsttime = True
    while True:
        if len(questions) > 1:
            found = False
            print('Enter the source from which the structure is coming or <q> to quit :')
            for key_character, qq in questions.items():
                print(' - <{}> for a structure from {}'.format(key_character, string_sources[qq]['string']))
            test = input(' ... ')
            if test == 'q':
                break
            if test not in list(questions.keys()):
                for key_character, qq in questions.items():
                    print(test)
                    print(test.__class__)
                    print(string_sources[qq]['regexp'].__class__)

                    if re.match(string_sources[qq]['regexp'], str(test)) is not None:
                        found = True
                        source_type = qq
                if not found:
                    print('Wrong key, try again ...')
                    continue
            else:
                source_type = questions[test]
        else:
            found = False
            source_type = list(questions.values())[0]
        if found and len(questions) > 1:
            input_source = test
        if source_type == 'cif':
            if not found:
                input_source = input('Enter path to cif file : ')
            cp = CifParser(input_source)
            structure = cp.get_structures()[0]
        elif source_type == 'mp':
            if not found:
                input_source = input('Enter materials project id (e.g. "mp-1902") : ')
            a = MPRester(chemenv_configuration.materials_project_api_key)
            structure = a.get_structure_by_material_id(input_source)
        lgf.setup_structure(structure)
        print('Computing environments for {} ... '.format(structure.composition.reduced_formula))
        se = lgf.compute_structure_environments_detailed_voronoi()
        print('Computing environments finished')
        while True:
            test = input('See list of environments determined for each (unequivalent) site ? '
                         '("y" or "n", "d" with details, "g" to see the grid) : ')
            strategy = default_strategy
            if test in ['y', 'd', 'g']:
                strategy.set_structure_environments(se)
                for eqslist in se.equivalent_sites:
                    site = eqslist[0]
                    isite = se.structure.index(site)
                    try:
                        if strategy.uniquely_determines_coordination_environments:
                            ces = strategy.get_site_coordination_environments(site)
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
                    if test in ['d', 'g'] and strategy.uniquely_determines_coordination_environments:
                        if thecg.mp_symbol != UNCLEAR_ENVIRONMENT_SYMBOL:
                            mystring += '  <Continuous symmetry measures>  '
                            mingeoms = se.ce_list[isite][thecg.coordination_number][0].minimum_geometries()
                            for mingeom in mingeoms:
                                mystring += '{} : {:.2f}       '.format(mingeom[0], mingeom[1]['symmetry_measure'])
                    print(mystring)
            if test == 'g':
                test = input('Enter index of site(s) for which you want to see the grid of parameters : ')
                indices = list(map(int, test.split()))
                print(indices)
                for isite in indices:
                    se.plot_environments(isite, additional_condition=se.AC.ONLY_ACB)
            if no_vis:
                test = input('Go to next structure ? ("y" to do so)')
                if test == 'y':
                    break
                continue
            test = input('View structure with environments ? ("y" for the unit cell or "m" for a supercell or "n") : ')
            if test in ['y', 'm']:
                if test == 'm':
                    mydeltas = []
                    test = input('Enter multiplicity (e.g. 3 2 2) : ')
                    nns = test.split()
                    for i0 in range(int(nns[0])):
                        for i1 in range(int(nns[1])):
                            for i2 in range(int(nns[2])):
                                mydeltas.append(np.array([1.0*i0, 1.0*i1, 1.0*i2], np.float))
                else:
                    mydeltas = [np.zeros(3, np.float)]
                if firsttime:
                    vis = StructureVis(show_polyhedron=False, show_unit_cell=True)
                    vis.show_help = False
                    firsttime = False
                vis.set_structure(se.structure)
                strategy.set_structure_environments(se)
                for isite, site in enumerate(se.structure):
                    try:
                        ces = strategy.get_site_coordination_environments(site)
                    except NeighborsNotComputedChemenvError:
                        continue
                    if len(ces) == 0:
                        continue
                    ce = strategy.get_site_coordination_environment(site)
                    if ce is not None and ce[0] != UNCLEAR_ENVIRONMENT_SYMBOL:
                        for mydelta in mydeltas:
                            psite = PeriodicSite(site._species, site._fcoords + mydelta, site._lattice,
                                                 properties=site._properties)
                            vis.add_site(psite)
                            neighbors = strategy.get_site_neighbors(psite)
                            draw_cg(vis, psite, neighbors, cg=lgf.cg.get_geometry_from_mp_symbol(ce[0]),
                                    perm=ce[1]['permutation'])
                vis.show()
            test = input('Go to next structure ? ("y" to do so) : ')
            if test == 'y':
                break
        print('')
