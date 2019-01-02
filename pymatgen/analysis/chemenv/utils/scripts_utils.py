# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.



from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser
try:
    import vtk
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
import re
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import UNCLEAR_ENVIRONMENT_SYMBOL
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import AbstractGeometry
from pymatgen.analysis.chemenv.utils.coordination_geometry_utils import rotateCoords
from pymatgen.analysis.chemenv.utils.defs_utils import chemenv_citations
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy
#from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimpleAbundanceChemenvStrategy
#from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import TargettedPenaltiedAbundanceChemenvStrategy
from pymatgen.core.structure import Molecule
from collections import OrderedDict
import numpy as np

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



strategies_class_lookup = OrderedDict()
strategies_class_lookup['SimplestChemenvStrategy'] = SimplestChemenvStrategy

#strategies_class_lookup['SimpleAbundanceChemenvStrategy'] = SimpleAbundanceChemenvStrategy
#strategies_class_lookup['TargettedPenaltiedAbundanceChemenvStrategy'] = TargettedPenaltiedAbundanceChemenvStrategy


def draw_cg(vis, site, neighbors, cg=None, perm=None, perfect2local_map=None,
            show_perfect=False, csm_info=None, symmetry_measure_type='csm_wcs_ctwcc', perfect_radius=0.1,
            show_distorted=True, faces_color_override=None):
    if show_perfect:
        if csm_info is None:
            raise ValueError('Not possible to show perfect environment without csm_info')
        csm_suffix = symmetry_measure_type[4:]
        perf_radius = (perfect_radius - 0.2) / 0.002
    if perm is not None and perfect2local_map is not None:
        raise ValueError('Only "perm" or "perfect2local_map" should be provided in draw_cg, not both')
    if show_distorted:
        vis.add_bonds(neighbors, site)
        for n in neighbors:
            vis.add_site(n)
    if len(neighbors) < 3:
        if show_distorted:
            vis.add_bonds(neighbors, site, color=[0.0, 1.0, 0.0], opacity=0.4, radius=0.175)
        if show_perfect:
            if len(neighbors) == 2:
                perfect_geometry = AbstractGeometry.from_cg(cg)
                trans = csm_info['other_symmetry_measures']['translation_vector_{}'.format(csm_suffix)]
                rot = csm_info['other_symmetry_measures']['rotation_matrix_{}'.format(csm_suffix)]
                scale = csm_info['other_symmetry_measures']['scaling_factor_{}'.format(csm_suffix)]
                points = perfect_geometry.points_wcs_ctwcc()
                rotated_points = rotateCoords(points, rot)
                points = [scale * pp + trans for pp in rotated_points]
                if 'wcs' in csm_suffix:
                    ef_points = points[1:]
                else:
                    ef_points = points
                edges = cg.edges(ef_points, input='coords')
                vis.add_edges(edges, color=[1.0, 0.0, 0.0])
                for point in points:
                    vis.add_partial_sphere(coords=point, radius=perf_radius, color=[0.0, 0.0, 0.0],
                                           start=0, end=360, opacity=1)
    else:
        if show_distorted:
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
            if faces_color_override:
                mycolor = faces_color_override
            else:
                mycolor = [float(i) / 255 for i in vis.el_color_mapping[symbol]]
            vis.add_faces(faces, mycolor, opacity=0.4)
            vis.add_edges(edges)
        if show_perfect:
            perfect_geometry = AbstractGeometry.from_cg(cg)
            trans = csm_info['other_symmetry_measures']['translation_vector_{}'.format(csm_suffix)]
            rot = csm_info['other_symmetry_measures']['rotation_matrix_{}'.format(csm_suffix)]
            scale = csm_info['other_symmetry_measures']['scaling_factor_{}'.format(csm_suffix)]
            points = perfect_geometry.points_wcs_ctwcc()
            rotated_points = rotateCoords(points, rot)
            points = [scale*pp + trans for pp in rotated_points]
            if 'wcs' in csm_suffix:
                ef_points = points[1:]
            else:
                ef_points = points
            edges = cg.edges(ef_points, input='coords')
            vis.add_edges(edges, color=[1.0, 0.0, 0.0])
            for point in points:
                vis.add_partial_sphere(coords=point, radius=perf_radius, color=[0.0, 0.0, 0.0],
                                       start=0, end=360, opacity=1)




# Visualizing a coordination geometry
def visualize(cg, zoom=None, vis=None, myfactor=1.0, view_index=True, faces_color_override=None):
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
    draw_cg(vis, site=structure[0], neighbors=structure[1:], cg=cg, faces_color_override=faces_color_override)
    if view_index:
        for ineighbor, neighbor in enumerate(structure[1:]):
            vis.add_text(neighbor.coords, '{}'.format(ineighbor), color=(0, 0, 0))
    if zoom is not None:
        vis.zoom(zoom)
    return vis


def welcome(chemenv_config):
    print('Chemical Environment package (ChemEnv)')
    print(chemenv_citations())
    print(chemenv_config.package_options_description())


def thankyou():
    print('Thank you for using the ChemEnv package')
    print(chemenv_citations())


def compute_environments(chemenv_configuration):
    string_sources = {'cif': {'string': 'a Cif file', 'regexp': r'.*\.cif$'},
                      'mp': {'string': 'the Materials Project database',
                             'regexp': r'mp-[0-9]+$'}}
    questions = {'c': 'cif'}
    questions['m'] = 'mp'
    lgf = LocalGeometryFinder()
    lgf.setup_parameters()
    allcg = AllCoordinationGeometries()
    strategy_class = strategies_class_lookup[chemenv_configuration.package_options['default_strategy']['strategy']]
    #TODO: Add the possibility to change the parameters and save them in the chemenv_configuration
    default_strategy = strategy_class()
    default_strategy.setup_options(chemenv_configuration.package_options['default_strategy']['strategy_options'])
    max_dist_factor = chemenv_configuration.package_options['default_max_distance_factor']
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
            a = MPRester()
            structure = a.get_structure_by_material_id(input_source)
        lgf.setup_structure(structure)
        print('Computing environments for {} ... '.format(structure.composition.reduced_formula))
        se = lgf.compute_structure_environments(maximum_distance_factor=max_dist_factor)
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
                            csm = ce[1]['other_symmetry_measures']['csm_wcs_ctwcc']
                            mystring += ' - {} ({}): {:.2f} % (csm : {:2f})\n'.format(cg.name, cg.mp_symbol,
                                                                                      100.0*ce[2],
                                                                                      csm)
                    if test in ['d', 'g'] and strategy.uniquely_determines_coordination_environments:
                        if thecg.mp_symbol != UNCLEAR_ENVIRONMENT_SYMBOL:
                            mystring += '  <Continuous symmetry measures>  '
                            mingeoms = se.ce_list[isite][thecg.coordination_number][0].minimum_geometries()
                            for mingeom in mingeoms:
                                csm = mingeom[1]['other_symmetry_measures']['csm_wcs_ctwcc']
                                mystring += '{} : {:.2f}       '.format(mingeom[0], csm)
                    print(mystring)
            if test == 'g':
                while True:
                    test = input('Enter index of site(s) (e.g. 0 1 2, separated by spaces) for which you want to see the grid of parameters : ')
                    try:
                         indices=[int(x) for x in test.split()]
                         print(str(indices))
                         for isite in indices:
                             if isite <0:
                                 raise IndexError
                             se.plot_environments(isite, additional_condition=se.AC.ONLY_ACB)
                         break
                    except ValueError:
                         print('This is not a valid site')
                    except IndexError:
                         print('This site is out of the site range')


            if no_vis:
                test = input('Go to next structure ? ("y" to do so)')
                if test == 'y':
                    break
                continue
            test = input('View structure with environments ? ("y" for the unit cell or "m" for a supercell or "n") : ')
            if test in ['y', 'm']:
                if test == 'm':
                    mydeltas = []
                    while True:
                        try:
                            test = input('Enter multiplicity (e.g. 3 2 2) : ')
                            nns = test.split()
                            for i0 in range(int(nns[0])):
                                for i1 in range(int(nns[1])):
                                    for i2 in range(int(nns[2])):
                                        mydeltas.append(np.array([1.0*i0, 1.0*i1, 1.0*i2], np.float))
                            break

                        except (ValueError,IndexError):
                            print('Not a valid multiplicity')
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
                            draw_cg(vis, psite, neighbors, cg=lgf.allcg.get_geometry_from_mp_symbol(ce[0]),
                                    perm=ce[1]['permutation'])
                vis.show()
            test = input('Go to next structure ? ("y" to do so) : ')
            if test == 'y':
                break
        print('')
