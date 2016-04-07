#!/usr/bin/env python

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.analysis.chemenv.utils.scripts_utils import visualize


if __name__ == '__main__':
    allcg = AllCoordinationGeometries()
    vis = None
    while True:
        cg_symbol = raw_input('Enter symbol of the geometry you want to see or "q" to quit : ')
        if cg_symbol == 'q':
            break
        try:
            cg = allcg[cg_symbol]
        except LookupError:
            print('Wrong geometry, try again ...')
            continue
        print(cg.name)
        for ipoint, point in enumerate(cg.points):
            print('Point #{:d} : {} {} {}'.format(ipoint, repr(point[0]), repr(point[1]), repr(point[2])))
        print('Algorithms used :')
        for ialgo, algo in enumerate(cg.algorithms):
            print('Algorithm #{:d} :'.format(ialgo))
            print(algo)
            print('')
        if vis is None:
            vis = visualize(cg=cg, zoom=1.0)
        else:
            vis = visualize(cg=cg, vis=vis)

