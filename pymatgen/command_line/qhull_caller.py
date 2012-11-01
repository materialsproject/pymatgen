#!/usr/bin/env python

'''
Interface with command line qhull.
Needs qhull installed. You can get it from http://www.qhull.org/.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ = "Sep 23, 2011"

import subprocess
import itertools
import math

import numpy as np


def run_qhull_command(command, data, proc_command=int, output_skip=1):
    """
    Helper function for actual qconvex and qvoronoi and qvertex commands.
    """
    assert len(data) > 0,"Data is empty"
    prep_str = str(len(data[0])) + "\n"
    prep_str += str(len(data)) + "\n"
    prep_str += "\n".join([' '.join([str(i) for i in row]) for row in data])
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE, close_fds=True)
    #print prep_str
    output = p.communicate(input=prep_str)[0]
    output = output.split("\n")
    for i in xrange(output_skip):
        output.pop(0)
    results = list()
    for row in output:
        cleanrow = row.strip()
        if cleanrow != "":
            results.append([proc_command(i) for i in cleanrow.split()])
    return results


def qconvex(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns the facets of the convex hull as a list of a list of integers.

    Args:
        data:
            Sequence of sequence of floats, e.g.[[0.0,0.1,0.2], [0.3,1,0.5],
            ...]

    Returns:
        facets as a list of list of int. e.g., [[1,2,3],[4,5,6],..]
    """
    return run_qhull_command(['qconvex', 'i', 'Qt'], data, int, 1)


def qvoronoi(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns voronoi results as a list of a list of integers.

    Args:
        data:
            Sequence of sequence of floats, e.g.[[0.0, 0.1, 0.2],
            [0.3, 1, 0.5], ...]

    Returns:
        facets as a list of list of int. e.g., [[1, 2, 3], [4, 5, 6], ..]
    """
    return run_qhull_command(['qvoronoi', 'Fv'], data, int, 1)


def qvertex(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns the facets of voronoi construction as a list of a list of float.

    Args:
        data:
            Sequence of sequence of floats, e.g.[[0.0,0.1,0.2], [0.3,1,0.5],
            ...]

    Returns:
        facets as a list of list of int. e.g., [[1.2,2.4,3.5],[4,5,6],..]
    """
    return run_qhull_command(['qvoronoi', 'p'], data, float, 2)


def qvertex_target(data, index):
    """
    Input data should be in the form of a list of a list of floats.
    index is the index of the targeted point
    Returns the vertices of the voronoi construction around this target point.
    """
    return run_qhull_command(['qvoronoi', 'p QV' + str(index)], data, float, 2)


def get_lines_voronoi(data):
    prep_str = str(len(data[0])) + "\n"
    prep_str += str(len(data)) + "\n"
    prep_str += "\n".join([' '.join([str(i) for i in row]) for row in data])
    p = subprocess.Popen(['qconvex', 'o'], stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE, close_fds=True)
    output = p.communicate(input=prep_str)[0]
    output = output.split("\n")

    nb_points = int(output[1].split(" ")[0])
    list_lines = []
    list_points = []
    for i in range(2, 2 + nb_points):
        list_points.append([float(c) for c in output[i].strip().split()])
    facets = []
    for i in range(2 + nb_points, len(output)):
        if output[i] != '':
            tmp = output[i].strip().split(" ")
            facets.append([int(tmp[j]) for j in range(1, len(tmp))])

    for i in range(len(facets)):
        for line in itertools.combinations(facets[i], 2):
            for j in range(len(facets)):
                if i != j and line[0] in facets[j] and line[1] in facets[j]:
                    #check if the two facets i and j are not coplanar
                    vector1 = np.array(list_points[facets[j][0]]) \
                        - np.array(list_points[facets[j][1]])
                    vector2 = np.array(list_points[facets[j][0]]) \
                        - np.array(list_points[facets[j][2]])
                    n1 = np.cross(vector1, vector2)
                    vector1 = np.array(list_points[facets[i][0]]) \
                        - np.array(list_points[facets[i][1]])
                    vector2 = np.array(list_points[facets[i][0]]) \
                        - np.array(list_points[facets[i][2]])
                    n2 = np.cross(vector1, vector2)

                    dot = math.fabs(np.dot(n1, n2) / (np.linalg.norm(n1)
                                                      * np.linalg.norm(n2)))
                    if(dot < 1.05 and dot > 0.95):
                        continue
                    list_lines.append({'start': list_points[line[0]],
                                       'end': list_points[line[1]]})
                    break
    return list_lines
