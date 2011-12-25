#!/usr/bin/env python

'''
Interface with command line qhull.
Right now only tested on Linux systems.
'''

from __future__ import division

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import subprocess
import re

def run_qhull_command(command, data, proc_command = int, output_skip=1):
    """
    Helper function for actual qconvex and qvoronoi and qvertex commands.
    """
    prep_str = str(len(data[0])) + "\n"
    prep_str += str(len(data)) +"\n"
    prep_str += "\n".join([' '.join([str(i) for i in row]) for row in data])
    p = subprocess.Popen(command,stdout=subprocess.PIPE,stdin=subprocess.PIPE, close_fds=True)
    output = p.communicate(input=prep_str)[0]
    output = re.split("\n", output)
    for i in xrange(output_skip):
        output.pop(0)
    results = list()
    for row in output:
        cleanrow = row.strip()
        if cleanrow != "":
            results.append([proc_command(i) for i in re.split("\s+",cleanrow)])
    return results

def qconvex(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns the facets of the convex hull as a list of a list of integers.
    """
    return run_qhull_command(['qconvex','i','Qt'], data, int, 1)

def qvoronoi(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns voronoi results as a list of a list of integers.
    """
    return run_qhull_command(['qvoronoi','Fv'], data, int, 1)

def qvertex(data):
    """
    Input data should be in the form of a list of a list of floats.
    Returns the facets of voronoi construction as a list of a list of float.
    """
    return run_qhull_command(['qvoronoi','p'], data, float, 2)

