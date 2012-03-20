#!/usr/bin/env python

'''
Interface with command line qhull.
Needs qhull installed. You can get it from http://www.qhull.org/.

As far as we know, no formal Python extension for higher dim hulls exist.
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
    p = subprocess.Popen(command, stdout = subprocess.PIPE, 
                         stdin = subprocess.PIPE, close_fds = True)
    #print prep_str
    output = p.communicate(input = prep_str)[0]
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

def qvertex_target(data, index):
    """
    Input data should be in the form of a list of a list of floats.
    index is the index of the targeted point
    Returns the vertices of the voronoi construction around this target point.
    """
    return run_qhull_command(['qvoronoi','p QV'+str(index)], data, float, 2)

def get_lines_voronoi(data):
    prep_str = str(len(data[0])) + "\n"
    prep_str += str(len(data)) +"\n"
    prep_str += "\n".join([' '.join([str(i) for i in row]) for row in data])
    #print prep_str
    p = subprocess.Popen(['qconvex', 'o'], stdout = subprocess.PIPE,stdin = subprocess.PIPE, close_fds = True)
    output = p.communicate(input = prep_str)[0]
    output = re.split("\n", output)
    #print output
    nb_points=int(output[1].split(" ")[0])
    points=[]
    list_lines=[]
    list_points=[]
    for i in range(2,2+nb_points):
        #if(not output[i]==''):
            #print output[i]
            #print [float(c) for c in output[i].strip().split(" ")]
            list_points.append([float(c) for c in output[i].strip().split()])
    #print len(list_points)
    facets=[]
    for i in range(2+nb_points,len(output)):
        if (output[i]==''):
            continue
        tmp=output[i].strip().split(" ")
        #if(len(tmp)<=4):
           # continue
        facets.append([int(tmp[j]) for j in range(1,len(tmp)) ])
        #for j in range(1,len(indices)-1):
            #check
            #list_lines.append({'start':list_points[indices[j]],'end':list_points[indices[j+1]]})
    #go through each combinations of two vertices and see if it is included in two facets
    import itertools
    import numpy
    import math
    #import pylab as plt
    #from mpl_toolkits.mplot3d import Axes3D
    
    for i in range(len(facets)):
        #print i
        #print "start"+str(facets[i])
        #vector1=numpy.array(list_points[facets[i][0]])
        #vector2=numpy.array(list_points[facets[i][1]])
        #n2=numpy.cross(vector1,vector2)
        #fig = plt.figure()
        #ax=Axes3D(fig)
        #ax.scatter([list_points[facets[i][j]][0] for j in range(len(facets[i]))],[list_points[facets[i][j]][1] for j in range(len(facets[i]))],[list_points[facets[i][j]][2] for j in range(len(facets[i]))],color='r')
        #plt.show()
        for line in itertools.combinations(facets[i],2):
            #print line[0]
            #print line[1]
            for j in range(len(facets)):
                if(i==j):
                    continue
                if(line[0] in facets[j] and line[1] in facets[j]):
                    #check if the two facets i and j are not coplanar
                    vector1=numpy.array(list_points[facets[j][0]])-numpy.array(list_points[facets[j][1]])
                    vector2=numpy.array(list_points[facets[j][0]])-numpy.array(list_points[facets[j][2]])
                    n1=numpy.cross(vector1,vector2)
                    vector1=numpy.array(list_points[facets[i][0]])-numpy.array(list_points[facets[i][1]])
                    vector2=numpy.array(list_points[facets[i][0]])-numpy.array(list_points[facets[i][2]])
                    n2=numpy.cross(vector1,vector2)
                    
                    dot=math.fabs(numpy.dot(n1,n2)/(numpy.linalg.norm(n1)*numpy.linalg.norm(n2)))
                    #print "common"+str(facets[j])
                    #print n1
                    #print n2
                    #print dot
                    if(dot<1.05 and dot>0.95):
                        continue
                    #print j
                   # print facets[j]
                    list_lines.append({'start':list_points[line[0]],'end':list_points[line[1]]})
                    break
    #print list_lines        
    return list_lines
    #list=run_qhull_command(['qconvex','o'], data, float, 1)
    #print list
    #return None
    
