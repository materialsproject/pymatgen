#!/usr/bin/env python

"""
This module provides classes for plotting PhaseDiagram objects.
"""

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="$Sep 23, 2011M$"

import math
import numpy as np
import re
import itertools

class PDPlotter(object):
    '''
    A plotter class for phase diagrams.
    '''

    def __init__(self, phasediagram):
        """
        Arguments:
            phasediagram - a PhaseDiagram object.
        """
        self._pd = phasediagram
        self._dim = len(self._pd.elements)
        self.lines = uniquelines(self._pd.facets)
        if self._dim < 2 or self._dim > 4:
            raise ValueError("Only 2-4 components supported!")

    @property
    def pd_plot_data(self):        
        '''
        Plot data for phase diagram.
        2-comp - Full hull with energies
        3/4-comp - Projection into 2D or 3D gibbs triangle.
        Returns:
            (lines, stable_entries, unstable_entries) 
                - lines is a list of list of coordinates for lines in the PD. 
                - stable_entries is a {coordinate : entry} for each stable node in the phase diagram. (Each coordinate can only have one stable phase)
                - unstable_entries is a {entry: coordinates} for all unstable nodes in the phase diagram.
        '''
        pd = self._pd
        entries = pd.qhull_entries
        data = np.array(pd.qhull_data)
        facetlines = self.lines
        lines = list()
        stable_entries = dict()
        for line in facetlines:
            entry1 = entries[line[0]]
            entry2 = entries[line[1]]
            if self._dim == 2:
                x = [data[line[0]][0], data[line[1]][0]]
                y = [pd.get_form_energy_per_atom(entry1), pd.get_form_energy_per_atom(entry2)]
                coord = [x,y]
            elif self._dim == 3:
                coord = triangular_coord(data[line, 0:2])
            else:
                coord = tet_coord(data[line, 0:3])
            lines.append(coord)
            labelcoord = list(zip(*coord))
            stable_entries[labelcoord[0]] = entry1
            stable_entries[labelcoord[1]] = entry2    

        allentries = pd.all_entries
        alldata = np.array(pd.all_entries_hulldata)
        unstable_entries = dict()
        stable = pd.stable_entries
        for i in xrange(0,len(allentries)):
            entry = allentries[i]
            if entry not in stable:
                if self._dim == 2:
                    x = [alldata[i][0], alldata[i][0]]
                    y = [pd.get_form_energy_per_atom(entry), pd.get_form_energy_per_atom(entry)]
                    coord = [x,y]
                elif self._dim == 3:
                    coord = triangular_coord([alldata[i, 0:2],alldata[i, 0:2]])
                else:
                    coord = tet_coord([alldata[i, 0:3],alldata[i, 0:3],alldata[i, 0:3]])
                labelcoord = list(zip(*coord))
                unstable_entries[entry] = labelcoord[0]

        return (lines, stable_entries, unstable_entries)

    def showplot(self):
        """
        Draws the plot using Matplotlib.
        """
        if self._dim <4:
            self._show_2d_plot()
        elif self._dim == 4:
            self._show_3d_plot()

    def _show_2d_plot(self):
        '''
        Shows the plot using pylab.  Usually I won't do imports in methods,
        but since plotting is a fairly expensive library to load and not all 
        machines have matplotlib installed, I have done it this way.
        '''
        import matplotlib.pyplot as plt
        from matplotlib.font_manager import FontProperties
        (lines, labels, unstable) = self.pd_plot_data
        for x, y in lines:
            plt.plot(x, y, 'bo-', linewidth=3, markeredgecolor='b', markerfacecolor='r', markersize=10)
        font = FontProperties()
        font.set_weight('bold')
        font.set_size(20)
        count = 1
                    
        if len(self._pd.elements) == 3:
            plt.axis('equal')
            plt.xlim((-0.1, 1.2))
            plt.ylim((-0.1, 1.0))
            plt.axis('off')
            legendstart = [1.0, 0.55]
        else:
            plt.xlim((-0.1, 1.4))
            legendstart = [1.1, 0.0]

        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = entry.name
            x = coords[0]
            if coords[0] >= math.sqrt(3) / 2:
                halign = 'left'
                x += 0.02
            else:
                halign = 'right'
                x += -0.02
            if coords[1] > 0:
                valign = 'bottom'
            else:
                valign = 'top'

            if len(entry.composition.elements) == 1:
                plt.text(x, coords[1], label, horizontalalignment=halign, verticalalignment=valign, fontproperties=font)
            else:
                plt.text(x, coords[1], str(count), horizontalalignment=halign, verticalalignment=valign, fontproperties=font)
                plt.text(legendstart[0], legendstart[1]-0.05 * count, str(count) + " : " + label, horizontalalignment='left', verticalalignment='top', fontproperties=font)
                count += 1

        for entry,coords in unstable.items():
            label = entry.name
            plt.plot(coords[0], coords[1], 'bx', linewidth=3, markeredgecolor='b', markerfacecolor='b', markersize=10)

        F = plt.gcf()
        F.set_size_inches((8, 6.4))
        plt.show()

    def _show_3d_plot(self):
        '''
        Shows the plot using pylab.  Usually I won't do imports in methods,
        but since plotting is a fairly expensive library to load and not all 
        machines have matplotlib installed, I have done it this way.
        '''
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as p3
        from matplotlib.font_manager import FontProperties
        fig=plt.figure()
        ax = p3.Axes3D(fig)
        font = FontProperties()
        font.set_weight('bold')
        font.set_size(20)
        (lines, labels, unstable) = self.pd_plot_data
        count = 1
        newlabels = list()
        for x, y, z in lines:
            ax.plot(x, y, z, 'bo-', linewidth=3, markeredgecolor='b', markerfacecolor='r', markersize=10)
        for coords in sorted(labels.keys()):
            entry = labels[coords]
            label = entry.name
            if len(entry.composition.elements) == 1:
                # commented out options are only for matplotlib 1.0.  Removed them so that most ppl can use this class.
                ax.text(coords[0], coords[1], coords[2], label)#, horizontalalignment=halign, verticalalignment=valign, fontproperties=font)
            else:
                ax.text(coords[0], coords[1], coords[2], str(count))#, horizontalalignment=halign, verticalalignment=valign, fontproperties=font)
                newlabels.append(str(count) + " : " + label)
                count += 1
        plt.figtext(0.01,0.01,'\n'.join(newlabels))
        ax.axis('off')
        plt.show()

    def write_image(self, stream, image_format = "svg"):
        '''
        Writes the phase diagram to an image in a stream.
        Arguments:
            stream - stream to write to. Can be a file stream or a StringIO stream.
            image_format - format for image. CAn be any of matplotlib supported formats. Defaults to svg for best results for vector graphics.
        '''
        (lines, labels, unstable) = self.pd_plot_data
        dim = len(self._pd.elements)
        elementref = re.compile("^[A-Z][a-z]*$")
        count = 1
        import matplotlib as mpl
        from matplotlib.font_manager import FontProperties
            
        # chose a non-GUI backend
        mpl.use( 'Agg' )
        import matplotlib.pyplot as plt
        font = FontProperties()
        font.set_weight('bold')
        font.set_size(20)
        
        if dim == 4:
            plt.clf()
            plt.cla()
            import mpl_toolkits.mplot3d.axes3d as p3
            fig=plt.figure()
            ax = p3.Axes3D(fig)
            
            newlabels = list()
            for x, y, z in lines:
                ax.plot(x, y, z, 'bo-', linewidth=4, markeredgecolor='b', markerfacecolor='r', markersize=12)
            for coords in sorted(labels.keys()):
                label = labels[coords].name
                if elementref.match(label):
                    ax.text(coords[0], coords[1], coords[2], label, fontproperties=font)
                else:
                    ax.text(coords[0], coords[1], coords[2], str(count), fontproperties=font)
                    newlabels.append(str(count) + " : " + label)
                    count += 1
            plt.figtext(0.01,0.01,'\n'.join(newlabels), fontproperties=font)
        
        elif dim < 4 and dim > 1:
            plt.clf()
            plt.cla()
            
            for x,y in lines:
                plt.plot(x,y,'bo-',linewidth=4,markeredgecolor='b',markerfacecolor='r',markersize=12)
            if dim == 3:
                plt.axis('equal')
                plt.xlim( (-0.02, 1.18) )
                plt.ylim( (-0.1, 1.0) )
                plt.axis('off')
                legendstart = [1.0, 1.0]
                legendspacing = 0.05
            else:
                plt.xlim((-0.1, 1.4))
                legendstart = [1.1, 0.0]
            ymin,ymax = plt.ylim()
            legendspacing = (ymax - ymin)/len(labels)
        
            for coords in sorted(labels.keys()):
                label = labels[coords].name
                x = coords[0]
                if coords[0] >= math.sqrt(3)/2:
                    halign = 'left'
                    x += 0.02
                else:
                    halign = 'right'
                    x += -0.02
                if coords[1] > 0:
                    valign = 'bottom'
                else:
                    valign = 'top'
        
                if elementref.match(label):
                    plt.text(x, coords[1], label,horizontalalignment=halign,verticalalignment=valign,fontproperties=font)
                else:
                    plt.text(x, coords[1], str(count),horizontalalignment=halign,verticalalignment=valign,fontproperties=font)
                    plt.text(legendstart[0], legendstart[1]-legendspacing*count, str(count) + " : "+label,horizontalalignment='left',verticalalignment='top',fontproperties=font)
                    count +=1
        f = plt.gcf()
        f.set_size_inches( (12, 10) )
        
        plt.savefig(stream, format=image_format)

def uniquelines(q):
    '''
    Given all the facets, convert it into a set of unique lines.  Specifically used for converting convex hull facets into line pairs of coordinates.
    Arguments:
        q - a 2-dim sequence, where each row represents a facet. E.g., [[1,2,3],[3,6,7],...]
    Returns:
        setoflines - a set of tuple of lines.  E.g., {(1,2), (1,3), (2,3), ....}
    '''
    setoflines = set()
    for facets in q:
        for line in itertools.combinations(facets,2):
            setoflines.add(tuple(line))
    return setoflines


def triangular_coord(coord):
    '''
    Convert a two component coordinate into a triangle based coordinate system
    for a prettier phase diagram.
    Arguments:
        coordinate - coordinate used in the convex hull computation.
    Returns:
        coordinates in a triangular-based coordinate system. 
    '''
    unitvec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])
    result = np.dot(np.array(coord),unitvec)
    return result.transpose()

def tet_coord(coord):
    '''
    Convert a four component coordinate into a tetrahedron based coordinate system
    for a prettier phase diagram.
    Arguments:
        coordinate - coordinate used in the convex hull computation.
    Returns:
        coordinates in a tetrahedron-based coordinate system. 
    '''
    unitvec = np.array([[1,0,0],[0.5,math.sqrt(3)/2,0],[0.5,1.0/3.0*math.sqrt(3)/2,math.sqrt(6)/3]])
    result = np.dot(np.array(coord),unitvec)
    return result.transpose()
