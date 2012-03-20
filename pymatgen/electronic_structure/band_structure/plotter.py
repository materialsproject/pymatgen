#!/usr/bin/env python

"""
This module provides classes and utilities to plot band structures
"""




__author__="Geoffroy Hautier, Shyue Ping Ong, Michael Kocher"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ ="March 14, 2012"

class BSPlotter(object):
    
    """
    class used to plot or get data to facilitate the plot of band structure line objects
    """
    
    def __init__(self, bs):
        """
        Arguments:
        bs - a Bandstructure_line object.
        """
        self._bs=bs
        self._nb_bands=bs._nb_bands
    
    @property
    def bs_plot_data(self):
        
        """
        get the data nicely formatted for a plot
        returns a dict:
            'ticks': a dictionnary with the 'distances' at which there is a kpoint (the x axis) and the labels (None if no label)
            'energy': an array (one element for each band) of energy for each kpoint
            'occup': similar to energy but giving occupations
        """
        
           
        energy=[]
        occup=[]
        distance=[self._bs._distance[j] for j in range(len(self._bs._kpoints))]
        ticks=self.get_ticks()
        for i in range(self._nb_bands):
            #pylab.plot([self._distance[j] for j in range(len(self._kpoints))],[self._bands[i]['energy'][j] for j in range(len(self._kpoints))],'b-',linewidth=5)
            energy.append([self._bs._bands[i]['energy'][j] for j in range(len(self._bs._kpoints))])
            occup.append([self._bs._bands[i]['occup'][j] for j in range(len(self._bs._kpoints))])
        
        return {'ticks':ticks,'distances':distance,'energy':energy,'occup':occup}
        
    def showplot(self):
        """
        plot the band structure.show it on the screen
        """
        import pylab
        from matplotlib import rc
        rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':20})
        rc('text', usetex=True)
        pylab.figure
        data=self.bs_plot_data
        energy=[]
        occup=[]
        distance=[self._bs._distance[j] for j in range(len(self._bs._kpoints))]
        ticks=self.get_ticks()
        for i in range(self._nb_bands):
            pylab.plot(data['distances'],data['energy'][i],'b-',linewidth=5)
            
        ticks=self.get_ticks()
        for i in range(len(ticks['label'])):
            if(ticks['label'][i]!=None):
                pylab.axvline(ticks['distance'][i],color='k')
        
        #pylab.axhline(self._efermi, color='r')
            
        pylab.gca().set_xticks(ticks['distance'])
        pylab.gca().set_xticklabels(ticks['label'])
        pylab.xlabel('Kpoints', fontsize = 'large')
        pylab.ylabel('Energy(eV)', fontsize = 'large')
        if(self._bs.is_metal()==False):
            vbm=self._bs.get_vbm()
            cbm=self._bs.get_cbm()
            if(cbm['kpoint'].label!=None):
                for i in range(len(self._bs._kpoints)):
                    if(self._bs._kpoints[i].label==cbm['kpoint'].label):
                        pylab.scatter(self._bs._distance[i],cbm['energy'],color='r',marker='o',s=100)
            else:
                pylab.scatter(self._bs._distance[cbm['kpoint_index']],cbm['energy'],color='r',marker='o',s=100)
                
            if(vbm['kpoint'].label!=None):
                for i in range(len(self._bs._kpoints)):
                    if(self._bs._kpoints[i].label==vbm['kpoint'].label):
                        pylab.scatter(self._bs._distance[i],vbm['energy'],color='G',marker='o',s=100)
            else:
                pylab.scatter(self._bs._distance[vbm['kpoint_index']],vbm['energy'],color='g',marker='o',s=100)
            
            pylab.ylim(vbm['energy']-4,cbm['energy']+4)
        
        else:
            pylab.axhline(self._bs._efermi, color='r')
            pylab.ylim(self._bs._efermi-4,self._bs._efermi+4)
            
        pylab.legend()
        pylab.show()
        
        
        
        
    def get_ticks(self):
        """
        get all ticks and labels for a band structure plot
        """
        tick_distance=[]
        tick_labels=[]
        previous_label=self._bs._kpoints[0].label
        previous_branch=self._bs.get_branch_name(0)
        for i in range(len(self._bs._kpoints)):
            c=self._bs._kpoints[i]
            if(c.label!=None):
                tick_distance.append(self._bs._distance[i])
                if(c.label!=previous_label and previous_branch!=self._bs.get_branch_name(i)):
                    label1=c.label
                    if(label1.startswith("\\") or label1.find("_")!=-1):
                        label1="$"+label1+"$"
                    label0=previous_label
                    if(label0.startswith("\\") or label0.find("_")!=-1):
                        label0="$"+label0+"$"
                    tick_labels.pop()
                    tick_distance.pop()
                    tick_labels.append(label0+"$|$"+label1)
                    #print label0+","+label1
                else:
                    if(c.label.startswith("\\") or c.label.find("_")!=-1):
                        tick_labels.append("$"+c.label+"$")
                    else:
                        tick_labels.append(c.label)
                previous_label=c.label
                previous_branch=self._bs.get_branch_name(i)
        return {'distance':tick_distance,'label':tick_labels} 
    
    def plot_compare(self,other_plotter):
        """
        plot two band structure for comparison.
        TODO: still a lot of work to do that nicely!
        """
        import pylab
        data=self.bs_plot_data
        data_other=other_plotter.bs_plot_data
        for i in range(self._nb_bands):
            pylab.plot(data['distances'],data['energy'][i],'b-',linewidth=3)
            
        for i in range(self._nb_bands):
            pylab.plot(data['distances'],data_other['energy'][i],'r--',linewidth=3)
        
        
        ticks=self.get_ticks()
        
        pylab.gca().set_xticks(ticks['distance'])
        pylab.gca().set_xticklabels(ticks['label'])
        pylab.xlabel('Kpoints', fontsize = 'large')
        pylab.ylabel('Energy(eV)', fontsize = 'large')
        #pylab.ylim(vbm-4,cbm+4)
        for i in range(len(ticks['label'])):
            if(ticks['label'][i]!=None):
                pylab.axvline(ticks['distance'][i],color='k')
        pylab.show()
        pylab.legend()
    
    def plot_brillouin(self):
        import pylab as plt
        import pymatgen.command_line.qhull_caller
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax=Axes3D(fig)
        vec1=self._bs._lattice_rec.matrix[0]
        vec2=self._bs._lattice_rec.matrix[1]
        vec3=self._bs._lattice_rec.matrix[2]
        #ax.plot([0,vec1[0]],[0,vec1[1]],[0,vec1[2]],color='k')
        #ax.plot([0,vec2[0]],[0,vec2[1]],[0,vec2[2]],color='k')
        #ax.plot([0,vec3[0]],[0,vec3[1]],[0,vec3[2]],color='k')
        
        #make the grid
        max_x=-1000
        max_y=-1000
        max_z=-1000
        min_x=1000
        min_y=1000
        min_z=1000
        list_k_points=[]
        for i in[-1,0,1]:
            for j in [-1,0,1]:
                for k in [-1,0,1]:
                    list_k_points.append(i*vec1+j*vec2+k*vec3)
                    if(list_k_points[-1][0]>max_x):
                        max_x=list_k_points[-1][0]
                    if(list_k_points[-1][1]>max_y):
                        max_y=list_k_points[-1][1]
                    if(list_k_points[-1][2]>max_z):
                        max_z=list_k_points[-1][0]
                    
                    if(list_k_points[-1][0]<min_x):
                        min_x=list_k_points[-1][0]
                    if(list_k_points[-1][1]<min_y):
                        min_y=list_k_points[-1][1]
                    if(list_k_points[-1][2]<min_z):
                        min_z=list_k_points[-1][0]    
                        
                    #ax.scatter([list_k_points[-1][0]],[list_k_points[-1][1]],[list_k_points[-1][2]])
        #plt.show()
        vertex=pymatgen.command_line.qhull_caller.qvertex_target(list_k_points,13)
        #print vertex
        lines=pymatgen.command_line.qhull_caller.get_lines_voronoi(vertex)
        #[vertex[i][0] for i in range(len(vertex))],[vertex[i][1] for i in range(len(vertex))]+" "+str(vertex[i][2])
        #ax.scatter([vertex[i][0] for i in range(len(vertex))],[vertex[i][1] for i in range(len(vertex))],[vertex[i][2] for i in range(len(vertex))],color='r')
        for i in range(len(lines)):
            vertex1=lines[i]['start']
            vertex2=lines[i]['end']
            ax.plot([vertex1[0],vertex2[0]],[vertex1[1],vertex2[1]],[vertex1[2],vertex2[2]],color='k')
        
        
        for b in self._bs._branches:
            vertex1=self._bs._kpoints[b['start_index']].cart_coords
            vertex2=self._bs._kpoints[b['end_index']].cart_coords
            ax.plot([vertex1[0],vertex2[0]],[vertex1[1],vertex2[1]],[vertex1[2],vertex2[2]],color='r',linewidth=3)
        #plot the labelled points    
        
        for k in self._bs._kpoints:
            if(not k.label==None):
                label=k.label
                if(k.label.startswith("\\") or k.label.find("_")!=-1):
                       label="$"+k.label+"$"
                off=0.01
                ax.text(k.cart_coords[0]+off,k.cart_coords[1]+off,k.cart_coords[2]+off,label,color='b')
                ax.scatter([k.cart_coords[0]],[k.cart_coords[1]],[k.cart_coords[2]],color='b')
        
        # make ticklabels and ticklines invisible
        for a in ax.w_xaxis.get_ticklines()+ax.w_xaxis.get_ticklabels():
            a.set_visible(False) 
        for a in ax.w_yaxis.get_ticklines()+ax.w_yaxis.get_ticklabels():
            a.set_visible(False)
        for a in ax.w_zaxis.get_ticklines()+ax.w_zaxis.get_ticklabels():
            a.set_visible(False) 
        
        #ax.set_xlim3d(0.5*min_x, 0.5*max_x) 
        #ax.set_ylim3d(0.5*min_y, 0.5*max_y) 
        #ax.set_zlim3d(0.5*min_z, 0.5*max_z) 
        ax.grid(False) 
        plt.show()
        