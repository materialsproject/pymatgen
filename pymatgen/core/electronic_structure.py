#!/usr/bin/env python

"""
This module provides classes to define electronic structure, such as the density of states, etc.
"""

from __future__ import division

__author__="Shyue Ping Ong, Vincent L Chevrier, Rickard Armiento, Geoffroy Hautier"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__status__ = "Production"
__date__ ="Sep 23, 2011"

import math
import numpy as np
import scipy.interpolate as spint
import pymatgen.command_line.aconvasp_caller

class _SpinImpl(object):
    """
    Internal representation of a Spin. 
    Do not use directly.
    """
       
    def __init__(self,name):
        self._name = name
            
    def __int__(self):
        return 1 if self._name == "up" else -1
        
    def __repr__(self):
        return self._name
        
    def __eq__(self, other):
        if other == None:
            return False
        return self._name == other._name
        
    def __hash__(self):
        return self.__int__()
           
    def __str__(self):
        return self._name

class Spin(object):
    """
    Enum type for Spin.  Only up and down.  Design follows somewhat the familiar Java syntax.
    """
    
    up = _SpinImpl("up")
    down = _SpinImpl("down")
    all_spins = (up, down)
    
    @staticmethod
    def from_int(i):
        if i == 1:
            return Spin.up
        elif i == -1:
            return Spin.down
        else:
            raise ValueError("Spin integers must be 1 or -1")
        
class _OrbitalImpl(object):
    """
    Internal representation of an orbital.  Do not use directly. 
    Use the Orbital class enum types.
    """
    
    def __init__(self,name, vasp_index):
        self._name = name
        self._vasp_index = vasp_index
            
    def __int__(self):
        return self._vasp_index
        
    def __repr__(self):
        return self._name
        
    def __eq__(self, other):
        if other == None:
            return False
        return self._name == other._name
        
    def __hash__(self):
        return self.__int__()
        
    @property
    def orbital_type(self):
        return self._name[0].upper()
        
    def __str__(self):
        return self._name

class Orbital(object):
    """
    Enum type for OrbitalType. Indices are basically the azimutal quantum number, l.
    Design follows somewhat the familiar Java syntax.
    """
    
    s = _OrbitalImpl("s",0)
    py = _OrbitalImpl("py",1)
    pz = _OrbitalImpl("pz",2)
    px = _OrbitalImpl("px",3)
    dxy = _OrbitalImpl("dxy",4)
    dyz = _OrbitalImpl("dyz",5)
    dz2 = _OrbitalImpl("dz2",6)
    dxz = _OrbitalImpl("dxz",7)
    dx2 = _OrbitalImpl("dx2",8)
    all_orbitals = (s, py, pz, px, dxy, dyz, dz2, dxz, dx2)
            
    @staticmethod
    def from_vasp_index(i):
        for orb in Orbital.all_orbitals:
            if int(orb) == i:
                return orb
    
    @staticmethod
    def from_string(orb_str):
        for orb in Orbital.all_orbitals:
            if str(orb) == orb_str:
                return orb
        raise ValueError("Illegal orbital definition!")

class Dos(object):

    def __init__(self, efermi, energies, densities):
        self._efermi = efermi
        self._energies = np.array(energies)
        self._dos = {k:np.array(d) for k,d in densities.items()} # should be a dict having as keys the Spin.up and Spin.down objects

    def get_densities(self, spin=None):
        if self._dos == None:
            result = None
        elif spin == None:
            if Spin.down in self._dos:
                result = self._dos[Spin.up] + self._dos[Spin.down]
            else:
                result = self._dos[Spin.up]
        else:
            result = self._dos[spin]
        return result

    @property
    def densities(self):
        return self._dos

    @property
    def efermi(self):
        return self._efermi

    @property
    def energies(self):
        return self._energies

    def __add__(self,other):
        """
        Adds two DOS together. Checks that energy scales are the same. Otherwise, a ValueError is thrown.
        """
        if not (self.energies == other.energies).all():
            raise ValueError("Energies of both DOS are not compatible!")
        densities = {spin: self._dos[spin] + other._dos[spin] for spin in self.densities.keys()}
        return Dos(self.efermi, self.energies, densities)

    def get_interpolated_value(self,energy):
        f = {}
        for spin in self._dos.keys():
            f[spin] = spint.interp1d(self._energies, self._dos[spin])(energy)
        return f
    
    def get_interpolated_gap(self,tol=0.001,abs_tol=False,spin=None):
        """
        Expects a DOS object and finds the gap
        
        Arguments:
            tol:
                tolerance in occupations for determining the gap
            abs_tol:
                tolerance an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the ap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            (gap, cbm, vbm) :  tuple of floats in eV corresponding to the gap, cbm and vbm
        """
        
        tdos = self.get_densities(spin)
        if abs_tol == False:
            tol = tol * tdos.sum() / tdos.shape[0]
        energies = self._energies
        below_fermi = [i for i in xrange(len(energies)) if energies[i] < self._efermi and tdos[i] > tol]
        above_fermi = [i for i in xrange(len(energies)) if energies[i] > self._efermi and tdos[i] > tol]
        vbm_start = max(below_fermi)
        cbm_start = min(above_fermi)
        
        if vbm_start == cbm_start:
            return 0.0,  self._efermi,self._efermi
        else:
            # Interpolate between adjacent values
            f = spint.interp1d(tdos[vbm_start:vbm_start+2][::-1], energies[vbm_start:vbm_start+2][::-1])
            start = f(tol)
            f = spint.interp1d(tdos[cbm_start-1:cbm_start+1], energies[cbm_start-1:cbm_start+1])
            end = f(tol)
            return end - start, end, start 

    def get_cbm_vbm(self, tol = 0.001, abs_tol = False, spin = None):
        """
        Expects a DOS object and finds the cbm and vbm.
        
        Args:
            tol: 
                tolerance in occupations for determining the gap
            abs_tol: 
                an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the gap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            (cbm, vbm): float in eV corresponding to the gap
        """
        #determine tolerance
        tdos = self.get_densities(spin)
        if abs_tol == False:
            tol = tol * tdos.sum() / tdos.shape[0]
    
        # find index of fermi energy
        i_fermi=0
        while (self._energies[i_fermi] <= self._efermi):
            i_fermi+=1
    
        # work backwards until tolerance is reached
        i_gap_start=i_fermi
        while i_gap_start-1 >=0 and tdos[i_gap_start-1] <= tol :
            i_gap_start -= 1
    
        # work forwards until tolerance is reached
        i_gap_end=i_gap_start
        while i_gap_end < tdos.shape[0] and tdos[i_gap_end] <= tol :
            i_gap_end += 1
        i_gap_end -= 1
        return (self._energies[i_gap_end],self._energies[i_gap_start])

    def get_gap(self, tol = 0.001, abs_tol = False, spin = None):
        """
        Expects a DOS object and finds the gap.
        
        Args:
            tol: 
                tolerance in occupations for determining the gap
            abs_tol: 
                an absolute tolerance (True) and a relative one (False)
            spin:
                Possible values are:
                    None - finds the gap in the summed densities
                    Up - finds the gap in the up spin channel
                    Down - finds the gap in teh down spin channel
        
        Returns:
            gap in eV
        """
        (cbm,vbm) = self.get_cbm_vbm(tol,abs_tol,spin)
        return max(cbm - vbm,0.0)

    def __str__(self):
        """
        Returns a string which can be easily plotted
        """
        if Spin.down in self._dos:
            stringarray = ["#%30s %30s %30s" %('Energy','DensityUp','DensityDown')]
            stringarray.extend(["%.5f %.5f %.5f" % (self._energies[i], self._dos[Spin.up][i], self._dos[Spin.down][i]) for i in range(len(self._energies))])
        else:
            stringarray = ["#%30s %30s" %('Energy','DensityUp')]
            stringarray.extend(["%.5f %.5f" % (self._energies[i], self._dos[Spin.up][i]) for i in range(len(self._energies))])
        return "\n".join(stringarray)

class PDos(Dos):
    """
    Projected DOS for a specific orbital. Extends the Dos object.
    """
    def __init__(self,efermi,energies,densities,orbital):
        Dos.__init__(self,efermi,energies,densities)
        self.orbital = orbital
        
    def __str__(self):
        return "#"+str(self.orbital) + "\n" + super(PDos,self).__str__()

class CompleteDos(Dos):
    """
    This wrapper class defines a total dos, and also provides a list of Pdos.
    Mainly used by pymatgen.io.vaspio.Vasprun to create a complete Dos from
    a vasprun.xml file. 
    """
    
    def __init__(self, structure, total_dos, pdoss):
        """
        Arguments:
            structure:
                Structure associated with this particular DOS.
            total_dos:
                total Dos for structure
            pdoss:
                a list of array of Pdos.  pdoss corresponds to site order in structure 
        """
        self._efermi = total_dos.efermi
        self._energies = total_dos.energies
        self._dos = total_dos.densities
        self._pdos = {structure[i]:{Orbital.from_vasp_index(j) : pdoss[i][j] for j in range(len(pdoss[i]))} for i in range(structure.num_sites)}
        self._structure = structure

    @property
    def structure(self):
        return self._structure

    def get_site_orbital_dos(self, site, orbital):
        return self._pdos[site][orbital]
        
    def get_site_dos(self, site):
        site_dos = None
        for pdos in self._pdos[site].values():
            if site_dos == None:
                site_dos = Dos(pdos.efermi, pdos.energies,pdos.densities)
            else:
                site_dos += pdos
        return site_dos

    def get_spd_dos(self):
        """
        Get orbital projected Dos.
        
        Returns:
            dict of {orbital: Dos}, e.g. {'s': Dos object, ...}
        """
        spd_dos = dict()
        for atom_dos in self._pdos.values():
            for pdos in atom_dos.values():
                orbital_type = pdos.orbital.orbital_type
                if orbital_type not in spd_dos:
                    spd_dos[orbital_type] = Dos(pdos.efermi, pdos.energies,pdos.densities)
                else:
                    spd_dos[orbital_type] += pdos
        return spd_dos

    def get_element_dos(self):
        """
        Get element projected Dos.
        
        Returns:
            dict of {Element: Dos}
        """
        
        el_dos = dict()
        for site, atom_dos in self._pdos.items():
            el = site.specie
            for pdos in atom_dos.values():
                if el not in el_dos:
                    el_dos[el] = Dos(pdos.efermi, pdos.energies,pdos.densities)
                else:
                    el_dos[el] += pdos
        return el_dos
    
    def __str__(self):
        return "Complete DOS for "+str(self._structure)
    
def plot_dos(dos_dict, zero_at_efermi = True):
    """
    Plots a series of Dos using matplotlib.
    
    Args:
        dos_dict:
            dict of {label: Dos}
        zero_at_efermi:
            Whether to shift all Dos to have zero energy at the fermi energy.
            Defaults to True.
    """
    import pylab
    color_order = ['r', 'b', 'g', 'y']
    count = 0
    for key,dos in dos_dict.items():
        energies = dos.energies - dos.efermi if zero_at_efermi else dos.energies
        densities = dos.densities
        if Spin.up in densities:
            pylab.plot(energies, densities[Spin.up], color_order[count % 4], label=str(key) + ' up')
        if Spin.down in densities:
            pylab.plot(energies, - densities[Spin.down], color_order[count % 4], label=str(key) + ' down')
        count += 1
    
    pylab.xlabel('Energies (eV)', fontsize = 'large')
    pylab.ylabel('Density of states', fontsize = 'large')
    pylab.legend()
    pylab.show()
    
class Bandstructure(object):
    
    """
    Created on Nov 2, 2011

    @author: geoffroy (geoffroy.hautier@uclouvain.be)
    
    The object stores a band structure along symmetry directions.
    
    the constructor takes:
    
    -kpoints: a list of kpoints in a numpy array format
    -eigenvals: a list of dictionnary (one per band), each dictionary has
    a list of float linked to an 'energy' key and a list of float linked to an 'occupation' key
    the order in which the list is set corresponds to the kpoints given in the kpoints list
    -labels_dict: a dictionnary of label associated with a kpoint
    -the structure as a Structure object
    -efermi the fermi energy
    
    Here is an example of call for the class (two kpoints, one band):
    
    kpoints=[np.array([0,0,0]),np.array([0.5,0.5,0.5])]
    eigenvals=[]
    eigenvals.append({'energy':[0,6.0],'occup':[1.0,1.0]})
    eigenvals.append({'energy':[2.0,7.0],'occup':[0.0,0.0]})
    labels_dict={'Gamma':np.array([0,0,0]),'X':np.array([0.5,0.5,0.5])
    structure= a structure object
    fermi=0.2
    
    bs_example=BandStructure(kpoints,eigenvals,labels_dict,rec_lattice,fermi)
    
    The kpoints are stored in self._kpoints as list of dictionary
    each element in the list stores 
    "kpoint": the actual kpoint in cartesian coordinates (not fractional and WITH the factor 2pi)
            stored as a numpy array.
    "distance": the distance along the x axes when band structures are plotted
    "label" the label of the kpoint (gamma, X etc...), None value if no label has been assigned
    "branch": the branch the kpoint belongs to (e.g., gamma-X)
    
    the energy bands are stored in self._bands, it's a list of dictionary, each dictionary relates
    to one band. The number of bands is stored also in self._nb_bands.
    each dictionary stores:
    "energy": a list of energy (ordered as the kpoints are ordered in self._kpoints), 
    "occup": a list of occupations for this band for each kpoint(ordered as the kpoints are ordered in self._kpoints)


    Many TODOs in this class:
    
    -add a way to enter k-points in fractional coords too
    -deal better with the need or not to pass a structure (compute the direct lattice from reciprocal?)
    -add spins
    -...

    """

    def __init__(self,kpoints,eigenvals,labels_dict,structure,efermi):
        
        self._efermi=efermi
        
        self._structure=structure
        
        """
        the structure as a Structure object
        """
        
        self._lattice_rec=structure.lattice.reciprocal_lattice
        
        """
        reciprocal lattice as a Lattice object
        """
        
        self._kpoints=[]
        self._kpoints=[{'kpoint':kpoints[i]} for i in range(len(kpoints))]
        """
        all kpoints, (order matter!)
        """
        
        
        self._labels_dict=labels_dict
        
        """
        the labels to the kpoints
        """
        
        self._branches=[]
        """
        list used to know what kpoints are in what branch
        """
        self._bands=eigenvals
        """
        all energy values for each band at the different kpoints
        """
        self._nb_bands=len(eigenvals)
        """
        the number of bands
        """
        self._branch_labels=set()
        """
        all branches labels (ex: Gamma-Z, etc...)
        """ 
        
        one_group=[]
        branches=[]
        #get labels and distance for each kpoint
        previous_kpoint=self._kpoints[0]['kpoint']
        previous_distance=0.0
        
        previous_label=self.get_label(self._kpoints[0]['kpoint'],self._labels_dict)
        for i in range(len(self._kpoints)):
            label=self.get_label(self._kpoints[i]['kpoint'],self._labels_dict)
            if(label!=None and previous_label!=None):
                self._kpoints[i]['distance']=previous_distance
            else:
                self._kpoints[i]['distance']=np.linalg.norm(self._kpoints[i]['kpoint']-previous_kpoint)+previous_distance
            previous_kpoint=self._kpoints[i]['kpoint']
            previous_distance=self._kpoints[i]['distance']
            label=self.get_label(self._kpoints[i]['kpoint'],self._labels_dict)
            self._kpoints[i]['label']=label
            
            if label!=None:
                if(previous_label!=None):
                    if(len(one_group)!=0):
                        branches.append(one_group)
                    one_group=[]
            previous_label=label
            one_group.append(i)
            
        if(len(one_group)!=0):
            branches.append(one_group)
        self._branches=branches
        
        #go through all the branches and assign their branch name to each kpoint
        for b in branches:
            label_start=self.get_label(self._kpoints[b[0]]['kpoint'],self._labels_dict)
            label_end=self.get_label(self._kpoints[b[-1]]['kpoint'],self._labels_dict)
            for c in b:
                self._kpoints[c]['branch']=label_start+"-"+label_end
                self._branch_labels.add(self._kpoints[c]['branch'])
          
            
    def get_label(self,kpoint,labels_dict):
        """
        simple method giving a label for any kpoint
        """
        for c in labels_dict:
            if(np.linalg.norm(kpoint-np.array(labels_dict[c]))<0.0001):
                return c
        return None
                
    def getVBM(self):
        """
        get the valence band minimum (VBM). returns a dictionnary with
        'band_index': a list of the indices of the band containing the VBM (please note that you can have several bands 
        sharing the VBM)
        'kpoint_index': the index in self._kpoints of the kpoint vbm
        'kpoint': the kpoint (in cartesian with 2pi factor)
        'energy': the energy of the VBM
        'label': the label of the vbm kpoint if any
        """
        max_tmp=-1000.0
        index=None
        for i in range(self._nb_bands):
            for j in range(len(self._kpoints)):
                if(self._bands[i]['energy'][j]<self._efermi):
                    if(self._bands[i]['energy'][j]>max_tmp):
                        max_tmp=self._bands[i]['energy'][j]
                        index=j
                        kpointvbm=self._kpoints[j]
        #get all other bands sharing the vbm
        list_index_band=[]
        for i in range(self._nb_bands):
            if(math.fabs(self._bands[i]['energy'][index]-max_tmp)<0.001):
                list_index_band.append(i)
        return {'band_index':list_index_band,'kpoint_index':index,'kpoint':kpointvbm['kpoint'],'energy':max_tmp,'label':kpointvbm['label']}
    
    def getCBM(self):
        """
        get the conduction band minimum (CBM). returns a dictionnary with
        'band_index': a list of the indices of the band containing the CBM (please note that you can have several bands 
        sharing the CBM)
        'kpoint_index': the index in self._kpoints of the kpoint cbm
        'kpoint': the kpoint (in cartesian with 2pi factor)
        'energy': the energy of the CBM
        'label': the label of the cbm kpoint if any
        """
        max_tmp=1000.0
        index=None
        for i in range(self._nb_bands):
            for j in range(len(self._kpoints)):
                if(self._bands[i]['energy'][j]>self._efermi):
                    if(self._bands[i]['energy'][j]<max_tmp):
                        max_tmp=self._bands[i]['energy'][j]
                        index=j
                        kpointvbm=self._kpoints[j]
        #get all other bands sharing the vbm
        list_index_band=[]
        for i in range(self._nb_bands):
            if(math.fabs(self._bands[i]['energy'][index]-max_tmp)<0.001):
                list_index_band.append(i)
        return {'band_index':list_index_band,'kpoint_index':index,'kpoint':kpointvbm['kpoint'],'energy':max_tmp,'label':kpointvbm['label']}
    
    
    


    def get_band_gap(self):
        """
        get the band gap 
        returns a dictionary with:
        'energy': the band gap energy in eV
        'direct': a boolean telling if the gap is direct (True) or not (False)
        'transition': the kpoint labels of the transition (e.g., \Gamma to X)
        
        TODO: not sure if the direct works, to test!
        
        """
        if(self.is_metal()==True):
            return {'energy':0.0,'direct':False,'transition':None}
        cbm=self.getCBM()
        vbm=self.getVBM()
        result={}
        result['energy']=cbm['energy']-vbm['energy']
        result['direct']=False
        if(cbm['label']==vbm['label'] or np.linalg.norm(cbm['kpoint']-vbm['kpoint'])<0.01):
            result['direct']=True
        result['transition']=str(vbm['label'])+"-"+str(cbm['label'])
        return result
    
    def get_direct_gap(self):
        """
        get the direct gap
        
        """
        if(self.is_metal()==True):
            return 0.0
        #go through each kpoint and find the vbm and cbm at this kpoint
        #find the minimum of this difference
        diffmin=10000
        for i in range(len(self._kpoints)):
            diff=0.0
            for j in range(self._nb_bands):
                if(self._bands[j]['energy'][i]>self._efermi):
                    cbm=self._bands[j]['energy'][i]
                    vbm=self._bands[j-1]['energy'][i]
                    #print str(i)+" "+str(cbm)+" "+str(vbm)
                    diff=cbm-vbm
                    #print str(diff)+" "+str(diffmin)
                    if(diff<diffmin):
                        diffmin=diff
                    break
        #print diffmin
        return diffmin
        
        
    def get_pg_matrices_rec(self):
        """
            get the point group matrices in the reciprocal lattice of the structure
            calls aconvasp to get this
        """
        return pymatgen.command_line.aconvasp_caller.get_point_group_rec(self._structure)
    
    def kpoints_from_label(self,label):
        """
            get the kpoint corresponding to a given label
        """
        to_return=[]
        count=0
        for c in self._kpoints:
            if(c['label']==label):
                tmp=c['kpoint_index']=count+1
                to_return.append(tmp)
            count=count+1
        return to_return
    
    def plot_kpoints_IBZ(self):
        """
        gives the irreducible brillouin zone and the kpoint in the band structure
        this is a first pass at this
        """
        
        import pylab as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print self._lattice_rec.abc
        ax.scatter([self._kpoints[i]['kpoint'][0] for i in range(len(self._kpoints))],[self._kpoints[i]['kpoint'][1] for i in range(len(self._kpoints))],[self._kpoints[i]['kpoint'][2] for i in range(len(self._kpoints))])
        
        for i in range(len(self._kpoints)):
            if(self._kpoints[i]['label']!=None):
            #ax.text(newkpoint[0],newkpoint[1],newkpoint[2],self._kpoints[i]['label'])
                ax.text(self._kpoints[i]['kpoint'][0],self._kpoints[i]['kpoint'][1],self._kpoints[i]['kpoint'][2],self._kpoints[i]['label'])
                
            
        plt.show()
        
        
    def get_principal_directions(self,kpoint):
        """
        using the point group symmetry figures out what are the principal directions
        for the effective mass tensor at kpoint
        of course, in the general case we can't know all principal directions
        """
        #find the sym. operation keeping ko
        listUc=self.get_pg_matrices_rec()
        list_return=[]
        list_matrix=[]
        for u in listUc:
            #print u
            newkpoint=np.dot(u['matrix'],kpoint.T)
            if(np.linalg.norm(newkpoint-kpoint)<0.001):
                list_matrix.append(u['matrix'])
                already=False
                for a in list_return:
                    if(np.linalg.norm(u['axis']-a)<0.001):
                        already=True
                if(already==False):
                    list_return.append(u['axis'])
            
        """    
        for i in range(len(list_return)):
            for j in range(len(list_return)):
                if(i==j):
                    continue
                fits=False
                for u in list_matrix:
                    newkpoint=np.dot(u,list_return[i].T)
                    kpoint=list_return[j].T
                    if(np.linalg.norm(newkpoint-kpoint)<0.001):    
                        fits=True
         """                   
            
        return list_return
        
        
    
    def get_effective_mass_average(self,target):
        """
        this method averages the effective masses of several branches in one average effective mass tensor
        the average is based on the occupation and heavier bands see a larger occupation
        the weight factor is based on the occupation of parabolic bands
        """
        tensors=self.get_effective_mass_smart(target)
        #fist we computed the mDOS=(mx*my*mz)^1/3 for each band
        mDOS=[]
        tot=0.0
        for t in tensors:
            m=(math.fabs(t[0]*t[1]*t[2]))**(1.0/3.0)
            mDOS.append(m)
            tot=tot+m**(3.0/2.0)
        #then we average with the mDOS^3/2 weights
        result=[0.0,0.0,0.0]
        for i in range(len(tensors)):
            result[0]=result[0]+mDOS[i]**(3.0/2.0)*np.float(tensors[i][0])/tot
            result[1]=result[1]+mDOS[i]**(3.0/2.0)*np.float(tensors[i][1])/tot
            result[2]=result[2]+mDOS[i]**(3.0/2.0)*np.float(tensors[i][2])/tot
        return result
    
    def get_effective_mass(self,target):
        """
        this gives all the effective masses for a given kpoint (even if there are several bands including the point)
        this gives a list of list of eigenvalues for one given kpoint
        """
        list_return=[]
        for c in target['band_index']:
            list_return.append(self.get_effective_mass_one_band(target['kpoint'], target['kpoint_index'], target['label'], c))
        return list_return
    
    def get_effective_mass_one_band(self, kpoint, kpoint_index, label, band_index):
        """
        this gives all the eigenvalues of the effective mass tensor that can be deducted from the band structure and ONE band
        It uses the symmetry of the Brillouin zone to do so
        """
        set_branches_nb=set()
        #look at all branches containing the point with label
        branches_of_interest=[]
        if(label!=None):
            for i in range(len(self._kpoints)):
                if(self._kpoints[i]['label']==label):
                    for j in range(len(self._branches)):
                        d=self._branches[j]
                        if i in d:
                            branches_of_interest.append({'indices':d,'vector':(self._kpoints[d[-1]]['kpoint']-self._kpoints[d[0]]['kpoint']),'index':j})
        else:
                for j in range(len(self._branches)):
                    d=self._branches[j]
                    if kpoint_index in d:
                        branches_of_interest.append({'indices':d,'vector':(self._kpoints[d[-1]]['kpoint']-self._kpoints[d[0]]['kpoint']),'index':j})
        #print branches_of_interest
        
        principal_directions=self.get_principal_directions(kpoint)
        
        #check if any of the branches are colinear with one of the all_eigenvectors
        for c in principal_directions:
            for d in branches_of_interest:
                if(math.fabs(np.dot(d['vector'],c)/(np.linalg.norm(d['vector'])*np.linalg.norm(c))<1.01 and math.fabs(np.dot(d['vector'],c)/(np.linalg.norm(d['vector']))*np.linalg.norm(c)))>0.99):
                    #print c
                    #print d['vector']
                    #print np.dot(d['vector'],c)/(np.linalg.norm(d['vector'])*np.linalg.norm(c))
                    set_branches_nb.add(d['index'])
        #print set_branches_nb       
        mass=[]
        #for the each branch get the mass
        
        for b in set_branches_nb:
            #print str(self._branches[b])+" "+str(self.get_effective_mass_along_line(kpoint_index, band_index, b)[0])
            mass.append(self.get_effective_mass_along_line(kpoint_index, band_index, b)[0])
        
        return mass
        
        #print self._branches
        #for c in self._branches:
            #if(target['kpoint_index'] in c):
                #print c
            
    def get_effective_mass_along_line(self, index_k, index_band, index_branch):
        """
            this is a very simple method giving the effective mass along a specific branch
            These do not have to correspond to one of the eigenvalue of the effective mass tensor
            This is based on a simple fit of a parabola on a few points
            TODO: change the way the points are picked
        """
        import scipy

        nb_sample=4
        mass=[]
        #get all points around (from -4 to 4) index_k on the band=
        
        
        #two cases: there a nb_sample points in the branch with higher indices or lower ones
        forward=False
        backward=False
        local_index=None
        for i in range(len(self._branches[index_branch])):
            if(np.linalg.norm(self._kpoints[self._branches[index_branch][i]]['kpoint']-self._kpoints[index_k]['kpoint'])<0.001):
                local_index=i 
                if(len(self._branches[index_branch])>i+4):
                    forward=True
                if(i-4>0):
                    backward=True
        x0=self._kpoints[self._branches[index_branch][local_index]]['distance']
        y0=self._bands[index_band]['energy'][index_k]
        if(forward):
            x=[self._kpoints[self._branches[index_branch][local_index+i]]['distance']-x0 for i in range(nb_sample)]
            y=[self._bands[index_band]['energy'][self._branches[index_branch][local_index+i]]-y0 for i in range(nb_sample)]
            tck = scipy.interpolate.splrep(x,y)
            yder = scipy.interpolate.splev(0,tck,der=2)
            #import pylab, numpy
            #pylab.plot(x,y,'x')
            #xnew=numpy.arange(0,x[3],0.01)
            #ynew=scipy.interpolate.splev(xnew,tck)
            #pylab.plot(xnew,ynew)
            #print yder/2.0
            #a=np.polyfit(x,y,2)
            #print a
            #pylab.plot(xnew,a[0]*xnew**2+a[1]*xnew+a[2])
            #pylab.show()
            
            mass.append(3.77/(yder/2.0))
            #mass.append(3.77/a[0])
         
        if(backward):
            x=[-1.0*(self._kpoints[self._branches[index_branch][local_index-i]]['distance']-x0) for i in range(nb_sample)]
            y=[self._bands[index_band]['energy'][self._branches[index_branch][local_index-i]]-y0 for i in range(nb_sample)]
            a=np.polyfit(x,y,2)
            tck = scipy.interpolate.splrep(x,y)
            yder = scipy.interpolate.splev(0,tck,der=2)
            mass.append(3.77/(yder/2.0))
            #mass.append(3.77/a[0])
               
        return mass
        
    
        
    def get_stationary_pg_for_kpoints(self,kpoint):
        listUc=self.get_pg_matrices_rec()
        list_return=[]
        for u in listUc:
                #print i
            newkpoint=np.dot(u,kpoint.T)

                
            if(np.linalg.norm(newkpoint-kpoint)<0.001):
                list_return.append(u)
        return list_return
    
    def get_eigen_vectors_for_pg(self):
        listUc=self.get_pg_matrices_rec()
        list_return=[]
        for u in listUc:
            print u
            print np.linalg.eig(u)
            list_return.append(np.linalg.eig(u))
        return list_return
    

    
    def is_metal(self):
        """
        check if the band structure indicates a metal by looking if the fermi level crosses a band
        
        """
        for i in range(self._nb_bands):
            below=False
            above=False
            for j in range(len(self._kpoints)):
                if(self._bands[i]['energy'][j]<self._efermi):
                    below=True
                if(self._bands[i]['energy'][j]>self._efermi):
                    above=True
            if(above==True and below==True):
                return True
        return False
                
    
    def plot_bands_compare(self, otherbands):
        """
        plot two band structure for comparison.
        TODO: still a lot of work to do that nicely!
        """
        import pylab
        for bl in self._branch_labels:
            for i in range(self._nb_bands):
                pylab.plot([self._kpoints[j]['distance'] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],[self._bands[i]['energy'][j] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],'-r')
        
        for bl in otherbands._branch_labels:
            for i in range(otherbands._nb_bands):
                pylab.plot([otherbands._kpoints[j]['distance'] for j in range(len(otherbands._kpoints)) if otherbands._kpoints[j]['branch']==bl],[otherbands._bands[i]['energy'][j] for j in range(len(otherbands._kpoints)) if otherbands._kpoints[j]['branch']==bl],'-k')
        
        ticks=self.get_ticks()
        
        pylab.gca().set_xticks(ticks['distance'])
        pylab.gca().set_xticklabels(ticks['label'])
        pylab.xlabel('Kpoints', fontsize = 'large')
        pylab.ylabel('Energy(eV)', fontsize = 'large')
        #pylab.ylim(vbm-4,cbm+4)
        pylab.show()
        pylab.legend()
        
    def plot_bands(self,band=None,kpoint_index=[]):
        """
        plot the band structure.
        band indicates the index number of the specific band to plot, None plots all bands
        TODO: add more options
        """
        import pylab
        pylab.figure
        if(band==None): 
            for bl in self._branch_labels:
                for i in range(self._nb_bands):
                    pylab.plot([self._kpoints[j]['distance'] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],[self._bands[i]['energy'][j] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],'b-',linewidth=5)
            for bl in self._branch_labels:
                for i in range(self._nb_bands):
                    pylab.plot([self._kpoints[j]['distance'] for j in kpoint_index if self._kpoints[j]['branch']==bl],[self._bands[i]['energy'][j] for j in kpoint_index if self._kpoints[j]['branch']==bl],'r^')
        else:
            for bl in self._branch_labels:
                i=band
                pylab.plot([self._kpoints[j]['distance'] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],[self._bands[i]['energy'][j] for j in range(len(self._kpoints)) if self._kpoints[j]['branch']==bl],'k-o')
            for bl in self._branch_labels:
                pylab.plot([self._kpoints[j]['distance'] for j in kpoint_index if self._kpoints[j]['branch']==bl],[self._bands[i]['energy'][j] for j in kpoint_index if self._kpoints[j]['branch']==bl],'r^')
        
        ticks=self.get_ticks()
        for i in range(len(ticks['label'])):
            if(ticks['label'][i]!=None):
                pylab.axvline(ticks['distance'][i],color='k')
            
        pylab.gca().set_xticks(ticks['distance'])
        pylab.gca().set_xticklabels(ticks['label'])
        pylab.xlabel('Kpoints', fontsize = 'large')
        pylab.ylabel('Energy(eV)', fontsize = 'large')
        vbm=self.getVBM()
        cbm=self.getCBM()
        if(band==None):
            pylab.ylim(vbm['energy']-4,cbm['energy']+4)
        pylab.show()
        pylab.legend()
        
        
    def get_ticks(self):
        """
        get all ticks and labels for a band structure plot
        """
        tick_distance=[]
        tick_labels=[]
        for c in self._kpoints:
            if(c['label']!=None):
                tick_distance.append(c['distance'])
                tick_labels.append(c['label'])
        return {'distance':tick_distance,'label':tick_labels}
    
    def get_difference(self,otherbs):
        diff=[]
        for n in range(self._nb_bands):
            diff.append([self._bands[n]['energy'][i]-otherbs._bands[n]['energy'][i] for i in range(len(self._kpoints))])
        return diff
    
    def to_dict(self):
        dictio={}
        dictio['lattice_rec']=self._lattice_rec.to_dict
        dictio['efermi']=self._efermi
        dictio['kpoints']=self._kpoints
        dictio['branches']=self._branches
        dictio['bands']=self._bands
        dictio['is_metal']=self.is_metal()
        dictio['VBM']=self.getVBM()
        dictio['CBM']=self.getCBM()
        dictio['band_gap']=self.get_band_gap()
        return dictio
    
       
def get_reconstructed_band_structure(list_bs):
    """
        this method takes a list of band structure (divided by branches)
        and reconstruct one band structure object from all of them
        TODO: check how to assign the fermi level
    """
    kpoints=[]
    eigenvals=[]
    labels_dict={}
    rec_lattice=list_bs[0]._lattice_rec
    structure=list_bs[0]._structure
    efermi=list_bs[0]._efermi
    nb_bands=list_bs[0]._nb_bands
    for bs in list_bs:
        for k in bs._kpoints:
            kpoints.append(k['kpoint'])
        for k, v in bs._labels_dict.iteritems():
            labels_dict[k]=v
        #eigenvals.append({'energy':[0,6.0],'occup':[1.0,1.0]})
    for i in range(nb_bands):
        eigenvals.append({'energy':[],'occup':[]})
        for bs in list_bs:
            for e in bs._bands[i]['energy']:
                eigenvals[i]['energy'].append(e)
            for u in bs._bands[i]['occup']:
                eigenvals[i]['occup'].append(u)
    return Bandstructure(kpoints,eigenvals,labels_dict,rec_lattice,structure,efermi) 

class PBandstructure(Bandstructure):
    """
    Projected Bandstructure for a specific orbital. Extends the Bandstructure object.
    For simplicity, I tried to follow the structure of the DOS implementation
    WARNING: I just started this... not really working yet
    """
    def __init__(self,kpoints,eigenvals,labels_dict,rec_lattice,structure, orbital):
        Bandstructure.__init__(self,kpoints,eigenvals,labels_dict,rec_lattice,structure, orbital)
        self.orbital = orbital
       
class CompleteBandStructure(Bandstructure):
    """
    This wrapper class defines a BandStructure, and also provides a list of PBandStructure (one for each site in the structure)
    WARNING: I just started this... not really working yet
    """
    
    def __init__(self, total_bandstructure, pbandstructure, structure):
        """
        """
        self._bandstructure = total_bandstructure
        self._pdos = {structure[i]:{Orbital.from_vasp_index(j) : pbandstructure[i][j] for j in range(len(pbandstructure[i]))} for i in range(structure.num_sites)}
        self._structure = structure        
    
    
    