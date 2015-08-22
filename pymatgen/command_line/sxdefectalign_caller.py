# coding: utf-8
from __future__ import division, unicode_literals

"""
TO DO: include compatibility with all the extra flags that sxdefectalign has?

This module implements an interface to the Freysoldt et al.'s excellent
code for calculating a correction for formation energy of a charged 
defect.

This module depends on a compiled sxdefectalign executable available in 
the path. Please download the executable and manual at 
http://sxlib.mpie.de/wiki/AddOns.

If you use this module, please cite the following:

Christoph Freysoldt, Jörg Neugebauer, and Chris G. Van de Walle, 
Phys. Rev. Lett. 102, 016402, 2009.
"""

__author__ = "Danny Broberg, Geoffroy Hautier, Bharat Medasani"
__version__ = "0.1"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com"
__status__ = "Beta"
__date__ = "8/5/15"

import os
import numpy as np
import subprocess
import shutil

from pymatgen.io.vaspio import Locpot
from monty.os.path import which
from monty.dev import requires
from monty.tempfile import ScratchDir

@requires(which("sxdefectalign"),
          "Charge defect correction requires the executable sxdefectalign to "
          "be in the path. Please download the executable at "
          "http://sxlib.mpie.de/wiki/AddOns.")
class FreysoldtCorrection(object):
    """
    Charge correction wraparound for sxdefectalign code written by Freysoldt

    """
    def __init__(self, locpot_ref, locpot_def, charge, epsilon, encut,
                 site_frac_coords, align=[0,0,0], lengths=None, name=''):
        """
        Initialize the sxdefectalign caller

        Args:
            locpot_ref: The path of the LOCPOT_ref (pure LOCPOT)
            locpot_def: The path of the LOCPOT_def (defect LOCPOT)
                note: these must have been pre-processed so that the 6th line
                of the VASP output for LOCPOT file is deleted (required for
                sxdefectalign code compatibility)
            charge: charge difference with respect to neutral defect cell
                (adding three electrons = charge of -3 )
            epsilon: dielectric constant to use for calculation
            encut: energy cutoff for Vasp (eV)
            site_frac_coords: fractional co-ordinates for defect site
            align: alignment constants for the three axes
                potential alignment from planar averages
            lengths: Optional: lengths of two lattice parameters
                note: these are optional because they can be figured out from
                the given Locpot's but this is time consuming to load these files

        """
	if os.path.exists(locpot_ref):
        	self._locpotref=os.path.abspath(str(locpot_ref))
	else:
		print 'Could not find Locpot_vref in specified path ' \
			'Double check path input for locpot_ref'
	if os.path.exists(locpot_def):
        	self._locpotdef=os.path.abspath(str(locpot_def))
	else:
		print 'Could not find Locpot_vdef in specified path ' \
			'Double check path input for locpot_ref'
        self._charge = charge
        self._epsilon = epsilon
        self._encut = encut
        self._frac_coords = site_frac_coords
        self._align=align
        if not lengths:
            struct=Locpot.from_file(locpot_ref) #maybe include an exception for if the Locpot can't be read
            self._lengths=struct.structure.lattice.abc
            print 'had to import lengths, if you want to speed up class ' \
                  'instantiation set lengths to: '+str(self._lengths)+'\n'
        else:
            self._lengths = lengths
	self.name=name
        self.PCen=None
        self.errors={}
        #check to see if sxdefectalign in path
        command=['sxdefectalign','-v']
        valid=['+-----------------------------------------------------------------------------\n',
               '| S/PHI/nX DFT package by S. Boeck, J. Neugebauer et al.\n', '| S/PHI/nX release 1.0\n',
               '| sxdefectalign by C. Freysoldt\n', '| version 1.3\n',
               '+-----------------------------------------------------------------------------\n']
        with ScratchDir('.'):
            #in case NERSC (Hopper) has issues with python subprocess can use hack
            #p = subprocess.Popen(command, stdout=subprocess.PIPE,
            #        stdin=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True)
            #output, err = p.communicate()
            #out = out.decode('utf-8')
            #err = err.decode('utf-8')

            #this is hack wrap-around for when subprocess doesn't work
            cmd = ' '.join(command)
            cmd += ' > tmptest'
            os.system(cmd)
            with open('tmptest') as f:
                out = f.readlines()
            if out!=valid:
                print 'Could not find sxdefectalign version 1.3 code! ' \
                      'Please make sure the correct version is in your path\n'
                self.errors['code']=True
            else:
                print 'Ready to run sxdefectalign code\n'
                self.errors['code']=None

    def runsxdefect(self):
        if not self._charge:
            print 'charge=0, so no charge correction needed'+'\n'
            self.PCen=[0,0,0]  #triplet for PC correction with respect to each axis
            self.potalign=[0,0,0]
            return

        with ScratchDir('.'):
            result=[]
            platy=[]
            relpos=",".join(str(i) for i in self._frac_coords)

            for axis in [0,1,2]:
                print 'do axis'+str(axis+1)
                command = ['~/sxdefectalign', '--vasp', '-a'+str(axis+1),
                    '--relative', '--pos', relpos,
                    '--charge', str(-self._charge),
                    '--ecut', str(self._encut/13.6057), #eV to Ry for sxdefect
                    '--eps', str(self._epsilon),
                    '-C', str(-float(self._align[axis])),
                    '--vref', self._locpotref,
                    '--vdef', self._locpotdef]
                print str(command)+'\n'

                ##standard way of running NERSC commands.
                #in case NERSC (Hopper) has issues with python subprocess can use hack
                #p = subprocess.Popen(command, stdout=subprocess.PIPE,
                #        stdin=subprocess.PIPE, close_fds=True)
                #out, err = p.communicate()
                #out = out.decode('utf-8')
                #err = err.decode('utf-8')
                #print 'output from sxdefectalign = ', str(out)
                #val=(float(out[0].split("\n")[12].split()[3])) #THIS MAY BE WRONG?

                #this is hack wrap-around for when subprocess doesn't work
                cmd = ' '.join(command)
                cmd += ' > tmpoutput'
                os.system(cmd)
                with open('tmpoutput') as f:
                    output = f.readlines()
                print 'output from sxdefectalign = '+str(output)+'\n'
                val =  output[-1].split()[3].strip()

                #return to normal code
                result.append(float(val))
                print "chg correction is "+str(result[-1])+'\n'
                os.remove('tmpoutput')

                #for alignment and creating easy to plot files
                x_lr, y_lr = [], []
                x, y = [], []
                x_diff, y_diff = [], []
                with open("vline-eV.dat",'r') as f_sr: #read in potential
                    for r in f_sr:
                        tmp = r.split("\t")
                        if(len(tmp)<3 and not r.startswith("&")):
                           x_lr.append(float(tmp[0])/1.889725989)   # to Angstrom
                           y_lr.append(float(tmp[1]))
                        if len(tmp) > 2:
                            x.append(float(tmp[0])/1.889725989)     # to Angstrom
                            x_diff.append(float(tmp[0])/1.889725989)# to Angstrom
                            y.append(float(tmp[2].rstrip("\n")))
                            y_diff.append(float(tmp[1]))

                #Do I have to move this out of the scratch dir in a different way?
                #keep output files for future plotting usage
		if not self.name:
                	os.rename("vline-eV.dat","../axis"+str(axis)+"vline-eV.dat")
		else:
                	os.rename("vline-eV.dat","../"+str(self.name)+"axis"+str(axis)+"vline-eV.dat")


                # Extract potential alignment term averaging window of +/- 1 Ang
                # around point halfway between neighboring defects
                latt_len = self._lengths[axis]
                if self._frac_coords[axis] >= 0.5:
                     platx = (self._frac_coords[axis]-0.5) * latt_len
                else:
                     platx = (self._frac_coords[axis]+0.5) * latt_len
                print "half way between defects is: ", platx

                xmin = latt_len - (1-platx) if platx < 1 else platx-1
                xmax = 1-(latt_len-platx) if platx > latt_len-1 else 1+platx
                print 'means sampling region is (', xmin, ',', xmax, ')'+'\n'

                tmpalign=[]
                if xmax < xmin:
                    print 'wrap around detected, special alignment needed'
                    for i in range(len(x)):
                        if x[i]<xmax or x[i]>xmin:
                            tmpalign.append(y[i])
                        else:
                            continue
                else:
                    for i in range(len(x)):
                        if (x[i]>xmin and x[i]<xmax):
                            tmpalign.append(y[i])
                        else:
                            continue

                print 'alignment is ', -np.mean(tmpalign),'\n'
                platy.append(-np.mean(tmpalign))
                flag = 0
                for i in tmpalign:  #check to see if alignment region varies too much
                    if np.abs(i-platy[-1])>0.2:
                        flag = 1
                    else:
                        continue
                if flag != 0:
                    print 'Warning: potential aligned region varied by more ' + \
                          'than 0.2eV (in range of halfway between defects ' + \
                          '+/-1 \Angstrom). Might have issues with Freidel ' + \
                          'oscillations or atomic relaxation\n'
                    self.errors['alignment'+str(axis)]='large osc (>0.2) in potential,' \
                                            ' check plots for further information'
        self.PCen=result    #triplet for PC correction with respect to each axis
        self.potalign=platy #triplet of potential alignment for each axis


    def get_full(self):
        """
        Output triplet of PC energies and a triplet of potential alignments for each axis
        """
        if not self.PCen:
            self.runsxdefect()
        return [self.PCen,self.potalign]


    def get_avgpotalign(self):
        """
        Gives average of potential alignment corrections for each axis
        """
        if not self.PCen:
            self.runsxdefect()
        avg=np.mean(self.potalign)
        print 'For potential alignment, average of '+str(self.potalign)+' is '+str(avg)+'\n'
        return avg

    def get_avgPCen(self):
        """
        Gives average of PC energies for each axis
        """
        if not self.PCen:
            self.runsxdefect()
        avg=np.mean(self.PCen)
        print 'For PC energy, average of '+str(self.PCen)+' is '+str(avg)+'\n'
        return avg

    def get_singlecorrection(self):
        """
        If one has no intention of plotting potentials and checking for localization (not recommended)
        This returns a single correction after one run of the sxdefectalign code.
        """
        if not self.PCen:
            self.runsxdefect()
        PC=self.get_avgPCen()
        pot=self.get_avgpotalign()
        correction=PC-float(self._charge)*pot
        print 'yields Freysoldt (sxdefectalign) correction energy of '+str(correction)+'\n'
        print 'note that the following error dictionary was created: errors='+str(self.errors)
        return correction



