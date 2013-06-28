#!/usr/bin/env python

"""
Interface with command line GULP. http://projects.ivec.org/gulp/help/manuals.html
WARNING: you need to have GULP in your path for this to work
"""

__author__ = "Bharat Medasani, Wenhao Sun"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "bkmedasani@lbl.gov,wenhao@mit.edu"
__status__ = "Production"
__date__ = "$Jun 22, 2013M$"

import subprocess
import os

from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.command_line.aconvasp_caller import run_aconvasp_command
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.finder import SymmetryFinder


_anions = set(map(Element, ["O","S","F","Cl","Br","N","P"]))
_cations = set(map(Element, [
    "Li","Na","K", # alkali metals
    "Be","Mg","Ca", # alkaline metals
    "Al","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","Au"
    ]))
_gulp_kw = {
        #Control of calculation type
        "angle","bond","cosmo","cosmic","cost","defect","distance", 
        "eem","efg","fit","free_energy","gasteiger","genetic", 
        "gradients","md","montecarlo","noautobond","noenergy","optimise",
        "pot","predict","preserve_Q","property","phonon","qeq","qbond", 
        "single","sm","static_first","torsion","transition_state",
        #Geometric variable specification
        "breathe","bulk_noopt","cellonly","conp","conv","isotropic", 
        "orthorhombic","nobreathe","noflgs","shell","unfix",
        #Algorithm
        "c6","dipole","fbfgs","fix_molecule","full","hill","kfull", 
        "marvinSE","madelung","minimum_image","molecule","molmec","molq", 
        "newda","noanisotropic_2b","nod2sym","nodsymmetry", 
        "noelectrostatics","noexclude","nofcentral","nofirst_point",
        "noksymmetry","nolist_md","nomcediff","nonanal","noquicksearch", 
        "noreal","norecip","norepulsive","nosasinitevery","nosderv", 
        "nozeropt","numerical","qiter","qok","spatial","storevectors", 
        "nomolecularinternalke","voight","zsisa",
        #Optimisation method
        "conjugate","dfp","lbfgs","numdiag","positive","rfo","unit",
        #Output control
        "average","broaden_dos","cartesian","compare","conserved",
        "dcharge","dynamical_matrix",
        "eigenvectors","global","hessian","hexagonal","intensity","linmin",
        "meanke","nodensity_out","nodpsym","nofirst_point","nofrequency",
        "nokpoints","operators","outcon","prt_eam","prt_two",
        "prt_regi_before","qsas","restore","save","terse",
        #Structure control
        "full", "hexagonal", "lower_symmetry", "nosymmetry",
        #PDF control
        "PDF","PDFcut","PDFbelow","PDFkeep","coreinfo","nowidth","nopartial",
        #Miscellaneous
        "nomodcoord","oldunits","zero_potential"
        }
    
class GulpCaller():
    """
    Class that defines the methods to generate gulp input, run gulp
    and process the output
    """
    def __init__(self,cmd='gulp'):
        """
        Initialize with the executable if not in the standard path
        """
        def is_exe(f):
            return os.path.isfile(f) and os.access(f,os.X_OK)

        fpath, fname = os.path.split(cmd)
        if fpath:
            if is_exe(cmd):
                self._gulp_cmd = cmd
                return
        else:
            for path in os.environ['PATH'].split(os.pathsep):
                path = path.strip('"')
                file = os.path.join(path,cmd)
                if is_exe(file):
                    self._gulp_cmd = file
                    return
        raise GulpError("Executable not found")

    def gulp_kw_line(*args):
        """
        Checks if the input args are proper gulp keywords and
        generates the 1st line of gulp input. 
        """
        if len(list(filter(lambda x: x in _gulp_kw, args))) != len(args):
            raise GulpError("Wrong keywords given")
        gin = " ".join(args)
        gin += "\n"
        return gin

    def gulp_structure(structure, cell_flg=True, frac_flg=True, 
            anion_shell_flg=True, cation_shell_flg=False, symm_flg=True):
        """
        Generates GULP input string for pymatgen structure 
        Args:
            structure:
                pymatgen Structure object
            cell_flg (default = True):
                If true, unit cell is written. 
            fractional_flg (default = True):
                If True, fractional coordinates are used.
                Else, cartesian coodinates in Angstroms are used. 
            ******
            GULP Convention is to use fractional coordinates for periodic 
            structures and cartesian coordinates for non-periodic structures.
            ******
            anion_shell_flg (default = True):
                If True, anions are considered polarizable.
            cation_shell_flg (default = False):
                If True, cations are considered polarizable.
            symm_flg (default = True):
                If True, symmetry information is also written.
        Returns:
            string containing structure for gulp input
        """
        if cell_flg:
            gin = "cell\n"
            l = structure.lattice
            lat_str = map(str, [l.a,l.b,l.c,l.alpha,l.beta,l.gamma])
            gin += " ".join(lat_str)+"\n"

        if frac_flg:
            gin += "frac\n"
            coord_attr = "frac_coords"
        else:
            gin += "cart\n"
            coord_attr = "coords"
        for site in structure.sites: 
            coord = map(str, list(getattr(site, coord_attr)))
            specie = site.specie
            core_site_desc = specie.symbol+" core "+" ".join(coord)+"\n"
            gin += core_site_desc
            if ((specie in _anions and anion_shell_flg) or 
                (specie in _cations and cation_shell_flg)):
                shel_site_desc = specie.symbol+" shel "+" ".join(coord)+"\n"
                gin += shel_site_desc
            else:
                pass

        if (symm_flg):
            gin += "space\n"
            gin += str(SymmetryFinder(structure).get_spacegroup_number())+"\n"
        return gin

    def gulp_specie_potential(structure, potential, **kwargs):
        """
        Generates GULP input specie and potential string for pymatgen structure
        Args:
            structure:
                pymatgen Structure object
            potential:
                String specifying the type of potential used
            kwargs:
                Additional parameters related to potential
                For potential == "buckingham":
                    anion_shell_flg (default = False):
                        If True, anions are considered polarizable.
                        anion_core_chrg=float
                        anion_shell_chrg=float
                    cation_shell_flg (default = False):
                        If True, cations are considered polarizable.
                        cation_core_chrg=float:
                        cation_shell_chrg=float:
        Returns:
            string containing specie and potential specification for gulp input
        """
        raise GulpError("Function gulp_specie_potential not yet implemented"+
                "\nUse gulp_lib instead")

    def gulp_lib(file_name):
        """
        Specifies GULP library file to read species and potential parameters.
        If using library don't specify species and potential 
        in the input file and vice versa.
        Args:
            file_name:
                Name of GULP library file
        Returns:
            GULP input string specifying library keyword
        """
        gulplib_set = lambda: 'GULP_LIB' in os.environ.keys()
        readable = lambda f: os.path.isfile(f) and os.access(f, os.R_OK)

        dirpath, fname = os.path.split(file_name)
        if dirpath:       #Full path specified
            if readable(file_name):
                gin = 'library '+file_name
            else:
                raise GulpError('GULP Library not found')
        else:
            fpath = os.path.join(os.getcwd(), file_name)  #Check current dir
            if readable(fpath):
                gin = 'library '+fpath
            elif gulplib_set():
                fpath = os.path.join(os.environ['GULP_LIB'],file_name)
                if readable(fpath):
                    gin = 'library '+file_name
                else:
                    raise GulpError('GULP Library not found')
        gin += "\n"
        return gin

    def run_gulp_input(ginfile):
        """
        Given a complete gulp.in file, run GULP
        """
        command=["gulp"]
        p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stdin=subprocess.PIPE,
                         stderr=subprocess.PIPE)
        output = p.communicate(ginfile)
        if "Error" in output[1] or "error" in output[1]:
            print ginfile
            print "----output_0---------"
            print output[0]
            print "----End of output_0------"
            print "----output_1--------"
            print output[1]
            print "----End of output_1------"
            raise GulpError(output[1])

        # We may not need this
        if "ERROR" in output[0]:
            raise GulpError(output[0])

        # Mott_Littleton method may fail to reach proper convergence
        conv_err_string = "Conditions for a minimum have not been satisfied"
        if conv_err_string in output[0]:
            raise GulpConvergenceError()

        gout_string = ""
        for line in output[0].split("\n"):
            gout_string = gout_string + line + "\n"
        return gout_string

    def gulp_energy(gulpout):
        energy = None
        for line in gulpout.split("\n"):
            if "Total lattice energy" in line and "eV" in line:
                energy = line.split()
        if energy:
            return energy[4]
        else:
            #print gulpout
            raise GulpError("Energy not found in Gulp output")

    def binaryoxide_tersoff_gulpinput(structure):
        '''
        Gets a GULP input for an oxide structure
        CURRENTLY ONLY WORKS FOR BINARY OXIDES WITH A SINGLE OXIDATION STATE
        Calculates oxidation state from the formula
        '''
        

        gin = self.gulp_kw_line("static", "noelectrostatics")
        gin += self.gulp_structure(structure, False, False, False, False)

        sites=structure.sites
        comp=structure.composition.get_el_amt_dict()

        ginput = gin
        ginput = ginput + "species \n"
        
        lastspecie=[]
        metaloxi=[]
        endstring=""
        for ii in sites:
            if lastspecie != ii.specie:
                lastspecie=ii.specie
                
                '''ONE DAY UPDATE THIS WITH BVANALYZER'''
                '''The only challenge is, I don't know how to deal with cation-cation potential'''
                
                if str(lastspecie) != "O":
                    nummet=comp[str(lastspecie)]
                    numoxi=comp["O"]
                    oxistate=numoxi*2/nummet
                    if oxistate%1 != 0:
                        raise SystemError("Oxide has mixed valence on metal")
                    oxidationstring=str(lastspecie)+" core "+str(oxistate)
                    metaloxi.append(lastspecie)
                else:
                    oxidationstring=str(lastspecie)+" core -2"
                    
                ginput=ginput+oxidationstring+ "\n"
                endstring=endstring+"qerfc \n"+str(lastspecie)+" "+str(lastspecie)+"  0.6000 10.0000 \n" 
        
        ginput=ginput+"# noelectrostatics \n Morse \n"
        
        for metal in metaloxi:
            metal=str(metal)+"("+str(int(oxistate))+")"
            MetOxiTers=Tersoff_pot().data[metal]
            ginput=ginput+MetOxiTers
        
        ginput=ginput+endstring
        
        return ginput

    def get_binoxi_gulp_energy_tersoff(structure):
        output=binaryoxide_tersoff_gulpinput(structure)
        output2=run_gulp_input(output)
        return float(_gulp_energy(output2))

def binaryoxide_tersoff_gulpinput(structure):
    '''
    Gets a GULP input for an oxide structure
    CURRENTLY ONLY WORKS FOR BINARY OXIDES WITH A SINGLE OXIDATION STATE
    Calculates oxidation state from the formula
    '''
    
    #gin=get_gulpinput(structure)
    gin=_gulp_structure(structure, False, False, False, False)
    gin="static noelectrostatics \n "+gin
    #print gin
    specs=structure.sites
    
    comp=structure.composition.get_el_amt_dict()

    ginput = gin
    #    ginput=""
    #     c=0
    #    for line in gin.split("\n"):
    #        c=c+1
    #        if c != 2 and c != 3 and c != 8:
    #            if c==1:
    #                ginput = ginput+line   
    #            elif c>=10 and c < 10+structure.num_sites:
    #                d=c-10
    #                ginput = ginput + "\n" + str(specs[d].specie) +" core "+line
    #            else:
    #                ginput = ginput+"\n"+line 
    ginput = ginput + "species \n"
    
    lastspecie=[]
    metaloxi=[]
    endstring=""
    for ii in specs:
        if lastspecie != ii.specie:
            lastspecie=ii.specie
            
            '''ONE DAY UPDATE THIS WITH BVANALYZER'''
            '''The only challenge is, I don't know how to deal with cation-cation potential'''
            
            if str(lastspecie) != "O":
                nummet=comp[str(lastspecie)]
                numoxi=comp["O"]
                oxistate=numoxi*2/nummet
                if oxistate%1 != 0:
                    raise SystemError("Oxide has mixed valence on metal")
                oxidationstring=str(lastspecie)+" core "+str(oxistate)
                metaloxi.append(lastspecie)
            else:
                oxidationstring=str(lastspecie)+" core -2"
                
            ginput=ginput+oxidationstring+ "\n"
            endstring=endstring+"qerfc \n"+str(lastspecie)+" "+str(lastspecie)+"  0.6000 10.0000 \n" 
    
    ginput=ginput+"# noelectrostatics \n Morse \n"
    
    for metal in metaloxi:
        metal=str(metal)+"("+str(int(oxistate))+")"
        MetOxiTers=Tersoff_pot().data[metal]
        ginput=ginput+MetOxiTers
    
    ginput=ginput+endstring
    
    return ginput

def binaryoxide_buckingham_gulpinput(structure, library):
    '''
    Gets a GULP input for an oxide structure and 
    buckingham potential from library
    '''
    
    #gin=get_gulpinput(structure)
    gin=_gulp_structure(structure, True, False, False, False)
    gin="single \n "+gin
    gin += _gulp_lib('lewis.lib')
    return gin
    #print gin

def gulpduplicatecheck(structure):
    '''
    Gets a GULP input for any structure
    Uses only an Au potential for all atoms.
    Not to actually get energies - just used to attribute some arbitrary energy to the structure.
    Identical structures will have the same 'energy', which can be used to pre-screen structurematcher for large structures. 
    
    '''
    
    #gin=get_gulpinput(structure)
    #print structure
    #print '------'
    gin=_gulp_structure(structure,False,False,False,False)
    gin="single \n "+gin
    #print gin
    #print '------'
    
    comp=structure.composition.get_el_amt_dict()
    #print comp

    ginput=""
    c=0
    for line in gin.split("\n"):
        c=c+1
        if c != 2 and c != 3:
            if c==1 or c==4:
                ginput += line+"\n"
            elif c>=5 and c < 5+structure.num_sites:
                splt_line = line.split(" ")
                splt_line[0] = 'Au'
                line = " ".join(splt_line)
                ginput += line+"\n"
            else:
                ginput += line + "\n"
    ginput = ginput + """lennard 12 6
Au core Au core 214180.2000 625.482 40.000 0 0"""
    #print ginput
    #print '-------'
    output2=run_gulp_input(ginput)
    return float(_gulp_energy(output2))

def get_binoxi_gulp_energy_tersoff(structure):
    output=binaryoxide_tersoff_gulpinput(structure)
    output2=run_gulp_input(output)
    return float(_gulp_energy(output2))

def get_binoxi_gulp_energy_buckingham(structure):
    gulp_inp=binaryoxide_buckingham_gulpinput(structure,'lewis.lib')
    output2=run_gulp_input(gulp_inp)
    return float(_gulp_energy(output2))

class GulpError(Exception):
    """
    Exception class for GULP.
    Raised when the GULP gives an error
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return "GulpError : " + self.msg

class GulpConvergenceError(Exception):
    """
    Exception class for GULP.
    Raised when proper convergence is not reached in Mott-Littleton
    defect energy optimisation procedure in GULP
    """
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return "GulpConvergenceError: " + self.msg

class Tersoff_pot(object):
    def __init__(self, verbose=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fid = open(os.path.join(module_dir, "OxideTersoffPotentials"), "rU")
        data = dict()
        for row in fid:
            metaloxi=row.split()[0]
            line=row.split(")")
            data[metaloxi]=line[1]
        fid.close()
        self.data=data
        
def get_gulpinput(structure):
    """
    From a structure, get a starting gulp.in file using aconvasp
    (Deprecated)
    """
    output=run_aconvasp_command(["aconvasp", "--gulp"], structure)
    gin_string = ""
    for line in output[0].split("\n"):
        if gin_string=="":
            gin_string = gin_string + line
        else:
            gin_string = gin_string +  "\n" +line
    
    return gin_string

