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
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.analysis.bond_valence import BVAnalyzer


_anions = set(map(Element, ["O","S","F","Cl","Br","N","P"]))
_cations = set(map(Element, [
    "Li","Na","K", # alkali metals
    "Be","Mg","Ca", # alkaline metals
    "Al","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ge","As",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi",
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"
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
    
class GulpIO:
    """
    Class that defines the methods to generate gulp input and process output
    """
    def keyword_line(self, *args):
        """
        Checks if the input args are proper gulp keywords and
        generates the 1st line of gulp input. 
        """
        if len(list(filter(lambda x: x in _gulp_kw, args))) != len(args):
            raise GulpError("Wrong keywords given")
        gin = " ".join(args)
        gin += "\n"
        return gin

    def structure_lines(self, structure, cell_flg=True, frac_flg=True, 
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
        gin = ""
        if cell_flg:
            gin += "cell\n"
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

    def specie_potential_lines(structure, potential, **kwargs):
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
        raise NotImplementedError("Function gulp_specie_potential not yet implemented"+
                "\nUse gulp_lib instead")

    def library_line(self, file_name):
        """
        Specifies GULP library file to read species and potential parameters.
        If using library don't specify species and potential 
        in the input file and vice versa. Make sure the elements of 
        structure are in the library file.
        Args:
            file_name:
                Name of GULP library file
        Returns:
            GULP input string specifying library option
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

    def buckingham_input(self, structure, keywords, library=None, 
            uc=True, valence_dict=None):
        '''
        Gets a GULP input for an oxide structure and 
        buckingham potential from library
        '''
        gin = self.keyword_line(*keywords)
        gin += self.structure_lines(structure, symm_flg=not uc)
        if not library:
            gin += self.buckingham_potential(structure, valence_dict)
        else:
            gin += self.library_line(library)
        return gin

    def buckingham_potential(self, structure, val_dict=None):
        '''
        Generate species, buckingham, and spring options for an oxide structure
        using the parameters in default libraries
        Ref: 1) G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys., 
                18, 1149-1161 (1985)
             2) T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,  
                J. Mater Chem., 4, 831-837 (1994)
        Args:
            structure:
                pymatgen.core.structure.Structure
            val_dict (Needed if structure is not charge neutral)
                El:valence dictionary, where El is element. 
        '''
        if not val_dict:
            bv = BVAnalyzer()
            el = [site.species_string for site in structure.sites]
            valences = bv.get_valences(structure)
            val_dict = dict(zip(el, valences))

        #Try bush library first
        bpb = BuckinghamPotBush()
        bpl = BuckinghamPotLewis()
        gin = ""
        for key in val_dict.keys():
            use_bush = True
            #if key != "O" and val_dict[key]%1 != 0:
            #    raise GulpError("Oxide has mixed valence on metal")
            if key not in bpb.species_dict.keys():
                use_bush = False
            elif val_dict[key] != bpb.species_dict[key]['oxi']:
                use_bush = False
            if use_bush:
                gin += "species \n"
                gin += bpb.species_dict[key]['inp_str']
                gin += "buckingham \n"
                gin += bpb.pot_dict[key]
                gin += "spring \n"
                gin += bpb.spring_dict[key]
                continue

            #Try lewis library next if element is not in bush
            #use_lewis = True
            if key != "O":  # For metals the key is "Metal_OxiState+"
                k = key+'_'+str(int(val_dict[key]))+'+'
                if k not in bpl.species_dict.keys():
                    #use_lewis = False
                    raise GulpError("Element {} not in library".format(k))
                gin += "species\n"
                gin += bpl.species_dict[k]
                gin += "buckingham\n"
                gin += bpl.pot_dict[k]
            else:
                gin += "species\n"
                k = "O_core"
                gin += bpl.species_dict[k]
                k = "O_shel"
                gin += bpl.species_dict[k]
                gin += "buckingham\n"
                gin += bpl.pot_dict[key]
                gin += 'spring\n'
                gin += bpl.spring_dict[key]
        return gin

    def tersoff_input(self, structure, periodic=False, uc=True, 
                                  *keywords):
        '''
        Gets a GULP input with Tersoff potential for an oxide structure 
        '''
        #gin="static noelectrostatics \n "
        gin = self.keyword_line(*keywords)
        gin += self.structure_lines(
                structure, cell_flg=periodic, frac_flg=periodic, 
                anion_shell_flg=False, cation_shell_flg=False, symm_flg=not uc
                )
        gin += self.tersoff_potential(structure)
        return gin

    def tersoff_potential(self, structure):
        '''
        Generate the species, tersoff potential lines for an oxide structure
        '''
        bv = BVAnalyzer()
        el = [site.species_string for site in structure.sites]
        valences = bv.get_valences(structure)
        el_val_dict = dict(zip(el, valences))

        gin = "species \n"
        qerfstring = "qerfc\n"
        
        for key in el_val_dict.keys():
            if key != "O" and el_val_dict[key]%1 != 0:
                raise SystemError("Oxide has mixed valence on metal")
            specie_string = key + " core " + str(el_val_dict[key]) + "\n"
            gin += specie_string
            qerfstring += key + " " + key + " 0.6000 10.0000 \n"
        
        gin += "# noelectrostatics \n Morse \n"
        met_oxi_ters = Tersoff_pot().data
        for key in el_val_dict.keys():
            if key != "O":
                metal = key +"(" + str(int(el_val_dict[key])) + ")"
                ters_pot_str = met_oxi_ters[metal]
                gin += ters_pot_str
        
        gin += qerfstring
        return gin

    def tersoff_potential_without_bv(self, structure):
        '''
        Gets a GULP Tersoff potential lines for an oxide structure
        CURRENTLY ONLY WORKS FOR BINARY OXIDES WITH A SINGLE OXIDATION STATE
        Calculates oxidation state from the formula
        ***Deprecated***
        '''
        sites=structure.sites
        comp=structure.composition.get_el_amt_dict()

        gin = "species \n"
        
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
                    
                gin=gin+oxidationstring+ "\n"
                endstring=endstring+"qerfc \n"+str(lastspecie)+" "+str(lastspecie)+"  0.6000 10.0000 \n" 
        
        gin=gin+"# noelectrostatics \n Morse \n"
        for metal in metaloxi:
            metal=str(metal)+"("+str(int(oxistate))+")"
            MetOxiTers=Tersoff_pot().data[metal]
            gin=gin+MetOxiTers
        gin=gin+endstring
        return gin

    def get_energy(self, gout):
        energy = None
        for line in gout.split("\n"):
            if "Total lattice energy" in line and "eV" in line:
                energy = line.split()
            elif "Non-primitive unit cell" in line and "eV" in line:
                energy = line.split()
        if energy:
            return float(energy[4])
        else:
            #print gout
            raise GulpError("Energy not found in Gulp output")

    def get_relaxed_structure(self, gout):
        #Find the structure lines
        structure_lines = []
        cell_param_lines = []
        #print gout
        output_lines = gout.split("\n")
        no_lines = len(output_lines)
        i = 0
        while i < no_lines:
            line = output_lines[i]
            if "Final fractional coordinates of atoms" in line:
                # read the site coordinates in the following lines
                i += 6
                line = output_lines[i]
                while line[0:2] != '--':
                    structure_lines.append(line)
                    i += 1
                    line = output_lines[i]
                # read the cell parameters 
                i += 12
                for del_i in range(6):
                    line = output_lines[i+del_i]
                    cell_param_lines.append(line)

                break
            else:
                i += 1

        #Process the structure lines
        if structure_lines:
            sp = []
            coords = []
            for line in structure_lines:
                fields = line.split()
                if fields[2] == 'c':
                    sp.append(fields[1])
                    coords.append(list(float(x) for x in fields[3:6]))
        else:
            raise IOError("No structure found")

        if cell_param_lines:
            a = float(cell_param_lines[0].split()[1])
            b = float(cell_param_lines[1].split()[1])
            c = float(cell_param_lines[2].split()[1])
            alpha = float(cell_param_lines[3].split()[1])
            beta = float(cell_param_lines[4].split()[1])
            gamma = float(cell_param_lines[5].split()[1])
            latt = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        else:
            raise IOError("No structure found")

        return Structure(latt, sp, coords)


class GulpCaller:
    """
    Class to run gulp from commandline
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

    def run(self, gin):
        """
        Run GULP using the gin as input
        Args:
            gin:
                GULP input string
        Returns:
            gout:
                GULP output string
        """
        #command=["gulp"]
        p = subprocess.Popen(
                self._gulp_cmd, stdout=subprocess.PIPE,
                stdin=subprocess.PIPE, stderr=subprocess.PIPE
                )
        output = p.communicate(gin)

        if "Error" in output[1] or "error" in output[1]:
            print gin
            print "----output_0---------"
            print output[0]
            print "----End of output_0------\n\n\n"
            print "----output_1--------"
            print output[1]
            print "----End of output_1------"
            raise GulpError(output[1])

        # We may not need this
        if "ERROR" in output[0]:
            raise GulpError(output[0])

        # Sometimes optimisation may fail to reach convergence
        conv_err_string = "Conditions for a minimum have not been satisfied"
        if conv_err_string in output[0]:
            raise GulpConvergenceError()

        gout = ""
        for line in output[0].split("\n"):
            gout = gout + line + "\n"
        return gout


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

def get_energy_tersoff(structure, gulp_cmd='gulp'):
    gio = GulpIO()
    gc = GulpCaller(gulp_cmd)
    gin = gio.tersoff_input(structure)
    #print gin
    gout = gc.run(gin)
    return gio.get_energy(gout)

def get_energy_buckingham(structure, gulp_cmd='gulp', 
        keywords=('optimise', 'conp'), valence_dict=None):
    gio = GulpIO()
    gc = GulpCaller(gulp_cmd)
    gin = gio.buckingham_input(
            structure, keywords, valence_dict=valence_dict
            )
    gout = gc.run(gin)
    #print gout
    return gio.get_energy(gout)

def get_energy_relax_structure_buckingham(structure, 
        gulp_cmd='gulp', keywords=('optimise', 'conp'), valence_dict=None):
    gio = GulpIO()
    gc = GulpCaller(gulp_cmd)
    gin = gio.buckingham_input(
            structure, keywords, valence_dict=valence_dict
            )
    gout = gc.run(gin)
    #print gout
    energy =  gio.get_energy(gout)
    #sp, coords = gio.get_relaxed_structure(gout)
    #print sp
    #print coords
    #print len(sp)
    #print len(coords)
    #relax_structure = Structure(structure.lattice, sp, coords, to_unit_cell=True)
    relax_structure = gio.get_relaxed_structure(gout) 
    return energy, relax_structure

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
    def __init__(self, msg=""):
        self.msg = msg
    def __str__(self):
        return self.msg


class BuckinghamPotLewis(object):
    """
    Generate the Buckingham Potential Table based on "Lewis & Catlow" and
    the metal oxidation state
    Ref: 1) G.V. Lewis and C.R.A. Catlow, J. Phys. C: Solid State Phys., 18,
            1149-1161 (1985)
    """
    def __init__(self, verbose=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fid = open(os.path.join(module_dir, 'lewis.lib'), 'rU')
        # In lewis.lib there is no shell for cation
        species_dict, pot_dict, spring_dict  = {}, {}, {}
        sp_flg, pot_flg, spring_flg = False, False, False
        for row in fid:
            if row[0] == "#":
                continue
            if row.split()[0] == "species":
                sp_flg, pot_flg, spring_flg = True, False, False
                continue
            if row.split()[0] == "buckingham":
                sp_flg, pot_flg, spring_flg = False, True, False
                continue
            if row.split()[0] == "spring":
                sp_flg, pot_flg, spring_flg = False, False, True
                continue

            metaloxi = row.split()[0]
            if sp_flg:
                if metaloxi == "O":
                    if row.split()[1] == "core":
                        species_dict["O_core"] = row
                        continue
                    if row.split()[1] == "shel":
                        species_dict["O_shel"] = row
                        continue
                metal = metaloxi.split('_')[0]
                #oxi_state = metaloxi.split('_')[1][0]
                species_dict[metaloxi] = metal + " core " + row.split()[2]+"\n"
                continue

            if pot_flg:
                if metaloxi == "O":
                    pot_dict["O"] = row
                metal = metaloxi.split('_')[0]
                #oxi_state = metaloxi.split('_')[1][0]
                pot_dict[metaloxi] = metal+" "+" ".join(row.split()[1:])+"\n"
                continue

            if spring_flg:
                spring_dict["O"] = row
                continue
        self.species_dict = species_dict
        self.pot_dict = pot_dict
        self.spring_dict = spring_dict


class BuckinghamPotBush(object):
    """
    Generate the Buckingham Potential Table from the bush.lib file included
    with GULP
    Ref: 1) T.S.Bush, J.D.Gale, C.R.A.Catlow and P.D. Battle,  J. Mater Chem.,
            4, 831-837 (1994)
    """

    def __init__(self, verbose=False):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        fid = open(os.path.join(module_dir, 'bush.lib'), 'rU')
        # In lewis.lib there is no shell for cation
        species_dict, pot_dict, spring_dict  = {}, {}, {}
        sp_flg, pot_flg, spring_flg = False, False, False
        for row in fid:
            if row[0] == "#":
                continue
            if row.split()[0] == "species":
                sp_flg, pot_flg, spring_flg = True, False, False
                continue
            if row.split()[0] == "buckingham":
                sp_flg, pot_flg, spring_flg = False, True, False
                continue
            if row.split()[0] == "spring":
                sp_flg, pot_flg, spring_flg = False, False, True
                continue

            met = row.split()[0]
            if sp_flg:
                if met not in species_dict.keys():
                    species_dict[met] = {'inp_str':'','oxi':0}
                species_dict[met]['inp_str'] += row
                species_dict[met]['oxi'] += float(row.split()[2])

            if pot_flg:
                pot_dict[met] = row

            if spring_flg:
                spring_dict[met] = row

        #Fill the null keys in spring dict with empty strings
        for key in pot_dict.keys():
            if key not in spring_dict.keys():
                spring_dict[key] = ""

        self.species_dict = species_dict
        self.pot_dict = pot_dict
        self.spring_dict = spring_dict


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
        
