#!/usr/bin/python
import getopt
import os
import re
import sys
from pymatgen.io.vaspio import Vasprun, Outcar
from pymatgen.util.string_utils import str_aligned

# A convenience script engine using VaspObjects to do all manner of simple outputs.
# Shyue

def get_energies(mydir):
    all_data = []
    for (parent, subdirs, files) in os.walk(mydir):
        for f in files:
            if re.match("vasprun.xml.*", f):
                fullpath = os.path.join(parent, f)
                xmlrun = Vasprun(fullpath)
                energy = xmlrun.final_energy
                num_atoms = len(xmlrun.final_structure)            
                all_data.append((fullpath.lstrip("./"), energy, energy/num_atoms))
    print str_aligned(all_data, ("Directory", "Energy", "Energy/Atom"))
    
def get_magnetizations(mydir, ionList):
    print "%10s | %7s" % ("Ion", "Magmoms")
    print "-" * 20
    for (parent, subdirs, files) in os.walk(mydir):
        for f in files:
            if re.match("OUTCAR*", f):
                fullpath = os.path.join(parent, f)
                outcar = Outcar(fullpath)
                mags = outcar.magnetization
                mags = [m['tot'] for m in mags]
                allIons = xrange(len(mags))
                print "%16s" % (fullpath.lstrip("./"))
                if len(ionList) > 0:
                    allIons = ionList
                for ion in allIons:
                    print "%10d | %3.4f" % (ion, mags[ion])
                print "-" * 20

def usage():
    print "Convenient Vasp run analyzer which can recursively go into a dir to search results.\n\
Author: Shyue Version: 0.1 Last updated: Jul 29 2010\n\
Usage is as follows: VaspAnalyzer.py [arguments and options] dirnames.\n\
Default is to just print energies.\n\
Arguments supported:\n\
-energies : output all final energies in the vasprun.xml\n\
-magnetization=i-j : output magnetizations for ions i-j using OUTCAR"

def main():
        
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hem:", ["help", "energies", "magnetization="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    #print args
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        for d in args:
            if os.path.isdir(d):
                if o in ("-e", "--energies"):
                    get_energies(d)
                elif o in ("-m", "--magnetization"):
                    ionList = list()
                    (start, end) = map(int, re.split("-", a))
                    ionList = range(start, end + 1)
                    get_magnetizations(d, ionList)

if __name__ == "__main__":
    main()
