#!/usr/bin/env python
import argparse
import os
import re
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

parser = argparse.ArgumentParser(description='''Convenient vasp run analyzer which can recursively go into a directory to search results.
Author: Shyue Ping Ong
Version: 1.0
Last updated: Oct 26 2011''')
parser.add_argument('directories', metavar='dir', default='.', type=str, nargs = '*', help='directory to process (default to .)')
parser.add_argument('-e', '--energies', dest='get_energies', action='store_const', const=True, help='print energies')
parser.add_argument('-m', '--mag', dest="ion_list", type=str, nargs = 1, help='print magmoms. ION LIST can be a range (e.g., 1-2) or the string "All" for all ions.')

args = parser.parse_args()
if args.get_energies:
    map(get_energies, args.directories)
if args.ion_list:
    ion_list = list()
    (start, end) = map(int, re.split("-", args.ion_list[0]))
    ion_list = range(start, end + 1)
    get_magnetizations(args.directories, ion_list)
