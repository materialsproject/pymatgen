#!/usr/bin/env python

'''
Created on Apr 29, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 29, 2012"

import os
import subprocess
import glob
import shlex
import shutil

sphinxcmd = "sphinx-apidoc -o . -f ../pymatgen"
args = shlex.split(sphinxcmd)
p = subprocess.Popen(args)
output = p.communicate()[0]

for f in glob.glob("*.rst"):
    if f.endswith('tests.rst'):
        os.remove(f)
    elif f.startswith('pymatgen') and f.endswith('rst'):
        newoutput = []
        suboutput = []
        subpackage = False
        with open(f, 'r') as fid:
            for line in fid:
                if line.strip() == "Subpackages":
                    subpackage = True
                if not subpackage:
                    newoutput.append(line)
                else:
                    suboutput.append(line)
                    if line.strip().startswith("pymatgen") and not line.strip().endswith("tests"):
                        newoutput.extend(suboutput)
                        subpackage = False
                        suboutput = []

        with open(f, 'w') as fid:
            fid.write("".join(newoutput))

p = subprocess.Popen(["make", "html"])
output = p.communicate()[0]

shutil.copyfile("nature_mp.css", os.path.join("..", "..", "docs", "pymatgen", "html", "static", "nature.css"))
shutil.copyfile("favicon.ico", os.path.join("..", "..", "docs", "pymatgen", "html", "static", "favicon.ico"))

#cp favicon.ico .. / .. / docs / pymatgen / html / static / favicon.ico
#cp pymatgen.png .. / .. / docs / pymatgen / html / static / pymatgen.png

