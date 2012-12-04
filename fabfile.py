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
import glob

from fabric.api import local, lcd
from fabric.state import env

def makedoc():
    with lcd("docs"):
        local("sphinx-apidoc -o . -f ../pymatgen")
        local("rm pymatgen.*.tests.rst")
        for f in glob.glob("docs/*.rst"):
            if f.startswith('docs/pymatgen') and f.endswith('rst'):
                newoutput = []
                suboutput = []
                subpackage = False
                with open(f, 'r') as fid:
                    for line in fid:
                        clean = line.strip()
                        if clean == "Subpackages":
                            subpackage = True
                        if not subpackage and not clean.endswith("tests"):
                            newoutput.append(line)
                        else:
                            if not clean.endswith("tests"):
                                suboutput.append(line)
                            if clean.startswith("pymatgen") and not clean.endswith("tests"):
                                newoutput.extend(suboutput)
                                subpackage = False
                                suboutput = []

                with open(f, 'w') as fid:
                    fid.write("".join(newoutput))

        local("make html")
        local("cp nature_mp.css ../../docs/pymatgen/html/static/nature.css")
        local("cp favicon.ico ../../docs/pymatgen/html/static/favicon.ico")

def publish():
    local("python setup.py release")

def test():
    local("nosetests")

def setver():
    from pymatgen import __version__
    local("sed s/version=.*,/version=\\\"{}\\\",/ setup.py > newsetup".format(__version__))
    local("mv newsetup setup.py")

def update_dev_doc():
    makedoc()
    with lcd("../docs/pymatgen/html/"):
        local("git add .")
        local("git commit -a -m \"Update dev docs\"")
        local("git push origin gh-pages")

def release():
    setver()
    test()
    makedoc()
    publish()
