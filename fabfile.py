#!/usr/bin/env python

"""
Deployment file to facilitate releases of pymatgen.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 29, 2012"

import glob
import os

from fabric.api import local, lcd
from pymatgen import __version__ as ver


def makedoc():
    with lcd("examples"):
        local("ipython nbconvert --to html *.ipynb")
        local("mv *.html ../docs/_static")
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
        local("cp _static/* _build/html/_static")

        #This makes sure pymatgen.org works to redirect to the Gihub page
        local("echo \"pymatgen.org\" > _build/html/CNAME")
        #Avoid ths use of jekyll so that _dir works as intended.
        local("touch _build/html/.nojekyll")


def publish():
    local("python setup.py release")


def test():
    local("nosetests")


def setver():
    local("sed s/version=.*,/version=\\\"{}\\\",/ setup.py > newsetup"
          .format(ver))
    local("mv newsetup setup.py")


def update_doc():
    makedoc()
    with lcd("docs/_build/html/"):
        local("git add .")
        local("git commit -a -m \"Update dev docs\"")
        local("git push origin gh-pages")


def log_ver():
    filepath = os.path.join(os.environ["HOME"], "Dropbox", "Public",
                            "pymatgen", ver)
    with open(filepath, "w") as f:
        f.write("Release")


def release():
    setver()
    test()
    publish()
    log_ver()
    update_doc()
