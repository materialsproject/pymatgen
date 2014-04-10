"""
Deployment file to facilitate releases of pymatgen.
"""

__author__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Mar 6, 2014"

import glob
import os
import webbrowser

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


def merge_stable():
    local("git checkout stable")
    local("git pull")
    local("git merge master")
    local("git push")
    local("git checkout master")


def log_ver():
    filepath = os.path.join(os.environ["HOME"], "Dropbox", "Public",
                            "pymatgen", ver)
    with open(filepath, "w") as f:
        f.write("Release")


def release(skip_test=False):
    setver()
    if not skip_test:
        local("nosetests")
    publish()
    log_ver()
    update_doc()
    merge_stable()


def opendoc():
    pth = os.path.abspath("docs/_build/html/index.html")
    webbrowser.open("file://" + pth)
