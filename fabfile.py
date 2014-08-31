"""
Deployment file to facilitate releases of pymatgen.
"""

__author__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Mar 6, 2014"

import glob
import os
import json
import webbrowser
import requests
import re

from fabric.api import local, lcd
from pymatgen import __version__ as ver


def make_doc():
    with open("CHANGES") as f:
        contents = f.read()

    toks = re.split("-+", contents)
    n = len(toks[0].split()[-1])

    changes = ("-" * n).join(toks[0:2])

    with open("LATEST_CHANGES", "w") as f:
        f.write(changes)

    with lcd("examples"):
        local("ipython nbconvert --to html *.ipynb")
        local("mv *.html ../docs/_static")
    with lcd("docs"):
        local("sphinx-apidoc -d 6 -o . -f ../pymatgen")
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

    os.remove("LATEST_CHANGES")

def publish():
    local("python setup.py release")


def setver():
    local("sed s/version=.*,/version=\\\"{}\\\",/ setup.py > newsetup"
          .format(ver))
    local("mv newsetup setup.py")


def update_doc():
    make_doc()
    with lcd("docs/_build/html/"):
        local("git add .")
        local("git commit -a -m \"Update dev docs\"")
        local("git push origin gh-pages")


def merge_stable():
    local("git commit -a -m \"v%s release\"" % ver)
    local("git push")
    local("git checkout stable")
    local("git pull")
    local("git merge master")
    local("git push")
    local("git checkout master")


def release_github():
    desc = []
    read = False
    with open("docs/index.rst") as f:
        for l in f:
            if l.strip() == "v" + ver:
                read = True
            elif l.strip() == "":
                read = False
            elif read:
                desc.append(l.rstrip())
    desc.pop(0)
    payload = {
        "tag_name": "v" + ver,
        "target_commitish": "master",
        "name": "v" + ver,
        "body": "\n".join(desc),
        "draft": False,
        "prerelease": False
    }

    response = requests.post(
        "https://api.github.com/repos/materialsproject/pymatgen/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_RELEASES_TOKEN"]})
    print response.text


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
    release_github()


def open_doc():
    pth = os.path.abspath("docs/_build/html/index.html")
    webbrowser.open("file://" + pth)
