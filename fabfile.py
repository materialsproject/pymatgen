# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Deployment file to facilitate releases of pymatgen.
Note that this file is meant to be run from the root directory of the pymatgen
repo.
"""

__author__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "Sep 1, 2014"

import glob
import os
import json
import webbrowser
import requests
import re
import subprocess

from fabric.api import local, lcd
from pymatgen import __version__ as ver


def make_doc():
    with open("CHANGES.rst") as f:
        contents = f.read()

    toks = re.split("\-{3,}", contents)
    n = len(toks[0].split()[-1])
    changes = [toks[0]]
    changes.append("\n" + "\n".join(toks[1].strip().split("\n")[0:-1]))
    changes = ("-" * n).join(changes)

    with open("docs/latest_changes.rst", "w") as f:
        f.write(changes)

    with lcd("examples"):
        local("jupyter nbconvert --to html *.ipynb")
        local("mv *.html ../docs/_static")
    with lcd("docs"):
        local("cp ../CHANGES.rst change_log.rst")
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


def publish():
    local("python setup.py release")


def setver():
    local("sed s/version=.*,/version=\\\"{}\\\",/ setup.py > newsetup"
          .format(ver))
    local("mv newsetup setup.py")


def update_doc():
    with lcd("docs/_build/html/"):
        local("git pull")
    make_doc()
    with lcd("docs/_build/html/"):
        local("git add .")
        local("git commit -a -m \"Update dev docs\"")
        local("git push origin gh-pages")


def update_coverage():
    with lcd("docs/_build/html/"):
        local("git pull")
    local("nosetests --config=nose.cfg --cover-html --cover-html-dir=docs/_build/html/coverage")
    update_doc()


def merge_stable():
    local("git commit -a -m \"v%s release\"" % ver)
    local("git push")
    local("git checkout stable")
    local("git pull")
    local("git merge master")
    local("git push")
    local("git checkout master")


def release_github():
    with open("CHANGES.rst") as f:
        contents = f.read()
    toks = re.split("\-+", contents)
    desc = toks[1].strip()
    toks = desc.split("\n")
    desc = "\n".join(toks[:-1]).strip()
    payload = {
        "tag_name": "v" + ver,
        "target_commitish": "master",
        "name": "v" + ver,
        "body": desc,
        "draft": False,
        "prerelease": False
    }
    response = requests.post(
        "https://api.github.com/repos/materialsproject/pymatgen/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_RELEASES_TOKEN"]})
    print response.text


def update_changelog():
    output = subprocess.check_output(["git", "log", "--pretty=format:%s",
                                      "v%s..HEAD" % ver])
    lines = ["* " + l for l in output.strip().split("\n")]
    with open("CHANGES.rst") as f:
        contents = f.read()
    l = "=========="
    toks = contents.split(l)
    toks.insert(-1, "\n\nvXXXX\n--------\n" + "\n".join(lines))
    with open("CHANGES.rst", "w") as f:
        f.write(toks[0] + l + "".join(toks[1:]))


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
