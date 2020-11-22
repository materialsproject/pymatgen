"""
Pyinvoke tasks.py file for automating releases and admin stuff.

Author: Shyue Ping Ong
"""

from invoke import task
import glob
import os
import json
import webbrowser
import requests
import re
import subprocess
import datetime

from monty.os import cd
from pymatgen import __version__ as CURRENT_VER

NEW_VER = datetime.datetime.today().strftime("%Y.%-m.%-d")


@task
def make_doc(ctx):
    """
    Generate API documentation + run Sphinx.

    :param ctx:
    """
    with open("CHANGES.rst") as f:
        contents = f.read()

    toks = re.split(r"\-{3,}", contents)
    n = len(toks[0].split()[-1])
    changes = [toks[0]]
    changes.append("\n" + "\n".join(toks[1].strip().split("\n")[0:-1]))
    changes = ("-" * n).join(changes)

    with open("docs_rst/latest_changes.rst", "w") as f:
        f.write(changes)

    with cd("docs_rst"):
        ctx.run("cp ../CHANGES.rst change_log.rst")
        ctx.run("rm pymatgen.*.rst", warn=True)
        ctx.run("sphinx-apidoc --separate -d 7 -o . -f ../pymatgen")
        ctx.run("rm pymatgen*.tests.*rst")
        for f in glob.glob("*.rst"):
            if f.startswith('pymatgen') and f.endswith('rst'):
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
        ctx.run("make html")

        ctx.run("cp _static/* ../docs/html/_static", warn=True)

    with cd("docs"):
        ctx.run("rm *.html", warn=True)
        ctx.run("cp -r html/* .", warn=True)
        ctx.run("rm -r html", warn=True)
        ctx.run("rm -r doctrees", warn=True)
        ctx.run("rm -r _sources", warn=True)
        ctx.run("rm -r _build", warn=True)

        # This makes sure pymatgen.org works to redirect to the Github page
        ctx.run("echo \"pymatgen.org\" > CNAME")
        # Avoid the use of jekyll so that _dir works as intended.
        ctx.run("touch .nojekyll")


@task
def make_dash(ctx):
    """
    Make customized doc version for Dash

    :param ctx:
    """
    ctx.run("cp docs_rst/conf-docset.py docs_rst/conf.py")
    make_doc(ctx)
    ctx.run('rm docs/_static/pymatgen.docset.tgz', warn=True)
    ctx.run('doc2dash docs -n pymatgen -i docs/_images/pymatgen.png -u https://pymatgen.org/')
    plist = "pymatgen.docset/Contents/Info.plist"
    xml = []
    with open(plist, "rt") as f:
        for l in f:
            xml.append(l.strip())
            if l.strip() == "<dict>":
                xml.append("<key>dashIndexFilePath</key>")
                xml.append("<string>index.html</string>")
    with open(plist, "wt") as f:
        f.write("\n".join(xml))
    ctx.run('tar --exclude=".DS_Store" -cvzf pymatgen.tgz pymatgen.docset')
    # xml = []
    # with open("docs/pymatgen.xml") as f:
    #     for l in f:
    #         l = l.strip()
    #         if l.startswith("<version>"):
    #             xml.append("<version>%s</version>" % NEW_VER)
    #         else:
    #             xml.append(l)
    # with open("docs/pymatgen.xml", "wt") as f:
    #     f.write("\n".join(xml))
    ctx.run('rm -r pymatgen.docset')
    ctx.run("cp docs_rst/conf-normal.py docs_rst/conf.py")


@task
def contribute_dash(ctx):
    make_dash(ctx)
    ctx.run('cp pymatgen.tgz ../Dash-User-Contributions/docsets/pymatgen/pymatgen.tgz')
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        with open("docset.json", "rt") as f:
            data = json.load(f)
            data["version"] = NEW_VER
        with open("docset.json", "wt") as f:
            json.dump(data, f, indent=4)
        ctx.run('git commit --no-verify -a -m "Update to v%s"' % NEW_VER)
        ctx.run('git push')
    ctx.run("rm pymatgen.tgz")


@task
def submit_dash_pr(ctx):
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        payload = {
            "title": "Update pymatgen docset to v%s" % NEW_VER,
            "body": "Update pymatgen docset to v%s" % NEW_VER,
            "head": "Dash-User-Contributions:master",
            "base": "master"
        }
        response = requests.post(
            "https://api.github.com/repos/materialsvirtuallab/Dash-User-Contributions/pulls",
            data=json.dumps(payload))
        print(response.text)


@task
def update_doc(ctx):
    """
    Update the web documentation.

    :param ctx:
    """
    ctx.run("cp docs_rst/conf-normal.py docs_rst/conf.py")
    make_doc(ctx)
    ctx.run("git add .")
    ctx.run("git commit -a -m \"Update docs\"")
    ctx.run("git push")


@task
def publish(ctx):
    """
    Upload release to Pypi using twine.

    :param ctx:
    """
    ctx.run("rm dist/*.*", warn=True)
    ctx.run("python setup.py sdist bdist_wheel")
    ctx.run("twine upload dist/*")


@task
def set_ver(ctx):
    with open("pymatgen/__init__.py", "rt") as f:
        contents = f.read()
        contents = re.sub(r"__version__ = .*\n", '__version__ = "%s"\n' % NEW_VER, contents)

    with open("pymatgen/__init__.py", "wt") as f:
        f.write(contents)

    with open("setup.py", "rt") as f:
        contents = f.read()
        contents = re.sub(r'version=([^,]+),', 'version="%s",' % NEW_VER, contents)

    with open("setup.py", "wt") as f:
        f.write(contents)


@task
def merge_stable(ctx):
    """
    Tag and merge into stable branch.

    :param ctx:
    """
    ctx.run("git checkout stable")
    ctx.run("git pull")
    ctx.run("git merge master")
    ctx.run("git push")
    ctx.run("git checkout master")


@task
def release_github(ctx):
    """
    Release to Github using Github API.

    :param ctx:
    """
    with open("CHANGES.rst") as f:
        contents = f.read()
    toks = re.split(r"\-+", contents)
    desc = toks[1].strip()
    toks = desc.split("\n")
    desc = "\n".join(toks[:-1]).strip()
    payload = {
        "tag_name": "v" + NEW_VER,
        "target_commitish": "master",
        "name": "v" + NEW_VER,
        "body": desc,
        "draft": False,
        "prerelease": False
    }
    response = requests.post(
        "https://api.github.com/repos/materialsproject/pymatgen/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_RELEASES_TOKEN"]})
    print(response.text)


@task
def post_discourse(ctx):
    """
    Post release announcement to http://discuss.matsci.org/c/pymatgen.

    :param ctx:
    """
    with open("CHANGES.rst") as f:
        contents = f.read()
    toks = re.split(r"\-+", contents)
    desc = toks[1].strip()
    toks = desc.split("\n")
    desc = "\n".join(toks[:-1]).strip()
    raw = "v" + NEW_VER + "\n\n" + desc
    payload = {
        "topic_id": 36,
        "raw": raw,
    }
    response = requests.post(
        "https://discuss.matsci.org/c/pymatgen/posts.json",
        data=payload,
        params={
            "api_username": os.environ["DISCOURSE_API_USERNAME"],
            "api_key": os.environ["DISCOURSE_API_KEY"]}
    )
    print(response.text)


@task
def update_changelog(ctx):
    """
    Create a preliminary change log using the git logs.

    :param ctx:
    """
    output = subprocess.check_output(["git", "log", "--pretty=format:%s",
                                      "v%s..HEAD" % CURRENT_VER])
    lines = ["* " + l for l in output.decode("utf-8").strip().split("\n")]
    with open("CHANGES.rst") as f:
        contents = f.read()
    l = "=========="
    toks = contents.split(l)
    head = "\n\nv%s\n" % NEW_VER + "-" * (len(NEW_VER) + 1) + "\n"
    toks.insert(-1, head + "\n".join(lines))
    with open("CHANGES.rst", "w") as f:
        f.write(toks[0] + l + "".join(toks[1:]))
    ctx.run("open CHANGES.rst")


@task
def release(ctx, nodoc=False):
    """
    Run full sequence for releasing pymatgen.

    :param ctx:
    :param notest: Whether to skip tests.
    :param nodoc: Whether to skip doc generation.
    """
    ctx.run("rm -r dist build pymatgen.egg-info", warn=True)
    set_ver(ctx)
    if not nodoc:
        make_doc(ctx)
        ctx.run("git add .")
        ctx.run("git commit -a -m \"Update docs\"")
        ctx.run("git push")
    merge_stable(ctx)
    release_github(ctx)
    post_discourse(ctx, warn=True)


@task
def open_doc(ctx):
    """
    Open local documentation in web browser.

    :param ctx:
    """
    pth = os.path.abspath("docs/_build/html/index.html")
    webbrowser.open("file://" + pth)


@task
def lint(ctx):
    for cmd in ["pycodestyle", "mypy", "flake8", "pydocstyle"]:
        ctx.run("%s pymatgen" % cmd)
