"""
Pyinvoke tasks.py file for automating releases and admin stuff.

Author: Shyue Ping Ong
"""
import datetime
import glob
import json
import os
import re
import subprocess
import webbrowser

import requests
from invoke import task
from monty.os import cd

from pymatgen.core import __version__ as CURRENT_VER


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
        ctx.run("sphinx-apidoc --implicit-namespaces --separate -d 7 -o . -f ../pymatgen")
        ctx.run("rm *.tests.*rst")
        for f in glob.glob("*.rst"):
            if f.startswith("pymatgen") and f.endswith("rst"):
                newoutput = []
                suboutput = []
                subpackage = False
                with open(f) as fid:
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

                with open(f, "w") as fid:
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
        ctx.run('echo "pymatgen.org" > CNAME')
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
    ctx.run("rm docs/_static/pymatgen.docset.tgz", warn=True)
    ctx.run("doc2dash docs -n pymatgen -i docs/_images/pymatgen.png -u https://pymatgen.org/")
    plist = "pymatgen.docset/Contents/Info.plist"
    xml = []
    with open(plist) as f:
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
    #             xml.append("<version>%s</version>" % version)
    #         else:
    #             xml.append(l)
    # with open("docs/pymatgen.xml", "wt") as f:
    #     f.write("\n".join(xml))
    ctx.run("rm -r pymatgen.docset")
    ctx.run("cp docs_rst/conf-normal.py docs_rst/conf.py")


@task
def contribute_dash(ctx, version):
    make_dash(ctx)
    ctx.run("cp pymatgen.tgz ../Dash-User-Contributions/docsets/pymatgen/pymatgen.tgz")
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        with open("docset.json") as f:
            data = json.load(f)
            data["version"] = version
        with open("docset.json", "wt") as f:
            json.dump(data, f, indent=4)
        ctx.run('git commit --no-verify -a -m "Update to v%s"' % version)
        ctx.run("git push")
    ctx.run("rm pymatgen.tgz")


@task
def submit_dash_pr(ctx, version):
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        payload = {
            "title": "Update pymatgen docset to v%s" % version,
            "body": "Update pymatgen docset to v%s" % version,
            "head": "Dash-User-Contributions:master",
            "base": "master",
        }
        response = requests.post(
            "https://api.github.com/repos/materialsvirtuallab/Dash-User-Contributions/pulls", data=json.dumps(payload)
        )
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
    ctx.run('git commit -a -m "Update docs"')
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
def set_ver(ctx, version):
    with open("pymatgen/core/__init__.py") as f:
        contents = f.read()
        contents = re.sub(r"__version__ = .*\n", '__version__ = "%s"\n' % version, contents)

    with open("pymatgen/core/__init__.py", "wt") as f:
        f.write(contents)

    with open("setup.py") as f:
        contents = f.read()
        contents = re.sub(r"version=([^,]+),", 'version="%s",' % version, contents)

    with open("setup.py", "wt") as f:
        f.write(contents)


@task
def release_github(ctx, version):
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
        "tag_name": "v" + version,
        "target_commitish": "master",
        "name": "v" + version,
        "body": desc,
        "draft": False,
        "prerelease": False,
    }
    response = requests.post(
        "https://api.github.com/repos/materialsproject/pymatgen/releases",
        data=json.dumps(payload),
        headers={"Authorization": "token " + os.environ["GITHUB_RELEASES_TOKEN"]},
    )
    print(response.text)


@task
def post_discourse(ctx, version):
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
    raw = "v" + version + "\n\n" + desc
    payload = {
        "topic_id": 36,
        "raw": raw,
    }
    response = requests.post(
        "https://discuss.matsci.org/c/pymatgen/posts.json",
        data=payload,
        params={"api_username": os.environ["DISCOURSE_API_USERNAME"], "api_key": os.environ["DISCOURSE_API_KEY"]},
    )
    print(response.text)


@task
def update_changelog(ctx, version=datetime.datetime.now().strftime("%Y.%-m.%-d"), sim=False):
    """
    Create a preliminary change log using the git logs.

    :param ctx:
    """
    output = subprocess.check_output(["git", "log", "--pretty=format:%s", "v%s..HEAD" % CURRENT_VER])
    lines = []
    misc = []
    for l in output.decode("utf-8").strip().split("\n"):
        m = re.match(r"Merge pull request \#(\d+) from (.*)", l)
        if m:
            pr_number = m.group(1)
            contrib, pr_name = m.group(2).split("/", 1)
            response = requests.get(f"https://api.github.com/repos/materialsproject/pymatgen/pulls/{pr_number}")
            lines.append(f"* PR #{pr_number} from @{contrib} {pr_name}")
            if "body" in response.json():
                for ll in response.json()["body"].split("\n"):
                    ll = ll.strip()
                    if ll in ["", "## Summary"]:
                        continue
                    elif ll.startswith("## Checklist") or ll.startswith("## TODO"):
                        break
                    lines.append(f"    {ll}")
        misc.append(l)
    with open("CHANGES.rst") as f:
        contents = f.read()
    l = "=========="
    toks = contents.split(l)
    head = "\n\nv%s\n" % version + "-" * (len(version) + 1) + "\n"
    toks.insert(-1, head + "\n".join(lines))
    if not sim:
        with open("CHANGES.rst", "w") as f:
            f.write(toks[0] + l + "".join(toks[1:]))
        ctx.run("open CHANGES.rst")
    else:
        print(toks[0] + l + "".join(toks[1:]))
    print("The following commit messages were not included...")
    print("\n".join(misc))


@task
def release(ctx, version=datetime.datetime.now().strftime("%Y.%-m.%-d"), nodoc=False):
    """
    Run full sequence for releasing pymatgen.

    :param ctx:
    :param nodoc: Whether to skip doc generation.
    """
    ctx.run("rm -r dist build pymatgen.egg-info", warn=True)
    set_ver(ctx, version)
    if not nodoc:
        make_doc(ctx)
        ctx.run("git add .")
        ctx.run('git commit -a -m "Update docs"')
        ctx.run("git push")
    release_github(ctx, version)
    ctx.run("rm -f dist/*.*", warn=True)
    ctx.run("python setup.py sdist bdist_wheel", warn=True)
    ctx.run("twine upload --skip-existing dist/*.whl", warn=True)
    ctx.run("twine upload --skip-existing dist/*.tar.gz", warn=True)
    # post_discourse(ctx, warn=True)


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
