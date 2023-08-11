"""
Pyinvoke tasks.py file for automating releases and admin stuff.

To cut a new pymatgen release, use `invoke update-changelog` followed by `invoke release`.

Author: Shyue Ping Ong
"""

from __future__ import annotations

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

from pymatgen.core import __version__


@task
def make_doc(ctx):
    """
    Generate API documentation + run Sphinx.

    :param ctx:
    """
    with cd("docs"):
        ctx.run("touch index.rst")
        ctx.run("rm pymatgen.*.rst", warn=True)
        ctx.run("rm pymatgen.*.md", warn=True)
        ctx.run("sphinx-apidoc --implicit-namespaces -P -M -d 7 -o . -f ../pymatgen ../**/tests/*")
        # ctx.run("rm pymatgen*.html", warn=True)
        # ctx.run("sphinx-build -b html . ../docs")  # HTML building.
        ctx.run("sphinx-build -M markdown . .")
        ctx.run("rm *.rst", warn=True)
        ctx.run("cp markdown/pymatgen*.md .")
        ctx.run("rm pymatgen*tests*.md", warn=True)
        for fn in glob.glob("pymatgen*.md"):
            with open(fn) as f:
                lines = [line.rstrip() for line in f if "Submodules" not in line]
            if fn == "pymatgen.md":
                preamble = ["---", "layout: default", "title: API Documentation", "nav_order: 6", "---", ""]
            else:
                preamble = ["---", "layout: default", "title: " + fn, "nav_exclude: true", "---", "",
                            "1. TOC", "{:toc}", ""]
            with open(fn, "w") as f:
                f.write("\n".join(preamble + lines))
        ctx.run("rm -r markdown", warn=True)
        # ctx.run("cp ../README.md index.md")
        ctx.run("cp ../CHANGES.md CHANGES.md")
        ctx.run("rm -rf doctrees", warn=True)


@task
def make_dash(ctx):
    """
    Make customized doc version for Dash.

    :param ctx:
    """
    ctx.run("cp docs_rst/conf-docset.py docs_rst/conf.py")
    make_doc(ctx)
    ctx.run("rm docs/_static/pymatgen.docset.tgz", warn=True)
    ctx.run("doc2dash docs -n pymatgen -i docs/_images/pymatgen.svg -u https://pymatgen.org/")
    plist = "pymatgen.docset/Contents/Info.plist"
    xml = []
    with open(plist) as file:
        for line in file:
            xml.append(line.strip())
            if line.strip() == "<dict>":
                xml.append("<key>dashIndexFilePath</key>")
                xml.append("<string>index.html</string>")
    with open(plist, "w") as file:
        file.write("\n".join(xml))
    ctx.run('tar --exclude=".DS_Store" -cvzf pymatgen.tgz pymatgen.docset')
    ctx.run("rm -r pymatgen.docset")
    ctx.run("cp docs_rst/conf-normal.py docs_rst/conf.py")


@task
def contribute_dash(ctx, version):
    make_dash(ctx)
    ctx.run("cp pymatgen.tgz ../Dash-User-Contributions/docsets/pymatgen/pymatgen.tgz")
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        with open("docset.json") as file:
            data = json.load(file)
            data["version"] = version
        with open("docset.json", "w") as file:
            json.dump(data, file, indent=4)
        ctx.run(f'git commit --no-verify -a -m "Update to v{version}"')
        ctx.run("git push")
    ctx.run("rm pymatgen.tgz")


@task
def submit_dash_pr(ctx, version):
    with cd("../Dash-User-Contributions/docsets/pymatgen"):
        payload = {
            "title": f"Update pymatgen docset to v{version}",
            "body": f"Update pymatgen docset to v{version}",
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
    with open("pymatgen/core/__init__.py") as file:
        contents = file.read()
        contents = re.sub(r"__version__ = .*\n", f"__version__ = {version!r}\n", contents)

    with open("pymatgen/core/__init__.py", "w") as file:
        file.write(contents)

    with open("setup.py") as file:
        contents = file.read()
        contents = re.sub(r"version=([^,]+),", f"version={version!r},", contents)

    with open("setup.py", "w") as file:
        file.write(contents)


@task
def release_github(ctx, version):
    """
    Release to Github using Github API.

    :param ctx:
    """
    with open("CHANGES.md") as file:
        contents = file.read()
    tokens = re.split(r"\-+", contents)
    desc = tokens[1].strip()
    tokens = desc.split("\n")
    desc = "\n".join(tokens[:-1]).strip()
    payload = {
        "tag_name": f"v{version}",
        "target_commitish": "master",
        "name": f"v{version}",
        "body": desc,
        "draft": False,
        "prerelease": False,
    }
    response = requests.post(
        "https://api.github.com/repos/materialsproject/pymatgen/releases",
        data=json.dumps(payload),
        headers={"Authorization": f"token {os.environ['GITHUB_RELEASES_TOKEN']}"},
    )
    print(response.text)


def post_discourse(version):
    """
    Post release announcement to http://discuss.matsci.org/c/pymatgen.

    :param ctx:
    """
    with open("CHANGES.rst") as file:
        contents = file.read()
    toks = re.split(r"\-+", contents)
    desc = toks[1].strip()
    toks = desc.split("\n")
    desc = "\n".join(toks[:-1]).strip()
    raw = f"v{version}\n\n{desc}"
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
def update_changelog(ctx, version=None, dry_run=False):
    """
    Create a preliminary change log using the git logs.

    :param ctx:
    """
    version = version or f"{datetime.datetime.now():%Y.%-m.%-d}"
    output = subprocess.check_output(["git", "log", "--pretty=format:%s", f"v{__version__}..HEAD"])
    lines = []
    ignored_commits = []
    for line in output.decode("utf-8").strip().split("\n"):
        re_match = re.match(r"Merge pull request \#(\d+) from (.*)", line)
        if re_match and "materialsproject/dependabot/pip" not in line:
            pr_number = re_match.group(1)
            contributor, pr_name = re_match.group(2).split("/", 1)
            response = requests.get(f"https://api.github.com/repos/materialsproject/pymatgen/pulls/{pr_number}")
            lines.append(f"* PR #{pr_number} from @{contributor} {pr_name}")
            json_resp = response.json()
            if body := json_resp["body"]:
                for ll in map(str.strip, body.split("\n")):
                    if ll in ["", "## Summary"]:
                        continue
                    if ll.startswith(("## Checklist", "## TODO")):
                        break
                    lines.append(f"    {ll}")
        ignored_commits.append(line)
    with open("CHANGES.md") as file:
        contents = file.read()
    delim = "##"
    tokens = contents.split(delim)
    tokens.insert(1, f"## v{version}\n\n" + "\n".join(lines) + "\n")
    if dry_run:
        print(tokens[0] + "##".join(tokens[1:]))
    else:
        with open("CHANGES.md", "w") as file:
            file.write(tokens[0] + "##".join(tokens[1:]))
        ctx.run("open CHANGES.md")
    print("The following commit messages were not included...")
    print("\n".join(ignored_commits))


@task
def release(ctx, version=None, nodoc=False):
    """
    Run full sequence for releasing pymatgen.

    :param ctx:
    :param nodoc: Whether to skip doc generation.
    """
    version = version or f"{datetime.datetime.now():%Y.%-m.%-d}"
    ctx.run("rm -r dist build pymatgen.egg-info", warn=True)
    set_ver(ctx, version)
    ctx.run("black setup.py")
    ctx.run("black pymatgen/core/__init__.py")
    if not nodoc:
        make_doc(ctx)
        ctx.run("git add .")
        ctx.run('git commit --no-verify -a -m "Update docs"')
        ctx.run("git push")
    release_github(ctx, version)

    ctx.run("rm -f dist/*.*", warn=True)
    ctx.run("python setup.py sdist bdist_wheel", warn=True)
    check_egg_sources_txt_for_completeness()
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
    webbrowser.open(f"file://{pth}")


@task
def lint(ctx):
    for cmd in ["ruff", "mypy", "black", "pylint"]:
        ctx.run(f"{cmd} pymatgen")


def check_egg_sources_txt_for_completeness():
    """Check that all source and data files in pymatgen/ are listed in pymatgen.egg-info/SOURCES.txt."""
    src_txt = "pymatgen.egg-info/SOURCES.txt"
    if not os.path.exists(src_txt):
        raise FileNotFoundError(f"{src_txt} not found. Run `pip install .` to create")

    with open(src_txt) as file:
        sources = file.read()

    for src_file in sources.splitlines():
        if not os.path.exists(src_file):
            raise ValueError(f"{src_file} does not exist!")

    for ext in ("py", "json", "json.gz", "yaml", "csv"):
        for filepath in glob(f"pymatgen/**/*.{ext}", recursive=True):
            if "/tests/" in filepath or "dao" in filepath:
                continue
            if filepath not in sources:
                raise ValueError(f"{filepath} not found in {src_txt}")
