"""
Pyinvoke tasks.py file for automating releases and admin stuff.

To cut a new pymatgen release, use `invoke update-changelog` followed by `invoke release`.

Author: Shyue Ping Ong
"""

from __future__ import annotations

import datetime
import json
import os
import re
import subprocess
import webbrowser
from typing import TYPE_CHECKING

import requests
from invoke import task
from monty.os import cd

from pymatgen.core import __version__

if TYPE_CHECKING:
    from invoke import Context


@task
def make_doc(ctx: Context) -> None:
    """
    Generate API documentation + run Sphinx.

    Args:
        ctx (Context): The context.
    """
    with cd("docs"):
        ctx.run("touch apidoc/index.rst", warn=True)
        ctx.run("rm pymatgen.*.rst", warn=True)
        # ctx.run("rm pymatgen.*.md", warn=True)
        ctx.run("sphinx-apidoc --implicit-namespaces -M -d 7 -o apidoc -f ../pymatgen ../**/tests/*")

        # Note: we use HTML building for the API docs to preserve search functionality.
        ctx.run("sphinx-build -b html apidoc html")  # HTML building.
        ctx.run("rm apidocs/*.rst", warn=True)
        ctx.run("mv html/pymatgen*.html .")
        ctx.run("mv html/modules.html .")

        # ctx.run("cp markdown/pymatgen*.md .")
        # ctx.run("rm pymatgen*tests*.md", warn=True)
        # ctx.run("rm pymatgen*.html", warn=True)
        # for filename in glob("pymatgen*.md"):
        #     with open(filename) as file:
        #         lines = [line.rstrip() for line in file if "Submodules" not in line]
        #     if filename == "pymatgen.md":
        #         preamble = ["---", "layout: default", "title: API Documentation", "nav_order: 6", "---", ""]
        #     else:
        #         preamble = [
        #             "---",
        #             "layout: default",
        #             f"title: {filename}",
        #             "nav_exclude: true",
        #             "---",
        #             "",
        #             "1. TOC",
        #             "{:toc}",
        #             "",
        #         ]
        #     with open(filename, mode="w") as file:
        #         file.write("\n".join(preamble + lines))
        ctx.run("rm -r markdown", warn=True)
        ctx.run("rm -r html", warn=True)
        ctx.run('sed -I "" "s/_static/assets/g" pymatgen*.html')
        ctx.run("rm -rf doctrees", warn=True)


@task
def publish(ctx: Context) -> None:
    """
    Upload release to Pypi using twine.

    Args:
        ctx (Context): The context.
    """
    ctx.run("rm dist/*.*", warn=True)
    ctx.run("python setup.py sdist bdist_wheel")
    ctx.run("twine upload dist/*")


@task
def set_ver(ctx: Context, version: str):
    with open("setup.py") as file:
        contents = file.read()
        contents = re.sub(r"version=([^,]+),", f"version={version!r},", contents)

    with open("setup.py", mode="w") as file:
        file.write(contents)


@task
def release_github(ctx: Context, version: str) -> None:
    """
    Release to Github using Github API.

    Args:
        version (str): The version.
    """
    with open("docs/CHANGES.md") as file:
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


def post_discourse(version: str) -> None:
    """
    Post release announcement to http://discuss.matsci.org/c/pymatgen.

    Args:
        version (str): The version.
    """
    with open("CHANGES.rst") as file:
        contents = file.read()
    tokens = re.split(r"\-+", contents)
    desc = tokens[1].strip()
    tokens = desc.split("\n")
    desc = "\n".join(tokens[:-1]).strip()
    raw = f"v{version}\n\n{desc}"
    payload = {
        "topic_id": 36,
        "raw": raw,
    }
    response = requests.post(
        "https://discuss.matsci.org/c/pymatgen/posts.json",
        data=payload,
        params={
            "api_username": os.environ["DISCOURSE_API_USERNAME"],
            "api_key": os.environ["DISCOURSE_API_KEY"],
        },
    )
    print(response.text)


@task
def update_changelog(ctx: Context, version: str | None = None, dry_run: bool = False) -> None:
    """Create a preliminary change log using the git logs.

    Args:
        ctx (invoke.Context): The context object.
        version (str, optional): The version to use for the change log. If not provided, it will
            use the current date in the format 'YYYY.M.D'. Defaults to None.
        dry_run (bool, optional): If True, the function will only print the changes without
            updating the actual change log file. Defaults to False.
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
            response = requests.get(
                f"https://api.github.com/repos/materialsproject/pymatgen/pulls/{pr_number}", timeout=600
            )
            lines += [f"* PR #{pr_number} from @{contributor} {pr_name}"]
            json_resp = response.json()
            if body := json_resp["body"]:
                for ll in map(str.strip, body.split("\n")):
                    if ll in ("", "## Summary"):
                        continue
                    if ll.startswith(("## Checklist", "## TODO")):
                        break
                    lines += [f"    {ll}"]
        ignored_commits += [line]
    with open("docs/CHANGES.md") as file:
        contents = file.read()
    delim = "##"
    tokens = contents.split(delim)
    tokens.insert(1, f"## v{version}\n\n" + "\n".join(lines) + "\n")
    if dry_run:
        print(tokens[0] + "##".join(tokens[1:]))
    else:
        with open("docs/docs/CHANGES.md", mode="w") as file:
            file.write(tokens[0] + "##".join(tokens[1:]))
        ctx.run("open docs/CHANGES.md")
    print("The following commit messages were not included...")
    print("\n".join(ignored_commits))


@task
def release(ctx: Context, version: str | None = None, nodoc: bool = False) -> None:
    """
    Run full sequence for releasing pymatgen.

    Args:
        ctx (invoke.Context): The context object.
        version (str, optional): The version to release.
        nodoc (bool, optional): Whether to skip documentation generation.
    """
    version = version or f"{datetime.datetime.now():%Y.%-m.%-d}"
    ctx.run("rm -r dist build pymatgen.egg-info", warn=True)
    set_ver(ctx, version)
    if not nodoc:
        make_doc(ctx)
        ctx.run("git add .")
        ctx.run('git commit --no-verify -a -m "Update docs"')
        ctx.run("git push")
    release_github(ctx, version)

    ctx.run("rm -f dist/*.*", warn=True)
    ctx.run("python setup.py sdist bdist_wheel", warn=True)
    ctx.run("twine upload --skip-existing dist/*.whl", warn=True)
    ctx.run("twine upload --skip-existing dist/*.tar.gz", warn=True)
    # post_discourse(ctx, warn=True)


@task
def open_doc(ctx: Context) -> None:
    """
    Open local documentation in web browser.
    """
    pth = os.path.abspath("docs/_build/html/index.html")
    webbrowser.open(f"file://{pth}")


@task
def lint(ctx: Context) -> None:
    """
    Run linting tools.

    Args:
        ctx (invoke.Context): The context object.
    """
    for cmd in ("ruff", "mypy", "ruff format"):
        ctx.run(f"{cmd} pymatgen")
