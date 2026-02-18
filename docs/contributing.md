---
layout: default
title: Contributing
nav_order: 6
---

# Recommended developer workflow

For developers interested in expanding `pymatgen` for their own purposes, we recommend forking `pymatgen` from the [`pymatgen` GitHub repo](https://github.com/materialsproject/pymatgen). Here's a recommended workflow (updated August 2025):

1. Create a free GitHub account (if you don't already have one) and perform the necessary setup (e.g., [setup SSH keys](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) etc.).

2. [Fork the `pymatgen` GitHub repo](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo), i.e., go to the main [`pymatgen` GitHub repo](https://github.com/materialsproject/pymatgen) and click "Fork" to create a copy of the `pymatgen` code base in your own GitHub account.

3. Install [`git`](https://git-scm.com/downloads) on your local machine (if you haven't already).

4. Clone *your forked repo* to your local machine. You will work mostly with your forked repo and only publish changes when they are ready to be reviewed and merged:

    ```sh
    git clone https://github.com/<username>/pymatgen
    cd pymatgen  # change into pymatgen directory
    ```

   Note that the entire GitHub repo is fairly large because of the presence of test files, but these are necessary for rigorous testing. Alternatively, you could make a [shallow clone](https://git-scm.com/docs/git-clone#Documentation/git-clone.txt---depthdepth):

    ```sh
    git clone https://github.com/<username>/pymatgen  --depth 1

    # Convert into a full clone
    git fetch --unshallow
    ```

5. Install the [uv package manager](https://docs.astral.sh/uv/getting-started/installation/):

6. Create a [virtual env](https://docs.astral.sh/uv/pip/environments/) for pymatgen:

    ```sh
    uv venv  # A virtual env will be created in the `.venv` directory in the repo.
    ```

7. Set up development environment via `uv`:

    ```sh
    uv sync
    # uv sync --extra optional  # Install a specific optional dependency
    # uv sync --all-extras      # Install all optional dependencies

    uv run pre-commit install  # Install pre-commit hook for linters.
    ```

8. Make a new branch for your contributions:

    ```sh
    git checkout -b my-new-fix-or-feature  # should be run from up-to-date `master` branch
    ```

9. Code (see [Coding Guidelines](#coding-guidelines)). Commit early and commit often. Keep your code up to date. You need to add the main repository to the list of your remotes.

    ```sh
    git remote add upstream https://github.com/materialsproject/pymatgen
    ```

   Make sure your repository is clean (no uncommitted changes) and is currently on the master branch. If not, commit or stash any changes and switch to the master:

    ```sh
    git checkout master
    ```

   Then you can pull all the new commits from the main line:

    ```sh
    git pull upstream master
    ```

   Remember, pull is a combination of the commands fetch and merge, so there may be merge conflicts to be manually resolved.

10. Publish your contributions. Assuming that you now have a couple of commits that you would like to contribute to the main repository. Please follow the following steps:

    1. If your change is based on a relatively old state of the main repository, then you should probably bring your repository up-to-date first to see if the change is not creating any merge conflicts.

    2. Check that everything compiles cleanly and passes all tests. The `pymatgen` repo comes with a complete set of tests for all modules. If you have written new modules or methods, you must write tests for the new code as well (see [Coding Guidelines](#coding-guidelines)).

    3. If everything is ok, publish the commits to your GitHub repository.

    ```sh
    git push origin master
    ```

11. Now that your commit is published, it doesn't mean that it has already been merged into the main repository. You should [create a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to `pymatgen` maintainers.

"Work-in-progress" pull requests are encouraged, especially if this is your first time contribution, and the maintainers will be happy to help or provide code review as necessary. Put "\[WIP\]" in the title of your pull request to indicate it's not ready to be merged.

## Coding Guidelines

Given that `pymatgen` is intended to be a long-term code base, we adopt very strict quality control and coding guidelines for all contributions to `pymatgen`. The following must be satisfied for your contributions to be accepted into `pymatgen`.

1. **Unit tests** are required. The only way to minimize code regression is to ensure that all code is well tested. Untested contributions will not be accepted. Please use the modern [`pytest`](https://docs.pytest.org/en/stable/how-to/index.html) framework for tests and avoid the legacy `unittest` style.

    Run `pytest` in your local repo and fix all errors before continuing further:
    ```sh
    # (Recommended) Run tests for a specific module
    uv run pytest tests/path

    # Run the full test suite
    uv run pytest tests
    ```

2. [**Python PEP 8** code style](https://python.org/dev/peps/pep-0008). We allow a few exceptions when they are well-justified (e.g., Element's atomic number is given a variable name of capital Z, in line with accepted scientific convention), but generally, PEP 8 must be observed. Code style will be automatically checked for all PRs and must pass before any PR is merged. To aid you, you can install and run the same set of formatters and linters that will run in CI using:

    ```sh
    uv run pre-commit install  # ensure linters are run prior to all future commits

    # (Optional) Run pre-commit manually
    uv run pre-commit run --files path/to/changed/files  # ensure certain file don't offend linters
    # or
    uv run pre-commit run --all-files  # ensure your entire codebase passes linters
    ```

3. **Python 3**. Check the `pyproject.toml` for the supported Python versions.

4. **Documentation** is required for all modules, classes and methods. In particular, the method doc strings should make clear the arguments expected and the return values. For complex algorithms (e.g., an Ewald summation), a summary of the algorithm should be provided and preferably with a link to a publication outlining the method in detail.

For the above, if in doubt, please refer to the `core` classes in `pymatgen` for examples of what is expected.
