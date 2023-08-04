---
layout: default
title: Contributing
nav_order: 6
---

# Collaborative Github Workflow

For developers interested in expanding `pymatgen` for their own purposes, we recommend forking `pymatgen` directly from the [`pymatgen` GitHub repo](https://github.com/materialsproject/pymatgen). Here's a typical workflow:

1. Create a free GitHub account (if you don't already have one) and perform the necessary setup (e.g., install SSH keys etc.).

1. Fork the `pymatgen` GitHub repo, i.e., go to the main [`pymatgen` GitHub repo](https://github.com/materialsproject/pymatgen) and click fork to create a copy of the `pymatgen` code base on your own GitHub account.

1. Install git on your local machine (if you don't already have it).

1. Clone *your forked repo* to your local machine. You will work mostly with your local repo and only publish changes when they are ready to be merged:

    ```sh
    git clone git@github.com:YOURNAME/pymatgen
    ```

    Note that the entire Github repo is fairly large because of the presence of test files, but these are absolutely necessary for rigorous testing of the code.

1. It is highly recommended you install all the optional dependencies as well.

1. Code (see [Coding Guidelines](#coding-guidelines)). Commit early and commit often. Keep your code up to date. You need to add the main repository to the list of your remotes.

    ```sh
    git remote add upstream git://github.com/materialsproject/`pymatgen`
    ```

    Make sure your repository is clean (no uncommitted changes) and is currently on the master branch. If not, commit or stash any changes and switch to the master.

    ```sh
    git checkout master
    ```

    Then you can pull all the new commits from the main line

    ```sh
    git pull upstream master
    ```

    Remember, pull is a combination of the commands fetch and merge, so there may be merge conflicts to be manually resolved.

1. Publish your contributions. Assuming that you now have a couple of commits that you would like to contribute to the main repository. Please follow the following steps:

    1. If your change is based on a relatively old state of the main repository, then you should probably bring your repository up-to-date first to see if the change is not creating any merge conflicts.

    1. Check that everything compiles cleanly and passes all tests. The `pymatgen` repo comes with a complete set of tests for all modules. If you have written new modules or methods, you must write tests for the new code as well (see [Coding Guidelines](#coding-guidelines)). Install and run `pytest` in your local repo directory and fix all errors before continuing further.

    1. If everything is ok, publish the commits to your GitHub repository.

    ```sh
    git push origin master
    ```

1. Now that your commit is published, it doesn't mean that it has already been merged into the main repository. You should issue a merge request to `pymatgen` maintainers. They will pull your commits and run their own tests before releasing.

"Work-in-progress" pull requests are encouraged, especially if this is your first time contributing to `pymatgen`, and the maintainers will be happy to help or provide code review as necessary. Put "\[WIP\]" in the title of your pull request to indicate it's not ready to be merged.

## Coding Guidelines

Given that `pymatgen` is intended to be a long-term code base, we adopt very strict quality control and coding guidelines for all contributions to `pymatgen`. The following must be satisfied for your contributions to be accepted into `pymatgen`.

1. **Unit tests** are required for all new modules and methods. The only way to minimize code regression is to ensure that all code is well tested. Untested contributions will not be accepted.
1. **Python PEP 8** [code style](https://python.org/dev/peps/pep-0008). We allow a few exceptions when they are well-justified (e.g., Element's atomic number is given a variable name of capital Z, in line with accepted scientific convention), but generally, PEP 8 must be observed. Code style will be automatically checked for all PRs and must pass before any PR is merged. To aid you, you can install and run the same set of formatters and linters that will run in CI using

   ```sh
   pre-commit install  # ensures linters are run prior to all future commits
   pre-commit run --files path/to/changed/files  # ensure your current uncommitted changes don't offend linters
   # or
    pre-commit run --all-files  # ensure your entire codebase passes linters
   ```

1. **Python 3**. We only support Python 3.8+.
1. **Documentation** is required for all modules, classes and methods. In particular, the method doc strings should make clear the arguments expected and the return values. For complex algorithms (e.g., an Ewald summation), a summary of the algorithm should be provided, and preferably with a link to a publication outlining the method in detail.

For the above, if in doubt, please refer to the core classes in `pymatgen` for examples of what is expected.
