Contributing
============

Pymatgen is a community project. We welcome any contributions to either improve or expand its functionality. There are
two main ways of contributing to pymatgen::

1. **Direct contributions to pymatgen main distribution.** For improvements to existing pymatgen code or new
   functionality that is expected to be of broad interest and usage to the materials science community, you may make
   direct contributions to the pymatgen main distribution.
2. **Add-ons to pymatgen**. With effect from v2022.0.3, pymatgen now supports the development of add-on packages under
   the pymatgen, pymatgen.analysis, pymatgen.ext and and pymatgen.io
   `namespaces <http://packaging.python.org/guides/packaging-namespace-packages/>`_. If you are developing new
   functionality that is likely to be used by a sub-community of materials scientists, this is probably the best route
   going forward. Examples include functionality related to a specific technological application or io for a new
   quantum chemistry software, etc.

We understand that it is not always clear-cut what contributions should be in the main pymatgen distribution and what
should be as an add-on. It is helpful to always submit a `GitHub Issue
<http://github.com/materialsproject/pymatgen/issues>`_ with your proposed feature to get feedback from the pymatgen
maintainers as well as the broader community. It is also possible that a package that initially was developed as an
add-on subsequently gain enough broad interest and traction such that it is incorporated into the main pymatgen
distribution.

Direct contributions to pymatgen main distribution
--------------------------------------------------

1. Create a free GitHub account (if you don't already have one) and perform the necessary setup (e.g., install SSH
   keys etc.).
2. Fork the pymatgen GitHub repo, i.e., go to the main `pymatgen GitHub repo`_ and click fork to create a copy of the
   pymatgen code base on your own Github account.
3. Install git on your local machine (if you don't already have it).
4. Clone *your forked repo* to your local machine. You will work mostly with your local repo and only publish changes
   when they are ready to be merged::

       git clone git@github.com:YOURNAME/pymatgen.git

   Note that the entire Github repo is fairly large because of the presence of test files, but these are absolutely
   necessary for rigorous testing of the code.
5. It is highly recommended you install all the optional dependencies as well::

      pip install -e '.[dev,optional]'

6. Code (see `Coding Guidelines`_). Commit early and commit often. Keep your code up to date. You need to add the main
   repository to the list of your remotes. Let's name the upstream repo as mpmaster (materialsproject master)::

       git remote add mpmaster git://github.com/materialsproject/pymatgen.git

   Make sure your repository is clean (no uncommitted changes) and is currently on the master branch. If not, commit or
   stash any changes and switch to the master::

      git checkout master

   Then you can pull all the new commits from the main line::

      git pull mpmaster master

   Remember, pull is a combination of the commands fetch and merge, so there may be merge conflicts that may need to be
   manually resolved.
7. Publish your contributions. Assuming that you now have a couple of commits that you would like to contribute to the
   main repository. Please follow the following steps:

   a. If your change is based on a relatively old state of the main repository, then you should probably bring your
      repository up-to-date first to see if the change is not creating any merge conflicts.
   b. Check that everything compiles cleanly and passes all tests.
      The pymatgen repo comes with a complete set of tests for all modules. If
      you have written new modules or methods, you must write tests for the new
      code as well (see `Coding Guidelines`_). Install and run pytest in your
      local repo directory and fix all errors before continuing further.
   c. If you have `pre-commit <https://pre-commit.com/>`_ installed you can use
      the provided :code:`.pre-commit-config.yaml` file to perform automatic style checks
      before publishing your code. The pre-commit hooks can be installed using::

            pre-commit install

   d. If everything is ok, publish the commits to your github repository::

         git push origin master

8. Now that your commit is published, it does not mean that it has already been merged into the main repository. You
   should submit a `pull request <https://github.com/materialsproject/pymatgen/pulls>`_ to the pymatgen main repository.
   A set of linting and unittests will be automatically performed. Fix all lint and test errors identified. Your PR
   will only be merged if these checks passed.

"Work-in-progress" pull requests are encouraged, especially if this is your first time contributing to pymatgen, and
the maintainers will be happy to help or provide code review as necessary. Put "[WIP]" in the title of your
pull request to indicate it's not ready to be merged.

Writing add-ons for pymatgen
----------------------------

Add-on packages can be developed under the pymatgen, pymatgen.analysis, pymatgen.ext and and pymatgen.io
`namespaces <http://packaging.python.org/guides/packaging-namespace-packages/>`_. We ask that most add-ons be developed
under the pymatgen.analysis, pymatgen.ext and and pymatgen.io namespaces and not the root pymatgen namespace. The
pymatgen root namespace is meant for development of broad classes of functionality. If in doubt, please consult with
the pymatgen maintainers. The benefits of writing an add-on for pymatgen are:

* You control the development and distribution of the add-on. You also get full recognition of your work.
* The add-on does not affect the main pymatgen distribution and end users have a choice of whether to install the
  add-on or not via `pip install pymatgen-analysis-addon`.
* Once installed, the add-on functions exactly like a part of pymatgen in that the imports are still via
  `from pymatgen.analysis.addon import *`.
* We will help your add-on gain recognition via our `listing of pymatgen add-ons </addons>`_. We have plans to develop
  this into a full-fledged searchable database of add-ons.

The namespaces provide an important clue what kind of contributions are suitable for add-ons.

* `pymatgen.analysis.*`: A new type of analysis, such as those for a specific application. E.g. superconductors, solar,
  etc. or an entire category of analysis, e.g., machine learning, diffusion, etc.
* `pymatgen.ext.*`: A high-level API access to a new external resource, for example, a new database of crystal
  structure, molecules and/or properties.
* `pymatgen.io.*`: Support for input/output from another code, e.g., some quantum chemistry software.

To help developers write add-ons, we have written a `pymatgen add-on template
<http://github.com/materialsproject/pymatgen-addon-template>`_ with detailed instructions. For a real-world
example using this template, check out Materials Virtual Lab's `pymatgen-analysis-diffusion
<http://pypi.org/project/pymatgen-analysis-diffusion/>`_.

It should be noted that while the pymatgen maintainers will attempt to help developers as far as possible, **we provide
no guarantees whatsoever on the quality or reliability of any code that is not part of the main pymatgen distribution**.
The add-on architecture therefore provides flexibility for broad expansion of scope in pymatgen functionality by the
community by loosening up the tight control in the main repository, which is bottlenecked by the small team maintaining
it.

Coding Guidelines
-----------------

Given that pymatgen is intended to be long-term code base, we adopt very strict
quality control and coding guidelines for all contributions to pymatgen. The
following must be satisfied for your contributions to be accepted into pymatgen.

1. **Unittests** are required for all new modules and methods. The only way to
   minimize code regression is to ensure that all code are well-tested. If the
   maintainer cannot test your code, the contribution will be rejected.
2. **Python PEP 8** `code style <http://www.python.org/dev/peps/pep-0008/>`_.
   We allow a few exceptions when they are well-justified (e.g., Element's
   atomic number is given a variable name of capital Z, in line with accepted
   scientific convention), but generally, PEP 8 must be observed. Code style
   will be automatically checked for all PRs and must pass before any PR is merged.
   To aid you, you can copy the example pre-commit hook into your .git/hooks
   directly. This will automatically run pycodestyle and other linting services
   prior to any commits. At the very least, copy pre-commit to .git/hooks/pre-push.
3. **Python 3**. We only support Python 3.8+.
4. **Documentation** required for all modules, classes and methods. In
   particular, the method docstrings should make clear the arguments expected
   and the return values. For complex algorithms (e.g., an Ewald summation), a
   summary of the algorithm should be provided, and preferably with a link to a
   publication outlining the method in detail.
5. **IDE**. We highly recommend the use of Pycharm. You should also set up
   pycodestyle and turn those on within the IDE setup. This will warn of any
   issues with coding styles. Many code style errors can be done by simply
   selecting the entire code and using the Code->Reformat Code within Pycharm.

For the above, if in doubt, please refer to the core classes in pymatgen for
examples of what is expected.

.. _`pymatgen's Google Groups page`: https://groups.google.com/forum/?fromgroups#!forum/pymatgen/
.. _`pymatgen GitHub repo`: https://github.com/materialsproject/pymatgen
