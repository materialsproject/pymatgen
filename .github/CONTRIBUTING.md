# Coding Guidelines

1. **Unittests** are required for all new modules and methods. The only way to
   minimize code regression is to ensure that all code are well-tested. If the
   maintainer cannot test your code, the contribution will be rejected.
2. **Python PEP 8** [code style](http://www.python.org/dev/peps/pep-0008/).
   We allow a few exceptions when they are well-justified (e.g., Element's
   atomic number is given a variable name of capital Z, in line with accepted
   scientific convention), but generally, PEP 8 should be observed.
3. **Python 3**. All code should seamless work with Python 2.7 and Python 3.x.
   See more details below.
4. **Documentation** required for all modules, classes and methods. In
   particular, the method docstrings should make clear the arguments expected
   and the return values. For complex algorithms, a summary of the alogirthm 
   should be provided, and preferably with a link to a publication outlining 
   the method in detail.

If in doubt, please refer to the core classes in pymatgen as well as 
associated unittests for examples of what is expected.

# Coding for Python 3 compatibility

With effect from version 3.0, all pymatgen code must be both Python 2.7+ and 3
compatible. Please read [Python's official guidelines](https://docs.python.org/3/howto/pyporting.html) 
on how to write Python 3.x compatible code. Specifically, we have adopted the
following practices throughout pymatgen.

1. **Unicode-always.** Unless you are absolutely sure you need byte literals
   (rare for pymatgen), always use unicode. In particular, the following should
   be the first line of all pymatgen modules::

        from __future__ import division, print_function, unicode_literals

   Future division means that 1/2 returns a float (0.5),
   which is more intuitive scientifically, instead of 0 (default integer
   division in Python 2). print_function ensures that print() is used instead
   of the print statement. And unicode_literals makes it such that all
   strings are treated as unicode by default. If you need to use bytes,
   those should be marked up explicitly as b'byte literal'.
2. **Use of the six package**. Where necessary, use the six package to handle
   interoperability between Python 2 and 3. Examples include the six.moves
   functions (common ones are zip, filter, map), and six.stringtypes (testing
   for string types, which should be rarely done).
3. **Python-modernize**. Use python-modernize to check your code for any
   potential changes that need to be made.
4. **Unit testing**. The entire pymatgen code base is continuously being
   tested in both Python 2.7 and >=3.4. If your code fails either of the
   tests, you need to fix it.
   
# Collaborative Github Workflow

We use the following workflow (adapted from
http://www.eqqon.com/index.php/Collaborative_Github_Workflow):

1. Create a free GitHub account (if you don't already have one) and perform the
   necessary setup (e.g., install SSH keys etc.).
2. Fork the pymatgen GitHub repo.
3. Install git on your local machine (if you don't already have it).
4. Clone *your forked repo* to your local machine. You will work mostly with
   your local repo and only publish changes when they are ready to be merged:
   ```
   git clone git@github.com:YOURNAME/pymatgen.git
   ```
5. It is highly recommended you install all the optional dependencies as well.
6. Code, commit early and commit often. Keep your code up to date. You need 
   to add the main repository to the list of your remotes as "upstream".
   ```
   git remote add upstream git://github.com/materialsproject/pymatgen.git
   ```
   Make sure your repository is clean (no uncommitted changes) and is currently
   on the master branch. If not, commit or stash any changes and switch to the
   master.
   ```
   git checkout master
   ```
   Then you can pull all the new commits from the main line
   ```
   git pull upstream master
   ```
   Remember, pull is a combination of the commands fetch and merge, so there may
   be merge conflicts to be manually resolved.
7. Publish your contributions. Assuming that you now have a couple of commits
   that you would like to contribute to the main repository. Please follow the
   following steps:
   a. If your change is based on a relatively old state of the main repository,
      then you should probably bring your repository up-to-date first to see if
      the change is not creating any merge conflicts.
   b. Check that everything compiles cleanly and passes all tests.
      The pymatgen repo comes with a complete set of tests for all modules. If
      you have written new modules or methods, you must write tests for the new
      code as well (see `Coding Guidelines`_). Install and run nosetest in your
      local repo directory and fix all errors before continuing further. There
      must be **no errors** for the nosetest.
   c. If everything is ok, publish the commits to your github repository.
   ```
   git push origin master
   ```
8. Now that your commit is published, it does not mean that it has already been
   merged into the main repository. You should issue a pull request to
   pymatgen' maintainers. They will run their own tests and checks, merge if 
   appropriate and release.
