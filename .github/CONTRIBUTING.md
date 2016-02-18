# Coding Guidelines

1. **Unittests** are required for all new modules and methods. The only way to
   minimize code regression is to ensure that all code are well-tested. If the
   maintainer cannot test your code, the contribution will be rejected.
2. **Python PEP 8** `code style <http://www.python.org/dev/peps/pep-0008/>`_.
   We allow a few exceptions when they are well-justified (e.g., Element's
   atomic number is given a variable name of capital Z, in line with accepted
   scientific convention), but generally, PEP 8 must be observed.
3. **Python 3**. All code should seamless work with Python 2.7 and Python 3.x.
   Please read [Python's official guidelines](https://docs.python.org/3/howto/pyporting.html) 
   on how to write Python 3.x compatible code, including the usage of the 
   "six" package. It is recommended that you install the "python-modernize" 
   package and run it before submitting any pull requests.
4. **Documentation** required for all modules, classes and methods. In
   particular, the method docstrings should make clear the arguments expected
   and the return values. For complex algorithms (e.g., an Ewald summation), a
   summary of the alogirthm should be provided, and preferably with a link to a
   publication outlining the method in detail.

If in doubt, please refer to the core classes in pymatgen for examples of 
what is expected.

# A word on coding for Python 3 compatibility

With effect from version 3.0, all pymatgen code must be both Python 2.7+ and 3
compatible. Specifically, we have adopted the following practices throughout
pymatgen.

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
   tested in both Python 2.7 and >=3.3. If your code fails either of the
   tests, you need to fix it.