**Official docs:** http://www.pymatgen.org

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. It currently powers the public Materials Project
(http://www.materialsproject.org), an initiative to make calculated
properties of all known inorganic materials available to materials
researchers. These are some of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive io capabilities to manipulate many VASP
   (http://cms.mpi.univie.ac.at/vasp/) and ABINIT (http://www.abinit.org/)
   input and output files and the crystallographic information file format.
   This includes generating Structure objects from vasp input and output.
   There is also support for Gaussian input files and XYZ file for molecules.
3. Comprehensive tool to generate and view compositional and grand canonical
   phase diagrams.
4. Electronic structure analyses (DOS and Bandstructure).
5. Integration with the Materials Project REST API.

Pymatgen, like all scientific research, will always be a work in progress.
While the development team will always strive to avoid backward incompatible
changes, they are sometimes unavoidable, and tough decisions have to be made
for the long term health of the code.

Pymatgen is free to use. However, we also welcome your help to improve this
library by making your own contributions.  These contributions can be in the
form of additional tools or modules you develop, or even simple things such
as bug reports. Please report any bugs and issues at pymatgen's `Github page
<https://github.com/materialsproject/pymatgen>`_. If you wish to be notified
of pymatgen releases, you may become a member of `pymatgen's Google Groups page
<https://groups.google.com/forum/?fromgroups#!forum/pymatgen/>`_.

Why use pymatgen?
=================

There are many materials analysis codes out there, both commerical and free.
So you might ask - why should I use pymatgen over others? Pymatgen offer
several advantages over other codes out there:

1. **It is (fairly) robust.** Pymatgen is used in the Materials Project. As
   such, the analysis it produces survives rigorous scrutiny every single
   day. Bugs tend to be found and corrected quickly. Furthermore,
   pymatgen uses `CircleCI <https://circleci.com>`_ for continuous
   integration, which ensures that all unittests pass with every commit.
2. **It is well documented.** A fairly comprehensive documentation has been
   written to help you get to grips with it quickly. That means more
   efficient research.
3. **It is open.** That means you are free to use it, and you can also
   contribute to it. It also means that pymatgen is continuously being
   improved. We have a policy of attributing any code you contribute to any
   publication you choose. Contributing to pymatgen means your research
   becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in pymatgen have been
   optimized by vectorizing in numpy. This means that coordinate
   manipulations are extremely fast and are in fact comparable to codes
   written in other languages. Pymatgen also comes with a complete system for
   handling periodic boundary conditions.

Getting pymatgen
================

Before installing pymatgen, you may need to first install a few critical
dependencies manually. Please refer to the official `pymatgen page`_ for
installation details and requirements.

1. Installation has been tested to be most successful with gcc,
   and several external C dependencies have issues with icc. Use gcc where
   possible and do "export CC=gcc" prior to installation.
2. Numpy's distutils is needed to compile the spglib and pyhull
   dependencies. This should be the first thing you install.
3. Pyhull and PyCifRW. The recent versions of pip does not allow the
   installation of externally hosted files. Furthermore,
   there are some issues with easy_install for these extensions. Install
   both these dependencies manually using "pip install <package>
   --allow-external <package> --allow-unverified <package>".

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. The easiest way to
install pymatgen on any system is to use easy_install or pip, as follows::

    easy_install pymatgen

or::

    pip install pymatgen

Some extra functionality (e.g., generation of POTCARs) do require additional
setup (please see the `official pymatgen page <http://pymatgen.org/>`_).

**Note for Windows users**: Given that pymatgen requires several Python C
extensions, it is generally recommended that you install it in a cygwin or
equivalent environment with the necessary compilers.

Developmental version
---------------------

The bleeding edge developmental version is at the pymatgen's `Github repo
<https://github.com/materialsproject/pymatgen>`_. The developmental
version is likely to be more buggy, but may contain new features. The
Github version include test files as well for complete unit testing. After
cloning the source, you can type::

    python setup.py install

or to install the package in developmental mode::

    python setup.py develop

To run the very comprehensive suite of unittests, make sure you have nose
installed and then just type::

    nosetests

in the pymatgen root directory.

Using pymatgen
==============

Please refer to the official `pymatgen page`_ for tutorials and examples.

How to cite pymatgen
====================

If you use pymatgen in your research, please consider citing the following
work:

    Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
    Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
    Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
    Open-Source Python Library for Materials Analysis.* Computational
    Materials Science, 2013, 68, 314-319. `doi:10.1016/j.commatsci.2012.10.028
    <http://dx.doi.org/10.1016/j.commatsci.2012.10.028>`_

In addition, some of pymatgen's functionality is based on scientific advances
/ principles developed by the computational materials scientists in our team.
Please refer to `pymatgen's documentation <http://pymatgen.org/>`_ on how to
cite them.

License
=======

Pymatgen is released under the MIT License. The terms of the license are as
follows::

    The MIT License (MIT)
    Copyright (c) 2011-2012 MIT & LBNL

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


.. _`pymatgen page` : http://www.pymatgen.org