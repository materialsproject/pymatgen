.. pymatgen documentation master file, created by
   sphinx-quickstart on Tue Nov 15 00:13:52 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/pymatgen.png
   :width: 300 px
   :alt: pymatgen
   :align: center

.. image:: https://circleci.com/gh/materialsproject/pymatgen/tree/stable.png?circle-token=:circle-token

Introduction
============

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
as bug reports. Please report any bugs and issues at pymatgen's `Github
page`_. If you wish to be notified of pymatgen releases, you may become a
member of `pymatgen's Google Groups page`_.

    *The code is mightier than the pen.*

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

Latest Change Log
=================

v2.9.6
------
1. Patch to allow 1D phase diagrams (essentially finding the lowest energy
   phase).

v2.9.5
------
1. Bug fix for linear assignment, which may sometimes affect Structure
   Matcher results.
2. Minor improvement to the way grand canonical PDs work.

v2.9.4
------
1. Bug fix for Pourbaix Maker (Sai).
2. Streamline use of scratch directories for various calls. Require monty >=
   0.1.2.
3. High accuracy mode for Zeo++ (Bharat Medasani).

v2.9.3
------
1. Bug fix release for printing TransformedStructures from Substitutor (Will
   Richards).
2. Misc improvements in BVAnalyzer, coord_utils and defects (Will Richards,
   David Waroquiers and Bharat Medasani).

v2.9.2
------
1. Bug fix release for DummySpecie, which failed when deserializing from
   json and had bad hash function.

v2.9.1
------
1. Structure/Molecule now supports Pythonic list-like API for replacing and
   removing sites. See :ref:`quick_start` for examples.

v2.9.0
------
1. Updates to support ABINIT 7.6.1 (by Matteo Giantomassi).
2. Vastly improved docs.
3. Refactoring to move commonly used Python utility functions to `Monty
   package <https://pypi.python.org/pypi/monty>`_, which is now a dependency
   for pymatgen.
4. Minor fixes and improvements to DiffusionAnalyzer.
5. Many bug fixes and improvements.

:doc:`Older versions </changelog>`

Getting pymatgen
================

Guided install
--------------

For users who intend to use pymatgen purely as an analysis library (without
developing on it), a user-friendly script has been written to guide users
through the installation process for 64-bit Linux and Mac users. This
installation script requires only basic *Python 2.7+, setuptools,
and a working version of gcc* as prerequisites. Click to download the
`pmg_install.py <_static/pmg_install.py>`_ script. Move the script to an
empty directory and then run::

    python pmg_install.py

Unless you are working in a virtual environment, you will probably need to
run the above command with admin privileges (e.g., sudo). This will install
pymatgen with all *basic dependencies*.

To include more optional dependencies, build the enumlib and bader
executables as well as a step-by-step initial setup for POTCARs and Materials
API usage, run::

    python pmg_install.py -f

The full installation requires a Fortran compiler (ifort or gfortran) to be in
the PATH, as well as X11 (`XQuartz <http://xquartz.macosforge.org/>`_ on Mac)
to be installed for matplotlib.

Stable version
--------------

.. note:: Preparation

    Before installing pymatgen, you may need to first install a few critical
    dependencies manually.

    1. Numpy's distutils is needed to compile the spglib and pyhull
       dependencies. This should be the first thing you install.
    2. Pyhull and PyCifRW. The recent versions of pip does not allow the
       installation of externally hosted files. Furthermore,
       there are some issues with easy_install for these extensions. Install
       both these dependencies manually using "pip install <package>
       --allow-external <package> --allow-unverified <package>".

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. The easiest way to
install pymatgen on any system is to use easy_install or pip, as follows::

    easy_install pymatgen

or::

    pip install pymatgen

Detailed installation instructions for various platforms (Mac and Windows)
are given on this :doc:`page </installation>`.

Some extra functionality (e.g., generation of POTCARs) do require additional
setup. Please see the following sections for further details on the
dependencies needed, where to get them and how to install them.

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

Running unittests
~~~~~~~~~~~~~~~~~

To run the very comprehensive suite of unittests included with the
developmental version, make sure you have nose installed and then just type::

    nosetests

in the pymatgen root directory.

Installation help
-----------------

.. toctree::
   :maxdepth: 1

   installation

Using pymatgen
==============

.. figure:: _static/overview.jpg
   :width: 100%
   :alt: pymatgen overview
   :align: center

   Overview of a typical workflow for pymatgen.

The figure above provides an overview of the functionality in pymatgen. A
typical workflow would involve a user converting data (structure, calculations,
etc.) from various sources (first principles calculations, crystallographic and
molecule input files, Materials Project, etc.) into Python objects using
pymatgen's io packages, which are then used to perform further structure
manipulation or analyses.

.. _quick_start:

Quick start
-----------

Useful aliases for commonly used objects are now provided. Supported objects
include Element, Composition, Structure, Molecule, Spin and Orbital. Here are
some quick examples of the core capabilities and objects:

.. code-block:: pycon

    >>> import pymatgen as mg
    >>>
    >>> si = mg.Element("Si")
    >>> si.atomic_mass
    28.0855
    >>> si.melting_point
    u'1687 K'
    >>>
    >>> comp = mg.Composition("Fe2O3")
    >>> comp.weight
    159.6882
    >>> #Note that Composition conveniently allows strings to be treated just
    >>> #like an Element object.
    >>> comp["Fe"]
    2.0
    >>> comp.get_atomic_fraction("Fe")
    0.4
    >>> lattice = mg.Lattice.cubic(4.2)
    >>> structure = mg.Structure(lattice, ["Cs", "Cl"],
    ...                          [[0, 0, 0], [0.5, 0.5, 0.5]])
    >>> structure.volume
    74.088000000000008
    >>> structure[0]
    PeriodicSite: Cs (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
    >>>
    >>> #Integrated symmetry tools from spglib.
    >>> from pymatgen.symmetry.finder import SymmetryFinder
    >>> finder = SymmetryFinder(structure)
    >>> finder.get_spacegroup_symbol()
    'Pm-3m'
    >>>
    >>> # Convenient IO to various formats. Format is intelligently determined
    >>> # from file name and extension.
    >>> mg.write_structure(structure, "POSCAR")
    >>> mg.write_structure(structure, "CsCl.cif")
    >>>
    >>> # Reading a structure from a file.
    >>> structure = mg.read_structure("POSCAR")
    >>>
    >>> # Reading and writing a molecule from a file. Supports XYZ and
    >>> # Gaussian input and output by default. Support for many other
    >>> # formats via the optional openbabel dependency (if installed).
    >>> methane = mg.read_mol("methane.xyz")
    >>> mg.write_mol(mol, "methane.gjf")
    >>>
    >>> # Pythonic API for editing Structures and Molecules (v2.9.1 onwards)
    >>> # Changing the specie of a site.
    >>> structure[1] = "F"
    >>> print structure
    Structure Summary (Cs1 F1)
    Reduced Formula: CsF
    abc   :   4.200000   4.200000   4.200000
    angles:  90.000000  90.000000  90.000000
    Sites (2)
    1 Cs     0.000000     0.000000     0.000000
    2 F     0.500000     0.500000     0.500000
    >>>
    >>> #Changes species and coordinates (fractional assumed for structures)
    >>> structure[1] = "Cl", [0.51, 0.51, 0.51]
    >>> print structure
    Structure Summary (Cs1 Cl1)
    Reduced Formula: CsCl
    abc   :   4.200000   4.200000   4.200000
    angles:  90.000000  90.000000  90.000000
    Sites (2)
    1 Cs     0.000000     0.000000     0.000000
    2 Cl     0.510000     0.510000     0.510000
    >>>
    >>> # Because structure is like a list, it supports most list-like methods
    >>> # such as sort, reverse, etc.
    >>> structure.reverse()
    >>> print structure
    Structure Summary (Cs1 Cl1)
    Reduced Formula: CsCl
    abc   :   4.200000   4.200000   4.200000
    angles:  90.000000  90.000000  90.000000
    Sites (2)
    1 Cl     0.510000     0.510000     0.510000
    2 Cs     0.000000     0.000000     0.000000
    >>>
    >>> # Molecules function similarly, but with Site and cartesian coords.
    >>> # The following changes the C in CH4 to an N and displaces it by 0.01A
    >>> # in the x-direction.
    >>> methane[0] = "N", [0.01, 0, 0]

The above illustrates only the most basic capabilities of pymatgen.

Examples
--------

A good way to explore the functionality of pymatgen is to look at examples.
Please check out the ipython notebooks at our :doc:`examples page </examples>`.

Usage guide
-----------

Users are also strongly encouraged to explore the :doc:`usage pages </usage>`
(toc given below).

.. toctree::
   :maxdepth: 2

   usage

API documentation
-----------------

For detailed documentation of all modules and classes, please refer to the
:doc:`API docs </modules>`.

matgenie.py - Command line tool
-------------------------------

To demonstrate the capabilities of pymatgen and to make it easy for users to
quickly use the functionality, pymatgen comes with a set of useful scripts
that utilize the library to perform all kinds of analyses. These are
installed to your path by default when you install pymatgen through the
typical installation routes.

Here, we will discuss the most versatile of these scripts, known as
matgenie.py. The typical usage of matgenie.py is::

    matgenie.py {analyze, plotdos, plotchgint, convert, symm, view, compare} additional_arguments

At any time, you can use "matgenie.py --help" or "matgenie.py subcommand
--help" to bring up a useful help message on how to use these subcommands.
Here are a few examples of typical usages::

    #Parses all vasp runs in a directory and display the basic energy
    #information. Saves the data in a file called vasp_data.gz for subsequent
    #reuse.

    matgenie.py analyze .

    #Plot the dos from the vasprun.xml file.

    matgenie.py plotdos vasprun.xml

    #Convert between file formats. The script attempts to intelligently
    #determine the file type. Input file types supported include CIF,
    #vasprun.xml, POSCAR, CSSR. You can force the script to assume certain file
    #types by specifying additional arguments. See matgenie.py convert -h.

    matgenie.py convert input_filename output_filename.

    #Obtain spacegroup information.

    matgenie.py symm -s filename1 filename2

    #Visualize a structure. Requires VTK to be installed.

    matgenie.py view filename

    #Compare two structures for similarity

    matgenie.py compare filename1 filename2

    #Generate a POTCAR with symbols Li_sv O and the PBE functional

    matgenie.py generate --potcar Li_sv O --functional PBE

ipmg - A Custom ipython shell
-----------------------------

From version 2.5.2, A custom ipython shell for pymatgen has been implemented.
Upon installing pymatgen in the usual manner, the "ipmg" script will be
installed. Running ipmg will bring users into a custom ipython environment
where the most commonly used pymatgen objects (see Aliases below) are
automatically loaded into the environment.

Add-ons
-------

Some add-ons are available for pymatgen today:

1. The `pymatgen-db <https://pypi.python.org/pypi/pymatgen-db>`_ add-on
   provides tools to create databases of calculated run data using pymatgen.
2. The `custodian <https://pypi.python.org/pypi/custodian>`_ pacakge provides
   a JIT job management and error correction for calculations, particularly
   VASP calculations.

Contributing
============

Pymatgen is developed by a team of volunteers. It is started by a team
comprising of MIT and Lawrence Berkeley National Laboratory staff to be a
robust toolkit for materials researchers to perform advanced manipulations of
structures and analyses.

For pymatgen to continue to grow in functionality and robustness, we rely on
other volunteers to develop new analyses and report and fix bugs. We welcome
anyone to use our code as-is, but if you could take a few moment to give back
to pymatgen in some small way, it would be greatly appreciated. A benefit of
contributing is that your code will now be used by other researchers who use
pymatgen, and we will include an acknowledgement to you (and any related
publications) in pymatgen.

Reporting bugs
--------------

A simple way that anyone can contribute is simply to report bugs and issues
to the developing team. You can either send an email to the `pymatgen's
Google Groups page`_ or even better, submit an Issue in our `Github page`_.

Developing for pymatgen
-----------------------

Another way to contribute is to submit new code/bugfixes to pymatgen. While
you can always zip your code and email it to the maintainer of pymatgen,
the best way for anyone to develop pymatgen is by adopting the collaborative
Github workflow (see section below).

.. toctree::
   :maxdepth: 2

   contributing

How to cite pymatgen
====================

If you use pymatgen in your research, please consider citing the following
work:

    Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
    Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
    Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
    Open-Source Python Library for Materials Analysis.* Computational
    Materials Science, 2013, 68, 314–319. `doi:10.1016/j.commatsci.2012.10.028
    <http://dx.doi.org/10.1016/j.commatsci.2012.10.028>`_

In addition, some of pymatgen's functionality is based on scientific advances
/ principles developed by various scientists. Please refer to the
:doc:`references page </references>` for citation info.

License
=======

Pymatgen is released under the MIT License. The terms of the license are as
follows::

    The MIT License (MIT)
    Copyright (c) 2011-2012 MIT & LBNL

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software")
    , to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`pymatgen's Google Groups page`: https://groups.google.com/forum/?fromgroups#!forum/pymatgen/
.. _`PyPI` : http://pypi.python.org/pypi/pymatgen
.. _`Github page`: https://github.com/materialsproject/pymatgen/issues
