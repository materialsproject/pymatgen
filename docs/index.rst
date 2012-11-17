.. pymatgen documentation master file, created by
   sphinx-quickstart on Tue Nov 15 00:13:52 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/pymatgen.png
   :width: 300 px
   :alt: pymatgen
   :align: center

Introduction
============

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. It currently powers the public Materials Project
(http://www.materialsproject.org), an initiative to make calculated properties
on a large number of materials available to materials researchers and designers.
These are some of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive io capabilities to manipulate many VASP input and output files
   (http://cms.mpi.univie.ac.at/vasp/) and the crystallographic information
   file format. This includes generating Structure objects from vasp input and
   output. There is also support for Gaussian input files and XYZ file for
   molecules.
3. Comprehensive tool to generate and view compositional and grand canonical
   phase diagrams.
4. Electronic structure analyses (DOS and Bandstructure).
5. Integration with the Materials Project REST API.

The pymatgen library is free (as in free beer) to download and to use. However,
we would also like you to help us improve this library by making your own
contributions as well.  These contributions can be in the form of additional
tools or modules you develop, or even simple things such as bug reports. Please
read the Contributing_ section or contact the maintainer of this library
(shyue@mit.edu) to find out how to include your contributions via github or for
bug reports.

Note that pymatgen, like all scientific research, will always be a work in
progress. While the development team will always strive to avoid backward
incompatible changes, they are sometimes unavoidable, and tough decisions have
to be made for the long term health of the code.

The most up-to-date documention is available at our github page
(http://materialsproject.github.com/pymatgen/), where you can also report any
bugs/issues. If you wish to be notified via email of pymatgen releases, you may
become a member of `pymatgen's Google Groups page`_.

   *The code is mightier than the pen.*

Latest Change Log (v2.2.5dev)
-----------------------------

1. Brand new *beta* bond valence analyzer based on a Maximum A Posteriori
   algo using data-mined ICSD data.
2. Speed up and improvements to core classes.
3. Improved structure fitter (credits to Geoffroy Hautier).
4. Brand new entry_tools module (pymatgen.entries.entry_tools).
5. Vastly improved Outcar parser based on reverse parsing that speeds up
   reading of OUTCAR files by orders of magnitude.
6. Miscellaneous bug fixes.

.. toctree::
   :maxdepth: 2

   changelog

Getting pymatgen
================

pymatgen is now in the Python Package Index (`PyPI`_). The version on
PyPI is always the latest stable release that will be hopefully, be relatively
bug-free. If you have  distutils installed, you can just type:

::

   easy_install pymatgen

or

::

   pip install pymatgen

if you have setuptools or pip installed to install pymatgen with most of the
dependencies set up. Otherwise, the latest stable source can be downloaded at
the `PyPI`_ site as well. Note that you may need to install numpy before
installing pymatgen as numpy's distutils is needed to compile the spglib
library.

Alternatively, the bleeding edge developmental version is at the public
pymatgen github repo at
https://github.com/materialsproject/pymatgen/tarball/master. These developmental
versions are likely to be more buggy, but may contain new features. Note that
the GitHub versions include test files as well for complete unittesting.

From the source, you can type::

   python setup.py install

With these basic steps, you should be able to use most of the basic
functionality of pymatgen. However, some extra functionality do require
additional setup. Please see the following sections for further details on the
dependencies needed, where to get them and how to install them.

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
manipulation or analyses. Users are strongly encouraged to explore the
detailed :doc:`usage pages </usage>` (toc given below), and
:doc:`the API docs </modules>`.

.. toctree::
   :maxdepth: 2

   usage

Command line - matgenie.py
--------------------------

To demonstrate the capabilities of pymatgen and to make it easy for users to
quickly use the functionality, pymatgen comes with a set of useful scripts
that utilize the library to perform all kinds of analyses. You can find these
scripts in `scripts directory of pymatgen's github repo
<https://github.com/materialsproject/pymatgen/tree/master/scripts>`_.

Here, we will discuss the most versatile of these scripts,
known as matgenie.py. The typical usage of matgenie.py is::

    matgenie.py {analyze,plot,convert,symm,view} additional_arguments

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
    #vasprun.xml, POSCAR, CSSR.You can force the script to assume certain file
    #types by specifying additional arguments. See matgenie.py convert -h.

    matgenie.py convert input_filename output_filename.

    #Obtain spacegroup information.

    matgenie.py symm -s filename1 filename2

    #Visualize a structure. Requires VTK to be installed.

    matgenie.py view filename

Aliases
-------

From version 2.0.0 of pymatgen, useful aliases for commonly used Objects are
now provided, similar in style to numpy. Supported objects include Element,
Composition, Structure, Molecule, Spin and Orbital. Here are some quick
examples of the core capabilities and objects::

   >>> from pymatgen import Element, Composition, Lattice, Structure
   >>>
   >>> si = Element("Si")
   >>> si.atomic_mass
   28.0855
   >>> si.melting_point
   u'1687 K'
   >>>
   >>> comp = Composition("Fe2O3")
   >>> comp.weight
   159.6882
   >>> comp[Element("Fe")]
   2.0
   >>> comp.get_atomic_fraction(Element("Fe"))
   0.4
   >>>
   >>> structure = Structure(Lattice.cubic(4.2), ["Cs", "Cl"],
   ...                               [[0, 0, 0], [0.5, 0.5, 0.5]])
   >>> structure.volume
   74.088000000000008
   >>> structure[0]
   Periodic Site
   abc : (0.0000, 0.0000, 0.0000)
   element    : Cs
   occupation : 1.00
   >>>
   >>> #Integrated symmetry tools from spglib.
   ... from pymatgen.symmetry.finder import SymmetryFinder
   >>> finder = SymmetryFinder(structure)
   >>> finder.get_spacegroup_symbol()
   'Pm-3m'
   >>>
   >>> #Writing out a POSCAR file for VASP calculations.
   ... from pymatgen.io.vaspio import Poscar
   >>> poscar = Poscar(structure)
   >>> poscar.write_file("POSCAR")

The above illustrates only the most basic capabilities of pymatgen.

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
publications) in pymatgen. Read on to find out about the various ways you can
contribute.

.. toctree::
   :maxdepth: 2

   contributing

API/Reference Docs
==================

The API documentation for pymatgen is provided at the link below.

.. toctree::
   :maxdepth: 1

   modules

The API docs are generated using Sphinx auto-doc and outlines the purpose of all
modules and classes, and the expected argument and returned objects for most
methods.

Citing pymatgen
===============

If you use pymatgen in your research, please consider citing the following
work:

   Shyue Ping Ong, William Davidson Richard, Anubhav Jain, Geoffroy Hautier,
   Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
   Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
   Open-Source Python Library for Materials Analysis.* - Submitted

In addition, some of pymatgen's functionality is based on scientific advances
/ principles developed by the computational materials scientists in our team.
If you use some of these functionality in your research, you may wish to
consider citing the following works:

pymatgen.io.vaspio_set module
-----------------------------

The MIT parameter sets, which are optimized for high-throughput computing, are
outlined the following work:

   A. Jain, G. Hautier, C. Moore, S. P. Ong, C. C. Fischer, T. Mueller,
   K. A. Persson, and G. Ceder. *A high-throughput infrastructure for density
   functional theory calculations.* Computational Materials Science, 2011,
   50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023

pymatgen.phasediagram package
-----------------------------

The phase diagram code, in particular the grand canonical phase diagram
analysis, is based on the work of Ong et al. and are used in following works:

   S. P. Ong, L. Wang, B. Kang, and G. Ceder. *Li-Fe-P-O2 Phase Diagram from
   First Principles Calculations.* Chemistry of Materials, 2008, 20(5),
   1798-1807. doi:10.1021/cm702327g

   S. P. Ong, A. Jain, G. Hautier, B. Kang, and G. Ceder. *Thermal stabilities
   of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first
   principles calculations.* Electrochemistry Communications, 2010, 12(3),
   427-430. doi:10.1016/j.elecom.2010.01.010

pymatgen.entries.compatibility module
-------------------------------------

The compatibility processing, which allows mixing of GGA and GGA+U runs that
have been calculated using the MaterialsProjectVaspInputSet or MITVaspInputSet,
is based on the following work:

   A. Jain, G. Hautier, S. P. Ong, C. Moore, C. C. Fischer, K. A. Persson, and
   G. Ceder. *Formation enthalpies by mixing GGA and GGA + U calculations.*
   Physical Review B, 2011, 84(4), 045115. doi:10.1103/PhysRevB.84.045115

pymatgen.symmetry
-----------------

The symmetry package is based on the excellent spglib developed by Atz Togo. For
more information, please refer to Atz Togo's site at
http://spglib.sourceforge.net/.

License
=======

Pymatgen is released under the MIT License. The terms of the license are as
follows::

   The MIT License (MIT)
   Copyright (c) 2011-2012 MIT & LBNL

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`pymatgen's Google Groups page`: https://groups.google.com/forum/?fromgroups#!forum/pymatgen/
.. _`PyPI` : http://pypi.python.org/pypi/pymatgen