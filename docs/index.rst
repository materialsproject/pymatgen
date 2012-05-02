.. pymatgen documentation master file, created by
   sphinx-quickstart on Tue Nov 15 00:13:52 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Introduction
============
   
Pymatgen (python materials genomics) is the python library that powers the 
Materials Project (http://www.materialsproject.org). These are some of the main 
features:

1. Highly flexible classes for the representation of Element, Site, Molecule, 
   Structure objects.
2. Extensive io capabilities to manipulate many VASP input and output files 
   (http://cms.mpi.univie.ac.at/vasp/) and the crystallographic information file 
   format. This includes generating Structure objects from vasp input and
   output. There is also support for Gaussian input files and XYZ file for
   molecules.
3. Comprehensive tool to generate and view compositional and grand canonical phase 
   diagrams.
4. Electronic structure analyses (DOS and Bandstructure).

The public version of pymatgen is free (as in free beer) to download and to use. 
However, we would also like you to help us improve this library by making your 
own contributions as well.  These contributions can be in the form of
additional tools or modules you develop, or even simple things such as bug
reports. Please read the Contributing_ section or contact the maintainer of
this library (shyue@mit.edu) to find out how to include your contributions
via github or for bug reports.

Note that pymatgen, like all scientific research, will always be a work in
progress. While the development team will always strive to avoid backward 
incompatible changes, they are sometimes unavoidable, and tough decisions have 
to be made for the long term health of the code.

The most up-to-date documention is available at our github page
(http://materialsproject.github.com/pymatgen/), where you can also report any
bugs/issues. If you wish to be notified via email of pymatgen releases, you may
become a member of `pymatgen's Google Groups page`_.

   *The code is mightier than the pen.*

Latest Change Log (v1.9.0)
--------------------------

1. Completely new json encoder and decoder that support serialization of almost
   all pymatgen objects. See `Using pymatgen`_ section for more details.
2. Simplification to Borg API utilizing the new json API.
3. Bandstructure classes now support spin-polarized runs.

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
   
to install pymatgen with most of the dependencies set up. Otherwise, the latest
stable source can be downloaded at the `PyPI`_ site as well.

Alternatively, the bleeding edge developmental version is at the public
pymatgen github repo at 
https://github.com/materialsproject/pymatgen/tarball/master. These developmental
versions are likely to be more buggy, but may contain new features. Note that
the github versions include test files as well for complete unittesting.

From the source, you can type:

::

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

.. toctree::
   :maxdepth: 2 
   
   usage

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

Some of pymatgen's functionality is based on scientific advances / principles
developed by the computational materials scientists in our team. If you 
use some of these functionality in your research, you may wish to consider
citing the following works:

pymatgen.io.vaspio_set module
-----------------------------

The parameter sets, which are optimized for high-throughput computing, are 
outlined the following work:
      
   A. Jain, G. Hautier, C. Moore, S. P. Ong, C. C. Fischer, T. Mueller, 
   K. A. Persson, and G. Ceder. A high-throughput infrastructure for density 
   functional theory calculations. Computational Materials Science, 2011, 
   50(8), 2295-2310. doi:10.1016/j.commatsci.2011.02.023
      
pymatgen.phasediagram package
-----------------------------

The phase diagram code, in particular the grand canonical phase diagram
analysis, is based on the work of Ong et al. and are used in following works:

   S. P. Ong, L. Wang, B. Kang, and G. Ceder. Li-Fe-P-O2 Phase Diagram from 
   First Principles Calculations. Chemistry of Materials, 2008, 20(5), 1798-1807.
   doi:10.1021/cm702327g
      
   S. P. Ong, A. Jain, G. Hautier, B. Kang, and G. Ceder. Thermal stabilities 
   of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first 
   principles calculations. Electrochemistry Communications, 2010, 12(3), 
   427-430. doi:10.1016/j.elecom.2010.01.010

pymatgen.entries.compatibility module
-------------------------------------

The compatibility processing, which allows mixing of GGA and GGA+U runs that 
have been calculated using the MaterialsProjectVaspInputSet or MITVaspInputSet,
is based on the following work:
      
   A. Jain, G. Hautier, S. P. Ong, C. Moore, C. C. Fischer, K. A. Persson, and 
   G. Ceder. Formation enthalpies by mixing GGA and GGA + U calculations. 
   Physical Review B, 2011, 84(4), 045115. doi:10.1103/PhysRevB.84.045115

pymatgen.symmetry
-----------------

The symmetry package is based on the excellent spglib developed by Atz Togo. For
more information, please refer to Atz Togo's site at
http://spglib.sourceforge.net/.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`pymatgen's Google Groups page`: https://groups.google.com/forum/?fromgroups#!forum/pymatgen/
.. _`PyPI` : http://pypi.python.org/pypi/pymatgen