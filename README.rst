.. image:: https://travis-ci.org/materialsproject/pymatgen.png

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. It currently powers the public Materials Project
(http://www.materialsproject.org), an initiative to make calculated
properties of all known inorganic materials available to materials
researchers. These are some of the main features:

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
   such, the analysis it produces survives rigourous scrutiny every single
   day. Bugs tend to be found and corrected quickly.
2. **It is well documented.** A fairly comprehensive documentation has been
   written to help you get to grips with it quickly. That means more
   efficient research.
3. **It is open.** That means you are free to use it, and you can also
   contribute to it. It also means that pymatgen is continuously being
   improved. We have a policy of attributing any code you contribute to any
   publication you choose. Contributing to pymatgen means your research
   becomes more visible, which translates to greater impact.

Getting pymatgen
================

Stable version
--------------

The version at the Python Package Index (PyPI) is always the latest stable
release that will be hopefully, be relatively bug-free. The easiest way to
install pymatgen on any system is to use easy_install or pip, as follows::

    easy_install pymatgen

or::

    pip install pymatgen

Some extra functionality (e.g., generation of POTCARs) do require additional
setup (please see `pymatgen's documentation
<http://packages.python.org/pymatgen>`_).

**Note**: You may need to install numpy before installing pymatgen as numpy's
distutils is needed to compile the spglib and pyhull dependencies.

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

The docs for the developmental version are available at pymatgen's `Github
pages <http://materialsproject.github.com/pymatgen/>`_.

Requirements
============

All required dependencies should be automatically taken care of if you
install pymatgen using easy_install or pip. Otherwise, these packages should
be available on `PyPI <http://pypi.python.org>`_.

1. Python 2.7+ required. New default modules such as json are used, as well as
   new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively
   by all core modules.
3. pyhull 1.3.6+: For generation of phase diagrams.
4. PyCifRW 3.3+: For reading and writing Crystallographic Information Format
   (CIF) files.
5. requests 1.0+: For the high-level interface to the Materials API.

Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

1. scipy 0.10+ (highly recommended): For use in Gaussian smearing.
2. matplotlib 1.1+ (highly recommended): For plotting (e.g., Phase Diagrams).
3. VTK with Python bindings 5.8+ (http://www.vtk.org/): For visualization of
   crystal structures using the pymatgen.vis package.
4. Atomistic Simulation Environment or ASE 3.6+: Required for the usage of the
   adapters in pymatgen.io.aseio between pymatgen's core Structure object and
   the Atoms object used by ASE. Get it at https://wiki.fysik.dtu.dk/ase/.
5. OpenBabel with Python bindings (http://openbabel.org): Required for the
   usage of the adapters in pymatgen.io.babelio between pymatgen's Molecule
   and OpenBabel's OBMol. Opens up input and output support for the very large
   number of input and output formats supported by OpenBabel.
6. nose - For complete unittesting.

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good python alternative exists at
the moment) required only for certain features:

1. ffmpeg: For generation of movies in structure_vtk.py. The executable ffmpeg
   must be in the path. Get it at http://www.ffmpeg.org.
2. enum: For the use of EnumerateStructureTransformation and the
   pymatgen.command_line.enumlib_caller module. This library by Gus Hart
   provides a robust way to enumerate derivative structures. It can be used to
   completely enumerate all symmetrically distinct ordered structures of
   disordered structures via the EnumerateStructureTransformation. The
   multienum.x and makestr.x executables must be in the path. Get it at
   http://enum.sourceforge.org and follow the instructions to compile
   multienum.x and makestr.x.

Using pymatgen
==============

.. figure:: http://packages.python.org/pymatgen/images/overview.jpg
   :width: 70%
   :alt: pymatgen overview
   :align: center

The figure above provides an overview of the functionality in pymatgen. A
typical workflow would involve a user converting data (structure, calculations,
etc.) from various sources (first principles calculations, crystallographic and
molecule input files, Materials Project, etc.) into Python objects using
pymatgen's io packages, which are then used to perform further structure
manipulation or analyses.

Command line - matgenie.py
--------------------------

To demonstrate the capabilities of pymatgen and to make it easy for users to
quickly use the functionality, pymatgen comes with a set of useful scripts
that utilize the library to perform all kinds of analyses. You can find these
scripts in `scripts directory of pymatgen's github repo
<https://github.com/materialsproject/pymatgen/tree/master/scripts>`_.

Here, we will discuss the most versatile of these scripts,
known as matgenie.py. The typical usage of matgenie.py is::

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


Basic usage
-----------

Useful aliases for commonly used objects are provided. Supported objects
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
    >>> comp[mg.Element("Fe")]
    2.0
    >>> comp.get_atomic_fraction(mg.Element("Fe"))
    0.4
    >>> lattice = mg.Lattice.cubic(4.2)
    >>> structure = mg.Structure(lattice, ["Cs", "Cl"],
    ...                       [[0, 0, 0], [0.5, 0.5, 0.5]])
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
    >>> #Writing out a POSCAR file for VASP calculations.
    >>> poscar = Poscar(structure)
    >>> mg.write_structure(structure, "POSCAR")
    >>>
    >>> #Reading a structure from a file. Supported files include CIF, POSCAR, etc.
    >>> structure = mg.read_structure("POSCAR")

Advanced Usage
--------------

Users are strongly encouraged to explore the detailed `usage pages
<http://packages.python.org/pymatgen/usage.html>`_ and `api docs
<http://packages.python.org/pymatgen/modules.html>`_.

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
Please refer to `pymatgen's documentation
<http://packages.python.org/pymatgen>`_ on how to cite them.
