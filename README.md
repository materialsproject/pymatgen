## Introduction ##
[![Build Status](https://travis-ci.org/materialsproject/pymatgen.png)](https://travis-ci.org/materialsproject/pymatgen)

Pymatgen (**py**thon **mat**erials **gen**omics) is the python library that
powers the Materials Project (http://www.materialsproject.org). These are some
of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive io capabilities to manipulate many VASP input and output files
   (http://cms.mpi.univie.ac.at/vasp/) and the crystallographic information file
   format. This includes generating Structure objects from vasp input and
   output. There is also support for Gaussian input files and XYZ file for
   molecules.
3. Comprehensive tool to generate and view compositional and grand canonical
   phase diagrams.
4. Electronic structure analyses (DOS and Bandstructure).

The public version of pymatgen is free (as in free beer) to download and to use.
However, we would also like you to help us improve this library by making your
own contributions as well.  These contributions can be in the form of
additional tools or modules you develop, or even simple things such as bug
reports. Please contact the maintainer of this library (shyue@mit.edu) to find
out how to include your contributions via github or for bug reports.

Note that pymatgen, like all scientific research, will always be a work in
progress. While the development team will always strive to avoid backward
incompatible changes, they are sometimes unavoidable, and tough decisions have
to be made for the long term health of the code.

For documentation and usage guide, please refer to the latest documentation at
our github page (http://materialsproject.github.com/pymatgen/). If you wish to
be notified via email of pymatgen releases, you may become a member of
pymatgen's Google Groups page
(https://groups.google.com/forum/?fromgroups#!forum/pymatgen/).

## Requirements ##

### Required for proper functioning of the code ###

1. Python 2.7+ required.  New default modules such as json are used, as well as
   new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively
   by all core modules.
3. pyhull 1.3.6+: For generation of phase diagrams.
4. PyCifRW 3.3+: For reading and writing Crystallographic Information Format
   (CIF) files.
5. requests 0.14+: For the high-level interface to the Materials API.

### Optional Python Libraries ###

Optional python libraries that are required if you need certain features

1. scipy 0.10+ (highly recommended) - For use in Gaussian smearing.
2. matplotlib (highly recommended): For plotting (e.g., Phase Diagrams).
3. VTK with Python bindings (http://www.vtk.org/): For visualization of crystal
   structures using the pymatgen.vis package.
4. Atomistic Simulation Environment or ASE : Required for the usage of the
   adapters in pymatgen.io.aseio between pymatgen's core Structure object and
   the Atoms object used by ASE. Get it at https://wiki.fysik.dtu.dk/ase/.
5. OpenBabel with Python bindings (http://openbabel.org). Required for the
   usage of the adapters in pymatgen.io.babelio between pymatgen's Molecule
   and OpenBabel's OBMol. Opens up input and output support for the very large
   number of input and output formats supported by OpenBabel.
6. nose - For complete unittesting. This is NOT optional for developers!

### Optional non-Python programs ###

Optional non-python libraries (because no good pythonic alternative exists at
the moment) required only for certain features.

1. [ffmpeg](http://www.http://ffmpeg.org//) : Needed for generation of movies
   (structure_vtk.py).  The executable ffmpeg must be in the path.
2. [enum](http://enum.sourceforge.org) : Needed for the use of
   EnumerateStructureTransformation and the pymatgen.command_line.enumlib_caller
   module. This library by Gus Hart provides a robust way to enumerate
   derivative structures. It can be used to completely enumerate all
   symmetrically distinct ordered structures of disordered structures via the
   EnumerateStructureTransformation. The multienum.x and makestr.x executables
   must be in the path.

## Basic Setup for Non-developers ##

If you are using pymatgen purely as a library and do not intend to contribute
code, you may install pymatgen either easy_install or python setup.py.

If you have easy_install or pip in your Python setup, the simplest way to get
the latest stable release of pymatgen is to do:

	easy_install pymatgen

or

    pip install pymatgen

If you don't have easy_install / pip, or you prefer to install the latest
development version of pymatgen, you can download the tarball and then do the
following:

	tar -zxvf pymatgen.tar.gz
	cd pymatgen
	python setup.py install

You may need to run the above with root privileges on your machine. In addition,
you may need to install additional python libraries and dependencies.

With these two basic steps, you should be able to use most of the pymatgen code.
I recommend that you start by reading some of the unittests in the tests
subdirectory for each package.  The unittests demonstrate the expected behavior
and functionality of the code.

However, some extra functionality do require additional setup, as outlined
below.

## Setup for Developers ##

There are two categories of developers.  General developers should follow the
procedures outlined in the pymatgen documentation on collaborative Github
workflow to fork a copy of pymatgen to their own Github accounts and cloning it
into their local machine.

Core developers who have write access to the main Github repo may clone the
pymatgen repo directly.

For both kinds of developers, it is recommended that after you clone the repo,
you run:

	python setup.py develop

which will install pymatgen in development mode and install some of the
necessary dependencies.

## Generating POTCARs ##

For pymatgen to generate POTCAR files, it needs to know where the VASP
pseudopotential files are.  We are not allowed to distribute these under the
VASP license. The good news is that we have included a setup script to help you
along.

After installation, do:

    potcar_setup.py

and follow the instructions. If you have done it correctly, you should get a
resources directory with the following directory structure:

	- psp_resources
	|- POT_GGA_PAW_PBE
	||- POTCAR.Ac_s.gz
	||- POTCAR.Ac.gz
	||- POTCAR.Ag.gz
	...
	|- POT_GGA_PAW_PW91
	...

After generating the resources directory, you should add a VASP_PSP_DIR
environment variable pointing to the generated directory and you should then be
able to generate POTCARs.
