Requirements
============

All required dependencies should be automatically taken care of if you
install pymatgen using easy_install or pip. Otherwise, these packages should
be available on `PyPI <http://pypi.python.org>`_.

Required for proper functioning of the code
-------------------------------------------

1. Python 2.7+ required.  New default modules such as json are used, as well as
   new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively
   by all core modules.
3. pyhull 1.1+: For generation of phase diagrams.
4. PyCifRW: For reading and writing Crystallographic Information Format (CIF)
   files.

Optional Python Libraries
-------------------------
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
6. nose - For complete unittesting. This is **not optional for developers**!

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good pythonic alternative exists at
the moment) required only for certain features.

1. ffmpeg : Needed for generation of movies (structure_vtk.py).  The executable
   ffmpeg must be in the path. Get it at http://www.ffmpeg.org.
2. enum : Needed for the use of EnumerateStructureTransformation and the
   pymatgen.command_line.enumlib_caller module. This library by Gus Hart
   provides a robust way to enumerate derivative structures. It can be used to
   completely enumerate all symmetrically distinct ordered structures of
   disordered structures via the EnumerateStructureTransformation. The
   multienum.x and makestr.x executables must be in the path. Get it at
   http://enum.sourceforge.org and follow the instructions to compile
   multienum.x and makestr.x.

POTCAR Setup for Users
======================

For the code to generate POTCAR files, it needs to know where the VASP
pseudopotential files are.  We are not allowed to distribute these under the
VASP license. The good news is that we have included a setup script to help you
along.

If you cloned the repo directly from GitHub, you should have a run_me_first.sh
file in the root directory of your local repo. Otherwise, you can get it
directly from our github site at http://github.com/materialsproject/pymatgen.
Run the shell script and follow the instructions. If you have done it
correctly, you should get a resources directory with the following directory
structure::

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

Alternatively, you can setup the above directly structure manually and set the
VASP_PSP_DIR environment variable accordingly.

Setup for Developers (using GitHub)
===================================

1. Clone the repo at http://github.com/materialsproject/pymatgen.

2. In your root pymatgen repo directory, type (you may need to do this with root
   privileges)::

      python setup.py develop

3. Install any missing python libraries that are necessary.

I recommend that you start by reading some of the unittests in the tests
subdirectory for each package. The unittests demonstrate the expected behavior
and functionality of the code.

Please read up on pymatgen's :doc:`coding guidelines </contributing>` before
you start coding. It will make integration much easier.

Installation tips for optional libaries
=======================================

This section provides a guide for installing various optional libraries used in
pymatgen.  Some of the python libraries are rather tricky to build in certain
operating systems, especially for users unfamiliar with building C/C++ code.
Please feel free to send in suggestions to update the instructions based on
your experiences. In all the instructions, it is assumed that you have standard
gcc and other compilers (e.g., Xcode on Macs) already installed.

Scipy (tested on v0.10.1)
-------------------------

Mac OS X 10.7 - 10.8
~~~~~~~~~~~~~~~~~~~~

Typical installation of Xcode with python setup.py install seems to work fine.
The pre-compiled binary for OSX 10.6 also seems to work.

Matplotlib (tested on v1.10)
----------------------------

Mac OS X 10.7 - 10.8
~~~~~~~~~~~~~~~~~~~~

This setup assumes you have the latest version of python (2.7 as of this is written)
and numpy already installed. You will need to set the compiler flags to build
matplotlib from source.

::

	export CFLAGS="-arch x86_64 -I/usr/X11/include -I/usr/X11/include/freetype2"
	export LDFLAGS="-arch x86_64 -L/usr/X11/lib"
	python setup.py build
	sudo python setup.py install

Solaris 10
~~~~~~~~~~

First install solstudio 12.2. Then put the following code in a shell script and
run it.

::

	#!/bin/bash
	PATH=/opt/solstudio12.2/bin:/usr/ccs/bin:/usr/bin:/usr/sfw/bin:/usr/sbin; export PATH
	ATLAS=None; export ATLAS
	BLAS=/opt/solstudio12.2/lib/libsunperf.so; export BLAS
	LAPACK=/opt/solstudio12.2/lib/libsunmath.so; export LAPACK
	python setup.py build
	python setup.py install

VTK (tested on v5.10.0)
-----------------------

Mac OS X 10.7 and 10.8
~~~~~~~~~~~~~~~~~~~~~~

The easiest is to install cmake from
http://cmake.org/cmake/resources/software.html.

Type the following:

::

	cd VTK (this is the directory you expanded VTK into)
	cmake -i (this uses cmake in an interactive manner)

For all options, use the defaults, EXCEPT for BUILD_SHARED_LIBS and
VTK_WRAP_PYTHON which must be set to ON. You may also need to modify the python
paths and library paths if they are in non-standard locations. For example, if
you have installed the official version of Python instead of using the
Mac-provided version, you will probably need to edit the CMakeCache Python
links. Example configuration for Python 2.7 is given below (only variables that
need to be modified are shown):

::

   //Path to a program.
   PYTHON_EXECUTABLE:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/bin/python

   //Path to a file.
   PYTHON_INCLUDE_DIR:PATH=/Library/Frameworks/Python.framework/Versions/2.7/Headers

   //Path to a library.
   PYTHON_LIBRARY:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib

   //Also delete the prefix settings for python, which typically links to the Mac python.

After the CMakeCache.txt file is generated, type:

::

	make
	sudo make install

With any luck, you should have vtk with the necessary python wrappers installed.

OpenBabel (tested on v2.3.2)
----------------------------

Mac OS X 10.7 - 10.8
~~~~~~~~~~~~~~~~~~~~

openbabel must be compiled with python bindings for integration with pymatgen.
Here are the steps that I took to make it work:

1. Install cmake from http://cmake.org/cmake/resources/software.html.
2. Download openbabel 2.3.2 *source code* from
   http://sourceforge.net/projects/openbabel/files/openbabel/2.3.0/.
3. Download Eigen version 3.0 from
   http://eigen.tuxfamily.org/index.php?title=Main_Page
4. Install pkg-config (easiest way is to install homebrew and do the following:

::
    brew install pkg-config

5. Extract your Eigen and openbabel source distributions:

::

   tar -zxvf openbabel-2.3.2.tar.gz
   tar -zxvf eigen3.tar.gz

5. Now you should have two directories. Assuming that your openbabel src is in
   a directory called "openbabel-2.3.2" and your eigen source is in a directory
   called "eigen3", do the following steps.

::

   mv openbabel-2.3.2 ob-src
   mkdir ob-build
   cd ob-build
   cmake -DPYTHON_BINDINGS=ON -DEIGEN3_INCLUDE_DIR=../eigen3 ../ob-src 2>&1 |
    tee cmake.out
   make -j2
   sudo make install

With any luck, you should have openbabel with python bindings installed. You can
test your installation by trying to import openbabel from the python command
line.

Enumlib (tested as of version of Jul 2012)
------------------------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

There does not seem to be any issues with installation as per the instructions
given by the author. For convenience, the steps are reproduced here:

::

   tar -zxvf enum.tar.gz

   #Compile the symmetry library. Go to the celib/trunk directory:
   cd celib/trunk

   #Set an environment variable to identify your fortran compiler
   export F90=gfortran

   make

   Next, make the enumeration library
   cd ../../enumlib/trunk
   make

   # Make the necessary standalone executables
   make multienum.x
   make makestr.x

After doing the above, make sure that the multienum.x and makestr.x executables
are available in your path.