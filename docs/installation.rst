Requirements
============

Required for proper functioning of the code.
--------------------------------------------

1. Python 2.7+ required.  New default modules such as json are used, as well as new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively by all core modules.
3. scipy 0.9+ - For interpolation, physical constants and other functions. In particular, scipy.spatial.Delaunay is used for phase diagram construction.
4. nose - For complete unittesting. This is NOT optional for developers!

Optional Python Libraries
-------------------------
Optional python libraries that are required if you need certain features

1. matplotlib : For plotting (e.g., Phase Diagrams).
2. PyCifRW (http://prdownload.berlios.de/pycifrw/PyCifRW-3.3.tar.gz) : For reading and writing Crystallographic Information Format (CIF) files [more info](http://pycifrw.berlios.de/)
3. pyspglib (http://spglib.sourceforge.net/) : For symmetry finding [more info](http://spglib.sourceforge.net/)
4. VTK with Python bindings (http://www.vtk.org/): For visualization of crystal structures using the pymatgen.vis package.
5. Atomistic Simulation Environment or ASE (https://wiki.fysik.dtu.dk/ase/): Required for the usage of the adapters in pymatgen.io.aseio between pymatgen's core Structure object and the Atoms object used by ASE. 

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good pythonic alternative exists at the moment) required only for certain features.

1. Qhull (http://www.qhull.org/) : Needed for bond length analysis (structure_analyzer.py).  The executable qconvex and qvoronoi must be in the path.

Basic Setup for Users
=====================

pymatgen is now on PyPI (http://pypi.python.org/pypi/pymatgen).  If you have distutils installed, you can now just type: 

::

	easy_install pymatgen
	
to install pymatgen with most of the dependencies set up. Alternatively, you can download the latest source from https://github.com/CederGroupMIT/pymatgen_repo/downloads and run 

::

	python setup.py install

With these basic steps, you should be able to use most of the pymatgen code. However, some extra functionality do require additional setup, as outlined below.


Generating POTCARs
------------------

For the code to generate POTCAR files, it needs to know where the VASP pseudopotential files are.  We are not allowed to distribute these under the VASP license. The good news is that we have included a setup script to help you along.

1. cd to the root directory of the repo where a file called run_me_first.sh is present.
2. Run the run_me_first.sh file, which will generate a resources directory in a location of your choosing. Please choose a location *outside* of the repo itself.  The script will also write a pymatgen.cfg file in the pymatgen subdir.

Basic Setup for Developers (using github)
=========================================

1. Clone the repo at http://github.com/CederGroupMIT/pymatgen_repo.
2. Install the necessary python libraries.
3. (Recommended) Add pymatgen to your PYTHONPATH.
4. (Recommended) Copy hooks from the example-hooks directory into the .git/hooks/ directory in your local repo.  

I recommend that you start by reading some of the unittests in the tests subdirectory for each package.  The unittests demonstrate the expected behavior and functionality of the code.

Installation of optional libaries
=================================

This section provides a guide for installing various optional libraries used in pymatgen.  Some of the python libraries are rather tricky to build in certain operating systems, especially for users unfamiliar with building C/C++ code. Please feel free to update the instructions based on your experiences.

Scipy (tested on v0.10.1)
-------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

Typical installation of Xcode with python setup.py install seems to work fine. The pre-compiled binary for OSX 10.6 also seems to work.

Matplotlib (tested on v1.10)
----------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

This setup assumes you have the latest version of python (2.7 as of this is written) and numpy already installed. 
You will need to set the compiler flags to build matplotlib from source.

:: 
	
	export CFLAGS="-arch x86_64 -I/usr/X11/include -I/usr/X11/include/freetype2" 
	export LDFLAGS="-arch x86_64 -L/usr/X11/lib" 
	python setup.py build 
	sudo python setup.py install


Solaris 10
~~~~~~~~~~

First install solstudio 12.2. Then put the following code in a shell script and run it.

	#!/bin/bash
	PATH=/opt/solstudio12.2/bin:/usr/ccs/bin:/usr/bin:/usr/sfw/bin:/usr/sbin; export PATH
	ATLAS=None; export ATLAS
	BLAS=/opt/solstudio12.2/lib/libsunperf.so; export BLAS
	LAPACK=/opt/solstudio12.2/lib/libsunmath.so; export LAPACK
	python setup.py build
	python setup.py install
	
Qhull (tested on v2012.1)
-------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

Typical installation with make fails with the following error:

	cc1plus: error: unrecognized command line option "-Wno-sign-conversion"

Simply removing "-Wno-sign-conversion" where it appears in the Makefile and then doing make followed by make install works fine.

VTK (tested on v5.8.0)
----------------------

Mac OS X 10.7
~~~~~~~~~~~~~

The easiest is to install cmake from http://cmake.org/cmake/resources/software.html

Type the following:

::

	cd VTK (this is the directory you expanded VTK into)
	cmake -i (this uses cmake in an interactive manner)

For all options, use the defaults, EXCEPT for BUILD_SHARED_LIBS and VTK_WRAP_PYTHON which must be set to ON. You may also need to modify the python paths and library paths if they are in non-standard locations.  After the CMakeCache.txt file is generated, type:

::

	make (note that this takes a while)
	sudo make install
	
With any luck, you should have vtk with the necessary python wrappers installed.

