Requirements
============

Required for proper functioning of the code
-------------------------------------------

1. Python 2.7+ required.  New default modules such as json are used, as well as 
   new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively 
   by all core modules.
3. scipy 0.10+ - For interpolation, physical constants and other functions. In 
   particular, scipy.spatial.Delaunay is used for phase diagram construction.

Optional Python Libraries
-------------------------
Optional python libraries that are required if you need certain features

1. matplotlib (highly recommended): For plotting (e.g., Phase Diagrams).
2. PyCifRW (highly recommended): For reading and writing Crystallographic 
   Information Format (CIF) files. Get it from http://pycifrw.berlios.de/.
3. pyspglib 1.2+ (highly recommended): For symmetry finding. Needed if you are
   using the pymatgen.symmetry, pymatgen.transformation and pymatgen.alchemy
   packages. Get it at http://spglib.sourceforge.net/
4. VTK with Python bindings (http://www.vtk.org/): For visualization of crystal 
   structures using the pymatgen.vis package.
5. Atomistic Simulation Environment or ASE : Required for the usage of the 
   adapters in pymatgen.io.aseio between pymatgen's core Structure object and 
   the Atoms object used by ASE. Get it at https://wiki.fysik.dtu.dk/ase/.
6. OpenBabel with Python bindings (http://openbabel.org). Required for the
   usage of the adapters in pymatgen.io.babelio between pymatgen's Molecule
   and OpenBabel's OBMol. Opens up input and output support for the very large
   number of input and output formats supported by OpenBabel.
7. nose - For complete unittesting. This is NOT optional for developers!

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good pythonic alternative exists at 
the moment) required only for certain features.

1. Qhull : Needed for bond length analysis (structure_analyzer.py). The executable 
   qconvex and qvoronoi must be in the path. Get it at http://www.qhull.org/.
2. ffmpeg : Needed for generation of movies (structure_vtk.py).  The executable 
   ffmpeg must be in the path. Get it at http://www.ffmpeg.org.

POTCAR Setup for Users
======================

For the code to generate POTCAR files, it needs to know where the VASP 
pseudopotential files are.  We are not allowed to distribute these under the 
VASP license. The good news is that we have included a setup script to help you along.

If you cloned the repo directly from github, you should have a run_me_first.sh 
file in the root directory of your local repo. Otherwise, you can get it directly 
from our github site at http://github.com/materialsproject/pymatgen. Run the 
shell script and follow the instructions. If you have done it correctly, you 
should get a resources directory with the following directory structure::

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

Setup for Developers (using github)
===================================

1. Clone the repo at http://github.com/materialsproject/pymatgen.
2. Install the necessary python libraries.
3. (Recommended) Add pymatgen to your PYTHONPATH.
4. (Recommended) Copy hooks from the example-hooks directory into the .git/hooks/ 
   directory in your local repo.  

I recommend that you start by reading some of the unittests in the tests 
subdirectory for each package.  The unittests demonstrate the expected behavior 
and functionality of the code.

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

Mac OS X 10.7
~~~~~~~~~~~~~

Typical installation of Xcode with python setup.py install seems to work fine. 
The pre-compiled binary for OSX 10.6 also seems to work.

Matplotlib (tested on v1.10)
----------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

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
	
Spglib (tested on v1.2)
-----------------------

Mac OS X 10.7
~~~~~~~~~~~~~

Download spglib from http://spglib.sourceforge.net/ and then enter the following 
commands:

::

	tar -zxvf spglib-1.1.2.tar.gz
	cd spglib-1.1.2
	./configure
	make
	sudo make install
	cd python/ase
	python setup.py install
	

Qhull (tested on v2012.1)
-------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

Typical installation with make fails with the following error:

	cc1plus: error: unrecognized command line option "-Wno-sign-conversion"

Simply removing "-Wno-sign-conversion" where it appears in the Makefile and then 
doing make followed by make install works fine.

VTK (tested on v5.8.0)
----------------------

Mac OS X 10.7
~~~~~~~~~~~~~

The easiest is to install cmake from
http://cmake.org/cmake/resources/software.html.

Type the following:

::

	cd VTK (this is the directory you expanded VTK into)
	cmake -i (this uses cmake in an interactive manner)

For all options, use the defaults, EXCEPT for BUILD_SHARED_LIBS and 
VTK_WRAP_PYTHON which must be set to ON. You may also need to modify the python 
paths and library paths if they are in non-standard locations.  After the 
CMakeCache.txt file is generated, type:

::

	make (note that this takes a while)
	sudo make install
	
With any luck, you should have vtk with the necessary python wrappers installed.

OpenBabel (tested on v2.3.0)
----------------------------

Mac OS X 10.7
~~~~~~~~~~~~~

openbabel must be compiled with python bindings for integration with pymatgen.
For some reason, openbabel v2.3.1 is harder to compile on Mac OS Lion than I
thought. But I managed to get v2.3.0 to work. Here are the steps that I took to
make it work:

1. Install cmake from http://cmake.org/cmake/resources/software.html.
2. Download openbabel 2.3.0 *source code* from
   http://sourceforge.net/projects/openbabel/files/openbabel/2.3.0/.
3. Download Eigen version 2.0 (newer versions will *not* work) from
   http://eigen.tuxfamily.org/index.php?title=Main_Page
4. Extract your Eigen and openbabel source distributions:

::

   tar -zxvf openbabel-2.3.0.tar.gz
   tar -zxvf eigen2.tar.gz 
   
5. Now you should have two directories. Assuming that your openbabel src is in 
   a directory called "openbabel-2.3.0" and your eigen source is in a directory
   called "eigen2", do the following steps.
   
::
   mv openbabel-2.3.0 ob-src
   mkdir ob-build
   cd ob-build
   cmake -DPYTHON_BINDINGS=ON -DEIGEN2_INCLUDE_DIR=../eigen2 ../ob-src 2>&1 | tee cmake.out
   make -j2
   sudo make install
   
With any luck, you should have openbabel with python bindings installed. You can
test your installation by trying to import openbabel from the python command
line.
