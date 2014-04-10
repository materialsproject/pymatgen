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
6. monty 0.1.1+: For some common complementary functions,
   design patterns (e.g., singleton) and decorators to the Python
   standard library.

Optional dependencies
---------------------

Optional libraries that are required if you need certain features:

1. scipy 0.10+ (highly recommended): For use in Gaussian smearing and faster
   Phase Diagrams.
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
6. nose - For unittesting. Not optional for developers.

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good python alternative exists at
the moment) required only for certain features:

1. ffmpeg: For generation of movies in structure_vtk.py. The executable ffmpeg
   must be in the path. Get it at http://www.ffmpeg.org.
2. enum: For the use of
   :class:`pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation`
   and :mod:`pymatgen.command_line.enumlib_caller` module. This library by Gus
   Hart provides a robust way to enumerate derivative structures. It can be
   used to completely enumerate all symmetrically distinct ordered structures
   of disordered structures via EnumerateStructureTransformation. Many other
   advanced transformations (e.g., MagOrderingTransformation) use
   EnumerateStructureTransformation. The multienum.x and makestr.x
   executables must be in the path. Get it at http://enum.sourceforge.net and
   follow the instructions to compile multienum.x and makestr.x.
3. bader: For use with :class:`pymatgen.command_line.bader.BaderAnalysis`.
   This library by Henkelmann et al. provides a robust way to calculate the
   Bader analysis from a CHGCAR. The bader executable must be in the path.
   Get it at http://theory.cm.utexas.edu/bader.
4. gulp: For use with :mod:`pymatgen.command_line.gulp_caller`,
   which is in turn used extensively by :mod:`pymatgen.analysis.defects` to
   compute empirical defect energies.
5. aconvasp: For use with the :mod:`pymatgen.command_line.aconvasp_caller`.

Detailed installation instructions
==================================

Mac OSX (tested on 10.6-10.9)
-----------------------------

For Macs, the initial installation steps can be a bit complicated because
pymatgen does require a number of extensions to be built. Here are some
recommended step-by-step instructions for a minimal setup necessary for
pymatgen usage.

1. Download and install the basic compilers needed:
    a. Xcode - This provides the gcc compiler. Get it from the App Store.
        i.  *OSX < 10.9*. After installation, start XCode,
            go to Preferences->Downloads and install the Command Line Tools.
            You may need to quit Xcode and reopen if the Command Line Tools
            option does not appear.
        ii. *OSX 10.9*. The command line tools for OSX Mavericks is no longer
            provided as an option under Xcode downloads. To install command
            line tools, type the following in a terminal::

               xcode-select --install

    b. Gfortran 4.6.2+ - Get an installer at
       http://gcc.gnu.org/wiki/GFortranBinaries#MacOS.
2. It is recommended that you install the latest copy of Python 2.7+ (not 3+),
   even though your Mac should already have a compatible version. This makes it
   easier for future upgrades and minimizes issues. Get it from the `Python
   home page <http://www.python.org>`_ and install.
3. Ensure that your terminal is running the correct version by typing::

    which python

   You should see something like
   “/Library/Frameworks/Python.framework/Versions/2.7/bin/python”. If you don’t
   get this (e.g., if you get /usr/bin/python), you may need to change your
   PATH.

4. Python setuptools make it easier to install subsequent programs via
   “easy_install”. If you want to, you can install pip as well using “sudo
   easy_install pip”. Pip has several advantages over easy_install. In a
   terminal, run::

    curl -o setuptools-0.6c11-py2.7.egg http://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg
    sudo sh setuptools-0.6c11-py2.7.egg
    sudo easy_install pip #optional

5. Install numpy and a few other dependencies::

    sudo pip install numpy
    sudo pip install pycifrw
    sudo pip install pyyaml

6. Install pymatgen, either in development mode or via pip.

Windows (tested on Win7 64-bit)
-------------------------------

The installation of pymatgen on Windows machines is particularly tricky,
especially for the non-developer, because many of the tools that come bundled
with Unix-based systems (such as gcc and python itself!) are not part of a
standard Windows installation.

The instructions below is a tested installation procedure for getting python
and pymatgen working on a Windows system via `Cygwin
<http://www.cygwin.com>`_, which is the easiest procedure I know. If anyone
has a step-by-step guide for a native installation without cygwin,
please send me the details.

1. Download and install `Cygwin`_. When you get to the part that asks for the
   packages to install, make sure that the following are selected:

    a. Devel - gcc4 (not gcc, which is gcc 3.*), git
    b. Libs - libmpfr4
    c. Python - python, python-numpy, python-setuptools
    d. Net - openssh

   In all cases, make sure that both binary and src is selected.
2. Start the Cygwin terminal.
3. easy_install pip (this makes it much easier to manage packages)::

    easy_install pip

4. Install some required packages which seem to have issues when installed as
   part of the pymatgen setup.py process::

    pip install pycifrw
    pip install pyyaml

5. Install pymatgen either using pip or the Github developer procedures
   below.
6. Test your installation by entering the python interactive prompt and doing
   a "import pymatgen as mg".

Linux
-----

If you are using a Linux system, it is generally assumed that you will have
python, numpy and the standard compilers already on your system. Standard
easy_install or pip install should work automatically. Even if there are some
minor compilation error messages, I generally assume Linux users are usually
able to diagnose and solve those. For users of Ubuntu, most of the dependencies
(including the optional ones) are most easily installed using apt-get.

POTCAR Setup
============

For the code to generate POTCAR files, it needs to know where the VASP
pseudopotential files are.  We are not allowed to distribute these under the
VASP license. The good news is that we have included a setup script to help you
along.

After installation, do::

    potcar_setup.py

and follow the instructions. If you have done it correctly, you should get a
resources directory with the following directory structure::

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

Installation tips for optional libraries
========================================

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

Openbabel must be compiled with python bindings for integration with pymatgen.
Here are the steps that I took to make it work:

1. Install cmake from http://cmake.org/cmake/resources/software.html.

2. Install pcre-8.33 from
   ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.33.tar.gz.

3. Install pkg-config-0.28 using MacPorts or from
   http://pkgconfig.freedesktop.org/releases/pkg-config-0.28.tar.gz.

4. Install SWIG from
   http://prdownloads.sourceforge.net/swig/swig-2.0.10.tar.gz.

5. Download openbabel 2.3.2 *source code* from
   https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/.

6. Download Eigen version 3.1.2 from
   http://bitbucket.org/eigen/eigen/get/3.1.2.tar.gz.

7. Extract your Eigen and openbabel source distributions::

    tar -zxvf openbabel-2.3.2.tar.gz
    tar -zxvf eigen3.tar.gz

8. Now you should have two directories. Assuming that your openbabel src is in
   a directory called "openbabel-2.3.2" and your eigen source is in a directory
   called "eigen3", do the following steps::

    mv openbabel-2.3.2 ob-src
    cd ob-src/scripts/python; rm openbabel.py openbabel-python.cpp; cd ../../..

9. Edit ob-src/scripts/CMakeLists.txt, jump to line 70, change “eigen2_define”
   to “eigen_define”.

10. Let's create a build directory::

        mkdir ob-build
        cd ob-build
        cmake -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DEIGEN3_INCLUDE_DIR=../eigen3 ../ob-src 2>&1 | tee cmake.out

11. Before proceeding further, similar to the VTK installation process in the
    previous section, you may also need to modify the CMakeCache.txt
    file by hand if your python paths and library paths if they are in
    non-standard locations. For example, if you have installed the official
    version of Python instead of using the Mac-provided version,
    you will probably need to edit the CMakeCache Python links. Example
    configuration for Python 2.7 is given below (only variables that need to
    be modified are shown)::

        //Path to a program.
        PYTHON_EXECUTABLE:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/bin/python

        //Path to a file.
        PYTHON_INCLUDE_DIR:PATH=/Library/Frameworks/Python.framework/Versions/2.7/Headers

        //Path to a library.
        PYTHON_LIBRARY:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib

12. If you are using Mavericks (OSX 10.9) and encounter errors relating to <tr1/memory>, you might also need to include the following flag in your CMakeCache.txt::

		CMAKE_CXX_FLAGS:STRING=-stdlib=libstdc++

13. Run make and install as follows::

        make -j2
        sudo make install

14. With any luck, you should have openbabel with python bindings installed.
    You can test your installation by trying to import openbabel from the
    python command line. Please note that despite best efforts,
    openbabel seems to install the python bindings into /usr/local/lib even
    if your Python is not the standard Mac version. In that case,
    you may need to add the following into your .bash_profile::

        export PYTHONPATH=/usr/local/lib:$PYTHONPATH

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