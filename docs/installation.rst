Installation instructions
====================================

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

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good pythonic alternative exists at the moment) required only for certain features.

1. Qhull (http://www.qhull.org/) : Needed for bond length analysis (structure_analyzer.py).  The executable qconvex and qvoronoi must be in the path.

Basic Setup for Users
=====================

pymatgen is now on PyPI (http://pypi.python.org/pypi/pymatgen).  You can now just do 

::

	easy_install pymatgen
	
to install pymatgen with most of the dependencies set up. Alternatively, you can download the latest source and run 

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