.. pymatgen documentation master file, created by
   sphinx-quickstart on Tue Nov 15 00:13:52 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation instructions
====================================

Requirements
============

Required for proper functioning of the code.
--------------------------------------------

1. Python 2.7+ required.  New default modules such as json are used, as well as new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively by all core modules.
3. scipy 0.9+ - For interpolation, physical constants and other functions. In particular, scipy.spatial.Delaunay is used for phase diagram construction.
4. PyYAML - For parsing of important PyYaml configuration files.
5. nose - For complete unittesting. This is NOT optional!

Optional Python Libraries
-------------------------
Optional python libraries that are required if you need certain features

1. matplotlib : For plotting (e.g., Phase Diagrams).
2. PyCifRW (http://prdownload.berlios.de/pycifrw/PyCifRW-3.3.tar.gz) : For reading and writing Crystallographic Information Format (CIF) files [more info](http://pycifrw.berlios.de/)

Optional non-Python programs
----------------------------

Optional non-python libraries (because no good pythonic alternative exists at the moment) required only for certain features.

1. Qhull (http://www.qhull.org/) : Needed for bond length analysis (structure_analyzer.py).  The executable qconvex and qvoronoi must be in the path.

Basic Setup
===========

1. Clone the repo.
2. Install the necessary python libraries.
3. (Recommended) Add pymatgen to your PYTHONPATH.
4. (Recommended for developers) Copy hooks from the example-hooks directory into the .git/hooks/ directory in your local repo.  

With these two basic steps, you should be able to use most of the pymatgen code.  I recommend that you start by reading some of the unittests in the tests subdirectory for each package.  The unittests demonstrate the expected behavior and functionality of the code.

However, some extra functionality do require additional setup, as outlined below.

Generating POTCARs
------------------

For the code to generate POTCAR files, it needs to know where the VASP pseudopotential files are.  We are not allowed to distribute these under the VASP license. The good news is that we have included a setup script to help you along.

1. cd to the root directory of the repo where a file called run_me_first.sh is present.
2. Run the run_me_first.sh file, which will generate a resources directory in a location of your choosing. Please choose a location *outside* of the repo itself.  The script will also write a pymatgen.cfg file in the pymatgen subdir.

Basic usage
===========

Some example scripts have been provided in the scripts directory. In general, most file format conversions, manipulations and io can be done with a few quick lines of code. For example, to read a POSCAR and write a cif:

::

	from pymatgen.io.vaspio import Poscar
	from pymatgen.io.cifio import CifWriter

	p = Poscar('POSCAR')
	w = CifWriter(p.struct)
	w.write_file('mystructure.cif')

For more examples, please take a look in the scripts directory. More examples will be added soon.