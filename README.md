## Introduction ##

Pymatgen (**py**thon **mat**erials **gen**omics) is the python library that powers 
the Materials Project (http://www.materialsproject.org). These are some of the 
key features:

1. Highly flexible classes for the representation of Element, Site, Structure 
   objects.
2. Powerful io capabilities to manipulate many VASP input and output files 
   (http://cms.mpi.univie.ac.at/vasp/) and the crystallographic information 
   file format.  This includes generating Structure objects from vasp input and 
   output.
3. A comprehensive tool to generate and view compositional and grand canonical 
   phase diagrams.
4. Electronic structure analyses (DOS and Bandstructure).

The public version of pymatgen is free (as in free beer) to download and to use. 
However, we would also like you to help us improve this library by making your 
own contributions as well.  These contributions can be in the form of additional 
tools or modules you develop, or even simple things such as bug reports.  Please 
contact the maintainer of this library (shyue@mit.edu) to find out how to include 
your contributions via github or for bug reports.

Note that pymatgen, like all scientific research, will always be a work in
progress. While the development team will always strive to avoid backward 
incompatible changes, sometimes those are unavoidable and tough decisions have 
to be made for the long term health of the code. 

For documentation, usage examples and change log, please read the documentation 
at http://materialsproject.github.com/pymatgen.

## Requirements ##

Required for proper functioning of the code.

1. Python 2.7+ required.  New default modules such as json are used, as well as 
   new unittest features in Python 2.7.
2. numpy - For array, matrix and other numerical manipulations. Used extensively 
   by all core modules.
3. scipy 0.9+ - For interpolation, physical constants and other functions. In 
   particular, scipy.spatial.Delaunay is used for phase diagram construction.
5. nose - For complete unittesting. This is NOT optional!

## Optional Python Libraries ##

Optional python libraries that are required if you need certain features.

1. matplotlib : For plotting (e.g., Phase Diagrams) using the pymatgen.phasediagrams 
   package.
2. [PyCifRW](http://prdownload.berlios.de/pycifrw/PyCifRW-3.3.tar.gz) : For 
   reading and writing Crystallographic Information Format (CIF) files using 
   the pymatgen.io.cifio module [more info](http://pycifrw.berlios.de/)
3. [Pyspglib](http://spglib.sourceforge.net/) : For symmetry finding using the 
   pymatgen.symmetry package.
4. VTK with Python bindings (http://www.vtk.org/): For visualization of crystal 
   structures using the pymatgen.vis package.
5. Atomistic Simulation Environment or ASE (https://wiki.fysik.dtu.dk/ase/): 
   Required for the usage of the adapters in pymatgen.io.aseio between pymatgen's 
   core Structure object and the Atoms object used by ASE. 

## Optional non-Python programs ##

Optional non-python libraries (because no good pythonic alternative exists at 
the moment) required only for certain features.

1. [Qhull](http://www.qhull.org/) : Needed for bond length analysis 
   (structure_analyzer.py).  The executable qconvex and qvoronoi must be in the path.
2. [ffmpeg](http://www.http://ffmpeg.org//) : Needed for generation of movies 
   (structure_vtk.py).  The executable ffmpeg must be in the path.

## Basic Setup ##

1. Clone the repo.
2. Install the necessary python libraries.
3. (Recommended) Add pymatgen to your PYTHONPATH.
4. (Recommended for developers) Copy hooks from the example-hooks directory into 
   the .git/hooks/ directory in your local repo.  

With these two basic steps, you should be able to use most of the pymatgen code.  
I recommend that you start by reading some of the unittests in the tests 
subdirectory for each package.  The unittests demonstrate the expected behavior 
and functionality of the code.

However, some extra functionality do require additional setup, as outlined below.

### Generating POTCARs ###

For the code to generate POTCAR files, it needs to know where the VASP 
pseudopotential files are.  We are not allowed to distribute these under the 
VASP license. The good news is that we have included a setup script to help you 
along.

If you cloned the repo directly from github, you should have a run_me_first.sh 
file in the root directory of your local repo. Otherwise, you can get it directly 
from our github site at http://github.com/materialsproject/pymatgen. Run the 
shell script and follow the instructions. If you have done it correctly, you 
should get a resources directory with the following directory structure:

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
