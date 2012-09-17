import os
import sys

from distribute_setup import use_setuptools
use_setuptools(version='0.6.10')
from setuptools import setup, find_packages, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
from pymatgen import __version__

long_description = """
Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. It currently powers the public Materials Project
(http://www.materialsproject.org), an initiative to make calculated properties
on a large number of materials available to materials researchers and
designers. These are some of the main features:
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
The pymatgen library is free (as in free beer) to download and to use.
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
pymatgen"s Google Groups page
(https://groups.google.com/forum/?fromgroups#!forum/pymatgen/).
"""

spgsrcdir = os.path.join("dependencies", "spglib-1.2.2", "src")
include_dirs = [spgsrcdir]
sources = ["cell.c", "debug.c", "hall_symbol.c", "kpoint.c", "lattice.c",
           "mathfunc.c", "pointgroup.c", "primitive.c", "refinement.c",
           "sitesym_database.c", "site_symmetry.c", "spacegroup.c", "spin.c",
           "spg_database.c", "spglib.c", "symmetry.c"]
sources = [os.path.join(spgsrcdir, srcfile) for srcfile in sources]

if os.name == "posix" and sys.platform in ("darwin", "sunos5"):
    print "Mac OS detected. Compiling without openmp..."
    extra_compile = []
    extra_link = []
else:
    extra_compile = ["-fopenmp"]
    extra_link = ["-lgomp"]

extension = Extension("pymatgen._spglib",
                      include_dirs=include_dirs + get_numpy_include_dirs(),
                      sources=[os.path.join("dependencies", "spglib-1.2.2",
                               "_spglib.c")] + sources,
                      extra_compile_args=extra_compile,
                      extra_link_args=extra_link
                      )

scripts = [os.path.join("scripts", f) for f in os.listdir("scripts")]

setup(name="pymatgen",
      packages=find_packages(),
      version=__version__,
      install_requires=["numpy>=1.5"],
      extras_require={"phasediagrams": ["scipy>=0.10"],
                      "plotting": ["matplotlib>=1.1"],
                      "ase_adaptor": ["ase>=3.3"],
                      "cif": ["PyCifRW>=3.3"]},
      package_data={"pymatgen.core": ["bond_lengths.json",
                                      "periodic_table.json"],
                    "pymatgen.io": ["*.cfg"],
                    "pymatgen.entries": ["*.cfg"],
                    "pymatgen.vis": ["ElementColorSchemes.cfg"]},
      author="Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier,"
      "Will Richards, Dan Gunter, Shreyas Cholia, Vincent L Chevrier, "
      "Rickard Armiento",
      author_email="shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, "
      "geoffroy.hautier@uclouvain.be, wrichard@mit.edu, dkgunter@lbl.gov, "
      "vincentchevrier@gmail.com, armiento@mit.edu",
      maintainer="Shyue Ping Ong",
      url="https://github.com/materialsproject/pymatgen/",
      license="MIT",
      description="pymatgen is the Python library powering the "
      "Materials Project (www.materialsproject.org).",
      long_description=long_description,
      keywords=["vasp", "materials", "project", "electronic", "structure"],
      classifiers=["Programming Language :: Python :: 2.7",
                   "Development Status :: 4 - Beta",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Information Analysis",
                   "Topic :: Scientific/Engineering :: Physics",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   "Topic :: Software Development :: Libraries :: Python Modules"],
      download_url="https://github.com/materialsproject/pymatgen/tarball/master",
      ext_modules=[extension],
      scripts=scripts
      )
