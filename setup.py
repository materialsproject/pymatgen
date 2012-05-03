import os
from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs


long_description = """
Pymatgen (python materials genomics) is the python library that powers the 
Materials Project (http://www.materialsproject.org). These are some of the main 
features:

1. Highly flexible classes for the representation of Element, Site, Molecule, 
   Structure objects.
2. Extensive io capabilities to manipulate many VASP input and output files 
   (http://cms.mpi.univie.ac.at/vasp/) and the crystallographic information file 
   format. This includes generating Structure objects from vasp input and
   output. There is also support for Gaussian input files and XYZ file for
   molecules.
3. Comprehensive tool to generate and view compositional and grand canonical phase 
   diagrams.
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
"""

spgsrcdir = os.path.join('extensions', 'spglib-1.1.2', 'src')

include_dirs = [spgsrcdir]
sources = ['cell.c', 'debug.c', 'hall_symbol.c', 'kpoint.c', 'lattice.c',
    'mathfunc.c', 'pointgroup.c', 'primitive.c', 'refinement.c',
    'sitesym_database.c', 'site_symmetry.c', 'spacegroup.c', 'spin.c',
    'spg_database.c', 'spglib.c', 'symmetry.c'
]

sources = [os.path.join(spgsrcdir, srcfile) for srcfile in sources]

extension = Extension('pymatgen._spglib',
                      include_dirs=include_dirs + get_numpy_include_dirs(),
                      sources=['_spglib.c'] + sources,
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-lgomp'],
                      )

setup (
  name='pymatgen',
  version='1.9.0',
  packages=['pymatgen', 'pymatgen.alchemy', 'pymatgen.analysis',
              'pymatgen.command_line', 'pymatgen.core', 'pymatgen.io',
              'pymatgen.phasediagram', 'pymatgen.symmetry',
              'pymatgen.transformations', 'pymatgen.util', 'pymatgen.vis',
              'pymatgen.alchemy.tests', 'pymatgen.analysis.tests',
              'pymatgen.command_line.tests', 'pymatgen.core.tests',
              'pymatgen.io.tests', 'pymatgen.phasediagram.tests',
              'pymatgen.symmetry.tests', 'pymatgen.transformations.tests',
              'pymatgen.util.tests'],
  install_requires=['numpy>=1.6', 'scipy>=0.10'],
  extras_require={
        'plotting':  ['matplotlib>=1.1'],
        'ase_adaptor': ['ase>=3.3'],
        'cif':['PyCifRW>=3.3']
  },
  package_data={'pymatgen.core': ['*.json'],
                  'pymatgen.io': ['*.cfg'],
                  'pymatgen.vis': ['ElementColorSchemes.cfg']},
  author='Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier, Will Richards, Dan Gunter, Vincent L Chevrier, Rickard Armiento',
  author_email='shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, geoffroy.hautier@uclouvain.be, wrichard@mit.edu, dkgunter@lbl.gov, vincentchevrier@gmail.com, armiento@mit.edu',
  maintainer='Shyue Ping Ong',
  url='https://github.com/materialsproject/pymatgen/',
  license='MIT',
  description="pymatgen is the Python library powering the Materials Project (www.materialsproject.org).",
  long_description=long_description,
  keywords=["vasp", "materials", "project", "electronic", "structure"],
  classifiers=[
        "Programming Language :: Python :: 2.7",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
  ],
  download_url="https://github.com/materialsproject/pymatgen/tarball/master",
  ext_modules=[extension]
)
