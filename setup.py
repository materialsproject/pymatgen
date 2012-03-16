import os
from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

README = os.path.join(os.path.dirname(__file__), 'README.md')
long_description = open(README).read() + '\n\n'

spgsrcdir = os.path.join('extensions', 'spglib-1.1.2', 'src')

include_dirs = [spgsrcdir]
sources = ['cell.c', 'debug.c', 'hall_symbol.c', 'kpoint.c', 'lattice.c',
    'mathfunc.c', 'pointgroup.c', 'primitive.c', 'refinement.c',
    'sitesym_database.c', 'site_symmetry.c', 'spacegroup.c', 'spin.c',
    'spg_database.c', 'spglib.c', 'symmetry.c'
]

sources = [os.path.join(spgsrcdir, srcfile) for srcfile in sources]

extension = Extension('pymatgen._spglib',
                      include_dirs = include_dirs + get_numpy_include_dirs(),
                      sources = ['_spglib.c'] + sources,
                      extra_compile_args = ['-fopenmp'],
                      extra_link_args = ['-lgomp'],
                      )

setup (
  name = 'pymatgen',
  version = '1.5.0',
  packages = ['pymatgen', 'pymatgen.alchemy', 'pymatgen.analysis',
              'pymatgen.command_line', 'pymatgen.core', 'pymatgen.io',
              'pymatgen.phasediagram', 'pymatgen.symmetry',
              'pymatgen.transformations', 'pymatgen.util', 'pymatgen.vis',
              'pymatgen.alchemy.tests', 'pymatgen.analysis.tests',
              'pymatgen.command_line.tests', 'pymatgen.core.tests',
              'pymatgen.io.tests', 'pymatgen.phasediagram.tests',
              'pymatgen.symmetry.tests', 'pymatgen.transformations.tests',
              'pymatgen.util.tests'],
  install_requires = ['numpy', 'scipy', 'PyCIFRW'],
  package_data = {'pymatgen.core': ['*.json'],
                  'pymatgen.io': ['*.cfg'],
                  'pymatgen.vis': ['ElementColorSchemes.cfg']},
  author = 'Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier, Will Richards, Dan Gunter, Vincent L Chevrier, Rickard Armiento',
  author_email = 'shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, geoffroy.hautier@uclouvain.be, wrichard@mit.edu, dkgunter@lbl.gov, vincentchevrier@gmail.com, armiento@mit.edu',
  maintainer = 'Shyue Ping Ong',
  url = 'https://github.com/CederGroupMIT/pymatgen_repo/',
  license = 'MIT',
  description = "pymatgen is the Python library powering the Materials Project (www.materialsproject.org).",
  long_description = long_description,
  keywords = ["vasp", "materials", "project", "electronic", "structure"],
  classifiers = [
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
  download_url = "https://github.com/CederGroupMIT/pymatgen_repo/tarball/master",
  ext_modules = [extension],
  test_suite = 'nose.collector',
  test_requires = ['nose']
)
