import os
from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

spgsrcdir = os.path.join('extensions', 'spglib-1.2', 'src')

include_dirs = [spgsrcdir]
sources = ['cell.c', 'debug.c', 'hall_symbol.c', 'kpoint.c', 'lattice.c',
    'mathfunc.c', 'pointgroup.c', 'primitive.c', 'refinement.c',
    'sitesym_database.c', 'site_symmetry.c', 'spacegroup.c', 'spin.c',
    'spg_database.c', 'spglib.c', 'symmetry.c'
]

sources = [os.path.join(spgsrcdir, srcfile) for srcfile in sources]

extension = Extension('pyspglib._spglib',
                      include_dirs=include_dirs + get_numpy_include_dirs(),
                      sources=['_spglib.c'] + sources,
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-lgomp'],
                      )

setup (
       name='spglib',
       version='1.2',
       description='This is the spglib module.',
       author='Atsushi Togo',
       author_email='atz.togo@gmail.com',
       url='http://spglib.sourceforge.net/',
       packages=['pyspglib'],
       ext_modules=[extension]
)
