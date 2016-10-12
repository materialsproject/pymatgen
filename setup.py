# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import glob
import os
from io import open
import sys
import platform

from setuptools import setup, find_packages, Extension

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Please install numpy"
          "first before install pymatgen...")
    sys.exit(-1)

SETUP_PTH = os.path.dirname(__file__)

extra_link_args = []
if sys.platform.startswith('win') and platform.machine().endswith('64'):
    extra_link_args.append('-Wl,--allow-multiple-definition')

with open(os.path.join(SETUP_PTH, "README.rst")) as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind + 1:]

setup(
    name="pymatgen",
    packages=find_packages(),
    version="4.4.1",
    install_requires=["numpy>=1.9", "six", "atomicfile", "requests",
                      "pybtex", "pyyaml", "monty>=0.9.5", "scipy>=0.14",
                      "tabulate", "enum34", "spglib"],
    extras_require={"plotting": ["matplotlib>=1.1", "prettyplotlib"],
                    "pourbaix diagrams, bandstructure": ["pyhull>=1.5.3"],
                    "ase_adaptor": ["ase>=3.3"],
                    "vis": ["vtk>=6.0.0"],
                    "abinit": ["pydispatcher>=2.0.3", "apscheduler==2.1.0"],
                    "chemenv": ["unittest2"]},
    package_data={"pymatgen.core": ["*.json"],
                  "pymatgen.analysis": ["*.yaml", "*.csv"],
                  "pymatgen.analysis.chemenv.coordination_environments.coordination_geometries_files": ["*.txt", "*.json"],
                  "pymatgen.analysis.chemenv.coordination_environments.strategy_files": ["*.json"],
                  "pymatgen.io.vasp": ["*.yaml"],
                  "pymatgen.io.feff": ["*.yaml"],
                  "pymatgen.symmetry": ["*.yaml", "*.json"],
                  "pymatgen.entries": ["*.yaml"],
                  "pymatgen.structure_prediction": ["data/*.json"],
                  "pymatgen.vis": ["ElementColorSchemes.yaml"],
                  "pymatgen.command_line": ["OxideTersoffPotentials"],
                  "pymatgen.analysis.defects": ["*.json"],
                  "pymatgen.analysis.diffraction": ["*.json"],
                  "pymatgen.util": ["structures/*.json"]},
    author="Pymatgen Development Team",
    author_email="pymatgen@googlegroups.com",
    maintainer="Shyue Ping Ong",
    maintainer_email="ongsp@eng.ucsd.edu",
    url="http://www.pymatgen.org",
    license="MIT",
    description="Python Materials Genomics is a robust materials "
                "analysis code that defines core object representations for "
                "structures and molecules with support for many electronic "
                "structure codes. It is currently the core analysis code "
                "powering the Materials Project "
                "(https://www.materialsproject.org).",
    long_description=long_desc,
    keywords=["VASP", "gaussian", "ABINIT", "nwchem", "materials", "project",
              "electronic", "structure", "analysis", "phase", "diagrams"],
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    ext_modules=[Extension("pymatgen.optimization.linear_assignment",
                           ["pymatgen/optimization/linear_assignment.c"],
                           include_dirs=get_numpy_include_dirs(),
                           extra_link_args=extra_link_args),
                 Extension("pymatgen.util.coord_utils_cython",
                           ["pymatgen/util/coord_utils_cython.c"],
                           include_dirs=get_numpy_include_dirs(),
                           extra_link_args=extra_link_args)],
    entry_points={
          'console_scripts': [
              'pmg = pymatgen.cli.pmg:main',
              'feff_input_generation = pymatgen.cli.feff_input_generation:main',
              'feff_plot_cross_section = pymatgen.cli.feff_plot_cross_section:main',
              'feff_plot_dos = pymatgen.cli.feff_plot_dos:main',
              'gaussian_analyzer = pymatgen.cli.gaussian_analyzer:main',
              'gen_potcar = pymatgen.cli.gen_potcar:main',
              'get_environment = pymatgen.cli.get_environment:main',
              'pydii = pymatgen.cli.pydii:main',
          ]
    },
    scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*"))
)
