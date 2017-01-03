# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
from io import open
import sys
import platform

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        if sys.version_info[0] >= 3:
            import builtins
            if hasattr(builtins, '__NUMPY_SETUP__'):
                del builtins.__NUMPY_SETUP__
            import importlib
            import numpy
            importlib.reload(numpy)
        else:
            import __builtin__
            if hasattr(__builtin__, '__NUMPY_SETUP__'):
                del __builtin__.__NUMPY_SETUP__
            import imp
            import numpy
            imp.reload(numpy)

        self.include_dirs.append(numpy.get_include())

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
    version="4.5.7",
    cmdclass={'build_ext': build_ext},
    setup_requires=['numpy', 'setuptools>=18.0'],
    install_requires=["numpy>=1.9", "six", "requests", "pyyaml>=3.11",
                      "monty>=0.9.6", "scipy>=0.14", "pydispatcher>=2.0.5",
                      "tabulate", "spglib>=1.9.8.7",
                      "matplotlib>=1.5", "palettable>=2.1.1"],
    extras_require={
        ':python_version == "2.7"': [
            'enum34',
            'pathlib2',
        ],
        "matproj.snl": ["pybtex"],
        "pourbaix diagrams, bandstructure": ["pyhull>=1.5.3"],
        "ase_adaptor": ["ase>=3.3"],
        "vis": ["vtk>=6.0.0"],
        "abinit": ["pydispatcher>=2.0.5", "apscheduler==2.1.0"]},
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
        "Programming Language :: Python :: 3.6",
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
                           extra_link_args=extra_link_args),
                 Extension("pymatgen.util.coord_utils_cython",
                           ["pymatgen/util/coord_utils_cython.c"],
                           extra_link_args=extra_link_args)],
    entry_points={
          'console_scripts': [
              'pmg = pymatgen.cli.pmg:main',
              'feff_input_generation = pymatgen.cli.feff_input_generation:main',
              'feff_plot_cross_section = pymatgen.cli.feff_plot_cross_section:main',
              'feff_plot_dos = pymatgen.cli.feff_plot_dos:main',
              'gaussian_analyzer = pymatgen.cli.gaussian_analyzer:main',
              'get_environment = pymatgen.cli.get_environment:main',
              'pydii = pymatgen.cli.pydii:main',
          ]
    }
)
