# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""Setup.py for pymatgen."""

import sys
import platform

from setuptools import setup, find_namespace_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    """Extension builder that checks for numpy before install."""

    def finalize_options(self):
        """Override finalize_options."""
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        import builtins

        if hasattr(builtins, "__NUMPY_SETUP__"):
            # pylint: disable=E1101
            del builtins.__NUMPY_SETUP__
        import importlib
        import numpy

        importlib.reload(numpy)
        self.include_dirs.append(numpy.get_include())


extra_link_args = []
if sys.platform.startswith("win") and platform.machine().endswith("64"):
    extra_link_args.append("-Wl,--allow-multiple-definition")


long_desc = """
Official docs: [http://pymatgen.org](http://pymatgen.org/)

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. These are some of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive input/output support, including support for
   [VASP](http://cms.mpi.univie.ac.at/vasp/), [ABINIT](http://www.abinit.org/),
   CIF, Gaussian, XYZ, and many other file formats.
3. Powerful analysis tools, including generation of phase diagrams, Pourbaix
   diagrams, diffusion analyses, reactions, etc.
4. Electronic structure analyses, such as density of states and band structure.
5. Integration with the Materials Project REST API.

Pymatgen is free to use. However, we also welcome your help to improve this
library by making your own contributions.  These contributions can be in the
form of additional tools or modules you develop, or feature requests and bug
reports. Please report any bugs and issues at pymatgen's [Github page]
(https://github.com/materialsproject/pymatgen). For help with any pymatgen
issues, please use the [Discourse page](https://discuss.matsci.org/c/pymatgen).

Why use pymatgen?
=================

There are many materials analysis codes out there, both commerical and free,
but pymatgen offer several advantages:

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers,
   and is the analysis code powering the [Materials Project](https://www.materialsproject.org).
   The analysis it produces survives rigorous scrutiny every single day. Bugs
   tend to be found and corrected quickly. Pymatgen also uses
   [CircleCI](https://circleci.com) and [Appveyor](https://www.appveyor.com/)
   for continuous integration on the Linux and Windows platforms,
   respectively, which ensures that every commit passes a comprehensive suite
   of unittests.
2. **It is well documented.** A fairly comprehensive documentation has been
   written to help you get to grips with it quickly.
3. **It is open.** You are free to use and contribute to pymatgen. It also means
   that pymatgen is continuously being improved. We will attribute any code you
   contribute to any publication you specify. Contributing to pymatgen means
   your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in pymatgen have been
   optimized by vectorizing in numpy/scipy. This means that coordinate
   manipulations are extremely fast and are in fact comparable to codes
   written in other languages. Pymatgen also comes with a complete system for
   handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in
   the well-established Materials Project. It is also actively being developed
   and maintained by the [Materials Virtual Lab](https://www.materialsvirtuallab.org),
   the ABINIT group and many other research groups.

With effect from version 2021.1.1, pymatgen only supports Python >3.7.
"""

setup(
    name="pymatgen",
    packages=find_namespace_packages(
        include=["pymatgen.*", "pymatgen.analysis.*", "pymatgen.io.*", "pymatgen.ext.*"],
        exclude=["pymatgen.*.tests", "pymatgen.*.*.tests", "pymatgen.*.*.*.tests"],
    ),
    version="2022.0.10",
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.20.1",
        "requests",
        "ruamel.yaml>=0.15.6",
        "monty>=3.0.2",
        "scipy>=1.5.0",
        "tabulate",
        "spglib>=1.9.9.44",
        "networkx>=2.2",
        "matplotlib>=1.5",
        "palettable>=3.1.1",
        "sympy",
        "pandas",
        "plotly>=4.5.0",
        "uncertainties>=3.1.4",
    ],
    extras_require={
        "provenance": ["pybtex"],
        "ase": ["ase>=3.3"],
        "vis": ["vtk>=6.0.0"],
        "abinit": ["netcdf4"],
        ':python_version < "3.8"': [
            "typing-extensions>=3.7.4.3",
        ],
    },
    # All package data has to be explicitly defined. Do not use automated codes like last time. It adds
    # all sorts of useless files like test files and is prone to path errors.
    package_data={
        "pymatgen.analysis": ["*.yaml", "*.json", "*.csv"],
        "pymatgen.analysis.chemenv": [
            "coordination_environments/coordination_geometries_files/*.json",
            "coordination_environments/coordination_geometries_files/*.txt",
            "coordination_environments/strategy_files/ImprovedConfidenceCutoffDefaultParameters.json",
        ],
        "pymatgen.analysis.structure_prediction": ["*.yaml", "data/*.json"],
        "pymatgen.analysis.diffraction": ["*.json"],
        "pymatgen.analysis.magnetism": ["default_magmoms.yaml"],
        "pymatgen.entries": ["py.typed", "*.json.gz", "*.yaml", "data/*.json"],
        "pymatgen.core": ["py.typed", "*.json"],
        "pymatgen.io.vasp": ["*.yaml", "*.json"],
        "pymatgen.io.feff": ["*.yaml"],
        "pymatgen.io.cp2k": ["*.yaml"],
        "pymatgen.io.lobster": ["lobster_basis/*.yaml"],
        "pymatgen.command_line": ["OxideTersoffPotentials"],
        "pymatgen.util": ["structures/*.json", "*.json"],
        "pymatgen.vis": ["*.yaml"],
        "pymatgen.io.lammps": ["CoeffsDataType.yaml"],
        "pymatgen.symmetry": ["*.yaml", "*.json", "*.sqlite"],
    },
    author="Pymatgen Development Team",
    author_email="ongsp@eng.ucsd.edu",
    maintainer="Shyue Ping Ong, Matthew Horton",
    maintainer_email="ongsp@eng.ucsd.edu, mkhorton@lbl.gov",
    url="http://www.pymatgen.org",
    license="MIT",
    description="Python Materials Genomics is a robust materials "
    "analysis code that defines core object representations for "
    "structures and molecules with support for many electronic "
    "structure codes. It is currently the core analysis code "
    "powering the Materials Project "
    "(https://www.materialsproject.org).",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    keywords=[
        "VASP",
        "gaussian",
        "ABINIT",
        "nwchem",
        "qchem",
        "materials",
        "science",
        "project",
        "electronic",
        "structure",
        "analysis",
        "phase",
        "diagrams",
        "crystal",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    ext_modules=[
        Extension(
            "pymatgen.optimization.linear_assignment",
            ["pymatgen/optimization/linear_assignment.c"],
            extra_link_args=extra_link_args,
        ),
        Extension("pymatgen.util.coord_cython", ["pymatgen/util/coord_cython.c"], extra_link_args=extra_link_args),
        Extension(
            "pymatgen.optimization.neighbors", ["pymatgen/optimization/neighbors.c"], extra_link_args=extra_link_args
        ),
    ],
    entry_points={
        "console_scripts": [
            "pmg = pymatgen.cli.pmg:main",
            "feff_plot_cross_section = pymatgen.cli.feff_plot_cross_section:main",
            "feff_plot_dos = pymatgen.cli.feff_plot_dos:main",
            "gaussian_analyzer = pymatgen.cli.gaussian_analyzer:main",
            "get_environment = pymatgen.cli.get_environment:main",
        ]
    },
)
