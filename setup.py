# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""Setup.py for pymatgen."""

from __future__ import annotations

import platform
import sys

import numpy
from setuptools import Extension, find_namespace_packages, setup

extra_link_args: list[str] = []
if sys.platform.startswith("win") and platform.machine().endswith("64"):
    extra_link_args.append("-Wl,--allow-multiple-definition")


long_desc = """
Official docs: [https://pymatgen.org](https://pymatgen.org/)

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. These are some of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive input/output support, including support for
   [VASP](https://www.vasp.at), [ABINIT](https://www.abinit.org/),
   CIF, Gaussian, XYZ, and many other file formats.
3. Powerful analysis tools, including generation of phase diagrams, Pourbaix
   diagrams, diffusion analyses, reactions, etc.
4. Electronic structure analyses, such as density of states and band structure.
5. Integration with the Materials Project REST API.

Pymatgen is free to use. However, we also welcome your help to improve this
library by making your own contributions. These contributions can be in the
form of additional tools or modules you develop, or feature requests and bug
reports. Please report any bugs and issues at pymatgen's [Github page]
(https://github.com/materialsproject/pymatgen). For help with any pymatgen
issues, please use the [Discourse page](https://discuss.matsci.org/c/pymatgen).

Why use pymatgen?
=================

There are many materials analysis codes out there, both commercial and free,
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
"""

setup(
    name="pymatgen",
    packages=find_namespace_packages(
        include=["pymatgen.*", "pymatgen.analysis.*", "pymatgen.io.*", "pymatgen.ext.*"],
        exclude=["pymatgen.*.tests", "pymatgen.*.*.tests", "pymatgen.*.*.*.tests"],
    ),
    version="2022.7.25",
    python_requires=">=3.8",
    setup_requires=[
        "Cython>=0.29.23",
    ],
    install_requires=[
        "numpy>=1.20.1",
        "requests",
        "ruamel.yaml>=0.17.0",
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
        "pybtex",
        "tqdm",
    ],
    extras_require={"ase": ["ase>=3.3"], "vis": ["vtk>=6.0.0"], "abinit": ["netcdf4"], "relaxation": ["m3gnet"]},
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
        "pymatgen.analysis.solar": ["am1.5G.dat"],
        "pymatgen.entries": ["py.typed", "*.json.gz", "*.yaml", "data/*.json"],
        "pymatgen.core": ["py.typed", "*.json"],
        "pymatgen.io.vasp": ["*.yaml", "*.json"],
        "pymatgen.io.feff": ["*.yaml"],
        "pymatgen.io.cp2k": ["*.yaml"],
        "pymatgen.io.lobster": ["lobster_basis/*.yaml"],
        "pymatgen.command_line": ["OxideTersoffPotentials"],
        "pymatgen.util": ["structures/*.json", "*.json"],
        "pymatgen.vis": ["*.yaml"],
        "pymatgen.io.lammps": ["CoeffsDataType.yaml", "templates/md.txt"],
        "pymatgen.symmetry": ["*.yaml", "*.json", "*.sqlite"],
    },
    author="Pymatgen Development Team",
    author_email="ongsp@eng.ucsd.edu",
    maintainer="Shyue Ping Ong, Matthew Horton",
    maintainer_email="ongsp@eng.ucsd.edu, mkhorton@lbl.gov",
    url="https://pymatgen.org",
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
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
            ["pymatgen/optimization/linear_assignment.pyx"],
            extra_link_args=extra_link_args,
        ),
        Extension("pymatgen.util.coord_cython", ["pymatgen/util/coord_cython.pyx"], extra_link_args=extra_link_args),
        Extension(
            "pymatgen.optimization.neighbors", ["pymatgen/optimization/neighbors.pyx"], extra_link_args=extra_link_args
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
    include_dirs=numpy.get_include(),
)
