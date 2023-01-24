# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import platform
import sys

import numpy
from setuptools import Extension, find_namespace_packages, setup

extra_link_args: list[str] = []
if sys.platform.startswith("win") and platform.machine().endswith("64"):
    extra_link_args = ["-Wl,--allow-multiple-definition"]


setup(
    name="pymatgen",
    packages=find_namespace_packages(
        include=["pymatgen.*", "pymatgen.analysis.*", "pymatgen.io.*", "pymatgen.ext.*", "cmd_line"],
        exclude=["pymatgen.*.tests", "pymatgen.*.*.tests", "pymatgen.*.*.*.tests"],
    ),
    version="2023.1.20",
    python_requires=">=3.8",
    install_requires=[
        "matplotlib>=1.5",
        "monty>=3.0.2",
        "mp-api>=0.27.3",
        "networkx>=2.2",
        "numpy>=1.20.1",
        "palettable>=3.1.1",
        "pandas",
        "plotly>=4.5.0",
        "pybtex",
        "requests",
        "ruamel.yaml>=0.17.0",
        "scipy>=1.5.0",
        "spglib>=2.0.2",
        "sympy",
        "tabulate",
        "tqdm",
        "uncertainties>=3.1.4",
    ],
    extras_require={
        "ase": ["ase>=3.3"],
        "vis": ["vtk>=6.0.0"],
        "abinit": ["netcdf4"],
        "relaxation": ["m3gnet"],
        "dev": [
            "black",
            "coverage",
            "coveralls",
            "flake8",
            "mypy==0.991",  # pinned due to long list of errors starting with mypy 0.990
            "pre-commit",
            "pydocstyle",
            "pylint",
            "pytest",
            "pytest-cov",
            "pytest-split",
        ],
        "docs": [
            "sphinx",
            "sphinx_rtd_theme",
            "doc2dash",
        ],
        "optional": [
            # "hiphive>=0.6",
            # "m3gnet>=0.0.8",
            "ase>=3.22.1",
            # https://peps.python.org/pep-0508/#environment-markers
            "BoltzTraP2>=22.3.2; platform_system!='Windows'",
            "chemview>=0.6",
            "f90nml>=1.1.2",
            "fdint>=2.0.2",
            "galore>=0.6.1",
            "h5py==3.6.0",  # pinned due to 3.7 crashing on windows
            # https://github.com/h5py/h5py/issues/2110
            "jarvis-tools>=2020.7.14",
            "netCDF4>=1.5.8",
            "phonopy>=2.4.2",
            "seekpath>=1.9.4",
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
        "pymatgen.analysis.solar": ["am1.5G.dat"],
        "pymatgen.entries": ["py.typed", "*.json.gz", "*.yaml", "data/*.json"],
        "pymatgen.core": ["py.typed", "*.json"],
        "pymatgen.io.vasp": ["*.yaml", "*.json"],
        "pymatgen.io.feff": ["*.yaml"],
        "pymatgen.io.cp2k": ["*.yaml"],
        "pymatgen.io.lobster": ["lobster_basis/*.yaml"],
        "pymatgen.command_line": ["*"],
        "pymatgen.util": ["structures/*.json", "*.json"],
        "pymatgen.vis": ["*.yaml"],
        "pymatgen.io.lammps": ["CoeffsDataType.yaml", "templates/md.txt"],
        "pymatgen.symmetry": ["*.yaml", "*.json", "*.sqlite"],
        "cmd_line": ["**/*"],
    },
    author="Pymatgen Development Team",
    author_email="ongsp@eng.ucsd.edu",
    maintainer="Shyue Ping Ong, Matthew Horton, Janosh Riebesell",
    maintainer_email="ongsp@eng.ucsd.edu, mkhorton@lbl.gov, janosh.riebesell@gmail.com",
    url="https://pymatgen.org",
    license="MIT",
    project_urls={
        "Docs": "https://pymatgen.org",
        "Package": "https://pypi.org/project/pymatgen",
        "Repo": "https://github.com/materialsproject/pymatgen",
    },
    description="Python Materials Genomics is a robust materials "
    "analysis code that defines core object representations for "
    "structures and molecules with support for many electronic "
    "structure codes. It is currently the core analysis code "
    "powering the Materials Project "
    "(https://materialsproject.org).",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    keywords=[
        "ABINIT",
        "analysis",
        "crystal",
        "diagrams",
        "electronic",
        "gaussian",
        "materials",
        "nwchem",
        "phase",
        "project",
        "qchem",
        "science",
        "structure",
        "VASP",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
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
    include_dirs=[numpy.get_include()],
)
