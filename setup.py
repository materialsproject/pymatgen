# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import glob
import os
from io import open
import sys

from setuptools import setup, find_packages, Extension

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Please install numpy"
          "first before install pymatgen...")
    sys.exit(-1)

SETUP_PTH = os.path.dirname(os.path.abspath(__file__))


def get_spglib_ext():
    """
    Set up spglib extension.
    """
    spglibs = glob.glob(os.path.join(SETUP_PTH, "dependencies", "spglib*"))
    if len(spglibs) != 1:
        raise ValueError("Incorrect number of spglib found in dependencies. "
                         "Expected 1, got %d" % len(spglibs))
    spglibdir = spglibs[0]

    # set rest of spglib
    spgsrcdir = os.path.join(spglibdir, "src")
    include_dirs = [spgsrcdir]
    sources = glob.glob(os.path.join(spgsrcdir, "*.c"))
    c_opt = [] if sys.version_info.major < 3 else [
        "-Wno-error=declaration-after-statement"]
    return Extension(
        "pymatgen._spglib",
        include_dirs=include_dirs + get_numpy_include_dirs(),
        sources=[os.path.join(spglibdir, "_spglib.c")] + sources,
        extra_compile_args=c_opt)


with open("README.rst") as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind + 1:]


setup(
    name="pymatgen",
    packages=find_packages(),
    version="3.2.10",
    install_requires=["numpy>=1.9", "six", "atomicfile", "requests",
                      "pybtex", "pyyaml", "monty>=0.7.0", "scipy>=0.14",
                      "tabulate"],
    extras_require={"plotting": ["matplotlib>=1.1", "prettyplotlib"],
                    "pourbaix diagrams, bandstructure": ["pyhull>=1.5.3"],
                    "ase_adaptor": ["ase>=3.3"],
                    "vis": ["vtk>=6.0.0"],
                    "abinit": ["pydispatcher>=2.0.3", "apscheduler==2.1.0"]},
    package_data={"pymatgen.core": ["*.json"],
                  "pymatgen.analysis": ["*.yaml", "*.csv"],
                  "pymatgen.io.vasp": ["*.yaml"],
                  "pymatgen.io.feff": ["*.yaml"],
                  "pymatgen.symmetry": ["*.yaml"],
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
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
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
    ext_modules=[get_spglib_ext()],
    scripts=glob.glob(os.path.join(SETUP_PTH, "scripts", "*"))
)
