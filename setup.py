import glob
import os
import subprocess

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages, Extension

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Attempting to "
          "install...")
    subprocess.call(["easy_install", "numpy"])
    from numpy.distutils.misc_util import get_numpy_include_dirs


def get_spglib_ext():
    """
    Set up spglib extension.
    """
    spglibs = glob.glob(os.path.join("dependencies", "spglib*"))
    if len(spglibs) == 0:
        raise ValueError("No spglib found in dependencies.")
    spglibdir = spglibs[0]

    # set rest of spglib
    spgsrcdir = os.path.join(spglibdir, "src")
    include_dirs = [spgsrcdir]
    sources = [os.path.join(spgsrcdir, srcfile) for srcfile in
        os.listdir(spgsrcdir) if srcfile.endswith(".c")]
    return Extension("pymatgen._spglib",
                     include_dirs=include_dirs + get_numpy_include_dirs(),
                     sources=[os.path.join(spglibdir, "_spglib.c")] + sources)

with open("README.rst") as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind + 1:]

setup(
    name="pymatgen",
    packages=find_packages(),
    version="2.9.11",
    install_requires=["numpy>=1.5", "pyhull>=1.4.3", "PyCifRW>=3.3",
                      "requests>=1.0", "pybtex>=0.16", "pyyaml>=3.0",
                      "monty>=0.2.2"],
    extras_require={"electronic_structure": ["scipy>=0.10"],
                    "plotting": ["matplotlib>=1.1"],
                    "ase_adaptor": ["ase>=3.3"],
                    "vis": ["vtk>=6.0.0"],
                    "abinitio": ["pydispatcher>=2.0", "apscheduler>=2.1.1"]},
    package_data={"pymatgen.core": ["*.json"],
                  "pymatgen.analysis": ["bvparam_1991.json", "icsd_bv.json"],
                  "pymatgen.io": ["*.cfg", "*.json"],
                  "pymatgen.entries": ["*.cfg"],
                  "pymatgen.structure_prediction": ["data/*.json"],
                  "pymatgen.vis": ["ElementColorSchemes.cfg"],
                  "pymatgen.command_line": ["OxideTersoffPotentials"],
                  "pymatgen.analysis.defects": ["*.json"],
                  "pymatgen.analysis.diffraction": ["*.json"]},
    author="Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier,"
    "William Davidson Richards, Stephen Dacek, Dan Gunter, Shreyas Cholia, "
    "Matteo Giantomassi, Vincent L Chevrier, Rickard Armiento",
    author_email="ongsp@ucsd.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, "
    "geoffroy.hautier@uclouvain.be, wrichard@mit.edu, sdacek@mit.edu, "
    "dkgunter@lbl.gov, scholia@lbl.gov, gmatteo@gmail.com, "
    "vincentchevrier@gmail.com, armiento@mit.edu",
    maintainer="Shyue Ping Ong",
    url="https://github.com/materialsproject/pymatgen/",
    license="MIT",
    description="Python Materials Genomics is a robust materials "
                "analysis code that defines core object representations for "
                "structures and molecules with support for many electronic "
                "structure codes. It is currently the core analysis code "
                "powering the Materials Project (www.materialsproject.org).",
    long_description=long_desc,
    keywords=["VASP", "gaussian", "ABINIT", "nwchem", "materials", "project",
              "electronic", "structure", "analysis", "phase", "diagrams"],
    classifiers=[
        "Programming Language :: Python :: 2.7",
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
    scripts=[os.path.join("scripts", f) for f in os.listdir("scripts")]
)
