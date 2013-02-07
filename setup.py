import glob
import os
import sys

from distribute_setup import use_setuptools
use_setuptools(version='0.6.10')
from setuptools import setup, find_packages, Extension

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
except ImportError:
    print "numpy.distutils.misc_util cannot be imported."
    print "numpy.distutils.misc_util is needed to build the spglib extension."
    print "Please install numpy first before retrying setup."
    sys.exit(-1)

def get_spglib_ext():
    # Get spglib
    spglibs = glob.glob(os.path.join("dependencies", "spglib*"))
    if len(spglibs) == 0:
        raise ValueError("No spglib found in dependencies/")
    spglibdir = spglibs[0]

    # set rest of spglib
    spgsrcdir = os.path.join(spglibdir, "src")
    include_dirs = [spgsrcdir]
    sources = ["cell.c", "debug.c", "hall_symbol.c", "kpoint.c", "lattice.c",
               "mathfunc.c", "pointgroup.c", "primitive.c", "refinement.c",
               "sitesym_database.c", "site_symmetry.c", "spacegroup.c", "spin.c",
               "spg_database.c", "spglib.c", "symmetry.c"]
    sources = [os.path.join(spgsrcdir, srcfile) for srcfile in sources]
    return Extension("pymatgen._spglib",
                     include_dirs=include_dirs + get_numpy_include_dirs(),
                     sources=[os.path.join(spglibdir, "_spglib.c")] + sources)

with open("README.rst") as f:
    long_desc = f.read()
    ind = long_desc.find("\n")
    long_desc = long_desc[ind+1:]

setup(
    name="pymatgen",
    packages=find_packages(),
    version="2.5.1",
    install_requires=["numpy>=1.5", "pyhull>=1.3.6", "PyCifRW>=3.3",
                      "requests>=1.0"],
    extras_require={"electronic_structure": ["scipy>=0.10"],
                    "plotting": ["matplotlib>=1.1"],
                    "ase_adaptor": ["ase>=3.3"]},
    package_data={"pymatgen.core": ["bond_lengths.json",
                                    "periodic_table.json"],
                  "pymatgen.analysis": ["bvparam_1991.json", "icsd_bv.json"],
                  "pymatgen.io": ["*.cfg"],
                  "pymatgen.entries": ["*.cfg"],
                  "pymatgen.structure_prediction": ["data/*.json"],
                  "pymatgen.vis": ["ElementColorSchemes.cfg"]},
    author="Shyue Ping Ong, Anubhav Jain, Michael Kocher, Geoffroy Hautier,"
    "William Davidson Richards, Stephen Dacek, Dan Gunter, Shreyas Cholia, "
    "Vincent L Chevrier, Rickard Armiento",
    author_email="shyue@mit.edu, anubhavj@mit.edu, mpkocher@lbnl.gov, "
    "geoffroy.hautier@uclouvain.be, wrichard@mit.edu, sdacek@mit.edu, "
    "dkgunter@lbl.gov, scholia@lbl.gov, vincentchevrier@gmail.com, "
    "armiento@mit.edu",
    maintainer="Shyue Ping Ong",
    url="https://github.com/materialsproject/pymatgen/",
    license="MIT",
    description="pymatgen is the Python materials analysis library powering "
                "the Materials Project (www.materialsproject.org).",
    long_description=long_desc,
    keywords=["vasp", "gaussian", "materials", "project",
              "electronic", "structure"],
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
    download_url="https://github.com/materialsproject/pymatgen/tarball/master",
    ext_modules=[get_spglib_ext()],
    scripts=[os.path.join("scripts", f) for f in os.listdir("scripts")]
)
