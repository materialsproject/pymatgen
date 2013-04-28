#!/usr/bin/env python

"""
This is a simple installation script for casual users of pymatgen who simply
plan to use pymatgen as a basic analysis library and is not planning to
develop on it. This script should work on most Linux and Mac systems that
have Python 2.7+ installed and setuptools installed. These are the only
required pre-requisites. Once those are installed, the script should take
care of the remainder of the installation process.

There are only a few options in this script. Please note that you probably
have to *run all commands with sudo* for the installation to proceed correctly.

Simply running:

    ./pmg_install

will install pymatgen with the basic dependencies.

Running:

    ./pmg_install -f

will install pymatgen with a few more optional packages and also start an
initial setup process that guides you through basic configuration
for POTCAR and Materials API support.

Report any issues or suggestions for this script to shyuep@gmail.com.
"""

__author__ = "Shyue Ping Ong"
__version__ = "1.0"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2013"

import sys
import subprocess


py_ver = sys.version_info
print "Detected Python version {}".format(".".join(map(str, py_ver)))


if py_ver < (2, 7) or py_ver >= (2, 8):
    print "Python version 2.7+ required. Download and install the necessary " \
          "python version from http://www.python.org/download/."
    sys.exit(-1)

try:
    import setuptools
    print "Detected setuptools version {}".format(setuptools.__version__)
except ImportError:
    print "setuptools not detected. Get it from https://pypi.python" \
          ".org/pypi/setuptools and follow the instructions to install first."
    sys.exit(-1)


try:
    gcc_ver = subprocess.Popen(["gcc", "--version"], stdout=subprocess.PIPE)\
        .communicate()[0]
except:
    print "gcc not found in PATH. gcc is needed for installation of numpy " \
          "and C extensions. For Mac users, please install Xcode and its " \
          "corresponding command-line tools first."
    sys.exit(-1)

try:
    import pip
    print "Detected pip version {}".format(pip.__version__)
except ImportError:
    print "pip not detected. Installing..."
    subprocess.call(["easy_install", "pip"])

try:
    from numpy.distutils.misc_util import get_numpy_include_dirs
except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Installing...")
    subprocess.call(["pip", "install", "numpy>=1.6.0"])

for pk in ["pyhull>=1.3.6", "PyCifRW>=3.3", "requests>=1.0", "pybtex>=0.16"]:
    subprocess.call(["pip", "install", pk])

subprocess.call(["pip", "install", "pymatgen"])

if "-f" in sys.argv:
    for pk in ["matplotlib>1.1", "scipy"]:
        try:
            subprocess.call(["pip", "install", pk])
        except:
            print "Unable to install {}. Skipping...".format(pk)
    try:
        subprocess.call(["pip", "install", "-Iv",
                         "https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.0"
                         ".2515.tar.gz"])
    except:
        print "Unable to install ASE. Skipping..."

    subprocess.call(["potcar_setup.py"])

    print "To use the Materials API, get your Materials API key at " \
          "https://www.materialsproject.org/profile and add it to your " \
          "environment"
    print "export MAPI_KEY=YOUR_API_KEY"
