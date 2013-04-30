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
import urllib
import os
import shutil


def build_enum():
    enumlib_url = "http://downloads.sourceforge.net/project/enum/enum/enum.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fenum%2Ffiles%2Fenum%2F&ts=1367333150&use_mirror=iweb"
    try:
        os.makedirs("enumlib")
        os.chdir("enumlib")
        urllib.urlretrieve(enumlib_url, "enum.tar.gz")
        subprocess.call(["tar", "-zxf", "enum.tar.gz"])
        os.chdir("celib")
        os.chdir("trunk")
        os.environ["F90"] = "gfortran"
        subprocess.call(["make"])
        os.chdir(os.path.join("..", ".."))
        enumpath = os.path.join("enumlib", "trunk")
        os.chdir(enumpath)
        subprocess.call(["make"])
        for f in ["multienum.x", "makestr.x"]:
            subprocess.call(["make", f])
            shutil.move(f, os.path.join("..", "..", ".."))
        os.chdir(os.path.join("..", "..", ".."))
        shutil.rmtree("enumlib")
        return True
    except Exception as ex:
        print str(ex)
        return False


def build_bader():
    bader_url = "http://theory.cm.utexas.edu/bader/download/bader.tar.gz"
    try:
        urllib.urlretrieve(bader_url, "bader.tar.gz")
        subprocess.call(["tar", "-zxf", "bader.tar.gz"])
        os.chdir("bader")
        subprocess.call(["cp", "makefile.osx_gfortran", "makefile"])
        subprocess.call(["make"])
        shutil.move("bader", os.path.join("..", "bader_exe"))
        os.chdir("..")
        shutil.rmtree("bader")
        shutil.move("bader_exe", "bader")
        return True
    except Exception as ex:
        print str(ex)
        return False

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

enum = False
bader = False

if "-f" in sys.argv:
    for pk in ["matplotlib>1.1", "scipy"]:
        try:
            subprocess.call(["pip", "install", pk])
        except:
            print "Unable to install {}. Skipping...".format(pk)
    try:
        subprocess.call([
            "pip", "install", "-Ivq",
            "https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.0.2515.tar.gz"])
    except:
        print "Unable to install ASE. Skipping..."

    enum = build_enum()
    bader = build_bader()

    try:
        subprocess.call(["potcar_setup.py"])
    except:
        print "Skipping POTCAR setup."

print "------------ Setup complete --------------"
print "You still need to perform a few manual changes."
print
if enum or bader:
    print "Please add {} to your PATH.".format(os.path.abspath("."))
print

print "To use the Materials API, get your Materials API key at " \
      "https://www.materialsproject.org/profile and add it to your " \
      "environment"
print "export MAPI_KEY=YOUR_API_KEY"
