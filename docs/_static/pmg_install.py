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

from __future__ import print_function

import sys
import subprocess
import urllib
import os
import shutil

__author__ = "Shyue Ping Ong"
__version__ = "1.1"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2013"




def build_enum(fortran_command="gfortran"):
    enumlib_url = "http://downloads.sourceforge.net/project/enum/enum/enum.tar.gz"
    currdir = os.getcwd()
    state = True
    try:
        os.makedirs("enumlib")
        os.chdir("enumlib")
        urllib.urlretrieve(enumlib_url, "enum.tar.gz")
        subprocess.call(["tar", "-zxf", "enum.tar.gz"])
        os.chdir("celib")
        os.chdir("trunk")
        os.environ["F90"] = fortran_command
        subprocess.call(["make"])
        os.chdir(os.path.join("..", ".."))
        enumpath = os.path.join("enumlib", "trunk")
        os.chdir(enumpath)
        subprocess.call(["make"])
        for f in ["multienum.x", "makestr.x"]:
            subprocess.call(["make", f])
            shutil.copy(f, os.path.join("..", "..", ".."))
    except Exception as ex:
        print(str(ex))
        state = False
    finally:
        os.chdir(currdir)
        shutil.rmtree("enumlib")
    return state


def build_bader(fortran_command="gfortran"):
    bader_url = "http://theory.cm.utexas.edu/henkelman/code/bader/download/bader.tar.gz"
    currdir = os.getcwd()
    state = True
    try:
        urllib.urlretrieve(bader_url, "bader.tar.gz")
        subprocess.call(["tar", "-zxf", "bader.tar.gz"])
        os.chdir("bader")
        subprocess.call(["cp", "makefile.osx_"+fortran_command, "makefile"])
        subprocess.call(["make"])
        shutil.copy("bader", os.path.join("..", "bader_exe"))
        os.chdir("..")
        shutil.rmtree("bader")
        os.remove("bader.tar.gz")
        shutil.move("bader_exe", "bader")
    except Exception as ex:
        print(str(ex))
        state = False
    finally:
        os.chdir(currdir)
    return state

try:
    py_ver = sys.version_info
    print("Detected Python version %s" % ".".join(["%s" % i for i in py_ver]))
    if py_ver < (2, 7) or py_ver >= (2, 8):
        print("Python version 2.7+ required. Download and install the necessary "
              "python version from http://www.python.org/download/.")
        sys.exit(-1)
except:
    print("Python version 2.7+ required. Download and install the necessary "
          "python version from http://www.python.org/download/.")
    sys.exit(-1)
    
try:
    import setuptools
    print("Detected setuptools version {}".format(setuptools.__version__))
except ImportError:
    print("setuptools not detected. Get it from https://pypi.python"
          ".org/pypi/setuptools and follow the instructions to install first.")
    sys.exit(-1)

try:
    gcc_ver = subprocess.Popen(["gcc", "--version"], stdout=subprocess.PIPE)\
        .communicate()[0]
except:
    print("gcc not found in PATH. gcc is needed for installation of numpy "
          "and C extensions. For Mac users, please install Xcode and its "
          "corresponding command-line tools first.")
    sys.exit(-1)

try:
    import pip
    print("Detected pip version {}".format(pip.__version__))
except (AttributeError, ImportError):
    print("pip not detected. Installing...")
    subprocess.call(["easy_install", "pip"])

try:
    import numpy
    from numpy.distutils.misc_util import get_numpy_include_dirs
    print("Detected numpy version {}".format(numpy.__version__))
except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Installing...")
    subprocess.call(["pip", "install", "-q", "numpy>=1.8.0"])
    from numpy.distutils.misc_util import get_numpy_include_dirs

if subprocess.call(["pip", "install", "pymatgen"]) != 0:
    print("Error installing pymatgen")
    sys.exit(-1)
print("")

enum = False
bader = False

if "-f" in sys.argv:
    for pk in ["matplotlib>1.5"]:
        if subprocess.call(["pip", "install", pk]) != 0:
            print("Unable to install {}. Skipping...".format(pk))

    fortran_command = None
    try:
        if subprocess.call(["ifort", "--version"]) == 0:
            print("Found ifort")
            fortran_command = "ifort"
        elif subprocess.call(["gfortran", "--version"]) == 0:
            print("Found gfortran")
            fortran_command = "gfortran"
    except:
        fortran_command = None

    if fortran_command is not None:
        print("Building enumlib")
        enum = build_enum(fortran_command)
        print("")
        print("Building bader")
        bader = build_bader(fortran_command)
        print("")
    else:
        print("No fortran compiler found. Skipping enumlib and bader build.")

    print("Performing POTCAR setup. Press Ctrl-C at any prompt to skip this "
          "step.")
    try:
        subprocess.call(["potcar_setup"])
    except:
        print("Skipping POTCAR setup.")
    print("")

print("------------ Setup complete --------------")
print("You still need to perform a few manual changes.")
print("")
if enum or bader:
    print("Please add {} to your PATH or move the executables multinum.x, "
          "makestr.x and bader to a location in your PATH."
          .format(os.path.abspath(".")))
    print("")

print("To use the Materials API, get your Materials API key at "
      "https://www.materialsproject.org/profile and add it to your "
      "environment")
print("export MAPI_KEY=YOUR_API_KEY")
print("")