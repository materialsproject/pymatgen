#!/usr/bin/env python
# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
A master convenience script with many tools for vasp and structure analysis.
"""


import os
import sys
import glob
import shutil
import subprocess
from monty.serialization import loadfn, dumpfn

from pymatgen import SETTINGS_FILE
from urllib.request import urlretrieve


def setup_potcars(args):
    """
    Setup POTCAR directirt,

    :param args: args from command.
    """
    pspdir, targetdir = [os.path.abspath(d) for d in args.potcar_dirs]
    try:
        os.makedirs(targetdir)
    except OSError:
        r = input("Destination directory exists. Continue (y/n)? ")
        if r != "y":
            print("Exiting ...")
            sys.exit(0)

    print("Generating pymatgen resources directory...")

    name_mappings = {
        "potpaw_PBE": "POT_GGA_PAW_PBE",
        "potpaw_PBE_52": "POT_GGA_PAW_PBE_52",
        "potpaw_PBE_54": "POT_GGA_PAW_PBE_54",
        "potpaw_PBE.52": "POT_GGA_PAW_PBE_52",
        "potpaw_PBE.54": "POT_GGA_PAW_PBE_54",
        "potpaw_LDA": "POT_LDA_PAW",
        "potpaw_LDA.52": "POT_LDA_PAW_52",
        "potpaw_LDA.54": "POT_LDA_PAW_54",
        "potpaw_LDA_52": "POT_LDA_PAW_52",
        "potpaw_LDA_54": "POT_LDA_PAW_54",
        "potUSPP_LDA": "POT_LDA_US",
        "potpaw_GGA": "POT_GGA_PAW_PW91",
        "potUSPP_GGA": "POT_GGA_US_PW91"
    }

    for (parent, subdirs, files) in os.walk(pspdir):
        basename = os.path.basename(parent)
        basename = name_mappings.get(basename, basename)
        for subdir in subdirs:
            filenames = glob.glob(os.path.join(parent, subdir, "POTCAR*"))
            if len(filenames) > 0:
                try:
                    basedir = os.path.join(targetdir, basename)
                    if not os.path.exists(basedir):
                        os.makedirs(basedir)
                    fname = filenames[0]
                    dest = os.path.join(basedir, os.path.basename(fname))
                    shutil.copy(fname, dest)
                    ext = fname.split(".")[-1]
                    if ext.upper() in ["Z", "GZ"]:
                        subprocess.Popen(["gunzip", dest]).communicate()
                    elif ext.upper() in ["BZ2"]:
                        subprocess.Popen(["bunzip2", dest]).communicate()
                    if subdir == "Osmium":
                        subdir = "Os"
                    dest = os.path.join(basedir, "POTCAR.{}".format(subdir))
                    shutil.move(os.path.join(basedir, "POTCAR"), dest)
                    subprocess.Popen(["gzip", "-f", dest]).communicate()
                except Exception as ex:
                    print("An error has occured. Message is %s. Trying to "
                          "continue... " % str(ex))

    print("")
    print("PSP resources directory generated. It is recommended that you "
          "run 'pmg config --add PMG_VASP_PSP_DIR %s'" % os.path.abspath(targetdir))
    print("Start a new terminal to ensure that your environment variables "
          "are properly set.")


def build_enum(fortran_command="gfortran"):
    """
    Build enum.

    :param fortran_command:
    """
    currdir = os.getcwd()
    state = True
    try:
        subprocess.call(["git", "clone", "--recursive",
                         "https://github.com/msg-byu/enumlib.git"])
        os.chdir(os.path.join(currdir, "enumlib", "symlib", "src"))
        os.environ["F90"] = fortran_command
        subprocess.call(["make"])
        enumpath = os.path.join(currdir, "enumlib", "src")
        os.chdir(enumpath)
        subprocess.call(["make"])
        for f in ["enum.x", "makestr.x"]:
            subprocess.call(["make", f])
            shutil.copy(f, os.path.join("..", ".."))
    except Exception as ex:
        print(str(ex))
        state = False
    finally:
        os.chdir(currdir)
        shutil.rmtree("enumlib")
    return state


def build_bader(fortran_command="gfortran"):
    """
    Build bader package.

    :param fortran_command:
    """
    bader_url = "http://theory.cm.utexas.edu/henkelman/code/bader/download/bader.tar.gz"
    currdir = os.getcwd()
    state = True
    try:
        urlretrieve(bader_url, "bader.tar.gz")
        subprocess.call(["tar", "-zxf", "bader.tar.gz"])
        os.chdir("bader")
        subprocess.call(
            ["cp", "makefile.osx_" + fortran_command, "makefile"])
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


def install_software(args):
    """
    Install all optional external software.

    :param args:
    """
    try:
        subprocess.call(["ifort", "--version"])
        print("Found ifort")
        fortran_command = "ifort"
    except Exception:
        try:
            subprocess.call(["gfortran", "--version"])
            print("Found gfortran")
            fortran_command = "gfortran"
        except Exception as ex:
            print(str(ex))
            print("No fortran compiler found.")
            sys.exit(-1)

    enum = None
    bader = None
    if args.install == "enumlib":
        print("Building enumlib")
        enum = build_enum(fortran_command)
        print("")
    elif args.install == "bader":
        print("Building bader")
        bader = build_bader(fortran_command)
        print("")
    if bader or enum:
        print("Please add {} to your PATH or move the executables multinum.x, "
              "makestr.x and/or bader to a location in your PATH."
              .format(os.path.abspath(".")))
        print("")


def add_config_var(args):
    """
    Add configuration args.

    :param args:
    """
    d = {}
    if os.path.exists(SETTINGS_FILE):
        shutil.copy(SETTINGS_FILE, SETTINGS_FILE + ".bak")
        print("Existing %s backed up to %s"
              % (SETTINGS_FILE, SETTINGS_FILE + ".bak"))
        d = loadfn(SETTINGS_FILE)
    toks = args.var_spec
    if len(toks) % 2 != 0:
        print("Bad variable specification!")
        sys.exit(-1)
    for i in range(int(len(toks) / 2)):
        d[toks[2 * i]] = toks[2 * i + 1]
    dumpfn(d, SETTINGS_FILE, default_flow_style=False)
    print("New %s written!" % (SETTINGS_FILE))


def configure_pmg(args):
    """
    Handle configure command.

    :param args:
    """
    if args.potcar_dirs:
        setup_potcars(args)
    elif args.install:
        install_software(args)
    elif args.var_spec:
        add_config_var(args)
