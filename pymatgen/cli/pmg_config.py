#!/usr/bin/env python
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Implementation for `pmg config` CLI.
"""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys
from argparse import Namespace
from typing import Literal, Dict
from urllib.request import urlretrieve

from monty.serialization import dumpfn, loadfn

from pymatgen.core import OLD_SETTINGS_FILE, SETTINGS_FILE


def setup_cp2k_data(cp2k_data_dirs: list[str]):
    """Setup cp2k basis and potential data directory"""
    data_dir, targetdir = (os.path.abspath(d) for d in cp2k_data_dirs)
    try:
        os.mkdir(targetdir)
    except OSError:
        r = input("Destination directory exists. Continue (y/n)?")
        if r != "y":
            print("Exiting ...")
            sys.exit(0)
    print("Generating pymatgen resource directory for CP2K...")

    import glob
    from ruamel import yaml
    from pymatgen.core import Element
    from pymatgen.io.cp2k.inputs import GthPotential, GaussianTypeOrbitalBasisSet
    from pymatgen.io.cp2k.utils import chunk
    from monty.json import jsanitize

    basis_files = glob.glob(os.path.join(data_dir, "*BASIS*"))
    potential_files = glob.glob(os.path.join(data_dir, "*POTENTIAL*"))

    settings: Dict[str, Dict] = {
        str(el): {
            "potentials": {},
            "basis_sets": {},
        }
        for el in Element
    }

    for potential_file in potential_files:
        print(f"Processing... {potential_file}")
        with open(potential_file, 'rt') as f:
            try:
                chunks = chunk(f.read())
            except:
                continue
        for c in chunks:
            try:
                potential = GthPotential.from_string(c)
                potential.filename = os.path.basename(potential_file)
                potential.version = None
                settings[potential.element.symbol]['potentials'][potential.get_hash()] = jsanitize(potential.dict())
            except:
                # Chunk was readable, but invalid. Mostly likely "N/A" for this potential
                pass

    for basis_file in basis_files:
        print(f"Processing... {basis_file}")
        with open(basis_file, 'rt') as f:
            try:
                chunks = chunk(f.read())
            except:
                continue
        for c in chunks:
            try:
                basis = GaussianTypeOrbitalBasisSet.from_string(c)
                basis.filename = os.path.basename(basis_file)
                settings[basis.element.symbol]['basis_sets'][basis.get_hash()] = jsanitize(basis.dict())
            except:
                # Chunk was readable, but invalid. Mostly likely "N/A" for this potential
                pass

    print("Done processing cp2k data files")

    for el in settings:
        print(f"Writting {el} settings file")
        with open(os.path.join(targetdir, el), "wt") as f:
            yaml.dump(settings.get(el), f, default_flow_style=False)

    print(
        "\n CP2K resource directory generated. It is recommended that you run:"
        f"\n  'pmg config --add PMG_CP2K_DATA_DIR {os.path.abspath(targetdir)}' "
    )
    print(
        "\n It is also recommended that you set the following (with example values):"
        "\n  'pmg config --add PMG_DEFAULT_CP2K_FUNCTIONAL PBE' "
        "\n  'pmg config --add PMG_DEFAULT_CP2K_BASIS_TYPE TZVP-MOLOPT' "
        "\n  'pmg config --add PMG_DEFAULT_CP2K_AUX_BASIS_TYPE pFIT' "
    )
    print("\n Start a new terminal to ensure that your environment variables are properly set.")

def setup_potcars(potcar_dirs: list[str]):
    """Setup POTCAR directories."""
    pspdir, targetdir = (os.path.abspath(d) for d in potcar_dirs)
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
        "potUSPP_GGA": "POT_GGA_US_PW91",
    }

    for parent, subdirs, _files in os.walk(pspdir):
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
                        with subprocess.Popen(["gunzip", dest]) as p:
                            p.communicate()
                    elif ext.upper() in ["BZ2"]:
                        with subprocess.Popen(["bunzip2", dest]) as p:
                            p.communicate()
                    if subdir == "Osmium":
                        subdir = "Os"
                    dest = os.path.join(basedir, f"POTCAR.{subdir}")
                    shutil.move(os.path.join(basedir, "POTCAR"), dest)
                    with subprocess.Popen(["gzip", "-f", dest]) as p:
                        p.communicate()
                except Exception as ex:
                    print(f"An error has occurred. Message is {str(ex)}. Trying to continue... ")

    print(
        "\nPSP resources directory generated. It is recommended that you "
        f"run 'pmg config --add PMG_VASP_PSP_DIR {os.path.abspath(targetdir)}'"
    )
    print("Start a new terminal to ensure that your environment variables are properly set.")


def build_enum(fortran_command: str = "gfortran") -> bool:
    """
    Build enum.

    :param fortran_command:
    """
    cwd = os.getcwd()
    state = True
    try:
        subprocess.call(["git", "clone", "--recursive", "https://github.com/msg-byu/enumlib.git"])
        os.chdir(os.path.join(cwd, "enumlib", "symlib", "src"))
        os.environ["F90"] = fortran_command
        subprocess.call(["make"])
        enumpath = os.path.join(cwd, "enumlib", "src")
        os.chdir(enumpath)
        subprocess.call(["make"])
        for f in ["enum.x", "makestr.x"]:
            subprocess.call(["make", f])
            shutil.copy(f, os.path.join("..", ".."))
    except Exception as ex:
        print(ex)
        state = False
    finally:
        os.chdir(cwd)
        shutil.rmtree("enumlib")
    return state


def build_bader(fortran_command="gfortran"):
    """
    Build bader package.

    :param fortran_command:
    """
    bader_url = "http://theory.cm.utexas.edu/henkelman/code/bader/download/bader.tar.gz"
    cwd = os.getcwd()
    state = True
    try:
        urlretrieve(bader_url, "bader.tar.gz")
        subprocess.call(["tar", "-zxf", "bader.tar.gz"])
        os.chdir("bader")
        subprocess.call(["cp", "makefile.osx_" + fortran_command, "makefile"])
        subprocess.call(["make"])
        shutil.copy("bader", os.path.join("..", "bader_exe"))
        os.chdir("..")
        shutil.rmtree("bader")
        os.remove("bader.tar.gz")
        shutil.move("bader_exe", "bader")
    except Exception as exc:
        print(str(exc))
        state = False
    finally:
        os.chdir(cwd)
    return state


def install_software(install: Literal["enumlib", "bader"]):
    """Install all optional external software."""
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
    if install == "enumlib":
        print("Building enumlib")
        enum = build_enum(fortran_command)
        print()
    elif install == "bader":
        print("Building bader")
        bader = build_bader(fortran_command)
        print()
    if bader or enum:
        print(
            f"Please add {os.path.abspath('.')} to your PATH or move the executables multinum.x, "
            "makestr.x and/or bader to a location in your PATH.\n"
        )


def add_config_var(tokens: list[str], backup_suffix: str) -> None:
    """Add/update keys in .pmgrc.yaml config file."""
    if os.path.exists(SETTINGS_FILE):
        # read and write new config file if exists
        fpath = SETTINGS_FILE
    elif os.path.exists(OLD_SETTINGS_FILE):
        # else use old config file if exists
        fpath = OLD_SETTINGS_FILE
    else:
        # if neither exists, create new config file
        fpath = SETTINGS_FILE
    d = {}
    if os.path.exists(fpath):
        if backup_suffix:
            shutil.copy(fpath, fpath + backup_suffix)
            print(f"Existing {fpath} backed up to {fpath}{backup_suffix}")
        d = loadfn(fpath)
    if len(tokens) % 2 != 0:
        raise ValueError(f"Uneven number {len(tokens)} of tokens passed to pmg config. Needs a value for every key.")
    for key, val in zip(tokens[0::2], tokens[1::2]):
        d[key] = val
    dumpfn(d, fpath)
    print(f"New {fpath} written!")


def configure_pmg(args: Namespace):
    """Handle configure command."""
    if args.potcar_dirs:
        setup_potcars(args.potcar_dirs)
    elif args.install:
        install_software(args.install)
    elif args.var_spec:
        add_config_var(args.var_spec, args.backup)
    elif args.cp2k_data_dirs:
        setup_cp2k_data(args.cp2k_data_dirs)
