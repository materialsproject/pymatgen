#!/usr/bin/env python

"""Implementation for `pmg config` CLI."""

from __future__ import annotations

import os
import shutil
import subprocess
from glob import glob
from typing import TYPE_CHECKING
from urllib.request import urlretrieve

from monty.json import jsanitize
from monty.serialization import dumpfn, loadfn

from pymatgen.core import OLD_SETTINGS_FILE, SETTINGS_FILE, Element
from pymatgen.io.cp2k.inputs import GaussianTypeOrbitalBasisSet, GthPotential
from pymatgen.io.cp2k.utils import chunk
from pymatgen.io.vasp.inputs import POTCAR_FUNCTIONAL_MAP

if TYPE_CHECKING:
    from argparse import Namespace
    from typing import Literal


def setup_cp2k_data(cp2k_data_dirs: list[str]) -> None:
    """Setup CP2K basis and potential data directory."""
    # this function used to use ruamel.yaml which underwent breaking changes. was easier to
    # migrate to PyYAML than fix
    import yaml  # type: ignore[import]

    data_dir, target_dir = (os.path.abspath(dirc) for dirc in cp2k_data_dirs)
    try:
        os.mkdir(target_dir)
    except OSError:
        reply = input("Destination directory exists. Continue (y/n)?")
        if reply != "y":
            raise SystemExit("Exiting ...")
    print("Generating pymatgen resource directory for CP2K...")

    basis_files = glob(f"{data_dir}/*BASIS*")
    potential_files = glob(f"{data_dir}/*POTENTIAL*")

    settings: dict[str, dict] = {str(el): {"potentials": {}, "basis_sets": {}} for el in Element}

    for potential_file in potential_files:
        print(f"Processing... {potential_file}")
        with open(potential_file, encoding="utf-8") as file:
            try:
                chunks = chunk(file.read())
            except IndexError:
                continue
        for chk in chunks:
            try:
                potential = GthPotential.from_str(chk)
                potential.filename = os.path.basename(potential_file)
                potential.version = None
                if potential.element is not None:
                    settings[potential.element.symbol]["potentials"][potential.get_hash()] = jsanitize(
                        potential, strict=True
                    )
            except ValueError:
                # Chunk was readable, but the element is not pmg recognized
                continue
            except IndexError:
                # Chunk was readable, but invalid. Mostly likely "N/A" for this potential
                continue

    for basis_file in basis_files:
        print(f"Processing... {basis_file}")
        with open(basis_file, encoding="utf-8") as file:
            try:
                chunks = chunk(file.read())
            except IndexError:
                continue
        for chk in chunks:
            try:
                basis = GaussianTypeOrbitalBasisSet.from_str(chk)
                basis.filename = os.path.basename(basis_file)
                settings[basis.element.symbol]["basis_sets"][basis.get_hash()] = jsanitize(basis, strict=True)  # type: ignore[union-attr]
            except ValueError:
                # Chunk was readable, but the element is not pmg recognized
                continue
            except IndexError:
                # Chunk was readable, but invalid. Mostly likely "N/A" for this potential
                continue

    print("Done processing cp2k data files")

    for el in settings:
        print(f"Writing {el} settings file")
        with open(os.path.join(target_dir, el), mode="w", encoding="utf-8") as file:
            yaml.dump(settings.get(el), file, default_flow_style=False)

    print(
        "\n CP2K resource directory generated. It is recommended that you run:"
        f"\n  'pmg config --add PMG_CP2K_DATA_DIR {os.path.abspath(target_dir)}' "
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
    psp_dir, target_dir = (os.path.abspath(d) for d in potcar_dirs)
    try:
        os.makedirs(target_dir)
    except OSError:
        reply = input("Destination directory exists. Continue (y/n)? ")
        if reply != "y":
            raise SystemExit("Exiting ...")

    print("Generating pymatgen resources directory...")

    # name_mappings ensures consistency with the POTCAR directory structure
    # expected by pymatgen.io.vasp.inputs.PotcarSingle
    name_mappings = {
        k: POTCAR_FUNCTIONAL_MAP[v]
        for k, v in {
            "potpaw_PBE": "PBE",
            "potpaw_PBE_52": "PBE_52",
            "potpaw_PBE.52": "PBE_52",
            "potpaw_PBE_54": "PBE_54",
            "potpaw_PBE.54": "PBE_54",
            "potpaw_PBE_64": "PBE_64",
            "potpaw_PBE.64": "PBE_64",
            "potpaw_LDA": "LDA",
            "potpaw_LDA.52": "LDA_52",
            "potpaw_LDA_52": "LDA_52",
            "potpaw_LDA.54": "LDA_54",
            "potpaw_LDA_54": "LDA_54",
            "potpaw_LDA_64": "LDA_64",
            "potpaw_LDA.64": "LDA_64",
            "potpaw_GGA": "PW91",
            "potUSPP_LDA": "LDA_US",
            "potUSPP_GGA": "PW91_US",
        }.items()
    }

    for parent, subdirs, _files in os.walk(psp_dir):
        basename = os.path.basename(parent)
        basename = name_mappings.get(basename, basename)
        for subdir in subdirs:
            filenames = glob(os.path.join(parent, subdir, "POTCAR*"))
            if len(filenames) > 0:
                try:
                    base_dir = os.path.join(target_dir, basename)
                    os.makedirs(base_dir, exist_ok=True)
                    fname = filenames[0]
                    dest = os.path.join(base_dir, os.path.basename(fname))
                    shutil.copy(fname, dest)
                    ext = fname.split(".")[-1]
                    if ext.upper() in {"Z", "GZ"}:
                        with subprocess.Popen(["gunzip", dest]) as process:
                            process.communicate()
                    elif ext.upper() == "BZ2":
                        with subprocess.Popen(["bunzip2", dest]) as process:
                            process.communicate()
                    if subdir == "Osmium":
                        subdir = "Os"
                    dest = os.path.join(base_dir, f"POTCAR.{subdir}")
                    shutil.move(f"{base_dir}/POTCAR", dest)
                    with subprocess.Popen(["gzip", "-f", dest]) as process:
                        process.communicate()
                except Exception as exc:
                    print(f"An error has occurred. Message is {exc}. Trying to continue... ")

    print(
        "\nPSP resources directory generated. It is recommended that you "
        f"run 'pmg config --add PMG_VASP_PSP_DIR {os.path.abspath(target_dir)}'"
    )
    print("Start a new terminal to ensure that your environment variables are properly set.")


def build_enum(fortran_command: str = "gfortran") -> bool:
    """Build enum.

    Args:
        fortran_command: The Fortran compiler command.
    """
    cwd = os.getcwd()
    state = True
    try:
        subprocess.call(["git", "clone", "--recursive", "https://github.com/msg-byu/enumlib"])
        os.chdir(f"{cwd}/enumlib/symlib/src")
        os.environ["F90"] = fortran_command
        subprocess.call(["make"])
        enum_path = f"{cwd}/enumlib/src"
        os.chdir(enum_path)
        subprocess.call(["make"])
        subprocess.call(["make", "enum.x"])
        shutil.copy("enum.x", os.path.join("..", ".."))
    except Exception as exc:
        print(exc)
        state = False
    finally:
        os.chdir(cwd)
        shutil.rmtree("enumlib")
    return state


def build_bader(fortran_command="gfortran"):
    """Build bader package.

    Args:
        fortran_command: The Fortran compiler command.
    """
    bader_url = "https://theory.cm.utexas.edu/henkelman/code/bader/download/bader.tar.gz"
    cwd = os.getcwd()
    state = True
    try:
        urlretrieve(bader_url, "bader.tar.gz")  # noqa: S310
        subprocess.call(["/usr/bin/tar", "-zxf", "bader.tar.gz"])
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
        except Exception as exc:
            print(str(exc))
            raise SystemExit("No fortran compiler found.")

    enum = bader = None
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
    if len(tokens) % 2 != 0:
        raise ValueError(f"Uneven number {len(tokens)} of tokens passed to pmg config. Needs a value for every key.")
    if os.path.isfile(SETTINGS_FILE):
        # read and write new config file if exists
        rc_path = SETTINGS_FILE
    elif os.path.isfile(OLD_SETTINGS_FILE):
        # else use old config file if exists
        rc_path = OLD_SETTINGS_FILE
    else:
        # if neither exists, create new config file
        rc_path = SETTINGS_FILE
    dct = {}
    if os.path.isfile(rc_path):
        if backup_suffix:
            shutil.copy(rc_path, rc_path + backup_suffix)
            print(f"Existing {rc_path} backed up to {rc_path}{backup_suffix}")
        dct = loadfn(rc_path)
    special_vals = {"true": True, "false": False, "none": None, "null": None}
    for key, val in zip(tokens[::2], tokens[1::2], strict=True):
        dct[key] = special_vals.get(val.lower(), val)
    dumpfn(dct, rc_path)
    print(f"New {rc_path} written!")


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
