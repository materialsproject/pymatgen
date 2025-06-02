#!/usr/bin/env python

"""Implementation for `pmg potcar` CLI."""

from __future__ import annotations

import os

from pymatgen.io.vasp import Potcar


def proc_dir(dirname, proc_file_function):
    """Process a directory.

    Args:
        dirname (str): Directory name.
        proc_file_function (callable): Callable to execute on directory.
    """
    for file in os.listdir(dirname):
        if os.path.isdir(os.path.join(dirname, file)):
            proc_dir(os.path.join(dirname, file), proc_file_function)
        else:
            proc_file_function(dirname, file)


def gen_potcar(dirname, filename):
    """Generate POTCAR from POTCAR.spec in directories.

    Args:
        dirname (str): Directory name.
        filename (str): Filename in directory.
    """
    if filename == "POTCAR.spec":
        fullpath = os.path.join(dirname, filename)
        with open(fullpath, encoding="utf-8") as file:
            elements = file.readlines()
        symbols = [el.strip() for el in elements if el.strip() != ""]
        potcar = Potcar(symbols)
        potcar.write_file(f"{dirname}/POTCAR")


def generate_potcar(args):
    """Generate POTCAR.

    Args:
        args (dict): Args from argparse.
    """
    if args.recursive:
        proc_dir(args.recursive, gen_potcar)
    elif args.symbols:
        try:
            p = Potcar(args.symbols, functional=args.functional)
            p.write_file("POTCAR")
        except Exception as exc:
            print(f"An error has occurred: {exc}")

    else:
        print("No valid options selected.")


if __name__ == "__main__":
    proc_dir(os.getcwd(), gen_potcar)
