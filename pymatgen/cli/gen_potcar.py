#!/usr/bin/env python
# coding: utf-8

import os

from pymatgen.io.vasp import Potcar


def proc_dir(dirname, procfilefunction):
    for f in os.listdir(dirname):
        if os.path.isdir(os.path.join(dirname, f)):
            proc_dir(os.path.join(dirname, f), procfilefunction)
        else:
            procfilefunction(dirname, f)


def gen_potcar(dirname, filename):
    if filename == "POTCAR.spec":
        fullpath = os.path.join(dirname, filename)
        f = open(fullpath, "r")
        elements = f.readlines()
        f.close()
        symbols = [el.strip() for el in elements if el.strip() != ""]
        potcar = Potcar(symbols)
        potcar.write_file(os.path.join(dirname, "POTCAR"))


if __name__ == "__main__":
    proc_dir(os.getcwd(), gen_potcar)
