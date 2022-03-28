#!/usr/bin/env python

from os.path import abspath, dirname

from pymatgen.io.vasp.outputs import Vasprun

run = Vasprun(f"{dirname(abspath(__file__))}/vasprun.xml")
