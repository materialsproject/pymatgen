"""
This package implements modules for input and output to and from VASP. It
imports the key classes form both vasp_input and vasp_output to allow most
classes to be simply called as pymatgen.io.vasp.Incar for example, to retain
backwards compatibility.
"""

from __future__ import annotations

from .inputs import Incar, Kpoints, Poscar, Potcar, PotcarSingle, VaspInput
from .outputs import (
    BSVasprun,
    Chgcar,
    Dynmat,
    Elfcar,
    Locpot,
    Oszicar,
    Outcar,
    Procar,
    Vasprun,
    VolumetricData,
    Wavecar,
    Waveder,
    Xdatcar,
)
