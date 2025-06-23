"""
This module defines classes for parsing the FEFF output files.

Currently supports the xmu.dat, ldos.dat output files are for non-spin case.
"""

from __future__ import annotations

import re
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core import Element
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.io.feff import Header, Potential, Tags

if TYPE_CHECKING:
    from typing_extensions import Self

__author__ = "Alan Dozier, Kiran Mathew, Chen Zheng"
__credits__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0.3"
__maintainer__ = "Alan Dozier"
__email__ = "adozier@uky.edu"
__status__ = "Beta"
__date__ = "April 7, 2013"


class LDos(MSONable):
    """Parser for ldos files ldos01, ldos02, ..."""

    def __init__(self, complete_dos, charge_transfer):
        """
        Args:
            complete_dos (CompleteDos): complete dos object
            charge_transfer (dict): computed charge transfer between atoms
                dictionary.
        """
        self.complete_dos = complete_dos
        self.charge_transfer = charge_transfer

    @classmethod
    def from_file(cls, feff_inp_file: str = "feff.inp", ldos_file: str = "ldos") -> Self:
        """
        Creates LDos object from raw Feff ldos files by
        by assuming they are numbered consecutively, i.e. ldos01.dat
        ldos02.dat...

        Args:
            feff_inp_file (str): input file of run to obtain structure
            ldos_file (str): output ldos file of run to obtain dos info, etc.
        """
        header_str = Header.header_string_from_file(feff_inp_file)
        header = Header.from_str(header_str)
        structure = header.struct
        n_sites = len(structure)
        parameters = Tags.from_file(feff_inp_file)

        if "RECIPROCAL" in parameters:
            pot_dict = {}
            pot_read_start = re.compile(".*iz.*lmaxsc.*xnatph.*xion.*folp.*")
            pot_read_end = re.compile(".*ExternalPot.*switch.*")
            pot_inp = re.sub(r"feff.inp", r"pot.inp", feff_inp_file)
            dos_index = 1
            begin = 0

            with zopen(pot_inp, mode="rt", encoding="utf-8") as potfile:
                line: str
                for line in potfile:  # type:ignore[assignment]
                    if len(pot_read_end.findall(line)) > 0:
                        break

                    if begin == 1:
                        begin += 1
                        continue

                    if begin == 2:
                        z_number = int(line.strip().split()[0])
                        ele_name = Element.from_Z(z_number).name
                        if ele_name not in pot_dict:
                            pot_dict[ele_name] = dos_index
                        else:
                            pot_dict[ele_name] = min(dos_index, pot_dict[ele_name])
                        dos_index += 1

                    if len(pot_read_start.findall(line)) > 0:
                        begin = 1
        else:
            pot_string = Potential.pot_string_from_file(feff_inp_file)
            dicts = Potential.pot_dict_from_str(pot_string)
            pot_dict = dicts[0]

        with zopen(f"{ldos_file}00.dat", mode="rt", encoding="utf-8") as file:
            lines = file.readlines()
        e_fermi = float(lines[0].split()[4])

        dos_energies = []
        ldos = {}

        for idx in range(1, len(pot_dict) + 1):
            if len(str(idx)) == 1:
                ldos[idx] = np.loadtxt(f"{ldos_file}0{idx}.dat")
            else:
                ldos[idx] = np.loadtxt(f"{ldos_file}{idx}.dat")

        for idx in range(len(ldos[1])):
            dos_energies.append(ldos[1][idx][0])

        all_pdos: list[dict] = []
        vorb = {"s": Orbital.s, "p": Orbital.py, "d": Orbital.dxy, "f": Orbital.f0}
        forb = {"s": 0, "p": 1, "d": 2, "f": 3}

        d_length = len(ldos[1])

        for idx in range(n_sites):
            pot_index = pot_dict[structure.species[idx].symbol]
            all_pdos.append(defaultdict(dict))
            for k, v in vorb.items():
                density = [ldos[pot_index][j][forb[k] + 1] for j in range(d_length)]
                up_dos = density
                if down_dos := None:
                    all_pdos[-1][v] = {Spin.up: up_dos, Spin.down: down_dos}
                else:
                    all_pdos[-1][v] = {Spin.up: up_dos}

        pdos = all_pdos
        vorb2 = {0: Orbital.s, 1: Orbital.py, 2: Orbital.dxy, 3: Orbital.f0}
        pdoss = {structure[i]: {v: pdos[i][v] for v in vorb2.values()} for i in range(len(pdos))}

        forb = {"s": 0, "p": 1, "d": 2, "f": 3}

        t_dos = [0] * d_length
        for idx in range(n_sites):
            pot_index = pot_dict[structure.species[idx].symbol]
            for forb_val in forb.values():
                density = [ldos[pot_index][j][forb_val + 1] for j in range(d_length)]
                for j in range(d_length):
                    t_dos[j] += density[j]
        _t_dos: dict = {Spin.up: t_dos}

        dos = Dos(e_fermi, dos_energies, _t_dos)
        complete_dos = CompleteDos(structure, dos, pdoss)  # type:ignore[arg-type]
        charge_transfer = LDos.charge_transfer_from_file(feff_inp_file, ldos_file)
        return cls(complete_dos, charge_transfer)

    @staticmethod
    def charge_transfer_from_file(feff_inp_file, ldos_file):
        """Get charge transfer from file.

        Args:
            feff_inp_file (str): name of feff.inp file for run
            ldos_file (str): ldos filename for run, assume consecutive order,
                i.e., ldos01.dat, ldos02.dat....

        Returns:
            dictionary of dictionaries in order of potential sites
            ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        """
        cht = {}
        parameters = Tags.from_file(feff_inp_file)

        if "RECIPROCAL" in parameters:
            dicts = [{}]
            pot_dict = {}
            dos_index = 1
            begin = 0
            pot_inp = re.sub(r"feff.inp", r"pot.inp", feff_inp_file)
            pot_readstart = re.compile(".*iz.*lmaxsc.*xnatph.*xion.*folp.*")
            pot_readend = re.compile(".*ExternalPot.*switch.*")
            with zopen(pot_inp, mode="rt", encoding="utf-8") as potfile:
                for line in potfile:
                    if len(pot_readend.findall(line)) > 0:
                        break
                    if begin == 1:
                        z_number = int(line.strip().split()[0])
                        ele_name = Element.from_Z(z_number).name
                        if len(pot_dict) == 0:
                            pot_dict[0] = ele_name
                        elif len(pot_dict) > 0:
                            pot_dict[max(pot_dict) + 1] = ele_name
                        begin += 1
                        continue
                    if begin == 2:
                        z_number = int(line.strip().split()[0])
                        ele_name = Element.from_Z(z_number).name
                        dicts[0][ele_name] = dos_index
                        dos_index += 1
                        if len(pot_dict) == 0:
                            pot_dict[0] = ele_name
                        elif len(pot_dict) > 0:
                            pot_dict[max(pot_dict) + 1] = ele_name
                    if len(pot_readstart.findall(line)) > 0:
                        begin = 1
        else:
            pot_string = Potential.pot_string_from_file(feff_inp_file)
            dicts = Potential.pot_dict_from_str(pot_string)
            pot_dict = dicts[1]

        for idx in range(len(dicts[0]) + 1):
            if len(str(idx)) == 1:
                with zopen(f"{ldos_file}0{idx}.dat", mode="rt", encoding="utf-8") as file:
                    lines = file.readlines()
                    s = float(lines[3].split()[2])
                    p = float(lines[4].split()[2])
                    d = float(lines[5].split()[2])
                    f1 = float(lines[6].split()[2])
                    tot = float(lines[1].split()[4])
                    cht[str(idx)] = {pot_dict[idx]: {"s": s, "p": p, "d": d, "f": f1, "tot": tot}}
            else:
                with zopen(f"{ldos_file}{idx}.dat", mode="rt", encoding="utf-8") as file:
                    lines = file.readlines()
                    s = float(lines[3].split()[2])
                    p = float(lines[4].split()[2])
                    d = float(lines[5].split()[2])
                    f1 = float(lines[6].split()[2])
                    tot = float(lines[1].split()[4])
                    cht[str(idx)] = {pot_dict[idx]: {"s": s, "p": p, "d": d, "f": f1, "tot": tot}}

        return cht

    def charge_transfer_to_str(self):
        """Get charge transfer as string."""
        ch = self.charge_transfer
        chts = ["\nCharge Transfer\n\nabsorbing atom"]
        for i in range(len(ch)):
            for atom, v2 in ch[str(i)].items():
                a = [
                    "\n",
                    atom,
                    "\n",
                    "s   ",
                    str(v2["s"]),
                    "\n",
                    "p   ",
                    str(v2["p"]),
                    "\n",
                    "d   ",
                    str(v2["d"]),
                    "\n",
                    "f   ",
                    str(v2["f"]),
                    "\n",
                    "tot ",
                    str(v2["tot"]),
                    "\n",
                ]
                chts.extend(a)
        return "".join(chts)


class Xmu(MSONable):
    r"""
    Parser for data in 'xmu.dat' file.
    The file 'xmu.dat' contains XANES, EXAFS or NRIXS data depending on the
    situation; \mu, \mu_0, and \chi = \chi * \mu_0 / \mu_0 / (edge+50eV) as
    functions of absolute energy E, relative energy E - E_f and wave number k.

    Default attributes:
        xmu: Photon absorption cross section of absorbing atom in material
        Energies: Energies of data point
        relative_energies: E - E_fermi
        wavenumber: k=\sqrt(E -E_fermi)
        mu: The total absorption cross-section.
        mu0: The embedded atomic background absorption.
        chi: fine structure.
        Edge: Absorption Edge
        Absorbing atom: Species of absorbing atom
        Material: Formula of material
        Source: Source of structure
        Calculation: Type of Feff calculation performed
    """

    def __init__(self, header, parameters, absorbing_atom, data):
        """
        Args:
            header: Header object
            parameters: Tags object
            absorbing_atom (str | int): absorbing atom symbol or index
            data (numpy.ndarray, Nx6): cross_sections.
        """
        self.header = header
        self.parameters = parameters
        self.absorbing_atom = absorbing_atom
        self.data = np.array(data)

    @classmethod
    def from_file(cls, xmu_dat_file: str = "xmu.dat", feff_inp_file: str = "feff.inp") -> Self:
        """Get Xmu from file.

        Args:
            xmu_dat_file (str): filename and path for xmu.dat
            feff_inp_file (str): filename and path of feff.inp input file

        Returns:
            Xmu object
        """
        data = np.loadtxt(xmu_dat_file)
        header = Header.from_file(feff_inp_file)
        parameters = Tags.from_file(feff_inp_file)
        pots = Potential.pot_string_from_file(feff_inp_file)

        # site index (Note: in feff it starts from 1)
        # else case is species symbol
        absorbing_atom = parameters["TARGET"] if "RECIPROCAL" in parameters else pots.splitlines()[3].split()[2]
        return cls(header, parameters, absorbing_atom, data)

    @property
    def energies(self):
        """The absolute energies in eV."""
        return self.data[:, 0]

    @property
    def relative_energies(self):
        """Energy with respect to the Fermi level E - E_f."""
        return self.data[:, 1]

    @property
    def wavenumber(self):
        r"""Get the wave number in units of \AA^-1. k=\sqrt(E - E_f) where E is
        the energy and E_f is the Fermi level computed from electron gas theory
        at the average interstitial charge density.
        """
        return self.data[:, 2]

    @property
    def mu(self):
        """The total absorption cross-section."""
        return self.data[:, 3]

    @property
    def mu0(self):
        """The embedded atomic background absorption."""
        return self.data[:, 4]

    @property
    def chi(self):
        """The normalized fine structure."""
        return self.data[:, 5]

    @property
    def e_fermi(self):
        """The Fermi level in eV."""
        return self.energies[0] - self.relative_energies[0]

    @property
    def source(self):
        """Source identification from Header file."""
        return self.header.source

    @property
    def calc(self):
        """Type of Feff calculation, XANES or EXAFS."""
        return "XANES" if "XANES" in self.parameters else "EXAFS"

    @property
    def material_formula(self):
        """Chemical formula of material from feff.inp file."""
        try:
            form = self.header.formula
        except IndexError:
            form = "No formula provided"
        return "".join(map(str, form))

    @property
    def edge(self):
        """Excitation edge."""
        return self.parameters["EDGE"]

    def as_dict(self):
        """Get dict representations of Xmu object."""
        dct = MSONable.as_dict(self)
        dct["data"] = self.data.tolist()
        return dct


class Eels(MSONable):
    """Parse eels.dat file."""

    def __init__(self, data):
        """
        Args:
            data (numpy.ndarray): data from eels.dat file
        """
        self.data = np.array(data)

    @property
    def energies(self):
        """The energies in eV."""
        return self.data[:, 0]

    @property
    def total_spectrum(self):
        """The total eels spectrum."""
        return self.data[:, 1]

    @property
    def atomic_background(self) -> np.ndarray:
        """The atomic background of EELS."""
        return self.data[:, 2]

    @property
    def fine_structure(self) -> np.ndarray:
        """The fine structure of EELS."""
        return self.data[:, 3]

    @classmethod
    def from_file(cls, eels_dat_file: str = "eels.dat") -> Self:
        """Parse eels spectrum.

        Args:
            eels_dat_file (str): filename and path for eels.dat

        Returns:
            Eels
        """
        data = np.loadtxt(eels_dat_file)
        return cls(data)

    def as_dict(self) -> dict:
        """Get dict representations of Xmu object."""
        dct = MSONable.as_dict(self)
        dct["data"] = self.data.tolist()
        return dct
