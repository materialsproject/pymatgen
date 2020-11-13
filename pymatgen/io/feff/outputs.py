# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module defines classes for parsing the FEFF output files.

Currently supports the xmu.dat, ldos.dat output files are for non-spin case.
"""


from collections import defaultdict, OrderedDict
import re

import numpy as np

from monty.io import zopen
from monty.json import MSONable

from pymatgen import Orbital, Spin, Element
from pymatgen.electronic_structure.dos import Dos, CompleteDos
from pymatgen.io.feff import Header, Potential, Tags

__author__ = "Alan Dozier, Kiran Mathew, Chen Zheng"
__credits__ = "Anubhav Jain, Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "1.0.3"
__maintainer__ = "Alan Dozier"
__email__ = "adozier@uky.edu"
__status__ = "Beta"
__date__ = "April 7, 2013"


class LDos(MSONable):
    """
    Parser for ldos files ldos01, ldos02, .....
    """

    def __init__(self, complete_dos, charge_transfer):
        """
        Args:
            complete_dos (CompleteDos): complete dos object
            charge_transfer (dict): computed charge transfer between atoms
                dictionary
        """
        self.complete_dos = complete_dos
        self.charge_transfer = charge_transfer

    @staticmethod
    def from_file(feff_inp_file='feff.inp', ldos_file='ldos'):
        """
        Creates LDos object from raw Feff ldos files by
        by assuming they are numbered consecutively, i.e. ldos01.dat
        ldos02.dat...

        Args:
            feff_inp_file (str): input file of run to obtain structure
            ldos_file (str): output ldos file of run to obtain dos info, etc.
        """
        header_str = Header.header_string_from_file(feff_inp_file)
        header = Header.from_string(header_str)
        structure = header.struct
        nsites = structure.num_sites
        parameters = Tags.from_file(feff_inp_file)

        if "RECIPROCAL" in parameters:
            pot_dict = dict()
            pot_readstart = re.compile('.*iz.*lmaxsc.*xnatph.*xion.*folp.*')
            pot_readend = re.compile('.*ExternalPot.*switch.*')
            pot_inp = re.sub(r'feff.inp', r'pot.inp', feff_inp_file)
            dos_index = 1
            begin = 0

            with zopen(pot_inp, "r") as potfile:
                for line in potfile:
                    if len(pot_readend.findall(line)) > 0:
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

                    if len(pot_readstart.findall(line)) > 0:
                        begin = 1
        else:
            pot_string = Potential.pot_string_from_file(feff_inp_file)
            dicts = Potential.pot_dict_from_string(pot_string)
            pot_dict = dicts[0]

        with zopen(ldos_file + "00.dat", "r") as fobject:
            f = fobject.readlines()
        efermi = float(f[0].split()[4])

        dos_energies = []
        ldos = {}

        for i in range(1, len(pot_dict) + 1):
            if len(str(i)) == 1:
                ldos[i] = np.loadtxt("{}0{}.dat".format(ldos_file, i))
            else:
                ldos[i] = np.loadtxt("{}{}.dat".format(ldos_file, i))

        for i in range(0, len(ldos[1])):
            dos_energies.append(ldos[1][i][0])

        all_pdos = []
        vorb = {"s": Orbital.s, "p": Orbital.py, "d": Orbital.dxy,
                "f": Orbital.f0}
        forb = {"s": 0, "p": 1, "d": 2, "f": 3}

        dlength = len(ldos[1])

        for i in range(nsites):
            pot_index = pot_dict[structure.species[i].symbol]
            all_pdos.append(defaultdict(dict))
            for k, v in vorb.items():
                density = [ldos[pot_index][j][forb[k] + 1]
                           for j in range(dlength)]
                updos = density
                downdos = None
                if downdos:
                    all_pdos[-1][v] = {Spin.up: updos, Spin.down: downdos}
                else:
                    all_pdos[-1][v] = {Spin.up: updos}

        pdos = all_pdos
        vorb2 = {0: Orbital.s, 1: Orbital.py, 2: Orbital.dxy, 3: Orbital.f0}
        pdoss = {structure[i]: {v: pdos[i][v]
                                for v in vorb2.values()}
                 for i in range(len(pdos))}

        forb = {"s": 0, "p": 1, "d": 2, "f": 3}

        tdos = [0] * dlength
        for i in range(nsites):
            pot_index = pot_dict[structure.species[i].symbol]
            for v in forb.values():
                density = [ldos[pot_index][j][v + 1] for j in range(dlength)]
                for j in range(dlength):
                    tdos[j] = tdos[j] + density[j]
        tdos = {Spin.up: tdos}

        dos = Dos(efermi, dos_energies, tdos)
        complete_dos = CompleteDos(structure, dos, pdoss)
        charge_transfer = LDos.charge_transfer_from_file(feff_inp_file,
                                                         ldos_file)
        return LDos(complete_dos, charge_transfer)

    @staticmethod
    def charge_transfer_from_file(feff_inp_file, ldos_file):
        """
        Get charge transfer from file.

        Args:
            feff_inp_file (str): name of feff.inp file for run
            ldos_file (str): ldos filename for run, assume consequetive order,
                i.e., ldos01.dat, ldos02.dat....

        Returns:
            dictionary of dictionaries in order of potential sites
            ({"p": 0.154, "s": 0.078, "d": 0.0, "tot": 0.232}, ...)
        """
        cht = OrderedDict()
        parameters = Tags.from_file(feff_inp_file)

        if 'RECIPROCAL' in parameters:
            dicts = [dict()]
            pot_dict = dict()
            dos_index = 1
            begin = 0
            pot_inp = re.sub(r'feff.inp', r'pot.inp', feff_inp_file)
            pot_readstart = re.compile('.*iz.*lmaxsc.*xnatph.*xion.*folp.*')
            pot_readend = re.compile('.*ExternalPot.*switch.*')
            with zopen(pot_inp, "r") as potfile:
                for line in potfile:
                    if len(pot_readend.findall(line)) > 0:
                        break
                    if begin == 1:
                        z_number = int(line.strip().split()[0])
                        ele_name = Element.from_Z(z_number).name
                        if len(pot_dict) == 0:
                            pot_dict[0] = ele_name
                        elif len(pot_dict) > 0:
                            pot_dict[max(pot_dict.keys()) + 1] = ele_name
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
                            pot_dict[max(pot_dict.keys()) + 1] = ele_name
                    if len(pot_readstart.findall(line)) > 0:
                        begin = 1
        else:
            pot_string = Potential.pot_string_from_file(feff_inp_file)
            dicts = Potential.pot_dict_from_string(pot_string)
            pot_dict = dicts[1]

        for i in range(0, len(dicts[0]) + 1):
            if len(str(i)) == 1:
                with zopen("{}0{}.dat".format(ldos_file, i), "rt") \
                        as fobject:
                    f = fobject.readlines()
                    s = float(f[3].split()[2])
                    p = float(f[4].split()[2])
                    d = float(f[5].split()[2])
                    f1 = float(f[6].split()[2])
                    tot = float(f[1].split()[4])
                    cht[str(i)] = {pot_dict[i]: {'s': s, 'p': p, 'd': d,
                                                 'f': f1,
                                                 'tot': tot}}
            else:
                with zopen(ldos_file + str(i) + ".dat", "rt") as fid:
                    f = fid.readlines()
                    s = float(f[3].split()[2])
                    p = float(f[4].split()[2])
                    d = float(f[5].split()[2])
                    f1 = float(f[6].split()[2])
                    tot = float(f[1].split()[4])
                    cht[str(i)] = {pot_dict[i]: {'s': s, 'p': p, 'd': d,
                                                 'f': f1,
                                                 'tot': tot}}

        return cht

    def charge_transfer_to_string(self):
        """returns shrage transfer as string"""
        ch = self.charge_transfer
        chts = ['\nCharge Transfer\n\nabsorbing atom']
        for i in range(len(ch)):
            for atom, v2 in ch[str(i)].items():
                a = ['\n', atom, '\n', 's   ', str(v2['s']), '\n',
                     'p   ', str(v2['p']), '\n',
                     'd   ', str(v2['d']), '\n',
                     'f   ', str(v2['f']), '\n',
                     'tot ', str(v2['tot']), '\n']
                chts.extend(a)
        return ''.join(chts)


class Xmu(MSONable):
    r"""
    Parser for data in 'xmu.dat' file.
    The file 'xmu.dat' contains XANES, EXAFS or NRIXS data depending on the
    situation; \\mu, \\mu_0, and \\chi = \\chi * \\mu_0/ \\mu_0/(edge+50eV) as
    functions of absolute energy E, relative energy E − E_f and wave number k.

    Default attributes:
        xmu: Photon absorption cross section of absorbing atom in material
        Energies: Energies of data point
        relative_energies: E - E_fermi
        wavenumber: k=\\sqrt(E −E_fermi)
        mu: The total absorption cross-section.
        mu0: The embedded atomic background absorption.
        chi: fine structure.
        Edge: Aborption Edge
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
            absorbing_atom (str/int): absorbing atom symbol or index
            data (numpy.ndarray, Nx6): cross_sections
        """
        self.header = header
        self.parameters = parameters
        self.absorbing_atom = absorbing_atom
        self.data = np.array(data)

    @staticmethod
    def from_file(xmu_dat_file="xmu.dat", feff_inp_file="feff.inp"):
        """
        Get Xmu from file.

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
        if "RECIPROCAL" in parameters:
            absorbing_atom = parameters["TARGET"]
        # species symbol
        else:
            absorbing_atom = pots.splitlines()[3].split()[2]
        return Xmu(header, parameters, absorbing_atom, data)

    @property
    def energies(self):
        """
        Returns the absolute energies in eV.
        """
        return self.data[:, 0]

    @property
    def relative_energies(self):
        """
        Returns energy with respect to the fermi level.
        E - E_f
        """
        return self.data[:, 1]

    @property
    def wavenumber(self):
        r"""
        Returns The wave number in units of \\AA^-1. k=\\sqrt(E −E_f) where E is
        the energy and E_f is the Fermi level computed from electron gas theory
        at the average interstitial charge density.
        """
        return self.data[:, 2]

    @property
    def mu(self):
        """
        Returns the total absorption cross-section.
        """
        return self.data[:, 3]

    @property
    def mu0(self):
        """
        Returns the embedded atomic background absorption.
        """
        return self.data[:, 4]

    @property
    def chi(self):
        """
        Returns the normalized fine structure.
        """
        return self.data[:, 5]

    @property
    def e_fermi(self):
        """
        Returns the fermi level in eV.
        """
        return self.energies[0] - self.relative_energies[0]

    @property
    def source(self):
        """
        Returns source identification from Header file
        """
        return self.header.source

    @property
    def calc(self):
        """
        Returns type of Feff calculation, XANES or EXAFS
        """
        return "XANES" if "XANES" in self.parameters else "EXAFS"

    @property
    def material_formula(self):
        """
        Returns chemical formula of material from feff.inp file
        """
        try:
            form = self.header.formula
        except IndexError:
            form = 'No formula provided'
        return "".join(map(str, form))

    @property
    def edge(self):
        """
        Returns excitation edge.
        """
        return self.parameters["EDGE"]

    def as_dict(self):
        """
        Returns dict representations of Xmu object
        """
        d = MSONable.as_dict(self)
        d["data"] = self.data.tolist()
        return d


class Eels(MSONable):
    """
    Parse'eels.dat' file.
    """

    def __init__(self, data):
        """
        Args:
            data (): Eels data.
        """
        self.data = np.array(data)

    @property
    def energies(self):
        """
        Returns the energies in eV.
        """
        return self.data[:, 0]

    @property
    def total_spectrum(self):
        """
        Returns the total eels spectrum.
        """
        return self.data[:, 1]

    @property
    def atomic_background(self):
        """
        Returns: atomic background.
        """
        return self.data[:, 2]

    @property
    def fine_structure(self):
        """
        Returns: Fine structure of EELS.
        """
        return self.data[:, 3]

    @staticmethod
    def from_file(eels_dat_file="eels.dat"):
        """
        Parse eels spectrum.

        Args:
            eels_dat_file (str): filename and path for eels.dat

        Returns:
             Eels object
        """
        data = np.loadtxt(eels_dat_file)
        return Eels(data)

    def as_dict(self):
        """
        Returns dict representations of Xmu object
        """
        d = MSONable.as_dict(self)
        d["data"] = self.data.tolist()
        return d
