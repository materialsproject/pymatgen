#!/usr/bin/env python

"""
This module provides classes to run and analyze boltztrap on pymatgen band
structure objects. Boltztrap is a software interpolating band structures and
computing materials properties from this band structure using Boltzmann
semi-classical transport theory.

Boltztrap has been developped by Georg Madsen.

http://www.icams.de/content/departments/ams/madsen/boltztrap.html

References are::

    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71
"""

from __future__ import division

__author__ = "Geoffroy Hautier"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "August 23, 2013"


import os
import math
import numpy as np
import tempfile
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.electronic_structure.dos import Dos, Spin
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.util.io_utils import which
from pymatgen.util.decorators import requires
from pymatgen.core.units import Energy, Length
from pymatgen.core.physical_constants import e, ELECTRON_MASS
import subprocess

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


class BoltztrapRunner():
    """
    This class is used to run Boltztrap on a band structure object
    """

    @requires(which('x_trans'),
              "BoltztrapRunner requires the executables 'x_trans' to be in "
              "the path. Please download the Boltztrap at "
              "http://www.icams.de/content/departments/ams/madsen/boltztrap"
              ".html and follow the instructions in the README to compile "
              "Bolztrap accordingly. Then add x_trans to your path")
    def __init__(self, bs, nelec, dos_type="HISTO", energy_grid=0.005,
                 lpfac=10, sym_nb=None):
        """
        Args:
            bs:
                A band structure object
            nelec:
                the number of electrons
            dos_type:
                two options here for the band structure integration: "HISTO"
                (histogram) or "TETRA" using the tetrahedon method. TETRA
                gives typically better results especially for DOSes but takes
                more time
            energy_grid:
                the energy steps used for the integration. in eV
            lpfac:
                the number of interpolation points in the real space. By
                default 10 gives 10 time more points in the real space than
                the number of kpoints given in reciprocal space
            sym_nb:
                ?
        """
        self.lpfac = lpfac
        self._bs = bs
        self._nelec = nelec
        self.dos_type = dos_type
        self.energy_grid = energy_grid
        self.sym_nb = sym_nb
        self.error = []

    def _make_energy_file(self, file_name):
        with open(file_name, 'w') as f:
            f.write("test\n")
            f.write(str(len(self._bs.kpoints))+"\n")
            for i in range(len(self._bs.kpoints)):
                tmp_eigs = []
                for spin in self._bs._bands:
                    for j in range(int(math.floor(self._bs._nb_bands * 0.9))):
                        tmp_eigs.append(Energy(self._bs._bands[spin][j][i] -
                                        self._bs.efermi, "eV").to("Ry"))
                tmp_eigs.sort()
                f.write("%12.8f %12.8f %12.8f %d\n"
                        % (self._bs.kpoints[i].frac_coords[0],
                           self._bs.kpoints[i].frac_coords[1],
                           self._bs.kpoints[i].frac_coords[2], len(tmp_eigs)))
                for j in range(len(tmp_eigs)):
                    f.write("%18.8f\n" % float(tmp_eigs[j]))

    def _make_struc_file(self, file_name):
        sym = SymmetryFinder(self._bs._structure, symprec=0.01)
        with open(file_name, 'w') as f:
            f.write(self._bs._structure.composition.formula+" " +
                    str(sym.get_spacegroup_symbol())+"\n")
            for i in range(3):
                line = ''
                for j in range(3):
                    line += "%12.5f" % (
                        Length(self._bs._structure.lattice.matrix[i][j],
                               "ang").to("bohr"))
                f.write(line+'\n')
            ops = sym.get_symmetry_dataset()['rotations']
            f.write(str(len(ops))+"\n")
            for c in ops:
                f.write('\n'.join([' '.join([str(int(i)) for i in row])
                                   for row in c]))
                f.write('\n')

    def _make_intrans_file(self, file_name,
                           doping=[1e15, 1e16, 1e17, 1e18, 1e19, 1e20]):
        with open(file_name, 'w') as fout:
            fout.write("GENE          # use generic interface\n")
            fout.write("1 0 0 0.0         # iskip (not presently used) idebug setgap shiftgap \n")
            fout.write(
                "0.0 %f 0.1 %6.1f     # Fermilevel (Ry),energygrid,energy span around Fermilevel, number of electrons\n"
                % (Energy(self.energy_grid, "eV").to("Ry"), self._nelec))
            fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
            fout.write("%d                        # lpfac, number of latt-points per k-point\n" % self.lpfac)
            fout.write("BOLTZ                     # run mode (only BOLTZ is supported)\n")
            fout.write(".15                       # (efcut) energy range of chemical potential\n")
            fout.write("800. 100.                  # Tmax, temperature grid\n")
            fout.write("-1.  # energyrange of bands given DOS output sig_xxx and dos_xxx (xxx is band number)\n")
            fout.write(self.dos_type+"\n")
            fout.write("1 0 0 -1\n")
            fout.write(str(2*len(doping))+"\n")
            for d in doping:
                fout.write(str(d)+"\n")
            for d in doping:
                fout.write(str(-d)+"\n")

    def _make_all_files(self, path):
        if self._bs.is_spin_polarized:
            self._make_energy_file(os.path.join(path, "boltztrap.energyso"))
        else:
            self._make_energy_file(os.path.join(path, "boltztrap.energy"))
        self._make_struc_file(os.path.join(path, "boltztrap.struct"))
        self._make_intrans_file(os.path.join(path, "boltztrap.intrans"))

    def run(self, prev_sigma=None, path_dir=None):
        dir_bz_name = "boltztrap"
        path_dir_orig = path_dir
        if path_dir is None:
            temp_dir = tempfile.mkdtemp()
            path_dir_orig = temp_dir
            path_dir = os.path.join(temp_dir, dir_bz_name)
        else:
            path_dir = os.path.join(path_dir_orig, dir_bz_name)
        if not os.path.exists(path_dir):
            os.mkdir(path_dir)
        else:
            for c in os.listdir(path_dir):
                os.remove(path_dir+"/"+c)
        os.chdir(path_dir)

        self._make_all_files(path_dir)
        if self._bs.is_spin_polarized:
            p = subprocess.Popen(["x_trans", "BoltzTraP", "-so"],
                                 stdout=subprocess.PIPE,
                                 stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            p.wait()
        else:
            p = subprocess.Popen(["x_trans", "BoltzTraP"],
                                 stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            p.wait()
        for c in p.communicate():
            if "STOP error in factorization" in c:
                raise BoltztrapError("STOP error in factorization")

        with open(os.path.join(path_dir, dir_bz_name+".outputtrans")) as f:
            warning = False
            for l in f:
                if "WARNING" in l:
                    warning = True
                    break
            if warning:
                print "There was a warning! Increase lpfac to " + \
                      str(self.lpfac * 2)
                self.lpfac *= 2
                self._make_intrans_file(os.path.join(path_dir,
                                                     dir_bz_name + ".intrans"))
                if self.lpfac > 100:
                    raise BoltztrapError("lpfac higher than 100 and still a warning")
                self.run(path_dir_orig)
        #here we check if the doping levels were well computed
        #sometimes boltztrap mess this up because of two small energy grids
        analyzer = BoltztrapAnalyzer.from_files(path_dir)
        doping_ok = True
        for doping in ['n', 'p']:
            for c in analyzer.mu_doping[doping]:
                if len(analyzer.mu_doping[doping][c])!=len(analyzer.doping[doping]):
                    doping_ok = False
                    break
                if doping == 'p' and \
                        sorted(analyzer.mu_doping[doping][c], reverse=True) != analyzer.mu_doping[doping][c]:
                    doping_ok = False
                    break
                if doping == 'n' and sorted(analyzer.mu_doping[doping][c]) != analyzer.mu_doping[doping][c]:
                    doping_ok=False
                    break
        if not doping_ok:
            self.energy_grid /= 10
            print "lowers energy grid to "+str(self.energy_grid)
            if self.energy_grid < 0.00005:
                raise BoltztrapError("energy grid lower than 0.00005 and still no good doping")
            self._make_intrans_file(path_dir + "/" + dir_bz_name + ".intrans")
            self.run(prev_sigma=None, path_dir=path_dir_orig)
        analyzer = BoltztrapAnalyzer.from_files(path_dir)
        #here, we test if a property (eff_mass tensor) converges
        if prev_sigma is None or \
                abs(sum(analyzer.get_eig_average_eff_mass_tensor()['n']) / 3
                        - prev_sigma)\
                / prev_sigma > 0.05:
            if prev_sigma is not None:
                print abs(sum(analyzer.get_eig_average_eff_mass_tensor()['n'])
                          / 3 - prev_sigma) / prev_sigma, \
                    self.lpfac, \
                    analyzer.get_average_eff_mass_tensor(300, 1e18)
            self.lpfac *= 2
            if self.lpfac > 100:
                raise BoltztrapError("lpfac higher than 100 and still not converged")
            self._make_intrans_file(path_dir + "/" + dir_bz_name + ".intrans")
            self.run(
                prev_sigma=sum(analyzer.get_eig_average_eff_mass_tensor()
                               ['n']) / 3, path_dir=path_dir_orig)
        else:
            print "converged", \
                abs(sum(analyzer.get_eig_average_eff_mass_tensor()['n'])/3
                    - prev_sigma) / prev_sigma, \
                self.lpfac, analyzer.get_average_eff_mass_tensor(300, 1e18)
        return path_dir


class BoltztrapError(Exception):
    """
    Exception class for boltztrap.
    Raised when the boltztrap gives an error
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return "BoltztrapError : " + self.msg


class BoltztrapAnalyzer():
    """
    class used to store all the data from a boltztrap run
    """

    def __init__(self, gap, mu_steps, cond, seebeck, kappa, hall, doping,
                 mu_doping, seebeck_doping, cond_doping, kappa_doping,
                 hall_doping, dos, carrier_conc, vol, warning):
        """
        Constructor taking directly all the data generated by Boltztrap. You
        won't probably use it directly but instead use the from_files and
        from_dict methods.

        Args:
            gap:
                The gap after interpolation in eV
            mu_steps:
                The steps of electron chemical potential (or Fermi level) in
                 eV.
            cond:
                The electronic conductivity tensor divided by a constant
                relaxation time (sigma/tau) at different temperature and
                fermi levels.
                The format is {temperature: [array of 3x3 tensors at each
                fermi level in mu_steps]}. The units are 1/(Ohm*m*s).
            seebeck:
                The Seebeck tensor at different temperatures and fermi levels
                The format is {temperature: [array of 3x3 tensors at each
                fermi level in mu_steps]}. The units are V/K
            kappa:
                The electronic thermal conductivity tensor divided by a
                constant relaxation time (kappa/tau) at different temperature
                and fermi levels. The format is {temperature: [array of 3x3
                tensors at each fermi level in mu_steps]}
                The units are W/(m*K*s)
            hall:
                The hall tensor at different temperature and fermi levels
                The format is {temperature: [array of 27 coefficients list at
                each fermi level in mu_steps]}
                The units are m^3/C
            doping:
                The different doping levels that have been given to Boltztrap
                The format is {'p':[],'n':[]} with an array of doping levels
                The units are cm^-3
            mu_doping:
                Gives the electron chemical potential (or Fermi level) for a
                given set of doping.
                Format is {'p':{temperature: [fermi levels],'n':{temperature:
                [fermi levels]}}
                the fermi level array is ordered according to the doping
                levels in doping units for doping are in cm^-3 and for Fermi
                level in eV
            seebeck_doping:
                The Seebeck tensor at different temperatures and doping levels
                The format is {'p': {temperature: [Seebeck tensors]},
                'n':{temperature: [Seebeck tensors]}}
                The [Seebeck tensors] array is ordered according to the
                doping levels in doping units for doping are in cm^-3 and for
                Seebeck in V/K
            cond_doping:
                The electronic conductivity tensor divided by a constant
                relaxation time (sigma/tau) at different temperatures and
                doping levels
                The format is {'p':{temperature: [conductivity tensors]},
                'n':{temperature: [conductivity tensors]}}
                The [conductivity tensors] array is ordered according to the
                doping levels in doping units for doping are in cm^-3 and for
                conductivity in 1/(Ohm*m*s)
            kappa_doping:
                The thermal conductivity tensor divided by a constant
                relaxation time (kappa/tau) at different temperatures and
                doping levels.
                The format is {'p':{temperature: [thermal conductivity
                tensors]},'n':{temperature: [thermal conductivity tensors]}}
                The [thermal conductivity tensors] array is ordered according
                to the doping levels in doping units for doping are in cm^-3
                and for thermal conductivity in W/(m*K*s)
            hall_doping:
                The Hall tensor at different temperatures and doping levels
                The format is {'p':{temperature: [Hall tensors]},
                'n':{temperature: [Hall tensors]}}
                The [Hall tensors] array is ordered according to the doping
                levels in doping and each Hall tensor is represented by a 27
                coefficients list.
                The units are m^3/C
            the concentration of carriers:
                in electron (or hole) per unit cell
            dos:
                The dos computed by Boltztrap
                given as a pymatgen Dos object
            vol:
                the volume of the unit cell in angstrom cube (A^3)
            warning:
                True if Boltztrap spitted out a warning
        """

        self.gap = gap
        self.mu_steps = mu_steps
        self.cond = cond
        self.seebeck = seebeck
        self.kappa = kappa
        self.hall = hall
        self.warning = warning
        self.doping = doping
        self.mu_doping = mu_doping
        self.seebeck_doping = seebeck_doping
        self.cond_doping = cond_doping
        self.kappa_doping = kappa_doping
        self.hall_doping = hall_doping
        self.carrier_conc = carrier_conc
        self.dos = dos
        self.vol = vol

    @staticmethod
    def _make_boltztrap_analyzer_from_data(
            data_full, data_hall, data_dos, temperature_steps, mu_steps,
            efermi, gap, doping, data_doping_full, data_doping_hall, vol,
            warning=False):
        """
        Make a BoltztrapAnalyzer object from raw data typically parse from
        files.
        """
        cond = {t: [] for t in temperature_steps}
        seebeck = {t: [] for t in temperature_steps}
        kappa = {t: [] for t in temperature_steps}
        hall = {t: [] for t in temperature_steps}
        carrier_conc = {t: [] for t in temperature_steps}
        dos_full = {'energy': [], 'density': []}
        warning = warning
        new_doping = {'p': [], 'n': []}
        for d in doping:
            if d > 0:
                new_doping['p'].append(d)
            else:
                new_doping['n'].append(-d)

        mu_doping = {'p': {t: [] for t in temperature_steps},
                     'n': {t: [] for t in temperature_steps}}
        seebeck_doping = {'p': {t: [] for t in temperature_steps},
                          'n': {t: [] for t in temperature_steps}}
        cond_doping = {'p': {t: [] for t in temperature_steps},
                       'n': {t: [] for t in temperature_steps}}
        kappa_doping = {'p': {t: [] for t in temperature_steps},
                        'n': {t: [] for t in temperature_steps}}
        hall_doping = {'p': {t: [] for t in temperature_steps},
                       'n': {t: [] for t in temperature_steps}}
        for d in data_full:
            carrier_conc[d[1]].append(d[2])
            tens_cond = [[d[3], d[4], d[5]],
                         [d[6], d[7], d[8]],
                         [d[9], d[10], d[11]]]
            cond[d[1]].append(tens_cond)
            tens_seebeck = [[d[12], d[13], d[14]],
                            [d[15], d[16], d[17]],
                            [d[18], d[19], d[20]]]
            seebeck[d[1]].append(tens_seebeck)
            tens_kappa = [[d[21], d[22], d[23]],
                          [d[24], d[25], d[26]],
                          [d[27], d[28], d[29]]]
            kappa[d[1]].append(tens_kappa)

        for d in data_hall:
            hall_tens = d[3:]
            hall[d[1]].append(hall_tens)

        for d in data_doping_full:
            tens_cond = [[d[2], d[3], d[4]],
                         [d[5], d[6], d[7]],
                         [d[8], d[9], d[10]]]
            tens_seebeck = [[d[11], d[12], d[13]],
                            [d[14], d[15], d[16]],
                            [d[17], d[18], d[19]]]
            tens_kappa = [[d[20], d[21], d[22]],
                          [d[23], d[24], d[25]],
                          [d[26], d[27], d[28]]]

            if d[1] < 0:
                mu_doping['n'][d[0]].append(Energy(d[-1], "Ry").to("eV"))
                cond_doping['n'][d[0]].append(tens_cond)
                seebeck_doping['n'][d[0]].append(tens_seebeck)
                kappa_doping['n'][d[0]].append(tens_kappa)
            else:
                mu_doping['p'][d[0]].append(Energy(d[-1], "Ry").to("eV"))
                cond_doping['p'][d[0]].append(tens_cond)
                seebeck_doping['p'][d[0]].append(tens_seebeck)
                kappa_doping['p'][d[0]].append(tens_kappa)

        for i in range(len(data_doping_hall)):
            hall_tens = [data_hall[i][j] for j in range(3, len(data_hall[i]))]
            if data_doping_hall[i][1] < 0:
                hall_doping['n'][data_doping_hall[i][0]].append(hall_tens)
            else:
                hall_doping['p'][data_doping_hall[i][0]].append(hall_tens)

        for t in data_dos:
            dos_full['energy'].append(t[0])
            dos_full['density'].append(t[1])

        dos = Dos(efermi, dos_full['energy'], {Spin.up: dos_full['density']})
        return BoltztrapAnalyzer(
            gap, mu_steps, cond, seebeck, kappa, hall, new_doping, mu_doping,
            seebeck_doping, cond_doping, kappa_doping, hall_doping, dos, carrier_conc,
            vol, warning)

    def get_mu_bounds(self, temp=300):
        return min(self.mu_doping['p'][temp]), max(self.mu_doping['n'][temp])

    def get_average_eff_mass_tensor(self, temperature=300, doping=1e18):
        """
        Gives the average effective mass tensor at a given temperature and
        doping level.
        The average effective mass tensor is defined as the integrated
        average of the second derivative
        This effective mass tensor takes into account:
        -non-parabolicity
        -multiple extrema
        -multiple bands

        Args:
            temperature:
                the temperature in K
            doping:
                the doping in cm^-3

        Returns:
            a dictionnary {'p':[[]],'n':[[]]}
            The arrays are 3x3 and represent the effective mass tensor
            The 'p' links to hole effective mass tensor and 'n' to electron
            effective mass tensor.
        """
        index = None
        import math
        results = {'p': [], 'n': []}
        for t in ['n', 'p']:
            for d in range(len(self.doping[t])):
                if math.fabs(self.doping[t][d]-doping) < 0.001:
                    index = d
            results[t] = np.linalg.inv(
                self.cond_doping[t][temperature][index]) * doping \
                * 10 ** 6 * e ** 2 / ELECTRON_MASS
        return results

    def get_eig_average_eff_mass_tensor(self, temperature=300, doping=1e18):
        """
        Gives the eigenvalues of the average effective mass tensor at a given
        temperature and doping level. The average effective mass tensor is
        defined as the integrated average of the second derivative
        This effective mass tensor takes into account:
        -non-parabolicity
        -multiple extrema
        -multiple bands

        Args:
            temperature:
                the temperature in K
            doping:
                the doping in cm^-3

        Returns:
            a dictionnary {'p':[],'n':[]}
            The list contains the sorted three eigenvalues of the symmetric
            tensor. The 'p' links to hole effective mass tensor and 'n' to
            electron effective mass tensor.
        """
        return {'p':
                sorted(np.linalg.eig(self.get_average_eff_mass_tensor(
                    temperature=temperature, doping=doping)['p'])[0]),
                'n':
                sorted(np.linalg.eig(self.get_average_eff_mass_tensor(
                    temperature=temperature, doping=doping)['n'])[0])}

    @staticmethod
    def from_files(path_dir):
        """
        get a BoltztrapAnalyzer object from a set of files

        Args:
            path_dir:
                directory where the boltztrap files are

        Returns:
            a BoltztrapAnalyzer object

        """
        t_steps = set()
        m_steps = set()
        gap = None
        doping = []
        data_doping_full = []
        data_doping_hall = []
        with open(os.path.join(path_dir, "boltztrap.condtens"), 'r') as f:
            data_full = []
            for line in f:
                if not line.startswith("#"):
                    t_steps.add(int(float(line.split()[1])))
                    m_steps.add(float(line.split()[0]))
                    data_full.append([float(c) for c in line.split()])

        with open(os.path.join(path_dir, "boltztrap.halltens"), 'r') as f:
            data_hall = []
            for line in f:
                if not line.startswith("#"):
                    data_hall.append([float(c) for c in line.split()])

        data_dos = []
        with open(os.path.join(path_dir, "boltztrap.transdos"), 'r') as f:
            for line in f:
                if not line.startswith(" #"):
                    data_dos.append([Energy(line.split()[0], "Ry").to("eV"),
                                     2 * Energy(line.split()[1],
                                                "eV").to("Ry")])

        with open(os.path.join(path_dir, "boltztrap.outputtrans"), 'r') as f:
            warning = False
            step = 0
            for line in f:
                if "WARNING" in line:
                    warning = True
                if line.startswith("VBM"):
                    efermi = Energy(line.split()[1], "Ry").to("eV")

                if step == 2:
                    l_tmp = line.split("-")[1:]
                    doping.extend([-float(d) for d in l_tmp])
                    step = 0

                if step == 1:
                    doping.extend([float(d) for d in line.split()])
                    step = 2

                if line.startswith("Doping levels to be output for") or \
                        line.startswith(" Doping levels to be output for"):
                    step = 1

                if line.startswith("Egap:"):
                    gap = float(line.split()[1])
        if len(doping) != 0:
            with open(os.path.join(path_dir, "fort.26"), 'r') as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:
                        data_doping_full.append([float(c)
                                                 for c in line.split()])

            with open(os.path.join(path_dir, "fort.27"), 'r') as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:

                        data_doping_hall.append([float(c) for c in line.split()])

        with open(os.path.join(path_dir, "boltztrap.struct"), 'r') as f:
            tokens = f.readlines()
            vol = Lattice([[Length(float(tokens[i].split()[j]), "bohr").to("ang")
                            for j in range(3)] for i in range(1, 4)]).volume
        return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
            data_full, data_hall, data_dos, sorted([t for t in t_steps]),
            sorted([Energy(m, "Ry").to("eV") for m in m_steps]), efermi, Energy(gap, "Ry").to("eV"),
            doping, data_doping_full, data_doping_hall, vol, warning)


    @property
    def to_dict(self):
        from pymatgen.util.io_utils import clean_json
        results = {'gap': self.gap,
                   'mu_steps': self.mu_steps,
                   'cond': self.cond,
                   'seebeck': self.seebeck,
                   'kappa': self.kappa,
                   'hall': self.hall,
                   'warning': self.warning, 'doping': self.doping,
                   'mu_doping': self.mu_doping,
                   'seebeck_doping': self.seebeck_doping,
                   'cond_doping': self.cond_doping,
                   'kappa_doping': self.kappa_doping,
                   'hall_doping': self.hall_doping,
                   'dos': self.dos.to_dict,
                   'carrier_conc': self.carrier_conc,
                   'vol': self.vol}
        return clean_json(results)

    @staticmethod
    def from_dict(data):
        def _make_float_array(a):
            res = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            for i in range(3):
                for j in range(3):
                    res[i][j] = float(a[i][j])
            return res

        def _make_float_hall(a):
            return [i for i in a[:27]]

        return BoltztrapAnalyzer(
            float(data['gap']), [float(d) for d in data['mu_steps']],
            {int(d): [_make_float_array(v) for v in data['cond'][d]]
             for d in data['cond']},
            {int(d): [_make_float_array(v) for v in data['seebeck'][d]]
             for d in data['seebeck']},
            {int(d): [_make_float_array(v) for v in data['kappa'][d]]
             for d in data['kappa']},
            {int(d): [_make_float_hall(v) for v in data['hall'][d]]
             for d in data['hall']},
            {'p': [float(d) for d in data['doping']['p']],
             'n': [float(d) for d in data['doping']['n']]},
            {'p': {int(d): [float(v) for v in data['mu_doping']['p'][d]]
                   for d in data['mu_doping']['p']},
             'n': {int(d): [float(v) for v in data['mu_doping']['n'][d]]
                   for d in data['mu_doping']['n']}},
            {'p': {int(d): [_make_float_array(v)
                            for v in data['seebeck_doping']['p'][d]]
                   for d in data['seebeck_doping']['p']},
             'n': {int(d): [_make_float_array(v)
                            for v in data['seebeck_doping']['n'][d]]
                   for d in data['seebeck_doping']['n']}},
            {'p': {int(d): [_make_float_array(v)
                            for v in data['cond_doping']['p'][d]]
                   for d in data['cond_doping']['p']},
             'n': {int(d): [_make_float_array(v)
                            for v in data['cond_doping']['n'][d]]
                   for d in data['cond_doping']['n']}},
            {'p': {int(d): [_make_float_array(v)
                            for v in data['kappa_doping']['p'][d]]
                   for d in data['kappa_doping']['p']},
             'n': {int(d): [_make_float_array(v)
                            for v in data['kappa_doping']['n'][d]]
                   for d in data['kappa_doping']['n']}},
            {'p': {int(d): [_make_float_hall(v)
                            for v in data['hall_doping']['p'][d]]
                   for d in data['hall_doping']['p']},
             'n': {int(d): [_make_float_hall(v)
                            for v in data['hall_doping']['n'][d]]
                   for d in data['hall_doping']['n']}},
            Dos.from_dict(data['dos']), data['carrier_conc'],
            data['vol'], str(data['warning']))


class BoltztrapPlotter():
    """
    class containing methods to plot the data from Boltztrap
    """

    def __init__(self, bz):
        """
        Args:
            bz:
                a BoltztrapAnalyzer object
        """
        self._bz = bz

    def _plot_doping(self, temp):
        limit = 2.21e15
        plt.axvline(self._bz.mu_doping['n'][temp][1], linewidth=3.0,
                    linestyle="--")
        plt.text(self._bz.mu_doping['n'][temp][1] + 0.01,
                 limit,
                 "$n$=10$^{" + str(math.log10(self._bz.doping['n'][1])) + "}$",
                 color='b')
        plt.axvline(self._bz.mu_doping['n'][temp][-1], linewidth=3.0,
                    linestyle="--")
        plt.text(self._bz.mu_doping['n'][temp][-1] + 0.01,
                 limit,
                 "$n$=10$^{" + str(math.log10(self._bz.doping['n'][-1]))
                 + "}$", color='b')
        plt.axvline(self._bz.mu_doping['p'][temp][1], linewidth=3.0,
                    linestyle="--")
        plt.text(self._bz.mu_doping['p'][temp][1] + 0.01,
                 limit,
                 "$p$=10$^{" + str(math.log10(self._bz.doping['p'][1])) + "}$",
                 color='b')
        plt.axvline(self._bz.mu_doping['p'][temp][-1], linewidth=3.0,
                    linestyle="--")
        plt.text(self._bz.mu_doping['p'][temp][-1] + 0.01,
                 limit, "$p$=10$^{" +
                        str(math.log10(self._bz.doping['p'][-1])) + "}$",
                 color='b')

    def _plot_BG_limits(self):
        plt.axvline(0.0, color='k', linewidth=3.0)
        plt.axvline(self._bz.gap, color='k', linewidth=3.0)

    def plot_seebeck(self, temp=300):
        plt.plot(self._bz.mu_steps, [np.linalg.eig(c)[0] * 1e6
                                     for c in self._bz.seebeck[temp]],
                 linewidth=3.0)
        self._plot_BG_limits()
        self._plot_doping(temp)
        plt.legend(['S$_1$', 'S$_2$', 'S$_3$'])
        plt.xlim(-0.5, self._bz.gap+0.5)
        plt.ylabel("Seebeck coefficient ($\mu$V/K)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.show()

    def plot_power_factor(self, temp=300):
        print "blah"
        matrix = []
        for i in range(len(self._bz.mu_steps)):
            matrix.append(np.dot(np.dot(self._bz.seebeck[temp][i],
                                        self._bz.seebeck[temp][i]),
                                 self._bz.cond[temp][i]))
        plt.plot(self._bz.mu_steps, [sorted(np.linalg.eig(c)[0] * 1e6 * 0.01)
                                     for c in matrix],
                 linewidth=3.0)
        self._plot_BG_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap + 0.5)
        plt.ylabel("power factor (S$^2$ * $\sigma$/${\\tau}$)", fontsize=30)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.show()

    def plot_conductivity(self, temp=300):
        plt.semilogy(self._bz.mu_steps, [sorted(np.linalg.eig(c)[0] * 0.01)
                                         for c in self._bz.cond[temp]],
                     linewidth=3.0)
        self._plot_BG_limits()
        self._plot_doping(temp)
        plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
        plt.xlim(-0.5, self._bz.gap + 0.5)
        plt.ylim([1e13, 1e20])
        plt.ylabel("conductivity, $\sigma$/${\\tau}$ (1/($\Omega$ m s))",
                   fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.show()

    def plot_dos(self):
        plotter = DosPlotter(sigma=0.05)
        plotter.add_dos("t", self._bz.dos)
        plotter.get_plot().show()

    def plot_carriers(self, temp=300):
        plt.semilogy(self._bz.mu_steps,
                       abs(self._bz.carrier_conc[temp]/(self._bz.vol*1e-24)),
                       linewidth=3.0, color='r')
        self._plot_BG_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap+0.5)
        plt.ylim(1e14,1e22)
        plt.ylabel("carrier concentration (cm-3)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.show()

