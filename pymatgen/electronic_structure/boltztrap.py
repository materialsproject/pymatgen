# coding: utf-8

from __future__ import division, unicode_literals, print_function

"""
This module provides classes to run and analyze boltztrap on pymatgen band
structure objects. Boltztrap is a software interpolating band structures and
computing materials properties from this band structure using Boltzmann
semi-classical transport theory.

Boltztrap has been developped by Georg Madsen.

http://www.icams.de/content/departments/ams/madsen/boltztrap.html

You need the version 1.2.3

References are::

    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71
"""

__author__ = "Geoffroy Hautier, Zachary Gibbs"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "August 23, 2013"

import os
import math
import numpy as np
import tempfile
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.dos import Dos, Spin, CompleteDos
from pymatgen.electronic_structure.core import Orbital
from pymatgen.electronic_structure.plotter import DosPlotter
from monty.os.path import which
from monty.dev import requires
from monty.json import jsanitize
from pymatgen.core.units import Energy, Length
from pymatgen.core.physical_constants import e, ELECTRON_MASS
import subprocess

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


class BoltztrapRunner(object):
    """
    This class is used to run Boltztrap on a band structure object.

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
            type:
                type of boltztrap usage by default BOLTZ to compute transport
                coefficients
                but you can have also "FERMI" to compute fermi surface or
                more correctly to
                get certain bands interpolated
            band_nb:
                indicates a band number. Used for Fermi Surface interpolation
                (type="SURFACE")
            tauref:
                reference relaxation time. Only set to a value different than
                zero if we want to model
                beyond the constant relaxation time.
            tauexp:
                exponent for the energy in the non-constant relaxation time
                approach
            tauen:
                reference energy for the non-constant relaxation time approach
            soc:
                results from spin-orbit coupling (soc) computations give
                typically non-polarized (no spin up or down)
                results but 1 electron occupations. If the band structure
                comes from a soc computation, you should set
                soc to True (default False)
            doping:
                the fixed doping levels you want to compute. Boltztrap provides both transport values
                depending on electron chemical potential (fermi energy) and for a series of fixed
                carrier concentrations. By default, this is set to 1e16, 1e17, 1e18, 1e19, 1e20 and 1e21
            energy_span_around_fermi:
                usually the interpolation is not needed on the entire energy range but on a specific range around
                the fermi level. This energy gives this range in eV. by default it is 1.5 eV. If you want to have a dos
                on the entire energy range, you will need to change this value.
            scissor:
                scissor to the band gap in eV. This applies a scissor operation moving the band edges without changing
                the band shape. This is useful to correct the often underestimated band gap in DFT. Default is 0.0 (no
                scissor)
    """

    @requires(which('x_trans'),
              "BoltztrapRunner requires the executables 'x_trans' to be in "
              "the path. Please download the Boltztrap at "
              "http://www.icams.de/content/departments/ams/madsen/boltztrap"
              ".html and follow the instructions in the README to compile "
              "Bolztrap accordingly. Then add x_trans to your path")
    def __init__(self, bs, nelec, dos_type="HISTO", energy_grid=0.005,
                 lpfac=10, run_type="BOLTZ", band_nb=None, tauref=0, tauexp=0, tauen=0, soc=False, doping=None,
                 energy_span_around_fermi=1.5, scissor=0.0):
        self.lpfac = lpfac
        self._bs = bs
        self._nelec = nelec
        self.dos_type = dos_type
        self.energy_grid = energy_grid
        self.error = []
        self.run_type = run_type
        self.band_nb = band_nb
        self.tauref = tauref
        self.tauexp = tauexp
        self.tauen = tauen
        self.soc = soc
        self.doping = doping or [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]
        self.energy_span_around_fermi = energy_span_around_fermi
        self.scissor = scissor

    def _make_energy_file(self, file_name):
        with open(file_name, 'w') as f:
            f.write("test\n")
            f.write(str(len(self._bs.kpoints)) + "\n")
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
        sym = SpacegroupAnalyzer(self._bs._structure, symprec=0.01)
        with open(file_name, 'w') as f:
            f.write(self._bs._structure.composition.formula + " " +
                    str(sym.get_spacegroup_symbol()) + "\n")
            for i in range(3):
                line = ''
                for j in range(3):
                    line += "%12.5f" % (
                        Length(self._bs._structure.lattice.matrix[i][j],
                               "ang").to("bohr"))
                f.write(line + '\n')
            ops = sym.get_symmetry_dataset()['rotations']
            f.write(str(len(ops)) + "\n")
            for c in ops:
                f.write('\n'.join([' '.join([str(int(i)) for i in row])
                                   for row in c]))
                f.write('\n')

    def _make_def_file(self, def_file_name):
        with open(def_file_name, 'w') as f:
            so = ""
            if self._bs.is_spin_polarized or self.soc:
                so = "so"
            f.write("5, 'boltztrap.intrans',      'old',    'formatted',0\n" +
                    "6,'boltztrap.outputtrans',      'unknown',    "
                    "'formatted',0\n" +
                    "20,'boltztrap.struct',         'old',    'formatted',0\n" +
                    "10,'boltztrap.energy" + so + "',         'old',    "
                                                  "'formatted',0\n" +
                    "48,'boltztrap.engre',         'unknown',    "
                    "'unformatted',0\n" +
                    "49,'boltztrap.transdos',        'unknown',    "
                    "'formatted',0\n" +
                    "50,'boltztrap.sigxx',        'unknown',    'formatted',"
                    "0\n" +
                    "51,'boltztrap.sigxxx',        'unknown',    'formatted',"
                    "0\n" +
                    "21,'boltztrap.trace',           'unknown',    "
                    "'formatted',0\n" +
                    "22,'boltztrap.condtens',           'unknown',    "
                    "'formatted',0\n" +
                    "24,'boltztrap.halltens',           'unknown',    "
                    "'formatted',0\n" +
                    "30,'boltztrap_BZ.cube',           'unknown',    "
                    "'formatted',0\n" +
                    "35,'boltztrap.banddat',           'unknown',    "
                    "'formatted',0\n" +
                    "36,'boltztrap_band.gpl',           'unknown',    "
                    "'formatted',0\n")

    def _make_proj_files(self, file_name, def_file_name):
        for o in Orbital.all_orbitals:
            for site_nb in range(0, len(self._bs._structure.sites)):
                if o in self._bs._projections[Spin.up][0][0]:
                    with open(file_name + "_" + str(site_nb) + "_" + str(o),
                              'w') as f:
                        f.write(self._bs._structure.composition.formula + "\n")
                        f.write(str(len(self._bs.kpoints)) + "\n")
                        for i in range(len(self._bs.kpoints)):
                            tmp_proj = []
                            for spin in self._bs._bands:
                                for j in range(int(
                                        math.floor(self._bs._nb_bands * 0.9))):
                                    tmp_proj.append(
                                        self._bs._projections[spin][j][i][o][
                                            site_nb])
                            # TODO deal with the sorting going on at the
                            # energy level!!!
                            f.write("%12.8f %12.8f %12.8f %d\n"
                                    % (self._bs.kpoints[i].frac_coords[0],
                                       self._bs.kpoints[i].frac_coords[1],
                                       self._bs.kpoints[i].frac_coords[2],
                                       len(tmp_proj)))
                            for j in range(len(tmp_proj)):
                                f.write("%18.8f\n" % float(tmp_proj[j]))
        with open(def_file_name, 'w') as f:
            so = ""
            if self._bs.is_spin_polarized:
                so = "so"
            f.write("5, 'boltztrap.intrans',      'old',    'formatted',0\n" +
                    "6,'boltztrap.outputtrans',      'unknown',    "
                    "'formatted',0\n" +
                    "20,'boltztrap.struct',         'old',    'formatted',0\n" +
                    "10,'boltztrap.energy" + so + "',         'old',    "
                                                  "'formatted',0\n" +
                    "48,'boltztrap.engre',         'unknown',    "
                    "'unformatted',0\n" +
                    "49,'boltztrap.transdos',        'unknown',    "
                    "'formatted',0\n" +
                    "50,'boltztrap.sigxx',        'unknown',    'formatted',"
                    "0\n" +
                    "51,'boltztrap.sigxxx',        'unknown',    'formatted',"
                    "0\n" +
                    "21,'boltztrap.trace',           'unknown',    "
                    "'formatted',0\n" +
                    "22,'boltztrap.condtens',           'unknown',    "
                    "'formatted',0\n" +
                    "24,'boltztrap.halltens',           'unknown',    "
                    "'formatted',0\n" +
                    "30,'boltztrap_BZ.cube',           'unknown',    "
                    "'formatted',0\n" +
                    "35,'boltztrap.banddat',           'unknown',    "
                    "'formatted',0\n" +
                    "36,'boltztrap_band.gpl',           'unknown',    "
                    "'formatted',0\n")
            i = 1000
            for o in Orbital.all_orbitals:
                for site_nb in range(0, len(self._bs._structure.sites)):
                    if o in self._bs._projections[Spin.up][0][0]:
                        f.write(str(i) + ",\'" + file_name + "_" + str(
                            site_nb) + "_" + str(o)
                                + "\' \'old\', \'formatted\',0\n")
                        i += 1

    def _make_intrans_file(self, file_name):
        if self.run_type == "BOLTZ":
            with open(file_name, 'w') as fout:
                fout.write("GENE          # use generic interface\n")
                setgap = 0
                if self.scissor > 0.0001:
                    setgap = 1
                fout.write("1 0 %d %f         # iskip (not presently used) idebug setgap shiftgap \n"
                           % (setgap, Energy(self.scissor,"eV").to("Ry")))
                fout.write(
                    "0.0 %f %f %6.1f     # Fermilevel (Ry),energygrid,energy span around Fermilevel, "
                    "number of electrons\n"
                    % (Energy(self.energy_grid, "eV").to("Ry"), Energy(self.energy_span_around_fermi, "eV").to("Ry"),
                       self._nelec))
                fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
                fout.write("%d                        # lpfac, number of latt-points per k-point\n" % self.lpfac)
                fout.write("BOLTZ                     # run mode (only BOLTZ is supported)\n")
                fout.write(".15                       # (efcut) energy range of chemical potential\n")
                fout.write("1300. 100.                  # Tmax, temperature grid\n")
                fout.write("-1.  # energyrange of bands given DOS output sig_xxx and dos_xxx (xxx is band number)\n")
                fout.write(self.dos_type+"\n")
                fout.write(str(self.tauref)+" "+str(self.tauexp)+" "+str(self.tauen)+" 0 0 0\n")
                fout.write(str(2*len(self.doping))+"\n")
                for d in self.doping:
                    fout.write(str(d) + "\n")
                for d in self.doping:
                    fout.write(str(-d) + "\n")
        elif self.run_type == "FERMI":
            with open(file_name, 'w') as fout:
                fout.write("GENE          # use generic interface\n")
                fout.write(
                    "1 0 0 0.0         # iskip (not presently used) idebug "
                    "setgap shiftgap \n")
                fout.write(
                    "0.0 %f 0.1 %6.1f     # Fermilevel (Ry),energygrid,"
                    "energy span around Fermilevel, "
                    "number of electrons\n"
                    % (Energy(self.energy_grid, "eV").to("Ry"), self._nelec))
                fout.write(
                    "CALC                    # CALC (calculate expansion "
                    "coeff), NOCALC read from file\n")
                fout.write(
                    "%d                        # lpfac, number of latt-points "
                    "per k-point\n" % self.lpfac)
                fout.write(
                    "FERMI                     # run mode (only BOLTZ is "
                    "supported)\n")
                fout.write(str(self.band_nb + 1))

    def _make_all_files(self, path):
        if self._bs.is_spin_polarized or self.soc:
            self._make_energy_file(os.path.join(path, "boltztrap.energyso"))
        else:
            self._make_energy_file(os.path.join(path, "boltztrap.energy"))
        self._make_struc_file(os.path.join(path, "boltztrap.struct"))
        self._make_intrans_file(os.path.join(path, "boltztrap.intrans"))
        self._make_def_file("BoltzTraP.def")
        if len(self._bs._projections) != 0:
            self._make_proj_files(os.path.join(path, "boltztrap.proj"),
                                  os.path.join(path, "BoltzTraP.def"))

    def run(self, prev_sigma=None, path_dir=None, convergence=True):
        if self.run_type == "FERMI":
            convergence = False
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
                os.remove(path_dir + "/" + c)
        os.chdir(path_dir)

        self._make_all_files(path_dir)
        if self._bs.is_spin_polarized or self.soc:
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

        with open(os.path.join(path_dir, dir_bz_name + ".outputtrans")) as f:
            warning = False
            for l in f:
                if "WARNING" in l:
                    warning = True
                    break
            if warning:
                print("There was a warning! Increase lpfac to " +
                      str(self.lpfac * 2))
                self.lpfac *= 2
                self._make_intrans_file(os.path.join(path_dir,
                                                     dir_bz_name + ".intrans"))
                if self.lpfac > 100:
                    raise BoltztrapError(
                        "lpfac higher than 100 and still a warning")
                self.run(path_dir_orig)
        # here we check if the doping levels were well computed
        # sometimes boltztrap mess this up because of two small energy grids
        analyzer = BoltztrapAnalyzer.from_files(path_dir)
        doping_ok = True
        print(analyzer.mu_doping, analyzer.doping)
        for doping in ['n', 'p']:
            for c in analyzer.mu_doping[doping]:
                if len(analyzer.mu_doping[doping][c]) != len(
                        analyzer.doping[doping]):
                    doping_ok = False
                    break
                if doping == 'p' and \
                                sorted(analyzer.mu_doping[doping][c],
                                       reverse=True) != \
                                analyzer.mu_doping[doping][c]:
                    doping_ok = False
                    break
                if doping == 'n' and sorted(analyzer.mu_doping[doping][c]) != \
                        analyzer.mu_doping[doping][c]:
                    doping_ok = False
                    break
        if not doping_ok:
            self.energy_grid /= 10
            print("lowers energy grid to " + str(self.energy_grid))
            if self.energy_grid < 0.00005:
                raise BoltztrapError(
                    "energy grid lower than 0.00005 and still no good doping")
            self._make_intrans_file(path_dir + "/" + dir_bz_name + ".intrans")
            self.run(prev_sigma=None, path_dir=path_dir_orig)
        analyzer = BoltztrapAnalyzer.from_files(path_dir)
        # here, we test if a property (eff_mass tensor) converges
        if convergence is False:
            return path_dir
        if prev_sigma is None or abs(sum(
                analyzer.get_average_eff_mass()['n'][300][
                    int(len(self.doping) / 2)]) / 3
                                             - prev_sigma) / prev_sigma > 0.05:
            if prev_sigma is not None:
                print((abs(sum(analyzer.get_average_eff_mass()['n'][300][
                                   int(len(self.doping) / 2)]) / 3
                           - prev_sigma) / prev_sigma, self.lpfac))
            self.lpfac *= 2
            if self.lpfac > 100:
                raise BoltztrapError(
                    "lpfac higher than 100 and still not converged")
            self._make_intrans_file(path_dir + "/" + dir_bz_name + ".intrans")
            self.run(
                prev_sigma=sum(analyzer.get_average_eff_mass()
                               ['n'][300][int(len(self.doping) / 2)]) / 3,
                path_dir=path_dir_orig)
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


class BoltztrapAnalyzer(object):
    """
    Class used to store all the data from a boltztrap run
    """

    def __init__(self, gap, mu_steps, cond, seebeck, kappa, hall, doping,
                 mu_doping, seebeck_doping, cond_doping, kappa_doping,
                 hall_doping, dos, dos_partial, carrier_conc, vol, warning):
        """
        Constructor taking directly all the data generated by Boltztrap. You
        won't probably use it directly but instead use the from_files and
        from_dict methods.

        Args:
            gap: The gap after interpolation in eV
            mu_steps: The steps of electron chemical potential (or Fermi
                level) in eV.
            cond: The electronic conductivity tensor divided by a constant
                relaxation time (sigma/tau) at different temperature and
                fermi levels.
                The format is {temperature: [array of 3x3 tensors at each
                fermi level in mu_steps]}. The units are 1/(Ohm*m*s).
            seebeck: The Seebeck tensor at different temperatures and fermi
                levels. The format is {temperature: [array of 3x3 tensors at
                each fermi level in mu_steps]}. The units are V/K
            kappa: The electronic thermal conductivity tensor divided by a
                constant relaxation time (kappa/tau) at different temperature
                and fermi levels. The format is {temperature: [array of 3x3
                tensors at each fermi level in mu_steps]}
                The units are W/(m*K*s)
            hall: The hall tensor at different temperature and fermi levels
                The format is {temperature: [array of 27 coefficients list at
                each fermi level in mu_steps]}
                The units are m^3/C
            doping: The different doping levels that have been given to
                Boltztrap. The format is {'p':[],'n':[]} with an array of
                doping levels. The units are cm^-3
            mu_doping: Gives the electron chemical potential (or Fermi level)
                for a given set of doping.
                Format is {'p':{temperature: [fermi levels],'n':{temperature:
                [fermi levels]}}
                the fermi level array is ordered according to the doping
                levels in doping units for doping are in cm^-3 and for Fermi
                level in eV
            seebeck_doping: The Seebeck tensor at different temperatures and
                doping levels. The format is {'p': {temperature: [Seebeck
                tensors]}, 'n':{temperature: [Seebeck tensors]}}
                The [Seebeck tensors] array is ordered according to the
                doping levels in doping units for doping are in cm^-3 and for
                Seebeck in V/K
            cond_doping: The electronic conductivity tensor divided by a
                constant relaxation time (sigma/tau) at different
                temperatures and doping levels
                The format is {'p':{temperature: [conductivity tensors]},
                'n':{temperature: [conductivity tensors]}}
                The [conductivity tensors] array is ordered according to the
                doping levels in doping units for doping are in cm^-3 and for
                conductivity in 1/(Ohm*m*s)
            kappa_doping: The thermal conductivity tensor divided by a constant
                relaxation time (kappa/tau) at different temperatures and
                doping levels.
                The format is {'p':{temperature: [thermal conductivity
                tensors]},'n':{temperature: [thermal conductivity tensors]}}
                The [thermal conductivity tensors] array is ordered according
                to the doping levels in doping units for doping are in cm^-3
                and for thermal conductivity in W/(m*K*s)
            hall_doping: The Hall tensor at different temperatures and doping
                levels.
                The format is {'p':{temperature: [Hall tensors]},
                'n':{temperature: [Hall tensors]}}
                The [Hall tensors] array is ordered according to the doping
                levels in doping and each Hall tensor is represented by a 27
                coefficients list.
                The units are m^3/C
            carrier_conc: The concentration of carriers in electron (or hole)
                per unit cell
            dos: The dos computed by Boltztrap given as a pymatgen Dos object
            dos_partial: Data for the partial DOS projected on sites and
                orbitals
            vol: Volume of the unit cell in angstrom cube (A^3)
            warning: True if Boltztrap spitted out a warning
        """
        self.gap = gap
        self.mu_steps = mu_steps
        self._cond = cond
        self._seebeck = seebeck
        self._kappa = kappa
        self._hall = hall
        self.warning = warning
        self.doping = doping
        self.mu_doping = mu_doping
        self._seebeck_doping = seebeck_doping
        self._cond_doping = cond_doping
        self._kappa_doping = kappa_doping
        self._hall_doping = hall_doping
        self._carrier_conc = carrier_conc
        self.dos = dos
        self.vol = vol
        self._dos_partial = dos_partial

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
            hall_tens = [[[d[3], d[4], d[5]],
                          [d[6], d[7], d[8]],
                          [d[9], d[10], d[11]]],
                         [[d[12], d[13], d[14]],
                          [d[15], d[16], d[17]],
                          [d[18], d[19], d[20]]],
                         [[d[21], d[22], d[23]],
                          [d[24], d[25], d[26]],
                          [d[27], d[28], d[29]]]]
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

        for d in data_doping_hall:
            hall_tens = [[[d[2], d[3], d[4]],
                          [d[5], d[6], d[7]],
                          [d[8], d[9], d[10]]],
                         [[d[11], d[12], d[13]],
                          [d[14], d[15], d[16]],
                          [d[17], d[18], d[19]]],
                         [[d[20], d[21], d[22]],
                          [d[23], d[24], d[25]],
                          [d[26], d[27], d[28]]]]
            if d[1] < 0:
                hall_doping['n'][d[0]].append(hall_tens)
            else:
                hall_doping['p'][d[0]].append(hall_tens)

        for t in data_dos['total']:
            dos_full['energy'].append(t[0])
            dos_full['density'].append(t[1])

        dos = Dos(efermi, dos_full['energy'], {Spin.up: dos_full['density']})
        dos_partial = data_dos['partial']

        return BoltztrapAnalyzer(
            gap, mu_steps, cond, seebeck, kappa, hall, new_doping, mu_doping,
            seebeck_doping, cond_doping, kappa_doping, hall_doping, dos,
            dos_partial, carrier_conc,
            vol, warning)

    def get_seebeck(self, output='eig', doping_levels=True):
        """
            Gives the seebeck coefficient in either a full 3x3 tensor form,
            as 3 eigenvalues, or as the average value
            (trace/3.0) If doping_levels=True, the results are given at
            different p and n doping
            levels (given by self.doping), otherwise it is given as a series
            of electron chemical potential values

            Args:
                output (string): the type of output. 'tensor' give the full
                3x3 tensor, 'eig' its 3 eigenvalues and
                'average' the average of the three eigenvalues
                doping_levels (boolean): True for the results to be given at
                different doping levels, False for results
                at different electron chemical potentials

            Returns:
                If doping_levels=True, a dictionnary {temp:{'p':[],'n':[]}}.
                The 'p' links to Seebeck at p-type doping
                and 'n' to the Seebeck at n-type doping. Otherwise, returns a
                {temp:[]} dictionary
                The result contains either the sorted three eigenvalues of
                the symmetric
                Seebeck tensor (output='eig') or a full tensor (3x3 array) (
                output='tensor') or as an average
                (output='average').

                units are microV/K
        """
        return BoltztrapAnalyzer._format_to_output(self._seebeck,
                                                   self._seebeck_doping, output,
                                                   doping_levels, 1e6)

    def get_conductivity(self, output='eig', doping_levels=True,
                         relaxation_time=1e-14):
        """
            Gives the conductivity in either a full 3x3 tensor form,
            as 3 eigenvalues, or as the average value
            (trace/3.0) If doping_levels=True, the results are given at
            different p and n doping
            levels (given by self.doping), otherwise it is given as a series
            of electron chemical potential values

            Args:
                output (string): the type of output. 'tensor' give the full
                3x3 tensor, 'eig' its 3 eigenvalues and
                'average' the average of the three eigenvalues
                doping_levels (boolean): True for the results to be given at
                different doping levels, False for results
                at different electron chemical potentials

            Returns:
                If doping_levels=True, a dictionnary {temp:{'p':[],'n':[]}}.
                The 'p' links to conductivity
                at p-type doping and 'n' to the conductivity at n-type
                doping. Otherwise,
                returns a {temp:[]} dictionary. The result contains either
                the sorted three eigenvalues of the symmetric
                conductivity tensor (format='eig') or a full tensor (3x3
                array) (output='tensor') or as an average
                (output='average').
                The result includes a given constant relaxation time

                units are 1/Ohm*m
        """
        return BoltztrapAnalyzer._format_to_output(self._cond,
                                                   self._cond_doping, output,
                                                   doping_levels,
                                                   relaxation_time)

    def get_power_factor(self, output='eig', doping_levels=True,
                         relaxation_time=1e-14):
        """
        Gives the power factor (Seebeck^2 * conductivity) in either a full
        3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values

        Args:
            output (string): the type of output. 'tensor' give the full 3x3
            tensor, 'eig' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (boolean): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials

        Returns:
            If doping_levels=True, a dictionnary {temp:{'p':[],'n':[]}}. The
            'p' links to power factor
            at p-type doping and 'n' to the conductivity at n-type doping.
            Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            power factor tensor (format='eig') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time

            units are microW/(m K^2)
        """
        result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]}
                         for doping in self._seebeck_doping}
        for doping in result_doping:
            for temp in result_doping[doping]:
                for i in range(len(self.doping[doping])):
                    full_tensor = np.dot(self._cond_doping[doping][temp][i],
                                         np.dot(
                                             self._seebeck_doping[doping][temp][
                                                 i],
                                             self._seebeck_doping[doping][temp][
                                                 i]))
                    result_doping[doping][temp].append(full_tensor)

        result = {temp: [] for temp in self._seebeck}
        for temp in result:
            for i in range(len(self.mu_steps)):
                full_tensor = np.dot(self._cond[temp][i],
                                     np.dot(self._seebeck[temp][i],
                                            self._seebeck[temp][i]))
                result[temp].append(full_tensor)
        return BoltztrapAnalyzer._format_to_output(result, result_doping,
                                                   output, doping_levels,
                                                   multi=1e6 * relaxation_time)

    def get_thermal_conductivity(self, output='eig', doping_levels=True,
                                 relaxation_time=1e-14):
        """
        Gives the electronic part of the thermal conductivity in either a
        full 3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values

        Args:
            output (string): the type of output. 'tensor' give the full 3x3
            tensor, 'eig' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (boolean): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}. The
            'p' links to thermal conductivity
            at p-type doping and 'n' to the thermal conductivity at n-type
            doping. Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            conductivity tensor (format='eig') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time

            units are W/mK
        """
        result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]}
                         for doping in self._seebeck_doping}
        for doping in result_doping:
            for temp in result_doping[doping]:
                for i in range(len(self.doping[doping])):
                    pf_tensor = np.dot(self._cond_doping[doping][temp][i],
                                       np.dot(
                                           self._seebeck_doping[doping][temp][
                                               i],
                                           self._seebeck_doping[doping][temp][
                                               i]))
                    result_doping[doping][temp].append((self._kappa_doping[
                                                            doping][temp][
                                                            i] - pf_tensor *
                                                        temp))

        result = {temp: [] for temp in self._seebeck}
        for temp in result:
            for i in range(len(self.mu_steps)):
                pf_tensor = np.dot(self._cond[temp][i],
                                   np.dot(self._seebeck[temp][i],
                                          self._seebeck[temp][i]))
                result[temp].append((self._kappa[temp][i] - pf_tensor * temp))

        return BoltztrapAnalyzer._format_to_output(result, result_doping,
                                                   output, doping_levels,
                                                   multi=relaxation_time)

    def get_zt(self, output='eig', doping_levels=True, relaxation_time=1e-14,
               kl=0.2):
        """
        Gives the ZT coefficient (S^2*cond*T/thermal cond) in either a full
        3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values. We assume a constant relaxation
        time and a constant
        lattice thermal conductivity

        Args:
            output (string): the type of output. 'tensor' give the full 3x3
            tensor, 'eig' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (boolean): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}. The
            'p' links to ZT
            at p-type doping and 'n' to the ZT at n-type doping. Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            ZT tensor (format='eig') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time and lattice
            thermal conductivity
        """
        result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]}
                         for doping in self._seebeck_doping}
        for doping in result_doping:
            for temp in result_doping[doping]:
                for i in range(len(self.doping[doping])):
                    pf_tensor = np.dot(self._cond_doping[doping][temp][i],
                                       np.dot(
                                           self._seebeck_doping[doping][temp][
                                               i],
                                           self._seebeck_doping[doping][temp][
                                               i]))
                    thermal_conduct = (self._kappa_doping[doping][temp][
                                           i] - pf_tensor * temp) * \
                                      relaxation_time
                    result_doping[doping][temp].append(
                        np.dot(pf_tensor * relaxation_time * temp,
                               np.linalg.inv(
                                   thermal_conduct + kl * np.eye(3, 3))))

        result = {temp: [] for temp in self._seebeck}
        for temp in result:
            for i in range(len(self.mu_steps)):
                pf_tensor = np.dot(self._cond[temp][i],
                                   np.dot(self._seebeck[temp][i],
                                          self._seebeck[temp][i]))
                thermal_conduct = (self._kappa[temp][
                                       i] - pf_tensor * temp) * relaxation_time
                result[temp].append(np.dot(pf_tensor * relaxation_time * temp,
                                           np.linalg.inv(
                                               thermal_conduct + kl * np.eye(3,
                                                                             3))))
        return BoltztrapAnalyzer._format_to_output(result, result_doping,
                                                   output, doping_levels)

    def get_average_eff_mass(self, output='eig'):
        """
        Gives the average effective mass tensor. We call it average because
        it takes into account all the bands
        and regons in the Brillouin zone. This is different than the standard
        textbook effective mass which relates
        often to only one (parabolic) band.
        The average effective mass tensor is defined as the integrated
        average of the second derivative of E(k)
        This effective mass tensor takes into account:
        -non-parabolicity
        -multiple extrema
        -multiple bands

        For more information about it. See:

        Hautier, G., Miglio, A., Waroquiers, D., Rignanese, G., & Gonze,
        X. (2014).
        How Does Chemistry Influence Electron Effective Mass in Oxides?
        A High-Throughput Computational Analysis. Chemistry of Materials,
        26(19), 5447–5458. doi:10.1021/cm404079a

        or

        Hautier, G., Miglio, A., Ceder, G., Rignanese, G.-M., & Gonze,
        X. (2013).
        Identification and design principles of low hole effective mass
        p-type transparent conducting oxides.
        Nature Communications, 4, 2292. doi:10.1038/ncomms3292

        Depending on the value of output, we have either the full 3x3
        effective mass tensor,
        its 3 eigenvalues or an average

        Args:
            output (string): 'eigs' for eigenvalues, 'tensor' for the full
            tensor and 'average' for an average (trace/3)

        Returns:
            a dictionary {'p':{temp:[]},'n':{temp:[]}} with an array of
            effective mass tensor, eigenvalues of average
            value (depending on output) for each temperature and for each
            doping level.
            The 'p' links to hole effective mass tensor and 'n' to electron
            effective mass tensor.
        """

        result_doping = {doping: {t: [] for t in self._cond_doping[doping]} for
                         doping in self.doping}
        for doping in result_doping:
            for temp in result_doping[doping]:
                for i in range(len(self.doping[doping])):
                    if output == 'tensor':
                        result_doping[doping][temp].append(np.linalg.inv(
                            np.array(self._cond_doping[doping][temp][i])) \
                                                           *
                                                           self.doping[doping][
                                                               i] * 10 ** 6 *
                                                           e ** 2 /
                                                           ELECTRON_MASS)
                    elif output == 'eig':
                        result_doping[doping][temp].append(
                            sorted(np.linalg.eigh(np.linalg.inv(
                                np.array(self._cond_doping[doping][temp][i])) *
                                                  self.doping[doping][
                                                      i] * 10 ** 6 * e ** 2 \
                                                  / ELECTRON_MASS)[0]))
                    else:
                        full_tensor = np.linalg.inv(
                            np.array(self._cond_doping[doping][temp][i])) \
                                      * self.doping[doping][
                                          i] * 10 ** 6 * e ** 2 / ELECTRON_MASS
                        result_doping[doping][temp].append((full_tensor[0][0] \
                                                            + full_tensor[1][
                                                                1] \
                                                            + full_tensor[2][
                                                                2]) / 3.0)
        return result_doping

    @staticmethod
    def _format_to_output(tensor, tensor_doping, output, doping_levels,
                          multi=1.0):
        if doping_levels:
            full_tensor = tensor_doping
            result = {doping: {t: [] for t in tensor_doping[doping]} for doping
                      in tensor_doping}
            for doping in full_tensor:
                for temp in full_tensor[doping]:
                    for i in range(len(full_tensor[doping][temp])):
                        if output == 'eig':
                            result[doping][temp].append(sorted(
                                np.linalg.eigh(full_tensor[doping][temp][i])[
                                    0] * multi))
                        elif output == 'tensor':
                            result[doping][temp].append(
                                np.array(full_tensor[doping][temp][i]) * multi)
                        else:
                            result[doping][temp].append(
                                (full_tensor[doping][temp][i][0][0] \
                                 + full_tensor[doping][temp][i][1][1] \
                                 + full_tensor[doping][temp][i][2][
                                     2]) * multi / 3.0)
        else:
            full_tensor = tensor
            result = {t: [] for t in tensor}
            for temp in full_tensor:
                for i in range(len(tensor[temp])):
                    if output == 'eig':
                        result[temp].append(sorted(
                            np.linalg.eigh(full_tensor[temp][i])[0] * multi))
                    elif output == 'tensor':
                        result[temp].append(
                            np.array(full_tensor[temp][i]) * multi)
                    else:
                        result[temp].append((full_tensor[temp][i][0][0] \
                                             + full_tensor[temp][i][1][1] \
                                             + full_tensor[temp][i][2][
                                                 2]) * multi / 3.0)
        return result

    def get_complete_dos(self, structure):
        """
        Gives a CompleteDos object with the DOS from the interpolated
        projected band structure
        Args:
            the structure (necessary to identify sites for projection)

        Returns:
            a CompleteDos object
        """
        pdoss = {}
        for s in self._dos_partial:
            if structure.sites[int(s)] not in pdoss:
                pdoss[structure.sites[int(s)]] = {}
            for o in self._dos_partial[s]:
                if Orbital.from_string(o) not in pdoss[structure.sites[int(s)]]:
                    pdoss[structure.sites[int(s)]][Orbital.from_string(o)] = {}
                pdoss[structure.sites[int(s)]][Orbital.from_string(o)][
                    Spin.up] = self._dos_partial[s][o]
        return CompleteDos(structure, total_dos=self.dos, pdoss=pdoss)

    def get_mu_bounds(self, temp=300):
        return min(self.mu_doping['p'][temp]), max(self.mu_doping['n'][temp])

    def get_carrier_concentration(self):
        """
        gives the carrier concentration (in cm^-3)

        Returns
            a dictionary {temp:[]} with an array of carrier concentration (in cm^-3) at each temperature
            The array relates to each step of electron chemical potential
        """

        return {temp: [1e24*i/self.vol for i in self._carrier_conc[temp]] for temp in self._carrier_conc}

    def get_hall_carrier_concentration(self):
        """
        gives the Hall carrier concentration (in cm^-3). This is the trace of the Hall tensor (see Boltztrap source
        code) Hall carrier concentration are not always exactly the same than carrier concentration.

        Returns
            a dictionary {temp:[]} with an array of Hall carrier concentration (in cm^-3) at each temperature
            The array relates to each step of electron chemical potential
        """
        result = {temp: [] for temp in self._hall}
        for temp in self._hall:
            for i in self._hall[temp]:
                trace = (i[1][2][0]+i[2][0][1]+i[0][1][2])/3.0
                if trace != 0.0:
                    result[temp].append(1e-6/(trace*e))
                else:
                    result[temp].append(0.0)
        return result

    @staticmethod
    def from_files(path_dir):
        """
        get a BoltztrapAnalyzer object from a set of files

        Args:
            path_dir: directory where the boltztrap files are

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

        data_dos = {'total': [], 'partial': {}}
        with open(os.path.join(path_dir, "boltztrap.transdos"), 'r') as f:
            count_series = 0
            for line in f:
                if not line.startswith(" #"):
                    data_dos['total'].append(
                        [Energy(float(line.split()[0]), "Ry").to("eV"),
                         float(line.split()[1])])
                    total_elec = float(line.split()[2])
                else:
                    count_series += 1
                if count_series > 1:
                    break
        data_dos['total'] = [
            [data_dos['total'][i][0], 2 * data_dos['total'][i][1] / total_elec]
            for i in range(len(data_dos['total']))]
        for file_name in os.listdir(path_dir):
            if file_name.endswith(
                    "transdos") and file_name != 'boltztrap.transdos':
                tokens = file_name.split(".")[1].split("_")
                with open(os.path.join(path_dir, file_name), 'r') as f:
                    for line in f:
                        if not line.startswith(" #"):
                            if tokens[1] not in data_dos['partial']:
                                data_dos['partial'][tokens[1]] = {}
                            if tokens[2] not in data_dos['partial'][tokens[1]]:
                                data_dos['partial'][tokens[1]][tokens[2]] = []
                            data_dos['partial'][tokens[1]][tokens[2]].append(
                                2 * float(line.split()[1]))

        with open(os.path.join(path_dir, "boltztrap.outputtrans"), 'r') as f:
            warning = False
            step = 0
            for line in f:
                if "WARNING" in line:
                    warning = True
                if line.startswith("VBM"):
                    efermi = Energy(line.split()[1], "Ry").to("eV")
                if line.startswith("Doping level number"):
                    doping.append(float(line.split()[6]))
                if line.startswith("Egap:"):
                    gap = float(line.split()[1])
        if len(doping) != 0:
            with open(os.path.join(path_dir, "boltztrap.condtens_fixdoping"),
                      'r') as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:
                        data_doping_full.append([float(c)
                                                 for c in line.split()])

            with open(os.path.join(path_dir, "boltztrap.halltens_fixdoping"),
                      'r') as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:
                        data_doping_hall.append(
                            [float(c) for c in line.split()])

        with open(os.path.join(path_dir, "boltztrap.struct"), 'r') as f:
            tokens = f.readlines()
            vol = Lattice([[Length(float(tokens[i].split()[j]), "bohr").to("ang")
                            for j in range(3)] for i in range(1, 4)]).volume

        data_dos = {'total': [], 'partial': {}}
        with open(os.path.join(path_dir, "boltztrap.transdos"), 'r') as f:
            count_series = 0
            normalization_factor = None
            for line in f:
                if not line.startswith(" #"):
                    data_dos['total'].append([Energy(float(line.split()[0]), "Ry").to("eV"),
                                              float(line.split()[1])])
                else:
                    count_series += 1
                if count_series > 1:
                    break
        for file_name in os.listdir(path_dir):
            if file_name.endswith("transdos") and file_name != 'boltztrap.transdos':
                tokens = file_name.split(".")[1].split("_")
                with open(os.path.join(path_dir, file_name), 'r') as f:
                    for line in f:
                        if not line.startswith(" #"):
                            if tokens[1] not in data_dos['partial']:
                                data_dos['partial'][tokens[1]] = {}
                            if tokens[2] not in data_dos['partial'][tokens[1]]:
                                data_dos['partial'][tokens[1]][tokens[2]] = []
                            data_dos['partial'][tokens[1]][tokens[2]].append(float(line.split()[1]))

        return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
            data_full, data_hall, data_dos, sorted([t for t in t_steps]),
            sorted([Energy(m, "Ry").to("eV") for m in m_steps]), efermi,
            Energy(gap, "Ry").to("eV"),
            doping, data_doping_full, data_doping_hall, vol, warning)

    def as_dict(self):

        results = {'gap': self.gap,
                   'mu_steps': self.mu_steps,
                   'cond': self._cond,
                   'seebeck': self._seebeck,
                   'kappa': self._kappa,
                   'hall': self._hall,
                   'warning': self.warning, 'doping': self.doping,
                   'mu_doping': self.mu_doping,
                   'seebeck_doping': self._seebeck_doping,
                   'cond_doping': self._cond_doping,
                   'kappa_doping': self._kappa_doping,
                   'hall_doping': self._hall_doping,
                   'dos': self.dos.as_dict(),
                   'dos_partial': self._dos_partial,
                   'carrier_conc': self._carrier_conc,
                   'vol': self.vol}
        return jsanitize(results)

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
            Dos.from_dict(data['dos']), data['dos_partial'],
            data['carrier_conc'],
            data['vol'], str(data['warning']))


class BoltztrapPlotter(object):
    """
    class containing methods to plot the data from Boltztrap.

    Args:
        bz: a BoltztrapAnalyzer object
    """

    def __init__(self, bz):
        self._bz = bz

    def _plot_doping(self, temp):
        if len(self._bz.doping) != 0:
            limit = 2.21e15
            plt.axvline(self._bz.mu_doping['n'][temp][0], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['n'][temp][0] + 0.01,
                     limit,
                     "$n$=10$^{" + str(
                         math.log10(self._bz.doping['n'][0])) + "}$",
                     color='b')
            plt.axvline(self._bz.mu_doping['n'][temp][-1], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['n'][temp][-1] + 0.01,
                     limit,
                     "$n$=10$^{" + str(math.log10(self._bz.doping['n'][-1]))
                     + "}$", color='b')
            plt.axvline(self._bz.mu_doping['p'][temp][0], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['p'][temp][0] + 0.01,
                     limit,
                     "$p$=10$^{" + str(
                         math.log10(self._bz.doping['p'][0])) + "}$",
                     color='b')
            plt.axvline(self._bz.mu_doping['p'][temp][-1], linewidth=3.0,
                        linestyle="--")
            plt.text(self._bz.mu_doping['p'][temp][-1] + 0.01,
                     limit, "$p$=10$^{" +
                     str(math.log10(self._bz.doping['p'][-1])) + "}$",
                     color='b')

    def _plot_bg_limits(self):
        plt.axvline(0.0, color='k', linewidth=3.0)
        plt.axvline(self._bz.gap, color='k', linewidth=3.0)

    def plot_seebeck_mu(self, temp=600, output='eig', xlim=None):
        """
        Plot the seebeck coefficient in function of Fermi level

        Args:
            temp:
                the temperature
            xlim:
                a list of min and max fermi energy by default (0, and band gap)
        Returns:
            a matplotlib object
        """
        seebeck = self._bz.get_seebeck(output=output, doping_levels=False)[temp]
        plt.plot(self._bz.mu_steps, seebeck,
                 linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['S$_1$', 'S$_2$', 'S$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim[0], xlim[1])
        plt.ylabel("Seebeck \n coefficient  ($\mu$V/K)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_conductivity_mu(self, temp=600, output='eig',
                             relaxation_time=1e-14, xlim=None):
        """
        Plot the conductivity in function of Fermi level. Semi-log plot

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        cond = self._bz.get_conductivity(relaxation_time=relaxation_time,
                                         output=output, doping_levels=False)[
            temp]
        plt.semilogy(self._bz.mu_steps, cond, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylim([1e13 * relaxation_time, 1e20 * relaxation_time])
        plt.ylabel("conductivity,\n $\sigma$ (1/($\Omega$ m))", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_power_factor_mu(self, temp=600, output='eig',
                             relaxation_time=1e-14, xlim=None):
        """
        Plot the power factor in function of Fermi level. Semi-log plot

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        pf = self._bz.get_power_factor(relaxation_time=relaxation_time,
                                       output=output, doping_levels=False)[temp]
        plt.semilogy(self._bz.mu_steps, pf, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['PF$_1$', 'PF$_2$', 'PF$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylabel("Power factor, ($\mu$W/(mK$^2$))", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_zt_mu(self, temp=600, output='eig', relaxation_time=1e-14,
                   xlim=None):
        """
        Plot the ZT in function of Fermi level.

        Args:
            temp: the temperature
            xlim: a list of min and max fermi energy by default (0, and band
                gap)
            tau: A relaxation time in s. By default none and the plot is by
               units of relaxation time

        Returns:
            a matplotlib object
        """
        zt = self._bz.get_zt(relaxation_time=relaxation_time, output=output,
                             doping_levels=False)[temp]
        plt.plot(self._bz.mu_steps, zt, linewidth=3.0)
        self._plot_bg_limits()
        self._plot_doping(temp)
        if output == 'eig':
            plt.legend(['ZT$_1$', 'ZT$_2$', 'ZT$_3$'])
        if xlim is None:
            plt.xlim(-0.5, self._bz.gap + 0.5)
        else:
            plt.xlim(xlim)
        plt.ylabel("ZT", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30.0)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_dos(self, sigma=0.05):
        """
        plot dos

        Args:
            sigma: a smearing

        Returns:
            a matplotlib object
        """
        plotter = DosPlotter(sigma=sigma)
        plotter.add_dos("t", self._bz.dos)
        return plotter.get_plot()

    def plot_carriers(self, temp=300):
        """
        Plot the carrier concentration in function of Fermi level

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        plt.semilogy(self._bz.mu_steps,
                     abs(self._bz.carrier_conc[temp] / (self._bz.vol * 1e-24)),
                     linewidth=3.0, color='r')
        self._plot_bg_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap + 0.5)
        plt.ylim(1e14, 1e22)
        plt.ylabel("carrier concentration (cm-3)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt

    def plot_hall_carriers(self, temp=300):
        """
        Plot the Hall carrier concentration in function of Fermi level

        Args:
            temp: the temperature

        Returns:
            a matplotlib object
        """
        hall_carriers = [abs(i) for i in self._bz.get_hall_carrier_concentration()[temp]]
        plt.semilogy(self._bz.mu_steps,
                     hall_carriers,
                     linewidth=3.0, color='r')
        self._plot_bg_limits()
        self._plot_doping(temp)
        plt.xlim(-0.5, self._bz.gap+0.5)
        plt.ylim(1e14, 1e22)
        plt.ylabel("Hall carrier concentration (cm-3)", fontsize=30.0)
        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        return plt
