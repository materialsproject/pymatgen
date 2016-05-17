# coding: utf-8

from __future__ import division, unicode_literals, print_function

import math
import os
import subprocess
import tempfile

import numpy as np
from monty.dev import requires
from monty.json import jsanitize
from monty.os import cd
from monty.os.path import which
from scipy.constants import e, m_e
from scipy.spatial import distance

from pymatgen.core.lattice import Lattice
from pymatgen.core.units import Energy, Length
from pymatgen.electronic_structure.bandstructure import \
    BandStructureSymmLine, Kpoint
from pymatgen.electronic_structure.core import Orbital
from pymatgen.electronic_structure.dos import Dos, Spin, CompleteDos
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

"""
This module provides classes to run and analyze boltztrap on pymatgen band
structure objects. Boltztrap is a software interpolating band structures and
computing materials properties from this band structure using Boltzmann
semi-classical transport theory.

Boltztrap has been developed by Georg Madsen.

http://www.icams.de/content/departments/ams/madsen/boltztrap.html

You need version 1.2.3 or higher

References are::

    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71
"""

__author__ = "Geoffroy Hautier, Zachary Gibbs, Francesco Ricci, Anubhav Jain"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "August 23, 2013"


class BoltztrapRunner(object):
    """
    This class is used to run Boltztrap on a band structure object.

    Args:
            bs:
                A band structure object
            nelec:
                the number of electrons
            dos_type:
                two options for the band structure integration: "HISTO"
                (histogram) or "TETRA" using the tetrahedon method. TETRA
                typically gives better results (especially for DOSes)
                but takes more time
            energy_grid:
                the energy steps used for the integration (eV)
            lpfac:
                the number of interpolation points in the real space. By
                default 10 gives 10 time more points in the real space than
                the number of kpoints given in reciprocal space
            run_type:
                type of boltztrap usage. by default
                - BOLTZ: (default) compute transport coefficients
                - BANDS: interpolate all bands contained in the energy range
                specified in energy_span_around_fermi variable, along specified
                k-points
                - DOS: compute total and partial dos (custom BoltzTraP code
                needed!)
                - FERMI: compute fermi surface or more correctly to
                get certain bands interpolated
            band_nb:
                indicates a band number. Used for Fermi Surface interpolation
                (run_type="FERMI")
            spin:
                specific spin component (1: up, -1: down) of the band selected
                in FERMI mode (mandatory).
            cond_band:
                if a conduction band is specified in FERMI mode,
                set this variable as True
            tauref:
                reference relaxation time. Only set to a value different than
                zero if we want to model beyond the constant relaxation time.
            tauexp:
                exponent for the energy in the non-constant relaxation time
                approach
            tauen:
                reference energy for the non-constant relaxation time approach
            soc:
                results from spin-orbit coupling (soc) computations give
                typically non-polarized (no spin up or down) results but single
                electron occupations. If the band structure comes from a soc
                computation, you should set soc to True (default False)
            doping:
                the fixed doping levels you want to compute. Boltztrap provides
                both transport values depending on electron chemical potential
                (fermi energy) and for a series of fixed carrier
                concentrations. By default, this is set to 1e16 to 1e22 in
                increments of factors of 10.
            energy_span_around_fermi:
                usually the interpolation is not needed on the entire energy
                range but on a specific range around the fermi level.
                This energy gives this range in eV. by default it is 1.5 eV.
                If DOS or BANDS type are selected, this range is automatically
                set to cover the entire energy range.
            scissor:
                scissor to apply to the band gap (eV). This applies a scissor
                operation moving the band edges without changing the band
                shape. This is useful to correct the often underestimated band
                gap in DFT. Default is 0.0 (no scissor)
            kpt_line:
                list/array of kpoints in fractional coordinates for BANDS mode
                calculation (standard path of high symmetry k-points is
                automatically set as default)
            tmax:
                Maximum temperature (K) for calculation (default=1300)
            tgrid:
                Temperature interval for calculation (default=100)

    """

    @requires(which('x_trans'),
              "BoltztrapRunner requires the executables 'x_trans' to be in "
              "the path. Please download the Boltztrap at "
              "http://www.icams.de/content/departments/ams/madsen/boltztrap"
              ".html and follow the instructions in the README to compile "
              "Bolztrap accordingly. Then add x_trans to your path")
    def __init__(self, bs, nelec, dos_type="HISTO", energy_grid=0.005,
                 lpfac=10, run_type="BOLTZ", band_nb=None, tauref=0, tauexp=0,
                 tauen=0, soc=False, doping=None, energy_span_around_fermi=1.5,
                 scissor=0.0, kpt_line=None, spin=None, cond_band=False,
                 tmax=1300, tgrid=100):
        self.lpfac = lpfac
        self._bs = bs
        self._nelec = nelec
        self.dos_type = dos_type
        self.energy_grid = energy_grid
        self.error = []
        self.run_type = run_type
        self.band_nb = band_nb
        self.spin = spin
        self.cond_band = cond_band
        self.tauref = tauref
        self.tauexp = tauexp
        self.tauen = tauen
        self.soc = soc
        self.kpt_line = kpt_line
        self.doping = doping or [1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22]
        self.energy_span_around_fermi = energy_span_around_fermi
        self.scissor = scissor
        self.tmax = tmax
        self.tgrid = tgrid

    def _make_energy_file(self, file_name):
        with open(file_name, 'w') as f:
            f.write("test\n")
            f.write(str(len(self._bs.kpoints)) + "\n")
            sign = 1.0
            if self.cond_band:
                sign = -1.0

            if self.run_type in ("DOS", "BANDS"):
                # automatically set energy span around Fermi
                const = Energy(2, "eV").to("Ry")
                emin_up = min([e_k[0] for e_k in self._bs.bands[Spin.up]])
                emax_up = max([e_k[0] for e_k in self._bs.bands[Spin.up]])

                if self._bs.is_spin_polarized:
                    emin_dw = min([e_k[0] for e_k in
                                   self._bs.bands[Spin.down]])
                    low_en_lim = Energy(min((emin_up, emin_dw)) -
                                        self._bs.efermi, "eV").to("Ry")

                    emax_dw = max([e_k[0] for e_k in
                                   self._bs.bands[Spin.down]])
                    high_en_lim = Energy(max((emax_up, emax_dw)) -
                                         self._bs.efermi, "eV").to("Ry")

                    self._ll = low_en_lim - const
                    self._hl = high_en_lim + const

                    en_range = Energy(max((abs(self._ll), abs(self._hl))),
                                      "Ry").to("eV")
                else:
                    en_range = Energy(
                        max((abs(emin_up - const), abs(emax_up + const))),
                        "Ry").to("eV")

                self.energy_span_around_fermi = en_range * 1.01
                print("energy_span_around_fermi = ",
                      self.energy_span_around_fermi)

            if self.run_type != "FERMI":
                for i in range(len(self._bs.kpoints)):
                    tmp_eigs = []
                    if self.run_type == "DOS":
                        spin_lst = [self.spin]
                    else:
                        spin_lst = self._bs.bands

                    for spin in spin_lst:
                        for j in range(
                                int(math.floor(self._bs.nb_bands * 0.9))):
                            tmp_eigs.append(
                                Energy(self._bs.bands[Spin(spin)][j][i] -
                                       self._bs.efermi, "eV").to("Ry"))
                    tmp_eigs.sort()

                    if self.run_type == "DOS" and self._bs.is_spin_polarized:
                        tmp_eigs.insert(0, self._ll)
                        tmp_eigs.append(self._hl)

                    f.write("%12.8f %12.8f %12.8f %d\n"
                            % (self._bs.kpoints[i].frac_coords[0],
                               self._bs.kpoints[i].frac_coords[1],
                               self._bs.kpoints[i].frac_coords[2],
                               len(tmp_eigs)))

                    for j in range(len(tmp_eigs)):
                        f.write("%18.8f\n" % (sign * float(tmp_eigs[j])))

            else:
                for i in range(len(self._bs.kpoints)):
                    tmp_eigs = []
                    tmp_eigs.append(Energy(
                        self._bs.bands[Spin(self.spin)][self.band_nb][i] -
                        self._bs.efermi, "eV").to("Ry"))
                    f.write("%12.8f %12.8f %12.8f %d\n"
                            % (self._bs.kpoints[i].frac_coords[0],
                               self._bs.kpoints[i].frac_coords[1],
                               self._bs.kpoints[i].frac_coords[2],
                               len(tmp_eigs)))
                    for j in range(len(tmp_eigs)):
                        f.write("%18.8f\n" % (sign * float(tmp_eigs[j])))

    def _make_struc_file(self, file_name):
        sym = SpacegroupAnalyzer(self._bs.structure, symprec=0.01)

        with open(file_name, 'w') as f:
            f.write("{} {}\n".format(self._bs.structure.composition.formula,
                    sym.get_spacegroup_symbol()))

            f.write("{}\n".format(self._bs.structure.lattice))

            ops = sym.get_symmetry_dataset()['rotations']
            f.write("{}\n".format(len(ops)))

            for c in ops:
                for row in c:
                    f.write("{}\n".format(" ".join(str(i) for i in row)))
                #f.write('\n'.join([' '.join([str(int(i)) for i in row])
                #                   for row in c]))
                #f.write('\n')

    # This function is useless in std version of BoltzTraP code
    # because x_trans script overwrite BoltzTraP.def
    def _make_def_file(self, def_file_name):
        with open(def_file_name, 'w') as f:
            so = ""
            if self._bs.is_spin_polarized or self.soc:
                so = "so"
            f.write("5, 'boltztrap.intrans',      'old',    'formatted',0\n" +
                    "6,'boltztrap.outputtrans',      'unknown',    "
                    "'formatted',0\n" +
                    "20,'boltztrap.struct',         'old',    'formatted',0\n"
                    + "10,'boltztrap.energy" + so + "',         'old',    "
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
                    "'formatted',0\n")

    # This function is useless in std version of BoltzTraP code
    # because x_trans script overwrite BoltzTraP.def
    def _make_proj_files(self, file_name, def_file_name):
        for o in Orbital:
            for site_nb in range(0, len(self._bs.structure.sites)):
                if o in self._bs._projections[Spin.up][0][0]:
                    with open(file_name + "_" + str(site_nb) + "_" + str(o),
                              'w') as f:
                        f.write(self._bs.structure.composition.formula + "\n")
                        f.write(str(len(self._bs.kpoints)) + "\n")
                        for i in range(len(self._bs.kpoints)):
                            tmp_proj = []
                            for j in range(
                                    int(math.floor(self._bs.nb_bands * 0.9))):
                                tmp_proj.append(
                                    self._bs._projections[Spin(self.spin)][j][
                                        i][o][site_nb])
                            # TODO deal with the sorting going on at
                            # the energy level!!!
                            # tmp_proj.sort()

                            if self.run_type == "DOS" and \
                                    self._bs.is_spin_polarized:
                                tmp_proj.insert(0, self._ll)
                                tmp_proj.append(self._hl)

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
                    "20,'boltztrap.struct',         'old',    'formatted',0\n"
                    + "10,'boltztrap.energy" + so + "',         'old',    "
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
                    "'formatted',0\n")
            i = 1000
            for o in Orbital:
                for site_nb in range(0, len(self._bs.structure.sites)):
                    if o in self._bs._projections[Spin.up][0][0]:
                        f.write(str(i) + ",\'" + "boltztrap.proj_" + str(
                            site_nb) + "_" + str(o.name) +
                                "\' \'old\', \'formatted\',0\n")
                        i += 1

    def _make_intrans_file(self, file_name):
        if self.run_type == "BOLTZ" or self.run_type == "DOS":
            with open(file_name, 'w') as fout:
                fout.write("GENE          # use generic interface\n")
                setgap = 0
                if self.scissor > 0.0001:
                    setgap = 1
                fout.write(
                    "1 0 %d %f         # iskip (not presently used) idebug "
                    "setgap shiftgap \n"
                    % (setgap, Energy(self.scissor, "eV").to("Ry")))
                fout.write(
                    "0.0 %f %f %6.1f     # Fermilevel (Ry),energygrid,energy "
                    "span around Fermilevel, number of electrons\n"
                    % (Energy(self.energy_grid, "eV").to("Ry"),
                       Energy(self.energy_span_around_fermi, "eV").to("Ry"),
                       self._nelec))
                fout.write(
                    "CALC                    # CALC (calculate expansion "
                    "coeff), NOCALC read from file\n")
                fout.write(
                    "%d                        # lpfac, number of latt-points "
                    "per k-point\n" % self.lpfac)
                fout.write(
                    "%s                     # run mode (only BOLTZ is "
                    "supported)\n" % self.run_type)
                fout.write(
                    ".15                       # (efcut) energy range of "
                    "chemical potential\n")
                fout.write(
                    "{} {}                  # Tmax, temperature grid\n".\
                    format(self.tmax, self.tgrid))
                fout.write(
                    "-1.  # energyrange of bands given DOS output sig_xxx and "
                    "dos_xxx (xxx is band number)\n")
                fout.write(self.dos_type + "\n")
                fout.write(
                    str(self.tauref) + " " + str(self.tauexp) + " " + str(
                        self.tauen) + " 0 0 0\n")
                fout.write(str(2 * len(self.doping)) + "\n")
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
                fout.write(str(
                    1) + "                        # actual band selected: " +
                           str(self.band_nb + 1) + " spin: " + str(self.spin))

        elif self.run_type == "BANDS":

            setgap = 0
            if self.scissor > 0.0001:
                setgap = 1
            if self.kpt_line is None:
                kpath = HighSymmKpath(self._bs.structure)
                self.kpt_line = [Kpoint(k, self._bs.structure.lattice) for k
                                 in
                                 kpath.get_kpoints(coords_are_cartesian=False)[
                                     0]]
                self.kpt_line = np.array(
                    [kp.frac_coords for kp in self.kpt_line])

            with open(file_name, 'w') as fout:
                fout.write("GENE          # use generic interface\n")
                fout.write(
                    "1 0 %d %f         # iskip (not presently used) idebug "
                    "setgap shiftgap \n"
                    % (setgap, Energy(self.scissor, "eV").to("Ry")))
                fout.write(
                    "0.0 %f %f %6.1f     # Fermilevel (Ry),energygrid,energy "
                    "span around Fermilevel, "
                    "number of electrons\n"
                    % (Energy(self.energy_grid, "eV").to("Ry"),
                       Energy(self.energy_span_around_fermi, "eV").to("Ry"),
                       self._nelec))
                fout.write(
                    "CALC                    # CALC (calculate expansion "
                    "coeff), NOCALC read from file\n")
                fout.write(
                    "%d                        # lpfac, number of latt-points "
                    "per k-point\n" % self.lpfac)
                fout.write(
                    "BANDS                     # run mode (only BOLTZ is "
                    "supported)\n")
                fout.write("P " + str(len(self.kpt_line)) + "\n")
                for kp in self.kpt_line:
                    fout.writelines([str(k) + ' ' for k in kp])
                    fout.write('\n')

    def _make_all_files(self, path):
        if self._bs.is_spin_polarized or self.soc:
            self._make_energy_file(os.path.join(path, "boltztrap.energyso"))
        else:
            self._make_energy_file(os.path.join(path, "boltztrap.energy"))

        self._make_struc_file(os.path.join(path, "boltztrap.struct"))
        self._make_intrans_file(os.path.join(path, "boltztrap.intrans"))
        self._make_def_file("BoltzTraP.def")
        if len(self._bs._projections) != 0 and self.run_type == "DOS":
            self._make_proj_files(os.path.join(path, "boltztrap.proj"),
                                  os.path.join(path, "BoltzTraP.def"))

    def run(self, path_dir=None, convergence=True):
        if self.run_type in ("BANDS", "DOS", "FERMI"):
            convergence = False

        if self.run_type == "BANDS" and self._bs.is_spin_polarized:
            print("Reminder: for run_type " + str(
                self.run_type) + " spin component are not separated!")

        if self.run_type in ("FERMI", "DOS") and self.spin is None:
            if self._bs.is_spin_polarized:
                raise BoltztrapError(
                    "Spin component must be specified for spin polarized "
                    "case!")
            else:
                self.spin = 1

        dir_bz_name = "boltztrap"
        path_dir_orig = path_dir
        if path_dir is None:
            temp_dir = tempfile.mkdtemp()
            path_dir = os.path.join(temp_dir, dir_bz_name)
        else:
            path_dir = os.path.abspath(
                os.path.join(path_dir_orig, dir_bz_name))
        if not os.path.exists(path_dir):
            os.mkdir(path_dir)
        else:
            for c in os.listdir(path_dir):
                os.remove(path_dir + "/" + c)

        with cd(path_dir):

            ######## convergence loop over energy_grid, lpfac and not on
            # eff_mass (as previously) ########################
            lpfac_start = self.lpfac
            converged = False

            while self.energy_grid > 0.00004:
                sigma_ratio = 1
                self.lpfac = lpfac_start

                print("lpfac, energy_grid: ", self.lpfac, self.energy_grid)

                while self.lpfac < 160:

                    self._make_all_files(path_dir)
                    if self._bs.is_spin_polarized or self.soc:
                        p = subprocess.Popen(["x_trans", "BoltzTraP", "-so"],
                                             stdout=subprocess.PIPE,
                                             stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                        p.wait()
                    else:
                        p = subprocess.Popen(["x_trans", "BoltzTraP"],
                                             stdout=subprocess.PIPE,
                                             stdin=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
                        p.wait()

                    for c in p.communicate():
                        print(c)
                        if "STOP error in factorization" in c:
                            raise BoltztrapError("STOP error in factorization")

                    with open(os.path.join(path_dir,
                                           dir_bz_name + ".outputtrans")) as f:
                        warning = False
                        for l in f:
                            if "WARNING" in l:
                                warning = True
                                break
                            if "Option unknown" in l:
                                raise BoltztrapError(
                                    "DOS mode needs a custom version of "
                                    "BoltzTraP code is needed")
                    if warning:
                        print("There was a warning! Increase lpfac to " +
                              str(self.lpfac + 10))
                        self.lpfac += 10
                        continue

                    if convergence:
                        analyzer = BoltztrapAnalyzer.from_files(path_dir)
                        doping_ok = True
                        for doping in ['n', 'p']:
                            for c in analyzer.mu_doping[doping]:
                                if len(analyzer.mu_doping[doping][c]) != len(
                                        analyzer.doping[doping]):
                                    doping_ok = False
                                    break
                                if doping == 'p' and \
                                                sorted(
                                                    analyzer.mu_doping[doping][
                                                        c], reverse=True) != \
                                                analyzer.mu_doping[doping][c]:
                                    doping_ok = False
                                    break
                                if doping == 'n' and sorted(
                                        analyzer.mu_doping[doping][c]) != \
                                        analyzer.mu_doping[doping][c]:
                                    doping_ok = False
                                    break

                        print('doping_ok', doping_ok)
                        if not doping_ok:
                            self.lpfac += 10
                            print("doping not ok, increase lpfac to " + str(
                                self.lpfac))
                            continue

                        converged = True
                        break

                    else:
                        converged = True
                        break

                if not converged:
                    self.energy_grid /= 10
                else:
                    break

            if not converged:
                raise BoltztrapError(
                    "Doping convergence not reached with lpfac=" + str(
                        self.lpfac)
                    + ", energy_grid=" + str(self.energy_grid))

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

    def __init__(self, gap=None, mu_steps=None, cond=None, seebeck=None,
                 kappa=None, hall=None, doping=None,
                 mu_doping=None, seebeck_doping=None, cond_doping=None,
                 kappa_doping=None,
                 hall_doping=None, dos=None, dos_partial=None,
                 carrier_conc=None, vol=None, warning=None,
                 bz_bands=None, bz_kpoints=None, fermi_surface_data=None):
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
            bz_bands: Data for interpolated bands on a k-point line
                (run_type=BANDS)
            bz_kpoints: k-point in reciprocal coordinates for interpolated
                bands (run_type=BANDS)
            fermi_surface_data: energy values in a 3D grid imported from the
                output .cube file using ase.io.cube.read_cube
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
        self._bz_bands = bz_bands
        self._bz_kpoints = bz_kpoints
        self.fermi_surface_data = fermi_surface_data

    @staticmethod
    def _make_boltztrap_analyzer_from_data(
            data_full=None, data_hall=None, data_dos=None,
            temperature_steps=None, mu_steps=None,
            efermi=None, gap=None, doping=None, data_doping_full=None,
            data_doping_hall=None, vol=None,
            warning=False, bz_bands=None, bz_kpoints=None, run_type="BOLTZ",
            dos_spin=None, fermi_surface_data=None):
        """
        Make a BoltztrapAnalyzer object from raw data typically parse from
        files.
        """
        if run_type == "BOLTZ":
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

            dos = Dos(efermi, dos_full['energy'],
                      {Spin.up: dos_full['density']})
            dos_partial = data_dos['partial']

            return BoltztrapAnalyzer(
                gap, mu_steps, cond, seebeck, kappa, hall, new_doping,
                mu_doping,
                seebeck_doping, cond_doping, kappa_doping, hall_doping, dos,
                dos_partial, carrier_conc, vol, warning)

        elif run_type == "BANDS":

            return BoltztrapAnalyzer(bz_bands=bz_bands, bz_kpoints=bz_kpoints)

        elif run_type == "DOS":
            if not dos_spin:
                print(
                    "Spin component set to up (dos_spin=1). Set dos_spin=-1 in from_files function if you want spin down")
                dos_spin = 1

            dos_full = {'energy': [], 'density': []}

            for t in data_dos['total']:
                dos_full['energy'].append(t[0])
                dos_full['density'].append(t[1])

            dos = Dos(efermi, dos_full['energy'],
                      {Spin(dos_spin): dos_full['density']})
            dos_partial = data_dos['partial']

            return BoltztrapAnalyzer(dos=dos, dos_partial=dos_partial)

        elif run_type == "FERMI":
            return BoltztrapAnalyzer(fermi_surface_data=fermi_surface_data)

    def get_symm_bands(self, structure, efermi, kpt_line=None,
                       labels_dict=None):
        """
            Function useful to read bands from Boltztap output and get a
            BandStructureSymmLine object comparable with that one from a DFT
            calculation (if the same kpt_line is provided). Default kpt_line
            and labels_dict is the standard path of high symmetry k-point for
            the specified structure. They could be extracted from the
            BandStructureSymmLine object that you want to compare with. efermi
            variable must be specified to create the BandStructureSymmLine
            object (usually it comes from DFT or Boltztrap calc)
        """
        try:
            if kpt_line is None:
                kpath = HighSymmKpath(structure)
                kpt_line = [Kpoint(k, structure.reciprocal_lattice) for k in
                            kpath.get_kpoints(coords_are_cartesian=False)[0]]
                labels_dict = {l: k for k, l in zip(
                    *kpath.get_kpoints(coords_are_cartesian=False)) if l}
                kpoints = [kp.frac_coords for kp in kpt_line]
            else:
                kpoints = [kp.frac_coords for kp in kpt_line]
                labels_dict = {k: labels_dict[k].frac_coords for k in
                               labels_dict}

            idx_list = []
            #       kpt_dense=np.array([kp for kp in self._bz_kpoints])
            for i, kp in enumerate(kpt_line):
                w = []
                prec = 1e-05
                while len(w) == 0:
                    w = np.where(np.all(
                        np.abs(kp.frac_coords - self._bz_kpoints) < [prec] * 3,
                        axis=1))[0]
                    prec *= 10

                # print( prec )
                idx_list.append([i, w[0]])

                # if len(w)>0:
                #     idx_list.append([i,w[0]])
                # else:
                #     w=np.where(np.all(np.abs(kp.frac_coords-self._bz_kpoints)
                # <[1e-04,1e-04,1e-04],axis=1))[0]
                #     idx_list.append([i,w[0]])

            idx_list = np.array(idx_list)
            # print( idx_list.shape )

            bands_dict = {Spin.up: (self._bz_bands * Energy(1, "Ry").to(
                "eV") + efermi).T[:, idx_list[:, 1]].tolist()}
            # bz_kpoints = bz_kpoints[idx_list[:,1]].tolist()

            sbs = BandStructureSymmLine(kpoints, bands_dict,
                                        structure.reciprocal_lattice, efermi,
                                        labels_dict=labels_dict)

            return sbs

        except:
            raise BoltztrapError(
                "Bands are not in output of BoltzTraP.\nBolztrapRunner have "
                "to be run with run_type=BANDS")

    def check_acc_bzt_bands(self, sbs_bz, sbs_ref, warn_thr=0.01):
        """
            Compare sbs_bz BandStructureSymmLine calculated with boltztrap with
            the sbs_ref BandStructureSymmLine as reference (from MP for
            instance), using correlation.
            See compare_sym_bands function doc
            Return a list of correlation values of the eigth bands with index
            that ranges from vbm_idx-3 to cbm_idx +3
            Return also two list of sum of relative errors for each branch of
            vbm and cbm and a boolean variable to signal the presence of not
            accurate bands around the gap. warn_thr is a threshold to get a
            warning in the accuracy of Boltztap interpolated bands.
        """
        if not sbs_ref.is_metal():
            vbm_idx = sbs_ref.get_vbm()['band_index'][Spin.up][-1]
            cbm_idx = sbs_ref.get_cbm()['band_index'][Spin.up][0]
            corr, werr_vbm = compare_sym_bands(sbs_bz, sbs_ref, vbm_idx)
            corr, werr_cbm = compare_sym_bands(sbs_bz, sbs_ref, cbm_idx)

            acc_err = False
            warn_thr = 0.01
            if any(corr[vbm_idx - 3:vbm_idx + 1] > warn_thr) or any(
                            corr[cbm_idx:cbm_idx + 4] > warn_thr):
                print("Warning! some bands around gap are not accurate")
                acc_err = True

            return corr, werr_vbm, werr_cbm, acc_err

        else:
            raise BoltztrapError(
                "band check implemented only for bandstructure with gap")

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
                                                   self._seebeck_doping,
                                                   output,
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
                                             self._seebeck_doping[doping][
                                                 temp][
                                                 i],
                                             self._seebeck_doping[doping][
                                                 temp][
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
                    thermal_conduct = (self._kappa_doping[doping][temp][i]
                                       - pf_tensor * temp) * relaxation_time
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
                                               thermal_conduct + kl *
                                               np.eye(3, 3))))
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
                                                           m_e)
                    elif output == 'eig':
                        result_doping[doping][temp].append(
                            sorted(np.linalg.eigh(np.linalg.inv(
                                np.array(self._cond_doping[doping][temp][i])) *
                                                  self.doping[doping][
                                                      i] * 10 ** 6 * e ** 2 \
                                                  / m_e)[0]))
                    else:
                        full_tensor = np.linalg.inv(
                            np.array(self._cond_doping[doping][temp][i])) \
                                      * self.doping[doping][
                                          i] * 10 ** 6 * e ** 2 / m_e
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
                        result[temp].append((full_tensor[temp][i][0][0]
                                             + full_tensor[temp][i][1][1]
                                             + full_tensor[temp][i][2][
                                                 2]) * multi / 3.0)
        return result

    def get_complete_dos(self, structure, analyzer_for_second_spin=None):
        """
            Gives a CompleteDos object with the DOS from the interpolated
            projected band structure
            Args:
                the structure (necessary to identify sites for projection)
                analyzer_for_second_spin must be specified to have a
                CompleteDos with both Spin components
            Returns:
                a CompleteDos object
            Example of use in case of spin polarized case:
            
                BoltztrapRunner(bs=bs,nelec=10,run_type="DOS",spin=1).run(path_dir='dos_up/')
                an_up=BoltztrapAnalyzer.from_files("dos_up/boltztrap/",dos_spin=1)

                BoltztrapRunner(bs=bs,nelec=10,run_type="DOS",spin=-1).run(path_dir='dos_dw/')
                an_dw=BoltztrapAnalyzer.from_files("dos_dw/boltztrap/",dos_spin=-1)

                cdos=an_up.get_complete_dos(bs.structure,an_dw)

        """
        pdoss = {}
        spin_1 = self.dos.densities.keys()[0]

        if analyzer_for_second_spin:
            if not np.all(self.dos.energies ==
                                  analyzer_for_second_spin.dos.energies):
                raise BoltztrapError(
                    "Dos merging error: energies of the two dos are different")

            spin_2 = analyzer_for_second_spin.dos.densities.keys()[0]
            if spin_1 == spin_2:
                raise BoltztrapError(
                    "Dos merging error: spin component are the same")

        for s in self._dos_partial:
            if structure.sites[int(s)] not in pdoss:
                pdoss[structure.sites[int(s)]] = {}
            for o in self._dos_partial[s]:
                if Orbital[o] not in pdoss[structure.sites[int(s)]]:
                    pdoss[structure.sites[int(s)]][Orbital[o]] = {}
                pdoss[structure.sites[int(s)]][Orbital[o]][
                    spin_1] = self._dos_partial[s][o]
                if analyzer_for_second_spin:
                    pdoss[structure.sites[int(s)]][Orbital[o]][
                        spin_2] = analyzer_for_second_spin._dos_partial[s][o]
        if analyzer_for_second_spin:
            tdos = Dos(self.dos.efermi, self.dos.energies,
                       {spin_1: self.dos.densities[spin_1],
                        spin_2: analyzer_for_second_spin.dos.densities[
                            spin_2]})
        else:
            tdos = self.dos

        return CompleteDos(structure, total_dos=tdos, pdoss=pdoss)

    def get_mu_bounds(self, temp=300):
        return min(self.mu_doping['p'][temp]), max(self.mu_doping['n'][temp])

    def get_carrier_concentration(self):
        """
        gives the carrier concentration (in cm^-3)

        Returns
            a dictionary {temp:[]} with an array of carrier concentration
            (in cm^-3) at each temperature
            The array relates to each step of electron chemical potential
        """

        return {temp: [1e24 * i / self.vol for i in self._carrier_conc[temp]]
                for temp in self._carrier_conc}

    def get_hall_carrier_concentration(self):
        """
        gives the Hall carrier concentration (in cm^-3). This is the trace of
        the Hall tensor (see Boltztrap source code) Hall carrier concentration
        are not always exactly the same than carrier concentration.

        Returns
            a dictionary {temp:[]} with an array of Hall carrier concentration
            (in cm^-3) at each temperature The array relates to each step of
            electron chemical potential
        """
        result = {temp: [] for temp in self._hall}
        for temp in self._hall:
            for i in self._hall[temp]:
                trace = (i[1][2][0] + i[2][0][1] + i[0][1][2]) / 3.0
                if trace != 0.0:
                    result[temp].append(1e-6 / (trace * e))
                else:
                    result[temp].append(0.0)
        return result

    @staticmethod
    def from_files(path_dir, dos_spin=None):
        """
        get a BoltztrapAnalyzer object from a set of files

        Args:
            path_dir: directory where the boltztrap files are

        Returns:
            a BoltztrapAnalyzer object

        """

        with open(os.path.join(path_dir, "boltztrap.outputtrans"), 'r') as f:
            for line in f:
                if "Calc type:" in line:
                    run_type = line.split()[-1]
                    print(run_type, " calc type found")
                    break

        if run_type == "BOLTZ":
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
                [data_dos['total'][i][0],
                 2 * data_dos['total'][i][1] / total_elec]
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
                                if tokens[2] not in data_dos['partial'][
                                    tokens[1]]:
                                    data_dos['partial'][tokens[1]][
                                        tokens[2]] = []
                                data_dos['partial'][tokens[1]][
                                    tokens[2]].append(
                                    2 * float(line.split()[1]))

            with open(os.path.join(path_dir, "boltztrap.outputtrans"),
                      'r') as f:
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
                with open(
                        os.path.join(path_dir, "boltztrap.condtens_fixdoping"),
                        'r') as f:
                    for line in f:
                        if not line.startswith("#") and len(line) > 2:
                            data_doping_full.append([float(c)
                                                     for c in line.split()])

                with open(
                        os.path.join(path_dir, "boltztrap.halltens_fixdoping"),
                        'r') as f:
                    for line in f:
                        if not line.startswith("#") and len(line) > 2:
                            data_doping_hall.append(
                                [float(c) for c in line.split()])

            with open(os.path.join(path_dir, "boltztrap.struct"), 'r') as f:
                tokens = f.readlines()
                vol = Lattice(
                    [[Length(float(tokens[i].split()[j]), "bohr").to("ang")
                      for j in range(3)] for i in range(1, 4)]).volume

            data_dos = {'total': [], 'partial': {}}
            with open(os.path.join(path_dir, "boltztrap.transdos"), 'r') as f:
                count_series = 0
                normalization_factor = None
                for line in f:
                    if not line.startswith(" #"):
                        data_dos['total'].append(
                            [Energy(float(line.split()[0]), "Ry").to("eV"),
                             float(line.split()[1])])
                    else:
                        count_series += 1
                    if count_series > 1:
                        break
            for file_name in os.listdir(path_dir):
                if file_name.endswith(
                        "transdos") and file_name != 'boltztrap.transdos':
                    tokens = file_name.split(".")[1].split("_")
                    with open(os.path.join(path_dir, file_name), 'r') as f:
                        for line in f:
                            if not line.startswith(" #"):
                                if tokens[1] not in data_dos['partial']:
                                    data_dos['partial'][tokens[1]] = {}
                                if tokens[2] not in data_dos['partial'][
                                    tokens[1]]:
                                    data_dos['partial'][tokens[1]][
                                        tokens[2]] = []
                                data_dos['partial'][tokens[1]][
                                    tokens[2]].append(float(line.split()[1]))

            return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
                data_full, data_hall, data_dos, sorted([t for t in t_steps]),
                sorted([Energy(m, "Ry").to("eV") for m in m_steps]), efermi,
                Energy(gap, "Ry").to("eV"),
                doping, data_doping_full, data_doping_hall, vol, warning)

        elif run_type == "FERMI":
            from ase.io.cube import read_cube
            if os.path.exists(os.path.join(path_dir, 'boltztrap_BZ.cube')):
                fs_data = read_cube(
                    str(os.path.join(path_dir, 'boltztrap_BZ.cube')),
                    read_data=True)
            if os.path.exists(os.path.join(path_dir, 'fort.30')):
                fs_data = read_cube(str(os.path.join(path_dir, 'fort.30')),
                                    read_data=True)
            else:
                raise BoltztrapError("No data file found for fermi surface")

            return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
                run_type="FERMI", fermi_surface_data=fs_data)

        elif run_type == "BANDS":
            bz_kpoints = np.loadtxt(
                os.path.join(path_dir, "boltztrap_band.dat"))[:, -3:]
            bz_bands = np.loadtxt(
                os.path.join(path_dir, "boltztrap_band.dat"))[:, 1:-6]
            return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
                run_type="BANDS",
                bz_bands=bz_bands, bz_kpoints=bz_kpoints)

        elif run_type == "DOS":
            data_dos = {'total': [], 'partial': {}}
            with open(os.path.join(path_dir, "boltztrap.transdos"), 'r') as f:
                count_series = 0
                for line in f:
                    if not line.startswith(" #"):
                        data_dos['total'].append(
                            [Energy(float(line.split()[0]), "Ry").to("eV"),
                             float(line.split()[1])])
                    else:
                        count_series += 1
                    if count_series > 1:
                        break
            print("here")
            tmp_data = np.array(data_dos['total'])
            tmp_den = np.trim_zeros(tmp_data[:, 1], 'f')[1:]
            lw_l = len(tmp_data[:, 1]) - len(tmp_den)
            tmp_ene = tmp_data[lw_l:, 0]
            tmp_den = np.trim_zeros(tmp_den, 'b')[:-1]
            hg_l = len(tmp_ene) - len(tmp_den)
            tmp_ene = tmp_ene[:-hg_l]
            tmp_data = np.vstack((tmp_ene, tmp_den)).T
            data_dos['total'] = tmp_data.tolist()

            for file_name in os.listdir(path_dir):
                if file_name.endswith(
                        "transdos") and file_name != 'boltztrap.transdos':
                    tokens = file_name.split(".")[1].split("_")
                    site = tokens[1]
                    orb = '_'.join(tokens[2:])
                    with open(os.path.join(path_dir, file_name), 'r') as f:
                        for line in f:
                            if not line.startswith(" #"):
                                if site not in data_dos['partial']:
                                    data_dos['partial'][site] = {}
                                if orb not in data_dos['partial'][site]:
                                    data_dos['partial'][site][orb] = []
                                data_dos['partial'][site][orb].append(
                                    float(line.split()[1]))
                    data_dos['partial'][site][orb] = data_dos['partial'][site][
                                                         orb][lw_l:-hg_l]

            with open(os.path.join(path_dir, "boltztrap.outputtrans"),
                      'r') as f:
                for line in f:
                    if line.startswith("VBM"):
                        efermi = Energy(line.split()[1], "Ry").to("eV")

            return BoltztrapAnalyzer._make_boltztrap_analyzer_from_data(
                efermi=efermi, run_type="DOS",
                data_dos=data_dos, dos_spin=dos_spin)

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


def compare_sym_bands(bands_obj, bands_ref_obj, nb=None):
    """
        Compute the mean of correlation between bzt and vasp bandstructure on
        sym line, for all bands and locally (for each branches) the difference
        squared (%) if nb is specified.
    """

    nkpt = len(bands_obj.kpoints)
    if bands_ref_obj.is_spin_polarized:
        nbands = min(bands_obj.nb_bands, 2 * bands_ref_obj.nb_bands)
    else:
        # TODO: why is this needed? Shouldn't pmg take care of nb_bands?
        nbands = min(len(bands_obj.bands[Spin.up]),
                     len(bands_ref_obj.bands[Spin.up]))
    # print(nbands)
    arr_bands = np.array(bands_obj.bands[Spin.up][:nbands])
    # arr_bands_lavg = (arr_bands-np.mean(arr_bands,axis=1).reshape(nbands,1))

    if bands_ref_obj.is_spin_polarized:
        arr_bands_ref_up = np.array(bands_ref_obj.bands[Spin.up])
        arr_bands_ref_dw = np.array(bands_ref_obj.bands[Spin.down])
        # print(arr_bands_ref_up.shape)
        arr_bands_ref = np.vstack((arr_bands_ref_up, arr_bands_ref_dw))
        arr_bands_ref = np.sort(arr_bands_ref, axis=0)[:nbands]
        # print(arr_bands_ref.shape)
    else:
        arr_bands_ref = np.array(bands_ref_obj.bands[Spin.up][:nbands])

    # arr_bands_ref_lavg =
    # (arr_bands_ref-np.mean(arr_bands_ref,axis=1).reshape(nbands,1))

    # err = np.sum((arr_bands_lavg-arr_bands_ref_lavg)**2,axis=1)/nkpt
    corr = np.array(
        [distance.correlation(arr_bands[idx], arr_bands_ref[idx]) for idx in
         range(nbands)])

    if type(nb) == int and nb < nbands:
        branches = [[s['start_index'], s['end_index'], s['name']] for s in
                    bands_ref_obj._branches]
        werr = {}
        for start, end, name in branches:
            # werr.append((sum((arr_bands_corr[nb][start:end+1] -
            # arr_bands_ref_corr[nb][start:end+1])**2)/(end+1-start)*100,name))
            werr[name] = np.sum(abs(
                arr_bands[nb][start:end + 1] - arr_bands_ref[nb][
                                               start:end + 1]) / abs(
                arr_bands_ref[nb][start:end + 1])) / (end + 1 - start) * 100
    else:
        werr = "No nb given"

    return corr, werr
