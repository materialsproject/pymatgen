"""This module provides classes to run and analyze boltztrap on pymatgen band
structure objects. Boltztrap is a software interpolating band structures and
computing materials properties from this band structure using Boltzmann
semi-classical transport theory.

Boltztrap has been developed by Georg Madsen.

http://www.icams.de/content/research/software-development/boltztrap/

You need version 1.2.3 or higher

References are::

    Madsen, G. K. H., and Singh, D. J. (2006).
    BoltzTraP. A code for calculating band-structure dependent quantities.
    Computer Physics Communications, 175, 67-71
"""

from __future__ import annotations

import logging
import math
import os
import subprocess
import tempfile
import time
from shutil import which
from typing import TYPE_CHECKING

import numpy as np
from monty.dev import requires
from monty.json import MSONable, jsanitize
from monty.os import cd
from scipy import constants
from scipy.spatial import distance

from pymatgen.core.lattice import Lattice
from pymatgen.core.units import Energy, Length
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine, Kpoint
from pymatgen.electronic_structure.core import Orbital
from pymatgen.electronic_structure.dos import CompleteDos, Dos, Spin
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from pymatgen.core.sites import PeriodicSite
    from pymatgen.core.structure import Structure

__author__ = "Geoffroy Hautier, Zachary Gibbs, Francesco Ricci, Anubhav Jain"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "1.1"
__maintainer__ = "Geoffroy Hautier"
__email__ = "geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "August 23, 2013"


class BoltztrapRunner(MSONable):
    """This class is used to run Boltztrap on a band structure object."""

    @requires(
        which("x_trans"),
        "BoltztrapRunner requires the executables 'x_trans' to be in PATH. Please download "
        "Boltztrap at http://www.icams.de/content/research/software-development/boltztrap/ "
        "and follow the instructions in the README to compile Bolztrap accordingly. "
        "Then add x_trans to your path",
    )
    def __init__(
        self,
        bs,
        nelec,
        dos_type="HISTO",
        energy_grid=0.005,
        lpfac=10,
        run_type="BOLTZ",
        band_nb=None,
        tauref=0,
        tauexp=0,
        tauen=0,
        soc=False,
        doping=None,
        energy_span_around_fermi=1.5,
        scissor=0.0,
        kpt_line=None,
        spin=None,
        cond_band=False,
        tmax=1300,
        tgrid=50,
        symprec=1e-3,
        cb_cut=10,
        timeout=7200,
    ):
        """
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
                range but on a specific range around the Fermi level.
                This energy gives this range in eV. by default it is 1.5 eV.
                If DOS or BANDS type are selected, this range is automatically
                set to cover the entire energy range.
            scissor:
                scissor to apply to the band gap (eV). This applies a scissor
                operation moving the band edges without changing the band
                shape. This is useful to correct the often underestimated band
                gap in DFT. Default is 0.0 (no scissor)
            kpt_line:
                list of fractional coordinates of kpoints as arrays or list of
                Kpoint objects for BANDS mode calculation (standard path of
                high symmetry k-points is automatically set as default)
            tmax:
                Maximum temperature (K) for calculation (default=1300)
            tgrid:
                Temperature interval for calculation (default=50)
            symprec: 1e-3 is the default in pymatgen. If the kmesh has been
                generated using a different symprec, it has to be specified
                to avoid a "factorization error" in BoltzTraP calculation.
                If a kmesh that spans the whole Brillouin zone has been used,
                or to disable all the symmetries, set symprec to None.
            cb_cut: by default 10% of the highest conduction bands are
                removed because they are often not accurate.
                Tune cb_cut to change the percentage (0-100) of bands
                that are removed.
            timeout: overall time limit (in seconds): mainly to avoid infinite
                loop when trying to find Fermi levels.
        """
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
        self.cb_cut = cb_cut / 100.0
        if isinstance(doping, list) and len(doping) > 0:
            self.doping = doping
        else:
            self.doping = []
            for d in [1e16, 1e17, 1e18, 1e19, 1e20, 1e21]:
                self.doping.extend([1 * d, 2.5 * d, 5 * d, 7.5 * d])
            self.doping.append(1e22)
        self.energy_span_around_fermi = energy_span_around_fermi
        self.scissor = scissor
        self.tmax = tmax
        self.tgrid = tgrid
        self._symprec = symprec
        if self.run_type in ("DOS", "BANDS"):
            self._auto_set_energy_range()
        self.timeout = timeout
        self.start_time = time.perf_counter()

    def _auto_set_energy_range(self):
        """Automatically determine the energy range as min/max eigenvalue
        minus/plus the buffer_in_ev.
        """
        emins = [min(e_k[0] for e_k in self._bs.bands[Spin.up])]
        emaxs = [max(e_k[0] for e_k in self._bs.bands[Spin.up])]

        if self._bs.is_spin_polarized:
            emins.append(min(e_k[0] for e_k in self._bs.bands[Spin.down]))

            emaxs.append(max(e_k[0] for e_k in self._bs.bands[Spin.down]))

        min_eigenval = Energy(min(emins) - self._bs.efermi, "eV").to("Ry")
        max_eigenval = Energy(max(emaxs) - self._bs.efermi, "eV").to("Ry")

        # set energy range to buffer around min/max EV
        # buffer does not increase CPU time but will help get equal
        # energies for spin up/down for band structure
        const = Energy(2, "eV").to("Ry")
        self._ll = min_eigenval - const
        self._hl = max_eigenval + const

        en_range = Energy(max((abs(self._ll), abs(self._hl))), "Ry").to("eV")

        self.energy_span_around_fermi = en_range * 1.01
        print("energy_span_around_fermi = ", self.energy_span_around_fermi)

    @property
    def bs(self):
        """The BandStructure."""
        return self._bs

    @property
    def nelec(self):
        """Number of electrons."""
        return self._nelec

    def write_energy(self, output_file):
        """Writes the energy to an output file.

        :param output_file: Filename
        """
        with open(output_file, "w") as f:
            f.write("test\n")
            f.write(f"{len(self._bs.kpoints)}\n")

            if self.run_type == "FERMI":
                sign = -1.0 if self.cond_band else 1.0
                for i, kpt in enumerate(self._bs.kpoints):
                    eigs = []
                    eigs.append(
                        Energy(
                            self._bs.bands[Spin(self.spin)][self.band_nb][i] - self._bs.efermi,
                            "eV",
                        ).to("Ry")
                    )
                    a, b, c = kpt.frac_coords
                    f.write(f"{a:12.8f} {b:12.8f} {c:12.8f}{len(eigs)}\n")
                    for e in eigs:
                        f.write(f"{sign * float(e):18.8f}\n")

            else:
                for i, kpt in enumerate(self._bs.kpoints):
                    eigs = []
                    spin_lst = [self.spin] if self.run_type == "DOS" else self._bs.bands

                    for spin in spin_lst:
                        # use 90% of bottom bands since highest eigenvalues
                        # are usually incorrect
                        # ask Geoffroy Hautier for more details
                        nb_bands = int(math.floor(self._bs.nb_bands * (1 - self.cb_cut)))
                        for j in range(nb_bands):
                            eigs.append(
                                Energy(
                                    self._bs.bands[Spin(spin)][j][i] - self._bs.efermi,
                                    "eV",
                                ).to("Ry")
                            )
                    eigs.sort()

                    if self.run_type == "DOS" and self._bs.is_spin_polarized:
                        eigs.insert(0, self._ll)
                        eigs.append(self._hl)
                    a, b, c = kpt.frac_coords
                    f.write(f"{a:12.8f} {b:12.8f} {c:12.8f} {len(eigs)}\n")

                    for e in eigs:
                        f.write(f"{float(e):18.8f}\n")

    def write_struct(self, output_file):
        """Writes the structure to an output file.

        :param output_file: Filename
        """
        if self._symprec is not None:
            sym = SpacegroupAnalyzer(self._bs.structure, symprec=self._symprec)
        elif self._symprec is None:
            pass

        with open(output_file, "w") as f:
            if self._symprec is not None:
                f.write(f"{self._bs.structure.composition.formula} {sym.get_space_group_symbol()}\n")
            elif self._symprec is None:
                f.write(f"{self._bs.structure.composition.formula} symmetries disabled\n")

            f.write(
                "\n".join(
                    " ".join(f"{Length(i, 'ang').to('bohr'):.5f}" for i in row)
                    for row in self._bs.structure.lattice.matrix
                )
                + "\n"
            )

            if self._symprec is not None:
                ops = sym.get_symmetry_dataset()["rotations"]
            elif self._symprec is None:
                ops = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
            f.write(f"{len(ops)}\n")

            for op in ops:
                for row in op:
                    f.write(f"{' '.join(map(str, row))}\n")

    def write_def(self, output_file):
        """Writes the def to an output file.

        :param output_file: Filename
        """
        # This function is useless in std version of BoltzTraP code
        # because x_trans script overwrite BoltzTraP.def
        with open(output_file, "w") as f:
            so = ""
            if self._bs.is_spin_polarized or self.soc:
                so = "so"
            f.write(
                "5, 'boltztrap.intrans',      'old',    'formatted',0\n"
                "6,'boltztrap.outputtrans',      'unknown',    "
                "'formatted',0\n"
                "20,'boltztrap.struct',         'old',    'formatted',0\n"
                "10,'boltztrap.energy" + so + "',         'old',    "
                "'formatted',0\n48,'boltztrap.engre',         'unknown',    "
                "'unformatted',0\n49,'boltztrap.transdos',        'unknown',    "
                "'formatted',0\n50,'boltztrap.sigxx',        'unknown',    'formatted',"
                "0\n51,'boltztrap.sigxxx',        'unknown',    'formatted',"
                "0\n21,'boltztrap.trace',           'unknown',    "
                "'formatted',0\n22,'boltztrap.condtens',           'unknown',    "
                "'formatted',0\n24,'boltztrap.halltens',           'unknown',    "
                "'formatted',0\n30,'boltztrap_BZ.cube',           'unknown',    "
                "'formatted',0\n"
            )

    def write_proj(self, output_file_proj, output_file_def):
        """Writes the projections to an output file.

        :param output_file: Filename
        """
        # This function is useless in std version of BoltzTraP code
        # because x_trans script overwrite BoltzTraP.def
        for oi, o in enumerate(Orbital):
            for site_nb in range(len(self._bs.structure)):
                if oi < len(self._bs.projections[Spin.up][0][0]):
                    with open(f"{output_file_proj}_{site_nb}_{o}", "w") as f:
                        f.write(self._bs.structure.composition.formula + "\n")
                        f.write(str(len(self._bs.kpoints)) + "\n")
                        for i, kpt in enumerate(self._bs.kpoints):
                            tmp_proj = []
                            for j in range(int(math.floor(self._bs.nb_bands * (1 - self.cb_cut)))):
                                tmp_proj.append(self._bs.projections[Spin(self.spin)][j][i][oi][site_nb])
                            # TODO deal with the sorting going on at
                            # the energy level!!!
                            # tmp_proj.sort()

                            if self.run_type == "DOS" and self._bs.is_spin_polarized:
                                tmp_proj.insert(0, self._ll)
                                tmp_proj.append(self._hl)
                            a, b, c = kpt.frac_coords
                            f.write(f"{a:12.8f} {b:12.8f} {c:12.8f} {len(tmp_proj)}\n")
                            for t in tmp_proj:
                                f.write(f"{float(t):18.8f}\n")
        with open(output_file_def, "w") as f:
            so = ""
            if self._bs.is_spin_polarized:
                so = "so"
            f.write(
                "5, 'boltztrap.intrans',      'old',    'formatted',0\n"
                "6,'boltztrap.outputtrans',      'unknown',    "
                "'formatted',0\n"
                "20,'boltztrap.struct',         'old',    'formatted',0\n"
                "10,'boltztrap.energy" + so + "',         'old',    "
                "'formatted',0\n48,'boltztrap.engre',         'unknown',    "
                "'unformatted',0\n49,'boltztrap.transdos',        'unknown',    "
                "'formatted',0\n50,'boltztrap.sigxx',        'unknown',    'formatted',"
                "0\n51,'boltztrap.sigxxx',        'unknown',    'formatted',"
                "0\n21,'boltztrap.trace',           'unknown',    "
                "'formatted',0\n22,'boltztrap.condtens',           'unknown',    "
                "'formatted',0\n24,'boltztrap.halltens',           'unknown',    "
                "'formatted',0\n30,'boltztrap_BZ.cube',           'unknown',    "
                "'formatted',0\n"
            )
            i = 1000
            for oi, o in enumerate(Orbital):
                for site_nb in range(len(self._bs.structure)):
                    if oi < len(self._bs.projections[Spin.up][0][0]):
                        f.write(f"{i},'boltztrap.proj_{site_nb}_{o.name}old', 'formatted',0\n")
                        i += 1

    def write_intrans(self, output_file):
        """Writes the intrans to an output file.

        :param output_file: Filename
        """
        setgap = 1 if self.scissor > 0.0001 else 0

        if self.run_type in ("BOLTZ", "DOS"):
            with open(output_file, "w") as fout:
                fout.write("GENE          # use generic interface\n")
                fout.write(
                    f"1 0 {setgap} {Energy(self.scissor, 'eV').to('Ry')}         "
                    "# iskip (not presently used) idebug setgap shiftgap \n"
                )
                fout.write(
                    f"0.0 {Energy(self.energy_grid, 'eV').to('Ry')} "
                    f"{Energy(self.energy_span_around_fermi, 'eV').to('Ry')} {self._nelec}.1f     "
                    f"# Fermilevel (Ry),energygrid,energy span around Fermilevel, number of electrons\n"
                )
                fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
                fout.write(f"{self.lpfac}                        # lpfac, number of latt-points per k-point\n")
                fout.write(f"{self.run_type}                     # run mode (only BOLTZ is supported)\n")
                fout.write(".15                       # (efcut) energy range of chemical potential\n")
                fout.write(f"{self.tmax} {self.tgrid}                  # Tmax, temperature grid\n")
                fout.write("-1.  # energyrange of bands given DOS output sig_xxx and dos_xxx (xxx is band number)\n")
                fout.write(self.dos_type + "\n")  # e.g., HISTO or TETRA
                fout.write(f"{self.tauref} {self.tauexp} {self.tauen} 0 0 0\n")
                fout.write(f"{2 * len(self.doping)}\n")

                for d in self.doping:
                    fout.write(str(d) + "\n")
                for d in self.doping:
                    fout.write(str(-d) + "\n")

        elif self.run_type == "FERMI":
            with open(output_file, "w") as fout:
                fout.write("GENE          # use generic interface\n")
                fout.write("1 0 0 0.0         # iskip (not presently used) idebug setgap shiftgap \n")
                fout.write(
                    f"0.0 {Energy(self.energy_grid, 'eV').to('Ry')} 0.1 {self._nelec:6.1f}     # Fermilevel "
                    "(Ry),energygrid,energy span around Fermilevel, number of electrons\n"
                )
                fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
                fout.write(f"{self.lpfac}                        # lpfac, number of latt-points per k-point\n")
                fout.write("FERMI                     # run mode (only BOLTZ is supported)\n")
                fout.write(f"1                        # actual band selected: {self.band_nb + 1} spin: {self.spin}")

        elif self.run_type == "BANDS":
            if self.kpt_line is None:
                kpath = HighSymmKpath(self._bs.structure)
                self.kpt_line = [
                    Kpoint(k, self._bs.structure.lattice) for k in kpath.get_kpoints(coords_are_cartesian=False)[0]
                ]
                self.kpt_line = [kp.frac_coords for kp in self.kpt_line]
            elif isinstance(self.kpt_line[0], Kpoint):
                self.kpt_line = [kp.frac_coords for kp in self.kpt_line]

            with open(output_file, "w") as fout:
                fout.write("GENE          # use generic interface\n")
                fout.write(
                    f"1 0 {setgap} {Energy(self.scissor, 'eV').to('Ry')}         # iskip (not presently used) "
                    "idebug setgap shiftgap \n"
                )
                fout.write(
                    f"0.0 {Energy(self.energy_grid, 'eV').to('Ry')} "
                    f"{Energy(self.energy_span_around_fermi, 'eV').to('Ry')} {self._nelec:.1f}     "
                    f"# Fermilevel (Ry),energygrid,energy span around Fermilevel, number of electrons\n"
                )
                fout.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
                fout.write(f"{self.lpfac}                        # lpfac, number of latt-points per k-point\n")
                fout.write("BANDS                     # run mode (only BOLTZ is supported)\n")
                fout.write(f"P {len(self.kpt_line)}\n")
                for kp in self.kpt_line:
                    fout.writelines([f"{k} " for k in kp])
                    fout.write("\n")

    def write_input(self, output_dir):
        """Writes the input files.

        :param output_dir: Directory to write the input files.
        """
        if self._bs.is_spin_polarized or self.soc:
            self.write_energy(f"{output_dir}/boltztrap.energyso")
        else:
            self.write_energy(f"{output_dir}/boltztrap.energy")

        self.write_struct(f"{output_dir}/boltztrap.struct")
        self.write_intrans(f"{output_dir}/boltztrap.intrans")
        self.write_def(f"{output_dir}/BoltzTraP.def")

        if len(self.bs.projections) != 0 and self.run_type == "DOS":
            self.write_proj(
                f"{output_dir}/boltztrap.proj",
                f"{output_dir}/BoltzTraP.def",
            )

    def run(
        self,
        path_dir=None,
        convergence=True,
        write_input=True,
        clear_dir=False,
        max_lpfac=150,
        min_egrid=0.00005,
    ):
        """Write inputs (optional), run BoltzTraP, and ensure
        convergence (optional).

        Args:
            path_dir (str): directory in which to run BoltzTraP
            convergence (bool): whether to check convergence and make
                corrections if needed
            write_input: (bool) whether to write input files before the run
                (required for convergence mode)
            clear_dir: (bool) whether to remove all files in the path_dir
                before starting
            max_lpfac: (float) maximum lpfac value to try before reducing egrid
                in convergence mode
            min_egrid: (float) minimum egrid value to try before giving up in
                convergence mode
        """
        # TODO: consider making this a part of custodian rather than pymatgen
        # A lot of this functionality (scratch dirs, handlers, monitors)
        # is built into custodian framework

        if convergence and not write_input:
            raise ValueError("Convergence mode requires write_input to be true")

        if self.run_type in ("BANDS", "DOS", "FERMI"):
            convergence = False
            if self.lpfac > max_lpfac:
                max_lpfac = self.lpfac

        if self.run_type == "BANDS" and self.bs.is_spin_polarized:
            print(
                f"Reminder: for run_type {self.run_type}, spin component are not separated! "
                "(you have a spin polarized band structure)"
            )

        if self.run_type in ("FERMI", "DOS") and self.spin is None:
            if self.bs.is_spin_polarized:
                raise BoltztrapError("Spin parameter must be specified for spin polarized band structures!")
            self.spin = 1

        dir_bz_name = "boltztrap"
        if path_dir is None:
            temp_dir = tempfile.mkdtemp()
            path_dir = os.path.join(temp_dir, dir_bz_name)
        else:
            path_dir = os.path.abspath(os.path.join(path_dir, dir_bz_name))

        if not os.path.exists(path_dir):
            os.mkdir(path_dir)
        elif clear_dir:
            for c in os.listdir(path_dir):
                os.remove(os.path.join(path_dir, c))

        FORMAT = "%(message)s"
        logging.basicConfig(
            level=logging.INFO,
            format=FORMAT,
            filename=f"{path_dir}/../boltztrap.out",
        )

        with cd(path_dir):
            lpfac_start = self.lpfac
            converged = False

            while self.energy_grid >= min_egrid and not converged:
                self.lpfac = lpfac_start
                if time.perf_counter() - self.start_time > self.timeout:
                    raise BoltztrapError(f"no doping convergence after timeout of {self.timeout} s")

                logging.info(f"lpfac, energy_grid: {self.lpfac} {self.energy_grid}")

                while self.lpfac <= max_lpfac and not converged:
                    if time.perf_counter() - self.start_time > self.timeout:
                        raise BoltztrapError(f"no doping convergence after timeout of {self.timeout} s")

                    if write_input:
                        self.write_input(path_dir)

                    bt_exe = ["x_trans", "BoltzTraP"]
                    if self._bs.is_spin_polarized or self.soc:
                        bt_exe.append("-so")

                    with subprocess.Popen(
                        bt_exe,
                        stdout=subprocess.PIPE,
                        stdin=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                    ) as p:
                        p.wait()

                        for c in p.communicate():
                            logging.info(c.decode())
                            if "error in factorization" in c.decode():
                                raise BoltztrapError("error in factorization")

                    warning = ""

                    with open(os.path.join(path_dir, dir_bz_name + ".outputtrans")) as file:
                        for line in file:
                            if "Option unknown" in line:
                                raise BoltztrapError("DOS mode needs a custom version of BoltzTraP code is needed")
                            if "WARNING" in line:
                                warning = line
                                break
                            if "Error - Fermi level was not found" in line:
                                warning = line
                                break

                    if not warning and convergence:
                        # check convergence for warning
                        analyzer = BoltztrapAnalyzer.from_files(path_dir)
                        for doping in ["n", "p"]:
                            for c in analyzer.mu_doping[doping]:
                                if len(analyzer.mu_doping[doping][c]) != len(analyzer.doping[doping]):
                                    warning = "length of mu_doping array is incorrect"
                                    break

                                if (
                                    doping == "p"
                                    and sorted(analyzer.mu_doping[doping][c], reverse=True)
                                    != analyzer.mu_doping[doping][c]
                                ):
                                    warning = "sorting of mu_doping array incorrect for p-type"
                                    break

                                # ensure n-type doping sorted correctly
                                if (
                                    doping == "n"
                                    and sorted(analyzer.mu_doping[doping][c]) != analyzer.mu_doping[doping][c]
                                ):
                                    warning = "sorting of mu_doping array incorrect for n-type"
                                    break

                    if warning:
                        self.lpfac += 10
                        logging.warning(f"Warning detected: {warning}! Increase lpfac to {self.lpfac}")

                    else:
                        converged = True

                if not converged:
                    self.energy_grid /= 10
                    logging.info(f"Could not converge with max lpfac; Decrease egrid to {self.energy_grid}")

            if not converged:
                lpfac, energy_grid = self.lpfac, self.energy_grid
                raise BoltztrapError(f"Doping convergence not reached with {lpfac=}, {energy_grid=}")

            return path_dir

    def as_dict(self):
        """MSONable dict."""
        results = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "lpfac": self.lpfac,
            "bs": self.bs.as_dict(),
            "nelec": self._nelec,
            "dos_type": self.dos_type,
            "run_type": self.run_type,
            "band_nb": self.band_nb,
            "spin": self.spin,
            "cond_band": self.cond_band,
            "tauref": self.tauref,
            "tauexp": self.tauexp,
            "tauen": self.tauen,
            "soc": self.soc,
            "kpt_line": self.kpt_line,
            "doping": self.doping,
            "energy_span_around_fermi": self.energy_span_around_fermi,
            "scissor": self.scissor,
            "tmax": self.tmax,
            "tgrid": self.tgrid,
            "symprec": self._symprec,
        }
        return jsanitize(results)


class BoltztrapError(Exception):
    """Exception class for boltztrap.
    Raised when the boltztrap gives an error.
    """


class BoltztrapAnalyzer:
    """Class used to store all the data from a boltztrap run."""

    def __init__(
        self,
        gap=None,
        mu_steps=None,
        cond=None,
        seebeck=None,
        kappa=None,
        hall=None,
        doping=None,
        mu_doping=None,
        seebeck_doping=None,
        cond_doping=None,
        kappa_doping=None,
        hall_doping=None,
        intrans=None,
        dos=None,
        dos_partial=None,
        carrier_conc=None,
        vol=None,
        warning=None,
        bz_bands=None,
        bz_kpoints=None,
        fermi_surface_data=None,
    ):
        """Constructor taking directly all the data generated by Boltztrap. You
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
                Fermi level in mu_steps]}. The units are 1/(Ohm*m*s).
            seebeck: The Seebeck tensor at different temperatures and fermi
                levels. The format is {temperature: [array of 3x3 tensors at
                each Fermi level in mu_steps]}. The units are V/K
            kappa: The electronic thermal conductivity tensor divided by a
                constant relaxation time (kappa/tau) at different temperature
                and fermi levels. The format is {temperature: [array of 3x3
                tensors at each Fermi level in mu_steps]}
                The units are W/(m*K*s)
            hall: The hall tensor at different temperature and fermi levels
                The format is {temperature: [array of 27 coefficients list at
                each Fermi level in mu_steps]}
                The units are m^3/C
            doping: The different doping levels that have been given to
                Boltztrap. The format is {'p':[],'n':[]} with an array of
                doping levels. The units are cm^-3
            mu_doping: Gives the electron chemical potential (or Fermi level)
                for a given set of doping.
                Format is {'p':{temperature: [fermi levels],'n':{temperature:
                [fermi levels]}}
                the Fermi level array is ordered according to the doping
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
            intrans: a dictionary of inputs e.g. {"scissor": 0.0}
            carrier_conc: The concentration of carriers in electron (or hole)
                per unit cell
            dos: The dos computed by Boltztrap given as a pymatgen Dos object
            dos_partial: Data for the partial DOS projected on sites and
                orbitals
            vol: Volume of the unit cell in angstrom cube (A^3)
            warning: string if BoltzTraP outputted a warning, else None
            bz_bands: Data for interpolated bands on a k-point line
                (run_type=BANDS)
            bz_kpoints: k-point in reciprocal coordinates for interpolated
                bands (run_type=BANDS)
            fermi_surface_data: energy values in a 3D grid imported from the
                output .cube file.
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
        self.intrans = intrans
        self._carrier_conc = carrier_conc
        self.dos = dos
        self.vol = vol
        self._dos_partial = dos_partial
        self._bz_bands = bz_bands
        self._bz_kpoints = bz_kpoints
        self.fermi_surface_data = fermi_surface_data

    def get_symm_bands(self, structure: Structure, efermi, kpt_line=None, labels_dict=None):
        """Function useful to read bands from Boltztrap output and get a
        BandStructureSymmLine object comparable with that one from a DFT
        calculation (if the same kpt_line is provided). Default kpt_line
        and labels_dict is the standard path of high symmetry k-point for
        the specified structure. They could be extracted from the
        BandStructureSymmLine object that you want to compare with. efermi
        variable must be specified to create the BandStructureSymmLine
        object (usually it comes from DFT or Boltztrap calc).
        """
        try:
            if kpt_line is None:
                kpath = HighSymmKpath(structure)
                kpt_line = [
                    Kpoint(kpt, structure.lattice.reciprocal_lattice)
                    for kpt in kpath.get_kpoints(coords_are_cartesian=False)[0]
                ]
                labels_dict = {
                    label: key for key, label in zip(*kpath.get_kpoints(coords_are_cartesian=False)) if label
                }
                kpt_line = [kp.frac_coords for kp in kpt_line]
            elif isinstance(kpt_line[0], Kpoint):
                kpt_line = [kp.frac_coords for kp in kpt_line]
                labels_dict = {k: labels_dict[k].frac_coords for k in labels_dict}

            _idx_list: list[tuple[int, ArrayLike]] = []
            for idx, kp in enumerate(kpt_line):
                w: list[bool] = []
                prec = 1e-5
                while len(w) == 0:
                    w = np.where(np.all(np.abs(kp - self._bz_kpoints) < [prec] * 3, axis=1))[0]  # type: ignore
                    prec *= 10

                _idx_list.append((idx, w[0]))

            idx_list = np.array(_idx_list)

            bz_bands_in_eV = (self._bz_bands * Energy(1, "Ry").to("eV") + efermi).T
            bands_dict = {Spin.up: bz_bands_in_eV[:, idx_list[:, 1]].tolist()}  # type: ignore

            return BandStructureSymmLine(
                kpt_line, bands_dict, structure.lattice.reciprocal_lattice, efermi, labels_dict=labels_dict
            )

        except Exception:
            raise BoltztrapError(
                "Bands are not in output of BoltzTraP.\nBolztrapRunner must be run with run_type=BANDS"
            )

    @staticmethod
    def check_acc_bzt_bands(sbs_bz, sbs_ref, warn_thr=(0.03, 0.03)):
        """Compare sbs_bz BandStructureSymmLine calculated with boltztrap with
        the sbs_ref BandStructureSymmLine as reference (from MP for
        instance), computing correlation and energy difference for eight bands
        around the gap (semiconductors) or Fermi level (metals).
        warn_thr is a threshold to get a warning in the accuracy of Boltztap
        interpolated bands.
        Return a dictionary with these keys:
        - "N": the index of the band compared; inside each there are:
            - "Corr": correlation coefficient for the 8 compared bands
            - "Dist": energy distance for the 8 compared bands
            - "branch_name": energy distance for that branch
        - "avg_corr": average of correlation coefficient over the 8 bands
        - "avg_dist": average of energy distance over the 8 bands
        - "nb_list": list of indexes of the 8 compared bands
        - "acc_thr": list of two float corresponding to the two warning
                     thresholds in input
        - "acc_err": list of two bools:
                     True if the avg_corr > warn_thr[0], and
                     True if the avg_dist > warn_thr[1]
        See also compare_sym_bands function doc.
        """
        if not sbs_ref.is_metal() and not sbs_bz.is_metal():
            vbm_idx = sbs_bz.get_vbm()["band_index"][Spin.up][-1]
            cbm_idx = sbs_bz.get_cbm()["band_index"][Spin.up][0]
            nb_list = range(vbm_idx - 3, cbm_idx + 4)

        else:
            bnd_around_efermi = []
            delta = 0
            spin = next(iter(sbs_bz.bands))
            while len(bnd_around_efermi) < 8 and delta < 100:
                delta += 0.1
                bnd_around_efermi = []
                for nb in range(len(sbs_bz.bands[spin])):
                    for kp in range(len(sbs_bz.bands[spin][nb])):
                        if abs(sbs_bz.bands[spin][nb][kp] - sbs_bz.efermi) < delta:
                            bnd_around_efermi.append(nb)
                            break
            if len(bnd_around_efermi) < 8:
                print(f"Warning! check performed on {len(bnd_around_efermi)}")
                nb_list = bnd_around_efermi
            else:
                nb_list = bnd_around_efermi[:8]

        bcheck = compare_sym_bands(sbs_bz, sbs_ref, nb_list)
        acc_err = [False, False]
        avg_corr = sum(item[1]["Corr"] for item in bcheck.items()) / 8
        avg_distance = sum(item[1]["Dist"] for item in bcheck.items()) / 8

        if avg_corr > warn_thr[0]:
            acc_err[0] = True
        if avg_distance > warn_thr[0]:
            acc_err[1] = True

        bcheck["avg_corr"] = avg_corr
        bcheck["avg_distance"] = avg_distance
        bcheck["acc_err"] = acc_err
        bcheck["acc_thr"] = warn_thr
        bcheck["nb_list"] = nb_list

        if True in acc_err:
            print("Warning! some bands around gap are not accurate")

        return bcheck

    def get_seebeck(self, output="eigs", doping_levels=True):
        """Gives the seebeck coefficient (microV/K) in either a
        full 3x3 tensor form, as 3 eigenvalues, or as the average value
        (trace/3.0) If doping_levels=True, the results are given at
        different p and n doping
        levels (given by self.doping), otherwise it is given as a series
        of electron chemical potential values.

        Args:
            output (str): the type of output. 'tensor' give the full
            3x3 tensor, 'eigs' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}.
            The 'p' links to Seebeck at p-type doping
            and 'n' to the Seebeck at n-type doping. Otherwise, returns a
            {temp:[]} dictionary
            The result contains either the sorted three eigenvalues of
            the symmetric
            Seebeck tensor (output='eigs') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').

            units are microV/K
        """
        return BoltztrapAnalyzer._format_to_output(self._seebeck, self._seebeck_doping, output, doping_levels, 1e6)

    def get_conductivity(self, output="eigs", doping_levels=True, relaxation_time=1e-14):
        """Gives the conductivity (1/Ohm*m) in either a full 3x3 tensor
        form, as 3 eigenvalues, or as the average value
        (trace/3.0) If doping_levels=True, the results are given at
        different p and n doping
        levels (given by self.doping), otherwise it is given as a series
        of electron chemical potential values.

        Args:
            output (str): the type of output. 'tensor' give the full
            3x3 tensor, 'eigs' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials
            relaxation_time (float): constant relaxation time in secs

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}.
            The 'p' links to conductivity
            at p-type doping and 'n' to the conductivity at n-type
            doping. Otherwise,
            returns a {temp:[]} dictionary. The result contains either
            the sorted three eigenvalues of the symmetric
            conductivity tensor (format='eigs') or a full tensor (3x3
            array) (output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time

            units are 1/Ohm*m
        """
        return BoltztrapAnalyzer._format_to_output(
            self._cond, self._cond_doping, output, doping_levels, relaxation_time
        )

    def get_power_factor(self, output="eigs", doping_levels=True, relaxation_time=1e-14):
        """Gives the power factor (Seebeck^2 * conductivity) in units
        microW/(m*K^2) in either a full 3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values.

        Args:
            output (str): the type of output. 'tensor' give the full 3x3
            tensor, 'eigs' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials
            relaxation_time (float): constant relaxation time in secs

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}. The
            'p' links to power factor
            at p-type doping and 'n' to the conductivity at n-type doping.
            Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            power factor tensor (format='eigs') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time

            units are microW/(m K^2)
        """
        result = result_doping = None
        if doping_levels:
            result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]} for doping in self._seebeck_doping}

            for doping in result_doping:
                for t in result_doping[doping]:
                    for i in range(len(self.doping[doping])):
                        full_tensor = np.dot(
                            self._cond_doping[doping][t][i],
                            np.dot(
                                self._seebeck_doping[doping][t][i],
                                self._seebeck_doping[doping][t][i],
                            ),
                        )
                        result_doping[doping][t].append(full_tensor)

        else:
            result = {t: [] for t in self._seebeck}
            for t in result:
                for i in range(len(self.mu_steps)):
                    full_tensor = np.dot(
                        self._cond[t][i],
                        np.dot(self._seebeck[t][i], self._seebeck[t][i]),
                    )
                    result[t].append(full_tensor)

        return BoltztrapAnalyzer._format_to_output(
            result, result_doping, output, doping_levels, multi=1e6 * relaxation_time
        )

    def get_thermal_conductivity(self, output="eigs", doping_levels=True, k_el=True, relaxation_time=1e-14):
        """Gives the electronic part of the thermal conductivity in either a
        full 3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values.

        Args:
            output (str): the type of output. 'tensor' give the full 3x3
            tensor, 'eigs' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials
            k_el (bool): True for k_0-PF*T, False for k_0
            relaxation_time (float): constant relaxation time in secs

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}. The
            'p' links to thermal conductivity
            at p-type doping and 'n' to the thermal conductivity at n-type
            doping. Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            conductivity tensor (format='eigs') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time

            units are W/mK
        """
        result = result_doping = None
        if doping_levels:
            result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]} for doping in self._seebeck_doping}
            for doping in result_doping:
                for t in result_doping[doping]:
                    for i in range(len(self.doping[doping])):
                        if k_el:
                            pf_tensor = np.dot(
                                self._cond_doping[doping][t][i],
                                np.dot(
                                    self._seebeck_doping[doping][t][i],
                                    self._seebeck_doping[doping][t][i],
                                ),
                            )
                            result_doping[doping][t].append(self._kappa_doping[doping][t][i] - pf_tensor * t)
                        else:
                            result_doping[doping][t].append(self._kappa_doping[doping][t][i])
        else:
            result = {t: [] for t in self._seebeck}
            for t in result:
                for i in range(len(self.mu_steps)):
                    if k_el:
                        pf_tensor = np.dot(
                            self._cond[t][i],
                            np.dot(self._seebeck[t][i], self._seebeck[t][i]),
                        )
                        result[t].append(self._kappa[t][i] - pf_tensor * t)
                    else:
                        result[t].append(self._kappa[t][i])

        return BoltztrapAnalyzer._format_to_output(result, result_doping, output, doping_levels, multi=relaxation_time)

    def get_zt(self, output="eigs", doping_levels=True, relaxation_time=1e-14, k_l=1):
        """Gives the ZT coefficient (S^2*cond*T/thermal cond) in either a full
        3x3 tensor form,
        as 3 eigenvalues, or as the average value (trace/3.0) If
        doping_levels=True, the results are given at
        different p and n doping levels (given by self.doping), otherwise it
        is given as a series of
        electron chemical potential values. We assume a constant relaxation
        time and a constant
        lattice thermal conductivity.

        Args:
            output (str): the type of output. 'tensor' give the full 3x3
            tensor, 'eigs' its 3 eigenvalues and
            'average' the average of the three eigenvalues
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials
            relaxation_time (float): constant relaxation time in secs
            k_l (float): lattice thermal cond in W/(m*K)

        Returns:
            If doping_levels=True, a dictionary {temp:{'p':[],'n':[]}}. The
            'p' links to ZT
            at p-type doping and 'n' to the ZT at n-type doping. Otherwise,
            returns a {temp:[]} dictionary. The result contains either the
            sorted three eigenvalues of the symmetric
            ZT tensor (format='eigs') or a full tensor (3x3 array) (
            output='tensor') or as an average
            (output='average').
            The result includes a given constant relaxation time and lattice
            thermal conductivity
        """
        result = result_doping = None
        if doping_levels:
            result_doping = {doping: {t: [] for t in self._seebeck_doping[doping]} for doping in self._seebeck_doping}

            for doping in result_doping:
                for t in result_doping[doping]:
                    for i in range(len(self.doping[doping])):
                        pf_tensor = np.dot(
                            self._cond_doping[doping][t][i],
                            np.dot(
                                self._seebeck_doping[doping][t][i],
                                self._seebeck_doping[doping][t][i],
                            ),
                        )
                        thermal_conduct = (self._kappa_doping[doping][t][i] - pf_tensor * t) * relaxation_time
                        result_doping[doping][t].append(
                            np.dot(
                                pf_tensor * relaxation_time * t,
                                np.linalg.inv(thermal_conduct + k_l * np.eye(3, 3)),
                            )
                        )
        else:
            result = {t: [] for t in self._seebeck}
            for t in result:
                for i in range(len(self.mu_steps)):
                    pf_tensor = np.dot(
                        self._cond[t][i],
                        np.dot(self._seebeck[t][i], self._seebeck[t][i]),
                    )
                    thermal_conduct = (self._kappa[t][i] - pf_tensor * t) * relaxation_time
                    result[t].append(
                        np.dot(
                            pf_tensor * relaxation_time * t,
                            np.linalg.inv(thermal_conduct + k_l * np.eye(3, 3)),
                        )
                    )

        return BoltztrapAnalyzer._format_to_output(result, result_doping, output, doping_levels)

    def get_average_eff_mass(self, output="eigs", doping_levels=True):
        """Gives the average effective mass tensor. We call it average because
        it takes into account all the bands
        and regions in the Brillouin zone. This is different than the standard
        textbook effective mass which relates
        often to only one (parabolic) band.
        The average effective mass tensor is defined as the integrated
        average of the second derivative of E(k)
        This effective mass tensor takes into account:
        -non-parabolicity
        -multiple extrema
        -multiple bands.

        For more information about it. See:

        Hautier, G., Miglio, A., Waroquiers, D., Rignanese, G., & Gonze,
        X. (2014).
        How Does Chemistry Influence Electron Effective Mass in Oxides?
        A High-Throughput Computational Analysis. Chemistry of Materials,
        26(19), 5447-5458. doi:10.1021/cm404079a

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
            output (str): 'eigs' for eigenvalues, 'tensor' for the full
            tensor and 'average' for an average (trace/3)
            doping_levels (bool): True for the results to be given at
            different doping levels, False for results
            at different electron chemical potentials

        Returns:
            If doping_levels=True,a dictionary {'p':{temp:[]},'n':{temp:[]}}
            with an array of effective mass tensor, eigenvalues of average
            value (depending on output) for each temperature and for each
            doping level.
            The 'p' links to hole effective mass tensor and 'n' to electron
            effective mass tensor.
        """
        result = result_doping = None
        conc = self.get_carrier_concentration()
        if doping_levels:
            result_doping = {doping: {t: [] for t in self._cond_doping[doping]} for doping in self.doping}
            for doping in result_doping:
                for temp in result_doping[doping]:
                    for i in range(len(self.doping[doping])):
                        try:
                            result_doping[doping][temp].append(
                                np.linalg.inv(np.array(self._cond_doping[doping][temp][i]))
                                * self.doping[doping][i]
                                * 10**6
                                * constants.e**2
                                / constants.m_e
                            )
                        except np.linalg.LinAlgError:
                            pass
        else:
            result = {t: [] for t in self._seebeck}
            for temp in result:
                for i in range(len(self.mu_steps)):
                    try:
                        cond_inv = np.linalg.inv(np.array(self._cond[temp][i]))
                    except np.linalg.LinAlgError:
                        pass
                    result[temp].append(cond_inv * conc[temp][i] * 10**6 * constants.e**2 / constants.m_e)

        return BoltztrapAnalyzer._format_to_output(result, result_doping, output, doping_levels)

    def get_seebeck_eff_mass(self, output="average", temp=300, doping_levels=False, Lambda=0.5):
        """Seebeck effective mass calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
        from ab initio band structure calculations.
        npj Computational Materials 3, 8 (2017).

        Args:
            output: 'average' returns the seebeck effective mass calculated using
                    the average of the three diagonal components of the seebeck tensor.
                    'tensor' returns the seebeck effective mass respect to the three
                    diagonal components of the seebeck tensor.
            doping_levels: False means that the seebeck effective mass is calculated
                           for every value of the chemical potential
                           True means that the seebeck effective mass is calculated
                           for every value of the doping levels for both n and p types
            temp:   temperature of calculated seebeck.
            Lambda: fitting parameter used to model the scattering (0.5 means constant
                    relaxation time).

        Returns:
            a list of values for the seebeck effective mass w.r.t the chemical potential,
            if doping_levels is set at False;
            a dict with n an p keys that contain a list of values for the seebeck effective
            mass w.r.t the doping levels, if doping_levels is set at True;
            if 'tensor' is selected, each element of the lists is a list containing
            the three components of the seebeck effective mass.
        """
        if doping_levels:
            sbk_mass = {}
            for dt in ("n", "p"):
                conc = self.doping[dt]
                seebeck = self.get_seebeck(output=output, doping_levels=True)[dt][temp]
                sbk_mass[dt] = []
                for i, c in enumerate(conc):
                    if output == "average":
                        sbk_mass[dt].append(seebeck_eff_mass_from_seebeck_carr(abs(seebeck[i]), c, temp, Lambda))
                    elif output == "tensor":
                        sbk_mass[dt].append([])
                        for j in range(3):
                            sbk_mass[dt][-1].append(
                                seebeck_eff_mass_from_seebeck_carr(abs(seebeck[i][j][j]), c, temp, Lambda)
                            )

        else:
            seebeck = self.get_seebeck(output=output, doping_levels=False)[temp]
            conc = self.get_carrier_concentration()[temp]
            sbk_mass = []
            for i, c in enumerate(conc):
                if output == "average":
                    sbk_mass.append(seebeck_eff_mass_from_seebeck_carr(abs(seebeck[i]), c, temp, Lambda))
                elif output == "tensor":
                    sbk_mass.append([])
                    for j in range(3):
                        sbk_mass[-1].append(seebeck_eff_mass_from_seebeck_carr(abs(seebeck[i][j][j]), c, temp, Lambda))
        return sbk_mass

    def get_complexity_factor(self, output="average", temp=300, doping_levels=False, Lambda=0.5):
        """Fermi surface complexity factor respect to calculated as explained in Ref.
        Gibbs, Z. M. et al., Effective mass and fermi surface complexity factor
        from ab initio band structure calculations.
        npj Computational Materials 3, 8 (2017).

        Args:
            output: 'average' returns the complexity factor calculated using the average
                    of the three diagonal components of the seebeck and conductivity tensors.
                    'tensor' returns the complexity factor respect to the three
                    diagonal components of seebeck and conductivity tensors.
            doping_levels: False means that the complexity factor is calculated
                           for every value of the chemical potential
                           True means that the complexity factor is calculated
                           for every value of the doping levels for both n and p types
            temp:   temperature of calculated seebeck and conductivity.
            Lambda: fitting parameter used to model the scattering (0.5 means constant
                    relaxation time).

        Returns:
            a list of values for the complexity factor w.r.t the chemical potential,
            if doping_levels is set at False;
            a dict with n an p keys that contain a list of values for the complexity factor
            w.r.t the doping levels, if doping_levels is set at True;
            if 'tensor' is selected, each element of the lists is a list containing
            the three components of the complexity factor.
        """
        if doping_levels:
            cmplx_fact = {}
            for dt in ("n", "p"):
                sbk_mass = self.get_seebeck_eff_mass(output, temp, doping_levels=True, Lambda=Lambda)[dt]
                cond_mass = self.get_average_eff_mass(output=output, doping_levels=True)[dt][temp]

                if output == "average":
                    cmplx_fact[dt] = [(m_s / abs(m_c)) ** 1.5 for m_s, m_c in zip(sbk_mass, cond_mass)]
                elif output == "tensor":
                    cmplx_fact[dt] = []
                    for i, sm in enumerate(sbk_mass):
                        cmplx_fact[dt].append([])
                        for j in range(3):
                            cmplx_fact[dt][-1].append((sm[j] / abs(cond_mass[i][j][j])) ** 1.5)

        else:
            sbk_mass = self.get_seebeck_eff_mass(output, temp, doping_levels=False, Lambda=Lambda)
            cond_mass = self.get_average_eff_mass(output=output, doping_levels=False)[temp]

            if output == "average":
                cmplx_fact = [(m_s / abs(m_c)) ** 1.5 for m_s, m_c in zip(sbk_mass, cond_mass)]
            elif output == "tensor":
                cmplx_fact = []
                for i, sm in enumerate(sbk_mass):
                    cmplx_fact.append([])
                    for j in range(3):
                        cmplx_fact[-1].append((sm[j] / abs(cond_mass[i][j][j])) ** 1.5)

        return cmplx_fact

    def get_extreme(
        self,
        target_prop,
        maximize=True,
        min_temp=None,
        max_temp=None,
        min_doping=None,
        max_doping=None,
        isotropy_tolerance=0.05,
        use_average=True,
    ):
        """This method takes in eigenvalues over a range of carriers,
        temperatures, and doping levels, and tells you what is the "best"
        value that can be achieved for the given target_property. Note that
        this method searches the doping dict only, not the full mu dict.

        Args:
            target_prop: target property, i.e. "seebeck", "power factor",
                         "conductivity", "kappa", or "zt"
            maximize: True to maximize, False to minimize (e.g. kappa)
            min_temp: minimum temperature allowed
            max_temp: maximum temperature allowed
            min_doping: minimum doping allowed (e.g., 1E18)
            max_doping: maximum doping allowed (e.g., 1E20)
            isotropy_tolerance: tolerance for isotropic (0.05 = 5%)
            use_average: True for avg of eigenval, False for max eigenval

        Returns:
            A dictionary with keys {"p", "n", "best"} with sub-keys:
            {"value", "temperature", "doping", "isotropic"}
        """

        def is_isotropic(x, isotropy_tolerance) -> bool:
            """Internal method to tell you if 3-vector "x" is isotropic.

            Args:
                x: the vector to determine isotropy for
                isotropy_tolerance: tolerance, e.g. 0.05 is 5%
            """
            if len(x) != 3:
                raise ValueError("Invalid input to is_isotropic!")

            st = sorted(x)
            return bool(
                all([st[0], st[1], st[2]])
                and (abs((st[1] - st[0]) / st[1]) <= isotropy_tolerance)
                and (abs(st[2] - st[0]) / st[2] <= isotropy_tolerance)
                and (abs((st[2] - st[1]) / st[2]) <= isotropy_tolerance)
            )

        if target_prop.lower() == "seebeck":
            d = self.get_seebeck(output="eigs", doping_levels=True)

        elif target_prop.lower() == "power factor":
            d = self.get_power_factor(output="eigs", doping_levels=True)

        elif target_prop.lower() == "conductivity":
            d = self.get_conductivity(output="eigs", doping_levels=True)

        elif target_prop.lower() == "kappa":
            d = self.get_thermal_conductivity(output="eigs", doping_levels=True)
        elif target_prop.lower() == "zt":
            d = self.get_zt(output="eigs", doping_levels=True)

        else:
            raise ValueError(f"Target property: {target_prop} not recognized!")

        absval = True  # take the absolute value of properties

        x_val = x_temp = x_doping = x_isotropic = None
        output = {}

        min_temp = min_temp or 0
        max_temp = max_temp or float("inf")
        min_doping = min_doping or 0
        max_doping = max_doping or float("inf")

        for pn in ("p", "n"):
            for t in d[pn]:  # temperatures
                if min_temp <= float(t) <= max_temp:
                    for didx, evs in enumerate(d[pn][t]):
                        doping_lvl = self.doping[pn][didx]
                        if min_doping <= doping_lvl <= max_doping:
                            isotropic = is_isotropic(evs, isotropy_tolerance)
                            if absval:
                                evs = [abs(x) for x in evs]
                            val = float(sum(evs)) / len(evs) if use_average else max(evs)
                            if x_val is None or (val > x_val and maximize) or (val < x_val and not maximize):
                                x_val = val
                                x_temp = t
                                x_doping = doping_lvl
                                x_isotropic = isotropic

            output[pn] = {
                "value": x_val,
                "temperature": x_temp,
                "doping": x_doping,
                "isotropic": x_isotropic,
            }
            x_val = None

        if maximize:
            max_type = "p" if output["p"]["value"] >= output["n"]["value"] else "n"
        else:
            max_type = "p" if output["p"]["value"] <= output["n"]["value"] else "n"

        output["best"] = output[max_type]
        output["best"]["carrier_type"] = max_type

        return output

    @staticmethod
    def _format_to_output(tensor, tensor_doping, output, doping_levels, multi=1.0):
        if doping_levels:
            full_tensor = tensor_doping
            result = {doping: {t: [] for t in tensor_doping[doping]} for doping in tensor_doping}
            for doping in full_tensor:
                for temp in full_tensor[doping]:
                    for i in range(len(full_tensor[doping][temp])):
                        if output in ["eig", "eigs"]:
                            result[doping][temp].append(sorted(np.linalg.eigh(full_tensor[doping][temp][i])[0] * multi))
                        elif output == "tensor":
                            result[doping][temp].append(np.array(full_tensor[doping][temp][i]) * multi)
                        elif output == "average":
                            result[doping][temp].append(
                                (
                                    full_tensor[doping][temp][i][0][0]
                                    + full_tensor[doping][temp][i][1][1]
                                    + full_tensor[doping][temp][i][2][2]
                                )
                                * multi
                                / 3.0
                            )
                        else:
                            raise ValueError(f"Unknown output format: {output}")
        else:
            full_tensor = tensor
            result = {t: [] for t in tensor}
            for temp in full_tensor:
                for i in range(len(tensor[temp])):
                    if output in ["eig", "eigs"]:
                        result[temp].append(sorted(np.linalg.eigh(full_tensor[temp][i])[0] * multi))
                    elif output == "tensor":
                        result[temp].append(np.array(full_tensor[temp][i]) * multi)
                    elif output == "average":
                        result[temp].append(
                            (full_tensor[temp][i][0][0] + full_tensor[temp][i][1][1] + full_tensor[temp][i][2][2])
                            * multi
                            / 3.0
                        )
                    else:
                        raise ValueError(f"Unknown output format: {output}")
        return result

    def get_complete_dos(self, structure: Structure, analyzer_for_second_spin=None):
        """Gives a CompleteDos object with the DOS from the interpolated projected band structure.

        Args:
            structure: necessary to identify sites for projection
            analyzer_for_second_spin: must be specified to have a CompleteDos with both Spin components

        Returns:
            a CompleteDos object

        Example of use in case of spin polarized case:

            BoltztrapRunner(bs=bs,nelec=10,run_type="DOS",spin=1).run(path_dir='dos_up/')
            an_up=BoltztrapAnalyzer.from_files("dos_up/boltztrap/",dos_spin=1)

            BoltztrapRunner(bs=bs,nelec=10,run_type="DOS",spin=-1).run(path_dir='dos_dw/')
            an_dw=BoltztrapAnalyzer.from_files("dos_dw/boltztrap/",dos_spin=-1)

            cdos=an_up.get_complete_dos(bs.structure,an_dw)
        """
        pdoss: dict[PeriodicSite, dict[Orbital, dict[Spin, ArrayLike]]] = {}
        spin_1 = next(iter(self.dos.densities))

        if analyzer_for_second_spin:
            if not np.all(self.dos.energies == analyzer_for_second_spin.dos.energies):
                raise BoltztrapError("Dos merging error: energies of the two dos are different")

            spin_2 = next(iter(analyzer_for_second_spin.dos.densities))
            if spin_1 == spin_2:
                raise BoltztrapError("Dos merging error: spin component are the same")

        for s in self._dos_partial:
            idx = int(s)
            if structure[idx] not in pdoss:
                pdoss[structure[idx]] = {}
            for o in self._dos_partial[s]:
                if Orbital[o] not in pdoss[structure[idx]]:
                    pdoss[structure[idx]][Orbital[o]] = {}
                pdoss[structure[idx]][Orbital[o]][spin_1] = self._dos_partial[s][o]
                if analyzer_for_second_spin:
                    pdoss[structure[idx]][Orbital[o]][spin_2] = analyzer_for_second_spin._dos_partial[s][o]
        if analyzer_for_second_spin:
            tdos = Dos(
                self.dos.efermi,
                self.dos.energies,
                {
                    spin_1: self.dos.densities[spin_1],
                    spin_2: analyzer_for_second_spin.dos.densities[spin_2],
                },
            )
        else:
            tdos = self.dos

        return CompleteDos(structure, total_dos=tdos, pdoss=pdoss)

    def get_mu_bounds(self, temp=300):
        """:param temp: Temperature.

        Returns:
            The chemical potential bounds at that temperature.
        """
        return min(self.mu_doping["p"][temp]), max(self.mu_doping["n"][temp])

    def get_carrier_concentration(self):
        """Gives the carrier concentration (in cm^-3).

        Returns:
            a dictionary {temp:[]} with an array of carrier concentration
            (in cm^-3) at each temperature
            The array relates to each step of electron chemical potential
        """
        return {temp: [1e24 * i / self.vol for i in self._carrier_conc[temp]] for temp in self._carrier_conc}

    def get_hall_carrier_concentration(self):
        """Gives the Hall carrier concentration (in cm^-3). This is the trace of
        the Hall tensor (see Boltztrap source code) Hall carrier concentration
        are not always exactly the same than carrier concentration.

        Returns:
            a dictionary {temp:[]} with an array of Hall carrier concentration
            (in cm^-3) at each temperature The array relates to each step of
            electron chemical potential
        """
        result = {temp: [] for temp in self._hall}
        for temp in self._hall:
            for i in self._hall[temp]:
                trace = (i[1][2][0] + i[2][0][1] + i[0][1][2]) / 3.0
                if trace != 0.0:
                    result[temp].append(1e-6 / (trace * constants.e))
                else:
                    result[temp].append(0.0)
        return result

    @staticmethod
    def parse_outputtrans(path_dir):
        """Parses .outputtrans file.

        Args:
            path_dir: dir containing boltztrap.outputtrans

        Returns:
            tuple - (run_type, warning, efermi, gap, doping_levels)
        """
        run_type = warning = efermi = gap = None
        doping_levels = []

        with open(f"{path_dir}/boltztrap.outputtrans") as f:
            for line in f:
                if "WARNING" in line:
                    warning = line
                elif "Calc type:" in line:
                    run_type = line.split()[-1]
                elif line.startswith("VBM"):
                    efermi = Energy(line.split()[1], "Ry").to("eV")
                elif line.startswith("Egap:"):
                    gap = Energy(float(line.split()[1]), "Ry").to("eV")
                elif line.startswith("Doping level number"):
                    doping_levels.append(float(line.split()[6]))

        return run_type, warning, efermi, gap, doping_levels

    @staticmethod
    def parse_transdos(path_dir, efermi, dos_spin=1, trim_dos=False):
        """Parses .transdos (total DOS) and .transdos_x_y (partial DOS) files.

        Args:
            path_dir: (str) dir containing DOS files
            efermi: (float) Fermi energy
            dos_spin: (int) -1 for spin down, +1 for spin up
            trim_dos: (bool) whether to post-process / trim DOS.

        Returns:
            tuple - (DOS, dict of partial DOS)
        """
        data_dos = {"total": [], "partial": {}}
        # parse the total DOS data
        # format is energy, DOS, integrated DOS
        with open(f"{path_dir}/boltztrap.transdos") as f:
            count_series = 0  # TODO: why is count_series needed?
            for line in f:
                if line.lstrip().startswith("#"):
                    count_series += 1
                    if count_series > 1:
                        break
                else:
                    data_dos["total"].append(
                        [
                            Energy(float(line.split()[0]), "Ry").to("eV"),
                            float(line.split()[1]),
                        ]
                    )

        lw_l = 0
        hg_l = -len(data_dos["total"])
        if trim_dos:
            # Francesco knows what this does
            # It has something to do with a trick of adding fake energies
            # at the endpoints of the DOS, and then re-trimming it. This is
            # to get the same energy scale for up and down spin DOS.
            tmp_data = np.array(data_dos["total"])
            tmp_den = np.trim_zeros(tmp_data[:, 1], "f")[1:]
            lw_l = len(tmp_data[:, 1]) - len(tmp_den)
            tmp_ene = tmp_data[lw_l:, 0]
            tmp_den = np.trim_zeros(tmp_den, "b")[:-1]
            hg_l = len(tmp_ene) - len(tmp_den)
            tmp_ene = tmp_ene[:-hg_l]
            tmp_data = np.vstack((tmp_ene, tmp_den)).T
            data_dos["total"] = tmp_data.tolist()

        # parse partial DOS data
        for file_name in os.listdir(path_dir):
            if file_name.endswith("transdos") and file_name != "boltztrap.transdos":
                tokens = file_name.split(".")[1].split("_")
                site = tokens[1]
                orb = "_".join(tokens[2:])
                with open(os.path.join(path_dir, file_name)) as f:
                    for line in f:
                        if not line.lstrip().startswith(" #"):
                            if site not in data_dos["partial"]:
                                data_dos["partial"][site] = {}
                            if orb not in data_dos["partial"][site]:
                                data_dos["partial"][site][orb] = []
                            data_dos["partial"][site][orb].append(float(line.split()[1]))
                data_dos["partial"][site][orb] = data_dos["partial"][site][orb][lw_l:-hg_l]

        dos_full = {"energy": [], "density": []}

        for t in data_dos["total"]:
            dos_full["energy"].append(t[0])
            dos_full["density"].append(t[1])

        dos = Dos(efermi, dos_full["energy"], {Spin(dos_spin): dos_full["density"]})
        dos_partial = data_dos["partial"]  # TODO: make this real DOS object?

        return dos, dos_partial

    @staticmethod
    def parse_intrans(path_dir):
        """Parses boltztrap.intrans mainly to extract the value of scissor applied
        to the bands or some other inputs.

        Args:
            path_dir: (str) dir containing the boltztrap.intrans file

        Returns:
            intrans (dict): a dictionary containing various inputs that had
                been used in the Boltztrap run.
        """
        intrans = {}
        with open(f"{path_dir}/boltztrap.intrans") as f:
            for line in f:
                if "iskip" in line:
                    intrans["scissor"] = Energy(float(line.split(" ")[3]), "Ry").to("eV")
                if "HISTO" in line or "TETRA" in line:
                    intrans["dos_type"] = line[:-1]
        return intrans

    @staticmethod
    def parse_struct(path_dir):
        """Parses boltztrap.struct file (only the volume).

        Args:
            path_dir: (str) dir containing the boltztrap.struct file

        Returns:
            (float) volume
        """
        with open(f"{path_dir}/boltztrap.struct") as f:
            tokens = f.readlines()
            return Lattice(
                [[Length(float(tokens[i].split()[j]), "bohr").to("ang") for j in range(3)] for i in range(1, 4)]
            ).volume

    @staticmethod
    def parse_cond_and_hall(path_dir, doping_levels=None):
        """Parses the conductivity and Hall tensors.

        Args:
            path_dir: Path containing .condtens / .halltens files
            doping_levels: ([float]) - doping lvls, parse outtrans to get this.

        Returns:
            mu_steps, cond, seebeck, kappa, hall, pn_doping_levels,
            mu_doping, seebeck_doping, cond_doping, kappa_doping,
            hall_doping, carrier_conc
        """
        # Step 1: parse raw data but do not convert to final format
        t_steps = set()
        mu_steps = set()
        data_full = []
        data_hall = []
        data_doping_full = []
        data_doping_hall = []
        doping_levels = doping_levels or []

        # parse the full conductivity/Seebeck/kappa0/etc data
        # also initialize t_steps and mu_steps
        with open(f"{path_dir}/boltztrap.condtens") as f:
            for line in f:
                if not line.startswith("#"):
                    mu_steps.add(float(line.split()[0]))
                    t_steps.add(int(float(line.split()[1])))
                    data_full.append([float(c) for c in line.split()])

        # parse the full Hall tensor
        with open(f"{path_dir}/boltztrap.halltens") as f:
            for line in f:
                if not line.startswith("#"):
                    data_hall.append([float(c) for c in line.split()])

        if len(doping_levels) != 0:
            # parse doping levels version of full cond. tensor, etc.
            with open(f"{path_dir}/boltztrap.condtens_fixdoping") as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:
                        data_doping_full.append([float(c) for c in line.split()])

            # parse doping levels version of full hall tensor
            with open(f"{path_dir}/boltztrap.halltens_fixdoping") as f:
                for line in f:
                    if not line.startswith("#") and len(line) > 2:
                        data_doping_hall.append([float(c) for c in line.split()])

        # Step 2: convert raw data to final format

        # sort t and mu_steps (b/c they are sets not lists)
        # and convert to correct energy
        t_steps = sorted(t_steps)
        mu_steps = sorted(Energy(m, "Ry").to("eV") for m in mu_steps)

        # initialize output variables - could use defaultdict instead
        # I am leaving things like this for clarity
        cond = {t: [] for t in t_steps}
        seebeck = {t: [] for t in t_steps}
        kappa = {t: [] for t in t_steps}
        hall = {t: [] for t in t_steps}
        carrier_conc = {t: [] for t in t_steps}

        mu_doping = {"p": {t: [] for t in t_steps}, "n": {t: [] for t in t_steps}}
        seebeck_doping = {"p": {t: [] for t in t_steps}, "n": {t: [] for t in t_steps}}
        cond_doping = {"p": {t: [] for t in t_steps}, "n": {t: [] for t in t_steps}}
        kappa_doping = {"p": {t: [] for t in t_steps}, "n": {t: [] for t in t_steps}}
        hall_doping = {"p": {t: [] for t in t_steps}, "n": {t: [] for t in t_steps}}

        # process doping levels
        pn_doping_levels = {"p": [], "n": []}
        for d in doping_levels:
            if d > 0:
                pn_doping_levels["p"].append(d)
            else:
                pn_doping_levels["n"].append(-d)

        # process raw conductivity data, etc.
        for d in data_full:
            temp, doping = d[1], d[2]
            carrier_conc[temp].append(doping)

            cond[temp].append(np.reshape(d[3:12], (3, 3)).tolist())
            seebeck[temp].append(np.reshape(d[12:21], (3, 3)).tolist())
            kappa[temp].append(np.reshape(d[21:30], (3, 3)).tolist())

        # process raw Hall data
        for d in data_hall:
            temp, doping = d[1], d[2]
            hall_tens = [
                np.reshape(d[3:12], (3, 3)).tolist(),
                np.reshape(d[12:21], (3, 3)).tolist(),
                np.reshape(d[21:30], (3, 3)).tolist(),
            ]
            hall[temp].append(hall_tens)

        # process doping conductivity data, etc.
        for d in data_doping_full:
            temp, doping, mu = d[0], d[1], d[-1]
            pn = "p" if doping > 0 else "n"
            mu_doping[pn][temp].append(Energy(mu, "Ry").to("eV"))
            cond_doping[pn][temp].append(np.reshape(d[2:11], (3, 3)).tolist())
            seebeck_doping[pn][temp].append(np.reshape(d[11:20], (3, 3)).tolist())
            kappa_doping[pn][temp].append(np.reshape(d[20:29], (3, 3)).tolist())

        # process doping Hall data
        for d in data_doping_hall:
            temp, doping, mu = d[0], d[1], d[-1]
            pn = "p" if doping > 0 else "n"
            hall_tens = [
                np.reshape(d[2:11], (3, 3)).tolist(),
                np.reshape(d[11:20], (3, 3)).tolist(),
                np.reshape(d[20:29], (3, 3)).tolist(),
            ]
            hall_doping[pn][temp].append(hall_tens)

        return (
            mu_steps,
            cond,
            seebeck,
            kappa,
            hall,
            pn_doping_levels,
            mu_doping,
            seebeck_doping,
            cond_doping,
            kappa_doping,
            hall_doping,
            carrier_conc,
        )

    @staticmethod
    def from_files(path_dir, dos_spin=1):
        """Get a BoltztrapAnalyzer object from a set of files.

        Args:
            path_dir: directory where the boltztrap files are
            dos_spin: in DOS mode, set to 1 for spin up and -1 for spin down

        Returns:
            a BoltztrapAnalyzer object
        """
        run_type, warning, efermi, gap, doping_levels = BoltztrapAnalyzer.parse_outputtrans(path_dir)

        vol = BoltztrapAnalyzer.parse_struct(path_dir)

        intrans = BoltztrapAnalyzer.parse_intrans(path_dir)

        if run_type == "BOLTZ":
            dos, pdos = BoltztrapAnalyzer.parse_transdos(path_dir, efermi, dos_spin=dos_spin, trim_dos=False)

            *cond_and_hall, carrier_conc = BoltztrapAnalyzer.parse_cond_and_hall(path_dir, doping_levels)

            return BoltztrapAnalyzer(gap, *cond_and_hall, intrans, dos, pdos, carrier_conc, vol, warning)

        if run_type == "DOS":
            trim = intrans["dos_type"] == "HISTO"
            dos, pdos = BoltztrapAnalyzer.parse_transdos(path_dir, efermi, dos_spin=dos_spin, trim_dos=trim)

            return BoltztrapAnalyzer(gap=gap, dos=dos, dos_partial=pdos, warning=warning, vol=vol)

        if run_type == "BANDS":
            bz_kpoints = np.loadtxt(f"{path_dir}/boltztrap_band.dat")[:, -3:]
            bz_bands = np.loadtxt(f"{path_dir}/boltztrap_band.dat")[:, 1:-6]
            return BoltztrapAnalyzer(bz_bands=bz_bands, bz_kpoints=bz_kpoints, warning=warning, vol=vol)

        if run_type == "FERMI":
            if os.path.exists(f"{path_dir}/boltztrap_BZ.cube"):
                fs_data = read_cube_file(f"{path_dir}/boltztrap_BZ.cube")
            elif os.path.exists(f"{path_dir}/fort.30"):
                fs_data = read_cube_file(f"{path_dir}/fort.30")
            else:
                raise BoltztrapError("No data file found for fermi surface")
            return BoltztrapAnalyzer(fermi_surface_data=fs_data)

        raise ValueError(f"{run_type=} not recognized!")

    def as_dict(self):
        """MSONable dict."""
        results = {
            "gap": self.gap,
            "mu_steps": self.mu_steps,
            "intrans": self.intrans,
            "cond": self._cond,
            "seebeck": self._seebeck,
            "kappa": self._kappa,
            "hall": self._hall,
            "doping": self.doping,
            "mu_doping": self.mu_doping,
            "seebeck_doping": self._seebeck_doping,
            "cond_doping": self._cond_doping,
            "kappa_doping": self._kappa_doping,
            "hall_doping": self._hall_doping,
            "dos": self.dos.as_dict(),
            "dos_partial": self._dos_partial,
            "carrier_conc": self._carrier_conc,
            "vol": self.vol,
            "warning": self.warning,
        }
        return jsanitize(results)

    @staticmethod
    def from_dict(data):
        """:param data: Dict representation.

        Returns:
            BoltztrapAnalyzer
        """

        def _make_float_array(a):
            res = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            for i in range(3):
                for j in range(3):
                    res[i][j] = float(a[i][j])
            return res

        def _make_float_hall(a):
            return list(a[:27])

        gap = data.get("gap")
        mu_steps = [float(d) for d in data["mu_steps"]] if "mu_steps" in data else None
        cond = (
            {int(d): [_make_float_array(v) for v in data["cond"][d]] for d in data["cond"]} if "cond" in data else None
        )
        seebeck = (
            {int(d): [_make_float_array(v) for v in data["seebeck"][d]] for d in data["seebeck"]}
            if "seebeck" in data
            else None
        )
        kappa = (
            {int(d): [_make_float_array(v) for v in data["kappa"][d]] for d in data["kappa"]}
            if "kappa" in data
            else None
        )
        hall = (
            {int(d): [_make_float_hall(v) for v in data["hall"][d]] for d in data["hall"]} if "hall" in data else None
        )
        doping = (
            {
                "p": [float(d) for d in data["doping"]["p"]],
                "n": [float(d) for d in data["doping"]["n"]],
            }
            if "doping" in data
            else None
        )

        mu_doping = (
            {
                "p": {int(d): [float(v) for v in data["mu_doping"]["p"][d]] for d in data["mu_doping"]["p"]},
                "n": {int(d): [float(v) for v in data["mu_doping"]["n"][d]] for d in data["mu_doping"]["n"]},
            }
            if "mu_doping" in data
            else None
        )

        seebeck_doping = (
            {
                "p": {
                    int(d): [_make_float_array(v) for v in data["seebeck_doping"]["p"][d]]
                    for d in data["seebeck_doping"]["p"]
                },
                "n": {
                    int(d): [_make_float_array(v) for v in data["seebeck_doping"]["n"][d]]
                    for d in data["seebeck_doping"]["n"]
                },
            }
            if "seebeck_doping" in data
            else None
        )

        cond_doping = (
            {
                "p": {
                    int(d): [_make_float_array(v) for v in data["cond_doping"]["p"][d]]
                    for d in data["cond_doping"]["p"]
                },
                "n": {
                    int(d): [_make_float_array(v) for v in data["cond_doping"]["n"][d]]
                    for d in data["cond_doping"]["n"]
                },
            }
            if "cond_doping" in data
            else None
        )

        kappa_doping = (
            {
                "p": {
                    int(d): [_make_float_array(v) for v in data["kappa_doping"]["p"][d]]
                    for d in data["kappa_doping"]["p"]
                },
                "n": {
                    int(d): [_make_float_array(v) for v in data["kappa_doping"]["n"][d]]
                    for d in data["kappa_doping"]["n"]
                },
            }
            if "kappa_doping" in data
            else None
        )

        hall_doping = (
            {
                "p": {
                    int(d): [_make_float_hall(v) for v in data["hall_doping"]["p"][d]] for d in data["hall_doping"]["p"]
                },
                "n": {
                    int(d): [_make_float_hall(v) for v in data["hall_doping"]["n"][d]] for d in data["hall_doping"]["n"]
                },
            }
            if "hall_doping" in data
            else None
        )

        dos = Dos.from_dict(data["dos"]) if "dos" in data else None
        dos_partial = data.get("dos_partial")
        carrier_conc = data.get("carrier_conc")
        vol = data.get("vol")
        warning = data.get("warning")

        return BoltztrapAnalyzer(
            gap=gap,
            mu_steps=mu_steps,
            cond=cond,
            seebeck=seebeck,
            kappa=kappa,
            hall=hall,
            doping=doping,
            mu_doping=mu_doping,
            seebeck_doping=seebeck_doping,
            cond_doping=cond_doping,
            kappa_doping=kappa_doping,
            hall_doping=hall_doping,
            dos=dos,
            dos_partial=dos_partial,
            carrier_conc=carrier_conc,
            vol=vol,
            warning=warning,
        )


def read_cube_file(filename):
    """:param filename: Cube filename

    Returns:
        Energy data.
    """
    with open(filename) as f:
        natoms = 0
        count_line = 0
        for line in f:
            line = line.rstrip("\n")
            if count_line == 0 and "CUBE" not in line:
                raise ValueError("CUBE file format not recognized")

            if count_line == 2:
                tokens = line.split()
                natoms = int(tokens[0])
            if count_line == 3:
                tokens = line.split()
                n1 = int(tokens[0])
            elif count_line == 4:
                tokens = line.split()
                n2 = int(tokens[0])
            elif count_line == 5:
                tokens = line.split()
                n3 = int(tokens[0])
            elif count_line > 5:
                break

            count_line += 1

    if "fort.30" in filename:
        energy_data = np.genfromtxt(filename, skip_header=natoms + 6, skip_footer=1)
        nlines_data = len(energy_data)
        last_line = np.genfromtxt(filename, skip_header=nlines_data + natoms + 6)
        energy_data = np.append(energy_data.flatten(), last_line).reshape(n1, n2, n3)  # pylint: disable=E1121
    elif "boltztrap_BZ.cube" in filename:
        energy_data = np.loadtxt(filename, skiprows=natoms + 6).reshape(n1, n2, n3)

    energy_data /= Energy(1, "eV").to("Ry")

    return energy_data


def compare_sym_bands(bands_obj, bands_ref_obj, nb=None):
    """Compute the mean of correlation between bzt and vasp bandstructure on
    sym line, for all bands and locally (for each branches) the difference
    squared (%) if nb is specified.
    """
    if bands_ref_obj.is_spin_polarized:
        nbands = min(bands_obj.nb_bands, 2 * bands_ref_obj.nb_bands)
    else:
        # TODO: why is this needed? Shouldn't pmg take care of nb_bands?
        nbands = min(len(bands_obj.bands[Spin.up]), len(bands_ref_obj.bands[Spin.up]))
    arr_bands = np.array(bands_obj.bands[Spin.up][:nbands])
    # arr_bands_lavg = (arr_bands-np.mean(arr_bands,axis=1).reshape(nbands,1))

    if bands_ref_obj.is_spin_polarized:
        arr_bands_ref_up = np.array(bands_ref_obj.bands[Spin.up])
        arr_bands_ref_dw = np.array(bands_ref_obj.bands[Spin.down])
        arr_bands_ref = np.vstack((arr_bands_ref_up, arr_bands_ref_dw))
        arr_bands_ref = np.sort(arr_bands_ref, axis=0)[:nbands]
    else:
        arr_bands_ref = np.array(bands_ref_obj.bands[Spin.up][:nbands])

    # arr_bands_ref_lavg =
    # (arr_bands_ref-np.mean(arr_bands_ref,axis=1).reshape(nbands,1))

    # err = np.sum((arr_bands_lavg-arr_bands_ref_lavg)**2,axis=1)/nkpt
    corr = np.array([distance.correlation(arr_bands[idx], arr_bands_ref[idx]) for idx in range(nbands)])

    if isinstance(nb, int):
        nb = [nb]

    bcheck = {}

    if max(nb) < nbands:
        branches = [[s["start_index"], s["end_index"], s["name"]] for s in bands_ref_obj.branches]

        if not bands_obj.is_metal() and not bands_ref_obj.is_metal():
            zero_ref = bands_ref_obj.get_vbm()["energy"]
            zero = bands_obj.get_vbm()["energy"]
            if not zero:
                vbm = bands_ref_obj.get_vbm()["band_index"][Spin.up][-1]
                zero = max(arr_bands[vbm])
        else:
            zero_ref = 0  # bands_ref_obj.efermi
            zero = 0  # bands_obj.efermi
            print(zero, zero_ref)

        for nbi in nb:
            bcheck[nbi] = {}

            bcheck[nbi]["Dist"] = np.mean(abs(arr_bands[nbi] - zero - arr_bands_ref[nbi] + zero_ref))
            bcheck[nbi]["Corr"] = corr[nbi]

            for start, end, name in branches:
                # werr.append((sum((arr_bands_corr[nb][start:end+1] -
                # arr_bands_ref_corr[nb][start:end+1])**2)/(end+1-start)*100,name))
                bcheck[nbi][name] = np.mean(
                    abs(arr_bands[nbi][start : end + 1] - zero - arr_bands_ref[nbi][start : end + 1] + zero_ref)
                )
    else:
        bcheck = "No nb given"

    return bcheck


def seebeck_spb(eta, Lambda=0.5):
    """Seebeck analytic formula in the single parabolic model."""
    try:
        from fdint import fdk
    except ImportError:
        raise BoltztrapError(
            "fdint module not found. Please, install it.\nIt is needed to calculate Fermi integral quickly."
        )

    return (
        constants.k
        / constants.e
        * ((2.0 + Lambda) * fdk(1.0 + Lambda, eta) / ((1.0 + Lambda) * fdk(Lambda, eta)) - eta)
        * 1e6
    )


def eta_from_seebeck(seeb, Lambda):
    """It takes a value of seebeck and adjusts the analytic seebeck until it's equal.

    Returns:
        float: eta where the two seebeck coefficients are equal (reduced chemical potential).
    """
    from scipy.optimize import fsolve

    out = fsolve(lambda x: (seebeck_spb(x, Lambda) - abs(seeb)) ** 2, 1.0, full_output=True)
    return out[0][0]


def seebeck_eff_mass_from_carr(eta, n, T, Lambda):
    """Calculate seebeck effective mass at a certain carrier concentration
    eta in kB*T units, n in cm-3, T in K, returns mass in m0 units.
    """
    try:
        from fdint import fdk
    except ImportError:
        raise BoltztrapError(
            "fdint module not found. Please, install it.\nIt is needed to calculate Fermi integral quickly."
        )

    return (2 * np.pi**2 * abs(n) * 10**6 / (fdk(0.5, eta))) ** (2.0 / 3) / (
        2 * constants.m_e * constants.k * T / (constants.h / 2 / np.pi) ** 2
    )


def seebeck_eff_mass_from_seebeck_carr(seeb, n, T, Lambda):
    """Find the chemical potential where analytic and calculated seebeck are identical
    and then calculate the seebeck effective mass at that chemical potential and
    a certain carrier concentration n.
    """
    eta = eta_from_seebeck(seeb, Lambda)
    return seebeck_eff_mass_from_carr(eta, n, T, Lambda)
