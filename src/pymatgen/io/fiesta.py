"""
This module implements input and output for Fiesta (http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html).

- Nwchem2Fiesta: create the input files needed for a Fiesta run
- FiestaRun: run gw_fiesta and bse_fiesta
- Localized BasisSetReader
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from string import Template
from typing import TYPE_CHECKING

from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.structure import Molecule

if TYPE_CHECKING:
    from pathlib import Path

    from typing_extensions import Self


__author__ = "ndardenne"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__email__ = "n.dardenne@uclouvain.be"
__date__ = "24/5/15"


class Nwchem2Fiesta(MSONable):
    """To run NWCHEM2FIESTA inside python:

    If nwchem.nw is the input, nwchem.out the output, and structure.movecs the
    "movecs" file, the syntax to run NWCHEM2FIESTA is: NWCHEM2FIESTA
    nwchem.nw nwchem.nwout structure.movecs > log_n2f
    """

    def __init__(self, folder, filename="nwchem", log_file="log_n2f"):
        """
        Args:
            folder: where are stored the nwchem
            filename: name of nwchem files read by NWCHEM2FIESTA (filename.nw, filename.nwout and filename.movecs)
            logfile: logfile of NWCHEM2FIESTA.
        """
        self.folder = folder
        self.filename = filename
        self.log_file = log_file

        self._NWCHEM2FIESTA_cmd = "NWCHEM2FIESTA"
        self._nwcheminput_fn = f"{filename}.nw"
        self._nwchemoutput_fn = f"{filename}.nwout"
        self._nwchemmovecs_fn = f"{filename}.movecs"

    def run(self):
        """Perform actual NWCHEM2FIESTA run."""
        init_folder = os.getcwd()
        os.chdir(self.folder)

        with zopen(self.log_file, mode="wt", encoding="utf-8") as fout:
            subprocess.call(
                [
                    self._NWCHEM2FIESTA_cmd,
                    self._nwcheminput_fn,
                    self._nwchemoutput_fn,
                    self._nwchemmovecs_fn,
                ],
                stdout=fout,
            )

        os.chdir(init_folder)

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "filename": self.filename,
            "folder": self.folder,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            Nwchem2Fiesta
        """
        return cls(folder=dct["folder"], filename=dct["filename"])


class FiestaRun(MSONable):
    """
    To run FIESTA inside python:
        if grid is [x,x] then bse runs
        if grid is [x,x,y] the fiesta(gw) runs
        otherwise it breaks.
    """

    def __init__(
        self,
        folder: str | None = None,
        grid: tuple[int, int, int] = (2, 2, 2),
        log_file: str = "log",
    ) -> None:
        """
        Args:
            folder: Folder to look for runs.
            grid:
            log_file: logfile of Fiesta.
        """
        self.folder = folder or os.getcwd()
        self.log_file = log_file
        self.grid = grid

    def run(self):
        """Perform FIESTA (gw) run."""
        if len(self.grid) == 3:
            self.mpi_procs = self.grid[0] * self.grid[1] * self.grid[2]
            self._gw_run()
        elif len(self.grid) == 2:
            self.mpi_procs = self.grid[0] * self.grid[1]
            self.bse_run()
        else:
            raise ValueError("Wrong grid size: must be [nrow, ncolumn, nslice] for gw of [nrow, nslice] for bse")

    def _gw_run(self):
        """Perform FIESTA (gw) run."""
        init_folder = os.getcwd()
        if self.folder != init_folder:
            os.chdir(self.folder)

        with zopen(self.log_file, mode="wt", encoding="utf-8") as fout:
            subprocess.call(
                [
                    "mpirun",
                    "-n",
                    str(self.mpi_procs),
                    "fiesta",
                    str(self.grid[0]),
                    str(self.grid[1]),
                    str(self.grid[2]),
                ],
                stdout=fout,
            )

        if self.folder != init_folder:
            os.chdir(init_folder)

    def bse_run(self):
        """Perform BSE run."""
        init_folder = os.getcwd()
        if self.folder != init_folder:
            os.chdir(self.folder)

        with zopen(self.log_file, mode="wt", encoding="utf-8") as fout:
            subprocess.call(
                [
                    "mpirun",
                    "-n",
                    str(self.mpi_procs),
                    "bse",
                    str(self.grid[0]),
                    str(self.grid[1]),
                ],
                stdout=fout,
            )

        if self.folder != init_folder:
            os.chdir(init_folder)

    def as_dict(self):
        """MSONable dict."""
        return {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "log_file": self.log_file,
            "grid": self.grid,
            "folder": self.folder,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            FiestaRun
        """
        return cls(folder=dct["folder"], grid=dct["grid"], log_file=dct["log_file"])


class BasisSetReader:
    """
    A basis set reader.
    Basis set are stored in data as a dict:
    :key l_zeta_ng for each nl orbitals which contain list of tuple (alpha, coef) for each of the ng gaussians
    in l_zeta orbital.
    """

    def __init__(self, filename):
        """
        Args:
            filename: Filename to read.
        """
        self.filename = filename

        with zopen(filename, mode="rt", encoding="utf-8") as file:
            basis_set = file.read()

        self.data = self._parse_file(basis_set)
        # compute the number of nlm orbitals per atom
        self.data.update(n_nlmo=self.set_n_nlmo())

    @staticmethod
    def _parse_file(lines):
        lmax_nnlo_patt = re.compile(r"\s* (\d+) \s+ (\d+) \s+ \# .* ", re.VERBOSE)

        nl_orbital_patt = re.compile(r"\s* (\d+) \s+ (\d+) \s+ (\d+) \s+ \# .* ", re.VERBOSE)

        coef_alpha_patt = re.compile(r"\s* ([-\d.\D]+) \s+ ([-\d.\D]+) \s* ", re.VERBOSE)

        preamble = []
        basis_set = {}
        parse_preamble = parse_lmax_nnlo = False
        parse_nl_orbital = False
        nnlo = lmax = None
        l_angular = zeta = ng = None

        for line in lines.split("\n"):
            if parse_nl_orbital:
                match_orb = nl_orbital_patt.search(line)
                match_alpha = coef_alpha_patt.search(line)
                if match_orb:
                    l_angular = match_orb[1]
                    zeta = match_orb[2]
                    ng = match_orb[3]
                    basis_set[f"{l_angular}_{zeta}_{ng}"] = []
                elif match_alpha:
                    alpha = match_alpha[1]
                    coef = match_alpha[2]
                    basis_set[f"{l_angular}_{zeta}_{ng}"].append((alpha, coef))
            elif parse_lmax_nnlo:
                if match_orb := lmax_nnlo_patt.search(line):
                    lmax = match_orb[1]
                    nnlo = match_orb[2]
                    parse_lmax_nnlo = False
                    parse_nl_orbital = True
            elif parse_preamble:
                preamble.append(line.strip())

            if line.find("</preamble>") != -1:
                parse_preamble = False
                parse_lmax_nnlo = True
            elif line.find("<preamble>") != -1:
                parse_preamble = True

        basis_set.update(lmax=lmax, n_nlo=nnlo, preamble=preamble)
        return basis_set

    def set_n_nlmo(self):
        """The number of nlm orbitals for the basis set."""
        n_nlm_orbs = 0

        data_tmp = self.data
        data_tmp.pop("lmax")
        data_tmp.pop("n_nlo")
        data_tmp.pop("preamble")

        for l_zeta_ng in data_tmp:
            n_l = l_zeta_ng.split("_")[0]
            n_nlm_orbs += 2 * int(n_l) + 1

        return str(n_nlm_orbs)

    def infos_on_basis_set(self):
        return (
            f"=========================================\n"
            f"Reading basis set:\n"
            f"\n"
            f"Basis set for {self.filename} atom \n"
            f"Maximum angular momentum = {self.data['lmax']}\n"
            f"Number of atomics orbitals = {self.data['n_nlo']}\n"
            f"Number of nlm orbitals = {self.data['n_nlmo']}\n"
            f"========================================="
        )


class FiestaInput(MSONable):
    """Input File for Fiesta called "cell.in" by default (mandatory in Fiesta for now)."""

    def __init__(
        self,
        mol,
        correlation_grid: dict[str, str] | None = None,
        exc_dft_option: dict[str, str] | None = None,
        cohsex_options: dict[str, str] | None = None,
        gw_options: dict[str, str] | None = None,
        bse_tddft_options: dict[str, str] | None = None,
    ):
        """
        Args:
            mol: pymatgen mol
            correlation_grid: dict
            Exc_DFT_option: dict
            COHSEX_options: dict
            GW_options: dict
            BSE_TDDFT_options: dict.
        """
        self._mol = mol
        self.correlation_grid = correlation_grid or {"dE_grid": "0.500", "n_grid": "14"}
        self.Exc_DFT_option = exc_dft_option or {"rdVxcpsi": "1"}
        self.cohsex_options = cohsex_options or {
            "eigMethod": "C",
            "mix_cohsex": "0.500",
            "nc_cohsex": "0",
            "nit_cohsex": "0",
            "nv_cohsex": "0",
            "resMethod": "V",
            "scf_cohsex_wf": "0",
        }
        self.GW_options = gw_options or {
            "nc_corr": "10",
            "nit_gw": "3",
            "nv_corr": "10",
        }
        self.bse_tddft_options = bse_tddft_options or {
            "do_bse": "1",
            "do_tddft": "0",
            "nc_bse": "382",
            "nit_bse": "50",
            "npsi_bse": "1",
            "nv_bse": "21",
        }

    def set_auxiliary_basis_set(self, folder, auxiliary_folder, auxiliary_basis_set_type="aug_cc_pvtz"):
        """
        copy in the desired folder the needed auxiliary basis set "X2.ion" where X is a specie.

        Args:
            auxiliary_folder: folder where the auxiliary basis sets are stored
            auxiliary_basis_set_type: type of basis set (string to be found in the extension of the file name; must
                be in lower case). ex: C2.ion_aug_cc_pvtz_RI_Weigend find "aug_cc_pvtz".
        """
        list_files = os.listdir(auxiliary_folder)

        for specie in self._mol.symbol_set:
            for file in list_files:
                if file.upper().find(f"{specie.upper()}2") != -1 and file.lower().find(auxiliary_basis_set_type) != -1:
                    shutil.copyfile(f"{auxiliary_folder}/{file}", f"{folder}/{specie}2.ion")

    def set_gw_options(self, nv_band=10, nc_band=10, n_iteration=5, n_grid=6, dE_grid=0.5):
        """Set parameters in cell.in for a GW computation.

        Args:
            nv__band: number of valence bands to correct with GW
            nc_band: number of conduction bands to correct with GW
            n_iteration: number of iteration
            n_grid and dE_grid: number of points and spacing in eV for correlation grid.
        """
        self.GW_options.update(nv_corr=nv_band, nc_corr=nc_band, nit_gw=n_iteration)
        self.correlation_grid.update(dE_grid=dE_grid, n_grid=n_grid)

    @staticmethod
    def make_full_bse_densities_folder(folder):
        """Mkdir "FULL_BSE_Densities" folder (needed for bse run) in the desired folder."""
        if os.path.isfile(f"{folder}/FULL_BSE_Densities"):
            return "FULL_BSE_Densities folder already exists"

        os.makedirs(f"{folder}/FULL_BSE_Densities")
        return "makedirs FULL_BSE_Densities folder"

    def set_bse_options(self, n_excitations=10, nit_bse=200):
        """Set parameters in cell.in for a BSE computation.

        Args:
            nv_bse: number of valence bands
            nc_bse: number of conduction bands
            n_excitations: number of excitations
            nit_bse: number of iterations.
        """
        self.bse_tddft_options.update(npsi_bse=n_excitations, nit_bse=nit_bse)

    def dump_bse_data_in_gw_run(self, BSE_dump=True):
        """Set the "do_bse" variable to 1 or 0 in cell.in.

        Args:
            BSE_dump (bool): Defaults to True.
        """
        self.bse_tddft_options.update(do_bse=int(BSE_dump), do_tddft=0)

    def dump_tddft_data_in_gw_run(self, tddft_dump: bool = True):
        """Set the do_tddft variable to 1 or 0 in cell.in.

        Args:
            tddft_dump (bool): Defaults to True.
        """
        self.bse_tddft_options.update(do_bse="0", do_tddft="1" if tddft_dump else "0")

    @property
    def infos_on_system(self):
        """Infos on initial parameters as in the log file of Fiesta."""
        lst = [
            "=========================================",
            "Reading infos on system:",
            "",
            f" Number of atoms = {self._mol.composition.num_atoms} ; number of species = {len(self._mol.symbol_set)}",
            f" Number of valence bands = {int(self._mol.nelectrons / 2)}",
            f" Sigma grid specs: n_grid = {self.correlation_grid['n_grid']} ;  "
            f"dE_grid = {self.correlation_grid['dE_grid']} (eV)",
        ]
        if int(self.Exc_DFT_option["rdVxcpsi"]) == 1:
            lst.append(" Exchange and correlation energy read from Vxcpsi.mat")
        elif int(self.Exc_DFT_option["rdVxcpsi"]) == 0:
            lst.append(" Exchange and correlation energy re-computed")

        if self.cohsex_options["eigMethod"] == "C":
            lst += [
                f" Correcting  {self.cohsex_options['nv_cohsex']} valence bands and  "
                f"{self.cohsex_options['nc_cohsex']} conduction bands at COHSEX level",
                f" Performing   {self.cohsex_options['nit_cohsex']} diagonal COHSEX iterations",
            ]
        elif self.cohsex_options["eigMethod"] == "HF":
            lst += [
                f" Correcting  {self.cohsex_options['nv_cohsex']} valence bands and  "
                f"{self.cohsex_options['nc_cohsex']} conduction bands at HF level",
                f" Performing   {self.cohsex_options['nit_cohsex']} diagonal HF iterations",
            ]

        lst += [
            f" Using resolution of identity : {self.cohsex_options['resMethod']}",
            f" Correcting  {self.GW_options['nv_corr']} valence bands and "
            f"{self.GW_options['nc_corr']} conduction bands at GW level",
            f" Performing   {self.GW_options['nit_gw']} GW iterations",
        ]

        if int(self.bse_tddft_options["do_bse"]) == 1:
            lst.append(" Dumping data for BSE treatment")

        if int(self.bse_tddft_options["do_tddft"]) == 1:
            lst.append(" Dumping data for TD-DFT treatment")
        lst.extend(("", " Atoms in cell cartesian A:"))
        symbols = list(self._mol.symbol_set)

        for site in self._mol:
            lst.append(f" {site.x} {site.y} {site.z} {int(symbols.index(site.specie.symbol)) + 1}")

        lst.append("=========================================")

        return str(lst)

    @property
    def molecule(self):
        """Molecule associated with this FiestaInput."""
        return self._mol

    def __str__(self):
        symbols = list(self._mol.symbol_set)

        geometry = [f" {site.x} {site.y} {site.z} {symbols.index(site.specie.symbol) + 1}" for site in self._mol]

        t = Template(
            """# number of atoms and species
   $nat    $nsp
# number of valence bands
    $nvbands
# number of points and spacing in eV for correlation grid
    $n_grid    $dE_grid
# relire=1 ou recalculer=0 Exc DFT
    $rdVxcpsi
# number of COHSEX corrected occp and unoccp bands: C=COHSEX  H=HF
    $nv_cohsex    $nc_cohsex   $eigMethod
# number of COHSEX iter, scf on wfns, mixing coeff; V=RI-V  I=RI-D
    $nit_cohsex   $resMethod       $scf_cohsex_wf       $mix_cohsex
# number of GW corrected occp and unoccp bands
   $nv_corr   $nc_corr
# number of GW iterations
    $nit_gw
# dumping for BSE and TDDFT
    $do_bse    $do_tddft
# number of occp. and virtual bands of BSE: nocore and up to 40 eVs
    $nv_bse   $nc_bse
# number of excitations needed and number of iterations
    $npsi_bse   $nit_bse
# list of symbols in order
$symbols
# scaling factor
    1.000
# atoms x,y,z cartesian .. will be multiplied by scale
$geometry
            """
        )

        return t.substitute(
            nat=int(self._mol.composition.num_atoms),
            nsp=len(self._mol.symbol_set),
            nvbands=int(self._mol.nelectrons / 2),
            n_grid=self.correlation_grid["n_grid"],
            dE_grid=self.correlation_grid["dE_grid"],
            rdVxcpsi=self.Exc_DFT_option["rdVxcpsi"],
            nv_cohsex=self.cohsex_options["nv_cohsex"],
            nc_cohsex=self.cohsex_options["nc_cohsex"],
            eigMethod=self.cohsex_options["eigMethod"],
            nit_cohsex=self.cohsex_options["nit_cohsex"],
            resMethod=self.cohsex_options["resMethod"],
            scf_cohsex_wf=self.cohsex_options["scf_cohsex_wf"],
            mix_cohsex=self.cohsex_options["mix_cohsex"],
            nv_corr=self.GW_options["nv_corr"],
            nc_corr=self.GW_options["nc_corr"],
            nit_gw=self.GW_options["nit_gw"],
            do_bse=self.bse_tddft_options["do_bse"],
            do_tddft=self.bse_tddft_options["do_tddft"],
            nv_bse=self.bse_tddft_options["nv_bse"],
            nc_bse=self.bse_tddft_options["nc_bse"],
            npsi_bse=self.bse_tddft_options["npsi_bse"],
            nit_bse=self.bse_tddft_options["nit_bse"],
            symbols="\n".join(symbols),
            geometry="\n".join(geometry),
        )

    def write_file(self, filename: str | Path) -> None:
        """Write FiestaInput to a file.

        Args:
            filename: Filename.
        """
        with zopen(filename, mode="wt", encoding="utf-8") as file:
            file.write(str(self))  # type:ignore[arg-type]

    def as_dict(self):
        """MSONable dict."""
        return {
            "mol": self._mol.as_dict(),
            "correlation_grid": self.correlation_grid,
            "Exc_DFT_option": self.Exc_DFT_option,
            "COHSEX_options": self.cohsex_options,
            "GW_options": self.GW_options,
            "BSE_TDDFT_options": self.bse_tddft_options,
        }

    @classmethod
    def from_dict(cls, dct: dict) -> Self:
        """
        Args:
            dct (dict): Dict representation.

        Returns:
            FiestaInput
        """
        return cls(
            mol=Molecule.from_dict(dct["mol"]),
            correlation_grid=dct["correlation_grid"],
            exc_dft_option=dct["Exc_DFT_option"],
            cohsex_options=dct["geometry_options"],
            gw_options=dct["symmetry_options"],
            bse_tddft_options=dct["memory_options"],
        )

    @classmethod
    def from_str(cls, string_input: str) -> FiestaInput:
        """
        Read an FiestaInput from a string. Currently tested to work with
        files generated from this class itself.

        Args:
            string_input: string_input to parse.

        Returns:
            FiestaInput object
        """
        correlation_grid = {}
        Exc_DFT_option = {}
        COHSEX_options = {}
        GW_options = {}
        BSE_TDDFT_options = {}

        lines = string_input.strip().split("\n")

        # number of atoms and species
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        nat = tokens[0]
        nsp = tokens[1]
        # number of valence bands
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()

        # correlation_grid
        # number of points and spacing in eV for correlation grid
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        correlation_grid["n_grid"] = tokens[0]
        correlation_grid["dE_grid"] = tokens[1]

        # Exc DFT
        # relire=1 ou recalculer=0 Exc DFT
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        Exc_DFT_option["rdVxcpsi"] = tokens[0]

        # COHSEX
        # number of COHSEX corrected occp and unoccp bands: C=COHSEX  H=HF
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        COHSEX_options["nv_cohsex"] = tokens[0]
        COHSEX_options["nc_cohsex"] = tokens[1]
        COHSEX_options["eigMethod"] = tokens[2]
        # number of COHSEX iter, scf on wfns, mixing coeff; V=RI-V  I=RI-D
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        COHSEX_options["nit_cohsex"] = tokens[0]
        COHSEX_options["resMethod"] = tokens[1]
        COHSEX_options["scf_cohsex_wf"] = tokens[2]
        COHSEX_options["mix_cohsex"] = tokens[3]

        # GW
        # number of GW corrected occp and unoccp bands
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        GW_options["nv_corr"] = tokens[0]
        GW_options["nc_corr"] = tokens[1]
        # number of GW iterations
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        GW_options["nit_gw"] = tokens[0]

        # BSE
        # dumping for BSE and TDDFT
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        BSE_TDDFT_options["do_bse"] = tokens[0]
        BSE_TDDFT_options["do_tddft"] = tokens[1]
        # number of occp. and virtual bands of BSE: nocore and up to 40 eVs
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        BSE_TDDFT_options["nv_bse"] = tokens[0]
        BSE_TDDFT_options["nc_bse"] = tokens[1]
        # number of excitations needed and number of iterations
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        BSE_TDDFT_options["npsi_bse"] = tokens[0]
        BSE_TDDFT_options["nit_bse"] = tokens[1]

        # Molecule
        # list of symbols in order
        lines.pop(0)
        atname = []
        i = int(nsp)
        while i != 0:
            line = lines.pop(0).strip()
            tokens = line.split()
            atname.append(tokens[0])
            i -= 1

        # scaling factor
        lines.pop(0)
        line = lines.pop(0).strip()
        tokens = line.split()
        # atoms x,y,z cartesian .. will be multiplied by scale
        lines.pop(0)
        # Parse geometry
        species = []
        coords = []
        i = int(nat)
        while i != 0:
            line = lines.pop(0).strip()
            tokens = line.split()
            coords.append([float(j) for j in tokens[:3]])
            species.append(atname[int(tokens[3]) - 1])
            i -= 1

        mol = Molecule(species, coords)

        return FiestaInput(
            mol=mol,
            correlation_grid=correlation_grid,
            exc_dft_option=Exc_DFT_option,
            cohsex_options=COHSEX_options,
            gw_options=GW_options,
            bse_tddft_options=BSE_TDDFT_options,
        )

    @classmethod
    def from_file(cls, filename: str | Path) -> FiestaInput:
        """
        Read an Fiesta input from a file. Currently tested to work with
        files generated from this class itself.

        Args:
            filename: Filename to parse.

        Returns:
            FiestaInput object
        """
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            return cls.from_str(file.read())  # type:ignore[arg-type]


class FiestaOutput:
    """
    A Fiesta output file parser.

    All energies are in eV.
    """

    def __init__(self, filename):
        """
        Args:
            filename: Filename to read.
        """
        self.filename = filename

        with zopen(filename, mode="rt", encoding="utf-8") as file:
            data = file.read()

        chunks = re.split(r"GW Driver iteration", data)

        # preamble: everything before the first GW Driver iteration
        chunks.pop(0)

        # self.job_info = self._parse_preamble(preamble)
        self.data = [self._parse_job(c) for c in chunks]

    @staticmethod
    def _parse_job(output):
        GW_BANDS_results_patt = re.compile(
            r"^<it.*  \| \s+ (\D+\d*) \s+ \| \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ \| "
            r" \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ \|"
            r" \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ",
            re.VERBOSE,
        )

        GW_GAPS_results_patt = re.compile(
            r"^<it.*  \| \s+ Egap_KS \s+ = \s+ ([-\d.]+) \s+ \| \s+ Egap_QP \s+ = \s+ ([-\d.]+) \s+ \| "
            r" \s+ Egap_QP \s+ = \s+ ([-\d.]+) \s+ \|",
            re.VERBOSE,
        )

        end_patt = re.compile(r"\s*program returned normally\s*")

        total_time_patt = re.compile(r"\s*total \s+ time: \s+  ([\d.]+) .*", re.VERBOSE)

        GW_results = {}
        parse_gw_results = parse_total_time = False

        for line in output.split("\n"):
            if parse_total_time:
                if match := end_patt.search(line):
                    GW_results.update(end_normally=True)

                if match := total_time_patt.search(line):
                    GW_results.update(total_time=match[1])

            if parse_gw_results:
                if line.find("Dumping eigen energies") != -1:
                    parse_total_time = True
                    parse_gw_results = False
                    continue

                if match := GW_BANDS_results_patt.search(line):
                    dct = {}
                    dct.update(
                        band=match[1].strip(),
                        eKS=match[2],
                        eXX=match[3],
                        eQP_old=match[4],
                        z=match[5],
                        sigma_c_Linear=match[6],
                        eQP_Linear=match[7],
                        sigma_c_SCF=match[8],
                        eQP_SCF=match[9],
                    )
                    GW_results[match[1].strip()] = dct

                if n := GW_GAPS_results_patt.search(line):
                    dct = {}
                    dct.update(
                        Egap_KS=n[1],
                        Egap_QP_Linear=n[2],
                        Egap_QP_SCF=n[3],
                    )
                    GW_results["Gaps"] = dct

            if line.find("GW Results") != -1:
                parse_gw_results = True

        return GW_results


class BSEOutput:
    """
    A bse output file parser. The start...

    All energies are in eV.
    """

    def __init__(self, filename):
        """
        Args:
            filename: Filename to read.
        """
        self.filename = filename

        with zopen(filename, mode="rt", encoding="utf-8") as file:
            log_bse = file.read()

        # self.job_info = self._parse_preamble(preamble)
        self.exiton = self._parse_job(log_bse)

    @staticmethod
    def _parse_job(output):
        BSE_exitons_patt = re.compile(
            r"^exiton \s+ (\d+)  : \s+  ([\d.]+) \( \s+ ([-\d.]+) \) \s+ \| .*  ",
            re.VERBOSE,
        )

        end_patt = re.compile(r"\s*program returned normally\s*")

        total_time_patt = re.compile(r"\s*total \s+ time: \s+  ([\d.]+) .*", re.VERBOSE)

        BSE_results = {}
        parse_BSE_results = parse_total_time = False

        for line in output.split("\n"):
            if parse_total_time:
                if match := end_patt.search(line):
                    BSE_results.update(end_normally=True)

                if match := total_time_patt.search(line):
                    BSE_results.update(total_time=match[1])

            if parse_BSE_results:
                if line.find("FULL BSE main valence -> conduction transitions weight:") != -1:
                    parse_total_time = True
                    parse_BSE_results = False
                    continue

                if match := BSE_exitons_patt.search(line):
                    dct = {}
                    dct.update(bse_eig=match[2], osc_strength=match[3])
                    BSE_results[str(match[1].strip())] = dct

            if line.find("FULL BSE eig.(eV), osc. strength and dipoles:") != -1:
                parse_BSE_results = True

        return BSE_results
