# coding: utf-8


import re
import os
import shutil
import subprocess

from string import Template
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core.structure import Molecule

"""
This module implements input and output processing for Fiesta (http://perso.neel.cnrs.fr/xavier.blase/fiesta/index.html).

and

-Nwchem2Fiesta class: to create the input files needed for a Fiesta run
-Fiesta_run: run gw_fiesta and bse_fiesta
-Localised Basis set reader
"""

__author__ = 'ndardenne'
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__email__ = "n.dardenne@uclouvain.be"
__date__ = "24/5/15"


class Nwchem2Fiesta(MSONable):
    """
    To run NWCHEM2FIESTA inside python:
    
    If nwchem.nw is the input, nwchem.out the output, and structure.movecs the
    "movecs" file, the syntax to run NWCHEM2FIESTA is: NWCHEM2FIESTA
    nwchem.nw  nwchem.nwout  structure.movecs > log_n2f
    """

    def __init__(self, folder, filename="nwchem", log_file="log_n2f"):
        """
        folder: where are stored the nwchem
        filename: name of nwchem files read by NWCHEM2FIESTA (filename.nw, filename.nwout and filename.movecs)
        logfile: logfile of NWCHEM2FIESTA

        the run method launchs NWCHEM2FIESTA

        """

        self.folder = folder
        self.filename = filename
        self.log_file = log_file

        self._NWCHEM2FIESTA_cmd = "NWCHEM2FIESTA"
        self._nwcheminput_fn = filename + ".nw"
        self._nwchemoutput_fn = filename + ".nwout"
        self._nwchemmovecs_fn = filename + ".movecs"

    def run(self):
        """
        Performs actual NWCHEM2FIESTA run
        """

        init_folder = os.getcwd()
        os.chdir(self.folder)

        with zopen(self.log_file, 'w') as fout:
            subprocess.call([self._NWCHEM2FIESTA_cmd, self._nwcheminput_fn,
                             self._nwchemoutput_fn, self._nwchemmovecs_fn],
                            stdout=fout)

        os.chdir(init_folder)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "filename": self.filename,
                "folder": self.folder}

    @classmethod
    def from_dict(cls, d):
        return Nwchem2Fiesta(folder=d["folder"], filename=d["filename"])


class Fiesta_run(MSONable):
    """
    To run FIESTA inside python:
     if grid is [x,x] then bse runs
     if grid is [x,x,y] the fiesta(gw) runs
     otherwise it breaks
    """

    def __init__(self, folder=os.getcwd(), grid=[2, 2, 2], log_file="log"):
        """
        folder:
        logfile: logfile of Fiesta
        """

        self.folder = folder
        self.log_file = log_file
        self.grid = grid

    def run(self):
        if len(self.grid) == 3:
            self.mpi_procs = self.grid[0] * self.grid[1] * self.grid[2]
            self.gw_run()
        elif len(self.grid) == 2:
            self.mpi_procs = self.grid[0] * self.grid[1]
            self.bse_run()
        else:
            raise ValueError(
                "Wrong grid size: must be [nrow, ncolumn, nslice] for gw of [nrow, nslice] for bse")

    def gw_run(self):
        """
        Performs FIESTA (gw) run
        """

        if self.folder != os.getcwd():
            init_folder = os.getcwd()
            os.chdir(self.folder)

        with zopen(self.log_file, 'w') as fout:
            subprocess.call(["mpirun", "-n", str(self.mpi_procs), "fiesta",
                             str(self.grid[0]), str(self.grid[1]),
                             str(self.grid[2])], stdout=fout)

        if self.folder != os.getcwd():
            os.chdir(init_folder)

    def bse_run(self):
        """
        Performs BSE run
        """

        if self.folder != os.getcwd():
            init_folder = os.getcwd()
            os.chdir(self.folder)

        with zopen(self.log_file, 'w') as fout:
            subprocess.call(
                ["mpirun", "-n", str(self.mpi_procs), "bse", str(self.grid[0]),
                 str(self.grid[1])], stdout=fout)

        if self.folder != os.getcwd():
            os.chdir(init_folder)

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "log_file": self.log_file,
                "grid": self.grid,
                "folder": self.folder}

    @classmethod
    def from_dict(cls, d):
        return Fiesta_run(folder=d["folder"], grid=d["grid"],
                          log_file=d['log_file'])


class Basis_set_reader:
    """
    A basis set reader.
    Args:
        filename: Filename to read.

        Basis set are stored in data as a dict:
         :key l_zeta_ng for each nl orbitals which contain list of tuple (alpha, coef) for each of the ng gaussians in l_zeta orbital
    """

    def __init__(self, filename):
        self.filename = filename

        with zopen(filename) as f:
            basis_set = f.read()

        self.data = self._parse_file(basis_set)
        # compute the number of nlm orbitals per atom
        self.data.update(n_nlmo=self.set_n_nlmo())

    def _parse_file(self, input):

        lmax_nnlo_patt = re.compile(r"\s* (\d+) \s+ (\d+) \s+ \# .* ",
                                    re.VERBOSE)

        nl_orbital_patt = re.compile(r"\s* (\d+) \s+ (\d+) \s+ (\d+) \s+ \# .* ",
                                     re.VERBOSE)

        coef_alpha_patt = re.compile(r"\s* ([-\d.\D]+) \s+ ([-\d.\D]+) \s* ",
                                     re.VERBOSE)

        preamble = []
        basis_set = {}
        parse_preamble = False
        parse_lmax_nnlo = False
        parse_nl_orbital = False

        for line in input.split("\n"):

            if parse_nl_orbital:
                m = nl_orbital_patt.search(line)
                n = coef_alpha_patt.search(line)
                if m:
                    l = m.group(1)
                    zeta = m.group(2)
                    ng = m.group(3)
                    basis_set[l + "_" + zeta + "_" + ng] = []
                elif n:
                    alpha = n.group(1)
                    coef = n.group(2)
                    basis_set[l + "_" + zeta + "_" + ng].append((alpha, coef))
            elif parse_lmax_nnlo:
                m = lmax_nnlo_patt.search(line)
                if m:
                    lmax = m.group(1)
                    nnlo = m.group(2)
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
        """
                :return: the number of nlm orbitals for the basis set
        """

        nnlmo = 0

        data_tmp = self.data
        data_tmp.pop('lmax')
        data_tmp.pop('n_nlo')
        data_tmp.pop('preamble')

        for l_zeta_ng in data_tmp:
            l = l_zeta_ng.split("_")[0]
            nnlmo = nnlmo + (2 * int(l) + 1)

        return str(nnlmo)

    def infos_on_basis_set(self):
        """
        infos on the basis set as in Fiesta log
        """
        o = []
        o.append("=========================================")
        o.append("Reading basis set:")
        o.append("")
        o.append(" Basis set for {} atom ".format(str(self.filename)))
        o.append(" Maximum angular momentum = {}".format(self.data['lmax']))
        o.append(" Number of atomics orbitals = {}".format(self.data['n_nlo']))
        o.append(" Number of nlm orbitals = {}".format(self.data['n_nlmo']))
        o.append("=========================================")

        return str(0)


class FiestaInput(MSONable):
    """
    Input File for Fiesta called "cell.in" by default (mandatory in Fiesta for now)
    """

    def __init__(self, mol,
                 correlation_grid={'dE_grid': u'0.500', 'n_grid': u'14'},
                 Exc_DFT_option={'rdVxcpsi': u'1'},
                 COHSEX_options={'eigMethod': u'C', 'mix_cohsex': u'0.500',
                                 'nc_cohsex': u'0', 'nit_cohsex': u'0',
                                 'nv_cohsex': u'0', 'resMethod': u'V',
                                 'scf_cohsex_wf': u'0'},
                 GW_options={'nc_corr': u'10', 'nit_gw': u'3',
                             'nv_corr': u'10'},
                 BSE_TDDFT_options={'do_bse': u'1', 'do_tddft': u'0',
                                    'nc_bse': u'382', 'nit_bse': u'50',
                                    'npsi_bse': u'1', 'nv_bse': u'21'}):
        """
        :param mol: pymatgen mol
        :param correlation_grid: dict
        :param Exc_DFT_option: dict
        :param COHSEX_options: dict
        :param GW_options: dict
        :param BSE_TDDFT_options: dict
        """

        self._mol = mol
        self.correlation_grid = correlation_grid
        self.Exc_DFT_option = Exc_DFT_option
        self.COHSEX_options = COHSEX_options
        self.GW_options = GW_options
        self.BSE_TDDFT_options = BSE_TDDFT_options

    def set_auxiliary_basis_set(self, folder, auxiliary_folder,
                                auxiliary_basis_set_type="aug_cc_pvtz"):
        """
        copy in the desired folder the needed auxiliary basis set "X2.ion" where X is a specie.
        :param auxiliary_folder: folder where the auxiliary basis sets are stored
        :param auxiliary_basis_set_type: type of basis set (string to be found in the extension of the file name; must be in lower case)
               ex: C2.ion_aug_cc_pvtz_RI_Weigend find "aug_cc_pvtz"
        """

        list_files = os.listdir(auxiliary_folder)

        for specie in self._mol.symbol_set:
            for file in list_files:
                if file.upper().find(
                                specie.upper() + "2") != -1 and file.lower().find(
                        auxiliary_basis_set_type) != -1:
                    shutil.copyfile(auxiliary_folder + "/" + file,
                                    folder + "/" + specie + "2.ion")

    def set_GW_options(self, nv_band=10, nc_band=10, n_iteration=5, n_grid=6,
                       dE_grid=0.5):
        """
        Set parameters in cell.in for a GW computation
        :param nv__band: number of valence bands to correct with GW
        :param nc_band: number of conduction bands to correct with GW
        :param n_iteration: number of iteration
        :param n_grid and dE_grid:: number of points and spacing in eV for correlation grid
        """

        self.GW_options.update(nv_corr=nv_band, nc_corr=nc_band,
                               nit_gw=n_iteration)
        self.correlation_grid.update(dE_grid=dE_grid, n_grid=n_grid)

    def make_FULL_BSE_Densities_folder(self, folder):
        """
        mkdir "FULL_BSE_Densities" folder (needed for bse run) in the desired folder
        """

        if os.path.exists(folder + "/FULL_BSE_Densities"):
            return "FULL_BSE_Densities folder already exists"
        else:
            os.makedirs(folder + "/FULL_BSE_Densities")
            return "makedirs FULL_BSE_Densities folder"

    def set_BSE_options(self, n_excitations=10, nit_bse=200):
        """
        Set parameters in cell.in for a BSE computation
        :param nv_bse: number of valence bands
        :param nc_bse: number of conduction bands
        :param n_excitations: number of excitations
        :param nit_bse: number of iterations
        """

        self.BSE_TDDFT_options.update(npsi_bse=n_excitations, nit_bse=nit_bse)

    def dump_BSE_data_in_GW_run(self, BSE_dump=True):
        """
        :param BSE_dump: boolean
        :return: set the "do_bse" variable to one in cell.in
        """

        if BSE_dump:
            self.BSE_TDDFT_options.update(do_bse=1, do_tddft=0)
        else:
            self.BSE_TDDFT_options.update(do_bse=0, do_tddft=0)

    def dump_TDDFT_data_in_GW_run(self, TDDFT_dump=True):
        """
        :param TDDFT_dump: boolen
        :return: set the do_tddft variable to one in cell.in
        """
        if TDDFT_dump == True:
            self.BSE_TDDFT_options.update(do_bse=0, do_tddft=1)
        else:
            self.BSE_TDDFT_options.update(do_bse=0, do_tddft=0)

    @property
    def infos_on_system(self):
        """
        Returns infos on initial parameters as in the log file of Fiesta
        """

        o = []
        o.append("=========================================")
        o.append("Reading infos on system:")
        o.append("")
        o.append(" Number of atoms = {} ; number of species = {}".format(
            int(self._mol.composition.num_atoms), len(self._mol.symbol_set)))
        o.append(" Number of valence bands = {}".format(
            int(self._mol.nelectrons / 2)))
        o.append(" Sigma grid specs: n_grid = {} ;  dE_grid = {} (eV)".format(
            self.correlation_grid['n_grid'], self.correlation_grid['dE_grid']))
        if int(self.Exc_DFT_option['rdVxcpsi']) == 1:
            o.append(" Exchange and correlation energy read from Vxcpsi.mat")
        elif int(self.Exc_DFT_option['rdVxcpsi']) == 0:
            o.append(" Exchange and correlation energy re-computed")

        if self.COHSEX_options['eigMethod'] == "C":
            o.append(
                " Correcting  {} valence bands and   {} conduction bands at COHSEX level".format(
                    self.COHSEX_options['nv_cohsex'],
                    self.COHSEX_options['nc_cohsex']))
            o.append(" Performing   {} diagonal COHSEX iterations".format(
                self.COHSEX_options['nit_cohsex']))
        elif self.COHSEX_options['eigMethod'] == "HF":
            o.append(
                " Correcting  {} valence bands and   {} conduction bands at HF level".format(
                    self.COHSEX_options['nv_cohsex'],
                    self.COHSEX_options['nc_cohsex']))
            o.append(" Performing   {} diagonal HF iterations".format(
                self.COHSEX_options['nit_cohsex']))

        o.append(" Using resolution of identity : {}".format(
            self.COHSEX_options['resMethod']))
        o.append(
            " Correcting  {} valence bands and  {} conduction bands at GW level".format(
                self.GW_options['nv_corr'], self.GW_options['nc_corr']))
        o.append(
            " Performing   {} GW iterations".format(self.GW_options['nit_gw']))

        if int(self.BSE_TDDFT_options['do_bse']) == 1:
            o.append(" Dumping data for BSE treatment")

        if int(self.BSE_TDDFT_options['do_tddft']) == 1:
            o.append(" Dumping data for TD-DFT treatment")
        o.append("")
        o.append(" Atoms in cell cartesian A:")
        symbols = []
        for syb in self._mol.symbol_set:
            symbols.append(syb)

        for site in self._mol:
            o.append(" {} {} {} {}".format(site.x, site.y,
                                           site.z, int(
                    symbols.index(site.specie.symbol)) + 1))

        o.append("=========================================")

        return str(o)

    @property
    def molecule(self):
        """
        Returns molecule associated with this FiestaInput.
        """
        return self._mol

    def __str__(self):

        symbols = []
        for syb in self._mol.symbol_set:
            symbols.append(syb)

        geometry = []
        for site in self._mol:
            geometry.append(" {} {} {} {}".format(site.x, site.y,
                                                  site.z, int(
                    symbols.index(site.specie.symbol)) + 1))

        t = Template("""# number of atoms and species
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
# number of occp. and virtual bands fo BSE: nocore and up to 40 eVs
    $nv_bse   $nc_bse
# number of excitations needed and number of iterations
    $npsi_bse   $nit_bse
# list of symbols in order
$symbols
# scaling factor
    1.000
# atoms x,y,z cartesian .. will be multiplied by scale
$geometry
            """)

        return t.substitute(nat=int(self._mol.composition.num_atoms),
                            nsp=len(self._mol.symbol_set),
                            nvbands=int(self._mol.nelectrons / 2),
                            n_grid=self.correlation_grid['n_grid'],
                            dE_grid=self.correlation_grid['dE_grid'],
                            rdVxcpsi=self.Exc_DFT_option['rdVxcpsi'],
                            nv_cohsex=self.COHSEX_options['nv_cohsex'],
                            nc_cohsex=self.COHSEX_options['nc_cohsex'],
                            eigMethod=self.COHSEX_options['eigMethod'],
                            nit_cohsex=self.COHSEX_options['nit_cohsex'],
                            resMethod=self.COHSEX_options['resMethod'],
                            scf_cohsex_wf=self.COHSEX_options['scf_cohsex_wf'],
                            mix_cohsex=self.COHSEX_options['mix_cohsex'],
                            nv_corr=self.GW_options['nv_corr'],
                            nc_corr=self.GW_options['nc_corr'],
                            nit_gw=self.GW_options['nit_gw'],
                            do_bse=self.BSE_TDDFT_options['do_bse'],
                            do_tddft=self.BSE_TDDFT_options['do_tddft'],
                            nv_bse=self.BSE_TDDFT_options['nv_bse'],
                            nc_bse=self.BSE_TDDFT_options['nc_bse'],
                            npsi_bse=self.BSE_TDDFT_options['npsi_bse'],
                            nit_bse=self.BSE_TDDFT_options['nit_bse'],
                            symbols="\n".join(symbols),
                            geometry="\n".join(geometry))

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())

    def as_dict(self):
        return {
            "mol": self._mol.as_dict(),
            "correlation_grid": self.correlation_grid,
            "Exc_DFT_option": self.Exc_DFT_option,
            "COHSEX_options": self.COHSEX_options,
            "GW_options": self.GW_options,
            "BSE_TDDFT_options": self.BSE_TDDFT_options
        }

    @classmethod
    def from_dict(cls, d):
        return FiestaInput(Molecule.from_dict(d["mol"]),
                           correlation_grid=d["correlation_grid"],
                           Exc_DFT_option=d["Exc_DFT_option"],
                           COHSEX_options=d["geometry_options"],
                           GW_options=d["symmetry_options"],
                           BSE_TDDFT_options=d["memory_options"])

    @classmethod
    def from_string(cls, string_input):
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
        l = lines.pop(0).strip()
        toks = l.split()
        nat = toks[0]
        nsp = toks[1]
        # number of valence bands
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        nvbands = toks[0]

        # correlation_grid
        # number of points and spacing in eV for correlation grid
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        correlation_grid['n_grid'] = toks[0]
        correlation_grid['dE_grid'] = toks[1]

        # Exc DFT
        # relire=1 ou recalculer=0 Exc DFT
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        Exc_DFT_option['rdVxcpsi'] = toks[0]

        # COHSEX
        # number of COHSEX corrected occp and unoccp bands: C=COHSEX  H=HF
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        COHSEX_options['nv_cohsex'] = toks[0]
        COHSEX_options['nc_cohsex'] = toks[1]
        COHSEX_options['eigMethod'] = toks[2]
        # number of COHSEX iter, scf on wfns, mixing coeff; V=RI-V  I=RI-D
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        COHSEX_options['nit_cohsex'] = toks[0]
        COHSEX_options['resMethod'] = toks[1]
        COHSEX_options['scf_cohsex_wf'] = toks[2]
        COHSEX_options['mix_cohsex'] = toks[3]

        # GW
        # number of GW corrected occp and unoccp bands
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        GW_options['nv_corr'] = toks[0]
        GW_options['nc_corr'] = toks[1]
        # number of GW iterations
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        GW_options['nit_gw'] = toks[0]

        # BSE
        # dumping for BSE and TDDFT
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        BSE_TDDFT_options['do_bse'] = toks[0]
        BSE_TDDFT_options['do_tddft'] = toks[1]
        # number of occp. and virtual bands fo BSE: nocore and up to 40 eVs
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        BSE_TDDFT_options['nv_bse'] = toks[0]
        BSE_TDDFT_options['nc_bse'] = toks[1]
        # number of excitations needed and number of iterations
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        BSE_TDDFT_options['npsi_bse'] = toks[0]
        BSE_TDDFT_options['nit_bse'] = toks[1]

        # Molecule
        # list of symbols in order
        lines.pop(0)
        atname = []
        i = int(nsp)
        while i != 0:
            l = lines.pop(0).strip()
            toks = l.split()
            atname.append(toks[0])
            i -= 1

        # scaling factor
        lines.pop(0)
        l = lines.pop(0).strip()
        toks = l.split()
        scale = toks[0]
        # atoms x,y,z cartesian .. will be multiplied by scale
        lines.pop(0)
        # Parse geometry
        species = []
        coords = []
        i = int(nat)
        while i != 0:
            l = lines.pop(0).strip()
            toks = l.split()
            coords.append([float(j) for j in toks[0:3]])
            species.append(atname[int(toks[3]) - 1])
            i -= 1

        mol = Molecule(species, coords)

        return FiestaInput(mol=mol, correlation_grid=correlation_grid,
                           Exc_DFT_option=Exc_DFT_option,
                           COHSEX_options=COHSEX_options,
                           GW_options=GW_options,
                           BSE_TDDFT_options=BSE_TDDFT_options)

    @classmethod
    def from_file(cls, filename):
        """
        Read an Fiesta input from a file. Currently tested to work with
        files generated from this class itself.

        Args:
            filename: Filename to parse.

        Returns:
            FiestaInput object
        """
        with zopen(filename) as f:
            return cls.from_string(f.read())


class FiestaOutput:
    """
    A Fiesta output file parser.

    All energies are in eV.

    Args:
        filename: Filename to read.
    """

    def __init__(self, filename):
        self.filename = filename

        with zopen(filename) as f:
            data = f.read()

        chunks = re.split(r"GW Driver iteration", data)

        # preamble: everything before the first GW Driver iteration
        preamble = chunks.pop(0)

        # self.job_info = self._parse_preamble(preamble)
        self.data = [self._parse_job(c) for c in chunks]

    def _parse_job(self, output):

        GW_BANDS_results_patt = re.compile(
            r"^<it.*  \| \s+ (\D+\d*) \s+ \| \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ \| "
            r" \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ \|"
            r" \s+ ([-\d.]+) \s+ ([-\d.]+) \s+ ", re.VERBOSE)

        GW_GAPS_results_patt = re.compile(
            r"^<it.*  \| \s+ Egap_KS \s+ = \s+ ([-\d.]+) \s+ \| \s+ Egap_QP \s+ = \s+ ([-\d.]+) \s+ \| "
            r" \s+ Egap_QP \s+ = \s+ ([-\d.]+) \s+ \|", re.VERBOSE)

        end_patt = re.compile(r"\s*program returned normally\s*")

        total_time_patt = re.compile(r"\s*total \s+ time: \s+  ([\d.]+) .*",
                                     re.VERBOSE)

        error_defs = {
            "calculations not reaching convergence": "Bad convergence"}

        GW_results = {}
        parse_gw_results = False
        parse_total_time = False

        for l in output.split("\n"):

            if parse_total_time:
                m = end_patt.search(l)
                if m:
                    GW_results.update(end_normally=True)

                m = total_time_patt.search(l)
                if m:
                    GW_results.update(total_time=m.group(1))

            if parse_gw_results:
                if l.find("Dumping eigen energies") != -1:
                    parse_total_time = True
                    parse_gw_results = False
                    continue
                else:
                    m = GW_BANDS_results_patt.search(l)
                    if m:
                        d = {}
                        d.update(band=m.group(1).strip(), eKS=m.group(2),
                                 eXX=m.group(3), eQP_old=m.group(4),
                                 z=m.group(5), sigma_c_Linear=m.group(6),
                                 eQP_Linear=m.group(7),
                                 sigma_c_SCF=m.group(8), eQP_SCF=m.group(9))
                        GW_results[m.group(1).strip()] = d

                    n = GW_GAPS_results_patt.search(l)
                    if n:
                        d = {}
                        d.update(Egap_KS=n.group(1), Egap_QP_Linear=n.group(2),
                                 Egap_QP_SCF=n.group(3))
                        GW_results["Gaps"] = d

            if l.find("GW Results") != -1:
                parse_gw_results = True

        return GW_results


class BSEOutput:
    """
    A bse output file parser. The start...

    All energies are in eV.

    Args:
        filename: Filename to read.
    """

    def __init__(self, filename):
        self.filename = filename

        with zopen(filename) as f:
            log_bse = f.read()

        # self.job_info = self._parse_preamble(preamble)
        self.exiton = self._parse_job(log_bse)

    def _parse_job(self, output):

        BSE_exitons_patt = re.compile(
            r"^exiton \s+ (\d+)  : \s+  ([\d.]+) \( \s+ ([-\d.]+) \) \s+ \| .*  ",
            re.VERBOSE)

        end_patt = re.compile(r"\s*program returned normally\s*")

        total_time_patt = re.compile(r"\s*total \s+ time: \s+  ([\d.]+) .*",
                                     re.VERBOSE)

        error_defs = {
            "calculations not reaching convergence": "Bad convergence"}

        BSE_results = {}
        parse_BSE_results = False
        parse_total_time = False

        for l in output.split("\n"):

            if parse_total_time:
                m = end_patt.search(l)
                if m:
                    BSE_results.update(end_normally=True)

                m = total_time_patt.search(l)
                if m:
                    BSE_results.update(total_time=m.group(1))

            if parse_BSE_results:
                if l.find(
                        "FULL BSE main valence -> conduction transitions weight:") != -1:
                    parse_total_time = True
                    parse_BSE_results = False
                    continue
                else:
                    m = BSE_exitons_patt.search(l)
                    if m:
                        d = {}
                        d.update(bse_eig=m.group(2), osc_strength=m.group(3))
                        BSE_results[str(m.group(1).strip())] = d

            if l.find("FULL BSE eig.(eV), osc. strength and dipoles:") != -1:
                parse_BSE_results = True

        return BSE_results
