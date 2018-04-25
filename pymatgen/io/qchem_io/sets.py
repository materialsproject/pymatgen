# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.inputs import QCInput

# Classes for reading/manipulating/writing QChem ouput files.

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QChemDictSet(QCInput):
    """
    Build a QCInput given all the various input parameters. Can be extended by standard implementations below.
    """

    def __init__(self, molecule, job_type, basis_set, scf_algorithm, dft_rung=4, pcm_dielectric=None,
                 max_scf_cycles=200, geom_opt_max_cycles=200):
        """
        Args:
            molecule (Pymatgen molecule object)
            job_type (str)
            basis_set (str)
            scf_algorithm (str)
            dft_rung (int)
            pcm_dielectric (str)
            max_scf_cycles (int)
            geom_opt_max_cycles (int)
        """
        self.molecule = molecule
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.dft_rung = dft_rung
        self.pcm_dielectric = pcm_dielectric
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles

        pcm_defaults = {"heavypoints": "194", "hpoints": "194",
                        "radii": "uff", "theory": "cpcm", "vdwscale": "1.1"}

        mypcm = {}
        mysolvent = {}
        myrem = {}
        myrem["job_type"] = job_type
        myrem["basis"] = self.basis_set
        myrem["max_scf_cycles"] = self.max_scf_cycles

        if self.dft_rung == 1:
            myrem["exchange"] = "B3LYP"
        elif self.dft_rung == 2:
            myrem["method"] = "B97-D3"
            myrem["dft_D"] = "D3_BJ"
        elif self.dft_rung == 3:
            myrem["method"] = "B97M-rV"
        elif self.dft_rung == 4:
            myrem["method"] = "wB97X-V"
        elif self.dft_rung == 5:
            myrem["method"] = "wB97M-V"
        else:
            raise ValueError("dft_rung should be between 1 and 5!")

        if self.job_type.lower() == "opt":
            myrem["geom_opt_max_cycles"] = self.geom_opt_max_cycles

        if self.pcm_dielectric != None:
            mypcm = pcm_defaults
            mysolvent["dielectric"] = self.pcm_dielectric
            myrem["solvent_method"] = 'pcm'

        super(QChemDictSet, self).__init__(self.molecule,
                                           rem=myrem, pcm=mypcm, solvent=mysolvent)


class OptSet(QChemDictSet):
    """
    QChemDictSet for a geometry optimization
    """

    def __init__(self, molecule, dft_rung=4, basis_set="6-311++G*", pcm_dielectric=None,
                 scf_algorithm="diis", max_scf_cycles=200, geom_opt_max_cycles=200):
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        super(OptSet, self).__init__(molecule=molecule, job_type="opt", dft_rung=dft_rung, pcm_dielectric=pcm_dielectric,
                                     basis_set=self.basis_set, scf_algorithm=self.scf_algorithm,
                                     max_scf_cycles=self.max_scf_cycles, geom_opt_max_cycles=self.geom_opt_max_cycles)


class SinglePointSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    def __init__(self, molecule, dft_rung=4, basis_set="6-311++G*",  pcm_dielectric=None, scf_algorithm="diis",
                 max_scf_cycles=200):
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super(SinglePointSet, self).__init__(molecule=molecule, job_type="sp", dft_rung=dft_rung,
                                             pcm_dielectric=pcm_dielectric, basis_set=self.basis_set,
                                             scf_algorithm=self.scf_algorithm, max_scf_cycles=self.max_scf_cycles)


class FreqSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    def __init__(self, molecule, dft_rung=4, basis_set="6-311++G*", pcm_dielectric=None, scf_algorithm="diis",
                 max_scf_cycles=200):
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super(FreqSet, self).__init__(molecule=molecule, job_type="freq", dft_rung=dft_rung,
                                      pcm_dielectric=pcm_dielectric, basis_set=self.basis_set,
                                      scf_algorithm=self.scf_algorithm, max_scf_cycles=self.max_scf_cycles)
