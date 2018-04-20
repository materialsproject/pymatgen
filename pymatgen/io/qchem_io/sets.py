# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging
from pymatgen.core import Molecule
from pymatgen.io.qchem_io.inputs import QCInput

# Classes for reading/manipulating/writing QChem ouput files.

__author__ = "Samuel Blau, Brandon Woods, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)

class QChemDictSet(QCInput):
    """
    Build a QCInput given all the various input parameters. Can be extended by standard implementations below.
    """
    def __init__(self, molecule, job_type, basis_set, SCF_algorithm, DFT_rung=4, PCM_solvent=None, max_scf_cycles=200, geom_opt_max_cycles=200):

        """
        Args:
            molecule (Pymatgen molecule)
            job_type (str)
            basis_set (str)
            SCF_algorithm (str)
            DFT_rung (int)
            PCM_solvent (str)
            max_scf_cycles (int)
            geom_opt_max_cycles (int)
        """
        if isinstance(molecule, Molecule):
            self.molecule = molecule
        else:
            raise ValueError('molecule must be a Pymatgen Molecule object!')
        self.job_type = job_type
        self.basis_set = basis_set
        self.SCF_algorithm = SCF_algorithm
        self.DFT_rung = DFT_rung
        self.PCM_solvent = PCM_solvent
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles

        myrem = {}
        myrem["job_type"] = job_type
        myrem["basis"] = self.basis_set
        myrem["max_scf_cycles"] = self.max_scf_cycles

        if self.DFT_rung == 1:
            myrem["exchange"] = "B3LYP"
        elif self.DFT_rung == 2:
            myrem["method"] = "B97-D3"
            myrem["DFT_D"] = "D3_BJ"
        elif self.DFT_rung == 3:
            myrem["method"] = "B97M-rV"
        elif self.DFT_rung == 4:
            myrem["method"] = "wB97X-V"
        elif self.DFT_rung == 5:
            myrem["method"] = "wB97M-V"
        else:
            print("DFT_rung should be between 1 and 5!")

        if self.job_type.lower() == "opt":
            myrem["geom_opt_max_cycles"] = self.geom_opt_max_cycles

        if self.PCM_solvent != None:
            print("Need PCM input implementation!")

        super(QChemDictSet, self).__init__(self.molecule, myrem)

class OptSet(QChemDictSet):
    """
    QChemDictSet for a geometry optimization
    """

    defaults = {"basis": "6-311++G*", "SCF_algorithm": "diis", "max_scf_cycles": 200, "geom_opt_max_cycles": 200}

    def __init__(self, molecule, DFT_rung=4, PCM_solvent=None):
        self.basis_set = defaults.get("basis")
        self.SCF_algorithm = defaults.get("SCF_algorithm")
        self.max_scf_cycles = defaults.get("max_scf_cycles")
        self.geom_opt_max_cycles = defaults.get("geom_opt_max_cycles")
        super(OptSet, self).__init__(molecule=molecule, job_type="opt", DFT_rung=DFT_rung, PCM_solvent=PCM_solvent, basis_set=self.basis_set, SCF_algorithm=self.SCF_algorithm, max_scf_cycles=self.max_scf_cycles, geom_opt_max_cycles=self.geom_opt_max_cycles)

class SinglePointSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    defaults = {"basis": "6-311++G*", "SCF_algorithm": "diis", "max_scf_cycles": 200}

    def __init__(self, molecule, DFT_rung=4, PCM_solvent=None):
        self.basis_set = defaults.get("basis")
        self.SCF_algorithm = defaults.get("SCF_algorithm")
        self.max_scf_cycles = defaults.get("max_scf_cycles")
        super(SinglePointSet, self).__init__(molecule=molecule, job_type="sp", DFT_rung=DFT_rung, PCM_solvent=PCM_solvent, basis_set=self.basis_set, SCF_algorithm=self.SCF_algorithm, max_scf_cycles=self.max_scf_cycles)

class FreqSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    defaults = {"basis": "6-311++G*", "SCF_algorithm": "diis", "max_scf_cycles": 200}

    def __init__(self, molecule, DFT_rung=4, PCM_solvent=None):
        self.basis_set = defaults.get("basis")
        self.SCF_algorithm = defaults.get("SCF_algorithm")
        self.max_scf_cycles = defaults.get("max_scf_cycles")
        super(FreqSet, self).__init__(molecule=molecule, job_type="freq", DFT_rung=DFT_rung, PCM_solvent=PCM_solvent, basis_set=self.basis_set, SCF_algorithm=self.SCF_algorithm, max_scf_cycles=self.max_scf_cycles)

