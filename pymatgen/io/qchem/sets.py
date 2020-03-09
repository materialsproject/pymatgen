# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for Qchem
"""

import logging
import os
from monty.io import zopen
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.utils import lower_and_check_unique


__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QChemDictSet(QCInput):
    """
    Build a QCInput given all the various input parameters. Can be extended by standard implementations below.
    """

    def __init__(self,
                 molecule,
                 job_type,
                 basis_set,
                 scf_algorithm,
                 dft_rung=4,
                 pcm_dielectric=None,
                 smd_solvent=None,
                 custom_smd=None,
                 max_scf_cycles=200,
                 geom_opt_max_cycles=200,
                 overwrite_inputs=None):
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
            overwrite_inputs (dict): This is dictionary of QChem input sections to add or overwrite variables,
            the available sections are currently rem, pcm, and solvent. So the accepted keys are rem, pcm, or solvent
            and the value is a dictionary of key value pairs relevant to the section. An example would be adding a
            new variable to the rem section that sets symmetry to false.
            ex. overwrite_inputs = {"rem": {"symmetry": "false"}}
            ***It should be noted that if something like basis is added to the rem dict it will overwrite
            the default basis.***
        """
        self.molecule = molecule
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.dft_rung = dft_rung
        self.pcm_dielectric = pcm_dielectric
        self.smd_solvent = smd_solvent
        self.custom_smd = custom_smd
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        self.overwrite_inputs = overwrite_inputs

        pcm_defaults = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1"
        }

        mypcm = {}
        mysolvent = {}
        mysmx = {}
        myrem = {}
        myrem["job_type"] = job_type
        myrem["basis"] = self.basis_set
        myrem["max_scf_cycles"] = self.max_scf_cycles
        myrem["gen_scfman"] = "true"
        myrem["xc_grid"] = "3"
        myrem["scf_algorithm"] = self.scf_algorithm
        myrem["resp_charges"] = "true"
        myrem["symmetry"] = "false"
        myrem["sym_ignore"] = "true"

        if self.dft_rung == 1:
            myrem["method"] = "b3lyp"
        elif self.dft_rung == 2:
            myrem["method"] = "b3lyp"
            myrem["dft_D"] = "D3_BJ"
        elif self.dft_rung == 3:
            myrem["method"] = "wb97xd"
        elif self.dft_rung == 4:
            myrem["method"] = "wb97xv"
        elif self.dft_rung == 5:
            myrem["method"] = "wb97mv"
        else:
            raise ValueError("dft_rung should be between 1 and 5!")

        if self.job_type.lower() == "opt":
            myrem["geom_opt_max_cycles"] = self.geom_opt_max_cycles

        if self.pcm_dielectric is not None and self.smd_solvent is not None:
            raise ValueError("Only one of pcm or smd may be used for solvation.")

        if self.pcm_dielectric is not None:
            mypcm = pcm_defaults
            mysolvent["dielectric"] = self.pcm_dielectric
            myrem["solvent_method"] = 'pcm'

        if self.smd_solvent is not None:
            if self.smd_solvent == "custom":
                mysmx["solvent"] = "other"
            else:
                mysmx["solvent"] = self.smd_solvent
            myrem["solvent_method"] = "smd"
            myrem["ideriv"] = "1"
            if self.smd_solvent == "custom" or self.smd_solvent == "other":
                if self.custom_smd is None:
                    raise ValueError(
                        'A user-defined SMD requires passing custom_smd, a string' +
                        ' of seven comma separated values in the following order:' +
                        ' dielectric, refractive index, acidity, basicity, surface' +
                        ' tension, aromaticity, electronegative halogenicity')

        if self.overwrite_inputs:
            for sec, sec_dict in self.overwrite_inputs.items():
                if sec == "rem":
                    temp_rem = lower_and_check_unique(sec_dict)
                    for k, v in temp_rem.items():
                        myrem[k] = v
                if sec == "pcm":
                    temp_pcm = lower_and_check_unique(sec_dict)
                    for k, v in temp_pcm.items():
                        mypcm[k] = v
                if sec == "solvent":
                    temp_solvent = lower_and_check_unique(sec_dict)
                    for k, v in temp_solvent.items():
                        mysolvent[k] = v
                if sec == "smx":
                    temp_smx = lower_and_check_unique(sec_dict)
                    for k, v in temp_smx.items():
                        mysmx[k] = v

        super().__init__(
            self.molecule, rem=myrem, pcm=mypcm, solvent=mysolvent, smx=mysmx)

    def write(self, input_file):
        """
        Args:
            input_file (str): Filename
        """
        self.write_file(input_file)
        if self.smd_solvent == "custom" or self.smd_solvent == "other":
            with zopen(os.path.join(os.path.dirname(input_file), "solvent_data"), 'wt') as f:
                f.write(self.custom_smd)


class OptSet(QChemDictSet):
    """
    QChemDictSet for a geometry optimization
    """

    def __init__(self,
                 molecule,
                 dft_rung=3,
                 basis_set="def2-tzvppd",
                 pcm_dielectric=None,
                 smd_solvent=None,
                 custom_smd=None,
                 scf_algorithm="diis",
                 max_scf_cycles=200,
                 geom_opt_max_cycles=200,
                 overwrite_inputs=None):
        """
        Args:
            molecule ():
            dft_rung ():
            basis_set ():
            pcm_dielectric ():
            smd_solvent ():
            custom_smd ():
            scf_algorithm ():
            max_scf_cycles ():
            geom_opt_max_cycles ():
            overwrite_inputs ():
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        super().__init__(
            molecule=molecule,
            job_type="opt",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            geom_opt_max_cycles=self.geom_opt_max_cycles,
            overwrite_inputs=overwrite_inputs)


class SinglePointSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    def __init__(self,
                 molecule,
                 dft_rung=3,
                 basis_set="def2-tzvppd",
                 pcm_dielectric=None,
                 smd_solvent=None,
                 custom_smd=None,
                 scf_algorithm="diis",
                 max_scf_cycles=200,
                 overwrite_inputs=None):
        """

        Args:
            molecule ():
            dft_rung ():
            basis_set ():
            pcm_dielectric ():
            smd_solvent ():
            custom_smd ():
            scf_algorithm ():
            max_scf_cycles ():
            overwrite_inputs ():
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="sp",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            overwrite_inputs=overwrite_inputs)


class FreqSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    def __init__(self,
                 molecule,
                 dft_rung=3,
                 basis_set="def2-tzvppd",
                 pcm_dielectric=None,
                 smd_solvent=None,
                 custom_smd=None,
                 scf_algorithm="diis",
                 max_scf_cycles=200,
                 overwrite_inputs=None):
        """
        Args:
            molecule ():
            dft_rung ():
            basis_set ():
            pcm_dielectric ():
            smd_solvent ():
            custom_smd ():
            scf_algorithm ():
            max_scf_cycles ():
            overwrite_inputs ():
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="freq",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            overwrite_inputs=overwrite_inputs)
