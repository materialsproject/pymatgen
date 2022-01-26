# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for Qchem
"""

import logging
import os
from typing import Dict, List, Literal, Optional

from monty.io import zopen

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.utils import lower_and_check_unique

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2021, The Materials Project"
__version__ = "0.1"

logger = logging.getLogger(__name__)


class QChemDictSet(QCInput):
    """
    Build a QCInput given all the various input parameters. Can be extended by standard implementations below.
    """

    def __init__(
        self,
        molecule: Molecule,
        job_type: str,
        basis_set: str,
        scf_algorithm: str,
        dft_rung: int = 4,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        opt_variables: Optional[Dict[str, List]] = None,
        scan_variables: Optional[Dict[str, List]] = None,
        max_scf_cycles: int = 200,
        geom_opt_max_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. For example, "def2-tzvpd".
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            opt_variables (dict): A dictionary of opt sections, where each opt section is a key
                and the corresponding values are a list of strings. Strings must be formatted
                as instructed by the QChem manual. The different opt sections are: CONSTRAINT, FIXED,
                DUMMY, and CONNECT.

                Ex. opt = {"CONSTRAINT": ["tors 2 3 4 5 25.0", "tors 2 5 7 9 80.0"], "FIXED": ["2 XY"]}
            scan_variables (dict): A dictionary of scan variables. Because two constraints of the
                same type are allowed (for instance, two torsions or two bond stretches), each TYPE of
                variable (stre, bend, tors) should be its own key in the dict, rather than each variable.
                Note that the total number of variable (sum of lengths of all lists) CANNOT be more than two.

                Ex. scan_variables = {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180 15"]}
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            nbo_params (list): A list of strings for the desired NBO params. If an empty list is passed,
                default NBO analysis will be performed. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
        """
        self.molecule = molecule
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.dft_rung = dft_rung
        self.pcm_dielectric = pcm_dielectric
        self.smd_solvent = smd_solvent
        self.custom_smd = custom_smd
        self.opt_variables = opt_variables
        self.scan_variables = scan_variables
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        self.plot_cubes = plot_cubes
        self.nbo_params = nbo_params
        self.overwrite_inputs = overwrite_inputs
        self.vdw_mode = vdw_mode

        pcm_defaults = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }

        plots_defaults = {"grid_spacing": "0.05", "total_density": "0"}

        if self.opt_variables is None:
            myopt = {}
        else:
            myopt = self.opt_variables

        if self.scan_variables is None:
            myscan = {}
        else:
            myscan = self.scan_variables

        mypcm = {}
        mysolvent = {}
        mysmx = {}
        myvdw = {}
        myplots = {}
        myrem = {}
        myrem["job_type"] = job_type
        myrem["basis"] = self.basis_set
        myrem["max_scf_cycles"] = str(self.max_scf_cycles)
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

        if self.job_type.lower() in ["opt", "ts", "pes_scan"]:
            myrem["geom_opt_max_cycles"] = str(self.geom_opt_max_cycles)

        if self.pcm_dielectric is not None and self.smd_solvent is not None:
            raise ValueError("Only one of pcm or smd may be used for solvation.")

        if self.pcm_dielectric is not None:
            mypcm = pcm_defaults
            mysolvent["dielectric"] = self.pcm_dielectric
            myrem["solvent_method"] = "pcm"

        if self.smd_solvent is not None:
            if self.smd_solvent == "custom":
                mysmx["solvent"] = "other"
            else:
                mysmx["solvent"] = self.smd_solvent
            myrem["solvent_method"] = "smd"
            myrem["ideriv"] = "1"
            if self.smd_solvent in ("custom", "other"):
                if self.custom_smd is None:
                    raise ValueError(
                        "A user-defined SMD requires passing custom_smd, a string"
                        + " of seven comma separated values in the following order:"
                        + " dielectric, refractive index, acidity, basicity, surface"
                        + " tension, aromaticity, electronegative halogenicity"
                    )

        if self.plot_cubes:
            myplots = plots_defaults
            myrem["plots"] = "true"
            myrem["make_cube_files"] = "true"

        if self.nbo_params is not None:
            myrem["nbo"] = "true"

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
                if sec == "scan":
                    temp_scan = lower_and_check_unique(sec_dict)
                    for k, v in temp_scan.items():
                        myscan[k] = v
                if sec == "van_der_waals":
                    temp_vdw = lower_and_check_unique(sec_dict)
                    for k, v in temp_vdw.items():
                        myvdw[k] = v
                    # set the PCM section to read custom radii
                    mypcm["radii"] = "read"
                if sec == "plots":
                    temp_plots = lower_and_check_unique(sec_dict)
                    for k, v in temp_plots.items():
                        myplots[k] = v
                if sec == "nbo":
                    temp_plots = lower_and_check_unique(sec_dict)
                    for k, v in temp_plots.items():
                        myplots[k] = v
                if sec == "opt":
                    temp_opts = lower_and_check_unique(sec_dict)
                    for k, v in temp_opts.items():
                        myopt[k] = v

        super().__init__(
            self.molecule,
            rem=myrem,
            opt=myopt,
            pcm=mypcm,
            solvent=mysolvent,
            smx=mysmx,
            scan=myscan,
            van_der_waals=myvdw,
            vdw_mode=self.vdw_mode,
            plots=myplots,
            nbo=self.nbo_params,
        )

    def write(self, input_file: str):
        """
        Args:
            input_file (str): Filename
        """
        self.write_file(input_file)
        if self.smd_solvent in ("custom", "other"):
            with zopen(os.path.join(os.path.dirname(input_file), "solvent_data"), "wt") as f:
                f.write(self.custom_smd)


class SinglePointSet(QChemDictSet):
    """
    QChemDictSet for a single point calculation
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
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
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class OptSet(QChemDictSet):
    """
    QChemDictSet for a geometry optimization
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        opt_variables: Optional[Dict[str, List]] = None,
        geom_opt_max_cycles: int = 200,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
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
            opt_variables=opt_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            geom_opt_max_cycles=self.geom_opt_max_cycles,
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class TransitionStateSet(QChemDictSet):
    """
    QChemDictSet for a transition-state search
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        opt_variables: Optional[Dict[str, List]] = None,
        geom_opt_max_cycles: int = 200,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode="atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        super().__init__(
            molecule=molecule,
            job_type="ts",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            opt_variables=opt_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            geom_opt_max_cycles=self.geom_opt_max_cycles,
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class ForceSet(QChemDictSet):
    """
    QChemDictSet for a force (gradient) calculation
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="force",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class FreqSet(QChemDictSet):
    """
    QChemDictSet for a frequency calculation
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
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
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class PESScanSet(QChemDictSet):
    """
    QChemDictSet for a potential energy surface scan (PES_SCAN) calculation,
    used primarily to identify possible transition states or to sample different
    geometries.
    Note: Because there are no defaults that can be used for a PES scan (the
    variables are completely dependent on the molecular structure), by default
    scan_variables = None. However, a PES Scan job should not be run with less
    than one variable (or more than two variables).
    """

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvppd",
        scf_algorithm: str = "diis",
        dft_rung: int = 3,
        pcm_dielectric: Optional[float] = None,
        smd_solvent: Optional[str] = None,
        custom_smd: Optional[str] = None,
        max_scf_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: Optional[Dict] = None,
        opt_variables: Optional[Dict[str, List]] = None,
        scan_variables: Optional[Dict[str, List]] = None,
        overwrite_inputs: Optional[Dict] = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            opt_variables (dict): A dictionary of opt sections, where each opt section is a key
                and the corresponding values are a list of strings. Strings must be formatted
                as instructed by the QChem manual. The different opt sections are: CONSTRAINT, FIXED,
                DUMMY, and CONNECT.

                Ex. opt = {"CONSTRAINT": ["tors 2 3 4 5 25.0", "tors 2 5 7 9 80.0"], "FIXED": ["2 XY"]}
            scan_variables (dict): A dictionary of scan variables. Because two constraints of the
                same type are allowed (for instance, two torsions or two bond stretches), each TYPE of
                variable (stre, bend, tors) should be its own key in the dict, rather than each variable.
                Note that the total number of variable (sum of lengths of all lists) CANNOT be more than two.

                Ex. scan_variables = {"stre": ["3 6 1.5 1.9 0.1"], "tors": ["1 2 3 4 -180 180 15"]}
            basis_set (str): Basis set to use. (Default: "def2-tzvppd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            dft_rung (int): Select the DFT functional among 5 recommended levels of theory,
                in order of increasing accuracy/cost. 1 = B3LYP, 2=B3lYP+D3, 3=ωB97X-D,
                4=ωB97X-V, 5=ωB97M-V. (Default: 3)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

                **Note that the "rungs" in this argument do NOT correspond to rungs on "Jacob's
                Ladder of Density Functional Approximations"**
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of smd_solvent and pcm_dielectric may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 200)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            overwrite_inputs (dict): Dictionary of QChem input sections to add or overwrite variables.
                The currently available sections (keys) are rem, pcm,
                solvent, smx, opt, scan, van_der_waals, and plots. The value of each key is a
                dictionary of key value pairs relevant to that section. For example, to add
                a new variable to the rem section that sets symmetry to false, use

                overwrite_inputs = {"rem": {"symmetry": "false"}}

                **Note that if something like basis is added to the rem dict it will overwrite
                the default basis.**

                **Note that supplying a van_der_waals section here will automatically modify
                the PCM "radii" setting to "read".**

                **Note that all keys must be given as strings, even when they are numbers!**
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies only if
                you are using overwrite_inputs to add a $van_der_waals section to the input. In 'atomic' mode
                (default), dict keys represent the atomic number associated with each radius (e.g., '12' = carbon).
                In 'sequential' mode, dict keys represent the sequential position of a single
                specific atom in the input structure.
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles

        if scan_variables is None:
            raise ValueError("Cannot run a pes_scan job without some variable to scan over!")

        super().__init__(
            molecule=molecule,
            job_type="pes_scan",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            smd_solvent=smd_solvent,
            custom_smd=custom_smd,
            opt_variables=opt_variables,
            scan_variables=scan_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )
