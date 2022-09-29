# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
Input sets for Qchem
"""

from __future__ import annotations

import logging
import os
from typing import Literal

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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        opt_variables: dict[str, list] | None = None,
        scan_variables: dict[str, list] | None = None,
        max_scf_cycles: int = 100,
        geom_opt_max_cycles: int = 200,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        new_geom_opt: dict | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            nbo_params (dict): A dict containing the desired NBO params. Note that a key:value pair of
                "version":7 will trigger NBO7 analysis. Otherwise, NBO5 analysis will be performed,
                including if an empty dict is passed. Besides a key of "version", all other key:value
                pairs will be written into the $nbo section of the QChem input file. (Default: False)
            new_geom_opt (dict): A dict containing parameters for the $geom_opt section of the QChem
                input file, which control the new geometry optimizer available starting in version 5.4.2.
                Note that the new optimizer remains under development and not officially released.
                Further note that even passig an empty dictionary will trigger the new optimizer.
                (Default: False)
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
        self.new_geom_opt = new_geom_opt
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
        myrem["thresh"] = "14"
        myrem["s2thresh"] = "16"
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

        mynbo = self.nbo_params
        if self.nbo_params is not None:
            myrem["nbo"] = "true"
            if "version" in self.nbo_params:
                if self.nbo_params["version"] == 7:
                    myrem["nbo_external"] = "true"
                else:
                    raise RuntimeError("nbo params version should only be set to 7! Exiting...")
            mynbo = {}
            for key in self.nbo_params:
                if key != "version":
                    mynbo[key] = self.nbo_params[key]

        my_geom_opt = self.new_geom_opt
        if self.new_geom_opt is not None:
            myrem["geom_opt2"] = "3"
            if "maxiter" in self.new_geom_opt:
                if self.new_geom_opt["maxiter"] != str(self.geom_opt_max_cycles):
                    raise RuntimeError("Max # of optimization cycles must be the same! Exiting...")
            else:
                self.new_geom_opt["maxiter"] = str(self.geom_opt_max_cycles)
            my_geom_opt = {}
            for key in self.new_geom_opt:
                my_geom_opt[key] = self.new_geom_opt[key]

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
                    if mynbo is None:
                        raise RuntimeError("Can't overwrite nbo params when NBO is not being run! Exiting...")
                    temp_nbo = lower_and_check_unique(sec_dict)
                    for k, v in temp_nbo.items():
                        mynbo[k] = v
                if sec == "geom_opt":
                    if my_geom_opt is None:
                        raise RuntimeError(
                            "Can't overwrite geom_opt params when not using the new optimizer! Exiting..."
                        )
                    temp_geomopt = lower_and_check_unique(sec_dict)
                    for k, v in temp_geomopt.items():
                        my_geom_opt[k] = v
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
            nbo=mynbo,
            geom_opt=my_geom_opt,
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        opt_variables: dict[str, list] | None = None,
        geom_opt_max_cycles: int = 200,
        new_geom_opt: dict | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
            new_geom_opt=new_geom_opt,
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        opt_variables: dict[str, list] | None = None,
        geom_opt_max_cycles: int = 200,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
        pcm_dielectric: float | None = None,
        smd_solvent: str | None = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        nbo_params: dict | None = None,
        opt_variables: dict[str, list] | None = None,
        scan_variables: dict[str, list] | None = None,
        overwrite_inputs: dict | None = None,
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
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
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
