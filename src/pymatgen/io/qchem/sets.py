"""Input sets for Qchem."""

from __future__ import annotations

import os
import warnings
from typing import TYPE_CHECKING

from monty.io import zopen

from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.utils import lower_and_check_unique

if TYPE_CHECKING:
    from typing import Any, Literal

    from pymatgen.core.structure import Molecule
    from pymatgen.util.typing import PathLike

__author__ = "Samuel Blau, Brandon Wood, Shyam Dwaraknath, Evan Spotte-Smith, Ryan Kingsbury"
__copyright__ = "Copyright 2018-2022, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"

# Note that in addition to the solvent-specific parameters, this dict contains
# dielectric constants for use with each solvent. The dielectric constants
# are used by the isodensity SS(V)PE electrostatic calculation part of CMIRS
# they are not part of the parameters tabulated by Q-Chem
# see https://manual.q-chem.com/latest/example_CMIRS-water.html
CMIRS_SETTINGS: dict[str, Any] = {
    "water": {
        "0.001": {
            "a": "-0.006736",
            "b": "0.032698",
            "c": "-1249.6",
            "d": "-21.405",
            "gamma": "3.7",
            "solvrho": "0.05",
        },
        "0.0005": {
            "a": "-0.006496",
            "b": "0.050833",
            "c": "-566.7",
            "d": "-30.503",
            "gamma": "3.2",
            "solvrho": "0.05",
        },
        "dielst": "78.39",
    },
    "benzene": {
        "0.001": {
            "a": "-0.00522",
            "b": "0.01294",
            "c": None,
            "d": None,
            "gamma": None,
            "solvrho": "0.0421",
        },
        "0.0005": {
            "a": "-0.00572",
            "b": "0.01116",
            "c": None,
            "d": None,
            "gamma": None,
            "solvrho": "0.0421",
        },
        "dielst": "2.28",
    },
    "cyclohexane": {
        "0.001": {
            "a": "-0.00938",
            "b": "0.03184",
            "c": None,
            "d": None,
            "gamma": None,
            "solvrho": "0.0396",
        },
        "0.0005": {
            "a": "-0.00721",
            "b": "0.05618",
            "c": None,
            "d": None,
            "gamma": None,
            "solvrho": "0.0396",
        },
        "dielst": "2.02",
    },
    "dimethyl sulfoxide": {
        "0.001": {
            "a": "-0.00951",
            "b": "0.044791",
            "c": None,
            "d": "-162.07",
            "gamma": "4.1",
            "solvrho": "0.05279",
        },
        "0.0005": {
            "a": "-0.002523",
            "b": "0.011757",
            "c": None,
            "d": "-817.93",
            "gamma": "4.3",
            "solvrho": "0.05279",
        },
        "dielst": "47",
    },
    "acetonitrile": {
        "0.001": {
            "a": "-0.008178",
            "b": "0.045278",
            "c": None,
            "d": "-0.33914",
            "gamma": "1.3",
            "solvrho": "0.03764",
        },
        "0.0005": {
            "a": "-0.003805",
            "b": "0.03223",
            "c": None,
            "d": "-0.44492",
            "gamma": "1.2",
            "solvrho": "0.03764",
        },
        "dielst": "36.64",
    },
}


class QChemDictSet(QCInput):
    """Build a QCInput given all the various input parameters. Can be extended by standard implementations below."""

    def __init__(
        self,
        molecule: Molecule,
        job_type: str,
        basis_set: str,
        scf_algorithm: str,
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        opt_variables: dict[str, list] | None = None,
        scan_variables: dict[str, list] | None = None,
        max_scf_cycles: int = 100,
        geom_opt_max_cycles: int = 200,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        geom_opt: dict | None = None,
        cdft_constraints: list[list[dict]] | None = None,
        almo_coupling_states: list[list[tuple[int, int]]] | None = None,
        overwrite_inputs: dict | None = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
        extra_scf_print: bool = False,
    ):
        """
        Args:
            molecule (Pymatgen Molecule object): Molecule to run QChem on.
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. For example, "def2-tzvpd".
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}

            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
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
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
                (Default: False)
            nbo_params (dict): A dict containing the desired NBO params. Note that a key:value pair of
                "version":7 will trigger NBO7 analysis. Otherwise, NBO5 analysis will be performed,
                including if an empty dict is passed. Besides a key of "version", all other key:value
                pairs will be written into the $nbo section of the QChem input file. (Default: False)
            geom_opt (dict): A dict containing parameters for the $geom_opt section of the Q-Chem input
                file, which control the new geometry optimizer available starting in version 5.4.2. The
                new optimizer remains under development but was officially released and became the default
                optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
                explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
                (Default: False)
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
                        cdft_constraints (list of lists of dicts):
                A list of lists of dictionaries, where each dictionary represents a charge
                constraint in the cdft section of the QChem input file.

                Each entry in the main list represents one state (allowing for multi-configuration
                calculations using constrained density functional theory - configuration interaction
                (CDFT-CI). Each state is represented by a list, which itself contains some number of
                constraints (dictionaries).

                Ex:

                1. For a single-state calculation with two constraints:
                 cdft_constraints=[[
                    {
                        "value": 1.0,
                        "coefficients": [1.0],
                        "first_atoms": [1],
                        "last_atoms": [2],
                        "types": [None]
                    },
                    {
                        "value": 2.0,
                        "coefficients": [1.0, -1.0],
                        "first_atoms": [1, 17],
                        "last_atoms": [3, 19],
                        "types": ["s"]
                    }
                ]]

                Note that a type of None will default to a charge constraint (which can also be
                accessed by requesting a type of "c" or "charge").

                2. For a CDFT-CI multi-reference calculation:
                cdft_constraints=[
                    [
                        {
                            "value": 1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ],
                    [
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": -1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ]
                ]
            cdft_constraints (list[list[dict]]): A list of lists of dictionaries, where each
            almo_coupling_states (list of lists of int 2-tuples):
                A list of lists of int 2-tuples used for calculations of diabatization and state
                coupling calculations relying on the absolutely localized molecular orbitals (ALMO)
                methodology. Each entry in the main list represents a single state (two states are
                included in an ALMO calculation). Within a single state, each 2-tuple represents the
                charge and spin multiplicity of a single fragment.
                ex: almo_coupling_states=[
                            [
                                (1, 2),
                                (0, 1)
                            ],
                            [
                                (0, 1),
                                (1, 2)
                            ]
                        ]
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
            extra_scf_print (bool): Whether to store extra information generated from the SCF
                cycle. If switched on, the Fock Matrix, coefficients of MO and the density matrix
                will be stored.
        """
        self.molecule = molecule
        self.job_type = job_type
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.qchem_version = qchem_version
        self.dft_rung = dft_rung
        self.pcm_dielectric = pcm_dielectric
        self.isosvp_dielectric = isosvp_dielectric
        self.smd_solvent = smd_solvent
        self.cmirs_solvent = cmirs_solvent
        self.custom_smd = custom_smd
        self.opt_variables = opt_variables
        self.scan_variables = scan_variables
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        self.plot_cubes = plot_cubes
        self.output_wavefunction = output_wavefunction
        self.nbo_params = nbo_params
        self.geom_opt = geom_opt
        self.cdft_constraints = cdft_constraints
        self.almo_coupling_states = almo_coupling_states
        self.overwrite_inputs = overwrite_inputs
        self.vdw_mode = vdw_mode
        self.extra_scf_print = extra_scf_print

        pcm_defaults = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1",
        }

        svp_defaults: dict[str, Any] = {
            "rhoiso": "0.001",
            "nptleb": "1202",
            "itrngr": "2",
            "irotgr": "2",
        }

        plots_defaults = {"grid_spacing": "0.05", "total_density": "0"}

        opt = {} if self.opt_variables is None else self.opt_variables

        scan = {} if self.scan_variables is None else self.scan_variables

        pcm: dict = {}
        solvent: dict = {}
        smx: dict = {}
        vdw: dict = {}
        plots: dict = {}
        svp: dict[str, Any] = {}
        pcm_nonels: dict[str, Any] = {}
        rem: dict[str, Any] = {
            "job_type": job_type,
            "basis": self.basis_set,
            "max_scf_cycles": str(self.max_scf_cycles),
            "gen_scfman": "true",
            "xc_grid": "3",
            "thresh": "14",
            "s2thresh": "16",
            "scf_algorithm": self.scf_algorithm,
            "resp_charges": "true",
            "symmetry": "false",
            "sym_ignore": "true",
        }

        if self.dft_rung == 1:
            rem["method"] = "spw92"
        elif self.dft_rung == 2:
            rem["method"] = "b97-d3"
            rem["dft_d"] = "d3_bj"
        elif self.dft_rung == 3:
            rem["method"] = "b97mv"
        elif self.dft_rung == 4:
            rem["method"] = "wb97mv"
        elif self.dft_rung == 5:
            rem["method"] = "wb97m(2)"
        else:
            raise ValueError("dft_rung should be between 1 and 5!")

        if self.job_type.lower() in ["opt", "ts", "pes_scan"]:
            rem["geom_opt_max_cycles"] = str(self.geom_opt_max_cycles)

        # To keep things simpler on the analysis side, don't give user option to change *.wfn file name
        if self.output_wavefunction:
            rem["write_wfn"] = "wavefunction"

        solvent_def = 0
        for a in [
            self.pcm_dielectric,
            self.isosvp_dielectric,
            self.smd_solvent,
            self.cmirs_solvent,
        ]:
            if a is not None:
                solvent_def += 1
        if solvent_def > 1:
            raise ValueError("Only one of PCM, ISOSVP, SMD, and CMIRSmay be used for solvation.")

        if self.pcm_dielectric is not None:
            pcm = pcm_defaults
            solvent["dielectric"] = str(self.pcm_dielectric)
            rem["solvent_method"] = "pcm"

        if self.isosvp_dielectric is not None:
            svp = svp_defaults
            svp["dielst"] = str(self.isosvp_dielectric)
            rem["solvent_method"] = "isosvp"
            rem["gen_scfman"] = "false"

        if self.smd_solvent is not None:
            if self.smd_solvent == "custom":
                smx["solvent"] = "other"
            else:
                smx["solvent"] = self.smd_solvent
            rem["solvent_method"] = "smd"
            rem["ideriv"] = "1"
            if self.smd_solvent in ("custom", "other"):
                if self.custom_smd is None:
                    raise ValueError(
                        "A user-defined SMD requires passing custom_smd, a string of seven comma separated values "
                        "in the following order: dielectric, refractive index, acidity, basicity, surface"
                        " tension, aromaticity, electronegative halogenicity"
                    )
                if self.qchem_version == 6:
                    custom_smd_vals = self.custom_smd.split(",")
                    smx["epsilon"] = custom_smd_vals[0]
                    smx["SolN"] = custom_smd_vals[1]
                    smx["SolA"] = custom_smd_vals[2]
                    smx["SolB"] = custom_smd_vals[3]
                    smx["SolG"] = custom_smd_vals[4]
                    smx["SolC"] = custom_smd_vals[5]
                    smx["SolH"] = custom_smd_vals[6]

        if self.cmirs_solvent is not None:
            # set up the ISOSVP calculation consistently with the CMIRS
            svp = svp_defaults
            rem["solvent_method"] = "isosvp"
            rem["gen_scfman"] = "false"
            svp["dielst"] = CMIRS_SETTINGS[self.cmirs_solvent]["dielst"]
            svp["idefesr"] = "1"  # this flag enables the CMIRS part
            svp["ipnrf"] = "1"  # this flag is also required for some undocumented reason
            pcm_nonels = CMIRS_SETTINGS[self.cmirs_solvent][svp["rhoiso"]]
            pcm_nonels["delta"] = "7"  # as recommended by Q-Chem. See manual.
            pcm_nonels["gaulag_n"] = "40"  # as recommended by Q-Chem. See manual.

        if self.plot_cubes:
            plots = plots_defaults
            rem["plots"] = "true"
            rem["make_cube_files"] = "true"

        nbo = self.nbo_params
        if self.nbo_params is not None:
            rem["nbo"] = "true"
            if "version" in self.nbo_params:
                if self.nbo_params["version"] == 7:
                    rem["nbo_external"] = "true"
                else:
                    raise RuntimeError("nbo params version should only be set to 7! Exiting...")
            nbo = {}
            for key in self.nbo_params:
                if key != "version":
                    nbo[key] = self.nbo_params[key]

        tmp_geom_opt = self.geom_opt
        geom_opt = self.geom_opt
        if (self.job_type.lower() in ["opt", "optimization"] and self.qchem_version == 6) or (
            self.qchem_version == 5 and self.geom_opt is not None
        ):
            if self.qchem_version == 5:
                rem["geom_opt2"] = "3"
            elif self.qchem_version == 6 and not self.geom_opt:
                tmp_geom_opt = {}
            if tmp_geom_opt is not None:
                if "maxiter" in tmp_geom_opt:
                    if tmp_geom_opt["maxiter"] != str(self.geom_opt_max_cycles):
                        raise RuntimeError("Max # of optimization cycles must be the same! Exiting...")
                else:
                    tmp_geom_opt["maxiter"] = str(self.geom_opt_max_cycles)
                if self.qchem_version == 6:
                    if "coordinates" not in tmp_geom_opt:
                        tmp_geom_opt["coordinates"] = "redundant"
                    if "max_displacement" not in tmp_geom_opt:
                        tmp_geom_opt["max_displacement"] = "0.1"
                    if "optimization_restart" not in tmp_geom_opt:
                        tmp_geom_opt["optimization_restart"] = "false"
                geom_opt = {}
                for key in tmp_geom_opt:
                    geom_opt[key] = tmp_geom_opt[key]

        if self.overwrite_inputs:
            for sec, sec_dict in self.overwrite_inputs.items():
                if sec == "rem":
                    rem |= lower_and_check_unique(sec_dict)
                if sec == "pcm":
                    pcm |= lower_and_check_unique(sec_dict)
                if sec == "solvent":
                    solvent |= lower_and_check_unique(sec_dict)
                    if rem["solvent_method"] != "pcm":
                        warnings.warn(
                            "The solvent section will be ignored unless solvent_method=pcm!",
                            stacklevel=2,
                        )
                if sec == "smx":
                    smx |= lower_and_check_unique(sec_dict)
                if sec == "scan":
                    scan |= lower_and_check_unique(sec_dict)
                if sec == "van_der_waals":
                    vdw |= lower_and_check_unique(sec_dict)
                    # set the PCM section to read custom radii
                    pcm["radii"] = "read"
                if sec == "plots":
                    plots |= lower_and_check_unique(sec_dict)
                if sec == "nbo":
                    raise RuntimeError("Set nbo parameters directly with nbo_params input! Exiting...")
                if sec == "geom_opt":
                    raise RuntimeError("Set geom_opt params directly with geom_opt input! Exiting...")
                if sec == "opt":
                    opt |= lower_and_check_unique(sec_dict)
                if sec == "svp":
                    temp_svp = lower_and_check_unique(sec_dict)
                    for k, v in temp_svp.items():
                        if k == "rhoiso" and self.cmirs_solvent is not None:
                            # must update both svp and pcm_nonels sections
                            if v not in ["0.001", "0.0005"]:
                                raise RuntimeError(
                                    "CMIRS is only parameterized for RHOISO values of 0.001 or 0.0005! Exiting..."
                                )
                            for k2 in pcm_nonels:
                                if CMIRS_SETTINGS[self.cmirs_solvent][v].get(k2):
                                    pcm_nonels[k2] = CMIRS_SETTINGS[self.cmirs_solvent][v].get(k2)
                        if k == "idefesr":
                            if self.cmirs_solvent is not None and v == "0":
                                warnings.warn(
                                    "Setting IDEFESR=0 will disable the CMIRS calculation you requested!",
                                    stacklevel=2,
                                )
                            if self.cmirs_solvent is None and v == "1":
                                warnings.warn(
                                    "Setting IDEFESR=1 will have no effect unless you specify a cmirs_solvent!",
                                    stacklevel=2,
                                )
                        if k == "dielst" and rem["solvent_method"] != "isosvp":
                            warnings.warn(
                                "Setting DIELST will have no effect unless you specify a solvent_method=isosvp!",
                                stacklevel=2,
                            )

                        svp[k] = v
                if sec == "pcm_nonels":
                    temp_pcm_nonels = lower_and_check_unique(sec_dict)
                    for k, v in temp_pcm_nonels.items():
                        pcm_nonels[k] = v

        if extra_scf_print:
            # Allow for the printing of the Fock matrix and the eigenvales
            rem["scf_final_print"] = "3"
            # If extra_scf_print is specified, make sure that the convergence of the
            # SCF cycle is at least 1e-8. Anything less than that might not be appropriate
            # for printing out the Fock Matrix and coefficients of the MO.
            if "scf_convergence" not in rem or int(rem["scf_convergence"]) < 8:
                rem["scf_convergence"] = "8"

        super().__init__(
            self.molecule,
            rem=rem,
            opt=opt,
            pcm=pcm,
            solvent=solvent,
            smx=smx,
            scan=scan,
            van_der_waals=vdw,
            vdw_mode=self.vdw_mode,
            plots=plots,
            nbo=nbo,
            geom_opt=geom_opt,
            cdft=self.cdft_constraints,
            almo_coupling=self.almo_coupling_states,
            svp=svp,
            pcm_nonels=pcm_nonels,
        )

    def write(self, input_file: PathLike) -> None:
        """
        Args:
            input_file (PathLike): Filename.
        """
        self.write_file(input_file)
        if self.smd_solvent in {"custom", "other"} and self.qchem_version == 5:
            with zopen(os.path.join(os.path.dirname(input_file), "solvent_data"), mode="wt", encoding="utf-8") as file:
                file.write(self.custom_smd)  # type:ignore[arg-type]


class SinglePointSet(QChemDictSet):
    """QChemDictSet for a single point calculation."""

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
        cdft_constraints: list[list[dict]] | None = None,
        almo_coupling_states: list[list[tuple[int, int]]] | None = None,
        extra_scf_print: bool = False,
        overwrite_inputs: dict | None = None,
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. (Default: "def2-tzvpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters. Note that due to limitations in Q-Chem, use of the ISOSVP
                or CMIRS solvent models will disable the GEN_SCFMAN algorithm, which may limit compatible choices
                for scf_algorithm.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics. Note also that due to limitations in Q-Chem, use of the ISOSVP
                or CMIRS solvent models will disable the GEN_SCFMAN algorithm, which may limit compatible choices
                for scf_algorithm.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
                (Default: False)
            cdft_constraints (list of lists of dicts):
                A list of lists of dictionaries, where each dictionary represents a charge
                constraint in the cdft section of the QChem input file.

                Each entry in the main list represents one state (allowing for multi-configuration
                calculations using constrained density functional theory - configuration interaction
                (CDFT-CI). Each state is represented by a list, which itself contains some number of
                constraints (dictionaries).

                Ex:

                1. For a single-state calculation with two constraints:
                 cdft_constraints=[[
                    {
                        "value": 1.0,
                        "coefficients": [1.0],
                        "first_atoms": [1],
                        "last_atoms": [2],
                        "types": [None]
                    },
                    {
                        "value": 2.0,
                        "coefficients": [1.0, -1.0],
                        "first_atoms": [1, 17],
                        "last_atoms": [3, 19],
                        "types": ["s"]
                    }
                ]]

                Note that a type of None will default to a charge constraint (which can also be
                accessed by requesting a type of "c" or "charge").

                2. For a CDFT-CI multi-reference calculation:
                cdft_constraints=[
                    [
                        {
                            "value": 1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ],
                    [
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": -1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ]
                ]
            almo_coupling_states (list of lists of int 2-tuples):
                A list of lists of int 2-tuples used for calculations of diabatization and state
                coupling calculations relying on the absolutely localized molecular orbitals (ALMO)
                methodology. Each entry in the main list represents a single state (two states are
                included in an ALMO calculation). Within a single state, each 2-tuple represents the
                charge and spin multiplicity of a single fragment.
                ex: almo_coupling_states=[
                            [
                                (1, 2),
                                (0, 1)
                            ],
                            [
                                (0, 1),
                                (1, 2)
                            ]
                        ]
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
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
            extra_scf_print (bool): Whether to store extra information generated from the SCF
                cycle. If switched on, the Fock Matrix, coefficients of MO and the density matrix
                will be stored.
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="sp",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            vdw_mode=vdw_mode,
            cdft_constraints=cdft_constraints,
            almo_coupling_states=almo_coupling_states,
            overwrite_inputs=overwrite_inputs,
            extra_scf_print=extra_scf_print,
        )


class OptSet(QChemDictSet):
    """QChemDictSet for a geometry optimization."""

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-svpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        opt_variables: dict[str, list] | None = None,
        geom_opt_max_cycles: int = 200,
        geom_opt: dict | None = None,
        cdft_constraints: list[list[dict]] | None = None,
        overwrite_inputs: dict | None = None,
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            job_type (str): QChem job type to run. Valid options are "opt" for optimization,
                "sp" for single point, "freq" for frequency calculation, or "force" for
                force evaluation.
            basis_set (str): Basis set to use. (Default: "def2-svpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            geom_opt (dict): A dict containing parameters for the $geom_opt section of the Q-Chem input
                file, which control the new geometry optimizer available starting in version 5.4.2. The
                new optimizer remains under development but was officially released and became the default
                optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
                explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
                (Default: False)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
                (Default: False)
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
            cdft_constraints (list of lists of dicts):
                A list of lists of dictionaries, where each dictionary represents a charge
                constraint in the cdft section of the QChem input file.

                Each entry in the main list represents one state (allowing for multi-configuration
                calculations using constrained density functional theory - configuration interaction
                (CDFT-CI). Each state is represented by a list, which itself contains some number of
                constraints (dictionaries).

                Ex:

                1. For a single-state calculation with two constraints:
                 cdft_constraints=[[
                    {
                        "value": 1.0,
                        "coefficients": [1.0],
                        "first_atoms": [1],
                        "last_atoms": [2],
                        "types": [None]
                    },
                    {
                        "value": 2.0,
                        "coefficients": [1.0, -1.0],
                        "first_atoms": [1, 17],
                        "last_atoms": [3, 19],
                        "types": ["s"]
                    }
                ]]

                Note that a type of None will default to a charge constraint (which can also be
                accessed by requesting a type of "c" or "charge").

                2. For a CDFT-CI multi-reference calculation:
                cdft_constraints=[
                    [
                        {
                            "value": 1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ],
                    [
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": -1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ]
                ]
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
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            opt_variables=opt_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            geom_opt_max_cycles=self.geom_opt_max_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            geom_opt=geom_opt,
            cdft_constraints=cdft_constraints,
            overwrite_inputs=overwrite_inputs,
        )


class TransitionStateSet(QChemDictSet):
    """QChemDictSet for a transition-state search."""

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-svpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        opt_variables: dict[str, list] | None = None,
        geom_opt_max_cycles: int = 200,
        geom_opt: dict | None = None,
        overwrite_inputs: dict | None = None,
        vdw_mode="atomic",
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-svpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            geom_opt_max_cycles (int): Maximum number of geometry optimization iterations. (Default: 200)
            geom_opt (dict): A dict containing parameters for the $geom_opt section of the Q-Chem input
                file, which control the new geometry optimizer available starting in version 5.4.2. The
                new optimizer remains under development but was officially released and became the default
                optimizer in Q-Chem version 6.0.0. Note that for version 5.4.2, the new optimizer must be
                explicitly requested by passing in a dictionary (empty or otherwise) for this input parameter.
                (Default: False)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
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
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        self.geom_opt_max_cycles = geom_opt_max_cycles
        super().__init__(
            molecule=molecule,
            job_type="ts",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            opt_variables=opt_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            geom_opt_max_cycles=self.geom_opt_max_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            geom_opt=geom_opt,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )


class ForceSet(QChemDictSet):
    """QChemDictSet for a force (gradient) calculation."""

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-tzvpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
        cdft_constraints: list[list[dict]] | None = None,
        overwrite_inputs: dict | None = None,
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-tzvpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
                (Default: False)
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
            cdft_constraints (list of lists of dicts):
                A list of lists of dictionaries, where each dictionary represents a charge
                constraint in the cdft section of the QChem input file.

                Each entry in the main list represents one state (allowing for multi-configuration
                calculations using constrained density functional theory - configuration interaction
                (CDFT-CI). Each state is represented by a list, which itself contains some number of
                constraints (dictionaries).

                Ex:

                1. For a single-state calculation with two constraints:
                 cdft_constraints=[[
                    {
                        "value": 1.0,
                        "coefficients": [1.0],
                        "first_atoms": [1],
                        "last_atoms": [2],
                        "types": [None]
                    },
                    {
                        "value": 2.0,
                        "coefficients": [1.0, -1.0],
                        "first_atoms": [1, 17],
                        "last_atoms": [3, 19],
                        "types": ["s"]
                    }
                ]]

                Note that a type of None will default to a charge constraint (which can also be
                accessed by requesting a type of "c" or "charge").

                2. For a CDFT-CI multi-reference calculation:
                cdft_constraints=[
                    [
                        {
                            "value": 1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ],
                    [
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": -1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ]
                ]
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
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="force",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            vdw_mode=vdw_mode,
            cdft_constraints=cdft_constraints,
            overwrite_inputs=overwrite_inputs,
        )


class FreqSet(QChemDictSet):
    """QChemDictSet for a frequency calculation."""

    def __init__(
        self,
        molecule: Molecule,
        basis_set: str = "def2-svpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
        nbo_params: dict | None = None,
        vdw_mode: Literal["atomic", "sequential"] = "atomic",
        cdft_constraints: list[list[dict]] | None = None,
        overwrite_inputs: dict | None = None,
    ):
        """
        Args:
            molecule (Pymatgen Molecule object)
            basis_set (str): Basis set to use. (Default: "def2-svpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2).

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
                (Default: False)
            vdw_mode ('atomic' | 'sequential'): Method of specifying custom van der Waals radii. Applies
                only if you are using overwrite_inputs to add a $van_der_waals section to the input.
                In 'atomic' mode (default), dict keys represent the atomic number associated with each
                radius (e.g., '12' = carbon). In 'sequential' mode, dict keys represent the sequential
                position of a single specific atom in the input structure.
            cdft_constraints (list of lists of dicts):
                A list of lists of dictionaries, where each dictionary represents a charge
                constraint in the cdft section of the QChem input file.

                Each entry in the main list represents one state (allowing for multi-configuration
                calculations using constrained density functional theory - configuration interaction
                (CDFT-CI). Each state is represented by a list, which itself contains some number of
                constraints (dictionaries).

                Ex:

                1. For a single-state calculation with two constraints:
                 cdft_constraints=[[
                    {
                        "value": 1.0,
                        "coefficients": [1.0],
                        "first_atoms": [1],
                        "last_atoms": [2],
                        "types": [None]
                    },
                    {
                        "value": 2.0,
                        "coefficients": [1.0, -1.0],
                        "first_atoms": [1, 17],
                        "last_atoms": [3, 19],
                        "types": ["s"]
                    }
                ]]

                Note that a type of None will default to a charge constraint (which can also be
                accessed by requesting a type of "c" or "charge").

                2. For a CDFT-CI multi-reference calculation:
                cdft_constraints=[
                    [
                        {
                            "value": 1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ],
                    [
                        {
                            "value": 0.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["c"]
                        },
                        {
                            "value": -1.0,
                            "coefficients": [1.0],
                            "first_atoms": [1],
                            "last_atoms": [27],
                            "types": ["s"]
                        },
                    ]
                ]
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
        """
        self.basis_set = basis_set
        self.scf_algorithm = scf_algorithm
        self.max_scf_cycles = max_scf_cycles
        super().__init__(
            molecule=molecule,
            job_type="freq",
            dft_rung=dft_rung,
            pcm_dielectric=pcm_dielectric,
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            vdw_mode=vdw_mode,
            cdft_constraints=cdft_constraints,
            overwrite_inputs=overwrite_inputs,
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
        basis_set: str = "def2-svpd",
        scf_algorithm: str = "diis",
        qchem_version: int = 5,
        dft_rung: int = 4,
        pcm_dielectric: float | None = None,
        isosvp_dielectric: float | None = None,
        smd_solvent: str | None = None,
        cmirs_solvent: (Literal["water", "acetonitrile", "dimethyl sulfoxide", "cyclohexane", "benzene"] | None) = None,
        custom_smd: str | None = None,
        max_scf_cycles: int = 100,
        plot_cubes: bool = False,
        output_wavefunction: bool = False,
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
            basis_set (str): Basis set to use. (Default: "def2-svpd")
            scf_algorithm (str): Algorithm to use for converging the SCF. Recommended choices are
                "DIIS", "GDM", and "DIIS_GDM". Other algorithms supported by Qchem's GEN_SCFMAN
                module will also likely perform well. Refer to the QChem manual for further details.
                (Default: "diis")
            qchem_version (int): Which major version of Q-Chem will be run. Supports 5 and 6. (Default: 5)
            dft_rung (int): Select the rung on "Jacob's Ladder of Density Functional Approximations" in
                order of increasing accuracy/cost. For each rung, we have prescribed one functional based
                on our experience, available benchmarks, and the suggestions of the Q-Chem manual:
                1 (LSDA) = SPW92
                2 (GGA) = B97-D3(BJ)
                3 (metaGGA) = B97M-V
                4 (hybrid metaGGA) = ωB97M-V
                5 (double hybrid metaGGA) = ωB97M-(2)

                (Default: 4)

                To set a functional not given by one of the above, set the overwrite_inputs
                argument to {"method":"<NAME OF FUNCTIONAL>"}
            pcm_dielectric (float): Dielectric constant to use for PCM implicit solvation model. (Default: None)
                If supplied, will set up the $pcm section of the input file for a C-PCM calculation.
                Other types of PCM calculations (e.g., IEF-PCM, SS(V)PE, etc.) may be requested by passing
                custom keywords to overwrite_inputs, e.g.
                overwrite_inputs = {"pcm": {"theory": "ssvpe"}}
                Refer to the QChem manual for further details on the models available.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            isosvp_dielectric (float): Dielectric constant to use for isodensity SS(V)PE implicit solvation model.
                (Default: None). If supplied, will set solvent_method to "isosvp" and populate the $svp section
                of the input file with appropriate parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            smd_solvent (str): Solvent to use for SMD implicit solvation model. (Default: None)
                Examples include "water", "ethanol", "methanol", and "acetonitrile". Refer to the QChem
                manual for a complete list of solvents available. To define a custom solvent, set this
                argument to "custom" and populate custom_smd with the necessary parameters.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            cmirs_solvent (str): Solvent to use for the CMIRS implicit solvation model. (Default: None).
                Only 5 solvents are presently available as of Q-Chem 6: "water", "benzene", "cyclohexane",
                "dimethyl sulfoxide", and "acetonitrile". Note that selection of a solvent here will also
                populate the iso SS(V)PE dielectric constant, because CMIRS uses the isodensity SS(V)PE model
                to compute electrostatics.

                **Note that only one of pcm_dielectric, isosvp_dielectric, smd_solvent, or cmirs_solvent may be set.**
            custom_smd (str): List of parameters to define a custom solvent in SMD. (Default: None)
                Must be given as a string of seven comma separated values in the following order:
                "dielectric, refractive index, acidity, basicity, surface tension, aromaticity,
                electronegative halogenicity"
                Refer to the QChem manual for further details.
            max_scf_cycles (int): Maximum number of SCF iterations. (Default: 100)
            plot_cubes (bool): Whether to write CUBE files of the electron density. (Default: False)
            output_wavefunction (bool): Whether to write a wavefunction (*.wfn) file of the electron density
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
            isosvp_dielectric=isosvp_dielectric,
            smd_solvent=smd_solvent,
            cmirs_solvent=cmirs_solvent,
            custom_smd=custom_smd,
            opt_variables=opt_variables,
            scan_variables=scan_variables,
            basis_set=self.basis_set,
            scf_algorithm=self.scf_algorithm,
            qchem_version=qchem_version,
            max_scf_cycles=self.max_scf_cycles,
            plot_cubes=plot_cubes,
            output_wavefunction=output_wavefunction,
            nbo_params=nbo_params,
            overwrite_inputs=overwrite_inputs,
            vdw_mode=vdw_mode,
        )
