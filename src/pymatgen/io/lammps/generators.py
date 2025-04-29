"""
Input set generators for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
(https://github.com/Matgenix/atomate2-lammps).
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from string import Template

import numpy as np
from monty.io import zopen
from monty.json import MSONable

from pymatgen.core import Structure
from pymatgen.io.core import InputGenerator
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile
from pymatgen.io.lammps.sets import LammpsInputSet
from pymatgen.util.typing import PathLike

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = f"{MODULE_DIR}/templates"
_BASE_LAMMPS_SETTINGS = {
    "units": "metal",
    "atom_style": "atomic",
    "dimension": 3,
    "boundary": ("p", "p", "p"),
    "pair_style": "lj/cut 10.0",
    "thermo": 100,
    "start_temp": 300,
    "end_temp": 300,
    "start_pressure": 0,
    "end_pressure": 0,
    "timestep": 0.001,
    "friction": 0.1,
    "log_interval": 100,
    "traj_interval": 100,
    "ensemble": "nvt",
    "thermostat": "nose-hoover",
    "barostat": None,
    "nsteps": 1000,
    "restart": "",
    "tol": 1e-6,
    "min_style": "cg",
}

FF_STYLE_KEYS = ["pair_style", "bond_style", "angle_style", "dihedral_style", "improper_style", "kspace_style"]
FF_COEFF_KEYS = ["pair_coeff", "bond_coeff", "angle_coeff", "dihedral_coeff", "improper_coeff", "kspace_coeff"]

LAMMPS_DEFINED_TYPES: dict[str, set[str | None]] = {
    "units": {"metal", "lj", "real", "si", "cgs", "electron", "micro", "nano"},
    "atom_style": {"atomic", "angle", "body", "bond", "charge", "electron", "full", "molecular"},
    "boundary": {"p", "f", "s", "m", "fs", "fm"},
    "ensemble": {"nve", "nvt", "npt", "nph", "minimize"},
    "thermostat": {"nose-hoover", "langevin", None},
    "barostat": {"nose-hoover", "berendsen", "langevin", None},
    "min_style": {"cg", "sd", "fire", "hftn", "quickmin", "spin", "spin/cg", "spin/lbfgs"},
}


@dataclass
class LammpsSettings(MSONable):
    """Define schema for LAMMPS convenience input class.

    The idea is to dynamically set the attributes of this class
    based on the LAMMPS input file template and the settings
    dictionary provided by the user. This way, standard inputs
    that do not typically vary (e.g., units, atom style, dimension,
    boundary conditions, styles) can be validated and set
    while allowing the user to specify the rest of the settings in
    a flexible manner (which won't be explicitly validated here).

    Args:
    units : str
        LAMMPS units style.
    atom_style : str
        LAMMPS atom style.
    dimension : int
        LAMMPS box dimension.
    boundary : tuple[str,str,str]
        Box boundary conditions. Each boundary can be 'p', 'f', 's', 'm', 'fs', 'fm'.
    pair_style : str
        FF pair style type. Default is 'lj/cut 10.0'.
    bond_style : str
        FF bond style type.
    angle_style : str
        FF angle style type.
    dihedral_style : str
        FF dihedral style type.
    start_temp : float
        Initial temperature. Default is 300 K. Initial velocities are generated based on this temperature.
    end_temp : float
        Final temperature. Default is 300 K.
    start_pressure : float | list | np.ndarray
        Initial pressure. Default is 0 atm. Can be a single value or a list for anisotropic pressure.
    end_pressure : float | list | np.ndarray
        Final pressure. Default is 0 atm. Can be a single value or a list for anisotropic pressure.
    timestep : float
        Simulation timestep. Default is 0.001 ps.
    friction : float
        Thermostat/Barostat friction coefficient. Default is 0.1 ps^-1.
    log_interval : int
        Interval for logging thermodynamic data (energy, temperature, pressure, volume).
        Default is 100 steps.
    traj_interval : int
        Interval for writing trajectory data (positions, velocities, forces) to a dump file.
        Default is 100 steps.
    ensemble : str
        Simulation ensemble. Default is 'nvt'.
    thermostat : str
        Thermostat type. Default is 'nose-hoover'.
    barostat : str
        Barostat type. Default is 'nose-hoover'.
    nsteps : int
        Number of simulation steps. Default is 1000. This is also the maximum number of
        steps for minimization simulations.
    """

    ''' units : Literal["metal", "lj", "real", "si", "cgs", "electron", "micro", "nano"] = "metal"
    atom_style : Literal[
        "atomic", "angle", "body", "bond", "charge", "electron", "full", "molecular"
    ] = "atomic"
    dimension : int = 3
    boundary : tuple[str,str,str] = ("p", "p", "p")
    pair_style : str = "lj/cut 10.0"
    bond_style : str = None
    angle_style : str = None
    dihedral_style : str = None
    start_temp : float = 300.
    end_temp : float = 300.
    start_pressure : float | list | np.ndarray = 0.
    end_pressure : float | list | np.ndarray = 0.
    timestep : float = 0.001
    friction : float = 0.1
    log_interval : int = 100
    traj_interval : int = 100
    ensemble : Literal["nve","nvt","npt","nph","minimize"] = "nvt"
    thermostat : Literal["nose-hoover", "langevin", None] = "nose-hoover"
    barostat : Literal["nose-hoover", "berendsen", "langevin", None] = "nose-hoover"
    nsteps : int = 1000
    restart : str = None
    tol : float = 1e-6
    min_style : Literal["cg", "sd", "fire", "hftn", "quickmin", "spin"] = "cg"'''

    validate_params: bool = field(default=True)

    def __init__(self, validate_params: bool = True, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        if validate_params:
            self.__post_init__()

    def __post_init__(self) -> None:
        """Validate input values."""

        for attr, accept_vals in LAMMPS_DEFINED_TYPES.items():
            curr_val = getattr(self, attr, None)
            if isinstance(curr_val, list | tuple):
                is_ok = all(v in accept_vals for v in curr_val)
            else:
                is_ok = curr_val in accept_vals
            if not is_ok:
                raise ValueError(f"Error validating key {attr}: set to {curr_val}, should be one of {accept_vals}.")

        if self.restart and not isinstance(self.restart, str):
            raise ValueError(
                f"restart should be the path to the restart file from the previous run, got {self.restart}."
            )

        if self.ensemble not in ["nve", "minimize"]:
            if isinstance(self.start_pressure, (list, np.ndarray)) and len(self.start_pressure) != 3:
                raise ValueError(f"start_pressure should be a list of 3 values, got {self.start_pressure}.")
            if isinstance(self.end_pressure, (list, np.ndarray)) and len(self.end_pressure) != 3:
                raise ValueError(f"end_pressure should be a list of 3 values, got {self.end_pressure}.")

            if (self.thermostat or self.barostat) and self.friction < self.timestep:
                warnings.warn(
                    f"Friction ({self.friction}) is smaller than the timestep ({self.timestep}). "
                    "This may lead to numerical instability.",
                    stacklevel=2,
                )

        if self.ensemble == "minimize":
            if self.nsteps < 1:
                raise ValueError(f"nsteps should be greater than 0 for minimization simulations, got {self.nsteps}.")
            if self.tol > 1e-4:
                warnings.warn(
                    f"Tolerance for minimization ({self.tol}) is larger than 1e-4. "
                    "This may lead to inaccurate results.",
                    stacklevel=2,
                )

    @property
    def dict(self) -> dict:
        return self.__dict__


@dataclass
class BaseLammpsSetGenerator(InputGenerator):
    """
    Base class for generating LAMMPS input sets.

    Args:
        inputfile : LammpsInputFile | str | Path
            Premade input file for the LAMMPS simulation.
            Useful if the user wants to use a custom input file
            (to make use of Lammps' flexibility).
            Default format based on the md.template file in the templates directory.
        settings : dict | LammpsInputSettings
            Settings for the LAMMPS simulation. Default settings are given in the
            _BASE_LAMMPS_SETTINGS object in 'metal' units for reference.
        force_field : Path | dict | ForceField
            Force field file or dictionary containing the force field parameters.
            Default is None.
        calc_type : str
            Type of calculation to be performed by LAMMPS.
        keep_stages : bool
            Whether to keep the stages of the input file or not. Default is True.
        override_updates : bool
            Whether to override the updates to the input file, i.e.,
            keep the input file as is. Default is False.
    """

    inputfile: LammpsInputFile | str | PathLike = field(default=None)
    settings: dict | LammpsSettings = field(default_factory=dict)
    force_field: Path | dict = field(default=None)
    calc_type: str = field(default="lammps")
    include_defaults: bool = field(default=True)
    validate_params: bool = field(default=True)
    keep_stages: bool = field(default=True)

    def __post_init__(self):
        """Determine correct template requested by user.

        If not inputfile, we assume the user is going to use the default template for
        one the defined flows (NVT/NVE/NPT, Minimize, MeltQuench, etc.)
        Instead, if the user/another flow specifies a template, we'll use that and an
        optional settings dict with it and not bother with validating the inputs in the
        inputfile due to the amount of flexibility LAMMPS allows
        If another flow is used that defines it's own inputfile template, we'll leave the
        validation of the user's inputs to that flow's init method.

        """
        if self.inputfile is None:
            ensemble = (
                self.settings.get("ensemble", "nvt") if isinstance(self.settings, dict) else self.settings.ensemble
            )

            file: str | None = None
            if ensemble in ["nve", "nvt", "npt", "nph"]:
                file = "md.template"
            elif ensemble in ["minimize"]:
                file = "minimization.template"
            if not file:
                raise ValueError(
                    f"Unknown {ensemble=}; acceptable values are\n{', '.join(LAMMPS_DEFINED_TYPES['ensemble'])}"
                )

            self.inputfile = LammpsInputFile.from_file(os.path.join(TEMPLATE_DIR, file), keep_stages=self.keep_stages)

        if isinstance(self.settings, dict):
            settings = self.settings.copy()
            if self.include_defaults:
                base_settings = _BASE_LAMMPS_SETTINGS.copy()
                for key, value in base_settings.items():
                    if key not in settings:
                        settings[key] = value
            self.settings = LammpsSettings(validate_params=self.validate_params, **settings)

    def update_settings(self, updates: dict):
        """
        Update the settings for the LammpsSettings object.
        Args:
            updates : dict
                Dictionary containing the settings to update.
        """
        if isinstance(self.settings, LammpsSettings):
            present_settings = self.settings.dict
            for k, v in updates.items():
                present_settings.update({k: v})
            self.settings = LammpsSettings(validate_params=self.validate_params, **present_settings)
        else:
            self.settings.update(updates)

    def get_input_set(
        self,
        data: Structure | LammpsData | CombinedData,
        additional_data: LammpsData | CombinedData | None = None,
        **kwargs,
    ) -> LammpsInputSet:
        """
        Generate a LAMMPS input set.
        Args:
            structure : Structure | LammpsData
                Structure or LammpsData object for the simulation.
            **kwargs : dict
                Additional keyword arguments to pass to the InputSet from pmg.
        """

        if not self.force_field:
            warnings.warn(
                "Force field not specified! Ensure you have the correct force field parameters "
                "in the data file/settings or will specify it manually using "
                "maker.input_set_generator.force_field.",
                stacklevel=2,
            )

        if isinstance(self.inputfile, PathLike):
            try:
                self.inputfile = LammpsInputFile.from_file(self.inputfile, keep_stages=self.keep_stages)
            except FileNotFoundError:
                try:
                    self.inputfile = LammpsInputFile.from_str(self.inputfile, keep_stages=self.keep_stages)
                except ValueError:
                    raise FileNotFoundError(
                        f"Input file {self.inputfile} not found. It was neither a path "
                        "nor a string repr of the inputfile. Please check your inputs!"
                    )

        settings_dict = self.settings.dict.copy() if isinstance(self.settings, LammpsSettings) else self.settings
        atom_style = settings_dict.get("atom_style", "full")
        print(f"Generating LAMMPS input set with settings: {settings_dict}")

        # Attempt to read data file and convert to LammpsData object
        if isinstance(data, Path):
            try:
                data = LammpsData.from_file(data, atom_style=atom_style)
            except FileNotFoundError:
                raise FileNotFoundError(f"Data file {data} not found. Please check the path.")

        if isinstance(data, str):
            try:
                data = LammpsData.from_str(data, atom_style=atom_style)
            except ValueError:
                raise ValueError(f"Data file {data} not recognized. Please check the format.")

        species = ""
        if isinstance(data, Structure):
            species = " ".join({s.symbol for s in data.species})
            data = LammpsData.from_structure(data, atom_style=atom_style)
            warnings.warn("Structure provided, converting to LammpsData object.", stacklevel=2)

        # Accounts for the restart file
        if settings_dict.get("restart", None):
            settings_dict.update(
                {"read_restart": f"{settings_dict['restart']}", "restart_flag": "read_restart", "read_data_flag": "#"}
            )

        # Housekeeping to fill up the default settings for the MD template
        settings_dict.update({f"{sys}_flag": "#" for sys in ["nve", "nvt", "npt", "nph", "restart", "extra_data"]})
        settings_dict.update({"read_data_flag": "read_data", "psymm": "iso"})
        # If the ensemble is not 'minimize', we set the read_data_flag to read_data
        # Convert start and end pressure to string if they are lists or arrays, and set psymm to accordingly
        if isinstance(settings_dict.get("start_pressure", None), (list, np.ndarray)):
            settings_dict.update(
                {"start_pressure": " ".join(map(str, settings_dict["start_pressure"])), "psymm": "aniso"}
            )
        if isinstance(settings_dict.get("end_pressure", None), (list, np.ndarray)):
            settings_dict.update({"end_pressure": " ".join(map(str, settings_dict["end_pressure"])), "psymm": "aniso"})

        # Loop over the LammpsSettings object and update the settings dictionary
        for attr, val in self.settings.dict.items():  # type: ignore[union-attr]
            if attr == "boundary":
                settings_dict.update({"boundary": " ".join(list(val))})

            elif attr == "ensemble":
                settings_dict.update({f"{val}_flag": "fix"})

            elif attr == "thermostat":
                if val == "langevin":
                    settings_dict.update({"nve_flag": "fix", "thermseed": 42, "thermostat": "langevin"})
                if val == "nose-hoover":
                    settings_dict.update({"thermostat": "nvt temp", "thermseed": ""})

            elif attr == "barostat":
                if val == "nose-hoover":
                    settings_dict.update(
                        {
                            "barostat": "npt temp",
                            "start_temp_barostat": settings_dict["start_temp"],
                            "end_temp_barostat": settings_dict["end_temp"],
                            "tfriction_barostat": settings_dict["friction"],
                        }
                    )

                if val in ["berendsen", "langevin"]:
                    settings_dict.update(
                        {
                            "barostat": "langevin" if val == "langevin" else "press/berendsen",
                            "nve_flag": "fix",
                            "nvt_flag": "fix",
                            "start_temp_barostat": "",
                            "end_temp_barostat": "",
                            "tfriction_barostat": "",
                            "thermostat": f"temp/{val}",
                        }
                    )
                    settings_dict.update({"thermoseed": 42 if val == "langevin" else ""})

            elif attr == "friction":
                settings_dict.update({"tfriction": val, "pfriction": val})

            else:
                settings_dict.update({attr: val})

        # Handle the force field input by writing a separate FF file
        # and making the necessary updates to the settings dict
        FF_string = ""
        if isinstance(self.force_field, str):
            FF_string += self.force_field
            settings_dict.update({f"{ff}_flag": "#" for ff in FF_STYLE_KEYS})

        if isinstance(self.force_field, dict):
            for key, value in self.force_field.items():
                if key in FF_STYLE_KEYS and value:
                    settings_dict.update({f"{key}": value, f"{key}_flag": f"{key}"})
                if key in FF_COEFF_KEYS and value:
                    FF_string += f"{key} {value}\n"
                if key in ["species"]:
                    # makes species specified in FF dict take precedence
                    species = " ".join(value) if isinstance(value, list) else value
                else:
                    warnings.warn(f"Force field key {key} not recognized, will be ignored.", stacklevel=2)

            for ff_key in FF_STYLE_KEYS:
                if ff_key not in self.settings.dict or not self.settings.dict[ff_key]:  # type: ignore[union-attr]
                    settings_dict.update({f"{ff_key}_flag": "#"})
                    warnings.warn(f"Force field key {ff_key} not found in the force field dictionary.", stacklevel=2)

        settings_dict.update({"dump_modify_flag": "dump_modify" if species else "#", "species": species})
        write_data = {"forcefield.lammps": FF_string}
        if additional_data:
            write_data.update({"extra.data": additional_data})
            settings_dict.update({"extra_data_flag": "include"})

        # Replace all variables
        input_str = Template(self.inputfile.get_str()).safe_substitute(**settings_dict)  # type: ignore[union-attr]
        lines = input_str.split("\n")
        # Filter out the lines where the substitution resulted in a line starting with '#'
        filtered_input_str = "\n".join([line for line in lines if not line.lstrip().startswith("#")])
        input_file = LammpsInputFile.from_str(filtered_input_str, keep_stages=self.keep_stages)

        return LammpsInputSet(
            inputfile=input_file, data=data, calc_type=self.calc_type, additional_data=write_data, **kwargs
        )


@dataclass
class BaseLammpsGenerator(InputGenerator):
    r"""
    Base class to generate LAMMPS input sets.
    Uses template files for the input. The variables that can be changed
    in the input template file are those starting with a $ sign, e.g. $nsteps.
    This generic class is specialized for each template in subclasses, e.g. LammpsMinimization.
    You can create a template for your own task following those present in pymatgen/io/lammps/templates.
    The parameters are then replaced based on the values found
    in the settings dictionary that you provide, e.g. `{"nsteps": 1000}`.

    Attributes:
        template: Path (string) to the template file used to create the InputFile for LAMMPS.
        calc_type: Human-readable string used to briefly describe the type of computations performed by LAMMPS.
        settings: Dictionary containing the values of the parameters to replace in the template.
        keep_stages: If True, the string is formatted in a block structure with stage names
        and newlines that differentiate commands in the respective stages of the InputFile.
        If False, stage names are not printed and all commands appear in a single block.

    /!\ This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    inputfile: LammpsInputFile | None = field(default=None)
    template: str | Path | LammpsInputFile = field(default_factory=str)
    data: LammpsData | CombinedData | None = field(default=None)
    settings: dict = field(default_factory=dict)
    calc_type: str = field(default="lammps")
    keep_stages: bool = field(default=True)

    def get_input_set(self, structure: Structure | LammpsData | CombinedData) -> LammpsInputSet:
        """Generate a LammpsInputSet from the structure/data, tailored to the template file."""

        data = (
            LammpsData.from_structure(structure, atom_style=self.settings.get("atom_style", "full"))
            if isinstance(structure, Structure)
            else structure
        )

        # Load the template
        if Path(self.template).is_file():
            with zopen(self.template, mode="rt", encoding="utf-8") as file:
                template_str = file.read()
        elif isinstance(self.template, LammpsInputFile):
            template_str = self.template.get_str()
        else:
            template_str = self.template

        # Replace all variables
        input_str = Template(template_str).safe_substitute(**self.settings)
        # Get LammpsInputFile
        input_file = LammpsInputFile.from_str(input_str, keep_stages=self.keep_stages)

        # Get the LammpsInputSet from the InputFile and data
        return LammpsInputSet(
            inputfile=input_file,
            data=data,
            calc_type=self.calc_type,
            template_file=self.template,
        )
