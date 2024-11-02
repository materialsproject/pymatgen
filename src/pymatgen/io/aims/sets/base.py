"""Module defining base FHI-aims input set and generator."""

from __future__ import annotations

import copy
import json
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any
from warnings import warn

import numpy as np
from monty.json import MontyDecoder, MontyEncoder

from pymatgen.core import Molecule, Structure
from pymatgen.io.aims.inputs import AimsControlIn, AimsGeometryIn
from pymatgen.io.aims.parsers import AimsParseError, read_aims_output
from pymatgen.io.core import InputFile, InputGenerator, InputSet

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from pymatgen.util.typing import PathLike

TMPDIR_NAME: str = "tmpdir"
OUTPUT_FILE_NAME: str = "aims.out"
CONTROL_FILE_NAME: str = "control.in"
PARAMS_JSON_FILE_NAME: str = "parameters.json"
GEOMETRY_FILE_NAME: str = "geometry.in"


DEFAULT_AIMS_PROPERTIES = (
    "energy",
    "free_energy",
    "forces",
    "stress",
    "stresses",
    "dipole",
    "magmom",
)


class AimsInputSet(InputSet):
    """A class to represent a set of Aims inputs."""

    def __init__(
        self,
        parameters: dict[str, Any],
        structure: Structure | Molecule,
        properties: Sequence[str] = ("energy", "free_energy"),
    ) -> None:
        """Construct the AimsInputSet.

        Args:
            parameters (dict[str, Any]): The ASE parameters object for the calculation
            structure (Structure or Molecule): The Structure/Molecule objects to
                create the inputs for
            properties (Sequence[str]): The properties to calculate for the calculation
        """
        self._parameters = parameters
        self._structure = structure
        self._properties = properties

        aims_control_in, aims_geometry_in = self.get_input_files()
        super().__init__(
            inputs={
                CONTROL_FILE_NAME: aims_control_in,
                GEOMETRY_FILE_NAME: aims_geometry_in,
                PARAMS_JSON_FILE_NAME: json.dumps(self._parameters, cls=MontyEncoder),
            }
        )

    def get_input_files(self) -> tuple[str, str]:
        """Get the input file contents for the calculation.

        Returns:
            tuple[str, str]: The contents of the control.in and geometry.in file
        """
        property_flags = {
            "forces": "compute_forces",
            "stress": "compute_analytical_stress",
            "stresses": "compute_heat_flux",
        }
        updated_params = dict(**self._parameters)
        for prop in self._properties:
            aims_name = property_flags.get(prop)
            if aims_name is not None:
                updated_params[aims_name] = True

        aims_geometry_in = AimsGeometryIn.from_structure(self._structure)
        aims_control_in = AimsControlIn(updated_params)

        return (
            aims_control_in.get_content(structure=self._structure),
            aims_geometry_in.content,
        )

    @property
    def control_in(self) -> str | slice | InputFile:
        """The control.in file contents."""
        return self[CONTROL_FILE_NAME]

    @property
    def geometry_in(self) -> str | slice | InputFile:
        """The geometry.in file contents."""
        return self[GEOMETRY_FILE_NAME]

    @property
    def params_json(self) -> str | slice | InputFile:
        """The JSON representation of the parameters dict."""
        return self[PARAMS_JSON_FILE_NAME]

    def set_parameters(self, *args, **kwargs) -> dict[str, Any]:
        """Set the parameters object for the AimsTemplate.

        This sets the parameters object that is passed to an AimsTemplate and
        resets the control.in file

        One can pass a dictionary mapping the aims variables to their values or
        the aims variables as keyword arguments. A combination of the two
        options is also allowed.

        Returns:
            dict[str, Any]: dictionary with the variables that have been added.
        """
        self._parameters.clear()
        for arg in args:
            self._parameters.update(arg)

        self._parameters.update(kwargs)

        aims_control_in, _ = self.get_input_files()

        self.inputs[CONTROL_FILE_NAME] = aims_control_in
        self.inputs[PARAMS_JSON_FILE_NAME] = json.dumps(self._parameters, cls=MontyEncoder)

        inputs = {str(key): val for key, val in self.inputs.items()}
        self.__dict__.update(inputs)

        return self._parameters

    def remove_parameters(
        self,
        keys: Iterable[str] | str,
        strict: bool = True,
    ) -> dict[str, Any]:
        """Remove the aims parameters listed in keys.

        This removes the aims variables from the parameters object.

        Args:
            keys (Iterable[str] or str): string or list of strings with the names of
                the aims parameters to be removed.
            strict (bool): whether to raise a KeyError if one of the aims parameters
                to be removed is not present.

        Returns:
            dict[str, Any]: Dictionary with the variables that have been removed.
        """
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            if key not in self._parameters:
                if strict:
                    raise ValueError(f"{key=} not in {list(self._parameters)=}")
                continue

            del self._parameters[key]

        return self.set_parameters(**self._parameters)

    def set_structure(self, structure: Structure | Molecule):
        """Set the structure object for this input set.

        Args:
            structure (Structure or Molecule): The new Structure or Molecule
                for the calculation
        """
        self._structure = structure

        aims_control_in, aims_geometry_in = self.get_input_files()
        self.inputs[GEOMETRY_FILE_NAME] = aims_geometry_in
        self.inputs[CONTROL_FILE_NAME] = aims_control_in
        inputs = {str(key): val for key, val in self.inputs.items()}
        self.__dict__.update(inputs)


@dataclass
class AimsInputGenerator(InputGenerator):
    """
    A class to generate Aims input sets.

    Attributes:
        user_params (dict[str, Any]): Updates the default
            parameters for the FHI-aims calculator
        user_kpoints_settings (dict[str, Any]):  The settings
            used to create the k-grid parameters for FHI-aims
        use_structure_charge (bool): If set to True, then the overall charge of the
            structure (structure.charge) is used to set the `charge` variable in the
            `control.in`. Default is False.
    """

    user_params: dict[str, Any] = field(default_factory=dict)
    user_kpoints_settings: dict[str, Any] = field(default_factory=dict)
    use_structure_charge: bool = False

    def get_input_set(
        self,
        structure: Structure | Molecule | None = None,
        prev_dir: PathLike | None = None,
        properties: list[str] | None = None,
    ) -> AimsInputSet:
        """Generate an AimsInputSet object.

        Args:
            structure (Structure or Molecule): Structure or
                Molecule to generate the input set for.
            prev_dir (str or Path): Path to the previous working directory
            properties (list[str]): System properties that are being calculated

        Returns:
            AimsInputSet: The input set for the calculation of structure
        """
        prev_structure, prev_parameters, _ = self._read_previous(prev_dir)

        structure = structure or prev_structure

        if structure is None:
            raise ValueError("No structure can be determined to generate the input set")

        parameters = self._get_input_parameters(structure, prev_parameters)
        properties = self._get_properties(properties, parameters)

        return AimsInputSet(parameters=parameters, structure=structure, properties=properties)

    @staticmethod
    def _read_previous(
        prev_dir: PathLike | None = None,
    ) -> tuple[Structure | Molecule | None, dict[str, Any], dict[str, Any]]:
        """Read in previous results.

        Args:
            prev_dir (str or Path): The previous directory for the calculation
        """
        prev_structure: Structure | Molecule | None = None
        prev_params = {}
        prev_results: dict[str, Any] = {}

        if prev_dir:
            # strip hostname from the directory (not good, works only with run_locally.
            # Should be checked with Fireworks, will not for sure work with
            # jobflow_remote)
            split_prev_dir = str(prev_dir).split(":")[-1]
            with open(f"{split_prev_dir}/parameters.json") as param_file:
                prev_params = json.load(param_file, cls=MontyDecoder)

            try:
                aims_output: Sequence[Structure | Molecule] = read_aims_output(
                    f"{split_prev_dir}/aims.out", index=slice(-1, None)
                )
                prev_structure = aims_output[0]

                prev_results = prev_structure.properties
                prev_results.update(prev_structure.site_properties)
            except (IndexError, AimsParseError):
                pass

        return prev_structure, prev_params, prev_results

    @staticmethod
    def _get_properties(
        properties: list[str] | None = None,
        parameters: dict[str, Any] | None = None,
    ) -> list[str]:
        """Get the properties to calculate.

        Args:
            properties (list[str]): The currently requested properties
            parameters (dict[str, Any]): The parameters for this calculation

        Returns:
            list[str]: The list of properties to calculate
        """
        if properties is None:
            properties = ["energy", "free_energy"]

        if parameters is None:
            return properties

        if "compute_forces" in parameters and "forces" not in properties:
            properties.append("forces")
        if "compute_heat_flux" in parameters and "stresses" not in properties:
            properties.append("stress")
            properties.append("stresses")
        if "stress" not in properties and (
            ("compute_analytical_stress" in parameters)
            or ("compute_numerical_stress" in parameters)
            or ("compute_heat_flux" in parameters)
        ):
            properties.append("stress")

        return properties

    def _get_input_parameters(
        self,
        structure: Structure | Molecule,
        prev_parameters: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Create the input parameters.

        Args:
            structure (Structure | Molecule): The structure
                or molecule for the system
            prev_parameters (dict[str, Any]): The previous
                calculation's calculation parameters

        Returns:
            dict: The input object
        """
        # Get the default config
        # FHI-aims recommends using their defaults so bare-bones default params
        params: dict[str, Any] = {"xc": "pbe", "relativistic": "atomic_zora scalar"}

        # Override default parameters with previous parameters
        prev_parameters = {} if prev_parameters is None else copy.deepcopy(prev_parameters)
        prev_parameters.pop("relax_geometry", None)
        prev_parameters.pop("relax_unit_cell", None)

        kpt_settings = copy.deepcopy(self.user_kpoints_settings)
        if isinstance(structure, Structure) and "k_grid" in prev_parameters:
            density = self.k2d(structure, prev_parameters.pop("k_grid"))
            if "density" not in kpt_settings:
                kpt_settings["density"] = density

        parameter_updates = self.get_parameter_updates(structure, prev_parameters)
        params = recursive_update(params, parameter_updates)
        # Add the structure charge (useful for defect workflows)
        if self.use_structure_charge:
            params["charge"] = structure.charge

        # Override default parameters with user_params
        params = recursive_update(params, self.user_params)
        if ("k_grid" in params) and ("density" in kpt_settings):
            warn(
                "WARNING: the k_grid is set in user_params and in the kpt_settings,"
                " using the one passed in user_params.",
                stacklevel=1,
            )
        elif isinstance(structure, Structure) and ("k_grid" not in params):
            density = kpt_settings.get("density", 5.0)
            even = kpt_settings.get("even", True)
            params["k_grid"] = self.d2k(structure, density, even)
        elif isinstance(structure, Molecule) and "k_grid" in params:
            warn("WARNING: removing unnecessary k_grid information", stacklevel=1)
            del params["k_grid"]

        return params

    def get_parameter_updates(
        self,
        structure: Structure | Molecule,
        prev_parameters: dict[str, Any],
    ) -> dict[str, Any]:
        """Update the parameters for a given calculation type.

        Args:
            structure (Structure or Molecule): The system to run
            prev_parameters (dict[str, Any]): Previous calculation parameters.

        Returns:
            dict: A dictionary of updates to apply.
        """
        return prev_parameters

    def d2k(
        self,
        structure: Structure,
        kpt_density: float | tuple[float, float, float] = 5.0,
        even: bool = True,
    ) -> Iterable[float]:
        """Convert k-point density to Monkhorst-Pack grid size.

        inspired by [ase.calculators.calculator.kptdensity2monkhorstpack]

        Args:
            structure (Structure): Contains unit cell and
                information about boundary conditions.
            kpt_density (float | list[float]): Required k-point
                density.  Default value is 5.0 point per Ang^-1.
            even (bool): Round up to even numbers.

        Returns:
            dict: Monkhorst-Pack grid size in all directions
        """
        recip_cell = structure.lattice.inv_matrix.transpose()
        return self.d2k_recip_cell(recip_cell, structure.lattice.pbc, kpt_density, even)

    def k2d(self, structure: Structure, k_grid: np.ndarray[int]):
        """Generate the kpoint density in each direction from given k_grid.

        Args:
            structure: Structure
                Contains unit cell and information about boundary conditions.
            k_grid: np.ndarray[int]
                k_grid that was used.

        Returns:
            dict: Density of kpoints in each direction. result.mean() computes average density
        """
        recip_cell = structure.lattice.inv_matrix.transpose()
        densities = k_grid / (2 * np.pi * np.sqrt((recip_cell**2).sum(axis=1)))
        return np.array(densities)

    @staticmethod
    def d2k_recip_cell(
        recip_cell: np.ndarray,
        pbc: Sequence[bool],
        kpt_density: float | tuple[float, float, float] = 5.0,
        even: bool = True,
    ) -> Sequence[int]:
        """Convert k-point density to Monkhorst-Pack grid size.

        Args:
            recip_cell (Cell): The reciprocal cell
            pbc (Sequence[bool]): If element of pbc is True
                then system is periodic in that direction
            kpt_density (float or list[floats]): Required k-point
                density. Default value is 5 points per Ang^-1.
        even(bool): Round up to even numbers.

        Returns:
            dict: Monkhorst-Pack grid size in all directions
        """
        if isinstance(kpt_density, float):
            kpt_density = (kpt_density, kpt_density, kpt_density)
        kpts: list[int] = []
        for i in range(3):
            if pbc[i]:
                k = 2 * np.pi * np.sqrt((recip_cell[i] ** 2).sum()) * float(kpt_density[i])
                if even:
                    kpts.append(2 * int(np.ceil(k / 2)))
                else:
                    kpts.append(int(np.ceil(k)))
            else:
                kpts.append(1)
        return kpts


def recursive_update(dct: dict, up: dict) -> dict:
    """
    Update a dictionary recursively and return it.

    Args:
        dct (dict): Input dictionary to modify
        up (dict): updates to apply

    Returns:
        dict: The updated dictionary.

    Example:
        d = {'activate_hybrid': {"hybrid_functional": "HSE06"}}
        u = {'activate_hybrid': {"cutoff_radius": 8}}

        yields {'activate_hybrid': {"hybrid_functional": "HSE06", "cutoff_radius": 8}}}
    """
    for key, val in up.items():
        if isinstance(val, dict):
            dct[key] = recursive_update(dct.get(key, {}), val)
        elif key == "output" and isinstance(val, list):  # for all other keys the list addition is not needed (I guess)
            old_v = dct.get(key, [])
            dct[key] = old_v + val
        else:
            dct[key] = val
    return dct
