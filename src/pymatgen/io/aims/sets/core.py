"""Module defining core FHI-aims input set generators."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from pymatgen.core import Structure
from pymatgen.io.aims.sets.base import AimsInputGenerator

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core import Molecule


_valid_dynamics: dict[str, tuple[str, ...]] = {
    "nve": ("", "4th_order", "damped"),
    "nvt": ("andersen", "berendsen", "parrinello", "nose-hoover"),
    "gle": ("thermostat",),
}


@dataclass
class StaticSetGenerator(AimsInputGenerator):
    """Common class for ground-state generators.

    Attributes:
        calc_type (str): The type of calculation
    """

    calc_type: str = "static"

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict[str, Any]:
        """Get the parameter updates for the calculation.

        Args:
            structure (Structure or Molecule): The structure to calculate the bands for
            prev_parameters (Dict[str, Any]): The previous parameters

        Returns:
            dict: The updated for the parameters for the output section of FHI-aims
        """
        return prev_parameters


@dataclass
class RelaxSetGenerator(AimsInputGenerator):
    """Generate FHI-aims relax sets for optimizing internal coordinates and lattice params.

    Attributes:
        calc_type (str): The type of calculation
        relax_cell (bool): If True then relax the unit cell from the structure
        max_force (float): Maximum allowed force in the calculation
        method (str): Method used for the geometry optimization
    """

    calc_type: str = "relaxation"
    relax_cell: bool = True
    max_force: float = 1e-3
    method: str = "trm"

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Args:
            structure (Structure or Molecule): The structure to calculate the bands for
        prev_parameters (Dict[str, Any]): The previous parameters

        Returns:
            dict: The updated for the parameters for the output section of FHI-aims
        """
        updates = {"relax_geometry": f"{self.method} {self.max_force:e}"}
        if isinstance(structure, Structure) and self.relax_cell:
            updates["relax_unit_cell"] = "full"
        elif isinstance(structure, Structure):
            updates["relax_unit_cell"] = "none"

        prev_parameters.update(updates)
        return prev_parameters


@dataclass
class SocketIOSetGenerator(AimsInputGenerator):
    """Generate FHI-aims input sets for running with the socket.

    Attributes:
        calc_type (str): The type of calculation
        host (str): The hostname for the server the socket is on
        port (int): The port the socket server is listening on
    """

    calc_type: str = "multi_scf"
    host: str = "localhost"
    port: int = 12345

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Args:
            structure (Structure or Molecule): The structure to calculate the bands for
            prev_parameters (Dict[str, Any]): The previous parameters

        Returns:
            dict: The updated for the parameters for the output section of FHI-aims
        """
        return {"use_pimd_wrapper": (self.host, self.port)}


@dataclass
class MDSetGenerator(AimsInputGenerator):
    """
    A class for generating FHI-aims input sets for molecular dynamics calculations.

    Parameters
    ----------
    ensemble
        Molecular dynamics ensemble to run. Options include `nvt`, `nve`, and `gle`.
        Default: `nve`
    ensemble_specs
        A dictionary containing the specifications of the molecular dynamics ensemble.
        Valid keys are `type` (the ensemble type, valid types are defined in
        `_valid_dynamics` dict), and `parameter` - the control parameter for the thermostat
        (not used for `nve` and `nve_4th_order`).
    temp
        Thermostat temperature. Default: None
    time
        Simulation time (in picoseconds). Negative value stands for indefinite run.
        Default: 5 ps
    time_step
        The time step (in picoseconds) for the simulation. default: 1 fs
    **kwargs
        Other keyword arguments that will be passed to :obj:`AimsInputGenerator`.
    """

    calc_type: str = "md"
    ensemble: str = "nve"
    ensemble_specs: dict[str, Any] = field(default_factory=dict)
    temp: float | None = None
    time: float = 5.0
    time_step: float = 0.001
    init_velocities: bool = True

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Parameters
        ----------
        structure (Structure or Molecule):
            The structure to calculate the bands for
        prev_parameters (Dict[str, Any]):
            The previous parameters

        Returns
        -------
        dict
            A dictionary of updates to apply.
        """
        updates: dict[str, Any] = {"MD_run": [self.time], "MD_time_step": self.time_step}

        # check for ensemble type validity
        default_ensemble_types = {"nve": "", "nvt": "parrinello", "gle": "thermostat"}
        if self.ensemble not in _valid_dynamics:
            raise ValueError(f"Ensemble {self.ensemble} not valid")
        ensemble_type = self.ensemble_specs.get("type", default_ensemble_types[self.ensemble])
        if ensemble_type not in _valid_dynamics[self.ensemble]:
            raise ValueError(
                f"Type {ensemble_type} is not valid for {self.ensemble} ensemble. "
                f"Valid types are: {' ,'.join(_valid_dynamics[self.ensemble])}"
            )
        ensemble_name = f"{self.ensemble.upper()}_{ensemble_type}" if ensemble_type else self.ensemble.upper()
        updates["MD_run"].append(ensemble_name)

        # add temperature
        if self.ensemble == "nve":
            if not self.init_velocities and "velocity" not in structure.site_properties:
                raise ValueError("Velocities must be initialized for NVE ensemble")
        else:
            if self.temp is None:
                raise ValueError(f"Temperature must be set for {ensemble_name} ensemble")
            updates["MD_run"].append(self.temp)

        # check for ensemble control parameter
        ensemble_parameter = self.ensemble_specs.get("parameter", None)
        if ensemble_name not in ("NVE", "NVE_4th_order"):
            if ensemble_parameter is None:
                raise ValueError(f"Ensemble {ensemble_name} parameter is not defined")
            updates["MD_run"].append(ensemble_parameter)

        # ...and put everything in the string
        updates["MD_run"] = " ".join(map(str, updates["MD_run"]))

        # initialize velocities
        if self.init_velocities:
            if self.temp is None:
                raise ValueError("Temperature must be set for velocity initialisation")
            updates["MD_MB_init"] = self.temp

        return updates
