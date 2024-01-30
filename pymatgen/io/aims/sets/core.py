"""Module defining core FHI-aims input set generators."""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from pymatgen.core import Structure
from pymatgen.io.aims.sets.base import AimsInputGenerator

if TYPE_CHECKING:
    from pymatgen.core import Molecule


@dataclass
class StaticSetGenerator(AimsInputGenerator):
    """Common class for ground-state generators.

    Parameters
    ----------
    calc_type: str
        The type of calculation
    """

    calc_type: str = "static"

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict[str, Any]:
        """Get the parameter updates for the calculation.

        Parameters
        ----------
        structure: Structure or Molecule
            The structure to calculate the bands for
        prev_parameters: Dict[str, Any]
            The previous parameters

        Returns
        -------
        The updated for the parameters for the output section of FHI-aims
        """
        return prev_parameters


@dataclass
class RelaxSetGenerator(AimsInputGenerator):
    """
    Class to generate FHI-aims relax sets.

    I.e., sets for optimization of internal and lattice coordinates.

    Parameters
    ----------
    calc_type: str
        The type of calculation
    relax_cell: bool
        If True then relax the unit cell from the structure
    max_force: float
        Maximum allowed force in the calculation
    method: str
        Method used for the geometry optimization
    """

    calc_type: str = "relaxation"
    relax_cell: bool = True
    max_force: float = 1e-3
    method: str = "trm"

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Parameters
        ----------
        structure: Structure or Molecule
            The structure to calculate the bands for
        prev_parameters: Dict[str, Any]
            The previous parameters

        Returns
        -------
        The updated for the parameters for the output section of FHI-aims
        """
        updates = {"relax_geometry": f"{self.method} {self.max_force:e}"}
        if isinstance(structure, Structure) and self.relax_cell:
            updates["relax_unit_cell"] = "full"
        elif isinstance(structure, Structure):
            updates["relax_unit_cell"] = "none"

        return updates


@dataclass
class SocketIOSetGenerator(AimsInputGenerator):
    """Class to generate FHI-aims input sets for running with the socket.

    Parameters
    ----------
    calc_type: str
        The type of calculation
    host: str
        The hostname for the server the socket is on
    port: int
        The port the socket server is listening on
    """

    calc_type: str = "multi_scf"
    host: str = "localhost"
    port: int = 12345

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Parameters
        ----------
        structure: Structure or Molecule
            The structure to calculate the bands for
        prev_parameters: Dict[str, Any]
            The previous parameters

        Returns
        -------
        The updated for the parameters for the output section of FHI-aims
        """
        return {"use_pimd_wrapper": (self.host, self.port)}
