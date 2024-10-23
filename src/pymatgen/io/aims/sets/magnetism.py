"""Define the InputSetGenerators for FHI-aims magnetism calculations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pymatgen.io.aims.sets.core import RelaxSetGenerator, StaticSetGenerator

if TYPE_CHECKING:
    from typing import Any

    from pymatgen.core.structure import Molecule, Structure


@dataclass
class MagneticStaticSetGenerator(StaticSetGenerator):
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
        updates = {
            "spin": "collinear",
            "output": [*prev_parameters.get("output", []), "mulliken"],
        }
        prev_parameters.update(updates)
        return prev_parameters


@dataclass
class MagneticRelaxSetGenerator(RelaxSetGenerator):
    """Generate FHI-aims relax sets for optimizing internal coordinates and lattice params.

    Attributes:
        calc_type (str): The type of calculation
        relax_cell (bool): If True then relax the unit cell from the structure
        max_force (float): Maximum allowed force in the calculation
        method (str): Method used for the geometry optimization
    """

    def get_parameter_updates(self, structure: Structure | Molecule, prev_parameters: dict[str, Any]) -> dict:
        """Get the parameter updates for the calculation.

        Args:
            structure (Structure or Molecule): The structure to calculate the bands for
        prev_parameters (Dict[str, Any]): The previous parameters

        Returns:
            dict: The updated for the parameters for the output section of FHI-aims
        """
        prev_parameters = super().get_parameter_updates(structure=structure, prev_parameters=prev_parameters)
        updates = {
            "spin": "collinear",
            "output": [*prev_parameters.get("output", []), "mulliken"],
        }
        prev_parameters.update(updates)
        return prev_parameters
