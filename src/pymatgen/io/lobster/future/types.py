from __future__ import annotations

from typing import TYPE_CHECKING, Literal, TypedDict

from numpy import complexfloating
from numpy.typing import NDArray

from pymatgen.electronic_structure.core import Spin

if TYPE_CHECKING:
    from typing import TypeAlias

    from numpy import floating, integer


class LobsterInteraction(TypedDict):
    """Dictionary representing a chemical interaction in LOBSTER.

    This dictionary stores information about a specific interaction between atoms in a structure.

    Attributes:
        index (int): The index of the interaction.
        centers (list[str]): List of strings representing the centers of the interaction (e.g., "Fe1", "O2").
        cells (list[list[int]] | NDArray[integer]): List of lists of integers representing the cells of the interaction
        (e.g., [0, 0, 0]).
        orbitals (list[str | None]): List of strings representing the orbitals involved in the interaction
        (e.g., "2s", "2p_x").
        length (float | None): The length of the interaction, representing the distance between the centers.
    """

    index: int
    centers: list[str]
    cells: list[list[int]] | NDArray[integer]
    orbitals: list[str | None]
    length: float | None


class LobsterInteractionData(LobsterInteraction, total=False):
    """Dictionary representing a chemical interaction in LOBSTER with additional COXX/ICOXX data.

    Extends `LobsterInteraction` by adding COXX and ICOXX values for each spin.

    Attributes:
        coxx (dict[Spin, NDArray[floating] | float]): COXX values for each spin.
        icoxx (dict[Spin, NDArray[floating] | float]): ICOXX values for each spin.
    """

    coxx: dict[Spin, NDArray[floating] | float]
    icoxx: dict[Spin, NDArray[floating] | float]


LobsterPopulations: TypeAlias = dict[
    str, dict[str, dict[Spin, dict[Literal["mulliken", "loewdin"], float]]]
]

LobsterMatrixData: TypeAlias = dict[int, dict[Spin | None, NDArray[complexfloating]]]


class LobsterBandOverlaps(TypedDict):
    """Dictionary representing band overlaps in LOBSTER.

    Attributes:
        k_points (dict[Spin, list[list[float]]]): List of k-points for each spin.
        matrices (dict[Spin, list[NDArray[floating]]]): List of matrices for each spin.
        max_deviations (dict[Spin, list[float]]): List of maximal deviations for each spin.
    """

    k_points: dict[Spin, list[list[float]]]
    matrices: dict[Spin, list[NDArray[floating]]]
    max_deviations: dict[Spin, list[float]]


class LobsterFatband(TypedDict):
    """Dictionary representing fatband data in LOBSTER.

    Attributes:
        center (str): Atom associated with the fatband.
        orbital (str): Orbital associated with the fatband.
        energies (dict[Spin, NDArray[floating]]): Energies at each k-point for each spin.
        projections (dict[Spin, NDArray[floating]]): Weights/projections at each k-point for each spin.
    """

    center: str
    orbital: str
    energies: dict[Spin, NDArray[floating]]
    projections: dict[Spin, NDArray[floating]]


class LobsterFatbands(TypedDict):
    """Dictionary representing multiple fatbands in LOBSTER.

    Attributes:
        k_points (list[list[float]]): List of k-points.
        bands (list[LobsterFatband]): List of fatband dictionaries.
    """

    k_points: list[list[float]]
    bands: list[LobsterFatband]
