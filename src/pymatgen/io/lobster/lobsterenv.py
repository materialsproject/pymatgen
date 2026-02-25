"""Deprecated class for analyzing NearNeighbors using ICOHPs/ICOOPs/ICOBIs."""

from __future__ import annotations

from monty.dev import deprecated


@deprecated(
    replacement="`pymatgen.analysis.lobster_env.LobsterNeighbors`", category=DeprecationWarning, deadline=(2026, 3, 31)
)
class LobsterNeighbors:
    """Deprecated LobsterNeighbors class for analyzing NearNeighbor interactions using ICOHPs/ICOOPs/ICOBIs."""
