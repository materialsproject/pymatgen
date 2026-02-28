"""Constants for LOBSTER input and output parsing."""

from __future__ import annotations

"""The version of LOBSTER that this code is compatible with. To be updated when necessary."""
LOBSTER_VERSION: str = "5.1.1"

"""A tuple of possible orbitals that LOBSTER can handle and their string representations in LOBSTER output files."""
LOBSTER_ORBITALS: tuple[str, ...] = (
    "s",
    "p_x",
    "p_y",
    "p_z",
    "d_xy",
    "d_xz",
    "d_yz",
    "d_z^2",
    "d_x^2-y^2",
    "f_xyz",
    "f_xz^2",
    "f_yz^2",
    "f_z^3",
    "f_x(x^2-3y^2)",
    "f_y(3x^2-y^2)",
    "f_z(x^2-y^2)",
)
