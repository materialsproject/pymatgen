"""Potential submodule defined for MW input"""
from __future__ import annotations

from monty.json import MSONable


class HarmonicBond(MSONable):
    """HarmonicBond class representing a harmonic bond potential.

    Attributes:
        site1 (str): Identifier of the first site involved in the harmonic bond.
        site2 (str): Identifier of the second site involved in the harmonic bond.
        k0 (float): Strength of the harmonic bond in atomic units (Eh/a0).
        r0 (float): Length of the harmonic bond.

    """

    def __init__(self, site1: str, site2: str, k0: float, r0: float):
        self.site1 = site1
        self.site2 = site2
        self.k0 = k0
        self.r0 = r0

    def as_dict(self):
        return {
            "type": "harmonic_bond",
            "site1": self.site1,
            "site2": self.site2,
            "k0": self.k0,
            "r0": self.r0,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(site1=d["site1"], site2=d["site2"], k0=d["k0"], r0=d["r0"])


class HarmonicAngle(MSONable):
    """HarmonicAngle class representing a harmonic angle potential.

    Attributes:
        site1 (str): Identifier of the first site involved in the harmonic angle.
        site2 (str): Identifier of the second site involved in the harmonic angle.
        site3 (str): Identifier of the third site involved in the harmonic angle.
        k0 (float): Strength of the harmonic angle in atomic units (Eh/rad).
        theta0 (float): Angle at equilibrium in radians.

    """

    def __init__(self, site1: str, site2: str, site3: str, k0: float, theta0: float):
        self.site1 = site1
        self.site2 = site2
        self.site3 = site3
        self.k0 = k0
        self.theta0 = theta0

    def as_dict(self):
        return {
            "type": "harmonic_angle",
            "site1": self.site1,
            "site2": self.site2,
            "site3": self.site3,
            "k0": self.k0,
            "theta0": self.theta0,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            site1=d["site1"],
            site2=d["site2"],
            site3=d["site3"],
            k0=d["k0"],
            theta0=d["theta0"],
        )


class Dihedral(MSONable):
    """Dihedral class representing a dihedral potential.

    Attributes:
        site1 (str): Identifier of the first site involved in the dihedral.
        site2 (str): Identifier of the second site involved in the dihedral.
        site3 (str): Identifier of the third site involved in the dihedral.
        site4 (str): Identifier of the fourth site involved in the dihedral.
        v1 (float): Parameter v1 of the dihedral potential in atomic units (Eh).
        v2 (float): Parameter v2 of the dihedral potential in atomic units (Eh).
        v3 (float): Parameter v3 of the dihedral potential in atomic units (Eh).
        v4 (float): Parameter v4 of the dihedral potential in atomic units (Eh).

    """

    def __init__(
        self,
        site1: str,
        site2: str,
        site3: str,
        site4: str,
        v1: float,
        v2: float,
        v3: float,
        v4: float,
    ):
        self.site1 = site1
        self.site2 = site2
        self.site3 = site3
        self.site4 = site4
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4

    def as_dict(self):
        return {
            "type": "dihedral",
            "site1": self.site1,
            "site2": self.site2,
            "site3": self.site3,
            "site4": self.site4,
            "v1": self.v1,
            "v2": self.v2,
            "v3": self.v3,
            "v4": self.v4,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            site1=d["site1"],
            site2=d["site2"],
            site3=d["site3"],
            site4=d["site4"],
            v1=d["v1"],
            v2=d["v2"],
            v3=d["v3"],
            v4=d["v4"],
        )


class Improper(MSONable):
    """Improper class representing an improper potential.

    Attributes:
        site1 (str): Identifier of the first site involved in the improper.
        site2 (str): Identifier of the second site involved in the improper.
        site3 (str): Identifier of the third site involved in the improper.
        site4 (str): Identifier of the fourth site involved in the improper.
        v1 (float): Parameter v1 of the improper potential in atomic units (Eh).
        v2 (float): Parameter v2 of the improper potential in atomic units (Eh).
        v3 (float): Parameter v3 of the improper potential in atomic units (Eh).
        v4 (float): Parameter v4 of the improper potential in atomic units (Eh).

    """

    def __init__(
        self,
        site1: str,
        site2: str,
        site3: str,
        site4: str,
        v1: float,
        v2: float,
        v3: float,
        v4: float,
    ):
        self.site1 = site1
        self.site2 = site2
        self.site3 = site3
        self.site4 = site4
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4

    def as_dict(self):
        return {
            "type": "improper",
            "site1": self.site1,
            "site2": self.site2,
            "site3": self.site3,
            "site4": self.site4,
            "v1": self.v1,
            "v2": self.v2,
            "v3": self.v3,
            "v4": self.v4,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            site1=d["site1"],
            site2=d["site2"],
            site3=d["site3"],
            site4=d["site4"],
            v1=d["v1"],
            v2=d["v2"],
            v3=d["v3"],
            v4=d["v4"],
        )
