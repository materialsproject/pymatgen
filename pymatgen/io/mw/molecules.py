"""Molecular submodule of Metalwalls input file"""

from __future__ import annotations

from monty.json import MSONable

from .potentials import Dihedral, HarmonicAngle, HarmonicBond, Improper


class Constraint:
    """Parent class for constraints."""

    @classmethod
    def from_dict(cls, d):
        """Returns corresponding class object based on input dict.

        Arguments:
            d (dict)

        Raises:
            ValueError: Raise error if unknown constraint type is given

        Returns:
            _description_
        """
        if d["type"] == "rigid":
            child_cls = RigidConstraint
        elif d["type"] == "rattle":
            child_cls = RattleConstraint
        else:
            raise ValueError(f"Unknown constraint type '{d['type']}'.")

        return child_cls.from_dict(d)

    def as_dict(self):
        raise NotImplementedError


class RigidConstraint(Constraint, MSONable):
    """RigidConstraint class representing a rigid constraint.

    Attributes:
        site1 (str): First site involved in the constraint.
        site2 (str): Second site involved in the constraint.
        distance (float): Distance for the rigid constraint.

    """

    def __init__(
        self,
        site1: str,
        site2: str,
        distance: float,
        tolerance: float = 1e-6,
        max_iterations: int = 100,
    ):
        self.site1 = site1
        self.site2 = site2
        self.distance = distance
        self.tolerance = tolerance
        self.max_iterations = max_iterations

    def as_dict(self):
        return {
            "type": "rigid",
            "site1": self.site1,
            "site2": self.site2,
            "distance": self.distance,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(site1=d["site1"], site2=d["site2"], distance=d["distance"])


class RattleConstraint(Constraint, MSONable):
    """RattleConstraint class representing a rattle constraint.

    Attributes:
        site1 (str): First site involved in the constraint.
        site2 (str): Second site involved in the constraint.
        distance (float): Distance for the rattle constraint.

    """

    def __init__(self, site1: str, site2: str, distance: float):
        self.site1 = site1
        self.site2 = site2
        self.distance = distance

    def as_dict(self):
        return {
            "type": "rattle",
            "site1": self.site1,
            "site2": self.site2,
            "distance": self.distance,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(site1=d["site1"], site2=d["site2"], distance=d["distance"])


class MoleculeType(MSONable):
    """MoleculeType class representing a molecule type.

    Attributes:
        name (str): Name of the molecule type.
        count (int): Number of molecules of this type in the system.
        sites (List[str]): List of site identifiers in the molecule.
        fourth_site (float): Fourth site parameter (default 0.0).

        constraints_algorithm (ConstraintAlgorithm): Constraints algorithm for the molecule type.
        constraints (List[Constraint]): List of constraints for the molecule type.

        harmonic_bonds (List[HarmonicBond]): List of harmonic bond potentials for the molecule type.
        harmonic_angles (List[HarmonicAngle]): List of harmonic angle potentials for the molecule type.
        dihedrals (List[Dihedral]): List of dihedral potentials for the molecule type.
        impropers (List[Improper]): List of improper potentials for the molecule type.

    """

    def __init__(
        self,
        name: str,
        count: int = 0,
        sites: list[str] | None = None,
        fourth_site: float = 0.0,
    ):
        self.name = name
        self.count = count
        self.sites = sites or []
        self.fourth_site = fourth_site

        # self.constraints_algorithm = None
        self.constraints: list[Constraint] = []

        self.harmonic_bonds: list[HarmonicBond] = []
        self.harmonic_angles: list[HarmonicAngle] = []
        self.dihedrals: list[Dihedral] = []
        self.impropers: list[Improper] = []

    def add_constraint(self, constraint: Constraint):
        """Add a constraint to the molecule type.

        Args:
            constraint (Constraint): Constraint object to add.

        """
        self.constraints.append(constraint)

    def add_harmonic_bond(self, harmonic_bond: HarmonicBond):
        """Add a harmonic bond potential to the molecule type.

        Args:
            harmonic_bond (HarmonicBond): HarmonicBond object to add.

        """
        self.harmonic_bonds.append(harmonic_bond)

    def add_harmonic_angle(self, harmonic_angle: HarmonicAngle):
        """Add a harmonic angle potential to the molecule type.

        Args:
            harmonic_angle (HarmonicAngle): HarmonicAngle object to add.

        """
        self.harmonic_angles.append(harmonic_angle)

    def add_dihedral(self, dihedral: Dihedral):
        """Add a dihedral potential to the molecule type.

        Args:
            dihedral (Dihedral): Dihedral object to add.

        """
        self.dihedrals.append(dihedral)

    def add_improper(self, improper: Improper):
        """Add an improper potential to the molecule type.

        Args:
            improper (Improper): Improper object to add.

        """
        self.impropers.append(improper)

    def as_dict(self):
        return {
            "name": self.name,
            "count": self.count,
            "sites": self.sites,
            "fourth_site": self.fourth_site,
            # "constraints_algorithm": self.constraints_algorithm.as_dict() if self.constraints_algorithm else None,
            "constraints": [constraint.as_dict() for constraint in self.constraints],
            "harmonic_bonds": [bond.as_dict() for bond in self.harmonic_bonds],
            "harmonic_angles": [angle.as_dict() for angle in self.harmonic_angles],
            "dihedrals": [dihedral.as_dict() for dihedral in self.dihedrals],
            "impropers": [improper.as_dict() for improper in self.impropers],
        }

    @classmethod
    def from_dict(cls, d):
        molecule_type = cls(
            name=d["name"],
            count=d["count"],
            sites=d["sites"],
            fourth_site=d["fourth_site"],
        )

        # if d.get("constraints_algorithm"):
        #     molecule_type.constraints_algorithm = ConstraintAlgorithm.from_dict(d["constraints_algorithm"])

        molecule_type.constraints = [Constraint.from_dict(c) for c in d["constraints"]]
        molecule_type.harmonic_bonds = [HarmonicBond.from_dict(b) for b in d["harmonic_bonds"]]
        molecule_type.harmonic_angles = [HarmonicAngle.from_dict(a) for a in d["harmonic_angles"]]
        molecule_type.dihedrals = [Dihedral.from_dict(dh) for dh in d["dihedrals"]]
        molecule_type.impropers = [Improper.from_dict(i) for i in d["impropers"]]

        return molecule_type
