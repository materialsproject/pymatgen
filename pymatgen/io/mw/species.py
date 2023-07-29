"""Species submodule of Metalwalls input file"""
from __future__ import annotations

from monty.json import MSONable


class SpeciesType(MSONable):
    """SpeciesType class representing the species_type block.

    Attributes:
        name (str): Unique identifier for this species. Up to 8 ASCII characters, no whitespaces, '#' or '!'.
        count (int): Number of species of this type in the system. Default is 0.
        mass (float): Mass of this species type in the unified atomic mass unit (u). Default is 0.0.
        charge (Union[str, Tuple[str, float], Tuple[str, float, float]]): Charge type and its parameters for this
                                                                          species type.
                     Possible values:
                       - 'neutral': Defines a species with no electrostatic interaction.
                       - ('point', magnitude): Defines a species with a point-charge distribution. 'magnitude' is the
                         magnitude of the charge in the atomic unit (e).
                       - ('gaussian', width, magnitude): Defines a species with a gaussian-charge distribution. 'width'
                         is the η width parameter, and 'magnitude' is the magnitude of the charge in the atomic unit
                         (e).
                     Default is 'neutral'.
        polarizability (float): Polarizability of the species as described in the force field section. Default is 0.0.
        deformability (Tuple[float, float, float]): Deformability parameters of the species as described in the force
                                                    field section.
                                                    Default is (0.0, 0.0, 0.0).
        mobile (bool): Determines if particles of this species type are mobile during the dynamics evolution. Default
                       is True.
        dump_lammps (bool): If True, information for the given species is dumped in the trajectories.lammpstrj file.
                            Frequency is set in the output block. Default is True.
        dump_xyz (bool): If True, particle positions of the given species are dumped in the trajectories.xyz file.
                         Frequency is set in the output block. Default is True.
        dump_trajectories (bool): If True, information for the given species is dumped in the trajectories.out file.
                                  Frequency is set in the output block. Default is True.
        dump_pdb (bool): If True, information for the given species is dumped in the trajectories.pdb file.
                         Frequency is set in the output block. Default is True.
        fourth_site_atom (bool): Flag that states that the atoms of this species are virtual atoms of four-site water
                                 models such as TIP4P or Dang-Chang water models, described in the force field section.
                                 If True, such four-site model is used. Default is False.

    """

    def __init__(
        self,
        name: str,
        count: int = 0,
        mass: float = 0.0,
        charge: str | tuple[str, float] | tuple[str, float, float] = "neutral",
        polarizability: float = 0.0,
        deformability: tuple[float, float, float] = (0.0, 0.0, 0.0),
        mobile: bool = True,
        dump_lammps: bool = True,
        dump_xyz: bool = True,
        dump_trajectories: bool = True,
        dump_pdb: bool = True,
        fourth_site_atom: bool = False,
    ):
        self.name = name
        self.count = count
        self.mass = mass
        self.charge = charge
        self.polarizability = polarizability
        self.deformability = deformability
        self.mobile = mobile
        self.dump_lammps = dump_lammps
        self.dump_xyz = dump_xyz
        self.dump_trajectories = dump_trajectories
        self.dump_pdb = dump_pdb
        self.fourth_site_atom = fourth_site_atom

    def as_dict(self):
        return {
            "name": self.name,
            "count": self.count,
            "mass": self.mass,
            "charge": self.charge,
            "polarizability": self.polarizability,
            "deformability": self.deformability,
            "mobile": self.mobile,
            "dump_lammps": self.dump_lammps,
            "dump_xyz": self.dump_xyz,
            "dump_trajectories": self.dump_trajectories,
            "dump_pdb": self.dump_pdb,
            "fourth_site_atom": self.fourth_site_atom,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            name=d["name"],
            count=d["count"],
            mass=d["mass"],
            charge=d["charge"],
            polarizability=d["polarizability"],
            deformability=d["deformability"],
            mobile=d["mobile"],
            dump_lammps=d["dump_lammps"],
            dump_xyz=d["dump_xyz"],
            dump_trajectories=d["dump_trajectories"],
            dump_pdb=d["dump_pdb"],
            fourth_site_atom=d["fourth_site_atom"],
        )


class DipolesMinimization(MSONable):
    """DipolesMinimization class representing the dipoles_minimization block in the simulation.

    Attributes:
        method (str): Method to compute the ionic dipoles. Available options:
                      - 'cg': Solves μ = aE for each ion using conjugate gradient minimization.
                              Additional parameters: tolerance (real), max_iterations (int).
                      - 'matrix_inversion': Solves μ = aE for each ion using direct matrix inversion.
                      - 'maze_iterative_shake': Solves mass-zero constrained equations using iterative SHAKE algorithm.
                              Additional parameters: tolerance (real), max_iterations (int), nblocks (Optional[int]).
                      - 'maze_inversion': Solves mass-zero constrained equations using direct matrix inversion.

        tolerance (float): Tolerance criterion for the convergence. Default is 1.0e-12.

        max_iterations (int): Maximum number of iterations allowed. Default is 100.

        nblocks (Optional[int]): Level of approximation for the SHAKE matrix. Optional parameter for the
                                 'maze_iterative_shake' method. Default is 0.

    """

    def __init__(
        self,
        method: str,
        tolerance: float = 1.0e-12,
        max_iterations: int = 100,
        nblocks: int | None = 0,
        # TODO: add preconditioner jacobi
    ):
        known_methods = [
            "cg",
            "matrix_inversion",
            "maze_iterative_shake",
            "maze_inversion",
        ]
        if method not in known_methods:
            raise ValueError(f"Unknown method '{method}'. Supported methods are: {', '.join(known_methods)}")
        self.method = method
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.nblocks = nblocks

    def as_dict(self):
        return {
            "method": self.method,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "nblocks": self.nblocks,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            method=d["method"],
            tolerance=d["tolerance"],
            max_iterations=d["max_iterations"],
            nblocks=d.get("nblocks"),
        )


class RadiusMinimization(MSONable):
    """RadiusMinimization class representing the radius_minimization block in the simulation.

    Attributes:
        method (str): Method to compute the radii. Available option:
                      - 'cg': Uses the nonlinear conjugate gradient algorithm for radii computation.
                              Additional parameters: tolerance (real), max_iterations (int).

        tolerance (float): Tolerance criterion for the convergence. Default is 1.0e-12.

        max_iterations (int): Maximum number of iterations allowed. Default is 100.

    """

    def __init__(
        self,
        method: str,
        tolerance: float = 1.0e-12,
        max_iterations: int = 100,
    ):
        known_methods = ["cg"]
        if method not in known_methods:
            raise ValueError(f"Unknown method '{method}'. Supported methods are: {', '.join(known_methods)}")

        self.method = method
        self.tolerance = tolerance
        self.max_iterations = max_iterations

    def as_dict(self):
        return {
            "method": self.method,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            method=d["method"],
            tolerance=d["tolerance"],
            max_iterations=d["max_iterations"],
        )
