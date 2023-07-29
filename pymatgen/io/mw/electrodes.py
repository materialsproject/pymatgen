"""Electrodes submodule of Metalwalls input file"""

from __future__ import annotations

from monty.json import MSONable


class ElectrodeType(MSONable):
    """Class to represent an electrode type in a simulation configuration.

    :param name: Unique identifier for this electrode species. Must be up to 8 ASCII characters and
                 cannot contain white spaces, '#' or '!'.
    :type name: str
    :param species: The name of the species used to fill this electrode.
    :type species: str
    :param potential: Constant potential for this electrode type. Default is 0.0.
    :type potential: float, optional
    :param piston: Parameters for NPT-like simulation. The first element is the target pressure,
                   and the second element is the direction of the applied force (should be ±1).
                   This option automatically sets compute_force to true. The piston can only be
                   applied with the conjugate gradient method. Default is (0.0, 1).
    :type piston: tuple(float, int), optional
    :param thomas_fermi_length: Thomas-Fermi length for this electrode. Default is 0.0.
    :type thomas_fermi_length: float, optional
    :param voronoi_volume: Voronoi volume for this electrode.
    :type voronoi_volume: float, optional
    """

    def __init__(
        self,
        name: str,
        species: str,
        potential: float = 0.0,
        piston: tuple[float, int] | None = (0.0, 1),
        thomas_fermi_length: float | None = 0.0,
        voronoi_volume: float | None = None,
    ):
        self.name = name
        self.species = species
        self.potential = potential
        self.piston = piston
        self.thomas_fermi_length = thomas_fermi_length
        self.voronoi_volume = voronoi_volume

    @classmethod
    def from_dict(cls, d: dict) -> ElectrodeType:
        """
        Create an ElectrodeType object from a dictionary.

        :param d: The dictionary containing the electrode type information.
        :type d: dict
        :return: An ElectrodeType object created from the dictionary.
        :rtype: ElectrodeType
        """
        return cls(
            name=d["name"],
            species=d["species"],
            potential=d.get("potential", 0.0),
            piston=d.get("piston", (0.0, 1)),
            thomas_fermi_length=d.get("thomas_fermi_length", 0.0),
            voronoi_volume=d.get("voronoi_volume"),
        )

    def as_dict(self) -> dict:
        """
        Convert the ElectrodeType object to a dictionary.

        :return: A dictionary representing the ElectrodeType object.
        :rtype: dict
        """
        return {
            "name": self.name,
            "species": self.species,
            "potential": self.potential,
            "piston": self.piston,
            "thomas_fermi_length": self.thomas_fermi_length,
            "voronoi_volume": self.voronoi_volume,
        }


class ElectrodeCharges(MSONable):
    """
    Class representing electrode charges configuration for a simulation.

    :param method: The method used to compute the charges on electrodes at each time step.
                   Supported methods: 'constant_charge', 'matrix_inversion', 'cg', 'maze_inversion',
                   'maze_iterative_shake'.
    :type method: str
    :param tolerance: The tolerance criterion for convergence. Default is 1.0e-12.
    :type tolerance: float, optional
    :param max_iterations: The maximum number of iterations allowed. Default is 100.
    :type max_iterations: int, optional
    :param preconditioner: The method used for preconditioning the conjugate gradient algorithm.
                           Only specified if method is 'cg'. Default is None.
    :type preconditioner: str, optional
    :param nblocks: The level of approximation for the SHAKE matrix. Only specified if method is
                    'maze_iterative_shake'.
                    Default is 0.
    :type nblocks: int, optional
    """

    def __init__(
        self,
        method: str,
        tolerance: float | None = 1.0e-12,
        max_iterations: int | None = 100,
        preconditioner: str | None | None = None,
        nblocks: int | None = 0,
    ):
        known_methods = [
            "constant_charge",
            "matrix_inversion",
            "cg",
            "maze_inversion",
            "maze_iterative_shake",
        ]
        if method not in known_methods:
            raise ValueError(f"Unknown method '{method}'. Supported methods are: {', '.join(known_methods)}")
        self.method = method
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.preconditioner = preconditioner
        self.nblocks = nblocks

    def as_dict(self):
        """
        Convert the ElectrodeCharges object to a dictionary.

        :return: A dictionary representing the ElectrodeCharges object.
        :rtype: dict
        """
        return {
            "method": self.method,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "preconditioner": self.preconditioner,
            "nblocks": self.nblocks,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Create an ElectrodeCharges object from a dictionary.

        :param d: The dictionary containing the ElectrodeCharges information.
        :type d: dict
        :return: An ElectrodeCharges object created from the dictionary.
        :rtype: ElectrodeCharges
        """
        return cls(
            method=d["method"],
            tolerance=d.get("tolerance", 1.0e-12),
            max_iterations=d.get("max_iterations", 100),
            preconditioner=d.get("preconditioner"),
            nblocks=d.get("nblocks", 0),
        )


class DipolesAndElectrodes(MSONable):
    """DipolesAndElectrodes class representing the dipoles_and_electrodes block in the simulation.

    Attributes:
        algorithm (str): Algorithm used to compute dipoles and electrode charges.
                         Available options:
                         - 'cg': Consistently computes dipoles and electrode charges using a conjugate-gradient method.
                                 Additional parameters: tolerance (real), max_iterations (int).
                         - 'cg_and_constant_charge': Keeps electrode charges constant while computing dipoles using cg.
                                 Additional parameters: tolerance (real), max_iterations (int).
        tolerance (float): Tolerance criterion for the convergence. Default is 1.0e-12.
        max_iterations (int): Maximum number of iterations allowed. Default is 100.
        preconditioner (str): Preconditioner method used for the conjugate-gradient algorithm.
                              Available options: 'jacobi'.

        charge_neutrality (bool): Indicates whether the charge neutrality constraint is used.
                                  Default is True.
        global_or_electrodes (str): Specifies the scope of the charge neutrality constraint.
                                    Available options: 'global', 'electrodes'.
                                    Default is 'global'.
    """

    def __init__(
        self,
        algorithm: str,
        tolerance: float = 1.0e-12,
        max_iterations: int = 100,
        preconditioner: str | None = None,
        charge_neutrality: bool = True,
        global_or_electrodes: str | None = "global",
    ):
        known_algorithms = ["cg", "cg_and_constant_charge"]
        if algorithm not in known_algorithms:
            raise ValueError(
                f"Unknown algorithm '{algorithm}'. Supported algorithms are: {', '.join(known_algorithms)}"
            )

        self.algorithm = algorithm
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.preconditioner = preconditioner
        self.charge_neutrality = charge_neutrality
        self.global_or_electrodes = global_or_electrodes

    def as_dict(self):
        return {
            "algorithm": self.algorithm,
            "tolerance": self.tolerance,
            "max_iterations": self.max_iterations,
            "preconditioner": self.preconditioner,
            "charge_neutrality": self.charge_neutrality,
            "global_or_electrodes": self.global_or_electrodes,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            algorithm=d["algorithm"],
            tolerance=d["tolerance"],
            max_iterations=d["max_iterations"],
            preconditioner=d.get("preconditioner"),
            charge_neutrality=d.get("charge_neutrality", True),
            global_or_electrodes=d.get("global_or_electrodes", "global"),
        )
