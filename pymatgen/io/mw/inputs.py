"""This module defines the classes and methods for the input of metalwalls
(https://gitlab.com/ampere2/metalwalls)
"""
from __future__ import annotations

from monty.json import MSONable

from .electrodes import ElectrodeCharges, ElectrodeType
from .interactions import Interactions
from .molecules import MoleculeType
from .species import DipolesMinimization, RadiusMinimization, SpeciesType


class GlobalParams(MSONable):
    """Global simulation parameter for MetalWalls

    Args:
        MSONable (num_steps, timestep, temperature, num_pbc): _description_
    """

    def __init__(
        self,
        num_steps: int,
        timestep: float,
        temperature: float,
        num_pbc: int = 3,
    ):
        self.num_steps = num_steps
        self.timestep = timestep
        self.temperature = temperature
        self.num_pbc = num_pbc

    def as_dict(self):
        return {
            "num_steps": self.num_steps,
            "timestep": self.timestep,
            "temperature": self.temperature,
            "num_pbc": self.num_pbc,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            num_steps=d["num_steps"],
            timestep=d["timestep"],
            temperature=d["temperature"],
        )


class Thermostat(MSONable):
    """_summary_

    Args:
        MSONable (_type_): _description_
    """

    def __init__(
        self,
        chain_length: int = 5,
        relaxation_time: float = 4134.1,
        max_iteration: int = 50,
        tolerance: float = 1e-15,
    ):
        self.chain_length = chain_length
        self.relaxation_time = relaxation_time
        self.tolerance = tolerance
        self.max_iteration = max_iteration

    def as_dict(self):
        return {
            "chain_length": self.chain_length,
            "relaxation_time": self.relaxation_time,
            "max_iteration": self.max_iteration,
            "tolerance": self.tolerance,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            chain_length=d["chain_length"],
            relaxation_time=d["relaxation_time"],
            max_iteration=d["max_iteration"],
            tolerance=d["tolerance"],
        )


class Barostat(MSONable):
    """_summary_

    Args:
        MSONable (_type_): _description_
    """

    def __init__(
        self,
        pressure: float = 0,
        chain_length: int = 5,
        relaxation_time: float = 20670.5,
    ):
        self.pressure = pressure
        self.chain_length = chain_length
        self.relaxation_time = relaxation_time

    def as_dict(self):
        return {
            "pressure": self.pressure,
            "chain_length": self.chain_length,
            "relaxation_time": self.relaxation_time,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            pressure=d["pressure"],
            chain_length=d["chain_length"],
            relaxation_time=d["relaxation_time"],
        )


class Velocity(MSONable):
    """Velocity class to define velocity-related actions in the simulation.

    Attributes:
        create (bool): If True, velocities are sampled from a pseudo-random
                       normal distribution at the given temperature.
                       If False, velocities are provided in the data.inpt
                       file and not resampled.
        com_threshold (float): A threshold value triggering a reset of the
                               centre of mass of all mobile particles if the
                               kinetic energy of the centre of mass exceeds
                               the threshold.
        scale (tuple): A tuple specifying the method used to scale velocities
                       and the scaling frequency.
                       The scale method keywords are:
                           - "global": All species are scaled together.
                           - "species": Species types are scaled independently from each other.
                           - "molecules_global": All species are scaled together, species belonging to the same molecule
                                                 are seen as one species type.
                           - "molecules_independent": Species types are scaled independently from each other, species
                                                      belonging to the same molecule are seen as one species type.
                       If scale keyword is not given, no scaling occurs except at the beginning of the simulation if
                       velocities are created. In this case, the default method used is "global" and the scaling
                       frequency is 100.

    """

    def __init__(
        self,
        create: bool = False,
        com_threshold: float = 1.0e6,
        scale: tuple[str, int] = ("global", 100),
    ):
        self.create = create
        self.com_threshold = com_threshold
        self.scale = scale

    def as_dict(self):
        return {
            "create": self.create,
            "com_threshold": self.com_threshold,
            "scale": self.scale,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            create=d["create"],
            com_threshold=d["com_threshold"],
            scale=d["scale"],
        )


class SpeciesBlock(MSONable):
    """SpeciesBlock class representing the species block in the simulation.

    Attributes:
        species (list): List of SpeciesType objects.
        dipoles_minimization (DipolesMinimization): Dipoles minimization configuration.
        radius_minimization (RadiusMinimization): Radius minimization configuration.
    """

    def __init__(
        self,
        species: list[SpeciesType] | None,
        dipoles_minimization: DipolesMinimization,
        radius_minimization: RadiusMinimization,
    ):
        self.species = species or []
        self.dipoles_minimization = dipoles_minimization
        self.radius_minimization = radius_minimization

    def as_dict(self):
        return {
            "species": [species_type.as_dict() for species_type in self.species],
            "dipoles_minimization": self.dipoles_minimization.as_dict(),
            "radius_minimization": self.radius_minimization.as_dict(),
        }

    @classmethod
    def from_dict(cls, d):
        species = [SpeciesType.from_dict(species_type_dict) for species_type_dict in d["species"]]
        dipoles_minimization = DipolesMinimization.from_dict(d["dipoles_minimization"])
        radius_minimization = RadiusMinimization.from_dict(d["radius_minimization"])

        return cls(species, dipoles_minimization, radius_minimization)


class MoleculesBlock(MSONable):
    """MoleculesBlock class representing the molecules block in the simulation.

    Attributes:
        molecule_types (list): List of MoleculeType objects.
    """

    def __init__(self, molecule_types: list[MoleculeType] | None):
        self.molecule_types = molecule_types or []

    def add_molecule_type(self, molecule_type: MoleculeType):
        """Add a MoleculeType to the block.

        Args:
            molecule_type (MoleculeType): MoleculeType object to add.

        """
        self.molecule_types.append(molecule_type)

    def as_dict(self):
        return {
            "molecule_types": [molecule_type.as_dict() for molecule_type in self.molecule_types],
        }

    @classmethod
    def from_dict(cls, d):
        molecule_types = (
            [MoleculeType.from_dict(molecule_type_dict) for molecule_type_dict in d["molecule_types"]]
            if d.get("molecule_types")
            else []
        )

        return cls(molecule_types)


class ElectrodesBlock(MSONable):
    """ElectrodesBlock class representing the electrodes block in the simulation.

    Attributes:
        electrodes (list): List of ElectrodeType objects.
        electrode_charges (ElectrodeCharges): Electrode charges configuration.
        charge_neutrality (bool): Indicates whether the charge neutrality constraint is used.
                                  Default is True.
        global_or_electrodes (str): Specifies the scope of the charge neutrality constraint.
                                    Available options: 'global', 'electrodes'.
                                    Default is 'global'.
    """

    def __init__(
        self,
        electrodes: list[ElectrodeType] | None,
        electrode_charges: ElectrodeCharges,
        charge_neutrality: bool = True,
        global_or_electrodes: str | tuple[str, float] | None = "global",
    ):
        self.electrodes = electrodes or []
        self.electrode_charges = electrode_charges
        self.charge_neutrality = charge_neutrality
        self.global_or_electrodes = global_or_electrodes

    def as_dict(self):
        return {
            "electrodes": [electrode.as_dict() for electrode in self.electrodes],
            "electrode_charges": self.electrode_charges.as_dict(),
            "charge_neutrality": self.charge_neutrality,
            "global_or_electrodes": self.global_or_electrodes,
        }

    @classmethod
    def from_dict(cls, d):
        electrodes = (
            [ElectrodeType.from_dict(electrode_dict) for electrode_dict in d["electrodes"]]
            if d.get("electrodes")
            else []
        )
        electrode_charges = ElectrodeCharges.from_dict(d["electrode_charges"])
        charge_neutrality = d["charge_neutrality"]
        global_or_electrodes = d["global_or_electrodes"]

        return cls(electrodes, electrode_charges, charge_neutrality, global_or_electrodes)


class DipolesElectrodesBlock(MSONable):
    """DipolesElectrodesBlock class representing the dipoles_and_electrodes block in the simulation.

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
        global_or_electrodes: str | tuple[str, float] | None = "global",
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


class Output:
    """Output class representing the output configuration for the simulation.

    Attributes:
        default (int): Default output frequency.
        step (int): Output frequency for general step data.
        restart (int): Output frequency for restart files.
        trajectories (int): Output frequency for atom trajectories.
        xyz (int): Output frequency for atom positions in xyz format.
        pdb (int): Output frequency for atom positions in pdb format.
        lammps (int): Output frequency for atom positions in lammps format.
        charges (int): Output frequency for atom charges.
        total_charges (int): Output frequency for total charges.
        energies (int): Output frequency for energy breakdown.
        forces (int): Output frequency for atom forces.
        elec_forces (int): Output frequency for electrostatic forces.
        dipoles (int): Output frequency for dipoles.
        radius (int): Output frequency for ion radius.
        polarization (int): Output frequency for polarization.
        potential_shift (int): Output frequency for potential shift.
        density (int): Output frequency for ion and charge density profiles.
        temperature (int): Output frequency for temperature.
        box_parameters (int): Output frequency for box parameters.
        pressure (int): Output frequency for pressure.
        stress_tensor (int): Output frequency for stress tensor.
    """

    def __init__(
        self,
        default: int,
        step: int,
        restart: int,
        trajectories: int,
        xyz: int,
        pdb: int,
        lammps: int,
        charges: int,
        total_charges: int,
        energies: int,
        forces: int,
        elec_forces: int,
        dipoles: int,
        radius: int,
        polarization: int,
        potential_shift: int,
        density: int,
        temperature: int,
        box_parameters: int,
        pressure: int,
        stress_tensor: int,
    ):
        self.default = default
        self.step = step
        self.restart = restart
        self.trajectories = trajectories
        self.xyz = xyz
        self.pdb = pdb
        self.lammps = lammps
        self.charges = charges
        self.total_charges = total_charges
        self.energies = energies
        self.forces = forces
        self.elec_forces = elec_forces
        self.dipoles = dipoles
        self.radius = radius
        self.polarization = polarization
        self.potential_shift = potential_shift
        self.density = density
        self.temperature = temperature
        self.box_parameters = box_parameters
        self.pressure = pressure
        self.stress_tensor = stress_tensor

    def as_dict(self):
        return {
            "default": self.default,
            "step": self.step,
            "restart": self.restart,
            "trajectories": self.trajectories,
            "xyz": self.xyz,
            "pdb": self.pdb,
            "lammps": self.lammps,
            "charges": self.charges,
            "total_charges": self.total_charges,
            "energies": self.energies,
            "forces": self.forces,
            "elec_forces": self.elec_forces,
            "dipoles": self.dipoles,
            "radius": self.radius,
            "polarization": self.polarization,
            "potential_shift": self.potential_shift,
            "density": self.density,
            "temperature": self.temperature,
            "box_parameters": self.box_parameters,
            "pressure": self.pressure,
            "stress_tensor": self.stress_tensor,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["default"],
            d["step"],
            d["restart"],
            d["trajectories"],
            d["xyz"],
            d["pdb"],
            d["lammps"],
            d["charges"],
            d["total_charges"],
            d["energies"],
            d["forces"],
            d["elec_forces"],
            d["dipoles"],
            d["radius"],
            d["polarization"],
            d["potential_shift"],
            d["density"],
            d["temperature"],
            d["box_parameters"],
            d["pressure"],
            d["stress_tensor"],
        )


class Plumed(MSONable):
    """Plumed class representing Plumed configuration for bias forces and more.

    Attributes:
        plumed_file (str): Name of the file that contains the Plumed instructions.
    """

    def __init__(self, plumed_file: str):
        self.plumed_file = plumed_file

    def as_dict(self):
        return {"plumed_file": self.plumed_file}

    @classmethod
    def from_dict(cls, d):
        return cls(d["plumed_file"])


class MWInput(MSONable):
    """MWInput class representing the overall structure of the simulation input.

    Attributes:
        global_params (GlobalParams): Global simulation parameters.
        thermostat (Thermostat): Thermostat configuration.
        barostat (Barostat): Barostat configuration.
        velocity (Velocity): Velocity section configuration.
        species_block (SpeciesBlock): Species block configuration.
        molecules_block (MoleculesBlock): Molecules block configuration.
        electrodes_block (ElectrodesBlock): Electrodes block configuration.
        dipoles_electrodes_block (DipolesElectrodesBlock): Dipoles and Electrodes block configuration.
        interactions (Interactions): Interatomic interactions block configuration.
        plumed (Plumed): Plumed configuration for bias forces and more.
        output (Output): Output configuration for the simulation.
    """

    def __init__(
        self,
        global_params: GlobalParams,
        thermostat: Thermostat,
        barostat: Barostat,
        velocity: Velocity,
        species_block: SpeciesBlock,
        molecules_block: MoleculesBlock,
        electrodes_block: ElectrodesBlock,
        dipoles_electrodes_block: DipolesElectrodesBlock,
        interactions: Interactions,
        plumed: Plumed,
        output: Output,
    ):
        self.global_params = global_params
        self.thermostat = thermostat
        self.barostat = barostat
        self.velocity = velocity
        self.species_block = species_block
        self.molecules_block = molecules_block
        self.electrodes_block = electrodes_block
        self.dipoles_electrodes_block = dipoles_electrodes_block
        self.interactions = interactions
        self.plumed = plumed
        self.output = output

    def as_dict(self):
        return {
            "global_params": self.global_params.as_dict(),
            "thermostat": self.thermostat.as_dict(),
            "barostat": self.barostat.as_dict(),
            "velocity": self.velocity.as_dict(),
            "species_block": self.species_block.as_dict(),
            "molecules_block": self.molecules_block.as_dict(),
            "electrodes_block": self.electrodes_block.as_dict(),
            "dipoles_electrodes_block": self.dipoles_electrodes_block.as_dict(),
            "interactions": self.interactions.as_dict(),
            "plumed": self.plumed.as_dict(),
            "output": self.output.as_dict(),
        }

    @classmethod
    def from_dict(cls, d):
        global_params = GlobalParams.from_dict(d["global_params"])
        thermostat = Thermostat.from_dict(d["thermostat"])
        barostat = Barostat.from_dict(d["barostat"])
        velocity = Velocity.from_dict(d["velocity"])
        species_block = SpeciesBlock.from_dict(d["species_block"])
        molecules_block = MoleculesBlock.from_dict(d["molecules_block"])
        electrodes_block = ElectrodesBlock.from_dict(d["electrodes_block"])
        dipoles_electrodes_block = DipolesElectrodesBlock.from_dict(d["dipoles_electrodes_block"])
        interactions = Interactions.from_dict(d["interactions"])
        plumed = Plumed.from_dict(d["plumed"])
        output = Output.from_dict(d["output"])

        return cls(
            global_params,
            thermostat,
            barostat,
            velocity,
            species_block,
            molecules_block,
            electrodes_block,
            dipoles_electrodes_block,
            interactions,
            plumed,
            output,
        )
