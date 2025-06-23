"""
This module provides conversion between the Atomic Simulation Environment (ASE)
Atoms object and pymatgen Structure objects.
"""

from __future__ import annotations

import warnings
from copy import deepcopy
from importlib.metadata import PackageNotFoundError
from typing import TYPE_CHECKING, TypeVar

import numpy as np
from monty.json import MontyDecoder, MSONable, jsanitize

from pymatgen.core.structure import IMolecule, IStructure, Lattice, Molecule, Structure

try:
    from ase.atoms import Atoms
    from ase.calculators.singlepoint import SinglePointDFTCalculator
    from ase.constraints import FixAtoms, FixCartesian
    from ase.io.jsonio import decode, encode
    from ase.spacegroup import Spacegroup

    NO_ASE_ERR = None

except ImportError:
    NO_ASE_ERR = PackageNotFoundError("AseAtomsAdaptor requires the ASE package. Use `pip install ase`")
    encode = decode = FixAtoms = FixCartesian = SinglePointDFTCalculator = Spacegroup = None

    class Atoms:  # type: ignore[no-redef]
        def __init__(self, *args, **kwargs):
            raise NO_ASE_ERR


if TYPE_CHECKING:
    from typing import Any

    from numpy.typing import ArrayLike
    from typing_extensions import Self

    from pymatgen.core.structure import SiteCollection

__author__ = "Shyue Ping Ong, Andrew S. Rosen"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 8, 2012"

StructT = TypeVar("StructT", bound=IStructure | IMolecule | Structure | Molecule)
MolT = TypeVar("MolT", bound=IMolecule)


class MSONAtoms(Atoms, MSONable):
    """A custom subclass of ASE Atoms that is MSONable, including `.as_dict()` and `.from_dict()` methods."""

    def as_dict(self: Atoms) -> dict[str, Any]:
        # Normally, we would want to this to be a wrapper around atoms.todict() with @module and
        # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() is not meant
        # to be used in a round-trip fashion and does not work properly with constraints.
        # See ASE issue #1387.
        atoms_no_info = self.copy()
        atoms_no_info.info = {}
        return {
            "@module": "pymatgen.io.ase",
            "@class": "MSONAtoms",
            "atoms_json": encode(atoms_no_info),
            "atoms_info": jsanitize(self.info, strict=True),
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        # Normally, we would want to this to be a wrapper around atoms.fromdict() with @module and
        # @class key-value pairs inserted. However, atoms.todict()/atoms.fromdict() is not meant
        # to be used in a round-trip fashion and does not work properly with constraints.
        # See ASE issue #1387.
        mson_atoms = cls(decode(dct["atoms_json"]))
        atoms_info = MontyDecoder().process_decoded(dct.get("atoms_info", {}))
        mson_atoms.info = atoms_info
        return mson_atoms


# NOTE: If making notable changes to this class, please ping @Andrew-S-Rosen on GitHub.
# There are some subtleties in here, particularly related to spins/charges.
class AseAtomsAdaptor:
    """Adaptor serves as a bridge between ASE Atoms and pymatgen objects."""

    @staticmethod
    def get_atoms(
        structure: SiteCollection,
        msonable: bool = True,
        **kwargs,
    ) -> MSONAtoms | Atoms:
        """Get ASE Atoms object from pymatgen structure or molecule.

        Args:
            structure (SiteCollection): pymatgen Structure or Molecule
            msonable (bool): Whether to return an MSONAtoms object, which is MSONable.
            **kwargs: passed to the ASE Atoms constructor

        Returns:
            Atoms: ASE Atoms object
        """
        if NO_ASE_ERR is not None:
            raise NO_ASE_ERR
        if not structure.is_ordered:
            raise ValueError("ASE Atoms only supports ordered structures")

        # Construct the base ASE Atoms object
        symbols = [str(site.specie.symbol) for site in structure]
        positions = [site.coords for site in structure]
        if hasattr(structure, "lattice"):
            pbc = getattr(structure.lattice, "pbc", True)
            cell = structure.lattice.matrix
        else:  # Molecule without lattice
            pbc = False
            cell = None

        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell, **kwargs)

        if msonable:
            atoms = MSONAtoms(atoms)

        if tags := structure.site_properties.get("tags"):
            atoms.set_tags(tags)

        # Set the site magmoms in the ASE Atoms object
        # Note: ASE distinguishes between initial and converged
        # magnetic moment site properties, whereas pymatgen does not. Therefore, we
        # have to distinguish between these two when constructing the Structure/Molecule.
        # Additionally, ASE does not have a global .charge or .spin_multiplicity property,
        # whereas Molecule does. Therefore, we have to patch these in to ensure they do not
        # get lost during interconversion.
        #
        # The mapping selected here is as follows:
        #
        # Site properties:
        # Atoms.get_initial_magnetic_moments() <--> Structure/Molecule.site_properties["magmom"]
        # Atoms.get_magnetic_moments() <--> Structure/Molecule.site_properties["final_magmom"]
        # Atoms.get_initial_charges() <--> Structure/Molecule.site_properties["charge"]
        # Atoms.get_charges() <--> Structure/Molecule.site_properties["final_charge"]
        #
        # Global properties:
        # Atoms.charge <--> Molecule.charge
        # Atoms.spin_multiplicity <--> Molecule.spin_multiplicity

        # Set Atoms initial magnetic moments and charges
        if "charge" in structure.site_properties:
            initial_charges = structure.site_properties["charge"]
            atoms.set_initial_charges(initial_charges)

        if "magmom" in structure.site_properties:
            initial_magmoms = structure.site_properties["magmom"]
            atoms.set_initial_magnetic_moments(initial_magmoms)

        # Set Atoms global charge and spin multiplicity.
        # This is patched into the Atoms object to ensure we don't lose the
        # Pymatgen global charge/spin_multiplicity when interconverting.
        if isinstance(structure, Molecule):
            atoms.charge = structure.charge
            atoms.spin_multiplicity = structure.spin_multiplicity

        # Get the oxidation states from the structure
        oxi_states: list[float | None] = [getattr(site.specie, "oxi_state", None) for site in structure]

        # Read in selective dynamics if present.
        # Note that FixCartesian class uses an opposite notion of
        # "fix" and "not fix" flags: in ASE True means fixed and False
        # means not fixed.
        fix_atoms: dict | None = None
        if "selective_dynamics" in structure.site_properties:
            fix_atoms = {
                str([xc, yc, zc]): ([xc, yc, zc], [])
                for xc in (True, False)
                for yc in (True, False)
                for zc in (True, False)
            }
            # [False, False, False] is free to move - no constraint in ASE.
            del fix_atoms[str([False, False, False])]
            for site in structure:
                selective_dynamics: ArrayLike = site.properties.get("selective_dynamics")  # type: ignore[assignment]
                for cmask_str in fix_atoms:
                    cmask_site = (~np.array(selective_dynamics)).tolist()
                    fix_atoms[cmask_str][1].append(cmask_str == str(cmask_site))
        else:
            fix_atoms = None

        # Set the selective dynamics with the FixCartesian class.
        if fix_atoms is not None:
            atoms.set_constraint(
                [
                    FixAtoms(indices) if cmask == [True, True, True] else FixCartesian(indices, mask=cmask)
                    for cmask, indices in fix_atoms.values()
                    # Do not add empty constraints
                    if any(indices)
                ]
            )

        # Add any remaining site properties to the ASE Atoms object
        for prop in structure.site_properties:
            if prop not in {
                "magmom",
                "charge",
                "final_magmom",
                "final_charge",
                "selective_dynamics",
            }:
                atoms.set_array(prop, np.array(structure.site_properties[prop]))
        if any(oxi_states):
            atoms.set_array("oxi_states", np.array(oxi_states))

        # Atoms.info <---> Structure.properties
        if properties := structure.properties:
            atoms.info = deepcopy(properties)

        # Regenerate Spacegroup object from `.todict()` representation
        if isinstance(atoms.info.get("spacegroup"), dict):
            atoms.info["spacegroup"] = Spacegroup(
                atoms.info["spacegroup"]["number"],
                setting=atoms.info["spacegroup"].get("setting", 1),
            )

        # Atoms.calc <---> Structure.calc
        if calc := getattr(structure, "calc", None):
            atoms.calc = calc
        else:
            # Set the Atoms final magnetic moments and charges if present.
            # This uses the SinglePointDFTCalculator as the dummy calculator
            # to store results.
            charges = structure.site_properties.get("final_charge")
            magmoms = structure.site_properties.get("final_magmom")
            if charges or magmoms:
                calc = SinglePointDFTCalculator(atoms, magmoms=magmoms, charges=charges)
                atoms.calc = calc

        return atoms

    @staticmethod
    def get_structure(
        atoms: Atoms,
        cls=Structure,
        **cls_kwargs,
    ) -> Structure | Molecule:
        """Get pymatgen structure from ASE Atoms.

        Args:
            atoms (Atoms): ASE Atoms object
            cls: The structure class to instantiate (defaults to pymatgen Structure)
            **cls_kwargs: Any additional kwargs to pass to the cls constructor

        Returns:
            (I)Structure/(I)Molecule: Equivalent pymatgen (I)Structure/(I)Molecule
        """
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        lattice = atoms.get_cell()

        # Get the tags
        tags = atoms.get_tags() if atoms.has("tags") else None

        # Get the (final) site magmoms and charges from the ASE Atoms object.
        if getattr(atoms, "calc", None) is not None and getattr(atoms.calc, "results", None) is not None:
            charges = atoms.calc.results.get("charges")
            magmoms = atoms.calc.results.get("magmoms")
        else:
            magmoms = charges = None

        # Get the initial magmoms and charges from the ASE Atoms object.
        initial_charges = atoms.get_initial_charges() if atoms.has("initial_charges") else None
        initial_magmoms = atoms.get_initial_magnetic_moments() if atoms.has("initial_magmoms") else None
        oxi_states = atoms.get_array("oxi_states") if atoms.has("oxi_states") else None

        # If the ASE Atoms object has constraints, make sure that they are of the
        # FixAtoms or FixCartesian kind, which are the only ones that
        # can be supported in Pymatgen.
        # By default, FixAtoms fixes all three (x, y, z) dimensions.
        if atoms.constraints:
            unsupported_constraint_type = False
            constraint_indices: dict = {
                str([xc, yc, zc]): ([xc, yc, zc], [])
                for xc in (True, False)
                for yc in (True, False)
                for zc in (True, False)
            }
            for constraint in atoms.constraints:
                if isinstance(constraint, FixAtoms):
                    constraint_indices[str([False] * 3)][1].extend(constraint.get_indices().tolist())
                elif isinstance(constraint, FixCartesian):
                    cmask = (~np.array(constraint.mask)).tolist()
                    constraint_indices[str(cmask)][1].extend(constraint.get_indices().tolist())
                else:
                    unsupported_constraint_type = True
            if unsupported_constraint_type:
                warnings.warn(
                    "Only FixAtoms and FixCartesian is supported by Pymatgen. Other constraints will not be set.",
                    stacklevel=2,
                )
            sel_dyn = []
            for atom in atoms:
                constrained = False
                for mask, indices in constraint_indices.values():
                    if atom.index in indices:
                        sel_dyn.append(mask)
                        constrained = True
                        break  # Assume no duplicates
                if not constrained:
                    sel_dyn.append([True] * 3)
        else:
            sel_dyn = None

        # Atoms.info <---> Structure.properties
        properties = deepcopy(getattr(atoms, "info", {}))
        # If present, convert Spacegroup object to JSON-serializable dict
        if properties.get("spacegroup") and isinstance(properties["spacegroup"], Spacegroup):
            properties["spacegroup"] = properties["spacegroup"].todict()

        # Return a (I)Molecule object if that was specifically requested;
        # otherwise return a (I)Structure object as expected
        if issubclass(cls, IMolecule):
            structure = cls(symbols, positions, properties=properties, **cls_kwargs)

        elif issubclass(cls, IStructure):
            structure = cls(
                Lattice(lattice, pbc=atoms.pbc),
                symbols,
                positions,
                coords_are_cartesian=True,
                properties=properties,
                **cls_kwargs,
            )

        else:
            raise TypeError(f"Unsupported {cls=}")

        # Atoms.calc <---> Structure.calc
        if calc := getattr(atoms, "calc", None):
            structure.calc = calc

        # Set the site magmoms in the Pymatgen structure object
        # Note: ASE distinguishes between initial and converged
        # magnetic moment site properties, whereas pymatgen does not. Therefore, we
        # have to distinguish between these two when constructing the Structure/Molecule.
        # Additionally, ASE does not have a global .charge or .spin_multiplicity property,
        # whereas Molecule does. Therefore, we have to patch these in to ensure they do not
        # get lost during interconversion.
        #
        # The mapping selected here is:
        #
        # Site properties:
        # Atoms.get_initial_magnetic_moments() <--> Structure/Molecule.site_properties["magmom"]
        # Atoms.get_magnetic_moments() <--> Structure/Molecule.site_properties["final_magmom"]
        # Atoms.get_initial_charges() <--> Structure/Molecule.site_properties["charge"]
        # Atoms.get_charges() <--> Structure/Molecule.site_properties["final_charge"]
        #
        # Global properties:
        # Atoms.charge <--> Molecule.charge
        # Atoms.spin_multiplicity <--> Molecule.spin_multiplicity

        if initial_charges is not None:
            structure.add_site_property("charge", initial_charges)
        if charges is not None:
            structure.add_site_property("final_charge", charges)
        if magmoms is not None:
            structure.add_site_property("final_magmom", magmoms)
        if initial_magmoms is not None:
            structure.add_site_property("magmom", initial_magmoms)
        if sel_dyn is not None and ~np.all(sel_dyn):
            structure.add_site_property("selective_dynamics", sel_dyn)
        if tags is not None:
            structure.add_site_property("tags", tags)

        # Add oxidation states by site
        if oxi_states is not None:
            structure.add_oxidation_state_by_site(oxi_states)

        # Add any remaining site properties to the Pymatgen structure object
        for prop in atoms.arrays:
            if prop not in {
                "numbers",
                "positions",
                "magmom",
                "initial_charges",
                "initial_magmoms",
                "final_magmom",
                "charge",
                "final_charge",
                "oxi_states",
            }:
                structure.add_site_property(prop, atoms.get_array(prop).tolist())

        return structure

    @staticmethod
    def get_molecule(atoms: Atoms, cls: type[MolT] = Molecule, **cls_kwargs) -> Molecule | IMolecule:  # type:ignore[assignment]
        """Get pymatgen molecule from ASE Atoms.

        Args:
            atoms (Atom): ASE Atoms object
            cls: The molecule class to instantiate (defaults to pymatgen Molecule)
            **cls_kwargs: Any additional kwargs to pass to the cls constructor

        Returns:
            (I)Molecule: Equivalent pymatgen (I)Molecule
        """
        molecule = AseAtomsAdaptor.get_structure(atoms, cls=cls, **cls_kwargs)

        # Set the global charge and spin multiplicity
        try:
            charge = atoms.charge
        except AttributeError:
            charge = round(np.sum(atoms.get_initial_charges())) if atoms.has("initial_charges") else 0

        try:
            spin_mult = atoms.spin_multiplicity
        except AttributeError:
            spin_mult = round(np.sum(atoms.get_initial_magnetic_moments())) + 1 if atoms.has("initial_magmoms") else 1

        molecule.set_charge_and_spin(charge, spin_multiplicity=spin_mult)

        return molecule  # type:ignore[return-value]
