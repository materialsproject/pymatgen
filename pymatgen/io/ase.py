"""
This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.
"""


from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from monty.json import jsanitize

from pymatgen.core.structure import Molecule, Structure

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from pymatgen.core.structure import SiteCollection

try:
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointDFTCalculator
    from ase.constraints import FixAtoms

    ase_loaded = True
except ImportError:
    ase_loaded = False

__author__ = "Shyue Ping Ong, Andrew S. Rosen"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 8, 2012"


# NOTE: If making notable changes to this class, please ping @Andrew-S-Rosen on GitHub.
# There are some subtleties in here, particularly related to spins/charges.
class AseAtomsAdaptor:
    """Adaptor serves as a bridge between ASE Atoms and pymatgen objects."""

    @staticmethod
    def get_atoms(structure: SiteCollection, **kwargs) -> Atoms:
        """
        Returns ASE Atoms object from pymatgen structure or molecule.

        Args:
            structure (SiteCollection): pymatgen Structure or Molecule
            **kwargs: passed to the ASE Atoms constructor

        Returns:
            Atoms: ASE Atoms object
        """
        if not ase_loaded:
            raise ImportError(
                "AseAtomsAdaptor requires the ASE package.\n"
                "Use `pip install ase` or `conda install ase -c conda-forge`"
            )
        if not structure.is_ordered:
            raise ValueError("ASE Atoms only supports ordered structures")

        # Construct the base ASE Atoms object
        symbols = [str(site.specie.symbol) for site in structure]
        positions = [site.coords for site in structure]
        if hasattr(structure, "lattice"):
            pbc = True
            cell = structure.lattice.matrix
        else:  # Molecule without lattice
            pbc = False
            cell = None

        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell, **kwargs)

        if "tags" in structure.site_properties:
            atoms.set_tags(structure.site_properties["tags"])

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

        # Set the Atoms final magnetic moments and charges if present.
        # This uses the SinglePointDFTCalculator as the dummy calculator
        # to store results.
        charges = structure.site_properties.get("final_charge")
        magmoms = structure.site_properties.get("final_magmom")
        if magmoms or charges:
            if magmoms and charges:
                calc = SinglePointDFTCalculator(atoms, magmoms=magmoms, charges=charges)
            elif magmoms:
                calc = SinglePointDFTCalculator(atoms, magmoms=magmoms)
            else:
                calc = SinglePointDFTCalculator(atoms, charges=charges)
            atoms.calc = calc

        # Get the oxidation states from the structure
        oxi_states: list[float | None] = [getattr(site.specie, "oxi_state", None) for site in structure]

        # Read in selective dynamics if present. Note that the ASE FixAtoms class fixes (x,y,z), so
        # here we make sure that [False, False, False] or [True, True, True] is set for the site selective
        # dynamics property. As a consequence, if only a subset of dimensions are fixed, this won't get passed to ASE.
        if "selective_dynamics" in structure.site_properties:
            fix_atoms = []
            for site in structure:
                selective_dynamics: ArrayLike = site.properties.get("selective_dynamics")  # type: ignore[assignment]
                if not (np.all(selective_dynamics) or not np.any(selective_dynamics)):
                    # should be [True, True, True] or [False, False, False]
                    raise ValueError(
                        "ASE FixAtoms constraint does not support selective dynamics in only some dimensions. "
                        f"Remove the {selective_dynamics=} and try again if you do not need them."
                    )
                is_fixed = bool(~np.all(site.properties["selective_dynamics"]))
                fix_atoms.append(is_fixed)

        else:
            fix_atoms = None

        # Set the selective dynamics with the FixAtoms class.
        if fix_atoms is not None:
            atoms.set_constraint(FixAtoms(mask=fix_atoms))

        # Add any remaining site properties to the ASE Atoms object
        for prop in structure.site_properties:
            if prop not in ["magmom", "charge", "final_magmom", "final_charge", "selective_dynamics"]:
                atoms.set_array(prop, np.array(structure.site_properties[prop]))
        if any(oxi_states):
            atoms.set_array("oxi_states", np.array(oxi_states))

        # Atoms.info <---> Structure.properties
        # Atoms.calc <---> Structure.calc
        if properties := getattr(structure, "properties"):  # noqa: B009
            atoms.info = properties
        if calc := getattr(structure, "calc", None):
            atoms.calc = calc

        return atoms

    @staticmethod
    def get_structure(atoms: Atoms, cls: type[Structure] = Structure, **cls_kwargs) -> Structure:
        """
        Returns pymatgen structure from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Structure class to instantiate (defaults to pymatgen Structure)
            **cls_kwargs: Any additional kwargs to pass to the cls

        Returns:
            Equivalent pymatgen.core.structure.Structure
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
        # kind FixAtoms, which are the only ones that can be supported in Pymatgen.
        # By default, FixAtoms fixes all three (x, y, z) dimensions.
        if atoms.constraints:
            unsupported_constraint_type = False
            constraint_indices = []
            for constraint in atoms.constraints:
                if isinstance(constraint, FixAtoms):
                    constraint_indices.extend(constraint.get_indices().tolist())
                else:
                    unsupported_constraint_type = True
            if unsupported_constraint_type:
                warnings.warn("Only FixAtoms is supported by Pymatgen. Other constraints will not be set.")
            sel_dyn = [[False] * 3 if atom.index in constraint_indices else [True] * 3 for atom in atoms]
        else:
            sel_dyn = None

        # Atoms.info <---> Structure.properties
        properties = jsanitize(getattr(atoms, "info", {}))

        # Return a Molecule object if that was specifically requested;
        # otherwise return a Structure object as expected
        if cls == Molecule:
            structure = cls(symbols, positions, properties=properties, **cls_kwargs)
        else:
            structure = cls(lattice, symbols, positions, coords_are_cartesian=True, properties=properties, **cls_kwargs)

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
            if prop not in [
                "numbers",
                "positions",
                "magmom",
                "initial_charges",
                "initial_magmoms",
                "final_magmom",
                "charge",
                "final_charge",
                "oxi_states",
            ]:
                structure.add_site_property(prop, atoms.get_array(prop).tolist())

        return structure

    @staticmethod
    def get_molecule(atoms: Atoms, cls: type[Molecule] = Molecule, **cls_kwargs) -> Molecule:
        """
        Returns pymatgen molecule from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Molecule class to instantiate (defaults to pymatgen molecule)
            **cls_kwargs: Any additional kwargs to pass to the cls

        Returns:
            Molecule: Equivalent pymatgen.core.structure.Molecule
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

        return molecule
