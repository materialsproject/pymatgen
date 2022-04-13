# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.
"""


__author__ = "Shyue Ping Ong, Andrew S. Rosen"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 8, 2012"

import warnings

import numpy as np

from pymatgen.core.structure import Molecule, Structure

try:
    from ase import Atoms

    ase_loaded = True
except ImportError:
    ase_loaded = False

if ase_loaded:
    from ase.calculators.singlepoint import SinglePointDFTCalculator
    from ase.constraints import FixAtoms


class AseAtomsAdaptor:
    """
    Adaptor serves as a bridge between ASE Atoms and pymatgen objects.
    """

    @staticmethod
    def get_atoms(structure, **kwargs):
        """
        Returns ASE Atoms object from pymatgen structure or molecule.

        Args:
            structure: pymatgen.core.structure.Structure or pymatgen.core.structure.Molecule
            **kwargs: other keyword args to pass into the ASE Atoms constructor

        Returns:
            ASE Atoms object
        """
        if not structure.is_ordered:
            raise ValueError("ASE Atoms only supports ordered structures")
        if not ase_loaded:
            raise ImportError(
                "AseAtomsAdaptor requires the ASE package.\nUse `pip install ase` or `conda install ase -c conda-forge`"
            )

        # Construct the base ASE Atoms object
        symbols = [str(site.specie.symbol) for site in structure]
        positions = [site.coords for site in structure]
        if hasattr(structure, "lattice"):
            cell = structure.lattice.matrix
            pbc = True
        else:
            cell = None
            pbc = None

        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell, **kwargs)

        # Set the site magmoms in the ASE Atoms object
        # Note: ASE distinguishes between initial and converged
        # magnetic moment site properties, whereas pymatgen does not. Therefore, we
        # have to distinguish between these two when constructing the Structure/Molecule.
        # The mapping selected here is:
        # ASE initial magmom <--> Pymatgen "magmom"
        # ASE final magmom <--> Pymatgen "final_magmom"
        # ASE initial charge <--> Pymatgen "charge"
        # ASE final charge <--> Pymatgen "final_charge"
        if "magmom" in structure.site_properties:
            initial_magmoms = structure.site_properties["magmom"]
            atoms.set_initial_magnetic_moments(initial_magmoms)
        if "charge" in structure.site_properties:
            initial_charges = structure.site_properties["charge"]
            atoms.set_initial_charges(initial_charges)

        if "final_magmom" in structure.site_properties:
            magmoms = structure.site_properties["final_magmom"]
        else:
            magmoms = None
        if "final_charge" in structure.site_properties:
            charges = structure.site_properties["final_charge"]
        else:
            charges = None
        if magmoms or charges:
            if magmoms and charges:
                calc = SinglePointDFTCalculator(atoms, **{"magmoms": magmoms, "charges": charges})
            elif magmoms:
                calc = SinglePointDFTCalculator(atoms, **{"magmoms": magmoms})
            elif charges:
                calc = SinglePointDFTCalculator(atoms, **{"charges": charges})
            atoms.calc = calc

        # Get the oxidation states from the structure
        oxi_states = [getattr(site.specie, "oxi_state", None) for site in structure]

        # Read in selective dynamics if present. Note that the ASE FixAtoms class fixes (x,y,z), so
        # here we make sure that [False, False, False] or [True, True, True] is set for the site selective
        # dynamics property. As a consequence, if only a subset of dimensions are fixed, this won't get passed to ASE.
        if "selective_dynamics" in structure.site_properties:
            fix_atoms = []
            for site in structure:
                site_prop = site.properties["selective_dynamics"]
                if site_prop not in [[True, True, True], [False, False, False]]:
                    raise ValueError(
                        "ASE FixAtoms constraint does not support selective dynamics in only some dimensions."
                        "Remove the selective dynamics and try again if you do not need them."
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
        if np.any(oxi_states):
            atoms.set_array("oxi_states", np.array(oxi_states))

        return atoms

    @staticmethod
    def get_structure(atoms, cls=None, **cls_kwargs):
        """
        Returns pymatgen structure from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Structure class to instantiate (defaults to pymatgen structure)
            **cls_kwargs: Any additional kwargs to pass to the cls

        Returns:
            Equivalent pymatgen.core.structure.Structure
        """

        cls = Structure if cls is None else cls

        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        lattice = atoms.get_cell()

        # Get the site magmoms from the ASE Atoms objects.
        if getattr(atoms, "calc", None) is not None and getattr(atoms.calc, "results", None) is not None:
            magmoms = atoms.calc.results.get("magmoms", None)
            charges = atoms.calc.results.get("charges", None)
        else:
            magmoms = None
            charges = None

        if atoms.has("initial_magmoms"):
            initial_magmoms = atoms.get_initial_magnetic_moments()
        else:
            initial_magmoms = None
        if atoms.has("initial_charges"):
            initial_charges = atoms.get_initial_charges()
        else:
            initial_charges = None
        if atoms.has("oxi_states"):
            oxi_states = atoms.get_array("oxi_states")
        else:
            oxi_states = None

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

        # Return a Molecule object if that was specifically requested;
        # otherwise return a Structure object as expected
        if cls == Molecule:
            structure = cls(symbols, positions, **cls_kwargs)
        else:
            structure = cls(lattice, symbols, positions, coords_are_cartesian=True, **cls_kwargs)

        # Set the site magmoms in the Pymatgen structure object
        # Note: ASE distinguishes between initial and converged
        # magnetic moment site properties, whereas pymatgen does not. Therefore, we
        # have to distinguish between these two when constructing the Structure/Molecule.
        # The mapping selected here is:
        # ASE initial magmom <--> Pymatgen "magmom"
        # ASE final magmom <--> Pymatgen "final_magmom"
        # ASE initial charge <--> Pymatgen "charge"
        # ASE final charge <--> Pymatgen "final_charge"
        if magmoms is not None:
            structure.add_site_property("final_magmom", magmoms)
        if initial_magmoms is not None:
            structure.add_site_property("magmom", initial_magmoms)
        if charges is not None:
            structure.add_site_property("final_charge", charges)
        if initial_charges is not None:
            structure.add_site_property("charge", initial_charges)
        if sel_dyn is not None and ~np.all(sel_dyn):
            structure.add_site_property("selective_dynamics", sel_dyn)

        # Add oxidation states by site
        if oxi_states is not None:
            structure.add_oxidation_state_by_site(oxi_states)

        # Add any remaining site properties to the Pymatgen structure object
        for prop in atoms.arrays:
            if prop not in [
                "numbers",
                "positions",
                "magmom",
                "final_magmom",
                "charge",
                "final_charge",
                "oxi_states",
            ]:
                structure.add_site_property(prop, atoms.get_array(prop).tolist())

        return structure

    @staticmethod
    def get_molecule(atoms, cls=None, **cls_kwargs):
        """
        Returns pymatgen molecule from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Molecule class to instantiate (defaults to pymatgen molecule)
            **cls_kwargs: Any additional kwargs to pass to the cls

        Returns:
            Equivalent pymatgen.core.structure.Molecule
        """

        cls = Molecule if cls is None else cls
        molecule = AseAtomsAdaptor.get_structure(atoms, cls=cls, **cls_kwargs)
        if atoms.has("initial_charges"):
            charge = round(np.sum(atoms.get_initial_charges()))
        else:
            charge = 0
        if atoms.has("initial_magmoms"):
            mult = round(np.sum(atoms.get_initial_magnetic_moments())) + 1
        else:
            mult = 1
        molecule.set_charge_and_spin(charge, spin_multiplicity=mult)

        return molecule
