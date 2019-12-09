# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides conversion between the Atomic Simulation Environment
Atoms object and pymatgen Structure objects.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 8, 2012"

from pymatgen.core.structure import Molecule, Structure

try:
    from ase import Atoms
    ase_loaded = True
except ImportError:
    ase_loaded = False


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
            raise ImportError("AseAtomsAdaptor requires ase package.\n"
                              "Use `pip install ase` or `conda install ase -c conda-forge`")
        symbols = [str(site.specie.symbol) for site in structure]
        positions = [site.coords for site in structure]
        if hasattr(structure, "lattice"):
            cell = structure.lattice.matrix
            pbc = True
        else:
            cell = None
            pbc = None
        return Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell,
                     **kwargs)

    @staticmethod
    def get_structure(atoms, cls=None):
        """
        Returns pymatgen structure from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Structure class to instantiate (defaults to pymatgen structure)

        Returns:
            Equivalent pymatgen.core.structure.Structure
        """
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        lattice = atoms.get_cell()

        cls = Structure if cls is None else cls
        return cls(lattice, symbols, positions,
                   coords_are_cartesian=True)

    @staticmethod
    def get_molecule(atoms, cls=None):
        """
        Returns pymatgen molecule from ASE Atoms.

        Args:
            atoms: ASE Atoms object
            cls: The Molecule class to instantiate (defaults to pymatgen molecule)

        Returns:
            Equivalent pymatgen.core.structure.Molecule
        """
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()

        cls = Molecule if cls is None else cls
        return cls(symbols, positions)
