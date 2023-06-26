"""
This module provides conversion between the JARVIS
Atoms object and pymatgen Structure objects.
"""

from __future__ import annotations

from pymatgen.core.structure import Structure

try:
    from jarvis.core.atoms import Atoms
except ImportError:
    Atoms = None


class JarvisAtomsAdaptor:
    """Adaptor serves as a bridge between JARVIS Atoms and pymatgen objects."""

    @staticmethod
    def get_atoms(structure):
        """
        Returns JARVIS Atoms object from pymatgen structure.

        Args:
            structure: pymatgen.core.structure.Structure

        Returns:
            JARVIS Atoms object
        """
        if not structure.is_ordered:
            raise ValueError("JARVIS Atoms only supports ordered structures")
        if Atoms is None:
            raise ImportError("JarvisAtomsAdaptor requires jarvis-tools package.\nUse `pip install -U jarvis-tools`")
        elements = [str(site.specie.symbol) for site in structure]
        coords = [site.frac_coords for site in structure]
        return Atoms(
            lattice_mat=structure.lattice.matrix,
            elements=elements,
            coords=coords,
            cartesian=False,
        )

    @staticmethod
    def get_structure(atoms):
        """
        Returns pymatgen structure from JARVIS Atoms.

        Args:
            atoms: JARVIS Atoms object

        Returns:
            Equivalent pymatgen.core.structure.Structure
        """
        return Structure(
            lattice=atoms.lattice_mat, species=atoms.elements, coords=atoms.frac_coords, coords_are_cartesian=False
        )
