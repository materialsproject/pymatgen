"""
This module provides conversion between the JARVIS
Atoms object and pymatgen Structure objects.
"""

from pymatgen.core.structure import Structure

try:
    from jarvis.core.atoms import Atoms

    jarvis_loaded = True
except ImportError:
    jarvis_loaded = False


class JarvisAtomsAdaptor:
    """
    Adaptor serves as a bridge between JARVIS Atoms and pymatgen objects.
    """

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
        if not jarvis_loaded:
            raise ImportError("JarvisAtomsAdaptor requires jarvis-tools package.\n" "Use `pip install -U jarvis-tools`")
        elements = [str(site.specie.symbol) for site in structure]
        coords = [site.frac_coords for site in structure]
        lattice_mat = structure.lattice.matrix
        return Atoms(
            lattice_mat=lattice_mat,
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
        elements = atoms.elements
        coords = atoms.frac_coords
        lattice_mat = atoms.lattice_mat

        return Structure(lattice_mat, elements, coords, coords_are_cartesian=False)
