"""
This module implements a MolecularOrbital class to represent band character in
solids. Useful for predicting PDOS character from structural information.
"""

from __future__ import annotations

from itertools import chain, combinations

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element


class MolecularOrbitals:
    """
    Represents the character of bands in a solid. The input is a chemical
    formula, since no structural characteristics are taken into account.

    The band character of a crystal emerges from the atomic orbitals of the
    constituent ions, hybridization/covalent bonds, and the spin-orbit
    interaction (ex: Fe2O3). Right now the orbitals are only built from
    the uncharged atomic species. Functionality can be improved by:
    1) calculate charged ion orbital energies
    2) incorporate the coordination environment to account for covalent bonds

    The atomic orbital energies are stored in pymatgen.core.periodic_table.JSON

    MOs = MolecularOrbitals('SrTiO3')
    MOs.band_edges
    # gives {'HOMO':['O','2p',-0.338381], 'LUMO':['Ti','3d',-0.17001], 'metal':False}
    """

    def __init__(self, formula):
        """
        Args:
            formula (str): Chemical formula. Must have integer subscripts. Ex: 'SrTiO3'.

        Attributes:
            composition: the composition as a dictionary. Ex: {'Sr': 1, 'Ti': 1, 'O', 3}
            elements: the dictionary keys for the composition
            elec_neg: the maximum pairwise electronegativity difference
            aos: the constituent atomic orbitals for each element as a dictionary
            band_edges: dictionary containing the highest occupied molecular orbital (HOMO),
                lowest unoccupied molecular orbital (LUMO), and whether the material is predicted
                to be a metal
        """
        self.composition = Composition(formula).as_dict()
        self.elements = list(self.composition)
        for subscript in self.composition.values():
            if not float(subscript).is_integer():
                raise ValueError("composition subscripts must be integers")

        self.elec_neg = self.max_electronegativity()
        self.aos = {
            str(el): [[str(el), k, v] for k, v in Element(el).atomic_orbitals.items()]  # pylint: disable=E1101
            for el in self.elements
        }
        self.band_edges = self.obtain_band_edges()

    def max_electronegativity(self):
        """
        Returns:
            The maximum pairwise electronegativity difference.
        """
        maximum = 0
        for e1, e2 in combinations(self.elements, 2):
            if abs(Element(e1).X - Element(e2).X) > maximum:
                maximum = abs(Element(e1).X - Element(e2).X)
        return maximum

    def aos_as_list(self):
        """
        Returns:
            A list of atomic orbitals, sorted from lowest to highest energy.

            The orbitals energies in eV are represented as
                [['O', '1s', -18.758245], ['O', '2s', -0.871362], ['O', '2p', -0.338381]]
            Data is obtained from
            https://www.nist.gov/pml/data/atomic-reference-data-electronic-structure-calculations
        """
        return sorted(
            chain.from_iterable([self.aos[el] * int(self.composition[el]) for el in self.elements]),
            key=lambda x: x[2],
        )

    def obtain_band_edges(self):
        """
        Fill up the atomic orbitals with available electrons.

        Returns:
            HOMO, LUMO, and whether it's a metal.
        """
        orbitals = self.aos_as_list()
        electrons = Composition(self.composition).total_electrons
        partial_filled = []
        for orbital in orbitals:
            if electrons <= 0:
                break
            if "s" in orbital[1]:
                electrons += -2
            elif "p" in orbital[1]:
                electrons += -6
            elif "d" in orbital[1]:
                electrons += -10
            elif "f" in orbital[1]:
                electrons += -14
            partial_filled.append(orbital)

        if electrons != 0:
            homo = partial_filled[-1]
            lumo = partial_filled[-1]
        else:
            homo = partial_filled[-1]
            try:
                lumo = orbitals[len(partial_filled)]
            except Exception:
                lumo = None

        return {"HOMO": homo, "LUMO": lumo, "metal": homo == lumo}
