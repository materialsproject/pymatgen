# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import warnings
import os

import numpy as np

from monty.serialization import loadfn
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
bond_params = loadfn(os.path.join(MODULE_DIR, 'DLS_bond_params.yaml'))


def is_ox(structure):
    comp = structure.composition
    for k in comp.keys():
        try:
            k.oxi_state
        except AttributeError:
            return False
    return True


class RLSVolumePredictor:
    """
    Reference lattice scaling (RLS) scheme that predicts the volume of a
    structure based on a known crystal structure.
    """
    def __init__(self, check_isostructural=True, radii_type="ionic-atomic",
                 use_bv=True):
        """
        Args:
            check_isostructural: Whether to test that the two structures are
                isostructural. This algo works best for isostructural compounds.
                Defaults to True.
            radii_type (str): Types of radii to use. You can specify "ionic"
                (only uses ionic radii), "atomic" (only uses atomic radii) or
                "ionic-atomic" (uses either ionic or atomic radii, with a
                preference for ionic where possible).
            use_bv (bool): Whether to use BVAnalyzer to determine oxidation
                states if not present.
        """
        self.check_isostructural = check_isostructural
        self.radii_type = radii_type
        self.use_bv = use_bv

    def predict(self, structure, ref_structure):
        """
        Given a structure, returns the predicted volume.
        Args:
            structure (Structure): structure w/unknown volume
            ref_structure (Structure): A reference structure with a similar
                structure but different species.
        Returns:
            a float value of the predicted volume
        """

        if self.check_isostructural:
            m = StructureMatcher()
            mapping = m.get_best_electronegativity_anonymous_mapping(
                structure, ref_structure)
            if mapping is None:
                raise ValueError("Input structures do not match!")

        if "ionic" in self.radii_type:
            try:
                # Use BV analyzer to determine oxidation states only if the
                # oxidation states are not already specified in the structure
                # and use_bv is true.
                if (not is_ox(structure)) and self.use_bv:
                    a = BVAnalyzer()
                    structure = a.get_oxi_state_decorated_structure(structure)
                if (not is_ox(ref_structure)) and self.use_bv:
                    a = BVAnalyzer()
                    ref_structure = a.get_oxi_state_decorated_structure(
                        ref_structure)

                comp = structure.composition
                ref_comp = ref_structure.composition

                # Check if all the associated ionic radii are available.
                if any([k.ionic_radius is None for k in list(comp.keys())]) or \
                        any([k.ionic_radius is None for k in
                             list(ref_comp.keys())]):
                    raise ValueError("Not all the ionic radii are available!")

                numerator = 0
                denominator = 0
                # Here, the 1/3 factor on the composition accounts for atomic
                # packing. We want the number per unit length.
                for k, v in comp.items():
                    numerator += k.ionic_radius * v ** (1 / 3)
                for k, v in ref_comp.items():
                    denominator += k.ionic_radius * v ** (1 / 3)

                return ref_structure.volume * (numerator / denominator) ** 3
            except Exception as ex:
                warnings.warn("Exception occured. Will attempt atomic radii.")
                # If error occurs during use of ionic radii scheme, pass
                # and see if we can resolve it using atomic radii.
                pass

        if "atomic" in self.radii_type:
            comp = structure.composition
            ref_comp = ref_structure.composition
            # Here, the 1/3 factor on the composition accounts for atomic
            # packing. We want the number per unit length.
            numerator = 0
            denominator = 0
            for k, v in comp.items():
                numerator += k.atomic_radius * v ** (1 / 3)
            for k, v in ref_comp.items():
                denominator += k.atomic_radius * v ** (1 / 3)
            return ref_structure.volume * (numerator / denominator) ** 3

        raise ValueError("Cannot find volume scaling based on radii choices "
                         "specified!")

    def get_predicted_structure(self, structure, ref_structure):
        """
        Given a structure, returns back the structure scaled to predicted
        volume.
        Args:
            structure (Structure): structure w/unknown volume
            ref_structure (Structure): A reference structure with a similar
                structure but different species.
        Returns:
            a Structure object with predicted volume
        """
        new_structure = structure.copy()
        new_structure.scale_lattice(self.predict(structure, ref_structure))

        return new_structure


class DLSVolumePredictor:
    """
    Data-mined lattice scaling (DLS) scheme that relies on data-mined bond
    lengths to predict the crystal volume of a given structure.

    As of 2/12/19, we suggest this method be used in conjunction with
    min_scaling and max_scaling to prevent instances of very large, unphysical
    predicted volumes found in a small subset of structures.
    """

    def __init__(self, cutoff=4.0, min_scaling=0.5, max_scaling=1.5):
        """
        Args:
            cutoff (float): cutoff radius added to site radius for finding
                site pairs. Necessary to increase only if your initial
                structure guess is extremely bad (atoms way too far apart). In
                all other instances, increasing cutoff gives same answer
                but takes more time.
            min_scaling (float): if not None, this will ensure that the new
                volume is at least this fraction of the original (preventing
                too-small volumes)
            max_scaling (float): if not None, this will ensure that the new
                volume is at most this fraction of the original (preventing
                too-large volumes)
        """
        self.cutoff = cutoff
        self.min_scaling = min_scaling
        self.max_scaling = max_scaling

    def predict(self, structure, icsd_vol=False):
        """
        Given a structure, returns the predicted volume.

        Args:
            structure (Structure) : a crystal structure with an unknown volume.
            icsd_vol (bool) : True if the input structure's volume comes from
                ICSD.

        Returns:
            a float value of the predicted volume.
        """

        # Get standard deviation of electronnegativity in the structure.
        std_x = np.std([site.specie.X for site in structure])
        # Sites that have atomic radii
        sub_sites = []
        # Record the "DLS estimated radius" from bond_params.
        bp_dict = {}

        for sp in list(structure.composition.keys()):
            if sp.atomic_radius:
                sub_sites.extend([site for site in structure
                                  if site.specie == sp])
            else:
                warnings.warn("VolumePredictor: no atomic radius data for "
                              "{}".format(sp))

            if sp.symbol not in bond_params:
                warnings.warn("VolumePredictor: bond parameters not found, "
                              "used atomic radii for {}".format(sp))
            else:
                r, k = bond_params[sp.symbol]["r"], bond_params[sp.symbol]["k"]
                bp_dict[sp] = float(r) + float(k) * std_x

        # Structure object that include only sites with known atomic radii.
        reduced_structure = Structure.from_sites(sub_sites)
        smallest_ratio = None

        for site1 in reduced_structure:
            sp1 = site1.specie
            neighbors = reduced_structure.get_neighbors(site1,
                                                        sp1.atomic_radius +
                                                        self.cutoff)

            for site2, dist in neighbors:
                sp2 = site2.specie

                if sp1 in bp_dict and sp2 in bp_dict:
                    expected_dist = bp_dict[sp1] + bp_dict[sp2]
                else:
                    expected_dist = sp1.atomic_radius + sp2.atomic_radius

                if not smallest_ratio or dist / expected_dist < smallest_ratio:
                    smallest_ratio = dist / expected_dist

        if not smallest_ratio:
            raise ValueError("Could not find any bonds within the given cutoff "
                             "in this structure.")

        volume_factor = (1 / smallest_ratio) ** 3

        # icsd volume fudge factor
        if icsd_vol:
            volume_factor *= 1.05

        if self.min_scaling:
            volume_factor = max(self.min_scaling, volume_factor)
        if self.max_scaling:
            volume_factor = min(self.max_scaling, volume_factor)

        return structure.volume * volume_factor

    def get_predicted_structure(self, structure, icsd_vol=False):
        """
        Given a structure, returns back the structure scaled to predicted
        volume.
        Args:
            structure (Structure): structure w/unknown volume

        Returns:
            a Structure object with predicted volume
        """
        new_structure = structure.copy()
        new_structure.scale_lattice(self.predict(structure, icsd_vol=icsd_vol))
        return new_structure
