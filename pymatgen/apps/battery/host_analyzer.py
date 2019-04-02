# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from collections import defaultdict
import math
from copy import deepcopy
import os
import logging
import sys
from monty.serialization import loadfn
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Composition
from pymatgen.analysis.structure_matcher import StructureMatcher, ElementComparator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Jimmy Shen"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Jimmy Shen"
__email__ = "jmmshn@lbl.gov"
__date__ = "April 1, 2019"

root = logging.getLogger()
root.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class MigrationPathAnalyzer():
    """
    Methods for analyzing the migration path of a material using the following scheme:
    - Map the relaxed sites of a material back to the empty host lattice
    - Apply symmetry operations of the empty lattice to obtain the other positions of the intercollated atom
    - Get the symmetry inequivalent hops
    - Get the migration barriers for each inequivalent hop

    """

    def __init__(self, base_entry, single_cat_entries, cation='Li', ltol=0.2, stol=0.3, angle_tol=5):
        """
        Pass in a entries for analysis
        Arguments:
            base_entry - the structure without a working ion for us to analyze the migration
            single_cat_entries - list of structures containing a single cation at different positions
            cation - a String symbol or Element for the cation. It must be positively charged, but can be 1+/2+/3+ etc.
        """

        self.sm = StructureMatcher(
            comparator=ElementComparator(),
            primitive_cell=False,
            ignored_species=[cation], ltol=ltol, stol=stol, angle_tol=angle_tol)


        logger.debug('See if the structures all match')
        for ent in single_cat_entries:
            assert(self.sm.fit(base_entry.structure, ent.structure))

        self.base_entry = base_entry
        self.single_cat_entries = single_cat_entries
        self.cation = cation

        self.translated_single_cat_entries = list(map(self.match_ent_to_base, self.single_cat_entries))

    def match_ent_to_base(self, ent):
        new_ent = deepcopy(ent)
        new_struct = self.sm.get_s2_like_s1(self.base_entry.structure, ent.structure)
        new_ent.structure = new_struct
        return new_ent

    def get_all_sym_sites(self, ent):
        sa = SpacegroupAnalyzer(self.base_entry.structure, symprec=0.3, angle_tolerance=10)
        host_allsites = self.base_entry.structure.copy()
        pos_Li = list(filter(lambda isite : isite.species_string=='Li', ent.structure.sites))

        host_allsites.insert(0, 'Li', pos_Li[0].frac_coords)

        for op in sa.get_space_group_operations():
            print("operation")
            struct_tmp = host_allsites.copy()
            struct_tmp.apply_operation(symmop=op,fractional=True)
            for isite in struct_tmp.sites:
                if isite.species_string == "Li":
                    host_allsites.insert(0, 'Li', isite.frac_coords)

        host_allsites.merge_sites(mode='delete')
        return host_allsites
    def get_full_structure(self):
        res = []
        for itr in mpa.translated_single_cat_entries:
            res.extend(mpa.get_all_sym_sites(itr).sites)
        return Structure.from_sites(res)

def main():
    test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            'test_files')
    test_ents = loadfn(test_dir+'/Mn6O5F7_cat_migration.json')

    mpa  = MigrationPathAnalyzer(test_ents['ent_base'], test_ents['one_cation'])
    mpa.get_all_sym_sites(mpa.translated_single_cat_entries[4])

if __name__ == "__main__":
    main()
