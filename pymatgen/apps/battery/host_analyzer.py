# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from pprint import pprint
from collections import defaultdict
import math
from copy import deepcopy
import os
import logging
import sys
from monty.serialization import loadfn
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Composition, Structure
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
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class MigrationPathAnalyzer():
    """
    Methods for analyzing the migration path of a material using the following scheme:
    - Map the relaxed sites of a material back to the empty host lattice
    - Apply symmetry operations of the empty lattice to obtain the other positions of the intercollated atom
    - Get the symmetry inequivalent hops
    - Get the migration barriers for each inequivalent hop

    """

    def __init__(self,
                 base_entry,
                 single_cat_entries,
                 cation='Li',
                 ltol=0.2,
                 stol=0.3,
                 angle_tol=5):
        """
        Pass in a entries for analysis

        :param base_entry: the structure without a working ion for us to analyze the migration
        :param single_cat_entries: list of structures containing a single cation at different positions
        :param cation: a String symbol or Element for the cation. It must be positively charged, but can be 1+/2+/3+ etc.
        :param ltol: parameter for StructureMatcher
        :param stol: parameter for StructureMa
        :param angle_tol: parameter for StructureMa

        """

        self.sm = StructureMatcher(
            comparator=ElementComparator(),
            primitive_cell=False,
            ignored_species=[cation],
            ltol=ltol,
            stol=stol,
            angle_tol=angle_tol)

        logger.debug('See if the structures all match')
        for ent in single_cat_entries:
            assert (self.sm.fit(base_entry.structure, ent.structure))

        self.base_entry = base_entry
        self.single_cat_entries = single_cat_entries
        self.cation = cation

        self.translated_single_cat_entries = list(
            map(self.match_ent_to_base, self.single_cat_entries))

        self.full_sites = self.get_full_sites()
        self.base_structure_full_sites = self.full_sites.copy()
        self.base_structure_full_sites.sites.extend(self.base_entry.structure.sites)

    def match_ent_to_base(self, ent):
        """
        Transform the structure of one entry to match the base structure

        :param ent: inserted structure with cation atoms
        :returns: entry with modified structure
        :rtype: ComputedStructureEntry

        """
        new_ent = deepcopy(ent)
        new_struct = self.sm.get_s2_like_s1(self.base_entry.structure,
                                            ent.structure)
        new_ent.structure = new_struct
        return new_ent

    def get_all_sym_sites(self, ent):
        """
        Return all of the symmetry equivalent sites

        :param ent: ComputedStructureEntry that contains cation
        :returns: Structure containing all of the symmetry equivalent sites
        :rtype: Structure

        """

        sa = SpacegroupAnalyzer(
            self.base_entry.structure, symprec=0.3, angle_tolerance=10)
        host_allsites = self.base_entry.structure.copy()
        host_allsites.remove_species(host_allsites.species)
        pos_Li = list(
            filter(lambda isite: isite.species_string == 'Li',
                   ent.structure.sites))

        for isite in pos_Li:
            host_allsites.insert(0, 'Li', isite.frac_coords, properties={'inserted_energy' : ent.energy})

        for op in sa.get_space_group_operations():
            struct_tmp = host_allsites.copy()
            struct_tmp.apply_operation(symmop=op, fractional=True)
            for isite in struct_tmp.sites:
                if isite.species_string == "Li":
                    host_allsites.insert(0, 'Li', isite.frac_coords, properties={'inserted_energy' : ent.energy})
        for isite in host_allsites.sites:
            logger.info(f"{isite.as_dict()['properties']}")

        host_allsites.merge_sites(mode='average')
        for isite in host_allsites.sites:
            logger.info(f"average {isite.as_dict()['properties']}")
        return host_allsites

    def get_full_sites(self):
        """
        Get each group of symmetry inequivalent sites and combine them

        :returns: Structure with all possible Li sites, the enregy of the structure is stored as a site property
        :rtype: Structure

        """
        res = []
        for itr in self.translated_single_cat_entries:
            res.extend(self.get_all_sym_sites(itr).sites)
        res = Structure.from_sites(res)
        res.merge_sites(tol=1.0, mode='average')
        return res

    def get_hops(self):
        # enumerate all of the hops in a systems
        pass

def main():
    test_dir = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", 'test_files')
    test_ents = loadfn(test_dir + '/Mn6O5F7_cat_migration.json')

    mpa = MigrationPathAnalyzer(test_ents['ent_base'], test_ents['one_cation'])
    print(mpa.full_sites)


if __name__ == "__main__":
    main()
