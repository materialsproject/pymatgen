#!/usr/bin/env python

'''
Created on Feb 24, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Feb 24, 2012"

import logging

from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.analysis.structure_fitter import StructureFitter

logger = logging.getLogger(__name__)


def _get_host(structure, species_to_remove):
    if species_to_remove:
        editor = StructureEditor(structure)
        editor.remove_species(species_to_remove)
        return editor.modified_structure
    else:
        return structure


def group_entries_by_structure(entries, species_to_remove=None,
                               anonymized=False, symmetry_tol=0,
                               fitting_accuracy=StructureFitter.FAST_FIT):
    """
    Given a sequence of ComputedStructureEntries, use structure fitter to group
    them by structural similarity.

    Args:
        entries:
            Sequence of entries.
        species_to_remove:
            Sometimes you want to compare a host framework (e.g., in Li-ion
            battery analysis). This allows you to specify species to remove
            before structural comparison.
        anonymized:
            Allow anonymized fitting.
        fitting_accuracy:
            The fitting accuracy to use for StructureFitter. Defaults to
            FAST_FIT.

    Returns:
        Sequence of sequence of entries by structural similarity. e.g,
        [[ entry1, entry2], [entry3, entry4, entry5]]
    """
    all_matches = []
    unmatched = [(entry, _get_host(entry.structure, species_to_remove))
                 for entry in entries]
    while len(unmatched) > 0:
        ref_host = unmatched[0][1]
        logger.info(
            "Reference tid = {}, formula = {}".format(unmatched[0][0].entry_id,
                                                      ref_host.formula)
        )
        ref_formula = ref_host.composition.anonymized_formula if anonymized \
            else ref_host.composition.reduced_formula
        logger.info("Reference host = {}".format(ref_formula))
        matches = [unmatched[0]]
        for i in xrange(1, len(unmatched)):
            test_host = unmatched[i][1]
            logger.info("Testing tid = {}, formula = {}"
                         .format(unmatched[i][0].entry_id, test_host.formula))
            test_formula = test_host.composition.anonymized_formula \
                if anonymized else test_host.composition.reduced_formula
            logger.info("Test host = {}".format(test_formula))
            if test_formula == ref_formula:
                fitter = StructureFitter(ref_host, test_host,
                                         anonymized=anonymized,
                                         symmetry_tol=symmetry_tol,
                                         fitting_accuracy=fitting_accuracy)
                if fitter.fit_found:
                    logger.info("Fit found")
                    matches.append(unmatched[i])
        all_matches.append([m[0] for m in matches])
        unmatched = filter(lambda x: x not in matches, unmatched)
        logger.info("{} unmatched remaining".format(len(unmatched)))
    return all_matches
