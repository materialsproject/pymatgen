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
import json
import datetime
import collections

from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.analysis.structure_fitter import StructureFitter
from pymatgen.serializers.json_coders import PMGJSONEncoder, PMGJSONDecoder

logger = logging.getLogger(__name__)


def _get_host(structure, species_to_remove):
    if species_to_remove:
        editor = StructureEditor(structure)
        editor.remove_species(species_to_remove)
        return editor.modified_structure
    else:
        return structure


def _perform_grouping(args):
    (entries_json, hosts_json, anonymized, symmetry_tol, fitting_accuracy,
     groups) = args
    entries = json.loads(entries_json, cls=PMGJSONDecoder)
    hosts = json.loads(hosts_json, cls=PMGJSONDecoder)
    unmatched = zip(entries, hosts)
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
        groups.append(json.dumps([m[0] for m in matches], cls=PMGJSONEncoder))
        unmatched = filter(lambda x: x not in matches, unmatched)
        logger.info("{} unmatched remaining".format(len(unmatched)))


def group_entries_by_structure(entries, species_to_remove=None,
                               anonymized=False, symmetry_tol=0,
                               fitting_accuracy=StructureFitter.FAST_FIT,
                               ncpus=None):
    """
    Given a sequence of ComputedStructureEntries, use structure fitter to group
    them by structural similarity.

    Args:
        entries:
            Sequence of ComputedStructureEntries.
        species_to_remove:
            Sometimes you want to compare a host framework (e.g., in Li-ion
            battery analysis). This allows you to specify species to remove
            before structural comparison.
        anonymized:
            Allow anonymized fitting.
        fitting_accuracy:
            The fitting accuracy to use for StructureFitter. Defaults to
            FAST_FIT.
        ncpus:
            Number of cpus to use. Use of multiple cpus can greatly improve
            fitting speed. Default of None means serial processing.

    Returns:
        Sequence of sequence of entries by structural similarity. e.g,
        [[ entry1, entry2], [entry3, entry4, entry5]]
    """

    start = datetime.datetime.now()
    logger.info("Started at {}".format(start))
    entries_host = [(entry, _get_host(entry.structure, species_to_remove))
                    for entry in entries]
    if ncpus:
        symm_entries = collections.defaultdict(list)
        for entry, host in entries_host:
            key = [host.composition.anonymized_formula
                   if anonymized
                   else host.composition.reduced_formula]
            if symmetry_tol:
                finder = SymmetryFinder(host, symmetry_tol)
                key.append(finder.get_spacegroup_number())
            symm_entries[tuple(key)].append((entry, host))
        import multiprocessing as mp
        logging.info("Using {} cpus".format(ncpus))
        manager = mp.Manager()
        groups = manager.list()
        p = mp.Pool(ncpus)
        #Parallel processing only supports Python primitives and not objects.
        p.map(_perform_grouping,
              [(json.dumps([e[0] for e in eh], cls=PMGJSONEncoder),
                json.dumps([e[1] for e in eh], cls=PMGJSONEncoder),
                anonymized,
                symmetry_tol, fitting_accuracy, groups)
               for eh in symm_entries.values()])
    else:
        groups = []
        hosts = [host for entry, host in entries_host]
        _perform_grouping((json.dumps(entries, cls=PMGJSONEncoder),
                           json.dumps(hosts, cls=PMGJSONEncoder),
                           anonymized,
                           symmetry_tol, fitting_accuracy, groups))
    entry_groups = []
    for g in groups:
        entry_groups.append(json.loads(g, cls=PMGJSONDecoder))
    logging.info("Finished at {}".format(datetime.datetime.now()))
    logging.info("Took {}".format(datetime.datetime.now() - start))
    return entry_groups
