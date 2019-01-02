# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module implements functions to perform various useful operations on
entries, such as grouping entries by structure.
"""



__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 24, 2012"

import logging
import json
import datetime
import collections

from monty.json import MontyEncoder, MontyDecoder

from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher, \
    SpeciesComparator

logger = logging.getLogger(__name__)


def _get_host(structure, species_to_remove):
    if species_to_remove:
        s = structure.copy()
        s.remove_species(species_to_remove)
        return s
    else:
        return structure


def _perform_grouping(args):
    (entries_json, hosts_json, ltol, stol, angle_tol,
     primitive_cell, scale, comparator, groups) = args

    entries = json.loads(entries_json, cls=MontyDecoder)
    hosts = json.loads(hosts_json, cls=MontyDecoder)
    unmatched = list(zip(entries, hosts))
    while len(unmatched) > 0:
        ref_host = unmatched[0][1]
        logger.info(
            "Reference tid = {}, formula = {}".format(unmatched[0][0].entry_id,
                                                      ref_host.formula)
        )
        ref_formula = ref_host.composition.reduced_formula
        logger.info("Reference host = {}".format(ref_formula))
        matches = [unmatched[0]]
        for i in range(1, len(unmatched)):
            test_host = unmatched[i][1]
            logger.info("Testing tid = {}, formula = {}"
                        .format(unmatched[i][0].entry_id, test_host.formula))
            test_formula = test_host.composition.reduced_formula
            logger.info("Test host = {}".format(test_formula))
            m = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol,
                                 primitive_cell=primitive_cell, scale=scale,
                                 comparator=comparator)
            if m.fit(ref_host, test_host):
                logger.info("Fit found")
                matches.append(unmatched[i])
        groups.append(json.dumps([m[0] for m in matches], cls=MontyEncoder))
        unmatched = list(filter(lambda x: x not in matches, unmatched))
        logger.info("{} unmatched remaining".format(len(unmatched)))


def group_entries_by_structure(entries, species_to_remove=None,
                               ltol=0.2, stol=.4, angle_tol=5,
                               primitive_cell=True, scale=True,
                               comparator=SpeciesComparator(),
                               ncpus=None):
    """
    Given a sequence of ComputedStructureEntries, use structure fitter to group
    them by structural similarity.

    Args:
        entries: Sequence of ComputedStructureEntries.
        species_to_remove: Sometimes you want to compare a host framework
            (e.g., in Li-ion battery analysis). This allows you to specify
            species to remove before structural comparison.
        ltol (float): Fractional length tolerance. Default is 0.2.
        stol (float): Site tolerance in Angstrom. Default is 0.4 Angstrom.
        angle_tol (float): Angle tolerance in degrees. Default is 5 degrees.
        primitive_cell (bool): If true: input structures will be reduced to
            primitive cells prior to matching. Defaults to True.
        scale: Input structures are scaled to equivalent volume if true;
            For exact matching, set to False.
        comparator: A comparator object implementing an equals method that
            declares equivalency of sites. Default is SpeciesComparator,
            which implies rigid species mapping.
        ncpus: Number of cpus to use. Use of multiple cpus can greatly improve
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
            symm_entries[comparator.get_structure_hash(host)].append((entry,
                                                                      host))
        import multiprocessing as mp
        logging.info("Using {} cpus".format(ncpus))
        manager = mp.Manager()
        groups = manager.list()
        p = mp.Pool(ncpus)
        #Parallel processing only supports Python primitives and not objects.
        p.map(_perform_grouping,
              [(json.dumps([e[0] for e in eh], cls=MontyEncoder),
                json.dumps([e[1] for e in eh], cls=MontyEncoder),
                ltol, stol, angle_tol, primitive_cell, scale,
                comparator, groups)
               for eh in symm_entries.values()])
    else:
        groups = []
        hosts = [host for entry, host in entries_host]
        _perform_grouping((json.dumps(entries, cls=MontyEncoder),
                           json.dumps(hosts, cls=MontyEncoder),
                           ltol, stol, angle_tol, primitive_cell, scale,
                           comparator, groups))
    entry_groups = []
    for g in groups:
        entry_groups.append(json.loads(g, cls=MontyDecoder))
    logging.info("Finished at {}".format(datetime.datetime.now()))
    logging.info("Took {}".format(datetime.datetime.now() - start))
    return entry_groups
