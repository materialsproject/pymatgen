# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module implements abstract base classes for post-processing entries.
Any class which modifies entries should inherit these classes.
"""

import six

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Oct 6, 2011"

import abc


class EntryPostProcessor(six.with_metaclass(abc.ABCMeta, object)):
    @abc.abstractmethod
    def process_entry(self, entry):
        """
        Process a single entry.

        Args:
            entry: An ComputedEntry object.

        Returns:
            An processed entry. None if entry is not compatible within the
            processing scheme.
        """
        return

    @abc.abstractmethod
    def process_entries(self, entries):
        """
        Process a sequence of entries.

        Args:
            entries: A sequence of ComputedEntries.

        Returns:
            An list of processed entries.  ComputedEntries in the original list
            which are not compatible with the processing scheme are excluded.
        """
        return

    @property
    @abc.abstractmethod
    def corrected_compound_formulas(self):
        """
        List of compound formulas that are corrected.
        """
        return
