# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
"""
This module implements defect equivalents of the basic ComputedEntry objects, 
which is the basic entity that can be used to perform defect thermodynamic analyses
"""

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils E. R. Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__status__ = "Development"
__date__ = "Mar 15, 2018"

import abc
import six
import numpy as np

from monty.json import MSONable

logger = logging.getLogger(__name__)


class DefectEntry(MSONable):
    """
    An lightweight DefectEntry object containing key computed data
    for many defect analysis.
    """

    def __init__(self, defect, uncorrected_energy, vbm=0, corrections={}, parameters={}, entry_id=None):
        """
        Args:
            defect:
                A Defect object from pymatgen.analysis.defects.core
            uncorrected_energy (float): Energy of the defect entry. Usually the difference between
                the final calculated energy for the defect supercell - the perfect 
                supercell energy
            vbm: 

            corrections ([Correction]):
                List of Correction classes (from pymatgen.analysis.defects.corrections)
                which correct energy due to charge (e.g. Freysoldt or Kumagai)
                or other factors (e.g. Shallow level shifts)
            parameters (dict): An optional dict of calculation parameters and data to 
                use with correction schemes
            entry_id (obj): An id to uniquely identify this defect, can be any MSONable
                type
        """
        self.defect = defect
        self.uncorrected_energy = uncorrected_energy
        self.corrections = corrections
        self.entry_id = entry_id
        self.parameter = parameters

    @property
    def site(self):
        return self.defect.site

    @property
    def multiplicty(self):
        return defect.multiplicty

    @property
    def charge(self):
        return defect.charge

    @property
    def energy(self):
        """
        Returns the *corrected* energy of the entry
        """
        return self.uncorrected_energy + np.sum(self.correction.values())

    def formation_energy(self, chemical_potentials, fermi_level=0):
        """
        Computes the formation energy for a defect taking into account a given chemical potential and fermi_level
        """
        chempot_correction = sum([
            chem_pot * (self.defec.structure.composition[el] - self.defect.defect_composition[el])
            for el, chem_pot in chemical_potentials
        ])

        formation_energy = self.energy + chempot_correction

        if "vbm" in self.parameters:
            formation_energy += self.charge * (self.parameters["vbm"] + fermi_level)

        return formation_energy

    def __repr__(self):
        """
        Human readable string representation of this entry
        """
        output = [
            "DefectEntry {} - {}".format(self.entry_id, self.defect.name), "Energy = {:.4f}".format(self.energy),
            "Correction = {:.4f}".format(np.sum(self.correction.values())), "Parameters:"
        ]
        for k, v in self.parameters.items():
            output.append("\t{} = {}".format(k, v))
        return "\n".join(output)

    def __str__(self):
        return self.__repr__()