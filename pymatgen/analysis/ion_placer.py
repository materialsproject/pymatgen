# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import logging

from monty.json import MSONable
from pymatgen import core
import numpy as np


__author__ = "Samuel Blau"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Samuel Blau"
__email__ = "samblau1@gmail.com"
__status__ = "Alpha"
__date__ = "11/6/18"


logger = logging.getLogger(__name__)


class IonPlacer(MSONable):

    def __init__(self, molecule, ion, stop_num):
        """
        Standard constructor for ion placement

        Args:
            molecule (Molecule)
            ion (str)
            stop_num (int)

        """
        self.mol = molecule
        self.vdwR = core.periodic_table.Specie(ion).van_der_waals_radius+0.0
        self.atR = core.periodic_table.Specie(ion).atomic_radius+0.0
        self.stop_num = stop_num
        self.accepted_points = []
        self._identify_reactive_sites()
        self._define_box()
        self._guess_until_stop()

    def _identify_reactive_sites(self):
        reactivity = []
        for site in self.mol:
            properties = site.properties
            if "charge" not in properties:
                raise KeyError("Each site in the molecule must have the charge property! Exiting...")
            if str(site.specie) == "H":
                reactivity.append(False)
            elif "spin" in properties:
                if float(properties["spin"]) > 0.05 or float(properties["charge"]) < 0.0:
                    reactivity.append(True)
                else:
                    reactivity.append(False)
            elif float(properties["charge"]) < 0.0:
                reactivity.append(True)
            else:
                reactivity.append(False)
        self.mol.add_site_property("reactivity",reactivity)

    def _define_box(self):
        self.min_vals = np.zeros(3)+100000000
        self.max_vals = np.zeros(3)-100000000
        for site in self.mol:
            if site.x < self.min_vals[0]:
                self.min_vals[0] = site.x
            elif site.x > self.max_vals[0]:
                self.max_vals[0] = site.x
            if site.y < self.min_vals[1]:
                self.min_vals[1] = site.y
            elif site.y > self.max_vals[1]:
                self.max_vals[1] = site.y
            if site.z < self.min_vals[2]:
                self.min_vals[2] = site.z
            elif site.z > self.max_vals[2]:
                self.max_vals[2] = site.z
        self.max_vals+=self.vdwR
        self.min_vals-=self.vdwR

    def _random_in_box(self):
        guess = np.random.random(3)
        for ii in range(3):
            guess[ii] = guess[ii]*(self.max_vals[ii]-self.min_vals[ii]) + self.min_vals[ii]
        return guess

    def _check_acceptance(self, point):
        reactive_in_range = 0
        for site in self.mol:
            dist = site.distance_from_point(point)
            if dist < self.atR:
                # Within atomic radius of an atom = reject
                return False
            if site.properties["reactivity"]:
                if dist < self.vdwR:
                    # In right range from a reactive site!
                    reactive_in_range += 1
            elif dist < self.vdwR*0.9:
                # Closer than vdw*0.9 of a nonreactive site = reject
                return False
        if reactive_in_range == 0:
            # No reactive sites in range = reject
            return False
        else:
            for good_pt in self.accepted_points:
                dist = np.linalg.norm(good_pt - point)
                if dist < self.atR and reactive_in_range == 1:
                    # Within an atomic radius of an accepted point
                    # and only one reactive site in range = reject
                    return False
                elif dist < self.atR*0.8:
                    # Within 0.8*atR of an accepted point and more
                    # than one reactive site in range = reject
                    return False
            return True

    def _guess_until_stop(self):
        num_reject_in_a_row = 0
        while num_reject_in_a_row < self.stop_num:
            guess = self._random_in_box()
            if self._check_acceptance(guess):
                self.accepted_points.append(guess)
                num_reject_in_a_row = 0
            else:
                num_reject_in_a_row +=1
            





        
