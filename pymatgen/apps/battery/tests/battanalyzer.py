# coding: utf-8

from __future__ import division, unicode_literals
from collections import defaultdict
from pymatgen.core.physical_constants import ELECTRON_VOLT, AVOGADROS_CONST

__author__ = "Anubhav Jain"
__copyright__ = "Copyright 2011, The Materials Project"
__credits__ = ["Shyue Ping Ong", "Geoffroy Hautier"]
__version__ = "1.0"
__maintainer__ = "Anubhav Jain"
__email__ = "ajain@lbl.gov"
__date__ = "Sep 20, 2011"

import math

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Composition


EV_PER_ATOM_TO_J_PER_MOL = ELECTRON_VOLT * AVOGADROS_CONST
ELECTRON_TO_AMPERE_HOURS = EV_PER_ATOM_TO_J_PER_MOL / 3600


class BatteryAnalyzer():
    """
    A suite of methods for starting with an oxidized structure and determining its potential as a battery
    """

    def __init__(self, struc_oxid, cation='Li'):
        """
        Pass in a structure for analysis
        Arguments:
            struc_oxid - a structure object. oxidation states *must* be assigned for this structure, else the analysis
            methods will NOT work! Disordered structures should be OK
            cation - a String symbol or Element for the cation. It must be positively charged, but can be 1+/2+/3+ etc.
        """
        self.struc_oxid = struc_oxid
        if (not isinstance(cation, Element)):
            self.cation = Element(cation)

        self.cation_charge = self.cation.max_oxidation_state

    @property
    def max_cation_removal(self):
        """
        Given a structure *with oxidation states (required)*, this returns the maximum number of cation A that can be removed
        while maintaining charge-balance (default = 'Li')
        Note that only cation removal is considered; insertion is considered in other methods
        The number of Li is not normalized in any way and will depend on cell size (this is an 'extrinsic' function!)
        """
        #how much 'spare charge' is left in the redox metals for oxidation?
        oxid_pot = sum([(Element(spec.symbol).max_oxidation_state - spec.oxi_state) * self.comp[spec] for spec in self.comp if is_redox_active_intercalation(Element(spec.symbol))])

        oxid_limit = oxid_pot / self.cation_charge

        #the number of A that exist in the structure for removal
        num_cation = self.comp[Specie(self.cation.symbol, self.cation_charge)]

        return min(oxid_limit, num_cation)


    @property
    def max_cation_insertion(self):
        """
        Given a structure *with oxidation states (required)*, this returns the maximum number of cation A that can be inserted
        while maintaining charge-balance (default = 'Li')
        Note that only cation insertion is considered; removal is considered in other methods
        No consideration is given to whether there (geometrically speaking) are Li sites to actually accommodate the extra Li
        The number of Li is not normalized in any way and will depend on cell size (this is an 'extrinsic' function!)
        """
        #how much 'spare charge' is left in the redox metals for reduction?
        lowest_oxid = defaultdict(lambda: 2, {'Cu':1})  # only Cu can go down to 1+
        oxid_pot = sum([(spec.oxi_state - min(e for e in Element(spec.symbol).oxidation_states if e>=lowest_oxid[spec.symbol])) * self.comp[spec] for spec in self.comp if is_redox_active_intercalation(Element(spec.symbol))])

        return oxid_pot / self.cation_charge


    def _get_max_cap_ah(self, remove, insert):
        """
        Give max capacity in mAh for inserting and removing a charged cation
        This method does not normalize the capacity and meant to be a helper method to prevent confusion
        """
        num_cations = 0
        if (remove):
            num_cations += self.max_cation_removal
        if (insert):
            num_cations += self.max_cation_insertion

        ampHours = num_cations * self.cation_charge * ELECTRON_TO_AMPERE_HOURS;
        return ampHours


    def get_max_capgrav(self, remove=True, insert=True):
        """
        Give max capacity in mAh/g for inserting and removing a charged cation
        Note that the weight is normalized to the most lithiated state:
            thus removal of 1 Li from LiFePO4 gives the same capacity as insertion of 1 Li into FePO4
        """
        weight = self.comp.weight
        if (insert):
            weight += self.max_cation_insertion * self.cation.atomic_mass
        return self._get_max_cap_ah(remove, insert) / (weight / 1000);


    def get_max_capvol(self, remove=True, insert=True, volume=None):
        """
        Give max capacity in mAh/cc for inserting and removing a charged cation
        Unfortunately the cell volume upon insertion/deinsertion cannot be known in advance:
            thus removal of 1 Li from LiFePO4 gives a different capacity than insertion of 1 Li into FePO4
            in one case the cell volume of LiFePO4 is used, in another the cell volume of FePO4 is used
        You can pass in an optional volume argument to normalize by in order to bypass this issue...
        """

        vol = volume if volume else self.struc_oxid.volume
        return self._get_max_cap_ah(remove, insert) * 1000 * 1E24 / (vol * AVOGADROS_CONST)

    def get_delith_int_oxid(self):
        """
        Returns a set of delithiation steps, e.g. set([1.0 2.0 4.0]) etc. where each element
        is a number of Li to remove from the raw, full cell (not normalized in any way)
            Thus, if you double the unit cell, your answers will be twice as large!
        The delithiation steps are chosen to be at integer oxidation states of the redox metals.
        If multiple redox metals are present, all combinations of reduction/oxidation are tested.
        Note that having more than 3 redox metals will likely slow down the algorithm; a more efficient
        implementation might solve this.

        Examples:
            LiFePO4 will return [1.0]
            Li4Fe3Mn1(PO4)4 will return [1.0, 2.0, 3.0, 4.0])
            Li6V4(PO4)6 will return [4.0, 6.0])  *note that this example is not normalized*
        """

        # the elements that can possibly be oxidized
        oxid_els = [Element(spec.symbol) for spec in self.comp if is_redox_active_intercalation(spec)]

        numa = set()
        for oxid_el in oxid_els:
            numa = numa.union(self._get_delith_intoxi(self.comp.copy(), oxid_el, oxid_els, numa))

        #convert from num a in structure to num a removed
        num_cation = self.comp[Specie(self.cation.symbol, self.cation_charge)]
        return set([num_cation - a for a in numa])

    def _get_delith_intoxi(self, spec_amts_oxi, oxid_el, oxid_els, numa):
        """
        This is a helper method for get_delith_int_oxid!
        Arguments:
            spec_amts_oxi - a dict of species to their amounts in the structure
            oxid_el - the element to oxidize
            oxid_els - the full list of elements that might be oxidized
            numa - a running set of numbers of A cation at integer oxidation steps, this is the recursive property to return
        Returns:
            a set of numbers A; steps for for oxidizing oxid_el first, then the other oxid_els in this list
        """

        #If Mn is the oxid_el, we have a mixture of Mn2+, Mn3+, determine the minimum oxidation state for Mn
        #this is the state we want to oxidize!
        oxid_old = min([spec.oxi_state for spec in spec_amts_oxi if spec.symbol == oxid_el.symbol])
        oxid_new = math.floor(oxid_old + 1)
        #if this is not a valid solution, break out of here and don't add anything to the list
        if (oxid_new > oxid_el.max_oxidation_state):
            return numa

        #update the spec_amts_oxi map to reflect that the oxidation took place
        spec_old = Specie(oxid_el.symbol, oxid_old)
        spec_new = Specie(oxid_el.symbol, oxid_new)
        specamt = spec_amts_oxi[spec_old]
        spec_amts_oxi = {sp:amt for sp, amt in spec_amts_oxi.items() if sp != spec_old}
        spec_amts_oxi[spec_new] = specamt
        spec_amts_oxi = Composition(spec_amts_oxi)
        #determine the amount of cation A in the structure needed for charge balance and add it to the list
        oxi_noA = sum([spec.oxi_state * spec_amts_oxi[spec] for spec in spec_amts_oxi if spec.symbol not in self.cation.symbol])
        a = max(0, -oxi_noA / self.cation_charge)
        numa = numa.union({a})
        #recursively try the other oxidation states
        if (a == 0):
            return numa
        else:
            for oxid_el in oxid_els:
                numa = numa.union(self._get_delith_intoxi(spec_amts_oxi.copy(), oxid_el, oxid_els, numa))
            return numa


def is_redox_active_intercalation (element):
        """
        True if element is redox active and interesting for intercalation materials
        """

        ns = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Nb', 'Mo', 'W', 'Sb', 'Sn', 'Bi']
        return element.symbol in ns

