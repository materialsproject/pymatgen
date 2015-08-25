# coding: utf-8

from __future__ import division, unicode_literals

"""
The class ChargedDefectsStructures aims at facilitating
setup of charged defect calculations with VASP.
"""

__author__ = "Bharat Medasani, Nils E. R. Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Bharat Medasani"
__email__ = "mbkumar@gmail.com,n.zimmermann@tuhh.de,geoffroy@uclouvain.be"
__status__ = "Development"
__date__ = "August 12, 2015"

import copy
from monty.string import str2unicode
from pymatgen.core.structure import PeriodicSite
from pymatgen.core.periodic_table import Specie, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.defects.point_defects import Vacancy
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator


def get_optimized_sc_scale(inp_struct, final_site_no):
    """
    Helper function to compute an integer vector for supercell scaling.
    Specifically, the function determines the supercell that maximizes
    the shortest separation of any point in the supercell
    to any of its 26 closest images.
    Intended usage is such that the supercell should be
    obtained by replicating the inp_struct structure
    along the crystallographic a-axis
    biggest[0] times, along the b-axis biggest[1] times, and along
    the c-axis biggest[2] times.

    Args:
        inp_struct (Structure): input structure (unit cell).
        final_site_no (integer): maximum number of sites in resulting
            supercell after scaling inp_struct with biggest.

    Returns:
        biggest (3x1 integer array): unit-cell scaling vector to be
            applied on inp_struct to obtain maximally final_site_no
            in supercell.
    """
    if inp_struct is None:
        raise ValueError('empty structure!')
    if final_site_no is None or final_site_no <= 0:
        raise ValueError('the target number of sites within the supercell'
                         ' has to be larger than zero!')

    # target_site = reference site in the original unit cell.
    target_site = inp_struct.sites[0]
    dictio = {}

    # At most, replicate unit cell 6 times
    # along each of the three crystallographic axes.
    # Should be safe for the systems anticipated so-far.
    for k1 in range(1,6):
        for k2 in range(1,6):
            for k3 in range(1,6):
                struct = inp_struct.copy()
                struct.make_supercell([k1, k2, k3])
                if len(struct.sites) > final_site_no:
                    continue

                # Find index of reference site target_site in supercell
                # with a distance criterion and a threshold of 0.001
                # Angstrom.
                index = None
                for i in range(struct.num_sites):
                    s = struct._sites[i]
                    if s.distance_from_point(target_site.coords) < 0.001:
                        index = i
                if index is None:
                    raise RuntimeError('could not find reference site in' \
                                       ' supercell structure!')

                # Find closest periodic image of reference site in the
                # supercell structure.  distance > 0.00001 makes sure
                # that we don't pick the the site in the original supercell
                # for which a=0, b=0, c=0.
                min_dist = 1000.0
                for a in range(-1,1):
                    for b in range(-1,1):
                        for c in range(-1,1):
                            distance = struct.get_distance(
                                    index, index, (a,b,c))
                            if  distance < min_dist and distance > 0.00001:
                                min_dist = distance
                min_dist = round(min_dist, 3)
                if dictio.has_key(min_dist):
                    if dictio[min_dist]['num_sites'] > struct.num_sites:
                    # Shouldn't this read:
                    # if dictio[min_dist]['num_sites'] < struct.num_sites:
                    # to ultimately get the supercell with the largest
                    # value of the shortest separation between periodic
                    # images while simultaneously trying to get as close to
                    # final_site_no as possible?!
                        dictio[min_dist]['num_sites'] = struct.num_sites
                        dictio[min_dist]['supercell'] = [k1, k2, k3]
                else:
                    dictio[min_dist] = {}
                    dictio[min_dist]['num_sites'] = struct.num_sites
                    dictio[min_dist]['supercell'] = [k1, k2, k3]
    if len(dictio) == 0:
        raise RuntimeError('no entries in dictionary containing candidates' \
                           ' of supercell scaling vectors!')
    min_dist = -1.0
    biggest = None
    for c in dictio:
        if c > min_dist:
            biggest = dictio[c]['supercell']
            min_dist = c
    if biggest is None or min_dist < 0.0:
        raise RuntimeError('could not find any supercell scaling vector!')
    return biggest


class ChargedDefectsStructures(object):
    """
    A class to generate charged defect structures for use in first 
    principles supercell formalism.  So-far antisites and vacancies
    are supported.  Automatic generation of interstitial defects
    are yet to come.
    """

    def __init__(self, structure, min_max_oxi={}, substitutions={}, 
                 oxi_states={}, cellmax=128, interstitial_sites=[],
                 antisites_flag=True, standardized=False, 
                 charge_states='liberal'):
        """
        Args:
            structure (Structure):
                bulk structure, on the basis of which the defective cells
                will be created.
            min_max_oxi (dict):
                The minimum and maximum oxidation state of each element as
                a dictionary; for example, {"O": (-2,0)}.  If not given,
                the oxidation states of pymatgen are considered.
            substitutions (dict):
                The allowed substitutions of elements as a dictionary.
                If not given, intrinsic defects (e.g., anti-sites)
                are computed only.  If given, extrinsic
                defects are considered as well, but they have to be
                explicitly specified.  Example: {"Co": ["Zn","Mn"]} means
                that Co sites will be substituted by Mn and Zn.
            oxi_states (dict):
                The oxidation state of the elements in the compound (e.g.,
                {"Fe": 2, "O": -2}.  If not given, the oxidation state of
                each site is computed with the bond-valence sum.
                WARNING: the bond-valence sum method can fail for
                mixed-valence compounds.
            cellmax (int):
                Maximum number of atoms allowed in the supercell.
            interstitials_sites ([PeriodicSite]):
                A list of periodic sites in which each corresponds to
                an interstitial atom that is to be inserted into the bulk
                structure.  Attention: the underlying lattice has to be
                consistent with the lattice of the bulk structure.
            antisites_flag (bool):
                If set to False, no anti-site configurations will be
                generated.
            standardized (bool):
                switch to indicate whether (True) or not (False) to use
                the standardized (i.e., primitive) structure rather than
                the unchanged input structure to generate the final
                bulk supercell, from which the defect cells will be
                generated.
            charge_states (string):
                Options are 'liberal' and 'conservative'.  If liberal is
                selected, more charge states are computed.
        """

        # Start setting up and error-checking.
        self.defects = []
        if cellmax <= 1:
            raise ValueError('maximum number of atoms in bulk supercell used' \
                             ' for setting up defect supercells must be'
                             ' larger than 1!')
        if min_max_oxi is not None:
            for key, value in min_max_oxi.items():
                if value[0] > value[1]:
                    raise ValueError('incorrect order of minimum and maximum' \
                                     ' oxidation state in dictionary for' \
                                     ' specie \"'+key+'\"!')
        self.cellmax = cellmax
        self.charge_states = charge_states
        if charge_states != 'liberal' and charge_states != 'conservative':
            raise ValueError('unrecognized charge-state switch!')

        # Create desired unit cell.
        spa = SpacegroupAnalyzer(structure,symprec=1e-2)
        prim_struct = spa.get_primitive_standard_structure()
        if standardized:
            self.struct = prim_struct
        else:
            self.struct = structure
        struct_species = self.struct.types_of_specie

        # Setup and error-check the substitution list.
        self.substitutions = {}
        for key, val in substitutions.items():
            self.substitutions[str2unicode(key)] = val
            found_key = False
            for specie in struct_species:
                if specie.symbol == str2unicode(key):
                    found_key = True
                    break
            if not found_key:
                raise ValueError('could not find any specie in' \
                                 ' input structure that matches key \"'+ \
                                 str2unicode(key)+'\" in list of targeted' \
                                 ' substitutions!')

        # Setup and error-check the oxidation states.
        if not oxi_states:
            if len(struct_species) == 1:
                oxi_states = {self.struct.types_of_specie[0].symbol: 0}
            else:
                vir = ValenceIonicRadiusEvaluator(self.struct)
                oxi_states = vir.valences
        else:
            if len(oxi_states) != len(struct_species):
                raise ValueError('number of oxidation states does not match' \
                                 ' number of distinct species in structure!')
        self.oxi_states = {}
        for key,val in oxi_states.items():
            strip_key = ''.join([s for s in key if s.isalpha()])
            self.oxi_states[str2unicode(strip_key)] = val
        for specie in self.struct.types_of_specie:
            found_oxi_state = False
            for key in self.oxi_states.keys():
                if specie.symbol == key:
                    found_oxi_state = True
                    break
            if not found_oxi_state:
                raise ValueError('could not find any oxidation state for' \
                                 ' specie \"'+specie.symbol+'\"!')
        if oxi_states and len(self.struct.types_of_specie) == 1:
            if self.oxi_states[self.struct.types_of_specie[0].symbol] != 0:
                raise ValueError('prohibited to assign an oxidation state'
                                 ' to an elemental structure!')

        # Generate bulk supercell of desired maximum size from unit cell.
        conv_prim_rat = int(self.struct.num_sites/prim_struct.num_sites)
        sc_scale = get_optimized_sc_scale(self.struct, cellmax)
        self.defects = {}
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        self.defects['bulk'] = {
                'name': 'bulk',
                'supercell': {'size': sc_scale, 'structure': sc}}

        # Setup and error-check ranges of oxidation states, which will be
        # used for setting up all charge states of each individual defect
        # structure.  If not explicitly specified during constructor call,
        # the ranges are based on the minimum and maximum value of
        # common oxidation state as implemented in the Element class.
        if not min_max_oxi:
            min_max_oxi = {}
            for s in struct_species:
                if isinstance(s, Specie):
                    el = s.element
                elif isinstance(s, Element):
                    el = s
                else:
                    raise TypeError('expected an object of Element or' \
                                    ' Specie type in input structure!')
                max_oxi = max(el.common_oxidation_states)
                min_oxi = min(el.common_oxidation_states)
                min_max_oxi[str2unicode(el.symbol)] = (min_oxi, max_oxi)
            for s, subspecies in self.substitutions.items():
                for subspecie in subspecies:
                    el = Element(subspecie)
                    max_oxi = max(el.common_oxidation_states)
                    min_oxi = min(el.common_oxidation_states)
                    min_max_oxi[str2unicode(el.symbol)] = (min_oxi, max_oxi)
        self.min_max_oxi = min_max_oxi

        vacancies = []
        as_defs = []
        sub_defs = []

        vac = Vacancy(self.struct, {}, {})
        vac_scs = vac.make_supercells_with_defects(sc_scale)

        for i in range(vac.defectsite_count()):
            vac_site = vac.get_defectsite(i)
            site_mult = vac.get_defectsite_multiplicity(i)
            site_mult = int(site_mult/conv_prim_rat)
            vac_specie = vac_site.specie
            vac_symbol = vac_site.specie.symbol
            vac_sc = vac_scs[i+1]
            vac_sc_site = list(set(vac_scs[0].sites) - set(vac_sc.sites))[0]

            list_charges=[]
            vac_oxi_state = self.oxi_states[str2unicode(vac_symbol)]
            if vac_oxi_state < 0:
                min_oxi = min(vac_oxi_state, self.min_max_oxi[vac_symbol][0])
                max_oxi = 0
            elif vac_oxi_state > 0:
                min_oxi = 0
                max_oxi = max(vac_oxi_state, self.min_max_oxi[vac_symbol][1])
            else:
                min_oxi = self.min_max_oxi[vac_symbol][0]
                max_oxi = self.min_max_oxi[vac_symbol][1]
            for c in range(min_oxi, max_oxi+1):
                list_charges.append(-c)

            vacancies.append({
                'name': "vac_{}_{}".format(i+1, vac_symbol),
                'unique_site': vac_site,
                'bulk_supercell_site': vac_sc_site,
                'defect_type': 'vacancy',
                'site_specie': vac_symbol,
                'site_multiplicity': site_mult,
                'supercell': {'size': sc_scale,'structure': vac_sc},
                'charges': list_charges })

            # Antisite defects generation
            if antisites_flag:
                for as_specie in set(struct_species)-set([vac_specie]):
                    as_symbol = as_specie.symbol
                    as_sc = vac_sc.copy()
                    as_sc.append(as_symbol, vac_sc_site.frac_coords)
                    if vac_oxi_state > 0:
                        oxi_max = max(self.min_max_oxi[as_symbol][1],0)
                        oxi_min = 0
                    else:
                        oxi_max = 0
                        oxi_min = min(self.min_max_oxi[as_symbol][0],0)
                    if self.charge_states=='liberal' and oxi_min==oxi_max:
                        if oxi_min - vac_oxi_state > 0:
                            charges = list(range(oxi_min-vac_oxi_state+1))
                        else:
                            charges = list(range(oxi_min-vac_oxi_state-1,1))
                    else:
                        charges = [c - vac_oxi_state for c in range(
                            oxi_min, oxi_max+1)]

                    as_defs.append({
                        'name': "as_{}_{}_on_{}".format(
                            i+1, as_symbol, vac_symbol),
                        'unique_site': vac_site,
                        'bulk_supercell_site': vac_sc_site,
                        'defect_type': 'antisite',
                        'site_specie': vac_symbol,
                        'substitution_specie': as_symbol,
                        'site_multiplicity': site_mult,
                        'supercell': {'size': sc_scale,'structure': as_sc},
                        'charges': charges})

            # Substitutional defects generation
            if vac_symbol in self.substitutions:
                for subspecie_symbol in self.substitutions[vac_symbol]:
                    sub_sc = vac_sc.copy()
                    sub_sc.append(subspecie_symbol, vac_sc_site.frac_coords)
                    if vac_oxi_state > 0:
                        oxi_max = max(self.min_max_oxi[subspecie_symbol][1],0)
                        oxi_min = 0
                    else:
                        oxi_max = 0
                        oxi_min = min(self.min_max_oxi[subspecie_symbol][0],0)
                    if self.charge_states=='liberal' and oxi_min==oxi_max:
                        if oxi_min - vac_oxi_state > 0:
                            charges = list(range(oxi_min-vac_oxi_state+1))
                        else:
                            charges = list(range(oxi_min-vac_oxi_state-1,1))
                    else:
                        charges = [c - vac_oxi_state for c in range(
                            oxi_min, oxi_max+1)]

                    sub_defs.append({
                        'name': "sub_{}_{}_on_{}".format(
                            i+1, subspecie_symbol, vac_symbol),
                        'unique_site': vac_site,
                        'bulk_supercell_site': vac_sc_site,
                        'defect_type':'antisite',
                        'site_specie':vac_symbol,
                        'substitution_specie':subspecie_symbol,
                        'site_multiplicity':site_mult,
                        'supercell':{'size':sc_scale,'structure':sub_sc},
                        'charges':charges})

        self.defects['vacancies'] = vacancies 
        self.defects['substitutions'] = sub_defs
        self.defects['substitutions'] += as_defs

        #interstitials
        interstitials = []
        for elt in self.struct.composition.elements:
            count = 1
            for frac_coord in interstitial_sites:
                site = PeriodicSite(elt, frac_coord, structure.lattice)
                interstitials.append({
                    'name': elt.symbol+str(count)+"_inter",
                    'unique_site': site,
                    'supercell': {'size': s_size,
                        'structure': self.make_interstitial(site, sc_scale)},
                    'charges': [c for c in range(
                        min_max_oxi[elt][0], min_max_oxi[elt][1]+1)]})
                count = count+1
        self.defects['interstitials'] = interstitials


    def make_interstitial(self, target_site, sc_scale):
        """
        Generates a supercell that contains a single interstitial defect
            from the input bulk structure.
        Args:
            target_site (PeridocSite): site to be included in the bulk
                structure to obtain a single interstitial defect.
            sc_scale ([int, int, int]): scaling vector to obtain the
                desired supercell from the input bulk structure.
        Return:
            sc (Structure): bulk structure scaled by sc_scale and
                containing the interstitial defect stored in target_site.
        """
        sc = self.struct.copy()
        sc.make_supercell(sc_scale)
        sc.append(target_site.specie, target_site.frac_coords)
        
        return sc


    def get_defects_data(self):
        """
        Get all data (i.e., structures) for the targeted charged-defect
            calculations.
        Return:
            defects (dict): all structures necessary for the targeted
                charged-defect calculations.
        """
        defects = self.defects
        return defects
