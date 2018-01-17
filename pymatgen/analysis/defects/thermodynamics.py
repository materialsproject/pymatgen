#!/usr/bin/env python

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Shyam Dwaraknath"
__email__ = "dbroberg@berkeley.edu, shyamd@lbl.gov"
__status__ = "Development"
__date__ = "January 11, 2018"

from math import exp
import numpy as np
from monty.json import MSONable

from pymatgen.analysis.defects.corrections import ChargeCorrection, OtherCorrection

class BulkEntry(object):
    """
    Similar to a SingleDefect object but for non-defective supercells
    Useful for parsing with DefectEntry, especially for storing information
        for formation energy corrections (electrostatic potentials etc.)
    """
    def __init__(self, structure, energy, supercell_size=(1, 1, 1), vbm=None, gap=None, bandedge_shifts=(0., 0.),
                 plnr_avg_esp=None, atomic_site_avg_esp=None, bulk_eigenvalues=None):
        self.structure = structure
        self.energy = energy
        self.supercell_size = supercell_size
        self.vbm = vbm
        self.band_gap = gap
        self.bandedge_shifts = bandedge_shifts
        self.plnr_avg_esp = plnr_avg_esp
        self.atomic_site_avg_esp = atomic_site_avg_esp
        self.bulk_eigenvalues = bulk_eigenvalues


class DefectEntry(object):
    """
    [Note very similar in form to ComputedDefect + DefectsAnalyzer objects of pycdt.core.defects_analyzer]

    Holds all the info concerning a defect computation
    Has ability to do formation energy analysis given chemical potential information
    """
    def __init__(self, defect, energy, bulk_entry, corrections = [], name=None):
        """
        Args:
            defect:
                A SingleDefect object from pymatgen.analysis.defects.point_defects
            energy (float): Energy of the defect entry. Usually the final calculated
                energy from VASP or other electronic structure codes.
            bulk_entry:
                A BulkEntry object
            corrections:
                List of Correction classes (from pymatgen.analysis.defects.corrections)
                which correct energy due to charge (e.g. Freysoldt or Kumagai)
                or other factors (e.g. Shallow level shifts)
            name:
                The name of the defect
        """

        self.defect = defect
        self.site = defect.defectsite
        self.multiplicity = defect.multiplicity
        self.supercell_size = defect.supercell_size
        self.charge = defect.charge
        self.energy = energy
        self.name = name
        if self.name:
            self.full_name = self.name + "_" + str(defect.charge)
        else:
            self.full_name = "defect_" + str(defect.charge)
        self.bulk = bulk_entry
        self.e_vbm = bulk_entry.vbm
        self.band_gap = bulk_entry.band_gap

        self.corrections = corrections # Can be added to after initialization

        self.e_fermi = None
        self.formation_energy = None

    def add_charge_correction(self, correction):
        """
        Manually change/add charge correction for defect
        Args:
            correction (float):
                New correction to be applied for defect
        """
        #first check if other charge correction already exists; delete all that exist
        c_corrindexlist = []
        for cindex, c in enumerate(self.corrections):
            if issubclass(type(c), ChargeCorrection):
                print('ChargeCorrection:',type(c),' already exists...will delete.')
                c_corrindexlist.append(cindex)

        for delcindex in sorted(c_corrindexlist, reverse=True):
            del c.corrections[delcindex]

        new_cc = ChargeCorrection(es_corr=correction) #note this does not contain information of correction...
        print('Adding Charge Correction: ',new_cc)
        c.corrections.append(new_cc)

    def add_other_correction(self, correction):
        """
        Manually change/add the other correction for defect
        Args:
            correction (float):
                New correction to be applied for defect
        """
        for cindex, c in enumerate(self.corrections):
            if issubclass(type(c), OtherCorrection):
                print('Note OtherCorrection exists:',type(c),'=',c.correction)

        new_oc = OtherCorrection(corr=correction) #note this does not contain information of correction...
        print('Adding Additional Other Correction: ',new_oc)
        c.corrections.append(new_oc)

    def get_formation_energy(self, mu_elts, ef=0.0):
        """
        compute the formation energy for a given dict of chemical potentials and fermi level
        """
        #TODO add flag that prohibits calculation of formens between different sized supercells
        mu_needed_coeffs = {}
        for elt in self.defect._structure.composition.elements:
            el_def_comp = self.defect.structure.composition[elt]
            el_blk_comp = self.bulk.structure.composition[elt]
            mu_needed_coeffs[elt] = el_blk_comp - el_def_comp

        sum_mus = 0.0
        for elt in mu_needed_coeffs:
            el = elt.symbol
            sum_mus += mu_needed_coeffs[elt] * mu_elts[el]

        total_corrections = 0.
        for c in self.corrections:
            total_corrections += c.correction

        self.e_fermi = ef
        self.formation_energy = self.energy - self.bulk.energy + \
                        sum_mus + self.charge*(self.e_vbm + ef) + \
                        total_corrections

    def add_other_correct_bg_simple(self, vbm_correct, cbm_correct):
        #TODO: note that this shift is available in the bulk entry that has been created
        # TODO: come up with better way to use this information for the corrections...
        """
        correct the band gap in the analyzer.
        We assume the defects level remain the same when moving the
        band edges
        Args:
            vbm_correct:
                The correction on the vbm as a positive number. e.g.,
                if the VBM goes 0.1 eV down vbm_correct=0.1
            cbm_correct:
                The correction on the cbm as a positive number. e.g.,
                if the CBM goes 0.1 eV up cbm_correct=0.1
        """
        self.band_gap = self.band_gap + cbm_correct + vbm_correct
        self.e_vbm = self._e_vbm - vbm_correct

    def get_defect_concentration(self, mu_elts, temp=300, ef=0.0):
        """
        Get the defect concentration for a temperature and Fermi level.
        Args:
            temp:
                the temperature in K
            Ef:
                the fermi level in eV (with respect to the VBM)
        Returns:
            defects concentration in m-3
        """
        struct = self.bulk.structure
        cell_multiplier = np.prod(self.supercell_size)
        n = self.multiplicity * cell_multiplier * 1e30 / struct.volume
        conc = n*exp( -self.get_formation_energy(mu_elts, ef=ef)/(kb*temp))

        return conc

class DefectPhaseDiagram(MSONable):
    """
    This is similar to a PhaseDiagram object in pymatgen, but has ability to do quick analysis of defect formation energies
    when fed DefectEntry objects

    uses many of the capabilities from PyCDT's DefectsAnalyzer class...

    This class is able to get:
        a) stability of charge states for a given defect,
        b) list of all formation ens
        c) transition levels in the gap
        d)

    Args:
        dentries ([DefectEntry]): A list of DefectEntry objects
    """
    def __init__(self, dentries):
        self.dentries = dentries
        self.band_gap = dentries[0].band_gap #TODO: run a check that all entries have same bandgaps and vbm values
        self._set_dpd_nochempots()

    def _set_dpd_nochempots(self):
        """
        Set the DefectPhaseDiagram attributes (which dont depend on chemical potential)
        """
        self.stable_charges = {} #keys are defect names, items are list of charge states that are stable
        self.finished_charges = {} #keys are defect names, items are list of charge states that are included in the phase diagram
        self.transition_levels = {} #keys are defect names, items are list of [fermi level for transition, previous q, next q] sets
        xlim = (-0.1, self.band_gap+.1)
        nb_steps = 10000
        x = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nb_steps)
        list_elts = [elt  for dfct in self.dentries  for elt in dfct._structure.composition.elements]
        no_chempots = {elt: 0. for elt in set(list_elts)} #zerod chemical potentials for calculating stable defects
        for t in self._get_all_defect_types():
            trans_level = []
            chg_type = []
            prev_min_q, cur_min_q = None, None
            for x_step in x:
                miny = 10000
                for dfct in self.dentries:
                    if dfct.name == t:
                        val = dfct.get_formation_energy(no_chempots, ef=x_step)
                        if val < miny:
                            miny = val
                            cur_min_q = dfct.charge

                if prev_min_q is not None:
                    if cur_min_q != prev_min_q:
                        trans_level.append((x_step, prev_min_q, cur_min_q))
                    if cur_min_q not in chg_type:
                        chg_type.append(cur_min_q)
                prev_min_q = cur_min_q

            self.stable_charges[dfct.name] = chg_type[:]
            self.finished_charges[dfct.name] = [e.charge for e in self.da._defects if e.name == t]
            self.transition_levels[dfct.name] = trans_level[:]

    def add_defect_entry(self, defect):
        """
        add a DefectEntry object to the entries list
        Args:
            defect:
                a DefectEntry object
        """
        self.dentries.append(defect)
        self._set_pd_nochempots()

    @property
    def all_defect_types(self):
        """
        List types of defects existing in the DefectPhaseDiagram
        """
        return self._get_all_defect_types()

    @property
    def all_defect_entries(self):
        """
        List all defect entries existing in the DefectPhaseDiagram
        """
        #TODO: return fulldefect type based on SingleDefect types, not the optional label
        return [e.full_name for e in self.dentries]

    @property
    def all_stable_entries(self):
        """
        List all stable entries (defect+charge) in the DefectPhaseDiagram
        """
        stable_entries = []
        for t, stabset in self.stable_charges.items():
            for c in stabset:
                stable_entries.append(t + "_" + str(c))
        return stable_entries

    @property
    def all_unstable_entries(self):
        """
        List all unstable entries (defect+charge) in the DefectPhaseDiagram
        """
        return [e for e in self.all_defect_entries if e not in self.all_stable_entries]

    @property
    def list_formation_energies(self, mu_elts, ef=0.):
        """
        Give list of all formation energies at specified efermi in the DefectPhaseDiagram
        args:
            mu_elts = {Element: number} is dictionary of chemical potentials to provide formation energies for
            ef: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy
        returns:
            list in format [ formation_energy, defect full-name ]
        """
        formation_energies = []
        for dfct in self.dentries:
            formation_energies.append( [dfct.get_formation_energy(mu_elts, ef=ef), dfct.full_name])

        return formation_energies

    def _get_all_defect_types(self):
        """
        return defect names from entry list. Note this works off of DefectEntry name field
        """
        #TODO: return defect type based on DefectEntry types, not the optional label
        to_return = []
        for d in self.dentries:
            if d.name not in to_return: to_return.append(d.name)
        return to_return

    def suggest_charges(self):
        """
        Based on entries, suggested follow up charge states to run
        (to make sure possibilities for charge states have been exhausted)
        """
        fullrecommendset = {}
        for t in self._get_all_defect_types():
            print('Consider recommendations for ',t)
            reccomendset = []
            allchgs = self.finished_charges[t]
            stablechgs = self.stable_charges[t]
            for followup_chg in range(min(stablechgs)-1,  max(stablechgs)+2):
                if followup_chg in allchgs:
                    continue
                else:
                    recflag = True
                    for tl in self.transition_levels[t]: #tl has list of [fermilev for trans, prev q, next q]
                        if tl[0] <= 0.1: #check if t.l. is within 0.1 eV of the VBM
                            morepos = int(tl[1])
                            if  (followup_chg > morepos):
                                print('Wont recommend:',followup_chg,'Because of this trans lev:', tl[1],'/',
                                      tl[2],' at ',tl[0])
                                recflag = False
                        if tl[0] >= (self.band_gap - 0.1): #check if t.l. is within 0.1 eV of CBM
                            moreneg = int(tl[2])
                            if  (followup_chg < moreneg):
                                print('Wont recommend:',followup_chg,'Because of this trans lev:', tl[1],'/',
                                      tl[2],' at ',tl[0], '(gap = ', self.band_gap,'eV)')
                                recflag = False
                    if recflag:
                        reccomendset.append(followup_chg)
            if len(reccomendset):
                print('charges recommending:',reccomendset)
            fullrecommendset[t] = reccomendset[:]

        return fullrecommendset