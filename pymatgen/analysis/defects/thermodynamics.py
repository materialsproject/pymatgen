#!/usr/bin/env python

__author__ = "Danny Broberg, Shyam Dwaraknath, Bharat Medasani, Nils Zimmermann, Geoffroy Hautier"
__copyright__ = "Copyright 2014, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Danny Broberg, Shyam Dwaraknath"
__email__ = "dbroberg@berkeley.edu, shyamd@lbl.gov"
__status__ = "Development"
__date__ = "January 11, 2018"


from monty.json import MSONable


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
        self.charge_correction = correction

    def add_other_correction(self, correction):
        """
        Manually change/add the other correction for defect
        Args:
            correction (float):
                New correction to be applied for defect
        """
        self.other_correction = correction

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

        self.e_fermi = ef
        self.formation_energy = self.energy - self.bulk.energy + \
                        sum_mus + self.charge*(self.e_vbm + ef) + \
                        self.charge_correction + self.other_correction

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
        self.compute_form_en()


class DefectPhaseDiagram(MSONable):
    """
    This is similar to a PhaseDiagram object in pymatgen, but has ability to do quick analysis of defect formation energies
    when fed ComputedDefectEntries ( or, for now, fw_metadata...

    uses many of the capabilities from PyCDT's DefectsAnalyzer class...

    SHould be able to get:
        a) stability of charge states for a given defect,
        b) list of all formation ens
        c) [later on] fermi energies...

    Args:
        fw_metadata (list): fw_metadata in update_spec after DefectAnalysisFireTask is run
    """
    #TODO: make this allow uploads of DefectEntry objects...
    def __init__(self, fw_metadata):
        tmpkey = ''
        tmpind = 0
        while not tmpkey:
            if fw_metadata.keys()[tmpind] in ['entry_bulk','_files_prev']:
                tmpind+=1
            else:
                tmpkey = fw_metadata.keys()[tmpind]
        basedat = fw_metadata[tmpkey]['corr_details']['base']
        #load up the defects analyzer objects
        self.fullset_da = {}
        self.da = None #only need one da for looking at stability...
        for mnom, muset in basedat['chem_lims'].items():
            self.fullset_da[mnom] = DefectsAnalyzer( fw_metadata['entry_bulk'], basedat['vbm'], muset,
                                                     basedat['gga_bandgap'])
            for dnom, def_entry in fw_metadata.items():
                if dnom not in ['entry_bulk', '_files_prev']:
                    print dnom,  def_entry['cd'].charge_correction, def_entry['cd'].other_correction
                    self.fullset_da[mnom].add_computed_defect(copy.copy(def_entry['cd']) )
            if not self.da:
                self.da = copy.copy(self.fullset_da[mnom])

        #now get stable defect sets
        self.stable_charges = {} #keys are defect names, items are list of charge states that are stable
        self.finished_charges = {} #keys are defect names, items are list of charge states that are included in the phase diagram
        self.transition_levels = {} #keys are defect names, items are list of [fermi level for transition, previous q, next q] sets
        xlim = (-0.1, self.da._band_gap+.1)
        nb_steps = 10000
        x = np.arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/nb_steps)
        for t in self.da._get_all_defect_types():
            print('hey dan:',t, 'tote defs=',len(self.da._defects))
            trans_level = []
            chg_type = []
            prev_min_q, cur_min_q = None, None
            for x_step in x:
                miny = 10000
                for i, dfct in enumerate(self.da._defects):
                    if dfct.name == t:
                        val = self.da._formation_energies[i] + \
                                dfct.charge*x_step
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


    @property
    def all_defect_types(self):
        """
        List types of defects existing in the DefectPhaseDiagram
        """
        return self.da._get_all_defect_types()

    @property
    def all_defect_entries(self):
        """
        List all defect entries existing in the DefectPhaseDiagram
        """
        return [e.full_name for e in self.da.defects]

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

    # @property
    # def formation_energies(self, efermi=0.):
    #     """
    #     Give dictionaries of all formation energies at specified efermi in the DefectPhaseDiagram
    #     Default efermi = 0 = VBM energy TODO: need to specify growth condition?
    #     """
    #     return self.da.get_formation_energies(ef=efermi)


    def suggest_charges(self):
        """
        Based on entries, suggested follow up charge states to run
        (to make sure possibilities for charge states have been exhausted)
        """
        fullrecommendset = {}
        for t in self.da._get_all_defect_types():
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
                            # tmpstrchg = tl[2].split('/') #str(prev_min_q)+'/'+str(cur_min_q)
                            # morepos = int(tmpstrchg[0]) #more positive charge state
                            if  (followup_chg > morepos):
                                print('Wont recommend:',followup_chg,'Because of this trans lev:', tl[1],'/',
                                      tl[2],' at ',tl[0])
                                recflag = False
                        if tl[0] >= (self.da._band_gap - 0.1): #check if t.l. is within 0.1 eV of CBM
                            moreneg = int(tl[2])
                            # tmpstrchg = tl[2].split('/') #str(prev_min_q)+'/'+str(cur_min_q)
                            # moreneg = int(tmpstrchg[1]) #more negative charge state
                            if  (followup_chg < moreneg):
                                print('Wont recommend:',followup_chg,'Because of this trans lev:', tl[1],'/',
                                      tl[2],' at ',tl[0], '(gap = ', self.da._band_gap,'eV)')
                                recflag = False
                    if recflag:
                        reccomendset.append(followup_chg)
            if len(reccomendset):
                print('charges recommending:',reccomendset)
            fullrecommendset[t] = reccomendset[:]

        return fullrecommendset