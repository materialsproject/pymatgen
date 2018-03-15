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
from scipy import integrate

from utils import kb

from pymatgen.analysis.defects.corrections import ChargeCorrection, OtherCorrection
from pymatgen.analysis.defects.chemical_potentials import ChemPotAnalyzer, MPChemPotAnalyzer

class BulkEntry(object):
    """
    Similar to a SingleDefect object but for non-defective supercells
    Useful for parsing with DefectEntry, especially for storing information
        for formation energy corrections (electrostatic potentials etc.)
    """
    def __init__(self, structure, energy, supercell_size=(1, 1, 1), vbm=None, gap=None, bandedge_shifts=(0., 0.),
                 plnr_avg_esp=None, atomic_site_avg_esp=None, bulk_eigenvalues=None):
        self.structure = structure
        self.composition = structure.composition
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
            self.finished_charges[dfct.name] = [e.charge for e in self.dentries if e.name == t]
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

    @property
    def list_defect_concentrations(self, mu_elts, temp=300, ef=0.):
        """
        Give list of all concentrations at specified efermi in the DefectPhaseDiagram
        args:
            mu_elts = {Element: number} is dictionary of chemical potentials to provide formation energies for
            temp = temperature to produce concentrations from
            ef: (float) is fermi level relative to valence band maximum
                Default efermi = 0 = VBM energy
        returns:
            list of dictionaries of defect concentrations
        """
        concentrations = []
        for dfct in self.dentries:
            concentrations.append( {'conc':dfct.get_defect_concentration(mu_elts, temp=temp, ef=ef),
                                        'name':dfct.name, 'charge': dfct.charge})

        return concentrations

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


class IntrinsicCarrier(object):
    """
    A class for producing free carrier concentrations from a pymatgen DOS object
    """
    def __init__(self, dos, exp_gap=None):
        self.dos = dos
        self.exp_gap = exp_gap

    def get_n_density(self, ef, T, ref='VBM', unitcell=True):
        """
        Obtain the free electron concentration as a function of Fermi level and
        temperature using Fermi-Dirac statistics and conduction band DOS
        Args:
            ef: Fermi energy
            dos: Dos object of pymatgen
            ref: Options 'VBM' or 'CBM'
                The second option is useful when ef > DFT gap but <
                experimental gap. In that case, ef close to CBM can be
                given as -ve value w.r.t. CBM.
            unitcell: Output units in terms of unitcell or not,
                    if True output is #/unitcell, if False output is #/cm^3
        Returns:
            Electron density
        """
        dos_gap = self.dos.get_gap()
        cbm, vbm = self.dos.get_cbm_vbm()
        gap = self.exp_gap if self.exp_gap else dos_gap

        if ref == 'CBM':
            ef += cbm
        else:
            ef += vbm + dos_gap - gap

        energies = self.dos.energies
        densities = self.dos.get_densities()
        if energies[-1] - ef < 3.0:
            print ("The upper limit of energy is within 3 eV of Fermi level. "
                   "Check for the convergence of electron concentration")
        i = np.searchsorted(energies, cbm)
        fd_stat = 1./(1 + np.exp((energies[i:] - ef) / (kb*T)))
        y = fd_stat * densities[i:]
        den = integrate.trapz(y, energies[i:])
        if not unitcell:
            den *= (1e24 / self.dos.structure.volume)
        return den

    def get_p_density(self, ef, T, ref='VBM', unitcell=True):
        """
        Obtain the hole concentration as a function of Fermi level and
        temperature using Fermi-Dirac statistics and conduction band DOS
        Args:
            ef: Fermi energy
            dos: Dos object of pymatgen
            ref: Options 'VBM' or 'CBM'
                    The second option is useful when ef > DFT gap
                    but < experimental gap
            unitcell: Output units in terms of unitcell or not,
                    if True output is #/unitcell, if False output is #/cm^3
        Returns:
            Hole density
        """
        dos_gap = self.dos.get_gap()
        cbm, vbm = self.dos.get_cbm_vbm()
        gap = self.exp_gap if self.exp_gap else dos_gap

        if ref == 'CBM':
            ef += cbm + gap - dos_gap
        else:
            ef += vbm

        energies = self.dos.energies
        densities = self.dos.get_densities()
        if ef - energies[0] < 3.0:
            print ("The lower limit of energy is within 3 eV of Fermi level. "
                   "Check for the convergence of hole concentration")
        i = np.searchsorted(energies, vbm) + 1
        fd_stat = 1./(1 + np.exp((ef - energies[:i]) / (kb*T)))
        y = fd_stat*densities[:i]
        den = integrate.trapz(y, energies[:i])
        if not unitcell:
            den *= (1e24 / self.dos.structure.volume)
        return  den


class GrandCanonicalDefectPhaseDiagram(DefectPhaseDiagram):
    """
    Adapts the DefectPhaseDiagram class to a grand canonical approach
    This allows for self consistently determining fermi levels,
    defect concentrations, and free carrier concentrations as a function
    of temperature.

    Args:
        mu_elts: {Element: float} is dictionary of chemical potentials, determined from a phase diagram
        T: Temperature (Kelvin)
        kwargs: Arguments to DefectsAnalyzer as keyword pair
    """
    def __init__(self, mu_elts, T=298.15,  **kwargs):
        super(self.__class__, self).__init__(**kwargs)
        self.mu_elts = mu_elts
        self.T = T
        self._ef = None
        self._all_possible_mu_elts = None

    @staticmethod
    def generate_gcdpd_from_pda(dentries, pda, T=298.15):
        """
        A static method for instantiating a GrandCanonicalDefectPhaseDiagram from a DefectEntry list
        and a PhaseDiagramAnalyzer object from MaterialsProject
        """
        #use bulk_entry object from first dentry to find the composition
        #    that one requires chemical potentials from...
        CPA = ChemPotAnalyzer(bulk_ce=dentries[0].bulk)
        all_possible_mu_elts = CPA.get_chempots_from_pda( pda)
        print('Generated following possible regions for defining '
              'chemical potentials: ',all_possible_mu_elts.keys())
        mu_region, mu_elts = all_possible_mu_elts.items()[0]
        print('Proceeding with chemical potentials from ',mu_region,'\n',mu_elts)

        gcdpd = GrandCanonicalDefectPhaseDiagram(mu_elts, T=T, dentries=dentries)
        gcdpd.all_possible_mu_elts = all_possible_mu_elts

        return gcdpd

    @staticmethod
    def generate_gcdpd_from_MP(dentries, mpid=None, composition=None, T=298.15, mapi_key=None):
        """
        A static method for instantiating a GrandCanonicalDefectPhaseDiagram
        from either a mpid string or a pymatgen Composition object
            [If both are provided, defaults to use of mpid]
            [If neither are provided, uses bulk_entry of first DefectEntry to generate the composition objectt]

        Uses the MaterialsProject database to create a phase diagram
        """
        all_species = set([elt  for dfct in dentries  for elt in dfct._structure.composition.elements])
        bulk_species = set(dentries[0].bulk.composition.elements)
        sub_species = all_species - bulk_species

        if mpid:
            MPcpa = MPChemPotAnalyzer(sub_species=sub_species, mpid=mpid, mapi_key=mapi_key)
            all_possible_mu_elts = MPcpa.analyze_GGA_chempots()
        else:
            if not composition:
                composition = dentries[0].bulk.composition
            MPcpa = MPChemPotAnalyzer(sub_species=sub_species, mapi_key=mapi_key)
            all_possible_mu_elts = MPcpa.get_chempots_from_composition(composition)

        print('Generated following possible regions for defining '
              'chemical potentials: ',all_possible_mu_elts.keys())
        mu_region, mu_elts = all_possible_mu_elts.items()[0]
        print('Proceeding with chemical potentials from ',mu_region,'\n',mu_elts)

        gcdpd = GrandCanonicalDefectPhaseDiagram(mu_elts, T=T, dentries=dentries)
        gcdpd.all_possible_mu_elts = all_possible_mu_elts

        return gcdpd


    def set_T(self, T):
        self.T = T

    @property
    def fermi_energy(self):
        return self._ef

    def _get_total_q(self, ef, gap, bulk_dos):
        qd_tot = 0
        for d in self.list_defect_concentrations(self.mu_elts, temp=self.T, ef=ef):
            qd_tot += d['charge'] * d['conc']

        q_h_cont = IntrinsicCarrier(bulk_dos, exp_gap=gap)
        ef_ref_cbm = ef - gap
        nden = q_h_cont.get_n_density(ef_ref_cbm, self.T, ref='CBM', unitcell=False)
        pden = q_h_cont.get_p_density(ef, self.T, unitcell=False)
        qd_tot += pden - nden
        return qd_tot

    def solve_for_fermi_energy(self, bulk_dos):
        """
        Solve for the Fermi energy self-consistently as a function of T
        and p_O2
        Observations are Defect concentrations, electron and hole conc
        Args:
            bulk_dos: bulk system dos (pymatgen Dos object)
            gap: Can be used to specify experimental gap.
                Will be useful if the self consistent Fermi level
                is > DFT gap
        Returns:
            Fermi energy
        """
        from scipy.optimize import bisect
        self._ef = bisect(lambda e: self._get_total_q(e, self.band_gap, bulk_dos), -1., self.band_gap+1.)

    def solve_non_eq_fermilevel(self, lowT, highT, bulk_dos, show_carriers=True, hold_htemp_ef = True):
        """
        Solve for the Fermi energy at a low temperature value, but with defect concentrations frozen in
            at some higher temperature value This means that only the fermi level and free charge carriers concentrations
            will be solved self consistently (simulating a 'frozen-in' defect concentration approach)

        Args:
            lowT: the temperature at which you want to solve for fermi level and free carrier concentrations
            highT: the fixed temperature that sets the defect concentrations at a higher temperature
            bulk_dos: bulk system dos (pymatgen Dos object)
            show_carriers: Dictates the outputs of this function
                Set to True if you want the new carrier concentrations in addition to the fermi levels
            hold_htemp_ef: Specifies whether the fermi level will be fixed for the defect concentrations at high T
                Set to False if you want the fermi level of the high temperature defects to be the same as the fermi level
                for the low temperature charge carriers (see below note on implementation for explanation of this flag)
        Returns:
            Fermi energy

        A NOTE ON IMPLEMENTATION: while testing this function, it became apparent that negative defect formation energies
            can sometimes cause stochastic results for the fermi level as a result of the defect concentrations not
            being allowed to change within the self consistent result. To circumvent this, this function has a flag ('hold_htemp_ef')
            to allow for the fermi level of high temperature defects to vary within the self consistent approach.
            While this is no longer the fully physical 'frozen-in' approach, it helps remove the stocastic results
            that result from numerical issues

        """
        from scipy.optimize import bisect
        q_h_cont = IntrinsicCarrier(bulk_dos, exp_gap=self.band_gap)

        if hold_htemp_ef:
            self.solve_for_fermi_energy(bulk_dos, gap=self.band_gap)
            defectlist = self.get_defects_concentration(highT, self._ef, unitcell=False)

        def noneq_total_q(ef):
            qd_tot = 0
            if hold_htemp_ef:
                for d in defectlist:
                    qd_tot += d['charge'] * d['conc']
            else:
                for d in self.get_defects_concentration(highT, ef, unitcell=False):
                    qd_tot += d['charge'] * d['conc']
            ef_ref_cbm = ef - self.band_gap
            nden = q_h_cont.get_n_density(ef_ref_cbm, lowT, ref='CBM', unitcell=False)
            pden = q_h_cont.get_p_density(ef, lowT, unitcell=False)
            qd_tot += pden - nden
            return qd_tot

        non_eq_ef = bisect(lambda e: noneq_total_q(e), -1., self.band_gap+1.)
        if show_carriers:
            ef_ref_cbm = non_eq_ef - self.band_gap
            nden = q_h_cont.get_n_density(ef_ref_cbm, lowT, ref='CBM', unitcell=False)
            pden = q_h_cont.get_p_density(non_eq_ef, lowT, unitcell=False)
            return non_eq_ef, nden, pden
        else:
            return non_eq_ef