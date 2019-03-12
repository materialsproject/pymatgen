# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest
import os
import json

from monty.serialization import loadfn

from pymatgen.core import Element
from pymatgen.analysis.defects.core import DefectEntry
from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram
from pymatgen.util.testing import PymatgenTest
from pymatgen.electronic_structure.dos import CompleteDos

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')

class DefectsThermodynamicsTest(PymatgenTest):
    def setUp(self):
        self.vbm_val = 2.6682
        self.gap = 2.
        self.entries = list(loadfn(os.path.join(os.path.dirname(__file__), "GaAs_test_defentries.json")).values())
        for entry in self.entries:
            entry.parameters.update( {'vbm': self.vbm_val})
        self.pd = DefectPhaseDiagram(self.entries, self.vbm_val, self.gap)

        # make Vac_As (q= -2) only defect test single-stable-charge exceptions
        self.extra_entry = DefectEntry(self.entries[5].defect.copy(), 100.)
        sep_entries = [ent for ent in self.entries if not (ent.name == 'Vac_As_mult4' and
                                                           ent.charge in [-2,-1,0,1,2])]
        sep_entries.append( self.extra_entry.copy())
        self.sep_pd = DefectPhaseDiagram( sep_entries, self.vbm_val, self.gap)

        # make Vac_As (q= -2) is incompatible for larger supercell
        ls_entries = self.entries[:]
        for entry in ls_entries:
            if entry.name == 'Vac_As_mult4' and entry.charge == -2.:
                entry.parameters['is_compatible'] = False
        self.pd_ls_fcTrue = DefectPhaseDiagram(ls_entries, self.vbm_val, self.gap, filter_compatible=True)
        self.pd_ls_fcFalse = DefectPhaseDiagram(ls_entries, self.vbm_val, self.gap, filter_compatible=False)

        # load complete dos for fermi energy solving
        with open(os.path.join(test_dir, "complete_dos.json"), "r") as f:
            dos_dict = json.load(f)
        tmp_dos = CompleteDos.from_dict(dos_dict)
        tmp_vbm = tmp_dos.get_cbm_vbm()[1]
        en_shift = self.vbm_val - tmp_vbm
        dos_dict["energies"] += en_shift
        dos_dict["efermi"] = self.vbm_val
        self.dos = CompleteDos.from_dict(dos_dict)


    def test_good_test_data(self):
        self.assertEqual(len(self.entries), 48)

    def test_suggest_charges(self):
        suggested_charges = self.pd.suggest_charges()
        for k in [
                    "Vac_As_mult4@0-1-2-3-4-5", "Sub_Ga_on_As_mult4@6-7-8-9-10-11", "Vac_Ga_mult4@12-13-14-15",
                    "Sub_As_on_Ga_mult4@16-17-18-19-20-21", "Int_Ga_mult1@22-23-24-25",
                    "Int_As_mult1@26-27-28-29-30-31-32-33-34", "Int_As_mult1@35-36-37-38-39-40-41-42-43",
                    "Int_Ga_mult1@44-45-46-47"
        ]:
            self.assertTrue(
                len(suggested_charges[k]) > 0, "Could not find any suggested charges for {} with band_gap of {}".format(
                    k, self.pd.band_gap))

        pd = DefectPhaseDiagram(self.entries, 2.6682, 1.0)
        suggested_charges = self.pd.suggest_charges()
        for k in ["Vac_As_mult4@0-1-2-3-4-5", "Vac_Ga_mult4@12-13-14-15"]:
            self.assertTrue(
                len(suggested_charges[k]) > 0, "Could not find any suggested charges for {} with band_gap of {}".format(
                    k, pd.band_gap))

        #test again but with only one charge state stable for Vac_As
        suggested_charges = self.sep_pd.suggest_charges()
        self.assertEqual( set(suggested_charges['Vac_As_mult4@0-43']), set([-4]))

    def test_suggest_larger_supercells(self):
        suggested_larger_cells = self.pd_ls_fcFalse.suggest_larger_supercells()
        self.assertEqual( suggested_larger_cells['Vac_As_mult4@0-1-2-3-4-5'], [-2])

        # raise error if filter_compatibile = True
        self.assertRaises( ValueError, self.pd_ls_fcTrue.suggest_larger_supercells)

    def test_entries(self):
        all_stable_entries = self.pd.all_stable_entries

        self.assertEqual(len(self.pd.defect_types), 8)
        self.assertEqual(len(all_stable_entries), sum([len(v) for v in self.pd.stable_charges.values()]))

        #test again but with only one charge state stable for Vac_As
        self.assertEqual( len(self.sep_pd.transition_level_map['Vac_As_mult4@0-43']), 0)
        self.assertEqual( len(self.sep_pd.stable_entries['Vac_As_mult4@0-43']), 1)
        self.assertEqual( len(self.sep_pd.finished_charges['Vac_As_mult4@0-43']), 2)

    def test_solve_for_fermi_energy(self):
        mu_elts = {Element("As"): -4.658070555, Element("Ga"): -3.7317319750000006}
        fermi_energy = self.pd.solve_for_fermi_energy( 100., mu_elts, self.dos)
        self.assertAlmostEqual( fermi_energy, 0.57224289)
        fermi_energy = self.pd.solve_for_fermi_energy( 1000., mu_elts, self.dos)
        self.assertAlmostEqual( fermi_energy, 0.68430342)

    def test_plot(self):
        #simple test that plot is produced
        p = self.pd.plot(saved=False)
        self.assertTrue( p)

if __name__ == "__main__":
    unittest.main()
