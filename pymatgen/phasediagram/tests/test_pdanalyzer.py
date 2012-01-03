import unittest
import os

from pymatgen.core.structure import Composition
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
from pymatgen.phasediagram.entries import PDEntryIO

class  PDAnalyzerTest(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements,entries) = PDEntryIO.from_csv(os.path.join(module_dir,"pdentries_test.csv"))
        self.pd = PhaseDiagram(entries)
        self.analyzer = PDAnalyzer(self.pd)

    def test_get_e_above_hull(self):
        for entry in self.pd.stable_entries:
            self.assertLess(self.analyzer.get_e_above_hull(entry), 1e-11, "Stable entries should have e above hull of zero!")
            
    def test_get_equilibrium_reaction_energy(self):
        for entry in self.pd.stable_entries:
            self.assertLessEqual(self.analyzer.get_equilibrium_reaction_energy(entry), 0, "Stable entries should have negative equilibrium reaction energy!")
    
    def test_get_decomposition(self):
        for entry in self.pd.stable_entries:
            self.assertEquals(len(self.analyzer.get_decomposition(entry.composition)), 1, "Stable composition should have only 1 decomposition!")
        dim = len(self.pd.elements)
        for entry in self.pd.all_entries:
            ndecomp = len(self.analyzer.get_decomposition(entry.composition))
            self.assertTrue(ndecomp > 0 and ndecomp <= dim, "The number of decomposition phases can at most be equal to the number of components.")
        
        #Just to test decomp for a ficitious composition
        ansdict = {entry.composition.formula: amt for entry, amt in self.analyzer.get_decomposition(Composition.from_formula("Li3Fe7O11")).items()}
        expected_ans = {"Fe2 O2" : 0.0952380952380949, "Li1 Fe1 O2": 0.5714285714285714, "Fe6 O8": 0.33333333333333393}
        for k, v in expected_ans.items():
            self.assertAlmostEqual(ansdict[k], v)

        
    def test_get_transition_chempots(self):
        for el in self.pd.elements:
            self.assertLessEqual(len(self.analyzer.get_transition_chempots(el)), len(self.pd.facets))

    def test_get_element_profile(self):
        for el in self.pd.elements:
            for entry in self.pd.stable_entries:
                if not (entry.composition.is_element):
                    self.assertLessEqual(len(self.analyzer.get_element_profile(el, entry.composition)), len(self.pd.facets))

if __name__ == '__main__':
    unittest.main()

