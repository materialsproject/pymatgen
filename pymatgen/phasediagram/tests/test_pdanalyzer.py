import unittest
import os
import random

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
    
    def test_get_decomposition(self):
        for entry in self.pd.stable_entries:
            self.assertEquals(len(self.analyzer.get_decomposition(entry.composition)), 1, "Stable composition should have only 1 decomposition!")
        dim = len(self.pd.elements)
        for entry in self.pd.all_entries:
            ndecomp = len(self.analyzer.get_decomposition(entry.composition))
            self.assertTrue(ndecomp > 0 and ndecomp <= dim, "The number of decomposition phases can at most be equal to the number of components.")

    def test_get_transition_chempots(self):
        el = self.pd.elements[random.randint(0,2)]
        self.assertLessEqual(len(self.analyzer.get_transition_chempots(el)), len(self.pd.facets))

    def test_get_element_profile(self):
        el = self.pd.elements[random.randint(0,2)]
        for entry in self.pd.stable_entries:
            if not (entry.composition.is_element):
                self.assertLessEqual(len(self.analyzer.get_element_profile(el, entry.composition)), len(self.pd.facets))

if __name__ == '__main__':
    unittest.main()

