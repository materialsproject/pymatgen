import unittest
import os

from pymatgen.analysis.pourbaix.maker import PourbaixDiagram
from pymatgen.analysis.pourbaix.entry import PourbaixEntryIO


class TestPourbaixDiagram(unittest.TestCase):

    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements, entries) = PourbaixEntryIO.from_csv(os.path.join(module_dir,
                                                    "test_entries.csv"))
        self._entries = entries
        self._pd = PourbaixDiagram(entries)
        
    def test_pourbaix_diagram(self):
        self.assertEqual(len(self._pd.facets), 6, "Incorrect number of facets")

#class TestPourbaixDiagramMultiElement(unittest.TestCase):
#    
#    def setUp(self):
        
            
if __name__ == '__main__':
    unittest.main()
