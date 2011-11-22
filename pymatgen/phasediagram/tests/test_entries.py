import unittest
import os

from pymatgen.phasediagram.entries import PDEntryIO, PDEntry, GrandPotPDEntry
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Composition
import json

class PDEntryTestCase (unittest.TestCase):
    '''
    Test all functions using a ficitious entry
    '''
    def setUp(self):
        comp = Composition({Element('Li'):1,Element('Fe'):1,Element('O'):2})
        self.entry = PDEntry(comp,53)
        self.gpentry = GrandPotPDEntry(self.entry,{Element('O'):1.5})

    def test_get_energy(self):
        self.assertEquals(self.entry.energy,53,"Wrong energy!")
        self.assertEquals(self.gpentry.energy,50,"Wrong energy!")

    def test_get_energy_per_atom(self):
        self.assertEquals(self.entry.energy_per_atom,53.0/4,"Wrong energy per atom!")
        self.assertEquals(self.gpentry.energy_per_atom,50.0/2,"Wrong energy per atom!")

    def test_get_name(self):
        self.assertEquals(self.entry.name,'LiFeO2',"Wrong name!")
        self.assertEquals(self.gpentry.name,'LiFeO2',"Wrong name!")

    def test_get_composition(self):
        comp = self.entry.composition
        expected_comp = Composition({Element('Fe'):1,Element('Li'):1,Element('O'):2})
        self.assertEquals(comp,expected_comp,"Wrong composition!")
        comp = self.gpentry.composition
        expected_comp = Composition({Element('Fe'):1,Element('Li'):1})
        self.assertEquals(comp,expected_comp,"Wrong composition!")

    def test_is_element(self):
        self.assertFalse(self.entry.is_element, "Wrongly assigned as element!")
        self.assertFalse(self.gpentry.is_element, "Wrongly assigned as element!")
 
        
class PDEntryIOTestCase(unittest.TestCase):

    def test_read_csv(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        (elements,entries) = PDEntryIO.from_csv(os.path.join(module_dir,"pdentries_test.csv"))
        self.assertEqual(elements,[Element('Li'),Element('Fe'),Element('O')], "Wrong elements!")
        self.assertEqual(len(entries),492,"Wrong number of entries!")

if __name__ == '__main__':
    unittest.main()

