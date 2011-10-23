import unittest
import os
import StringIO

from pymatgen.phasediagram.entries import PDEntryIO
from pymatgen.phasediagram.pdmaker import PhaseDiagram
from pymatgen.phasediagram.plotter import PDPlotter, uniquelines

module_dir = os.path.dirname(os.path.abspath(__file__))

class PDPlotterTest(unittest.TestCase):

    def setUp(self):        
        (elements, entries) = PDEntryIO.from_csv(os.path.join(module_dir,"pdentries_test.csv"))
        self.pd = PhaseDiagram(entries)
        self.plotter = PDPlotter(self.pd)
        
    def test_get_pd_plot_data(self):
        (lines, labels, unstable_entries) = self.plotter.pd_plot_data
        self.assertEqual(len(lines), 22)
        self.assertEqual(len(labels), len(self.pd.stable_entries), "Incorrect number of lines generated!")
        self.assertEqual(len(unstable_entries), len(self.pd.all_entries) - len(self.pd.stable_entries), "Incorrect number of lines generated!")
        
    def test_write_image(self):
        
        output = StringIO.StringIO()
        self.plotter.write_image(output, "svg")
        self.assertNotEqual(output.getvalue(), "")
        
        
class UtilityFunctionTest(unittest.TestCase):
    
    def test_unique_lines(self):
        testdata = [[5, 53, 353], [399, 20, 52], [399, 400, 20], [13, 399, 52], [21, 400, 353], [393, 5, 353], [400, 393, 353], [393, 400, 399], [393, 13, 5], [13, 393, 399], [400, 17, 20], [21, 17, 400]]
        expected_ans = set([(399, 20), (21, 353), (17, 20), (400, 353), (21, 400), (17, 400), (400, 17), (393, 399), (20, 52), (400, 393), (5, 53), (400, 399), (5, 353), (13, 399), (393, 400), (13, 393), (13, 52), (393, 13), (393, 353), (399, 52), (53, 353), (13, 5), (400, 20), (393, 5), (21, 17), (399, 400)])
        self.assertEqual(uniquelines(testdata), expected_ans)
        

if __name__ == '__main__':
    unittest.main()

