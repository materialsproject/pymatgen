import os
import unittest
import sys

from pymatgen.core.periodic_table import Element
from pymatgen.util.testing import PymatgenTest
from pymatgen.io.phonopy import *

from monty.tempfile import ScratchDir

if sys.version_info >= (3, 0):
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
        from phonopy.file_IO import write_disp_yaml
    except ImportError:
        Phonopy = None
else:
    Phonopy = None

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "phonopy")


class PhonopyParserTest(PymatgenTest):
    def test_get_ph_bs(self):
        ph_bs = get_ph_bs_symm_line(os.path.join(test_dir, 'NaCl_band.yaml'),
                                    has_nac=True)

        self.assertAlmostEqual(ph_bs.bands[1][10], 0.7753555184)
        self.assertAlmostEqual(ph_bs.bands[5][100], 5.2548379776)
        self.assertArrayEqual(ph_bs.bands.shape, (6, 204))
        self.assertArrayEqual(ph_bs.eigendisplacements.shape, (6, 204, 2, 3))
        self.assertArrayAlmostEqual(ph_bs.eigendisplacements[3][50][0],
                                    [0. + 0.j, 0.14166569 + 0.04098339j,
                                     -0.14166569 - 0.04098339j])
        self.assertTrue(ph_bs.has_eigendisplacements, True)
        self.assertArrayEqual(ph_bs.min_freq()[0].frac_coords, [0, 0, 0])
        self.assertAlmostEqual(ph_bs.min_freq()[1], -0.03700895020)
        self.assertTrue(ph_bs.has_imaginary_freq())
        self.assertFalse(ph_bs.has_imaginary_freq(tol=0.5))
        self.assertArrayAlmostEqual(ph_bs.asr_breaking(),
                                    [-0.0370089502, -0.0370089502,
                                     -0.0221388897])
        self.assertEqual(ph_bs.nb_bands, 6)
        self.assertEqual(ph_bs.nb_qpoints, 204)
        self.assertArrayAlmostEqual(ph_bs.qpoints[1].frac_coords, [0.01, 0, 0])
        self.assertTrue(ph_bs.has_nac)
        self.assertAlmostEqual(
            ph_bs.get_nac_frequencies_along_dir([1, 1, 0])[3], 4.6084532143)
        self.assertIsNone(ph_bs.get_nac_frequencies_along_dir([1, 1, 1]))
        self.assertArrayAlmostEqual(
            ph_bs.get_nac_eigendisplacements_along_dir([1, 1, 0])[3][1],
            [(0.1063906409128248 + 0j), 0j, 0j])
        self.assertIsNone(ph_bs.get_nac_eigendisplacements_along_dir([1, 1, 1]))

    def test_get_ph_dos(self):
        dos = get_ph_dos(os.path.join(test_dir, 'NaCl_total_dos.dat'))

        self.assertAlmostEqual(dos.densities[15], 0.0001665998)
        self.assertAlmostEqual(dos.frequencies[20], 0.0894965119)
        self.assertAlmostEqual(dos.get_interpolated_value(3.),
                               1.2915532670115628)
        self.assertEqual(len(dos.frequencies), 201)
        self.assertEqual(len(dos.densities), 201)

    def test_get_complete_dos(self):
        cdos = get_complete_ph_dos(
            os.path.join(test_dir, 'NaCl_partial_dos.dat'),
            os.path.join(test_dir, 'NaCl_phonopy.yaml'))
        site_Na = cdos.structure[0]
        site_Cl = cdos.structure[1]

        self.assertEqual(len(cdos.frequencies), 201)
        self.assertAlmostEqual(cdos.pdos[site_Na][30], 0.008058208)
        self.assertAlmostEqual(cdos.pdos[site_Cl][30], 0.0119040783)

        self.assertIn(Element.Na, cdos.get_element_dos())
        self.assertIn(Element.Cl, cdos.get_element_dos())


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class StructureConversionTest(PymatgenTest):
    def test_structure_conversion(self):
        s_pmg = PymatgenTest.get_structure("LiFePO4")
        s_ph = get_phonopy_structure(s_pmg)
        s_pmg2 = get_pmg_structure(s_ph)

        coords_ph = s_ph.get_scaled_positions()
        symbols_pmg = set([e.symbol for e in s_pmg.composition.keys()])
        symbols_pmg2 = set([e.symbol for e in s_pmg2.composition.keys()])

        self.assertAlmostEqual(s_ph.get_cell()[1, 1],
                               s_pmg.lattice._matrix[1, 1], 7)
        self.assertAlmostEqual(s_pmg.lattice._matrix[1, 1],
                               s_pmg2.lattice._matrix[1, 1], 7)
        self.assertEqual(symbols_pmg, set(s_ph.symbols))
        self.assertEqual(symbols_pmg, symbols_pmg2)
        self.assertArrayAlmostEqual(coords_ph[3], s_pmg.frac_coords[3])
        self.assertArrayAlmostEqual(s_pmg.frac_coords[3], s_pmg2.frac_coords[3])
        self.assertEqual(s_ph.get_number_of_atoms(), s_pmg.num_sites)
        self.assertEqual(s_pmg.num_sites, s_pmg2.num_sites)


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class GetDisplacedStructuresTest(PymatgenTest):
    def test_get_displaced_structures(self):
        pmg_s = Structure.from_file(os.path.join(test_dir, "POSCAR-unitcell"),
                                    False)
        supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 2]]
        structures = get_displaced_structures(pmg_structure=pmg_s,
                                              atom_disp=0.01,
                                              supercell_matrix=supercell_matrix)

        self.assertEqual(len(structures), 49)
        self.assertArrayAlmostEqual(structures[4].frac_coords[0],
                                    np.array(
                                        [0.10872682, 0.21783039, 0.12595286]),
                                    7)
        self.assertArrayAlmostEqual(structures[-1].frac_coords[9],
                                    np.array(
                                        [0.89127318, 0.78130015, 0.37404715]),
                                    7)
        self.assertEqual(structures[0].num_sites, 128)
        self.assertEqual(structures[10].num_sites, 128)
        self.assertArrayAlmostEqual(structures[0].lattice._matrix,
                                    structures[8].lattice._matrix, 8)

        # test writing output
        with ScratchDir("."):
            structures = get_displaced_structures(pmg_structure=pmg_s,
                                                  atom_disp=0.01,
                                                  supercell_matrix=supercell_matrix,
                                                  yaml_fname="test.yaml")
            self.assertTrue(os.path.exists("test.yaml"))


if __name__ == '__main__':
    unittest.main()
