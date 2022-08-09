import os
import unittest
from pathlib import Path

import numpy as np
from monty.tempfile import ScratchDir

from pymatgen.core.periodic_table import Element
from pymatgen.io.phonopy import (
    CompletePhononDos,
    PhononBandStructure,
    PhononBandStructureSymmLine,
    Structure,
    get_complete_ph_dos,
    get_displaced_structures,
    get_gruneisen_ph_bs_symm_line,
    get_gruneisenparameter,
    get_ph_bs_symm_line,
    get_ph_dos,
    get_phonon_band_structure_from_fc,
    get_phonon_band_structure_symm_line_from_fc,
    get_phonon_dos_from_fc,
    get_phonopy_structure,
    get_pmg_structure,
    get_thermal_displacement_matrices,
)
from pymatgen.util.testing import PymatgenTest

try:
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS
except ImportError as ex:
    print(ex)
    Phonopy = None

test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "phonopy")


class PhonopyParserTest(PymatgenTest):
    def test_get_ph_bs(self):
        ph_bs = get_ph_bs_symm_line(os.path.join(test_dir, "NaCl_band.yaml"), has_nac=True)

        self.assertAlmostEqual(ph_bs.bands[1][10], 0.7753555184)
        self.assertAlmostEqual(ph_bs.bands[5][100], 5.2548379776)
        self.assertArrayEqual(ph_bs.bands.shape, (6, 204))
        self.assertArrayEqual(ph_bs.eigendisplacements.shape, (6, 204, 2, 3))
        self.assertArrayAlmostEqual(
            ph_bs.eigendisplacements[3][50][0],
            [0.0 + 0.0j, 0.14166569 + 0.04098339j, -0.14166569 - 0.04098339j],
        )
        self.assertTrue(ph_bs.has_eigendisplacements, True)
        self.assertArrayEqual(ph_bs.min_freq()[0].frac_coords, [0, 0, 0])
        self.assertAlmostEqual(ph_bs.min_freq()[1], -0.03700895020)
        self.assertTrue(ph_bs.has_imaginary_freq())
        self.assertFalse(ph_bs.has_imaginary_freq(tol=0.5))
        self.assertArrayAlmostEqual(ph_bs.asr_breaking(), [-0.0370089502, -0.0370089502, -0.0221388897])
        self.assertEqual(ph_bs.nb_bands, 6)
        self.assertEqual(ph_bs.nb_qpoints, 204)
        self.assertArrayAlmostEqual(ph_bs.qpoints[1].frac_coords, [0.01, 0, 0])
        self.assertTrue(ph_bs.has_nac)
        self.assertAlmostEqual(ph_bs.get_nac_frequencies_along_dir([1, 1, 0])[3], 4.6084532143)
        self.assertIsNone(ph_bs.get_nac_frequencies_along_dir([1, 0, 1]))
        self.assertArrayAlmostEqual(
            ph_bs.get_nac_eigendisplacements_along_dir([1, 1, 0])[3][1],
            [(0.1063906409128248 + 0j), 0j, 0j],
        )
        self.assertIsNone(ph_bs.get_nac_eigendisplacements_along_dir([1, 0, 1]))

    def test_get_ph_dos(self):
        dos = get_ph_dos(os.path.join(test_dir, "NaCl_total_dos.dat"))

        self.assertAlmostEqual(dos.densities[15], 0.0001665998)
        self.assertAlmostEqual(dos.frequencies[20], 0.0894965119)
        self.assertAlmostEqual(dos.get_interpolated_value(3.0), 1.2915532670115628)
        self.assertEqual(len(dos.frequencies), 201)
        self.assertEqual(len(dos.densities), 201)

    def test_get_complete_dos(self):
        cdos = get_complete_ph_dos(
            os.path.join(test_dir, "NaCl_partial_dos.dat"),
            os.path.join(test_dir, "NaCl_phonopy.yaml"),
        )
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
        symbols_pmg = {e.symbol for e in s_pmg.composition}
        symbols_pmg2 = {e.symbol for e in s_pmg2.composition}

        self.assertAlmostEqual(s_ph.get_cell()[1, 1], s_pmg.lattice._matrix[1, 1], 7)
        self.assertAlmostEqual(s_pmg.lattice._matrix[1, 1], s_pmg2.lattice._matrix[1, 1], 7)
        self.assertEqual(symbols_pmg, set(s_ph.symbols))
        self.assertEqual(symbols_pmg, symbols_pmg2)
        self.assertArrayAlmostEqual(coords_ph[3], s_pmg.frac_coords[3])
        self.assertArrayAlmostEqual(s_pmg.frac_coords[3], s_pmg2.frac_coords[3])
        self.assertEqual(s_ph.get_number_of_atoms(), s_pmg.num_sites)
        self.assertEqual(s_pmg.num_sites, s_pmg2.num_sites)


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class GetDisplacedStructuresTest(PymatgenTest):
    def test_get_displaced_structures(self):
        pmg_s = Structure.from_file(os.path.join(test_dir, "POSCAR-unitcell"), False)
        supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 2]]
        structures = get_displaced_structures(pmg_structure=pmg_s, atom_disp=0.01, supercell_matrix=supercell_matrix)

        self.assertEqual(len(structures), 49)
        self.assertArrayAlmostEqual(
            structures[4].frac_coords[0],
            np.array([0.10872682, 0.21783039, 0.12595286]),
            7,
        )
        self.assertArrayAlmostEqual(
            structures[-1].frac_coords[9],
            np.array([0.89127318, 0.78130015, 0.37404715]),
            7,
        )
        self.assertEqual(structures[0].num_sites, 128)
        self.assertEqual(structures[10].num_sites, 128)
        self.assertArrayAlmostEqual(structures[0].lattice._matrix, structures[8].lattice._matrix, 8)

        # test writing output
        with ScratchDir("."):
            structures = get_displaced_structures(
                pmg_structure=pmg_s,
                atom_disp=0.01,
                supercell_matrix=supercell_matrix,
                yaml_fname="test.yaml",
            )
            self.assertTrue(os.path.exists("test.yaml"))


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class TestPhonopyFromForceConstants(unittest.TestCase):
    def setUp(self) -> None:
        test_path = Path(test_dir)
        structure_file = test_path / "POSCAR-NaCl"
        fc_file = test_path / "FORCE_CONSTANTS"

        self.structure = Structure.from_file(structure_file)
        self.supercell_matrix = np.eye(3) * 2
        self.force_constants = parse_FORCE_CONSTANTS(fc_file)

    def test_get_phonon_dos_from_fc(self):
        dos = get_phonon_dos_from_fc(
            self.structure,
            self.supercell_matrix,
            self.force_constants,
            mesh_density=10.0,
        )

        self.assertTrue(dos, CompletePhononDos)
        self.assertEqual(len(dos.frequencies), 201)
        self.assertIn(Element.Na, dos.get_element_dos())
        self.assertIn(Element.Cl, dos.get_element_dos())

    def test_get_phonon_band_structure_from_fc(self):
        bs = get_phonon_band_structure_from_fc(
            self.structure,
            self.supercell_matrix,
            self.force_constants,
            mesh_density=10.0,
        )

        self.assertTrue(bs, PhononBandStructure)
        self.assertEqual(bs.nb_bands, 8)
        self.assertEqual(bs.nb_qpoints, 8)
        self.assertAlmostEqual(bs.bands[2][10], 3.887125285018674)

    def test_get_phonon_band_structure_symm_line_from_fc(self):
        bs = get_phonon_band_structure_symm_line_from_fc(
            self.structure,
            self.supercell_matrix,
            self.force_constants,
            line_density=5.0,
        )

        self.assertTrue(bs, PhononBandStructureSymmLine)
        self.assertEqual(bs.nb_bands, 24)
        self.assertEqual(bs.nb_qpoints, 48)
        self.assertAlmostEqual(bs.bands[2][10], 2.869229797603161)


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class TestGruneisen(unittest.TestCase):
    def test_ph_bs_symm_line(self):
        self.bs_symm_line_1 = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_band_Si.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_Si"),
            fit=True,
        )
        self.bs_symm_line_2 = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_band_Si.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_Si"),
            fit=False,
        )

        # check if a bit of the gruneisen parameters happens

        self.assertNotEqual(self.bs_symm_line_1.gruneisen[0][0], self.bs_symm_line_2.gruneisen[0][0])
        with self.assertRaises(ValueError):
            self.bs_symm_line_2 = get_gruneisen_ph_bs_symm_line(
                gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_eq_plus_minus_InP.yaml")
            )

    def test_gruneisen_parameter(self):
        self.gruneisenobject_Si = get_gruneisenparameter(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_Si.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/eq/POSCAR_Si"),
        )

        self.assertAlmostEqual(self.gruneisenobject_Si.frequencies[0][0], 0.2523831291)
        self.assertAlmostEqual(self.gruneisenobject_Si.gruneisen[0][0], -0.1190736091)

        # catch the exception when no structure is present
        with self.assertRaises(ValueError):
            get_gruneisenparameter(
                os.path.join(PymatgenTest.TEST_FILES_DIR, "gruneisen/gruneisen_mesh_InP_without_struct.yaml")
            )


@unittest.skipIf(Phonopy is None, "Phonopy not present")
class TestThermalDisplacementMatrices(PymatgenTest):
    def test_get_thermal_displacement_matrix(self):
        list_matrices = get_thermal_displacement_matrices(
            os.path.join(
                PymatgenTest.TEST_FILES_DIR, "thermal_displacement_matrices", "thermal_displacement_matrices.yaml"
            ),
            os.path.join(PymatgenTest.TEST_FILES_DIR, "thermal_displacement_matrices", "POSCAR"),
        )

        self.assertArrayAlmostEqual(
            list(list_matrices[0].thermal_displacement_matrix_cart[0]),
            [0.00516, 0.00613, 0.00415, -0.00011, -0.00158, -0.00081],
        )
        self.assertArrayAlmostEqual(
            list(list_matrices[0].thermal_displacement_matrix_cart[21]),
            [0.00488, 0.00497, 0.00397, -0.00070, -0.00070, 0.00144],
        )

        self.assertArrayAlmostEqual(
            list(list_matrices[0].thermal_displacement_matrix_cif[0]),
            [0.00457, 0.00613, 0.00415, -0.00011, -0.00081, -0.00082],
        )
        self.assertArrayAlmostEqual(
            list(list_matrices[0].thermal_displacement_matrix_cif[21]),
            [0.00461, 0.00497, 0.00397, -0.00070, 0.00002, 0.00129],
        )

        # check if correct number of temperatures has been read
        self.assertEqual(len(list_matrices), 31)

        self.assertAlmostEqual(list_matrices[-1].temperature, 300.0)
        self.assertAlmostEqual(list_matrices[0].temperature, 0.0)


if __name__ == "__main__":
    unittest.main()
