import os

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.util.testing import PymatgenTest


class TrajectoryTest(PymatgenTest):
    def setUp(self):
        xdatcar = Xdatcar(os.path.join(PymatgenTest.TEST_FILES_DIR, "Traj_XDATCAR"))
        self.traj = Trajectory.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "Traj_XDATCAR"))
        self.structures = xdatcar.structures

    def _check_traj_equality(self, traj_1, traj_2):
        if np.sum(np.square(np.subtract(traj_1.lattice, traj_2.lattice))) > 0.0001:
            return False

        if traj_1.species != traj_2.species:
            return False

        return all(i == j for i, j in zip(self.traj, traj_2))

    def test_single_index_slice(self):
        self.assertTrue(all([self.traj[i] == self.structures[i] for i in range(0, len(self.structures), 19)]))

    def test_slice(self):
        sliced_traj = self.traj[2:99:3]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[2:99:3])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            self.assertTrue(all([sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj))]))
        else:
            self.assertTrue(False)

        sliced_traj = self.traj[:-4:2]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[:-4:2])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            self.assertTrue(all([sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj))]))
        else:
            self.assertTrue(False)

    def test_list_slice(self):
        sliced_traj = self.traj[[10, 30, 70]]
        sliced_traj_from_structs = Trajectory.from_structures([self.structures[i] for i in [10, 30, 70]])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            self.assertTrue(all([sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj))]))
        else:
            self.assertTrue(False)

    def test_conversion(self):
        # Convert to displacements and back. Check structures
        self.traj.to_displacements()
        self.traj.to_positions()

        self.assertTrue(all([struct == self.structures[i] for i, struct in enumerate(self.traj)]))

    def test_copy(self):
        traj_copy = self.traj.copy()
        self.assertTrue(all([i == j for i, j in zip(self.traj, traj_copy)]))

    def test_site_properties(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]
        site_properties = [
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [6, 6],
            },
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        traj = Trajectory(lattice, species, frac_coords, site_properties=site_properties)

        # compare the overall site properties list
        self.assertEqual(traj.site_properties, site_properties)
        # # compare the site properties after slicing
        self.assertEqual(traj[0].site_properties, site_properties[0])
        self.assertEqual(traj[1:].site_properties, site_properties[1:])

    def test_frame_properties(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]
        site_properties = [
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [6, 6],
            },
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        frame_properties = {"energy_per_atom": [-3.0001, -3.0971, -3.0465]}

        traj = Trajectory(
            lattice,
            species,
            frac_coords,
            site_properties=site_properties,
            frame_properties=frame_properties,
        )
        # compare the overall site properties list
        self.assertEqual(traj.frame_properties, frame_properties)
        # compare the site properties after slicing
        expected_output = {"energy_per_atom": [-3.0971, -3.0465]}
        self.assertEqual(traj[1:].frame_properties, expected_output)

    def test_extend(self):
        traj = self.traj.copy()

        # Case of compatible trajectories
        compatible_traj = Trajectory.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "Traj_Combine_Test_XDATCAR_1"))
        traj.extend(compatible_traj)

        full_traj = Trajectory.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "Traj_Combine_Test_XDATCAR_Full"))
        compatible_success = self._check_traj_equality(self.traj, full_traj)

        # Case of incompatible trajectories
        traj = self.traj.copy()
        incompatible_traj = Trajectory.from_file(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "Traj_Combine_Test_XDATCAR_2")
        )
        incompatible_test_success = False
        try:
            traj.extend(incompatible_traj)
        except Exception:
            incompatible_test_success = True

        self.assertTrue(compatible_success and incompatible_test_success)

    def test_extend_no_site_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        # Trajectory with no site_properties
        traj_1 = Trajectory(lattice, species, frac_coords)
        traj_2 = Trajectory(lattice, species, frac_coords)

        # Test combining two trajectories with no site properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        self.assertEqual(traj_combined.site_properties, None)

    def test_extend_equivalent_site_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        # Trajectories with constant site properties
        site_properties_1 = [
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            }
        ]
        traj_1 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_1)

        site_properties_2 = [
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            }
        ]
        traj_2 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_2)

        # Test combining two trajectories with similar site_properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        self.assertEqual(traj_combined.site_properties, site_properties_1)

    def test_extend_inequivalent_site_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        # Trajectories with constant site properties
        site_properties_1 = [
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [5, 5],
            }
        ]
        traj_1 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_1)

        site_properties_2 = [
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            }
        ]
        traj_2 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_2)

        # Test combining two trajectories with similar site_properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        expected_site_props = [
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[False, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # Trajectory with const site_properties and trajectory with changing site properties
        site_properties_1 = [
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            }
        ]
        traj_1 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_1)

        site_properties_2 = [
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        traj_2 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_2)

        # Test combining two trajectories with similar site_properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        expected_site_props = [
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # The other way around
        traj_combined = traj_2.copy()
        traj_combined.extend(traj_1)
        expected_site_props = [
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, False], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # Trajectory with no and trajectory with changing site properties
        site_properties_1 = None
        traj_1 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_1)

        site_properties_2 = [
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        traj_2 = Trajectory(lattice, species, frac_coords, site_properties=site_properties_2)

        # Test combining two trajectories with similar site_properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        expected_site_props = [
            None,
            None,
            None,
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
        ]
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # The other way around
        traj_combined = traj_2.copy()
        traj_combined.extend(traj_1)
        expected_site_props = [
            {
                "selective_dynamics": [[False, True, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, False, True], [False, False, False]],
                "magmom": [5, 5],
            },
            {
                "selective_dynamics": [[True, True, False], [False, False, False]],
                "magmom": [5, 5],
            },
            None,
            None,
            None,
        ]
        self.assertEqual(traj_combined.site_properties, expected_site_props)

    def test_extend_no_frame_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        # Trajectory with no site_properties
        traj_1 = Trajectory(lattice, species, frac_coords)
        traj_2 = Trajectory(lattice, species, frac_coords)

        # Test combining two trajectories with no site properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        self.assertEqual(traj_combined.frame_properties, None)

    def test_extend_frame_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        # Trajectories with constant site properties
        frame_properties_1 = {"energy": [-3, -3.9, -4.1]}
        traj_1 = Trajectory(lattice, species, frac_coords, frame_properties=frame_properties_1)

        frame_properties_2 = {"energy": [-4.2, -4.25, -4.3]}
        traj_2 = Trajectory(lattice, species, frac_coords, frame_properties=frame_properties_2)

        # Test combining two trajectories with similar site_properties
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_2)
        expected_frame_properties = {"energy": [-3, -3.9, -4.1, -4.2, -4.25, -4.3]}
        self.assertEqual(traj_combined.frame_properties, expected_frame_properties)

        # Mismatched frame propertied
        frame_properties_3 = {"energy": [-4.2, -4.25, -4.3], "pressure": [2, 2.5, 2.5]}
        traj_3 = Trajectory(lattice, species, frac_coords, frame_properties=frame_properties_3)
        traj_combined = traj_1.copy()
        traj_combined.extend(traj_3)
        expected_frame_properties = {
            "energy": [-3, -3.9, -4.1, -4.2, -4.25, -4.3],
            "pressure": [None, None, None, 2, 2.5, 2.5],
        }
        self.assertEqual(traj_combined.frame_properties, expected_frame_properties)

    def test_length(self):
        self.assertTrue(len(self.traj) == len(self.structures))

    def test_displacements(self):
        poscar = Poscar.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR, "POSCAR"))
        structures = [poscar.structure]
        displacements = np.zeros((11, *np.shape(structures[-1].frac_coords)))
        for i in range(10):
            displacement = np.random.random_sample(np.shape(structures[-1].frac_coords)) / 20
            new_coords = displacement + structures[-1].frac_coords
            structures.append(Structure(structures[-1].lattice, structures[-1].species, new_coords))
            displacements[i + 1, :, :] = displacement

        traj = Trajectory.from_structures(structures, constant_lattice=True)
        traj.to_displacements()

        self.assertTrue(np.allclose(traj.frac_coords, displacements))

    def test_variable_lattice(self):
        structure = self.structures[0]

        # Generate structures with different lattices
        structures = []
        for i in range(10):
            new_lattice = np.dot(structure.lattice.matrix, np.diag(1 + np.random.random_sample(3) / 20))
            temp_struct = structure.copy()
            temp_struct.lattice = Lattice(new_lattice)
            structures.append(temp_struct)

        traj = Trajectory.from_structures(structures, constant_lattice=False)

        # Check if lattices were properly stored
        self.assertTrue(
            all(np.allclose(struct.lattice.matrix, structures[i].lattice.matrix) for i, struct in enumerate(traj))
        )

    def test_to_from_dict(self):
        d = self.traj.as_dict()
        traj = Trajectory.from_dict(d)
        self.assertEqual(type(traj), Trajectory)

    def test_xdatcar_write(self):
        self.traj.write_Xdatcar(filename="traj_test_XDATCAR")
        # Load trajectory from written xdatcar and compare to original
        written_traj = Trajectory.from_file("traj_test_XDATCAR")
        self._check_traj_equality(self.traj, written_traj)
        os.remove("traj_test_XDATCAR")


if __name__ == "__main__":
    import unittest

    unittest.main()
