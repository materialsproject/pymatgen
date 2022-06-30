import os
from copy import deepcopy

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

    def test_extend_site_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        num_frames = len(frac_coords)

        props_1 = {
            "selective_dynamics": [[False, False, False], [False, False, False]],
            "magmom": [5, 5],
        }
        props_2 = {
            "selective_dynamics": [[True, True, True], [False, False, False]],
            "magmom": [5, 5],
        }
        props_3 = [
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

        traj_1 = Trajectory(lattice, species, frac_coords, site_properties=props_1)
        traj_2 = Trajectory(lattice, species, frac_coords, site_properties=props_2)
        traj_3 = Trajectory(lattice, species, frac_coords, site_properties=props_3)
        traj_4 = Trajectory(lattice, species, frac_coords, site_properties=None)

        # const & const (both constant and the same site properties)
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_1)
        expected_site_props = props_1
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # const & const (both constant but different site properties)
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_2)
        expected_site_props = [props_1] * num_frames + [props_2] * num_frames
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # const & changing
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_3)
        expected_site_props = [props_1] * num_frames + props_3
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # const & none
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_4)
        expected_site_props = [props_1] * num_frames + [None] * num_frames
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # changing & const
        traj_combined = deepcopy(traj_3)
        traj_combined.extend(traj_1)
        expected_site_props = props_3 + [props_1] * num_frames
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # changing & changing
        traj_combined = deepcopy(traj_3)
        traj_combined.extend(traj_3)
        expected_site_props = props_3 + props_3
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # changing & none
        traj_combined = deepcopy(traj_3)
        traj_combined.extend(traj_4)
        expected_site_props = props_3 + [None] * num_frames
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # none & const
        traj_combined = deepcopy(traj_4)
        traj_combined.extend(traj_1)
        expected_site_props = [None] * num_frames + [props_1] * num_frames
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # none & changing
        traj_combined = deepcopy(traj_4)
        traj_combined.extend(traj_3)
        expected_site_props = [None] * num_frames + props_3
        self.assertEqual(traj_combined.site_properties, expected_site_props)

        # none & none
        traj_combined = deepcopy(traj_4)
        traj_combined.extend(traj_4)
        expected_site_props = None
        self.assertEqual(traj_combined.site_properties, expected_site_props)

    def test_extend_frame_props(self):
        lattice = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        species = ["Si", "Si"]
        frac_coords = [
            [[0, 0, 0], [0.5, 0.5, 0.5]],
            [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
            [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
        ]

        energy_1 = [-3, -3.9, -4.1]
        energy_2 = [-4.2, -4.25, -4.3]
        pressure_2 = [2, 2.5, 2.5]

        # energy only properties
        props_1 = [{"energy": e} for e in energy_1]
        traj_1 = Trajectory(lattice, species, frac_coords, frame_properties=props_1)

        # energy and pressure properties
        props_2 = [{"energy": e, "pressure": p} for e, p in zip(energy_2, pressure_2)]
        traj_2 = Trajectory(lattice, species, frac_coords, frame_properties=props_2)

        # no properties
        traj_3 = Trajectory(lattice, species, frac_coords, frame_properties=None)

        # test combining two with different properties
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_2)
        expected_props = props_1 + props_2
        self.assertEqual(traj_combined.frame_properties, expected_props)

        # test combining two where one has properties and the other does not
        traj_combined = deepcopy(traj_1)
        traj_combined.extend(traj_3)
        expected_props = props_1 + [None] * len(frac_coords)
        self.assertEqual(traj_combined.frame_properties, expected_props)

        # test combining two both of which have no properties
        traj_combined = deepcopy(traj_3)
        traj_combined.extend(traj_3)
        self.assertEqual(traj_combined.frame_properties, None)

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

        # Check if the file is written correctly when lattice is not constant.
        traj.write_Xdatcar(filename="traj_test_XDATCAR")
        # Load trajectory from written xdatcar and compare to original
        written_traj = Trajectory.from_file("traj_test_XDATCAR", constant_lattice=False)
        self._check_traj_equality(traj, written_traj)
        os.remove("traj_test_XDATCAR")

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
