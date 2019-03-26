from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.trajectory import Trajectory
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import numpy as np
import os

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class TrajectoryTest(PymatgenTest):

    def setUp(self):
        xdatcar = Xdatcar(os.path.join(test_dir, "Traj_XDATCAR"))
        self.traj = Trajectory.from_file(os.path.join(test_dir, "Traj_XDATCAR"))
        self.structures = xdatcar.structures

    def test_single_index_slice(self):
        self.assertTrue(all([self.traj[i] == self.structures[i] for i in range(0, len(self.structures), 19)]))

    def test_slice(self):
        sliced_traj = self.traj[2:99:3]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[2:99:3])

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

    def test_extend(self):
        traj = self.traj.copy()

        # Case of compatible trajectories
        compatible_traj = Trajectory.from_file(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_1"))
        traj.extend(compatible_traj)

        full_traj = Trajectory.from_file(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_Full"))
        compatible_success = self._check_traj_equality(self.traj, full_traj)

        # Case of incompatible trajectories
        traj = self.traj.copy()
        incompatible_traj = Trajectory.from_file(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_2"))
        incompatible_test_success=False
        try:
            traj.extend(incompatible_traj)
        except:
            incompatible_test_success=True

        self.assertTrue(compatible_success and incompatible_test_success)

    def test_length(self):
        self.assertTrue(len(self.traj) == len(self.structures))

    def test_displacements(self):
        poscar = Poscar.from_file(os.path.join(test_dir, "POSCAR"))
        structures = [poscar.structure]
        displacements = np.zeros((11, *np.shape(structures[-1].frac_coords)))
        for i in range(10):
            displacement = np.random.random_sample(np.shape(structures[-1].frac_coords)) / 20
            new_coords = displacement + structures[-1].frac_coords
            structures.append(Structure(structures[-1].lattice, structures[-1].species, new_coords))
            displacements[i+1, :, :] = displacement

        traj = Trajectory.from_structures(structures, constant_lattice=True)
        traj.to_displacements()

        self.assertTrue(np.allclose(traj.frac_coords, displacements))

    def test_changing_lattice(self):
        structure = self.structures[0]

        # Generate structures with different lattices
        structures = []
        for i in range(10):
            new_lattice = np.dot(structure.lattice.matrix, np.diag(1 + np.random.random_sample(3)/20))
            temp_struct = structure.copy()
            temp_struct.lattice = Lattice(new_lattice)
            structures.append(temp_struct)

        traj = Trajectory.from_structures(structures, constant_lattice=False)

        # Check if lattices were properly stored
        self.assertTrue(
            all([np.allclose(struct.lattice.matrix, structures[i].lattice.matrix) for i, struct in enumerate(traj)]))

    def test_to_from_dict(self):
        d = self.traj.as_dict()
        traj = Trajectory.from_dict(d)
        self.assertEqual(type(traj), Trajectory)

    def _check_traj_equality(self, traj_1, traj_2):
        if np.sum(np.square(np.subtract(traj_1.lattice, traj_2.lattice))) > 0.0001:
            return False

        if traj_1.species != traj_2.species:
            return False

        return all([i == j for i, j in zip(self.traj, traj_2)])


if __name__ == '__main__':
    import unittest
    unittest.main()
