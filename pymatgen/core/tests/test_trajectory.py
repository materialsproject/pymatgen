from pymatgen.util.testing import PymatgenTest
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.trajectory import Trajectory
from pymatgen.core.structure import Structure
import numpy as np
import os

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files')


class TrajectoryTest(PymatgenTest):

    def setUp(self):
        xdatcar = Xdatcar(os.path.join(test_dir, "Traj_XDATCAR"))
        self.traj = Trajectory.from_structures(xdatcar.structures, const_lattice=True)
        self.structures = xdatcar.structures

    def testSingleIndexSlice(self):
        self.assertTrue(all([self.traj[i] == self.structures[i] for i in range(0, len(self.structures), 19)]))

    def testSlice(self):
        sliced_traj = self.traj[2:99:3]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[2:99:3])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            self.assertTrue(all([sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj))]))
        else:
            self.assertTrue(False)

    def testConversion(self):
        # Convert to displacements and back. Check structures
        self.traj.to_displacements()
        self.traj.to_positions()

        self.assertTrue(all([struct == self.structures[i] for i, struct in enumerate(self.traj)]))

    def testCopy(self):
        traj_copy = self.traj.copy()
        self.assertTrue(all([i==j for i, j in zip(self.traj, traj_copy)]))

    def testExtend(self):
        traj = self.traj.copy()

        # Case of compatible trajectories
        xdatcar = Xdatcar(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_1"))
        compatible_traj = Trajectory.from_structures(xdatcar.structures)
        traj.extend(compatible_traj)

        xdatcar = Xdatcar(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_Full"))
        full_traj = Trajectory.from_structures(xdatcar.structures)
        compatible_success = self.check_traj_equality(self.traj, full_traj)

        # Case off incompatible trajectories
        traj = self.traj.copy()
        xdatcar = Xdatcar(os.path.join(test_dir, "Traj_Combine_Test_XDATCAR_1"))
        incompatible_traj = Trajectory.from_structures(xdatcar.structures)
        traj.extend(incompatible_traj)
        incompatible_success = self.check_traj_equality(self.traj, traj)

        self.assertTrue(compatible_success and incompatible_success)

    def testLength(self):
        self.assertTrue(len(self.traj) == len(self.structures))

    def testDisplacements(self):
        poscar = Poscar.from_file(os.path.join(test_dir, "POSCAR"))
        structures = [poscar.structure]
        displacements = np.zeros((11, *np.shape(structures[-1].frac_coords)))
        for i in range(10):
            displacement = np.random.random_sample(np.shape(structures[-1].frac_coords)) / 20
            new_coords = displacement + structures[-1].frac_coords
            structures.append(Structure(structures[-1].lattice, structures[-1].species, new_coords))
            displacements[i+1, :, :] = displacement

        traj = Trajectory.from_structures(structures, const_lattice=True)
        traj.to_displacements()

        self.assertTrue(np.allclose(traj.frac_coords, displacements))

    def check_traj_equality(self, traj_1, traj_2):
        if np.sum(np.square(np.subtract(traj_1.lattice, traj_2.lattice))) > 0.0001:
            return False

        if traj_1.species != traj_2.species:
            return False

        return all([i==j for i, j in zip(self.traj, traj_2)])


if __name__ == '__main__':
    import unittest
    unittest.main()