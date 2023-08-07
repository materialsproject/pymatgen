from __future__ import annotations

import copy
import os

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Molecule, Structure
from pymatgen.core.trajectory import Trajectory
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Xdatcar
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestTrajectory(PymatgenTest):
    def setUp(self):
        xdatcar = Xdatcar(f"{TEST_FILES_DIR}/Traj_XDATCAR")
        self.traj = Trajectory.from_file(f"{TEST_FILES_DIR}/Traj_XDATCAR")
        self.structures = xdatcar.structures

        out = QCOutput(f"{TEST_FILES_DIR}/molecules/new_qchem_files/ts.out")
        last_mol = out.data["molecule_from_last_geometry"]
        species = last_mol.species
        coords = out.data["geometries"]

        self.molecules = []
        for c in coords:
            mol = Molecule(species, c, charge=int(last_mol.charge), spin_multiplicity=int(last_mol.spin_multiplicity))
            self.molecules.append(mol)

        self.traj_mols = Trajectory(
            species=species,
            coords=coords,
            charge=int(last_mol.charge),
            spin_multiplicity=int(last_mol.spin_multiplicity),
        )

    def _check_traj_equality(self, traj_1, traj_2):
        if traj_1.lattice is not None and not np.allclose(traj_1.lattice, traj_2.lattice):
            return False

        if traj_1.species != traj_2.species:
            return False

        return all(i == j for i, j in zip(self.traj, traj_2))

    def _get_lattice_species_and_coords(self):
        lattice = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        species = ["Si", "Si"]
        coords = np.asarray(
            [
                [[0, 0, 0], [0.5, 0.5, 0.5]],
                [[0.1, 0.1, 0.1], [0.6, 0.6, 0.6]],
                [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
            ]
        )
        return lattice, species, coords

    def _get_species_and_coords(self):
        species = ["C", "O"]
        coords = np.asarray(
            [
                [[1.5709474478, -0.16099953, 0.0], [1.9291378639, -1.2161950538, 0.0]],
                [[1.5688628148, -0.1548583957, 0.0], [1.9312224969, -1.2223361881, 0.0]],
                [[1.5690858055, -0.1555153055, 0.0], [1.9309995062, -1.2216792783, 0.0]],
            ]
        )
        return species, coords, 0, 1

    def test_single_index_slice(self):
        assert all(self.traj[i] == self.structures[i] for i in range(0, len(self.structures), 19))
        assert all(self.traj_mols[i] == self.molecules[i] for i in range(0, len(self.molecules)))

    def test_slice(self):
        sliced_traj = self.traj[2:99:3]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[2:99:3])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            assert all(sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

        sliced_traj = self.traj[:-4:2]
        sliced_traj_from_structs = Trajectory.from_structures(self.structures[:-4:2])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            assert all(sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

        sliced_traj = self.traj_mols[0:2]
        sliced_traj_from_mols = Trajectory.from_molecules(self.molecules[0:2])

        if len(sliced_traj) == len(sliced_traj_from_mols):
            assert all(sliced_traj[i] == sliced_traj_from_mols[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

        sliced_traj = self.traj_mols[:-2]
        sliced_traj_from_mols = Trajectory.from_molecules(self.molecules[:-2])

        if len(sliced_traj) == len(sliced_traj_from_mols):
            assert all(sliced_traj[i] == sliced_traj_from_mols[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

    def test_list_slice(self):
        sliced_traj = self.traj[[10, 30, 70]]
        sliced_traj_from_structs = Trajectory.from_structures([self.structures[i] for i in [10, 30, 70]])

        if len(sliced_traj) == len(sliced_traj_from_structs):
            assert all(sliced_traj[i] == sliced_traj_from_structs[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

        sliced_traj = self.traj_mols[[1, 3]]
        sliced_traj_from_mols = Trajectory.from_molecules([self.molecules[i] for i in [1, 3]])

        if len(sliced_traj) == len(sliced_traj_from_mols):
            assert all(sliced_traj[i] == sliced_traj_from_mols[i] for i in range(len(sliced_traj)))
        else:
            raise AssertionError

    def test_conversion(self):
        # Convert to displacements and back, and then check structures.
        self.traj.to_displacements()
        self.traj.to_positions()

        assert all(struct == self.structures[i] for i, struct in enumerate(self.traj))

        self.traj_mols.to_displacements()
        self.traj_mols.to_positions()

        assert all(mol == self.molecules[i] for i, mol in enumerate(self.traj_mols))

    def test_site_properties(self):
        lattice, species, coords = self._get_lattice_species_and_coords()

        props = [
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
        traj = Trajectory(lattice=lattice, species=species, coords=coords, site_properties=props)

        # compare the overall site properties list
        assert traj.site_properties == props

        # compare the site properties after slicing
        assert traj[0].site_properties == props[0]
        assert traj[1:].site_properties == props[1:]

        species, coords, charge, spin = self._get_species_and_coords()

        props = [
            {"test": [[True, True, True], [False, False, False]]},
            {"test": [[False, False, False], [False, False, False]]},
            {"test": [[True, True, True], [False, False, False]]},
        ]
        traj = Trajectory(
            species=species,
            coords=coords,
            charge=charge,
            spin_multiplicity=spin,
            site_properties=props,
        )

        # compare the overall site properties list
        assert traj.site_properties == props

        # compare the site properties after slicing
        assert traj[0].site_properties == props[0]
        assert traj[1:].site_properties == props[1:]

    def test_frame_properties(self):
        lattice, species, coords = self._get_lattice_species_and_coords()

        props = [{"energy_per_atom": e} for e in [-3.0001, -3.0971, -3.0465]]

        traj = Trajectory(lattice=lattice, species=species, coords=coords, frame_properties=props)

        # compare the overall site properties
        assert traj.frame_properties == props

        # compare the site properties after slicing
        expected = props[1:]
        assert traj[1:].frame_properties == expected

        species, coords, charge, spin = self._get_species_and_coords()

        props = [{"SCF_energy_in_the_final_basis_set": e} for e in [-113.3256885788, -113.3260019471, -113.326006415]]

        traj = Trajectory(
            species=species,
            coords=coords,
            charge=charge,
            spin_multiplicity=spin,
            frame_properties=props,
        )

        # compare the overall site properties
        assert traj.frame_properties == props

        # compare the site properties after slicing
        expected = props[1:]
        assert traj[1:].frame_properties == expected

    def test_extend(self):
        traj = copy.deepcopy(self.traj)

        # Case of compatible trajectories
        compatible_traj = Trajectory.from_file(f"{TEST_FILES_DIR}/Traj_Combine_Test_XDATCAR_1")
        traj.extend(compatible_traj)

        full_traj = Trajectory.from_file(f"{TEST_FILES_DIR}/Traj_Combine_Test_XDATCAR_Full")
        compatible_success = self._check_traj_equality(self.traj, full_traj)

        # Case of incompatible trajectories
        traj = copy.deepcopy(self.traj)
        incompatible_traj = Trajectory.from_file(f"{TEST_FILES_DIR}/Traj_Combine_Test_XDATCAR_2")
        incompatible_test_success = False
        try:
            traj.extend(incompatible_traj)
        except Exception:
            incompatible_test_success = True

        assert compatible_success
        assert incompatible_test_success

        traj = copy.deepcopy(self.traj_mols)

        # Case of compatible trajectories
        compatible_traj = Trajectory(
            species=traj.species,
            coords=[
                [
                    [-1.46958173, -0.47370158, -0.03391061],
                    [-0.79757102, 0.48588802, 0.94508206],
                    [0.50256405, 0.8947604, 0.47698504],
                    [1.56101382, 0.13356272, 0.79931048],
                    [1.43897567, -0.8642765, 1.56363034],
                    [2.66882238, 0.48431336, 0.30635727],
                    [-2.72606146, -0.81552889, 0.39696593],
                    [3.307822, -1.01132269, 1.26654957],
                    [-0.81092724, -1.35590014, -0.1458541],
                    [-1.48634516, 0.02121279, -1.02465009],
                    [-0.71212347, 0.03008471, 1.93272477],
                    [-1.37888759, 1.40819443, 1.02143913],
                    [-4.79241099, 0.80275103, -0.39852432],
                    [-4.28509927, -1.03484764, 0.86348452],
                ]
            ],
            charge=0,
            spin_multiplicity=2,
        )
        traj.extend(compatible_traj)

        assert len(traj) == 5
        assert traj[-2] == traj[-1]

        # Case of incompatible trajectories
        species, coords, charge, spin = self._get_species_and_coords()

        traj = copy.deepcopy(self.traj)
        incompatible_traj = Trajectory(species=species, coords=coords, charge=charge, spin_multiplicity=spin)
        incompatible_test_success = False
        try:
            traj.extend(incompatible_traj)
        except Exception:
            incompatible_test_success = True

        assert incompatible_test_success

    def test_extend_site_props(self):
        lattice, species, coords = self._get_lattice_species_and_coords()
        num_frames = len(coords)

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

        traj_1 = Trajectory(lattice=lattice, species=species, coords=coords, site_properties=props_1)

        traj_2 = Trajectory(lattice=lattice, species=species, coords=coords, site_properties=props_2)

        traj_3 = Trajectory(lattice=lattice, species=species, coords=coords, site_properties=props_3)

        traj_4 = Trajectory(lattice=lattice, species=species, coords=coords, site_properties=None)

        # const & const (both constant and the same site properties)
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_1)
        expected_site_props = props_1
        assert traj_combined.site_properties == expected_site_props

        # const & const (both constant but different site properties)
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_2)
        expected_site_props = [props_1] * num_frames + [props_2] * num_frames
        assert traj_combined.site_properties == expected_site_props

        # const & changing
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_3)
        expected_site_props = [props_1] * num_frames + props_3
        assert traj_combined.site_properties == expected_site_props

        # const & none
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_4)
        expected_site_props = [props_1] * num_frames + [None] * num_frames
        assert traj_combined.site_properties == expected_site_props

        # changing & const
        traj_combined = copy.deepcopy(traj_3)
        traj_combined.extend(traj_1)
        expected_site_props = props_3 + [props_1] * num_frames
        assert traj_combined.site_properties == expected_site_props

        # changing & changing
        traj_combined = copy.deepcopy(traj_3)
        traj_combined.extend(traj_3)
        expected_site_props = props_3 + props_3
        assert traj_combined.site_properties == expected_site_props

        # changing & none
        traj_combined = copy.deepcopy(traj_3)
        traj_combined.extend(traj_4)
        expected_site_props = props_3 + [None] * num_frames
        assert traj_combined.site_properties == expected_site_props

        # none & const
        traj_combined = copy.deepcopy(traj_4)
        traj_combined.extend(traj_1)
        expected_site_props = [None] * num_frames + [props_1] * num_frames
        assert traj_combined.site_properties == expected_site_props

        # none & changing
        traj_combined = copy.deepcopy(traj_4)
        traj_combined.extend(traj_3)
        expected_site_props = [None] * num_frames + props_3
        assert traj_combined.site_properties == expected_site_props

        # none & none
        traj_combined = copy.deepcopy(traj_4)
        traj_combined.extend(traj_4)
        expected_site_props = None
        assert traj_combined.site_properties == expected_site_props

    def test_extend_frame_props(self):
        lattice, species, coords = self._get_lattice_species_and_coords()

        energy_1 = [-3, -3.9, -4.1]
        energy_2 = [-4.2, -4.25, -4.3]
        pressure_2 = [2, 2.5, 2.5]

        # energy only properties
        props_1 = [{"energy": e} for e in energy_1]
        traj_1 = Trajectory(lattice=lattice, species=species, coords=coords, frame_properties=props_1)

        # energy and pressure properties
        props_2 = [{"energy": e, "pressure": p} for e, p in zip(energy_2, pressure_2)]
        traj_2 = Trajectory(lattice=lattice, species=species, coords=coords, frame_properties=props_2)

        # no properties
        traj_3 = Trajectory(lattice=lattice, species=species, coords=coords, frame_properties=None)

        # test combining two with different properties
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_2)
        expected_props = props_1 + props_2
        assert traj_combined.frame_properties == expected_props

        # test combining two where one has properties and the other does not
        traj_combined = copy.deepcopy(traj_1)
        traj_combined.extend(traj_3)
        expected_props = props_1 + [None] * len(coords)
        assert traj_combined.frame_properties == expected_props

        # test combining two both of which have no properties
        traj_combined = copy.deepcopy(traj_3)
        traj_combined.extend(traj_3)
        assert traj_combined.frame_properties is None

    def test_length(self):
        assert len(self.traj) == len(self.structures)
        assert len(self.traj_mols) == len(self.molecules)

    def test_displacements(self):
        poscar = Poscar.from_file(f"{TEST_FILES_DIR}/POSCAR")
        structures = [poscar.structure]
        displacements = np.zeros((11, *np.shape(structures[-1].frac_coords)))

        for i in range(10):
            displacement = np.random.random_sample(np.shape(structures[-1].frac_coords)) / 20
            new_coords = displacement + structures[-1].frac_coords
            structures.append(Structure(structures[-1].lattice, structures[-1].species, new_coords))
            displacements[i + 1, :, :] = displacement

        traj = Trajectory.from_structures(structures, constant_lattice=True)
        traj.to_displacements()

        assert np.allclose(traj.coords, displacements)

    def test_variable_lattice(self):
        structure = self.structures[0]

        # Generate structures with different lattices
        structures = []
        for _ in range(10):
            new_lattice = np.dot(structure.lattice.matrix, np.diag(1 + np.random.random_sample(3) / 20))
            temp_struct = structure.copy()
            temp_struct.lattice = Lattice(new_lattice)
            structures.append(temp_struct)

        traj = Trajectory.from_structures(structures, constant_lattice=False)

        # Check if lattices were properly stored
        assert all(np.allclose(struct.lattice.matrix, structures[i].lattice.matrix) for i, struct in enumerate(traj))

        # Check if the file is written correctly when lattice is not constant.
        traj.write_Xdatcar(filename="traj_test_XDATCAR")

        # Load trajectory from written xdatcar and compare to original
        written_traj = Trajectory.from_file("traj_test_XDATCAR", constant_lattice=False)
        self._check_traj_equality(traj, written_traj)
        os.remove("traj_test_XDATCAR")

    def test_to_from_dict(self):
        d = self.traj.as_dict()
        traj = Trajectory.from_dict(d)
        assert isinstance(traj, Trajectory)

        d = self.traj_mols.as_dict()
        traj = Trajectory.from_dict(d)
        assert isinstance(traj, Trajectory)

    def test_xdatcar_write(self):
        self.traj.write_Xdatcar(filename="traj_test_XDATCAR")

        # Load trajectory from written xdatcar and compare to original
        written_traj = Trajectory.from_file("traj_test_XDATCAR")
        self._check_traj_equality(self.traj, written_traj)
        os.remove("traj_test_XDATCAR")


if __name__ == "__main__":
    import unittest

    unittest.main()
