# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import unittest

import numpy as np

from pymatgen.core.interface import Interface
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest


class InterfaceTest(PymatgenTest):
    def setUp(self):
        self.interface: Interface = self.get_structure("Si_SiO2_Interface")

    def test_basic_props(self):
        interface = self.interface
        assert isinstance(interface, Interface)

        assert len(interface.substrate_indices) == 14
        assert len(interface.film_indices) == 36
        assert len(interface.film_sites) == len(interface.film_indices)
        assert len(interface.substrate_sites) == len(interface.substrate_indices)
        assert interface.gap == 2.0
        assert np.allclose(interface.in_plane_offset, [0, 0])
        assert interface.vacuum_over_film == 20.0
        assert interface.film_termination == "O2_P6/mmm_4"
        assert interface.substrate_termination == "Si_P6/mmm_7"
        assert interface.film_layers == 6
        assert interface.substrate_layers == 2

        iface_dict = interface.as_dict()
        for k in [
            "lattice",
            "sites",
            "in_plane_offset",
            "gap",
            "vacuum_over_film",
            "interface_properties",
        ]:
            assert k in iface_dict
        assert isinstance(interface.from_dict(iface_dict), Interface)

    def test_gap_setter(self):
        interface = self.interface

        assert np.allclose(interface.gap, 2.0)

        max_sub_c = np.max(np.array([s.frac_coords for s in interface.substrate])[:, 2])
        min_film_c = np.min(np.array([f.frac_coords for f in interface.film])[:, 2])
        gap = (min_film_c - max_sub_c) * interface.lattice.c
        assert np.allclose(interface.gap, gap)

        interface.gap += 1

        assert np.allclose(interface.gap, 3.0)

        max_sub_c = np.max(np.array([s.frac_coords for s in interface.substrate])[:, 2])
        min_film_c = np.min(np.array([f.frac_coords for f in interface.film])[:, 2])
        gap = (min_film_c - max_sub_c) * interface.lattice.c
        assert np.allclose(interface.gap, gap)

    def test_in_plane_offset_setter(self):

        interface = self.interface
        init_coords = np.array(self.interface.frac_coords)
        interface.in_plane_offset = np.array([0.2, 0.2])

        assert np.allclose(interface.in_plane_offset, np.array([0.2, 0.2]))

        test_coords = np.array(init_coords)
        for i in interface.film_indices:
            test_coords[i] += [0.2, 0.2, 0]
        assert np.allclose(np.mod(test_coords, 1.0), np.mod(interface.frac_coords, 1.0))

    def test_vacuum_over_film_setter(self):
        interface = self.interface
        init_coords = self.interface.cart_coords

        assert interface.vacuum_over_film == 20

        interface.vacuum_over_film += 10

        assert interface.vacuum_over_film == 30
        assert np.allclose(init_coords, interface.cart_coords)

    def test_get_shifts_based_on_adsorbate_sites(self):
        # Only testing two tolerances as there appears to be significant numerical noise in this method
        assert len(self.interface.get_shifts_based_on_adsorbate_sites()) == 42
        assert len(self.interface.get_shifts_based_on_adsorbate_sites(tolerance=20.0)) == 1

    def test_from_slabs(self):
        si_struct = self.get_structure("Si")
        sio2_struct = self.get_structure("SiO2")
        si_conventional = SpacegroupAnalyzer(si_struct).get_conventional_standard_structure()
        sio2_conventional = SpacegroupAnalyzer(sio2_struct).get_conventional_standard_structure()

        si_slab = SlabGenerator(si_conventional, (1, 1, 1), 5, 10, reorient_lattice=True).get_slab()
        sio2_slab = SlabGenerator(sio2_conventional, (1, 0, 0), 5, 10, reorient_lattice=True).get_slab()

        interface = Interface.from_slabs(film_slab=si_slab, substrate_slab=sio2_slab)
        assert isinstance(interface, Interface)


if __name__ == "__main__":
    unittest.main()
