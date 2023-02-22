from __future__ import annotations

import json
import os
import unittest

from pytest import approx

from pymatgen.analysis.wulff import WulffShape
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord import in_coord_list
from pymatgen.util.testing import PymatgenTest

__author__ = "Zihan Xu, Richard Tran, Balachandran Radhakrishnan"
__copyright__ = "Copyright 2013, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Zihan Xu"
__email__ = "zix009@eng.ucsd.edu"
__date__ = "May 05 2016"


class WulffShapeTest(PymatgenTest):
    def setUp(self):
        module_dir = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(module_dir, "surface_samples.json")) as data_file:
            surface_properties = json.load(data_file)

        surface_energies, miller_indices = {}, {}
        for mpid in surface_properties:
            e_surf_list, miller_list = [], []
            for surface in surface_properties[mpid]["surfaces"]:
                e_surf_list.append(surface["surface_energy"])
                miller_list.append(surface["miller_index"])
            surface_energies[mpid] = e_surf_list
            miller_indices[mpid] = miller_list

        # In the case of a high anisotropy material
        # Nb: mp-8636
        latt_Nb = Lattice.cubic(2.992)
        # In the case of an fcc material
        # Ir: mp-101
        latt_Ir = Lattice.cubic(3.8312)
        # In the case of a hcp material
        # Ti: mp-72
        latt_Ti = Lattice.hexagonal(4.6000, 2.8200)
        self.ucell_Nb = Structure(
            latt_Nb,
            ["Nb", "Nb", "Nb", "Nb"],
            [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
        )
        self.wulff_Nb = WulffShape(latt_Nb, miller_indices["mp-8636"], surface_energies["mp-8636"])

        self.ucell_Ir = Structure(
            latt_Nb,
            ["Ir", "Ir", "Ir", "Ir"],
            [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
        )
        self.wulff_Ir = WulffShape(latt_Ir, miller_indices["mp-101"], surface_energies["mp-101"])

        self.ucell_Ti = Structure(
            latt_Ti,
            ["Ti", "Ti", "Ti"],
            [[0, 0, 0], [0.333333, 0.666667, 0.5], [0.666667, 0.333333, 0.5]],
        )
        self.wulff_Ti = WulffShape(latt_Ti, miller_indices["mp-72"], surface_energies["mp-72"])
        self.cube = WulffShape(Lattice.cubic(1), [(1, 0, 0)], [1])
        self.hex_prism = WulffShape(Lattice.hexagonal(2.63, 5.21), [(0, 0, 1), (1, 0, 0)], [0.35, 0.53])

        self.surface_properties = surface_properties

    @unittest.skipIf("DISPLAY" not in os.environ, "Need display")
    def test_get_plot(self):
        # Basic test, not really a unittest.
        self.wulff_Ti.get_plot()
        self.wulff_Nb.get_plot()
        self.wulff_Ir.get_plot()

    @unittest.skipIf("DISPLAY" not in os.environ, "Need display")
    def test_get_plotly(self):
        # Basic test, not really a unittest.
        self.wulff_Ti.get_plotly()
        self.wulff_Nb.get_plotly()
        self.wulff_Ir.get_plotly()

    def symm_check(self, ucell, wulff_vertices):
        """
        # Checks if the point group of the Wulff shape matches
        # the point group of its conventional unit cell

        Args:
            ucell (str): Unit cell that the Wulff shape is based on.
            wulff_vertices (list): List of all vertices on the Wulff
                shape. Use wulff.wulff_pt_list to obtain the list
                (see wulff_generator.py).

        return (bool)
        """
        space_group_analyzer = SpacegroupAnalyzer(ucell)
        symm_ops = space_group_analyzer.get_point_group_operations(cartesian=True)
        for point in wulff_vertices:
            for op in symm_ops:
                symm_point = op.operate(point)
                if in_coord_list(wulff_vertices, symm_point):
                    continue
                else:
                    return False
        return True

    def consistency_tests(self):
        # For a set of given values, these tests will
        # ensure that the general result given by the
        # algorithm does not change as the code is edited

        # For fcc Ir, make sure the (111) direction
        # is the most dominant facet on the Wulff shape

        fractional_areas = self.wulff_Ir.area_fraction_dict
        miller_list = list(fractional_areas)
        area_list = [fractional_areas[hkl] for hkl in fractional_areas]
        assert miller_list[area_list.index(max(area_list))] == (1, 1, 1)

        # Overall weighted surface energy of fcc Nb should be
        # equal to the energy of the (310) surface, ie. fcc Nb
        # is anisotropic, the (310) surface is so low in energy,
        # its the only facet that exists in the Wulff shape

        Nb_area_fraction_dict = self.wulff_Nb.area_fraction_dict
        for hkl in Nb_area_fraction_dict:
            if hkl == (3, 1, 0):
                assert Nb_area_fraction_dict[hkl] == 1
            else:
                assert Nb_area_fraction_dict[hkl] == 0

        assert self.wulff_Nb.miller_energy_dict[(3, 1, 0)] == self.wulff_Nb.weighted_surface_energy

    def symmetry_test(self):
        # Maintains that all wulff shapes have the same point
        # groups as the conventional unit cell they were
        # derived from. This test should pass for all subsequent
        # updates of the surface_properties collection

        check_symmetry_Nb = self.symm_check(self.ucell_Nb, self.wulff_Nb.wulff_pt_list)
        check_symmetry_Ir = self.symm_check(self.ucell_Ir, self.wulff_Ir.wulff_pt_list)
        check_symmetry_Ti = self.symm_check(self.ucell_Ti, self.wulff_Ti.wulff_pt_list)
        assert check_symmetry_Nb
        assert check_symmetry_Ir
        assert check_symmetry_Ti

    def test_get_azimuth_elev(self):
        # Test out the viewing of the Wulff shape from Miller indices.
        azim, elev = self.wulff_Ir._get_azimuth_elev((0, 0, 1))
        assert azim == 0
        assert elev == 90
        azim, elev = self.wulff_Ir._get_azimuth_elev((1, 1, 1))
        assert azim == approx(45)

    def test_properties(self):
        # Simple test to check if the values of some
        # properties are consistent with what we already have

        wulff_shapes = {
            "mp-8636": self.wulff_Nb,
            "mp-72": self.wulff_Ti,
            "mp-101": self.wulff_Ir,
        }
        for mp_id, wulff in wulff_shapes.items():
            properties = self.surface_properties[mp_id]
            assert round(wulff.weighted_surface_energy, 3) == round(properties["weighted_surface_energy"], 3)
            assert round(wulff.shape_factor, 3) == round(properties["shape_factor"], 3)
            assert round(wulff.anisotropy, 3) == round(properties["surface_anisotropy"], 3)

    def test_corner_and_edges(self):
        # Test if it is returning the correct number of corner and edges
        self.assertArrayEqual(self.cube.tot_corner_sites, 8)
        self.assertArrayEqual(self.cube.tot_edges, 12)
        self.assertArrayEqual(self.hex_prism.tot_corner_sites, 12)
        self.assertArrayEqual(self.hex_prism.tot_edges, 18)


if __name__ == "__main__":
    unittest.main()
