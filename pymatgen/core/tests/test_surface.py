#!/usr/bin/python


from __future__ import annotations

import json
import os
import random
import unittest

import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.surface import (
    ReconstructionGenerator,
    Slab,
    SlabGenerator,
    generate_all_slabs,
    get_d,
    get_slab_regions,
    get_symmetrically_distinct_miller_indices,
    get_symmetrically_equivalent_miller_indices,
    miller_index_from_sites,
)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.util.testing import PymatgenTest


def get_path(path_str):
    path = os.path.join(PymatgenTest.TEST_FILES_DIR, "surface_tests", path_str)
    return path


class SlabTest(PymatgenTest):
    def setUp(self):
        zno1 = Structure.from_file(get_path("ZnO-wz.cif"), primitive=False)
        zno55 = SlabGenerator(zno1, [1, 0, 0], 5, 5, lll_reduce=False, center_slab=False).get_slab()

        Ti = Structure(
            Lattice.hexagonal(4.6, 2.82),
            ["Ti", "Ti", "Ti"],
            [
                [0.000000, 0.000000, 0.000000],
                [0.333333, 0.666667, 0.500000],
                [0.666667, 0.333333, 0.500000],
            ],
        )

        Ag_fcc = Structure(
            Lattice.cubic(4.06),
            ["Ag", "Ag", "Ag", "Ag"],
            [
                [0.000000, 0.000000, 0.000000],
                [0.000000, 0.500000, 0.500000],
                [0.500000, 0.000000, 0.500000],
                [0.500000, 0.500000, 0.000000],
            ],
        )

        m = [[3.913449, 0, 0], [0, 3.913449, 0], [0, 0, 5.842644]]
        latt = Lattice(m)
        fcoords = [[0.5, 0, 0.222518], [0, 0.5, 0.777482], [0, 0, 0], [0, 0, 0.5], [0.5, 0.5, 0]]
        non_laue = Structure(latt, ["Nb", "Nb", "N", "N", "N"], fcoords)

        self.ti = Ti
        self.agfcc = Ag_fcc
        self.zno1 = zno1
        self.zno55 = zno55
        self.nonlaue = non_laue
        self.h = Structure(Lattice.cubic(3), ["H"], [[0, 0, 0]])
        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_init(self):
        zno_slab = Slab(
            self.zno55.lattice,
            self.zno55.species,
            self.zno55.frac_coords,
            self.zno55.miller_index,
            self.zno55.oriented_unit_cell,
            0,
            self.zno55.scale_factor,
        )
        m = self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        assert round(abs(zno_slab.surface_area - area), 7) == 0
        assert zno_slab.lattice.parameters == self.zno55.lattice.parameters
        assert zno_slab.oriented_unit_cell.composition == self.zno1.composition
        assert len(zno_slab) == 8

        # check reorient_lattice. get a slab not oriented and check that orientation
        # works even with Cartesian coordinates.
        zno_not_or = SlabGenerator(
            self.zno1,
            [1, 0, 0],
            5,
            5,
            lll_reduce=False,
            center_slab=False,
            reorient_lattice=False,
        ).get_slab()
        zno_slab_cart = Slab(
            zno_not_or.lattice,
            zno_not_or.species,
            zno_not_or.cart_coords,
            zno_not_or.miller_index,
            zno_not_or.oriented_unit_cell,
            0,
            zno_not_or.scale_factor,
            coords_are_cartesian=True,
            reorient_lattice=True,
        )
        self.assertArrayAlmostEqual(zno_slab.frac_coords, zno_slab_cart.frac_coords)
        c = zno_slab_cart.lattice.matrix[2]
        self.assertArrayAlmostEqual([0, 0, np.linalg.norm(c)], c)

    def test_add_adsorbate_atom(self):
        zno_slab = Slab(
            self.zno55.lattice,
            self.zno55.species,
            self.zno55.frac_coords,
            self.zno55.miller_index,
            self.zno55.oriented_unit_cell,
            0,
            self.zno55.scale_factor,
        )
        zno_slab.add_adsorbate_atom([1], "H", 1)

        assert len(zno_slab) == 9
        assert str(zno_slab[8].specie) == "H"
        assert round(abs(zno_slab.get_distance(1, 8) - 1.0), 7) == 0
        assert zno_slab[8].c > zno_slab[0].c
        m = self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        assert round(abs(zno_slab.surface_area - area), 7) == 0
        assert zno_slab.lattice.parameters == self.zno55.lattice.parameters

    def test_get_sorted_structure(self):
        species = [str(site.specie) for site in self.zno55.get_sorted_structure()]
        assert species == ["Zn2+"] * 4 + ["O2-"] * 4

    def test_methods(self):
        # Test various structure methods
        self.zno55.get_primitive_structure()

    def test_as_from_dict(self):
        d = self.zno55.as_dict()
        obj = Slab.from_dict(d)
        assert obj.miller_index == (1, 0, 0)

    def test_dipole_and_is_polar(self):
        self.assertArrayAlmostEqual(self.zno55.dipole, [0, 0, 0])
        assert not self.zno55.is_polar()
        cscl = self.get_structure("CsCl")
        cscl.add_oxidation_state_by_element({"Cs": 1, "Cl": -1})
        slab = SlabGenerator(
            cscl,
            [1, 0, 0],
            5,
            5,
            reorient_lattice=False,
            lll_reduce=False,
            center_slab=False,
        ).get_slab()
        self.assertArrayAlmostEqual(slab.dipole, [-4.209, 0, 0])
        assert slab.is_polar()

    def test_surface_sites_and_symmetry(self):
        # test if surfaces are equivalent by using
        # Laue symmetry and surface site equivalence

        for bool in [True, False]:
            # We will also set the slab to be centered and
            # off centered in order to test the center of mass
            slabgen = SlabGenerator(self.agfcc, (3, 1, 0), 10, 10, center_slab=bool)
            slab = slabgen.get_slabs()[0]
            surf_sites_dict = slab.get_surface_sites()
            assert len(surf_sites_dict["top"]) == len(surf_sites_dict["bottom"])
            total_surf_sites = sum(len(surf_sites_dict[key]) for key in surf_sites_dict)
            assert slab.is_symmetric()
            assert total_surf_sites / 2 == 4

            # Test if the ratio of surface sites per area is
            # constant, ie are the surface energies the same
            r1 = total_surf_sites / (2 * slab.surface_area)
            slabgen = SlabGenerator(self.agfcc, (3, 1, 0), 10, 10, primitive=False)
            slab = slabgen.get_slabs()[0]
            surf_sites_dict = slab.get_surface_sites()
            total_surf_sites = sum(len(surf_sites_dict[key]) for key in surf_sites_dict)
            r2 = total_surf_sites / (2 * slab.surface_area)
            self.assertArrayAlmostEqual(r1, r2)

    def test_symmetrization(self):
        # Restricted to primitive_elemental materials due to the risk of
        # broken stoichiometry. For compound materials, use is_polar()

        # Get all slabs for P6/mmm Ti and Fm-3m Ag up to index of 2

        all_Ti_slabs = generate_all_slabs(
            self.ti,
            2,
            10,
            10,
            bonds=None,
            tol=1e-3,
            max_broken_bonds=0,
            lll_reduce=False,
            center_slab=False,
            primitive=True,
            max_normal_search=2,
            symmetrize=True,
        )

        all_Ag_fcc_slabs = generate_all_slabs(
            self.agfcc,
            2,
            10,
            10,
            bonds=None,
            tol=1e-3,
            max_broken_bonds=0,
            lll_reduce=False,
            center_slab=False,
            primitive=True,
            max_normal_search=2,
            symmetrize=True,
        )

        all_slabs = [all_Ti_slabs, all_Ag_fcc_slabs]

        for slabs in all_slabs:
            asymmetric_count = 0
            symmetric_count = 0

            for slab in slabs:
                sg = SpacegroupAnalyzer(slab)

                # Check if a slab is symmetric
                if not sg.is_laue():
                    asymmetric_count += 1
                else:
                    symmetric_count += 1

            # Check if slabs are all symmetric
            assert asymmetric_count == 0
            assert symmetric_count == len(slabs)

        # Check if we can generate symmetric slabs from bulk with no inversion
        all_non_laue_slabs = generate_all_slabs(self.nonlaue, 1, 15, 15, symmetrize=True)
        assert len(all_non_laue_slabs) > 0

    def test_get_symmetric_sites(self):
        # Check to see if we get an equivalent site on one
        # surface if we add a new site to the other surface

        all_Ti_slabs = generate_all_slabs(
            self.ti,
            2,
            10,
            10,
            bonds=None,
            tol=1e-3,
            max_broken_bonds=0,
            lll_reduce=False,
            center_slab=False,
            primitive=True,
            max_normal_search=2,
            symmetrize=True,
        )

        for slab in all_Ti_slabs:
            sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])
            site = sorted_sites[-1]
            point = np.array(site.frac_coords)
            point[2] = point[2] + 0.1
            point2 = slab.get_symmetric_site(point)
            slab.append("O", point)
            slab.append("O", point2)

            # Check if slab is all symmetric
            sg = SpacegroupAnalyzer(slab)
            assert sg.is_laue()

    def test_oriented_unit_cell(self):
        # Check to see if we get the fully reduced oriented unit
        # cell. This will also ensure that the constrain_latt
        # parameter for get_primitive_structure is working properly

        def surface_area(s):
            m = s.lattice.matrix
            return np.linalg.norm(np.cross(m[0], m[1]))

        all_slabs = generate_all_slabs(self.agfcc, 2, 10, 10, max_normal_search=3)
        for slab in all_slabs:
            ouc = slab.oriented_unit_cell

            assert round(abs(surface_area(slab) - surface_area(ouc)), 7) == 0
            assert len(slab) >= len(ouc)

    def test_get_slab_regions(self):
        # If a slab layer in the slab cell is not completely inside
        # the cell (noncontiguous), check that get_slab_regions will
        # be able to identify where the slab layers are located

        s = self.get_structure("LiFePO4")
        slabgen = SlabGenerator(s, (0, 0, 1), 15, 15)
        slab = slabgen.get_slabs()[0]
        slab.translate_sites([i for i, site in enumerate(slab)], [0, 0, -0.25])
        bottom_c, top_c = [], []
        for site in slab:
            if site.frac_coords[2] < 0.5:
                bottom_c.append(site.frac_coords[2])
            else:
                top_c.append(site.frac_coords[2])
        ranges = get_slab_regions(slab)
        assert tuple(ranges[0]) == (0, max(bottom_c))
        assert tuple(ranges[1]) == (min(top_c), 1)

    def test_as_dict(self):
        slabs = generate_all_slabs(
            self.ti,
            1,
            10,
            10,
            bonds=None,
            tol=1e-3,
            max_broken_bonds=0,
            lll_reduce=False,
            center_slab=False,
            primitive=True,
        )
        slab = slabs[0]
        s = json.dumps(slab.as_dict())
        d = json.loads(s)
        assert slab == Slab.from_dict(d)

        # test initialising with a list scale_factor
        slab = Slab(
            self.zno55.lattice,
            self.zno55.species,
            self.zno55.frac_coords,
            self.zno55.miller_index,
            self.zno55.oriented_unit_cell,
            0,
            self.zno55.scale_factor.tolist(),
        )
        s = json.dumps(slab.as_dict())
        d = json.loads(s)
        assert slab == Slab.from_dict(d)


class SlabGeneratorTest(PymatgenTest):
    def setUp(self):
        lattice = Lattice.cubic(3.010)
        frac_coords = [
            [0.00000, 0.00000, 0.00000],
            [0.00000, 0.50000, 0.50000],
            [0.50000, 0.00000, 0.50000],
            [0.50000, 0.50000, 0.00000],
            [0.50000, 0.00000, 0.00000],
            [0.50000, 0.50000, 0.50000],
            [0.00000, 0.00000, 0.50000],
            [0.00000, 0.50000, 0.00000],
        ]
        species = ["Mg", "Mg", "Mg", "Mg", "O", "O", "O", "O"]
        self.MgO = Structure(lattice, species, frac_coords)
        self.MgO.add_oxidation_state_by_element({"Mg": 2, "O": -6})

        lattice_Dy = Lattice.hexagonal(3.58, 25.61)
        frac_coords_Dy = [
            [0.00000, 0.00000, 0.00000],
            [0.66667, 0.33333, 0.11133],
            [0.00000, 0.00000, 0.222],
            [0.66667, 0.33333, 0.33333],
            [0.33333, 0.66666, 0.44467],
            [0.66667, 0.33333, 0.55533],
            [0.33333, 0.66667, 0.66667],
            [0.00000, 0.00000, 0.778],
            [0.33333, 0.66667, 0.88867],
        ]
        species_Dy = ["Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy", "Dy"]
        self.Dy = Structure(lattice_Dy, species_Dy, frac_coords_Dy)

    def test_get_slab(self):
        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        s = gen.get_slab(0.25)
        assert round(abs(s.lattice.abc[2] - 20.820740000000001), 7) == 0

        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10, max_normal_search=1)
        slab = gen.get_slab()
        assert len(slab) == 6
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10, primitive=False, max_normal_search=1)
        slab_non_prim = gen.get_slab()
        assert len(slab_non_prim) == len(slab) * 4

        # Some randomized testing of cell vectors
        for i in range(1, 231):
            i = random.randint(1, 230)
            sg = SpaceGroup.from_int_number(i)
            if sg.crystal_system == "hexagonal" or (
                sg.crystal_system == "trigonal"
                and (
                    sg.symbol.endswith("H")
                    or sg.int_number
                    in [
                        143,
                        144,
                        145,
                        147,
                        149,
                        150,
                        151,
                        152,
                        153,
                        154,
                        156,
                        157,
                        158,
                        159,
                        162,
                        163,
                        164,
                        165,
                    ]
                )
            ):
                latt = Lattice.hexagonal(5, 10)
            else:
                # Cubic lattice is compatible with all other space groups.
                latt = Lattice.cubic(5)
            s = Structure.from_spacegroup(i, latt, ["H"], [[0, 0, 0]])
            miller = (0, 0, 0)
            while miller == (0, 0, 0):
                miller = (
                    random.randint(0, 6),
                    random.randint(0, 6),
                    random.randint(0, 6),
                )
            gen = SlabGenerator(s, miller, 10, 10)
            a, b, c = gen.oriented_unit_cell.lattice.matrix
            assert round(abs(np.dot(a, gen._normal) - 0), 7) == 0
            assert round(abs(np.dot(b, gen._normal) - 0), 7) == 0

    def test_normal_search(self):
        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        for miller in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
            gen = SlabGenerator(fcc, miller, 10, 10)
            gen_normal = SlabGenerator(fcc, miller, 10, 10, max_normal_search=max(miller))
            slab = gen_normal.get_slab()
            assert round(abs(slab.lattice.alpha - 90), 7) == 0
            assert round(abs(slab.lattice.beta - 90), 7) == 0
            assert len(gen_normal.oriented_unit_cell) >= len(gen.oriented_unit_cell)

        graphite = self.get_structure("Graphite")
        for miller in [(1, 0, 0), (1, 1, 0), (0, 0, 1), (2, 1, 1)]:
            gen = SlabGenerator(graphite, miller, 10, 10)
            gen_normal = SlabGenerator(graphite, miller, 10, 10, max_normal_search=max(miller))
            assert len(gen_normal.oriented_unit_cell) >= len(gen.oriented_unit_cell)

        sc = Structure(
            Lattice.hexagonal(3.32, 5.15),
            ["Sc", "Sc"],
            [[1 / 3, 2 / 3, 0.25], [2 / 3, 1 / 3, 0.75]],
        )
        gen = SlabGenerator(sc, (1, 1, 1), 10, 10, max_normal_search=1)
        assert round(abs(gen.oriented_unit_cell.lattice.angles[1] - 90), 7) == 0

    def test_get_slabs(self):
        gen = SlabGenerator(self.get_structure("CsCl"), [0, 0, 1], 10, 10)

        # Test orthogonality of some internal variables.
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        assert round(abs(np.dot(a, gen._normal) - 0), 7) == 0
        assert round(abs(np.dot(b, gen._normal) - 0), 7) == 0

        assert len(gen.get_slabs()) == 1

        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        assert len(gen.get_slabs()) == 5

        assert len(gen.get_slabs(bonds={("P", "O"): 3})) == 2

        # There are no slabs in LFP that does not break either P-O or Fe-O
        # bonds for a miller index of [0, 0, 1].
        assert len(gen.get_slabs(bonds={("P", "O"): 3, ("Fe", "O"): 3})) == 0

        # If we allow some broken bonds, there are a few slabs.
        assert len(gen.get_slabs(bonds={("P", "O"): 3, ("Fe", "O"): 3}, max_broken_bonds=2)) == 2

        # At this threshold, only the origin and center Li results in
        # clustering. All other sites are non-clustered. So the of
        # slabs is of sites in LiFePO4 unit cell - 2 + 1.
        assert len(gen.get_slabs(tol=1e-4, ftol=1e-4)) == 15

        LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"), primitive=False)
        gen = SlabGenerator(LiCoO2, [0, 0, 1], 10, 10)
        lco = gen.get_slabs(bonds={("Co", "O"): 3})
        assert len(lco) == 1
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        assert round(abs(np.dot(a, gen._normal) - 0), 7) == 0
        assert round(abs(np.dot(b, gen._normal) - 0), 7) == 0

        scc = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Fe"], [[0, 0, 0]])
        gen = SlabGenerator(scc, [0, 0, 1], 10, 10)
        slabs = gen.get_slabs()
        assert len(slabs) == 1
        gen = SlabGenerator(scc, [1, 1, 1], 10, 10, max_normal_search=1)
        slabs = gen.get_slabs()
        assert len(slabs) == 1

        # Test whether using units of hkl planes instead of Angstroms for
        # min_slab_size and min_vac_size will give us the same number of atoms
        natoms = []
        for a in [1, 1.4, 2.5, 3.6]:
            s = Structure.from_spacegroup("Im-3m", Lattice.cubic(a), ["Fe"], [[0, 0, 0]])
            slabgen = SlabGenerator(s, (1, 1, 1), 10, 10, in_unit_planes=True, max_normal_search=2)
            natoms.append(len(slabgen.get_slab()))
        n = natoms[0]
        for i in natoms:
            assert n == i

    def test_triclinic_TeI(self):
        # Test case for a triclinic structure of TeI. Only these three
        # Miller indices are used because it is easier to identify which
        # atoms should be in a surface together. The closeness of the sites
        # in other Miller indices can cause some ambiguity when choosing a
        # higher tolerance.
        numb_slabs = {(0, 0, 1): 5, (0, 1, 0): 3, (1, 0, 0): 7}
        TeI = Structure.from_file(get_path("icsd_TeI.cif"), primitive=False)
        for k, v in numb_slabs.items():
            trclnc_TeI = SlabGenerator(TeI, k, 10, 10)
            TeI_slabs = trclnc_TeI.get_slabs()
            assert v == len(TeI_slabs)

    def test_get_orthogonal_c_slab(self):
        TeI = Structure.from_file(get_path("icsd_TeI.cif"), primitive=False)
        trclnc_TeI = SlabGenerator(TeI, (0, 0, 1), 10, 10)
        TeI_slabs = trclnc_TeI.get_slabs()
        slab = TeI_slabs[0]
        norm_slab = slab.get_orthogonal_c_slab()
        assert round(abs(norm_slab.lattice.angles[0] - 90), 7) == 0
        assert round(abs(norm_slab.lattice.angles[1] - 90), 7) == 0

    def test_get_orthogonal_c_slab_site_props(self):
        TeI = Structure.from_file(get_path("icsd_TeI.cif"), primitive=False)
        trclnc_TeI = SlabGenerator(TeI, (0, 0, 1), 10, 10)
        TeI_slabs = trclnc_TeI.get_slabs()
        slab = TeI_slabs[0]
        # Add site property to slab
        sd_list = [[True, True, True] for site in slab.sites]
        new_sp = slab.site_properties
        new_sp["selective_dynamics"] = sd_list
        slab_with_site_props = slab.copy(site_properties=new_sp)

        # Get orthogonal slab
        norm_slab = slab_with_site_props.get_orthogonal_c_slab()

        # Check if site properties is consistent (or kept)
        assert slab_with_site_props.site_properties == norm_slab.site_properties

    def test_get_tasker2_slabs(self):
        # The uneven distribution of ions on the (111) facets of Halite type
        # slabs are typical examples of Tasker 3 structures. We will test
        # this algo to generate a Tasker 2 structure instead
        slabgen = SlabGenerator(self.MgO, (1, 1, 1), 10, 10, max_normal_search=1)
        # We generate the Tasker 3 structure first
        slab = slabgen.get_slabs()[0]
        assert not slab.is_symmetric()
        assert slab.is_polar()
        # Now to generate the Tasker 2 structure, we must
        # ensure there are enough ions on top to move around
        slab.make_supercell([2, 1, 1])
        slabs = slab.get_tasker2_slabs()
        # Check if our Tasker 2 slab is nonpolar and symmetric
        for slab in slabs:
            assert slab.is_symmetric()
            assert not slab.is_polar()

    def test_nonstoichiometric_symmetrized_slab(self):
        # For the (111) halite slab, sometimes a non-stoichiometric
        # system is preferred over the stoichiometric Tasker 2.
        slabgen = SlabGenerator(self.MgO, (1, 1, 1), 10, 10, max_normal_search=1)
        slabs = slabgen.get_slabs(symmetrize=True)

        # We should end up with two terminations, one with
        # an Mg rich surface and another O rich surface
        assert len(slabs) == 2
        for slab in slabs:
            assert slab.is_symmetric()

        # For a low symmetry primitive_elemental system such as
        # R-3m, there should be some non-symmetric slabs
        # without using non-stoichiometric_symmetrized_slab
        slabs = generate_all_slabs(self.Dy, 1, 30, 30, center_slab=True, symmetrize=True)
        for s in slabs:
            assert s.is_symmetric()
            assert len(s) > len(self.Dy)

    def test_move_to_other_side(self):
        # Tests to see if sites are added to opposite side
        s = self.get_structure("LiFePO4")
        slabgen = SlabGenerator(s, (0, 0, 1), 10, 10, center_slab=True)
        slab = slabgen.get_slab()
        surface_sites = slab.get_surface_sites()

        # check if top sites are moved to the bottom
        top_index = [ss[1] for ss in surface_sites["top"]]
        slab = slabgen.move_to_other_side(slab, top_index)
        all_bottom = [slab[i].frac_coords[2] < slab.center_of_mass[2] for i in top_index]
        assert all(all_bottom)

        # check if bottom sites are moved to the top
        bottom_index = [ss[1] for ss in surface_sites["bottom"]]
        slab = slabgen.move_to_other_side(slab, bottom_index)
        all_top = [slab[i].frac_coords[2] > slab.center_of_mass[2] for i in bottom_index]
        assert all(all_top)

    def test_bonds_broken(self):
        # Querying the Materials Project database for Si
        s = self.get_structure("Si")
        # Conventional unit cell is supplied to ensure miller indices
        # correspond to usual crystallographic definitions
        conv_bulk = SpacegroupAnalyzer(s).get_conventional_standard_structure()
        slabgen = SlabGenerator(conv_bulk, [1, 1, 1], 10, 10, center_slab=True)
        # Setting a generous estimate for max_broken_bonds
        # so that all terminations are generated. These slabs
        # are ordered by ascending number of bonds broken
        # which is assigned to Slab.energy
        slabs = slabgen.get_slabs(bonds={("Si", "Si"): 2.40}, max_broken_bonds=30)
        # Looking at the two slabs generated in VESTA, we
        # expect 2 and 6 bonds broken so we check for this.
        # Number of broken bonds are floats due to primitive
        # flag check and subsequent transformation of slabs.
        assert slabs[0].energy, 2.0
        assert slabs[1].energy, 6.0


class ReconstructionGeneratorTests(PymatgenTest):
    def setUp(self):
        latt = Lattice.cubic(3.51)
        species = ["Ni"]
        coords = [[0, 0, 0]]
        self.Ni = Structure.from_spacegroup("Fm-3m", latt, species, coords)
        latt = Lattice.cubic(2.819000)
        species = ["Fe"]
        coords = [[0, 0, 0]]
        self.Fe = Structure.from_spacegroup("Im-3m", latt, species, coords)
        self.Si = Structure.from_spacegroup("Fd-3m", Lattice.cubic(5.430500), ["Si"], [(0, 0, 0.5)])

        with open(
            os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                "..",
                "reconstructions_archive.json",
            )
        ) as data_file:
            self.rec_archive = json.load(data_file)

    def test_build_slab(self):
        # First lets test a reconstruction where we only remove atoms
        recon = ReconstructionGenerator(self.Ni, 10, 10, "fcc_110_missing_row_1x2")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        assert recon_slab.reconstruction
        assert len(slab) == len(recon_slab) + 2
        assert recon_slab.is_symmetric()

        # Test if the ouc corresponds to the reconstructed slab
        recon_ouc = recon_slab.oriented_unit_cell
        ouc = slab.oriented_unit_cell
        assert ouc.lattice.b * 2 == recon_ouc.lattice.b
        assert len(ouc) * 2 == len(recon_ouc)

        # Test a reconstruction where we simply add atoms
        recon = ReconstructionGenerator(self.Ni, 10, 10, "fcc_111_adatom_t_1x1")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        assert len(slab) == len(recon_slab) - 2
        assert recon_slab.is_symmetric()

        # If a slab references another slab,
        # make sure it is properly generated
        recon = ReconstructionGenerator(self.Ni, 10, 10, "fcc_111_adatom_ft_1x1")
        slab = recon.build_slabs()[0]
        assert slab.is_symmetric

        # Test a reconstruction where it works on a specific
        # termination (Fd-3m (111))
        recon = ReconstructionGenerator(self.Si, 10, 10, "diamond_111_1x2")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        assert len(slab) == len(recon_slab) - 8
        assert recon_slab.is_symmetric()

        # Test a reconstruction where terminations give
        # different reconstructions with a non-primitive_elemental system

    def test_get_d(self):
        # Ensure that regardless of the size of the vacuum or slab
        # layer, the spacing between atomic layers should be the same

        recon = ReconstructionGenerator(self.Si, 10, 10, "diamond_100_2x1")

        recon2 = ReconstructionGenerator(self.Si, 20, 10, "diamond_100_2x1")
        s1 = recon.get_unreconstructed_slabs()[0]
        s2 = recon2.get_unreconstructed_slabs()[0]
        assert round(abs(get_d(s1) - get_d(s2)), 7) == 0

    @unittest.skip("This test relies on neighbor orders and is hard coded. Disable temporarily")
    def test_previous_reconstructions(self):
        # Test to see if we generated all reconstruction types correctly and nothing changes

        m = StructureMatcher()
        for n in self.rec_archive:
            if "base_reconstruction" in self.rec_archive[n]:
                arch = self.rec_archive[self.rec_archive[n]["base_reconstruction"]]
                sg = arch["spacegroup"]["symbol"]
            else:
                sg = self.rec_archive[n]["spacegroup"]["symbol"]
            if sg == "Fm-3m":
                rec = ReconstructionGenerator(self.Ni, 20, 20, n)
                el = self.Ni[0].species_string
            elif sg == "Im-3m":
                rec = ReconstructionGenerator(self.Fe, 20, 20, n)
                el = self.Fe[0].species_string
            elif sg == "Fd-3m":
                rec = ReconstructionGenerator(self.Si, 20, 20, n)
                el = self.Si[0].species_string

            slabs = rec.build_slabs()
            s = Structure.from_file(get_path(os.path.join("reconstructions", el + "_" + n + ".cif")))
            assert any(len(m.group_structures([s, slab])) == 1 for slab in slabs)


class MillerIndexFinderTests(PymatgenTest):
    def setUp(self):
        self.cscl = Structure.from_spacegroup("Pm-3m", Lattice.cubic(4.2), ["Cs", "Cl"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.Fe = Structure.from_spacegroup("Im-3m", Lattice.cubic(2.82), ["Fe"], [[0, 0, 0]])
        mglatt = Lattice.from_parameters(3.2, 3.2, 5.13, 90, 90, 120)
        self.Mg = Structure(mglatt, ["Mg", "Mg"], [[1 / 3, 2 / 3, 1 / 4], [2 / 3, 1 / 3, 3 / 4]])
        self.lifepo4 = self.get_structure("LiFePO4")
        self.tei = Structure.from_file(get_path("icsd_TeI.cif"), primitive=False)
        self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"), primitive=False)

        self.p1 = Structure(
            Lattice.from_parameters(3, 4, 5, 31, 43, 50),
            ["H", "He"],
            [[0, 0, 0], [0.1, 0.2, 0.3]],
        )
        self.graphite = self.get_structure("Graphite")
        self.trigBi = Structure(
            Lattice.from_parameters(3, 3, 10, 90, 90, 120),
            ["Bi", "Bi", "Bi", "Bi", "Bi", "Bi"],
            [
                [0.3333, 0.6666, 0.39945113],
                [0.0000, 0.0000, 0.26721554],
                [0.0000, 0.0000, 0.73278446],
                [0.6666, 0.3333, 0.60054887],
                [0.6666, 0.3333, 0.06611779],
                [0.3333, 0.6666, 0.93388221],
            ],
        )

    def test_get_symmetrically_distinct_miller_indices(self):
        # Tests to see if the function obtains the known number of unique slabs

        indices = get_symmetrically_distinct_miller_indices(self.cscl, 1)
        assert len(indices) == 3
        indices = get_symmetrically_distinct_miller_indices(self.cscl, 2)
        assert len(indices) == 6

        assert len(get_symmetrically_distinct_miller_indices(self.lifepo4, 1)) == 7

        # The TeI P-1 structure should have 13 unique millers (only inversion
        # symmetry eliminates pairs)
        indices = get_symmetrically_distinct_miller_indices(self.tei, 1)
        assert len(indices) == 13

        # P1 and P-1 should have the same # of miller indices since surfaces
        # always have inversion symmetry.
        indices = get_symmetrically_distinct_miller_indices(self.p1, 1)
        assert len(indices) == 13

        indices = get_symmetrically_distinct_miller_indices(self.graphite, 2)
        assert len(indices) == 12

        # Now try a trigonal system.
        indices = get_symmetrically_distinct_miller_indices(self.trigBi, 2, return_hkil=True)
        assert len(indices) == 17
        assert all(len(hkl) == 4 for hkl in indices)

    def test_get_symmetrically_equivalent_miller_indices(self):
        # Tests to see if the function obtains all equivalent hkl for cubic (100)
        indices001 = [
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (0, 0, -1),
            (0, -1, 0),
            (-1, 0, 0),
        ]
        indices = get_symmetrically_equivalent_miller_indices(self.cscl, (1, 0, 0))
        assert all(hkl in indices for hkl in indices001)

        # Tests to see if it captures expanded Miller indices in the family e.g. (001) == (002)
        hcp_indices_100 = get_symmetrically_equivalent_miller_indices(self.Mg, (1, 0, 0))
        hcp_indices_200 = get_symmetrically_equivalent_miller_indices(self.Mg, (2, 0, 0))
        assert len(hcp_indices_100) * 2 == len(hcp_indices_200)
        assert len(hcp_indices_100) == 6
        assert all(len(hkl) == 4 for hkl in hcp_indices_100)

    def test_generate_all_slabs(self):
        slabs = generate_all_slabs(self.cscl, 1, 10, 10)
        # Only three possible slabs, one each in (100), (110) and (111).
        assert len(slabs) == 3

        # make sure it generates reconstructions
        slabs = generate_all_slabs(self.Fe, 1, 10, 10, include_reconstructions=True)

        # Four possible slabs, (100), (110), (111) and the zigzag (100).
        assert len(slabs) == 4

        slabs = generate_all_slabs(self.cscl, 1, 10, 10, bonds={("Cs", "Cl"): 4})
        # No slabs if we don't allow broken Cs-Cl
        assert len(slabs) == 0

        slabs = generate_all_slabs(self.cscl, 1, 10, 10, bonds={("Cs", "Cl"): 4}, max_broken_bonds=100)
        assert len(slabs) == 3

        slabs2 = generate_all_slabs(self.lifepo4, 1, 10, 10, bonds={("P", "O"): 3, ("Fe", "O"): 3})
        assert len(slabs2) == 0

        # There should be only one possible stable surfaces, all of which are
        # in the (001) oriented unit cell
        slabs3 = generate_all_slabs(self.LiCoO2, 1, 10, 10, bonds={("Co", "O"): 3})
        assert len(slabs3) == 1
        mill = (0, 0, 1)
        for s in slabs3:
            assert s.miller_index == mill

        slabs1 = generate_all_slabs(self.lifepo4, 1, 10, 10, tol=0.1, bonds={("P", "O"): 3})
        assert len(slabs1) == 4

        # Now we test this out for repair_broken_bonds()
        slabs1_repair = generate_all_slabs(self.lifepo4, 1, 10, 10, tol=0.1, bonds={("P", "O"): 3}, repair=True)
        assert len(slabs1_repair) > len(slabs1)

        # Lets see if there are no broken PO4 polyhedrons
        miller_list = get_symmetrically_distinct_miller_indices(self.lifepo4, 1)
        all_miller_list = []
        for slab in slabs1_repair:
            hkl = tuple(slab.miller_index)
            if hkl not in all_miller_list:
                all_miller_list.append(hkl)
            broken = []
            for site in slab:
                if site.species_string == "P":
                    neighbors = slab.get_neighbors(site, 3)
                    cn = 0
                    for nn in neighbors:
                        cn += 1 if nn[0].species_string == "O" else 0
                    broken.append(cn != 4)
            assert not any(broken)

        # check if we were able to produce at least one
        # termination for each distinct Miller _index
        assert len(miller_list) == len(all_miller_list)

    def test_miller_index_from_sites(self):
        """Test surface miller index convenience function"""
        # test on a cubic system
        m = Lattice.cubic(1)
        s1 = np.array([0.5, -1.5, 3])
        s2 = np.array([0.5, 3.0, -1.5])
        s3 = np.array([2.5, 1.5, -4.0])
        assert miller_index_from_sites(m, [s1, s2, s3]) == (2, 1, 1)

        # test casting from matrix to Lattice
        m = [[2.319, -4.01662582, 0.0], [2.319, 4.01662582, 0.0], [0.0, 0.0, 7.252]]

        s1 = np.array([2.319, 1.33887527, 6.3455])
        s2 = np.array([1.1595, 0.66943764, 4.5325])
        s3 = np.array([1.1595, 0.66943764, 0.9065])
        hkl = miller_index_from_sites(m, [s1, s2, s3])
        assert hkl == (2, -1, 0)


if __name__ == "__main__":
    unittest.main()
