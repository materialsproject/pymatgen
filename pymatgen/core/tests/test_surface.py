#!/usr/bin/python


import unittest
import os
import random
import json

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import Slab, SlabGenerator, generate_all_slabs, \
    get_symmetrically_distinct_miller_indices, ReconstructionGenerator, \
    miller_index_from_sites, get_d
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.structure_matcher import StructureMatcher


def get_path(path_str):
    cwd = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(cwd, "..", "..", "..", "test_files", "surface_tests",
                        path_str)
    return path


class SlabTest(PymatgenTest):

    def setUp(self):
        zno1 = Structure.from_file(get_path("ZnO-wz.cif"), primitive=False)
        zno55 = SlabGenerator(zno1, [1, 0, 0], 5, 5, lll_reduce=False,
                              center_slab=False).get_slab()

        Ti = Structure(Lattice.hexagonal(4.6, 2.82), ["Ti", "Ti", "Ti"],
                       [[0.000000, 0.000000, 0.000000],
                       [0.333333, 0.666667, 0.500000],
                       [0.666667, 0.333333, 0.500000]])

        Ag_fcc = Structure(Lattice.cubic(4.06), ["Ag", "Ag", "Ag", "Ag"],
                           [[0.000000, 0.000000, 0.000000],
                           [0.000000, 0.500000, 0.500000],
                           [0.500000, 0.000000, 0.500000],
                           [0.500000, 0.500000, 0.000000]])

        self.ti = Ti
        self.agfcc = Ag_fcc
        self.zno1 = zno1
        self.zno55 = zno55
        self.h = Structure(Lattice.cubic(3), ["H"],
                            [[0, 0, 0]])
        self.libcc = Structure(Lattice.cubic(3.51004), ["Li", "Li"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

    def test_init(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species,
                        self.zno55.frac_coords,
                        self.zno55.miller_index,
                        self.zno55.oriented_unit_cell,
                        0, self.zno55.scale_factor)
        m =self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)
        self.assertEqual(zno_slab.oriented_unit_cell.composition,
                         self.zno1.composition)
        self.assertEqual(len(zno_slab), 8)

    def test_add_adsorbate_atom(self):
        zno_slab = Slab(self.zno55.lattice, self.zno55.species,
                        self.zno55.frac_coords,
                        self.zno55.miller_index,
                        self.zno55.oriented_unit_cell,
                        0, self.zno55.scale_factor)
        zno_slab.add_adsorbate_atom([1], 'H', 1)

        self.assertEqual(len(zno_slab), 9)
        self.assertEqual(str(zno_slab[8].specie), 'H')
        self.assertAlmostEqual(zno_slab.get_distance(1, 8), 1.0)
        self.assertTrue(zno_slab[8].c > zno_slab[0].c)
        m = self.zno55.lattice.matrix
        area = np.linalg.norm(np.cross(m[0], m[1]))
        self.assertAlmostEqual(zno_slab.surface_area, area)
        self.assertEqual(zno_slab.lattice.lengths_and_angles,
                         self.zno55.lattice.lengths_and_angles)

    def test_get_sorted_structure(self):
        species = [str(site.specie) for site in
                   self.zno55.get_sorted_structure()]
        self.assertEqual(species, ["Zn2+"] * 4 + ["O2-"] * 4)

    def test_methods(self):
        #Test various structure methods
        self.zno55.get_primitive_structure()

    def test_as_from_dict(self):
        d = self.zno55.as_dict()
        obj = Slab.from_dict(d)
        self.assertEqual(obj.miller_index, (1, 0, 0))

    def test_dipole_and_is_polar(self):
        self.assertArrayAlmostEqual(self.zno55.dipole, [0, 0, 0])
        self.assertFalse(self.zno55.is_polar())
        cscl = self.get_structure("CsCl")
        cscl.add_oxidation_state_by_element({"Cs": 1, "Cl": -1})
        slab = SlabGenerator(cscl, [1, 0, 0], 5, 5, reorient_lattice=False,
                             lll_reduce=False, center_slab=False).get_slab()
        self.assertArrayAlmostEqual(slab.dipole, [-4.209, 0, 0])
        self.assertTrue(slab.is_polar())

    def test_surface_sites_and_symmetry(self):
        # test if surfaces are equivalent by using
        # Laue symmetry and surface site equivalence

        for bool in [True, False]:
            # We will also set the slab to be centered and
            # off centered in order to test the center of mass
            slabgen = SlabGenerator(self.agfcc, (3, 1, 0), 10, 10, center_slab=bool)
            slab = slabgen.get_slabs()[0]
            surf_sites_dict = slab.get_surface_sites()
            self.assertEqual(len(surf_sites_dict["top"]), len(surf_sites_dict["bottom"]))
            total_surf_sites = sum([len(surf_sites_dict[key])
                                    for key in surf_sites_dict.keys()])
            self.assertTrue(slab.is_symmetric())
            self.assertEqual(total_surf_sites/2, 4)
            self.assertTrue(slab.have_equivalent_surfaces())

            # Test if the ratio of surface sites per area is
            # constant, ie are the surface energies the same
            r1 = total_surf_sites/(2*slab.surface_area)
            slabgen = SlabGenerator(self.agfcc, (3, 1, 0), 10, 10, primitive=False)
            slab = slabgen.get_slabs()[0]
            surf_sites_dict = slab.get_surface_sites()
            total_surf_sites = sum([len(surf_sites_dict[key])
                                    for key in surf_sites_dict.keys()])
            r2 = total_surf_sites/(2*slab.surface_area)
            self.assertArrayEqual(r1, r2)

    def test_symmetrization(self):

        # Restricted to primitive_elemental materials due to the risk of
        # broken stoichiometry. For compound materials, use is_polar()

        # Get all slabs for P6/mmm Ti and Fm-3m Ag up to index of 2

        all_Ti_slabs = generate_all_slabs(self.ti, 2, 10, 10, bonds=None,
                                          tol=1e-3, max_broken_bonds=0,
                                          lll_reduce=False, center_slab=False,
                                          primitive=True, max_normal_search=2,
                                          symmetrize=True)

        all_Ag_fcc_slabs = generate_all_slabs(self.agfcc, 2, 10, 10, bonds=None,
                                              tol=1e-3, max_broken_bonds=0,
                                              lll_reduce=False, center_slab=False,
                                              primitive=True, max_normal_search=2,
                                              symmetrize=True)

        all_slabs = [all_Ti_slabs, all_Ag_fcc_slabs]

        for i, slabs in enumerate(all_slabs):

            assymetric_count = 0
            symmetric_count = 0

            for i, slab in enumerate(slabs):
                sg = SpacegroupAnalyzer(slab)

                # Check if a slab is symmetric
                if not sg.is_laue():
                    assymetric_count += 1
                else:
                    symmetric_count += 1

            # Check if slabs are all symmetric
            self.assertEqual(assymetric_count, 0)
            self.assertEqual(symmetric_count, len(slabs))

    def test_get_symmetric_sites(self):

        # Check to see if we get an equivalent site on one
        # surface if we add a new site to the other surface

        all_Ti_slabs = generate_all_slabs(self.ti, 2, 10, 10, bonds=None,
                                          tol=1e-3, max_broken_bonds=0,
                                          lll_reduce=False, center_slab=False,
                                          primitive=True, max_normal_search=2,
                                          symmetrize=True)

        for slab in all_Ti_slabs:
            sorted_sites = sorted(slab, key=lambda site: site.frac_coords[2])
            site = sorted_sites[-1]
            point = site.frac_coords
            point[2] = point[2]+0.1
            point2 = slab.get_symmetric_site(point)
            slab.append("O", point)
            slab.append("O", point2)

            # Check if slab is all symmetric
            sg = SpacegroupAnalyzer(slab)
            self.assertTrue(sg.is_laue())

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
            self.assertAlmostEqual(surface_area(slab), surface_area(ouc))
            self.assertGreaterEqual(len(slab), len(ouc))


class SlabGeneratorTest(PymatgenTest):

    def setUp(self):

        lattice = Lattice.cubic(3.010)
        frac_coords = [[0.00000, 0.00000, 0.00000],
                       [0.00000, 0.50000, 0.50000],
                       [0.50000, 0.00000, 0.50000],
                       [0.50000, 0.50000, 0.00000],
                       [0.50000, 0.00000, 0.00000],
                       [0.50000, 0.50000, 0.50000],
                       [0.00000, 0.00000, 0.50000],
                       [0.00000, 0.50000, 0.00000]]
        species = ['Mg', 'Mg', 'Mg', 'Mg', 'O', 'O', 'O', 'O']
        self.MgO = Structure(lattice, species, frac_coords)
        self.MgO.add_oxidation_state_by_element({"Mg": 2, "O": -6})

        lattice_Dy = Lattice.hexagonal(3.58, 25.61)
        frac_coords_Dy = [[0.00000, 0.00000, 0.00000],
                       [0.66667, 0.33333, 0.11133],
                       [0.00000, 0.00000, 0.222],
                       [0.66667, 0.33333, 0.33333],
                       [0.33333, 0.66666, 0.44467],
                       [0.66667, 0.33333, 0.55533],
                       [0.33333, 0.66667, 0.66667],
                       [0.00000, 0.00000, 0.778],
                       [0.33333, 0.66667, 0.88867]]
        species_Dy = ['Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy', 'Dy']
        self.Dy = Structure(lattice_Dy, species_Dy, frac_coords_Dy)

    def test_get_slab(self):
        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        s = gen.get_slab(0.25)
        self.assertAlmostEqual(s.lattice.abc[2], 20.820740000000001)

        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10)
        slab = gen.get_slab()
        gen = SlabGenerator(fcc, [1, 1, 1], 10, 10, primitive=False)
        slab_non_prim = gen.get_slab()
        self.assertEqual(len(slab), 6)
        self.assertEqual(len(slab_non_prim), len(slab) * 4)

        # Some randomized testing of cell vectors
        for i in range(1, 231):
            i = random.randint(1, 230)
            sg = SpaceGroup.from_int_number(i)
            if sg.crystal_system == "hexagonal" or (sg.crystal_system == \
                    "trigonal" and (sg.symbol.endswith("H") or
                    sg.int_number in [143, 144, 145, 147, 149, 150, 151, 152,
                                      153, 154, 156, 157, 158, 159, 162, 163,
                                      164, 165])):
                latt = Lattice.hexagonal(5, 10)
            else:
                # Cubic lattice is compatible with all other space groups.
                latt = Lattice.cubic(5)
            s = Structure.from_spacegroup(i, latt, ["H"], [[0, 0, 0]])
            miller = (0, 0, 0)
            while miller == (0, 0, 0):
                miller = (random.randint(0, 6), random.randint(0, 6),
                          random.randint(0, 6))
            gen = SlabGenerator(s, miller, 10, 10)
            a, b, c = gen.oriented_unit_cell.lattice.matrix
            self.assertAlmostEqual(np.dot(a, gen._normal), 0)
            self.assertAlmostEqual(np.dot(b, gen._normal), 0)

    def test_normal_search(self):
        fcc = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        for miller in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
            gen = SlabGenerator(fcc, miller, 10, 10)
            gen_normal = SlabGenerator(fcc, miller, 10, 10,
                                       max_normal_search=max(miller))
            slab = gen_normal.get_slab()
            self.assertAlmostEqual(slab.lattice.alpha, 90)
            self.assertAlmostEqual(slab.lattice.beta, 90)
            self.assertGreaterEqual(len(gen_normal.oriented_unit_cell),
                                    len(gen.oriented_unit_cell))

        graphite = self.get_structure("Graphite")
        for miller in [(1, 0, 0), (1, 1, 0), (0, 0, 1), (2, 1, 1)]:
            gen = SlabGenerator(graphite, miller, 10, 10)
            gen_normal = SlabGenerator(graphite, miller, 10, 10,
                                       max_normal_search=max(miller))
            self.assertGreaterEqual(len(gen_normal.oriented_unit_cell),
                                    len(gen.oriented_unit_cell))

        sc = Structure(Lattice.hexagonal(3.32, 5.15), ["Sc", "Sc"],
                       [[1/3, 2/3, 0.25], [2/3, 1/3, 0.75]])
        gen = SlabGenerator(sc, (1, 1, 1), 10, 10, max_normal_search=1)
        self.assertAlmostEqual(gen.oriented_unit_cell.lattice.angles[1], 90)

    def test_get_slabs(self):
        gen = SlabGenerator(self.get_structure("CsCl"), [0, 0, 1], 10, 10)

        #Test orthogonality of some internal variables.
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        self.assertAlmostEqual(np.dot(a, gen._normal), 0)
        self.assertAlmostEqual(np.dot(b, gen._normal), 0)

        self.assertEqual(len(gen.get_slabs()), 1)

        s = self.get_structure("LiFePO4")
        gen = SlabGenerator(s, [0, 0, 1], 10, 10)
        self.assertEqual(len(gen.get_slabs()), 5)

        self.assertEqual(len(gen.get_slabs(bonds={("P", "O"): 3})), 2)

        # There are no slabs in LFP that does not break either P-O or Fe-O
        # bonds for a miller index of [0, 0, 1].
        self.assertEqual(len(gen.get_slabs(
            bonds={("P", "O"): 3, ("Fe", "O"): 3})), 0)

        #If we allow some broken bonds, there are a few slabs.
        self.assertEqual(len(gen.get_slabs(
            bonds={("P", "O"): 3, ("Fe", "O"): 3},
            max_broken_bonds=2)), 2)

        # At this threshold, only the origin and center Li results in
        # clustering. All other sites are non-clustered. So the of
        # slabs is of sites in LiFePO4 unit cell - 2 + 1.
        self.assertEqual(len(gen.get_slabs(tol=1e-4)), 15)

        LiCoO2=Structure.from_file(get_path("icsd_LiCoO2.cif"),
                                          primitive=False)
        gen = SlabGenerator(LiCoO2, [0, 0, 1], 10, 10)
        lco = gen.get_slabs(bonds={("Co", "O"): 3})
        self.assertEqual(len(lco), 1)
        a, b, c = gen.oriented_unit_cell.lattice.matrix
        self.assertAlmostEqual(np.dot(a, gen._normal), 0)
        self.assertAlmostEqual(np.dot(b, gen._normal), 0)

        scc = Structure.from_spacegroup("Pm-3m", Lattice.cubic(3), ["Fe"],
                                        [[0, 0, 0]])
        gen = SlabGenerator(scc, [0, 0, 1], 10, 10)
        slabs = gen.get_slabs()
        self.assertEqual(len(slabs), 1)
        gen = SlabGenerator(scc, [1, 1, 1], 10, 10, max_normal_search=1)
        slabs = gen.get_slabs()
        self.assertEqual(len(slabs), 1)

        # Test whether using units of hkl planes instead of Angstroms for
        # min_slab_size and min_vac_size will give us the same number of atoms
        natoms = []
        for a in [1, 1.4, 2.5, 3.6]:
            s = Structure.from_spacegroup("Im-3m", Lattice.cubic(a), ["Fe"], [[0,0,0]])
            slabgen = SlabGenerator(s, (1,1,1), 10, 10, in_unit_planes=True,
                                    max_normal_search=2)
            natoms.append(len(slabgen.get_slab()))
        n = natoms[0]
        for i in natoms:
            self.assertEqual(n, i)

    def test_triclinic_TeI(self):
        # Test case for a triclinic structure of TeI. Only these three
        # Miller indices are used because it is easier to identify which
        # atoms should be in a surface together. The closeness of the sites
        # in other Miller indices can cause some ambiguity when choosing a
        # higher tolerance.
        numb_slabs = {(0, 0, 1): 5, (0, 1, 0): 3, (1, 0, 0): 7}
        TeI = Structure.from_file(get_path("icsd_TeI.cif"),
                                  primitive=False)
        for k, v in numb_slabs.items():
            trclnc_TeI = SlabGenerator(TeI, k, 10, 10)
            TeI_slabs = trclnc_TeI.get_slabs()
            self.assertEqual(v, len(TeI_slabs))

    def test_get_orthogonal_c_slab(self):
        TeI = Structure.from_file(get_path("icsd_TeI.cif"),
                                  primitive=False)
        trclnc_TeI = SlabGenerator(TeI, (0, 0, 1), 10, 10)
        TeI_slabs = trclnc_TeI.get_slabs()
        slab = TeI_slabs[0]
        norm_slab = slab.get_orthogonal_c_slab()
        self.assertAlmostEqual(norm_slab.lattice.angles[0], 90)
        self.assertAlmostEqual(norm_slab.lattice.angles[1], 90)

    def test_get_tasker2_slabs(self):
        # The uneven distribution of ions on the (111) facets of Halite
        # type slabs are typical examples of Tasker 3 structures. We
        # will test this algo to generate a Tasker 2 structure instead
        slabgen = SlabGenerator(self.MgO, (1,1,1), 10, 10,
                                max_normal_search=1)
        # We generate the Tasker 3 structure first
        slab = slabgen.get_slabs()[0]
        self.assertFalse(slab.is_symmetric())
        self.assertTrue(slab.is_polar())
        # Now to generate the Tasker 2 structure, we must
        # ensure there are enough ions on top to move around
        slab.make_supercell([2,1,1])
        slabs = slab.get_tasker2_slabs()
        # Check if our Tasker 2 slab is nonpolar and symmetric
        for slab in slabs:
            self.assertTrue(slab.is_symmetric())
            self.assertFalse(slab.is_polar())

    def test_nonstoichiometric_symmetrized_slab(self):
        # For the (111) halite slab, sometimes a nonstoichiometric
        # system is preferred over the stoichiometric Tasker 2.
        slabgen = SlabGenerator(self.MgO, (1,1,1), 10, 10,
                                max_normal_search=1)
        slabs = slabgen.get_slabs(symmetrize=True)

        # We should end up with two terminations, one with
        # an Mg rich surface and another O rich surface
        self.assertEqual(len(slabs), 2)
        for slab in slabs:
            self.assertTrue(slab.is_symmetric())

        # For a low symmetry primitive_elemental system such as
        # R-3m, there should be some nonsymmetric slabs
        # without using nonstoichiometric_symmetrized_slab
        slabs = generate_all_slabs(self.Dy, 1, 30, 30,
                                   center_slab=True, symmetrize=True)
        for s in slabs:
            self.assertTrue(s.is_symmetric())
            self.assertGreater(len(s), len(self.Dy))

    def test_move_to_other_side(self):

        # Tests to see if sites are added to opposite side
        s = self.get_structure("LiFePO4")
        slabgen = SlabGenerator(s, (0,0,1), 10, 10, center_slab=True)
        slab = slabgen.get_slab()
        surface_sites = slab.get_surface_sites()

        # check if top sites are moved to the bottom
        top_index = [ss[1] for ss in surface_sites["top"]]
        slab = slabgen.move_to_other_side(slab, top_index)
        all_bottom = [slab[i].frac_coords[2] < slab.center_of_mass[2]
                      for i in top_index]
        self.assertTrue(all(all_bottom))

        # check if bottom sites are moved to the top
        bottom_index = [ss[1] for ss in surface_sites["bottom"]]
        slab = slabgen.move_to_other_side(slab, bottom_index)
        all_top = [slab[i].frac_coords[2] > slab.center_of_mass[2]
                      for i in bottom_index]
        self.assertTrue(all(all_top))


class ReconstructionGeneratorTests(PymatgenTest):

    def setUp(self):

        l = Lattice.cubic(3.51)
        species = ["Ni"]
        coords = [[0,0,0]]
        self.Ni = Structure.from_spacegroup("Fm-3m", l, species, coords)
        l = Lattice.cubic(2.819000)
        species = ["Fe"]
        coords = [[0,0,0]]
        self.Fe = Structure.from_spacegroup("Im-3m", l, species, coords)
        self.Si = Structure.from_spacegroup("Fd-3m", Lattice.cubic(5.430500),
                                            ["Si"], [(0, 0, 0.5)])

        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "..",
                               "reconstructions_archive.json")) as data_file:
            self.rec_archive = json.load(data_file)

    def test_build_slab(self):

        # First lets test a reconstruction where we only remove atoms
        recon = ReconstructionGenerator(self.Ni, 10, 10,
                                        "fcc_110_missing_row_1x2")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        self.assertTrue(recon_slab.reconstruction)
        self.assertEqual(len(slab), len(recon_slab)+2)
        self.assertTrue(recon_slab.is_symmetric())

        # Test if the ouc corresponds to the reconstructed slab
        recon_ouc = recon_slab.oriented_unit_cell
        ouc = slab.oriented_unit_cell
        self.assertEqual(ouc.lattice.b*2, recon_ouc.lattice.b)
        self.assertEqual(len(ouc)*2, len(recon_ouc))

        # Test a reconstruction where we simply add atoms
        recon = ReconstructionGenerator(self.Ni, 10, 10,
                                        "fcc_111_adatom_t_1x1")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        self.assertEqual(len(slab), len(recon_slab)-2)
        self.assertTrue(recon_slab.is_symmetric())

        # If a slab references another slab,
        # make sure it is properly generated
        recon = ReconstructionGenerator(self.Ni, 10, 10,
                                        "fcc_111_adatom_ft_1x1")
        slab = recon.build_slabs()[0]
        self.assertTrue(slab.is_symmetric)

        # Test a reconstruction where it works on a specific
        # termination (Fd-3m (111))
        recon = ReconstructionGenerator(self.Si, 10, 10,
                                        "diamond_111_1x2")
        slab = recon.get_unreconstructed_slabs()[0]
        recon_slab = recon.build_slabs()[0]
        self.assertEqual(len(slab), len(recon_slab)-8)
        self.assertTrue(recon_slab.is_symmetric())

        # Test a reconstruction where terminations give
        # different reconstructions with a non-primitive_elemental system

    def test_get_d(self):

        # Ensure that regardles of the size of the vacuum or slab
        # layer, the spacing between atomic layers should be the same

        recon = ReconstructionGenerator(self.Si, 10, 10,
                                        "diamond_100_2x1")

        recon2 = ReconstructionGenerator(self.Si, 20, 10,
                                         "diamond_100_2x1")
        s1 = recon.get_unreconstructed_slabs()[0]
        s2 = recon2.get_unreconstructed_slabs()[0]
        self.assertAlmostEqual(get_d(s1), get_d(s2))

    def test_previous_reconstructions(self):

        # Test to see if we generated all reconstruction
        # types correctly and nothing changes

        m = StructureMatcher()
        for n in self.rec_archive.keys():
            if "base_reconstruction" in self.rec_archive[n].keys():
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
            s = Structure.from_file(get_path(os.path.join("reconstructions",
                                                          el+"_"+n+".cif")))
            self.assertTrue(any([len(m.group_structures([s, slab]))==1 for slab in slabs]))


class MillerIndexFinderTests(PymatgenTest):

    def setUp(self):
        self.cscl = Structure.from_spacegroup(
            "Pm-3m", Lattice.cubic(4.2), ["Cs", "Cl"],
            [[0, 0, 0], [0.5, 0.5, 0.5]])
        self.Fe = Structure.from_spacegroup(\
            "Im-3m", Lattice.cubic(2.82), ["Fe"],
            [[0, 0, 0]])
        self.lifepo4 = self.get_structure("LiFePO4")
        self.tei = Structure.from_file(get_path("icsd_TeI.cif"),
                                       primitive=False)
        self.LiCoO2 = Structure.from_file(get_path("icsd_LiCoO2.cif"),
                                          primitive=False)

        self.p1 = Structure(Lattice.from_parameters(3, 4, 5, 31, 43, 50),
                            ["H", "He"], [[0, 0, 0], [0.1, 0.2, 0.3]])
        self.graphite = self.get_structure("Graphite")

    def test_get_symmetrically_distinct_miller_indices(self):

        # Tests to see if the function obtains the known number of unique slabs

        indices = get_symmetrically_distinct_miller_indices(self.cscl, 1)
        self.assertEqual(len(indices), 3)
        indices = get_symmetrically_distinct_miller_indices(self.cscl, 2)
        self.assertEqual(len(indices), 6)

        self.assertEqual(len(get_symmetrically_distinct_miller_indices(self.lifepo4, 1)), 7)

        # The TeI P-1 structure should have 13 unique millers (only inversion
        # symmetry eliminates pairs)
        indices = get_symmetrically_distinct_miller_indices(self.tei, 1)
        self.assertEqual(len(indices), 13)

        # P1 and P-1 should have the same # of miller indices since surfaces
        # always have inversion symmetry.
        indices = get_symmetrically_distinct_miller_indices(self.p1, 1)
        self.assertEqual(len(indices), 13)

        indices = get_symmetrically_distinct_miller_indices(self.graphite, 2)
        self.assertEqual(len(indices), 12)

    def test_generate_all_slabs(self):

        slabs = generate_all_slabs(self.cscl, 1, 10, 10)
        # Only three possible slabs, one each in (100), (110) and (111).
        self.assertEqual(len(slabs), 3)

        # make sure it generates reconstructions
        slabs = generate_all_slabs(self.Fe, 1, 10, 10,
                                   include_reconstructions = True)

        # Four possible slabs, (100), (110), (111) and the zigzag (100).
        self.assertEqual(len(slabs), 4)

        slabs = generate_all_slabs(self.cscl, 1, 10, 10,
                                   bonds={("Cs", "Cl"): 4})
        # No slabs if we don't allow broken Cs-Cl
        self.assertEqual(len(slabs), 0)

        slabs = generate_all_slabs(self.cscl, 1, 10, 10,
                                   bonds={("Cs", "Cl"): 4},
                                   max_broken_bonds=100)
        self.assertEqual(len(slabs), 3)

        slabs2 = generate_all_slabs(self.lifepo4, 1, 10, 10,
                                    bonds={("P", "O"): 3, ("Fe", "O"): 3})
        self.assertEqual(len(slabs2), 0)

        # There should be only one possible stable surfaces, all of which are
        # in the (001) oriented unit cell
        slabs3 = generate_all_slabs(self.LiCoO2, 1, 10, 10,
                                    bonds={("Co", "O"): 3})
        self.assertEqual(len(slabs3), 1)
        mill = (0, 0, 1)
        for s in slabs3:
            self.assertEqual(s.miller_index, mill)

        slabs1 = generate_all_slabs(self.lifepo4, 1, 10, 10, tol=0.1,
                                    bonds={("P", "O"): 3})
        self.assertEqual(len(slabs1), 4)

        # Now we test this out for repair_broken_bonds()
        slabs1_repair = generate_all_slabs(self.lifepo4, 1, 10, 10, tol=0.1,
                                    bonds={("P", "O"): 3}, repair=True)
        self.assertGreater(len(slabs1_repair), len(slabs1))

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
            self.assertFalse(any(broken))

        # check if we were able to produce at least one
        # termination for each distinct Miller _index
        self.assertEqual(len(miller_list), len(all_miller_list))

    def test_miller_index_from_sites(self):
        # test on a cubic system
        m = Lattice.cubic(1).matrix
        s1 = np.array([0.5, -1.5, 3])
        s2 = np.array([0.5, 3.,-1.5])
        s3 = np.array([2.5, 1.5,-4.])
        self.assertEqual(tuple(miller_index_from_sites(m, [s1, s2, s3])),
                         (-2,-1,-1))

        # test on a hexagonal system
        m = np.array([[2.319, -4.01662582, 0.],
                      [2.319, 4.01662582, 0.],
                      [0., 0., 7.252]])

        s1 = np.array([2.319, 1.33887527, 6.3455])
        s2 = np.array([1.1595, 0.66943764, 4.5325])
        s3 = np.array([1.1595, 0.66943764, 0.9065])
        hkl = [np.round(i, 6) for i in miller_index_from_sites(m, [s1, s2, s3])]
        self.assertEqual(tuple(hkl), (2, -1, 0))


if __name__ == "__main__":
    unittest.main()
