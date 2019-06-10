# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

__author__ = "Kyle Bystrom"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Kyle Bystrom"
__email__ = "kylebystrom@gmail.com"
__date__ = "5/29/2019"

import unittest
from pymatgen.analysis.interface import Interface, InterfaceBuilder
from pymatgen.analysis.substrate_analyzer import ZSLGenerator
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure
import numpy as np 

class InterfaceTest(PymatgenTest):

    @classmethod
    def setUpClass(cls):

        si_struct = cls.get_structure('Si')
        sio2_struct = cls.get_structure('SiO2')

        sga = SpacegroupAnalyzer(si_struct)
        si_conventional = sga.get_conventional_standard_structure()
        sga = SpacegroupAnalyzer(sio2_struct)
        sio2_conventional = sga.get_conventional_standard_structure()
        cls.ib = InterfaceBuilder(si_conventional, sio2_conventional)
        cls.ib.generate_interfaces(substrate_millers=[[1, 0, 0]], film_layers=3, substrate_layers=3)

    def test_offset_vector(self):
        interface = self.ib.interfaces[0]
        init_lattice = interface.lattice.matrix.copy()
        self.assertArrayAlmostEqual(interface.offset_vector, np.array([0,0,2.5]))
        init_film = interface.film
        init_sub = interface.substrate
        tst = Structure.from_sites(interface.sites)
        interface.z_shift += 1
        self.assertArrayAlmostEqual(interface.offset_vector, np.array([0,0,3.5]))
        tdm, idm = tst.distance_matrix, interface.distance_matrix
        diff = tdm - idm
        assert (tdm <= idm + 1e-10).all()
        assert (tdm + 0.5 < idm).any()
        interface.z_shift -= 1
        self.assertArrayAlmostEqual(interface.offset_vector, np.array([0,0,2.5]))
        idm = interface.distance_matrix
        assert (np.abs(tdm - idm) < 1e-10).all()
        interface.ab_shift += np.array([0.2,0.2])
        self.assertArrayAlmostEqual(interface.ab_shift, np.array([0.2,0.2]))
        idm = interface.distance_matrix
        assert (np.abs(tdm - idm) > 0.9).any()
        self.assertArrayAlmostEqual(init_lattice, interface.lattice.matrix)

        self.assertArrayAlmostEqual(init_film.distance_matrix, interface.film.distance_matrix)
        self.assertArrayAlmostEqual(init_sub.distance_matrix, interface.substrate.distance_matrix)

    def test_labels(self):
        interface = self.ib.interfaces[0]
        film = interface.film
        substrate = interface.substrate
        film_sites = [site for i, site in enumerate(interface)\
                        if 'film' in interface.site_properties['interface_label'][i]]
        substrate_sites = [site for i, site in enumerate(interface)\
                        if 'substrate' in interface.site_properties['interface_label'][i]]
        assert film.sites == film_sites
        assert substrate.sites == substrate_sites
        assert len(film) == len(interface.modified_film_structure)
        assert len(substrate) == len(interface.modified_sub_structure)

    def test_vertical_spacing(self):
        interface = self.ib.interfaces[0]
        self.assertAlmostEqual(interface.z_shift,
            np.min(interface.film.cart_coords[:,2]) - np.max(interface.substrate.cart_coords[:,2]))
        self.assertAlmostEqual(interface.lattice.c, interface.vacuum_thickness + interface.z_shift\
                                + np.max(interface.film.cart_coords[:,2])
                                - np.min(interface.film.cart_coords[:,2])\
                                + np.max(interface.substrate.cart_coords[:,2])\
                                - np.min(interface.substrate.cart_coords[:,2]))

    def test_inplane_spacing(self):
        pass

class InterfaceBuilderTest(PymatgenTest):

    @classmethod
    def setUpClass(cls):

        si_struct = cls.get_structure('Si')
        sio2_struct = cls.get_structure('SiO2')

        sga = SpacegroupAnalyzer(si_struct)
        si_conventional = sga.get_conventional_standard_structure()
        sga = SpacegroupAnalyzer(sio2_struct)
        sio2_conventional = sga.get_conventional_standard_structure()

        cls.ibs = []
        cls.ibs.append(cls.make_ib(si_conventional, sio2_conventional, [1,0,0]))
        cls.ibs.append(cls.make_ib(sio2_conventional, si_conventional, [1,0,0]))
        cls.ibs.append(cls.make_ib(si_struct, sio2_struct, [1,0,0]))
        cls.ibs.append(cls.make_ib(sio2_struct, si_struct, [1,0,0]))
        cls.ibs.append(cls.make_ib(si_struct, sio2_struct, [1,1,1]))
        cls.ibs.append(cls.make_ib(sio2_struct, si_struct, [1,1,1]))

    @staticmethod
    def make_ib(substrate, film, miller):
        ib = InterfaceBuilder(substrate, film)
        ib.generate_interfaces(substrate_millers=[miller])
        return ib

    def test_film_and_substrate_sites(self):
        for ib in self.ibs:
            interface = ib.interfaces[0]
            assert len(interface.film) == len(interface.modified_film_structure)
            assert len(interface.substrate) == len(interface.modified_sub_structure)

    def test_lattice(self):
        zsl = ZSLGenerator()
        for ib in self.ibs:
            interface = ib.interfaces[0]
            assert zsl.is_same_vectors(interface.modified_sub_structure.lattice.matrix[:2],
                                        interface.modified_film_structure.lattice.matrix[:2])

    def test_structure_preservation(self):
        for ib in self.ibs:
            for interface in ib.interfaces[:2]:
                assert interface.film.is_valid(tol=1.6)
                assert interface.substrate.is_valid(tol=1.6)
