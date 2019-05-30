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
from pymatgen import MPRester
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

    def test_shift_vector(self):
        interface = self.ib.interfaces[0]
        tst = Structure.from_sites(interface.sites)
        interface.vacuum_thickness += 1
        tdm, idm = tst.distance_matrix, interface.distance_matrix
        diff = tdm - idm
        assert (tdm <= idm + 1e-10).all()
        assert (tdm + 0.5 < idm).any()
        interface.vacuum_thickness -= 1
        idm = interface.distance_matrix
        assert (np.abs(tdm - idm) < 1e-10).all()
        interface.ab_shift += np.array([0.2,0.2])
        idm = interface.distance_matrix
        assert (np.abs(tdm - idm) > 0.5).any()

class InterfaceBuilderTest(PymatgenTest):

    @classmethod
    def setUpClass(cls):
        mpr = MPRester()

        si_struct = mpr.get_structure_by_material_id('mp-149')
        sio2_struct = mpr.get_structure_by_material_id('mp-6930')

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
            assert len(interface.film) == len(interface.oriented_film_cell)
            assert len(interface.substrate) == len(interface.oriented_sub_cell)

    def test_lattice(self):
        zsl = ZSLGenerator()
        for ib in self.ibs:
            interface = ib.interfaces[4]
            assert zsl.is_same_vectors(interface.oriented_sub_cell.lattice.matrix[:2],
                                        interface.oriented_film_cell.lattice.matrix[:2])

    def test_structure_preservation(self):
        for ib in self.ibs:
            for interface in ib.interfaces[:2]:
                assert interface.film.is_valid(tol=1.6)
                assert interface.substrate.is_valid(tol=1.6)
