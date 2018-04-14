# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals

import unittest
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.analysis.defects.core import DefectEntry, Vacancy
from pymatgen.analysis.defects.corrections import FreysoldtCorrection

class DefectsCorrectionsTest(PymatgenTest):
    def test_freysoldt(self):
        #a trivial test for freysoldt correction
        struc = PymatgenTest.get_structure("VO2")
        abc = struc.lattice.abc
        axisdata = [np.arange(0., lattval, 0.01) for lattval in abc]
        bldata = [np.array([1. for u in np.arange(0., lattval, 0.01)])  for lattval in abc]
        dldata = [np.array([(-1-np.cos(2*np.pi*u/lattval)) for u in np.arange(0., lattval, 0.01)])  for lattval in abc]
        params = {'axis_grid':axisdata, 'bulk_planar_averages': bldata, 'defect_planar_averages':dldata}
        vac = Vacancy(struc, struc.sites[0], charge=-3)
        fc = FreysoldtCorrection(15)

        #test electrostatic correction
        es_corr = fc.perform_es_corr(struc.lattice, 15., -3)
        self.assertEqual(es_corr, 2.929238)

        #test potential alignment method
        pot_corr = fc.perform_pot_corr( axisdata[0], bldata[0], dldata[0], struc.lattice, 15., -3, vac.site.coords, 0)
        self.assertEqual(pot_corr, 2.917553311042547)

        #test entry full correction method
        de = DefectEntry(vac, 0., corrections={}, parameters=params, entry_id=None)
        val = fc.get_correction(de)
        self.assertEqual(val['freysoldt_electrostatic'], 2.929238)
        self.assertEqual(val['freysoldt_potential_alignment'], 4.579907163820233)


if __name__ == "__main__":
    unittest.main()
