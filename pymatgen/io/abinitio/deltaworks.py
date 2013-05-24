#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import os.path
import collections
import numpy as np

from pymatgen.io.abinitio.abiobjects import Smearing
from pymatgen.io.abinitio.workflow import DeltaTest

try:
    from pymatgen.io.abinitio.deltafactor.dataset import DeltaFactorDataset
except ImportError:
    pass

##########################################################################################

class DeltaFactoryError(Exception):
    "Base Error class"

class CIFNotFoundError(DeltaFactoryError):
    "CIF file not found in CIFs directory"

class DeltaFactory(object):
    """
    Factory class producing work objects for the computation of the delta factor.
    """
    Error = DeltaFactoryError

    def __init__(self):
        self.delta_data = DeltaFactorDataset()

    def work_for_pseudo(self, workdir, runmode, pseudo, accuracy="normal", kppa=6750, 
        ecut=None, smearing="fermi_dirac:0.0005"):
        """
        Returns a Work object from the given pseudopotential.

        :Note: 0.001 Rydberg is the value used with WIEN2K
        """
        try:
            cif_path = self.delta_data.cif_paths[pseudo.symbol]
        except KeyError:
            raise CIFNotFoundError("%s: cannot find CIF file for pseudo" % pseudo.name)

        # Include spin polarization for O, Cr and Mn (antiferromagnetic) and Fe, Co, and Ni (ferromagnetic). 
        spin_mode = "unpolarized"

        if pseudo.symbol in ["Fe", "Co", "Ni"]: spin_mode = "polarized"
        if pseudo.symbol in ["O", "Cr", "Mn"]: spin_mode = "afm"

        work = DeltaTest(workdir, runmode, cif_path, pseudo, kppa,
                         spin_mode=spin_mode, smearing=smearing, accuracy=accuracy, ecut=ecut, ecutsm=0.05,
                        )
        return work

##########################################################################################

class DeltaDrone(object):

    def assimilate(self, top):
        "Assimilate the results"
        self.results = d = {}
        self.errors = {}

        for dirpath, dirnames, filenames in os.walk(top):
            for fname in filenames:

                if fname == "deltadata.txt":
                    path = os.path.join(dirpath, fname)
                    ppname, x = os.path.split(path)
                    ppname = os.path.basename(ppname)

                    try:
                        d[ppname] = self._analyze(path)
                    except Exception as exc:
                        self.errors[ppname] = str(exc)

    def dump_results(self, stream=sys.stdout, title=""):
        s = stream
        if title:
            s.write("# " + title + "\n")

        for (ppname, res) in self.results.items():
            idx = ppname.find("-")
            symbol = ppname[:idx]
            s.write("%s %.5f %.5f %.5f # %s\n" % (symbol, res[0], res[1], res[2], ppname))

    @staticmethod
    def _analyze(path):
        import eosfit
        data = np.loadtxt(path)
        volume, bulk_modulus, bulk_deriv, residuals = eosfit.BM(data)
        echarge = 1.60217733e-19
        return volume, bulk_modulus * echarge * 1.0e21, bulk_deriv

##########################################################################################
