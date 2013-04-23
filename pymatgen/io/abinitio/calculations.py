"Factory functions producing ABINIT workflows. Entry points for client code (high-level interface)"
from __future__ import division, print_function

import os

#from pprint import pprint
#
from pymatgen.io.abinitio.abiobjects import Smearing, KSampling #, AbiStructure, Screening, PPModel, SelfEnergy, 
from pymatgen.io.abinitio.strategies import ScfStrategy, NscfStrategy
from pymatgen.io.abinitio.workflow import BandStructure, PseudoIterativeConvergence, PseudoConvergence #,GW_Workflow, #, Relaxation, 

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

#__all__ = [
#]

##########################################################################################

class PPConvergenceFactory(object):
    "Factory object"

    def work_for_pseudo(self, workdir, pseudo, ecut_range, runmode="sequential", atols_mev=(10, 1, 0.1), 
                        spin_mode="polarized", acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV",):
        """
        Return a Work object given the pseudopotential pseudo

        Args:
            workdir: 
                Working directory.
            pseudo: 
                Pseudo object.
            ecut_range:
                range of cutoff energies in Ha units.
            runmode: 
            atols_mev:
                Tolerances in meV for accuracy in ["low", "normal", "high"]
            spin_mode:
                Spin polarization.
            acell: 
                Length of the real space lattice (Bohr units)
            smearing: 
                Defines the smearing technique.
        """
        workdir = os.path.abspath(workdir)

        smearing = Smearing.assmearing(smearing)

        if isinstance(ecut_range, slice):
            workflow = PseudoIterativeConvergence(workdir, pseudo, ecut_range, atols_mev,
                                              runmode    = runmode,
                                              spin_mode  = spin_mode,
                                              acell      = acell, 
                                              smearing   = smearing,
                                             )

        else:
            workflow = PseudoConvergence(workdir, pseudo, ecut_range, atols_mev,
                                     runmode    = runmode,
                                     spin_mode  = spin_mode,
                                     acell      = acell, 
                                     smearing   = smearing,
                                    )

        return workflow

##########################################################################################

def bandstructure(workdir, runmode, structure, pseudos, scf_kppa, nscf_nband, ndivsm, 
                  accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV", 
                  charge=0.0, scf_solver=None, dos_kppa=None,
                 ):

    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling, accuracy=accuracy, spin_mode=spin_mode, 
                               smearing=smearing, charge=charge, scf_solver=None)

    nscf_ksampling = KSampling.path_from_structure(ndivsm, structure)

    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, nscf_solver=None)

    dos_strategy = None

    if dos_kppa is not None:
        raise NotImplementedError("DOS must be tested")
        dos_ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0)
        dos_strategy = NscfStrategy(scf_strategy, dos_ksampling, nscf_nband, nscf_solver=None)

    return BandStructure(workdir, runmode, scf_strategy, nscf_strategy, dos_strategy=dos_strategy)

##########################################################################################

def g0w0_with_ppmodel(workdir, runmode, structure, pseudos, scf_kppa, nscf_nband,
                      accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV", 
                      ppmodel="godby", charge=0.0, scf_solver=None,
                     ):
    ecuteps, scr_nband, sigma_nband, ecutsigx = 2, 10, 10, 4

    screening = Screening(ecuteps, scr_nband)
    #print(screening.to_abivars())

    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening)
    print(self_energy.to_abivars())

    scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening)
    #print(scr_strategy.make_input())

    sigma_strategy = SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy)
    print(sigma_strategy.make_input())

    return GW_Workflow(workdir, runmode, scf_strategy, nscf_strategy, scr_strategy, sigma_strategy)

##########################################################################################
