# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Factory functions producing ABINIT Works.
Works are packed together in a flow. A flow can be ran using abirun (abipy)
Entry points for client code (high-level interface)
"""

import os

from .abiobjects import KSampling, Screening, SelfEnergy, ExcHamiltonian, HilbertTransform
from .works import BandStructureWork, G0W0Work, BseMdfWork


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"



def g0w0_extended_work(structure, pseudos, kppa, nscf_nband, ecuteps, ecutsigx, scf_nband, accuracy="normal",
                       spin_mode="polarized", smearing="fermi_dirac:0.1 eV", response_models=["godby"], charge=0.0,
                       inclvkb=2, scr_nband=None, sigma_nband=None, workdir=None, manager=None, gamma=True, nksmall=20,
                       work_class=None, **extra_abivars):
    """
    Returns a :class:`Work` object that performs G0W0 calculations for the given the material.

    Args:
        structure: Pymatgen structure.
        pseudos: List of `Pseudo` objects.
        scf_ Defines the sampling used for the SCF run.
        nscf_nband: Number of bands included in the NSCF run.
        ecuteps: Cutoff energy [Ha] for the screening matrix.
        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        ppmodel: Plasmonpole technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        inclvkb: Treatment of the dipole matrix elements (see abinit variable).
        scr_nband: Number of bands used to compute the screening (default is nscf_nband)
        sigma_nband: Number of bands used to compute the self-energy (default is nscf_nband)
        workdir: Working directory.
        manager: :class:`TaskManager` instance.
        nksamll: if not None, a DFT bandstucture calculation will be added after the sc run
        extra_abivars: Dictionary with extra variables passed to ABINIT.
    """
    # TODO: Cannot use istwfk != 1.

    # all these too many options are for development only the current idea for the final version is
    #if gamma:
    #    scf_ksampling = KSampling.automatic_density(structure=structure, kppa=10000, chksymbreak=0, shifts=(0, 0, 0))
    #    nscf_ksampling = KSampling.gamma_centered(kpts=(2, 2, 2))
    #    if kppa <= 13:
    #        nscf_ksampling = KSampling.gamma_centered(kpts=(scf_kppa, scf_kppa, scf_kppa))
    #    else:
    #        nscf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0, shifts=(0, 0, 0))
    #else:
    #    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
    #    nscf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    if gamma:
        if kppa == 1:
            scf_ksampling = KSampling.gamma_centered(kpts=(1, 1, 1))
            nscf_ksampling = KSampling.gamma_centered(kpts=(1, 1, 1))
        elif kppa == 2:
            scf_ksampling = KSampling.gamma_centered(kpts=(2, 2, 2))
            nscf_ksampling = KSampling.gamma_centered(kpts=(2, 2, 2))
        elif kppa < 0:
            scf_ksampling = KSampling.gamma_centered(kpts=(-kppa, -kppa, -kppa))
            nscf_ksampling = KSampling.gamma_centered(kpts=(2, 2, 2))
        elif kppa <= 13:
            scf_ksampling = KSampling.gamma_centered(kpts=(kppa, kppa, kppa))
            nscf_ksampling = KSampling.gamma_centered(kpts=(kppa, kppa, kppa))
        else:
            scf_ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=(0, 0, 0))
            nscf_ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=(0, 0, 0))
    else:
        #this is the original behaviour before the devellopment of the gwwrapper
        scf_ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0)
        nscf_ksampling = KSampling.automatic_density(structure, kppa, chksymbreak=0)

    print(scf_ksampling)
    print(nscf_ksampling)

    if "istwfk" not in extra_abivars:
        extra_abivars["istwfk"] = "*1"

    scf_inputs = []
    to_add = {}
    #scf_nband = min(nscf_nband)
    #print(scf_nband)
    extra_abivars.update(to_add)

    for k in extra_abivars.keys():
        if k[-2:] == '_s':
            var = k[:len(k)-2]
            values = extra_abivars.pop(k)
            to_add.update({k: values[-1]})
            for value in values:
                extra_abivars[var] = value
                extra_abivars['pawecutdg'] = extra_abivars['ecut']*2
                scf_inputs.append(ScfStrategy(structure, pseudos, scf_ksampling, accuracy=accuracy,
                                                spin_mode=spin_mode, smearing=smearing, charge=charge,
                                                scf_algorithm=None, nband=scf_nband, **extra_abivars))

    #temporary for testing a new approach ...
    spread_scr = False if os.path.isfile('no_spread_scr') else True

    if len(scf_strategy) == 0:
        scf_strategy.append(ScfStrategy(structure, pseudos, scf_ksampling, accuracy=accuracy, spin_mode=spin_mode,
                                        smearing=smearing, charge=charge, scf_algorithm=None, nband=scf_nband,
                                        **extra_abivars))


    nscf_strategy = NscfStrategy(scf_strategy[-1], nscf_ksampling, int(max(nscf_nband)*1.1)+1,
                                 nbdbuf=int(0.1*max(nscf_nband)), nstep=200, **extra_abivars)

    if scr_nband is None:
        scr_nband = nscf_nband
    if sigma_nband is None:
        sigma_nband = nscf_nband

    if ecutsigx < max(ecuteps):
        ecutsigx = max(ecuteps)

    sigma_strategy = []

    if 'cd' in response_models:
        hilbert = HilbertTransform(nomegasf=100, domegasf=None, spmeth=1, nfreqre=None, freqremax=None, nfreqim=None,
                                   freqremin=None)

    for response_model in response_models:
        for ecuteps_v in ecuteps:
            for nscf_nband_v in nscf_nband:
                scr_nband = nscf_nband_v
                sigma_nband = nscf_nband_v
                if response_model == 'cd':
                    screening = Screening(ecuteps_v, scr_nband, w_type="RPA", sc_mode="one_shot", hilbert=hilbert,
                                          ecutwfn=None, inclvkb=inclvkb)
                    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening, hilbert=hilbert)
                else:
                    ppmodel = response_model
                    screening = Screening(ecuteps_v, scr_nband, w_type="RPA", sc_mode="one_shot", ecutwfn=None,
                                          inclvkb=inclvkb)
                    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening, ppmodel=ppmodel,
                                             gw_qprange=1)
                scr_strategy = ScreeningStrategy(scf_strategy[-1], nscf_strategy, screening, **extra_abivars)
                sigma_strategy.append(SelfEnergyStrategy(scf_strategy[-1], nscf_strategy, scr_strategy, self_energy,
                                                         **extra_abivars))

    if work_class is None: work_class = G0W0Work
    print(work_class)

    return work_class(scf_strategy, nscf_strategy, scr_strategy, sigma_strategy, workdir=workdir, manager=manager,
                      spread_scr=spread_scr, nksmall=nksmall)
