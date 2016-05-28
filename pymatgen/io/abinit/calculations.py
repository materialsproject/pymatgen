# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
Factory functions producing ABINIT Works.
Works are packed together in a flow. A flow can be ran using abirun (abipy)
Entry points for client code (high-level interface)
"""
from __future__ import unicode_literals, division, print_function

import os

from .abiobjects import KSampling, Screening, SelfEnergy, ExcHamiltonian, HilbertTransform
#from .strategies import ScfStrategy, NscfStrategy, ScreeningStrategy, SelfEnergyStrategy, MdfBse_Strategy
from .works import BandStructureWork, G0W0Work, BseMdfWork


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"


#def bandstructure_work(structure, pseudos, scf_kppa, nscf_nband,
#                       ndivsm, accuracy="normal", spin_mode="polarized",
#                       smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
#                       dos_kppa=None, workdir=None, manager=None, work_class=None, **extra_abivars):
#    """
#    Returns a :class:`Work` for bandstructure calculations.
#
#    Args:
#        structure: Pymatgen structure.
#        pseudos: List of `Pseudo` objects.
#        scf_kppa: Defines the sampling used for the SCF run.
#        nscf_nband: Number of bands included in the NSCF run.
#        ndivs: Number of divisions used to sample the smallest segment of the k-path.
#        accuracy: Accuracy of the calculation.
#        spin_mode: Spin polarization.
#        smearing: Smearing technique.
#        charge: Electronic charge added to the unit cell.
#        scf_algorithm: Algorithm used for solving of the SCF cycle.
#        dos_kppa: Defines the k-point sampling used for the computation of the DOS
#            (None if DOS is not wanted).
#        workdir: Working directory.
#        manager: :class:`TaskManager` instance.
#        extra_abivars: Dictionary with extra variables passed to ABINIT.
#    """
#    #multi = MultiDataset(structure, pseudos, ndtset=2 if dos_kppa is None else 2 + len(dos_kppa))
#
#    # SCF calculation.
#    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
#
#    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
#                               accuracy=accuracy, spin_mode=spin_mode,
#                               smearing=smearing, charge=charge,
#                               scf_algorithm=scf_algorithm, **extra_abivars)
#
#    #scf_electrons = Electrons(spin_mode=spin_mode, smearing=smearing, algorithm=scf_algorithm,
#    #                          charge=charge, nband=scf_nband, fband=None)
#    #multi[0].set_vars(scf_ksampling.to_abivars())
#    #multi[0].set_vars(scf_electrons.to_abivars())
#
#    # Band structure calculation.
#    nscf_ksampling = KSampling.path_from_structure(ndivsm, structure)
#
#    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)
#
#    # DOS calculation.
#    dos_strategy = None
#    if dos_kppa is not None:
#        dos_ksampling = KSampling.automatic_density(structure, dos_kppa, chksymbreak=0)
#        #dos_ksampling = KSampling.monkhorst(dos_ngkpt, shiftk=dos_shiftk, chksymbreak=0)
#        dos_strategy = NscfStrategy(scf_strategy, dos_ksampling, nscf_nband, nscf_solver=None, **extra_abivars)
#        #dos_electrons = aobj.Electrons(spin_mode=spin_mode, smearing=smearing, algorithm={"iscf": -2},
#        #                               charge=charge, nband=nscf_nband)
#
#        #dt = 2 + i
#        #multi[dt].set_vars(dos_ksampling.to_abivars())
#        #multi[dt].set_vars(dos_electrons.to_abivars())
#        #multi[dt].set_vars(_stopping_criterion("nscf", accuracy))
#
#    if work_class is None: work_class = BandStructureWork
#    return work_class(scf_strategy, nscf_strategy, dos_inputs=dos_strategy, workdir=workdir, manager=manager)
#
#
#def g0w0_with_ppmodel_work(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx,
#                           accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
#                           ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None,
#                           sigma_nband=None, gw_qprange=1, workdir=None, manager=None, work_class=None, **extra_abivars):
#    """
#    Returns a :class:`Work` object that performs G0W0 calculations for the given the material.
#
#    Args:
#        structure: Pymatgen structure.
#        pseudos: List of `Pseudo` objects.
#        scf_kppa: Defines the sampling used for the SCF run.
#        nscf_nband: Number of bands included in the NSCF run.
#        ecuteps: Cutoff energy [Ha] for the screening matrix.
#        ecutsigx: Cutoff energy [Ha] for the exchange part of the self-energy.
#        accuracy: Accuracy of the calculation.
#        spin_mode: Spin polarization.
#        smearing: Smearing technique.
#        ppmodel: Plasmonpole technique.
#        charge: Electronic charge added to the unit cell.
#        scf_algorithm: Algorithm used for solving of the SCF cycle.
#        inclvkb: Treatment of the dipole matrix elements (see abinit variable).
#        scr_nband: Number of bands used to compute the screening (default is nscf_nband)
#        sigma_nband: Number of bands used to compute the self-energy (default is nscf_nband)
#        gw_qprange: Option for the automatic selection of k-points and bands for GW corrections.
#            See Abinit docs for more detail. The default value makes the code compute the
#            QP energies for all the point in the IBZ and one band above and one band below the Fermi level.
#        workdir: Working directory.
#        manager: :class:`TaskManager` instance.
#        extra_abivars: Dictionary with extra variables passed to ABINIT.
#    """
#    # TODO: Cannot use istwfk != 1.
#    if "istwfk" not in extra_abivars:
#        extra_abivars["istwfk"] = "*1"
#
#    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
#
#    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
#                               accuracy=accuracy, spin_mode=spin_mode,
#                               smearing=smearing, charge=charge,
#                               scf_algorithm=scf_algorithm, **extra_abivars)
#
#    nscf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
#
#    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)
#
#    if scr_nband is None: scr_nband = nscf_nband
#    if sigma_nband is None: sigma_nband = nscf_nband
#
#    screening = Screening(ecuteps, scr_nband, w_type="RPA", sc_mode="one_shot",
#                          hilbert=None, ecutwfn=None, inclvkb=inclvkb)
#
#    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
#                             gw_qprange=gw_qprange, ppmodel=ppmodel)
#
#    scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening, **extra_abivars)
#
#    sigma_strategy = SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy,
#                                        **extra_abivars)
#
#    if work_class is None: work_class = G0W0Work
#    return work_class(scf_strategy, nscf_strategy, scr_strategy, sigma_strategy, workdir=workdir, manager=manager)


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


#def bse_with_mdf_work(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk,
#                      ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf,
#                      exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized",
#                      smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, workdir=None, manager=None,
#                      work_class=None, **extra_abivars):
#    """
#    Returns a :class:`Work` object that performs a GS + NSCF + Bethe-Salpeter calculation.
#    The self-energy corrections are approximated with the scissors operator.
#    The screening in modeled by the model dielectric function.
#
#    Args:
#        structure: :class:`Structure` object.
#        pseudos: List of `Pseudo` objects.
#        scf_kppa: Defines the sampling used for the SCF run.
#        nscf_nband: Number of bands included in the NSCF run.
#        nscf_ngkpt: Divisions of the k-mesh used for the NSCF and the BSE run.
#        nscf_shiftk: Shifts used for the NSCF and the BSE run.
#        ecuteps: Cutoff energy [Ha] for the screening matrix.
#        bs_loband: Index of the first occupied band included the e-h basis set
#            (ABINIT convention i.e. first band starts at 1).
#            Can be scalar or array of shape (nsppol,)
#        bs_nband: Highest band idex used for the construction of the e-h basis set.
#        soenergy: Scissor energy in Hartree.
#        mdf_epsinf: Value of the macroscopic dielectric function used in expression for the model dielectric function.
#        exc_type: Approximation used for the BSE Hamiltonian (Tamm-Dancoff or coupling).
#        bs_algo: Algorith for the computatio of the macroscopic dielectric function.
#        accuracy: Accuracy of the calculation.
#        spin_mode: Spin polarization.
#        smearing: Smearing technique.
#        charge: Electronic charge added to the unit cell.
#        scf_algorithm: Algorithm used for solving the SCF cycle.
#        workdir: Working directory.
#        manager: :class:`TaskManger` instance.
#        extra_abivars: Dictionary with extra variables passed to ABINIT.
#    """
#    # TODO: Cannot use istwfk != 1.
#    if "istwfk" not in extra_abivars:
#        extra_abivars["istwfk"] = "*1"
#
#    # Ground-state strategy.
#    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
#
#    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
#                               accuracy=accuracy, spin_mode=spin_mode,
#                               smearing=smearing, charge=charge, scf_algorithm=None, **extra_abivars)
#
#    # NSCF calculation with the randomly-shifted k-mesh.
#    nscf_ksampling = KSampling.monkhorst(nscf_ngkpt, shiftk=nscf_shiftk, chksymbreak=0)
#
#    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)
#
#    # Strategy for the BSE calculation.
#    exc_ham = ExcHamiltonian(bs_loband, bs_nband, soenergy, coulomb_mode="model_df", ecuteps=ecuteps,
#                             spin_mode=spin_mode, mdf_epsinf=mdf_epsinf, exc_type=exc_type, algo=bs_algo,
#                             bs_freq_mesh=None, with_lf=True, zcut=None)
#
#    bse_strategy = MdfBse_Strategy(scf_strategy, nscf_strategy, exc_ham, **extra_abivars)
#
#    if work_class is None: work_class = BseMdfWork
#    return work_class(scf_strategy, nscf_strategy, bse_strategy, workdir=workdir, manager=manager)
