"""
Factory functions producing ABINIT workflows. Entry points for client code
(high-level interface)
"""
from __future__ import division, print_function

import os

from pymatgen.io.abinitio.abiobjects import (Smearing, KSampling, Screening,
    SelfEnergy, ExcHamiltonian)

from pymatgen.io.abinitio.strategies import (ScfStrategy, NscfStrategy,
    ScreeningStrategy, SelfEnergyStrategy, MDFBSE_Strategy)

from pymatgen.io.abinitio.workflows import (PseudoIterativeConvergence, 
    PseudoConvergence, BandStructureWorkflow, G0W0_Workflow, BSEMDF_Workflow)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"



class PPConvergenceFactory(object):
    """
    Factory object that constructs workflows for analyzing the converge of
    pseudopotentials.
    """
    def work_for_pseudo(self, workdir, manager, pseudo, ecut_range, 
                        toldfe=1.e-8, atols_mev=(10, 1, 0.1), spin_mode="polarized",
                        acell=(8, 9, 10), smearing="fermi_dirac:0.1 eV",):
        """
        Return a `Workflow` object given the pseudopotential pseudo.

        Args:
            workdir:
                Working directory.
            pseudo:
                Pseudo object.
            ecut_range:
                range of cutoff energies in Ha units.
            manager:
                `TaskManager` object.
            toldfe:
                Tolerance on the total energy (Ha).
            atols_mev:
                Tolerances in meV for accuracy in ["low", "normal", "high"]
            spin_mode:
                Spin polarization.
            acell:
                Length of the real space lattice (Bohr units)
            smearing:
                Smearing technique.
        """
        workdir = os.path.abspath(workdir)

        smearing = Smearing.assmearing(smearing)

        if isinstance(ecut_range, slice):
            workflow = PseudoIterativeConvergence(
                workdir, manager, pseudo, ecut_range, atols_mev, 
                toldfe=toldfe, spin_mode=spin_mode, 
                acell=acell, smearing=smearing)

        else:
            workflow = PseudoConvergence(
                workdir, manager, pseudo, ecut_range, atols_mev, 
                toldfe=toldfe, spin_mode=spin_mode, 
                acell=acell, smearing=smearing)

        return workflow


def bandstructure(structure, pseudos, scf_kppa, nscf_nband,
                  ndivsm, accuracy="normal", spin_mode="polarized",
                  smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None,
                  dos_kppa=None, workdir=None, manager=None, **extra_abivars):
    """
    Returns a Work object that computes that bandstructure of the material.

    Args:
        structure:
            Pymatgen structure.
        pseudos:
            List of `Pseudo` objects.
        scf_kppa:
            Defines the sampling used for the SCF run.
        nscf_nband:
            Number of bands included in the NSCF run.
        ndivs:
            Number of divisions used to sample the smallest segment of the
            k-path.
        accuracy:
            Accuracy of the calculation.
        spin_mode:
            Spin polarization.
        smearing:
            Smearing technique.
        charge:
            Electronic charge added to the unit cell.
        scf_algorithm:
            Algorithm used for solving of the SCF cycle.
        dos_kppa:
            Defines the k-point sampling used for the computation of the DOS 
            (None if DOS is not wanted).
        workdir:
            Working directory.
        manager:
            `TaskManager` instance.
        extra_abivars:
            Dictionary with extra variables passed to ABINIT.
    """
    # SCF calculation.
    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
                               accuracy=accuracy, spin_mode=spin_mode,
                               smearing=smearing, charge=charge,
                               scf_algorithm=scf_algorithm, **extra_abivars)

    # Band structure calculation.
    nscf_ksampling = KSampling.path_from_structure(ndivsm, structure)

    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)

    # DOS calculation.
    dos_strategy = None
    if dos_kppa is not None:
        raise NotImplementedError("DOS must be tested")
        dos_ksampling = KSampling.automatic_density(structure, dos_kppa, chksymbreak=0)
        #dos_ksampling = KSampling.monkhorst(dos_ngkpt, shiftk=dos_shiftk, chksymbreak=0)

        dos_strategy = NscfStrategy(scf_strategy, dos_ksampling, nscf_nband, nscf_solver=None, **extra_abivars)

    return BandStructureWorkflow(scf_strategy, nscf_strategy, dos_inputs=dos_strategy, 
                                 workdir=workdir, manager=manager)



#def relaxation(workdir, manager, structure, pseudos, scf_kppa,
#               accuracy="normal", spin_mode="polarized",
#               smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, **extra_abivars):
#    """
#    Returns a Work object that performs structural relaxations.
#
#    Args:
#        workdir:
#            Working directory.
#        manager:
#            `TaskManager` object.
#        structure:
#            Pymatgen structure.
#        pseudos:
#            List of `Pseudo` objects.
#        scf_kppa:
#            Defines the sampling used for the SCF run.
#        accuracy:
#            Accuracy of the calculation.
#        spin_mode:
#            Spin polarization.
#        smearing:
#            Smearing technique.
#        charge:
#            Electronic charge added to the unit cell.
#        scf_algorithm:
#            Algorithm used for solving the SCF cycle.
#    """
#    # SCF calculation.
#    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)
#    relax_algo = 
#
#    relax_strategy = RelaxStrategy(structure, pseudos, scf_ksampling, relax_algo, 
#                                   accuracy=accuracy, spin_mode=spin_mode, smearing=smearing, 
#                                   charge=charge, scf_algorithm=scf_algorithm)
#
#    #return Relaxation(relax_strategy, workdir=workdir, manager=manager)


def g0w0_with_ppmodel(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx, 
                      accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
                      ppmodel="godby", charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None, 
                      sigma_nband=None, gw_qprange=1, workdir=None, manager=None, **extra_abivars):
    """
    Returns a Work object that performs G0W0 calculations for the given the material.

    Args:
        structure:
            Pymatgen structure.
        pseudos:
            List of `Pseudo` objects.
        scf_kppa:
            Defines the sampling used for the SCF run.
        nscf_nband:
            Number of bands included in the NSCF run.
        ecuteps:
            Cutoff energy [Ha] for the screening matrix.
        ecutsigx:
            Cutoff energy [Ha] for the exchange part of the self-energy.
        accuracy:
            Accuracy of the calculation.
        spin_mode:
            Spin polarization.
        smearing:
            Smearing technique.
        ppmodel:
            Plasmonpole technique.
        charge:
            Electronic charge added to the unit cell.
        scf_algorithm:
            Algorithm used for solving of the SCF cycle.
        inclvkb:
            Treatment of the dipole matrix elements (see abinit variable).
        scr_nband:
            Number of bands used to compute the screening (default is nscf_nband)
        sigma_nband:
            Number of bands used to compute the self-energy (default is nscf_nband)
        gw_qprange:
            Option for the automatic selection of k-points and bands for GW corrections.
            See Abinit docs for more detail. The default value makes the code computie the 
            QP energies for all the point in the IBZ and one band above and one band below the Fermi level.
        workdir:
            Working directory.
        manager:
            `TaskManager` instance.
        extra_abivars
            Dictionary with extra variables passed to ABINIT.
    """
    # TODO: Cannot use istwfk != 1.
    if "istwfk" not in extra_abivars:
        extra_abivars["istwfk"] = "*1"

    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
                               accuracy=accuracy, spin_mode=spin_mode,
                               smearing=smearing, charge=charge,
                               scf_algorithm=None, **extra_abivars)

    nscf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)

    if scr_nband is None: scr_nband = nscf_nband
    if sigma_nband is None: sigma_nband = nscf_nband

    screening = Screening(ecuteps, scr_nband, w_type="RPA", sc_mode="one_shot",
                          hilbert=None, ecutwfn=None, inclvkb=inclvkb)

    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
                             ppmodel=ppmodel)

    scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening, **extra_abivars)

    sigma_strategy = SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy,
                                        **extra_abivars)

    return G0W0_Workflow(scf_strategy, nscf_strategy, scr_strategy, sigma_strategy, 
                         workdir=workdir, manager=manager)


#def g0w0_with_cd(structure, pseudos, scf_kppa, nscf_nband, ecuteps, ecutsigx, hilbert,
#                 accuracy="normal", spin_mode="polarized", smearing="fermi_dirac:0.1 eV",
#                 charge=0.0, scf_algorithm=None, inclvkb=2, scr_nband=None, 
#                 sigma_nband=None, workdir=None, manager=None, **extra_abivars):
#    """
#    Returns a Work object that performs G0W0 calculations for the given the material.
#
#    Args:
#        structure:
#            Pymatgen structure.
#        pseudos:
#            List of `Pseudo` objects.
#        scf_kppa:
#            Defines the sampling used for the SCF run.
#        nscf_nband:
#            Number of bands included in the NSCF run.
#        ecuteps:
#            Cutoff energy [Ha] for the screening matrix.
#        ecutsigx:
#            Cutoff energy [Ha] for the exchange part of the self-energy.
#        hilbert:
#            `HilbertTransform` object with the parameters defining the frequency mesh 
#            used for the spectral function and the frequency mesh used for the polarizability
#        accuracy:
#            Accuracy of the calculation.
#        spin_mode:
#            Spin polarization.
#        smearing:
#            Smearing technique.
#        charge:
#            Electronic charge added to the unit cell.
#        scf_algorithm:
#            Algorithm used for solving of the SCF cycle.
#        inclvkb:
#            Treatment of the dipole matrix elements (see abinit variable).
#        scr_nband:
#            Number of bands used to compute the screening (default is nscf_nband)
#        sigma_nband:
#            Number of bands used to compute the self-energy (default is nscf_nband)
#        workdir:
#            Working directory.
#        manager:
#            `TaskManager` instance.
#        extra_abivars
#            Dictionary with extra variables passed to ABINIT.
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
#                               scf_algorithm=None, **extra_abivars)
#
#    nscf_ksampling = KSampling.automatic_density(structure, 1, chksymbreak=0)
#
#    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)
#
#    if scr_nband is None: scr_nband = nscf_nband
#    if sigma_nband is None: sigma_nband = nscf_nband
#
#    screening = Screening(ecuteps, scr_nband, w_type="RPA", sc_mode="one_shot",
#                          hilbert=hilbert, ecutwfn=None, inclvkb=inclvkb)
#
#    self_energy = SelfEnergy("gw", "one_shot", sigma_nband, ecutsigx, screening,
#                             hilbert=hilbert)
#
#    scr_strategy = ScreeningStrategy(scf_strategy, nscf_strategy, screening, 
#                                     **extra_abivars)
#
#    sigma_strategy = SelfEnergyStrategy(scf_strategy, nscf_strategy, scr_strategy, self_energy,
#                                        **extra_abivars)
#
#    return G0W0_Workflow(scf_strategy, nscf_strategy, scr_strategy, sigma_strategy, 
#                         workdir=workdir, manager=manager)


def bse_with_mdf(structure, pseudos, scf_kppa, nscf_nband, nscf_ngkpt, nscf_shiftk, 
                 ecuteps, bs_loband, bs_nband, soenergy, mdf_epsinf, 
                 exc_type="TDA", bs_algo="haydock", accuracy="normal", spin_mode="polarized", 
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, workdir=None, manager=None, 
                 **extra_abivars):
    """
    Returns a `Workflow` object that performs a GS + NSCF + Bethe-Salpeter calculation.
    The self-energy corrections are approximated with the scissors operator. The screening
    in modeled by the model dielectric function.

    Args:
        structure:
            `Structure` object.
        pseudos:
            List of `Pseudo` objects.
        scf_kppa:
            Defines the sampling used for the SCF run.
        nscf_nband:
            Number of bands included in the NSCF run.
        nscf_ngkpt:
            Division of the k-mesh used for the NSCF and the BSE run.
        nscf_shiftk:
            Shifts used for the NSCF and the BSE run.
        ecuteps:
            Cutoff energy [Ha] for the screening matrix.
        bs_loband:
            Index of the first occupied band included the e-h basis set
            (ABINIT convention i.e. first band starts at 1).
            Can be scalar or array of shape (nsppol,)
        bs_nband:
            Highest band idex used for the construction of the e-h basis set.
        soenergy:
            Scissor energy in Hartree.
        mdf_epsinf:
            Value of the macroscopic dielectric function used in expression for the model dielectric function.
        exc_type:
            Approximation used for the BSE Hamiltonian (Tamm-Dancoff or coupling).
        bs_algo:
            Algorith for the computatio of the macroscopic dielectric function.
        accuracy:
            Accuracy of the calculation.
        spin_mode:
            Spin polarization.
        smearing:
            Smearing technique.
        charge:
            Electronic charge added to the unit cell.
        scf_algorithm:
            Algorithm used for solving the SCF cycle.
        workdir:
            Working directory.
        manager:
            `TaskManger` instance.
        extra_abivars:
            Dictionary with extra variables passed to ABINIT.
    """
    # TODO: Cannot use istwfk != 1.
    if "istwfk" not in extra_abivars:
        extra_abivars["istwfk"] = "*1"

    # Ground-state strategy.
    scf_ksampling = KSampling.automatic_density(structure, scf_kppa, chksymbreak=0)

    scf_strategy = ScfStrategy(structure, pseudos, scf_ksampling,
                               accuracy=accuracy, spin_mode=spin_mode,
                               smearing=smearing, charge=charge, scf_algorithm=None, **extra_abivars)

    # NSCF calculation with the randomly-shifted k-mesh.
    nscf_ksampling = KSampling.monkhorst(nscf_ngkpt, shiftk=nscf_shiftk, chksymbreak=0)

    nscf_strategy = NscfStrategy(scf_strategy, nscf_ksampling, nscf_nband, **extra_abivars)

    # Strategy for the BSE calculation.
    exc_ham = ExcHamiltonian(bs_loband, bs_nband, soenergy, coulomb_mode="model_df", ecuteps=ecuteps, 
                             spin_mode=spin_mode, mdf_epsinf=mdf_epsinf, exc_type=exc_type, algo=bs_algo,
                             bs_freq_mesh=None, with_lf=True, zcut=None)

    bse_strategy = MDFBSE_Strategy(scf_strategy, nscf_strategy, exc_ham, **extra_abivars)
    raise NotImplementedError("")

    return BSEMDF_Workflow(scf_strategy, nscf_strategy, bse_strategy, workdir=workdir, manager=manager)
