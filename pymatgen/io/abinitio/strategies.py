"""
Strategy objects used to create Abinit input files for particular types of calculations.
"""
from __future__ import division, print_function

import abc
from pprint import pprint, pformat

from .abiobjects import SpinMode, Smearing
from .pseudos import PseudoTable
from .input import SystemCard, ElectronsCard, ControlCard, KpointsCard, Input

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

##########################################################################################

# TODO
#class StrategyAbstractFactory(object):
#    """
#    Abstract factory for strategies. Returns a concrete factory
#
#    Example: abst_factory = GSStrategyAbstractFactory()
#        factory = abs_factory()
#        strategy = factory.make_strategy()
#    """
#    __metaclass__ = abc.ABCMeta
#        
#    #def __init__(self, mode, accuracy):
#        # Dict with Ground-State strategy subclasses
#        gs_subclasses = {}
#        for klass in SCFStrategy.__subclasses__():
#            gs_subclasses[klass.__name__] = klass

##########################################################################################

class Strategy(object):
    """
    A Strategy object generates the abinit input file used for a particular type of calculation
    e.g. ground-state runs, structural relaxations, self-energy calculations ...

    A Strategy can absorb data (e.g. data produced in the previous steps of a workflow) and 
    can use this piece of information to generate/optimize the input variables.
    Strategy object provides the method make_input that builds and returns the abinit input file.

    Attributes:
        accuracy:
            Accuracy of the calculation used to define basic parameters of the run.
            such as tolerances, basis set truncation ... 
    """
    __metaclass__ = abc.ABCMeta

    def __str__(self):
        return "<%s at %s, accuracy = %s>" % (self.__class__.__name__, id(self), self.accuracy)

    def learn(self, **data):
        "Update the data stored in self"
        if not hasattr(self, "_data"):
            self._data = dict(data)
        else:
            if [k in self._data for k in data].count(True) != 0:
                raise ValueError("Keys %s are already present in data" % str([k for k in data]))
            self._data.update(data)

    @property
    def accuracy(self):
        "Accuracy used by the strategy"
        try:
            return self._accuracy
        except AttributeError:
            self.set_accuracy("normal")
            return self._accuracy

    def set_accuracy(self, accuracy):
        "Accuracy setter"
        if hasattr(self, "_accuracy"):
            raise RuntimeError("object already has accuracy %s " % self._accuracy)

        assert accuracy in ["low", "normal", "high",]
        self._accuracy = accuracy

    @property
    def data(self):
        try:
            return self. _data
        except AttributeError:
            return {}

    @abc.abstractmethod
    def make_input(self, *args, **kwargs):
        "Returns an Input instance"

##########################################################################################

class SCFStrategy(Strategy):
    """
    Strategy for ground-state SCF calculations.
    """
    def __init__(self, structure, pseudos, ksampling, accuracy="normal", spin_mode="polarized", 
                 smearing="fermi_dirac:0.1 eV", charge=0, scf_algorithm=None, **kwargs):
        """
        Args:
            structure:
                pymatgen structure
            pseudos:
                List of pseudopotentials.
            ksampling:
                Ksampling object defining the sampling of the BZ.
            accuracy:
                Accuracy of the calculation.
            spin_mode: 
                Flag defining the spin polarization (nsppol, nspden, nspinor). Defaults to "polarized"
            smearing: 
                String or Smearing instance. 
            charge:
                Total charge of the system. Default is 0.
            scf_algorithm:
                ElectronsAlgorithm instance.
            **kwargs:
                Extra variables that will be directly added to the input file
        """
        super(SCFStrategy, self).__init__()

        self.set_accuracy(accuracy)

        self.structure  = structure
        self.pseudos    = PseudoTable.astable(pseudos)
        self.ksampling  = ksampling
        self.spin_mode  = SpinMode.asspinmode(spin_mode) 
        self.smearing   = Smearing.assmearing(smearing)
        self.charge     = charge
        self.scf_algo   = scf_algorithm if scf_algorithm else {}

        self._kwargs = kwargs

    def _make_systemcard(self):
        return SystemCard(self.structure, self.pseudos)

    def _make_electronscard(self):
        return ElectronsCard(spin_mode=self.spin_mode, smearing=self.smearing)

    def _make_kpointscard(self):
        # FIXME
        return KpointsCard()

    def make_input(self):
        # Initialize the system section from structure.
        system = self._make_systemcard()

        # Variables for electrons.
        electrons = self._make_electronscard()

        # K-point sampling.
        kpoints = self._make_kpointscard()

        scf_control = ControlCard(system, electrons, kpoints, **self._kwargs)    
                                                                         
        return Input(system, electrons, kpoints, scf_control)              

##########################################################################################

class NSCFStrategy(Strategy):
    """
    Strategy for non-self-consiste calculations.
    """
    def __init__(self, scf_strategy, nscf_nband, ksampling, nscf_algo=None, **kwargs):
        """
        Args:
            scf_strategy:
                SCFStrategy
            ksampling:
                Ksampling object defining the sampling of the BZ.
            nscf_nband:
                Number of bands to compute.
            nscf_algo
                ElectronsAlgorithm instance.
            **kwargs:
                Extra variables that will be directly added to the input file
        """
        super(NSCFStrategy, self).__init__()

        self.set_accuracy(scf_strategy.accuracy)

        self.scf_strategy = scf_strategy

        self.nscf_nband = nscf_nband
        self.ksampling  = ksampling

        if nscf_algo:
            self.nscf_algo = nscf_algo 
        else:
            self.nscf_algo = {"iscf": -3}

        self._kwargs = kwargs

    def _make_systemcard(self):
        scf_strategy = self.scf_strategy
        return SystemCard(scf_strategy.structure, scf_strategy.pseudos)

    def _make_electronscard(self):
        scf_strategy = self.scf_strategy
        #FIXME algo
        return ElectronsCard(spin_mode=scf_strategy.spin_mode, smearing=scf_strategy.smearing, charge=scf_strategy.charge, iscf=-3)

    def _make_kpointscard(self):
        # FIXME
        return KpointsCard()

    def make_input(self):
        # Initialize the system section from structure.
        system = self._make_systemcard()

        # Variables for electrons.
        electrons = self._make_electronscard()

        # K-point variables
        kpoints = self._make_kpointscard()

        scf_control = ControlCard(system, electrons, kpoints, **self._kwargs)    
                                                                         
        return Input(system, electrons, kpoints, scf_control)              

##########################################################################################

class RelaxStrategy(SCFStrategy):

    def __init__(structure, pseudos, ksampling, accuracy="normal", spin_mode="polarized", 
                 smearing="fermi_dirac:0.1 eV", charge=0, scf_algorithm=None, **kwargs):
        """
        Args:
            structure:
                pymatgen structure
            pseudos:
                List of pseudopotentials.
            ksampling:
                Ksampling object defining the sampling of the BZ.
            accuracy:
                Accuracy of the calculation.
            spin_mode: 
                Flag defining the spin polarization (nsppol, nspden, nspinor). Defaults to "polarized"
            smearing: 
                String or Smearing instance. 
            charge:
                Total charge of the system. Default is 0.
            scf_algorithm:
                ElectronsAlgorithm instance.
            **kwargs:
                Extra variables that will be directly added to the input file
        """
        super(RelaxStrategy, self).__init__(structure, pseudos, ksampling, 
                 accuracy=accuracy, spin_mode=spin_mode, 
                 smearing=smearing, charge=charge, scf_algorithm=scf_algorithm, **kwargs)

        self._kwargs = kwargs

    def _make_relaxcard(self):
        raise NotImplementedError("")

    def make_input(self):
        # Initialize the system section from structure.
        system = self._make_systemcard()
                                                                                        
        # Variables for electrons.
        electrons = self._make_electronscard()
                                                                                        
        kpoints = self._make_kpointscard()

        relax = self._make_relaxcard()
                                                                                        
        control = ControlCard(system, electrons, kmesh, relax, **self._kwargs)
                                                                         
        return Input(system, electrons, kpoints, relax, control)

##########################################################################################

class ScreeningStrategy(Strategy):

    def __init__(self, scf_strategy, nscf_strategy, screening, **kwargs):
        """
        Constructor for screening calculations.
                                                                                                       
        Args:
            scf_strategy:
                Strategy used for the ground-state calculation
            nscf_strategy:
                Strategy used for the non-self consistent calculation
            screening:
                Screening instance
            **kwargs:
                Extra variables added directly to the input file
        """
        super(ScreeningStrategy, self).__init__()

        self.scf_strategy = scf_strategy
        self.nscf_strategy = nscf_strategy
        self.screening = screening
        self._kwargs = kwargs

    def _make_screeningcard(self):
        raise NotImplementedError("")
        scr_cards.append(ScreeningCard(scr_strategy, comment="Generated from NSCF input"))
                                                                                         
    def make_input(self):
        scr_cards = nscf_input.copy_cards(exclude=["ElectronsCard",])

        nband_screening = scr_strategy.get_varvalue("nband")

        scr_cards.append(ElectronsCard(spin_mode=nscf_input.spin_mode, nband=nband_screening, smearing=smearing))

        scr_cards.append(ScreeningCard(scr_strategy, comment="Generated from NSCF input"))

        scr_cards.append(ControlCard(*scr_cards, **kwargs))

        return Input(*scr_cards)

##########################################################################################

class SelfEnergyStrategy(Strategy):

    def __init__(self, scf_strategy, nscf_strategy, scr_strategy, **kwargs):
        """
        Constructor for screening calculations.
                                                                                                       
        Args:
            scf_strategy:
                Strategy used for the ground-state calculation
            nscf_strategy:
                Strategy used for the non-self consistent calculation
            scr_strategy
                Strategy used for the screening calculation
            **kwargs:
                Extra variables added directly to the input file
        """
        # TODO Add consistency check between SCR and SIGMA strategies

        super(ScreeningStrategy, self).__init__()

        self.scf_strategy = scf_strategy
        self.nscf_strategy = nscf_strategy
        self.scr_strategy = scr_strategy
        self._kwargs = kwargs

    def _make_selfenergycard(self):
        raise NotImplementedError("")
        se_card = SelfEnergyCard(sigma_strategy, comment="Generated from SCR input")
                                                                                         
    def make_input(self):
        smearing = Smearing.assmearing(smearing)

        if (smearing != scr_input.smearing):
            raise ValueError("Input smearing differs from the one used in the SCR run")

        sigma_nband = sigma_strategy.get_varvalue("nband")

        se_card = SelfEnergyCard(sigma_strategy, comment="Generated from SCR input")

        sigma_cards = scr_input.copy_cards(exclude=["ScreeningCard", "ElectronsCard",])

        sigma_cards.append(ElectronsCard(spin_mode=scr_input.spin_mode, nband=sigma_nband, smearing=smearing))

        sigma_cards.append(se_card)

        sigma_cards.append(ControlCard(*sigma_cards, **kwargs))

        return Input(*sigma_cards)

##########################################################################################
