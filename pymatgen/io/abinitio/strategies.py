"""
Strategy objects used to create Abinit input files for particular types of calculations.
"""
from __future__ import division, print_function

from pprint import pprint, pformat
import abc

from .pseudos import PseudoTable

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
#        for klass in GroundStateStrategy.__subclasses__():
#            gs_subclasses[klass.__name__] = klass

##########################################################################################

class Strategy(object):
    """
    A Strategy object generates the abinit input files used for a particular type of calculation
    e.g. ground-state runs, structural relaxations, self-energy calculations ...

    A Strategy can absorb data (e.g. data produced in the previous steps of a workflow) and 
    can use this piece of information to generate/optimize the input variables.
    Client code can use the method get_input to obtain an object that describes the abinit input file.

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

    #def _validate_input(self, input):
    #    # Sanity check
    #    for (k, v) in self.vars.items():
    #        if v is MANDATORY:
    #            raise ValueError("%s: key %s must have a value specified by the user" % (self.__class__.__name__, k))

    #def get_input(self):
    #    """Return the list of abinit variables."""
    #    input = self.make_input()
    #    
    #    if self.has_data:
    #        self.optimize_input(input)
    #                                                                                                          
    #    self._validate_abivars(input)
    #                                                                                                          
    #    return input                                                                                       

    @abc.abstractmethod
    def make_input(self, *args, **kwargs):
        "Returns an Input instance"

##########################################################################################

class GroundStateStrategy(Strategy)
    """
    Strategy for ground-state calculations.
    """
    def __init__(self, structure, pseudos, ksampling, accuracy="normal", spin_mode="polarized", 
                 smearing="fermi_dirac:0.1 ev", charge=0, scf_solver=None, **kwargs):
        """
        Args:
            structure:
                pymatgen structure
            pseudos:
                List of pseudopotentials.
            ksamplng:
                Ksampling object defining the sampling of the BZ.
            accuracy:
                Accuracy of the calculation.
            spin_mode: 
                Flag defining the spin polarization (nsppol, nspden, nspinor). Defaults to "polarized"
            smearing: 
                String or Smearing instance. 
            charge:
                Total charge of the system. Default is 0.
            scf_solver:
                SCFSolver instance.
            **kwargs:
                Extra variables that will be directly added to the input file

        Returns:
            AbinitInput instance.
        """
        self.set_accuracy(accuracy)

        self.structure  = structure
        self.pseudos    = PseudoTable.astable(pseudos)
        self.ksampling  = ksampling
        self.spin_mode  = SpinMode.asspinmode(spin_mode) 
        self.smearing   = Smearing.assmearing(smearing)
        self.charge     = self.charge
        self.scf_solver = self.scf_solver if scf_solver else {}

        self._kwargs = kwargs

    def make_scf_input(self):
        # Initialize the system section from structure.
        system = SystemCard(self.structure, self.pseudos)

        # Variables for electrons.
        electrons = ElectronsCard(spin_mode=self.spin_mode, smearing=self.smearing)

        # K-point sampling.
        if [obj is not None for obj in [ngkpt, kppa,]].count(True) != 1:
            raise ValueError("only one variable among ngkpt, kppa should be supplied")

        if ngkpt is not None:
            kmesh = KpointsCard.monkhorst_automatic(structure, ngkpt)
        elif kppa is not None:
            kmesh = KpointsCard.automatic_density(structure, kppa)       

        scf_control = ControlCard(system, electrons, kmesh, **self._kwargs)    
                                                                         
        return Input(system, electrons, kmesh, scf_control)              

##########################################################################################
                                                                                           
if __name__ == "__main__":
    import sys
    gs_strategy = GS_strategy()
