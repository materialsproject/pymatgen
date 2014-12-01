# coding: utf-8
"""Strategy objects for creating ABINIT calculations."""
from __future__ import unicode_literals, division, print_function

import sys
import os
import abc
import collections
import copy
import six
import numpy as np

from six.moves import map, zip
from monty.string import is_string
from pymatgen.util.string_utils import str_aligned, str_delimited
from .abiobjects import Electrons
from .pseudos import PseudoTable

import logging
logger = logging.getLogger(__name__)


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"


def select_pseudos(pseudos, structure, ret_table=True):
    """
    Given a list of pseudos and a pymatgen structure, extract the pseudopotentials
    for the calculation (useful when we receive an entire periodic table).
    If ret_table is True, the function will return a PseudoTable instead of
    a list of `Pseudo` objects

    Raises:
        ValueError if no pseudo is found or multiple occurrences are found.
    """
    table = PseudoTable.as_table(pseudos)

    pseudos = []
    for symbol in structure.types_of_specie:
        # Get the list of pseudopotentials in table from atom symbol.
        pseudos_for_type = table.pseudos_with_symbol(symbol)

        if not pseudos_for_type:
            raise ValueError("Cannot find pseudo for symbol %s" % symbol)

        if len(pseudos_for_type) > 1:
            raise ValueError("Find multiple pseudos for symbol %s" % symbol)

        pseudos.append(pseudos_for_type[0])

    if ret_table:
        return PseudoTable(pseudos)
    else:
        return pseudos


def order_pseudos(pseudos, structure):
    #logger.info('calling order pseudos')
    return select_pseudos(pseudos, structure) #, ret_table=False)


def num_valence_electrons(pseudos, structure):
    """
    Compute the number of valence electrons from
    a list of pseudopotentials and the crystalline structure.

    Args:
        pseudos: List of strings, list of of pseudos or `PseudoTable` instance.
        structure: Pymatgen structure.

    Raises:
        ValueError if cannot find a pseudo in the input pseudos or if the
        input list contains more than one pseudo for the chemical symbols appearing in structure.
    """
    table = PseudoTable.as_table(pseudos)

    valence = 0.0
    for site in structure:
        entries = table.pseudos_with_symbol(site.specie.symbol)
        if len(entries) != 1:
            raise ValueError("Found %d entries for symbol %s" % (len(entries), site.specie.symbol))
        valence += entries[0].Z_val

    return valence


class AbstractStrategy(six.with_metaclass(abc.ABCMeta, object)):
    """
    A Strategy object generates the ABINIT input file used for a particular type of calculation
    e.g. ground-state runs, structural relaxations, self-energy calculations ...

    A Strategy can absorb data (e.g. data produced in the previous steps of a workflow) and
    can use this piece of information to generate/optimize the input variables.
    Strategy objects must provide the method make_input that builds and returns the abinit input file.

    Attributes:

        pseudos: List of pseudopotentials.
    """

    #@abc.abstractproperty
    #def pseudos(self):

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return self.pseudos.allnc

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return self.pseudos.allpaw

    def num_valence_electrons(self):
        """Number of valence electrons computed from the pseudos and the structure."""
        return num_valence_electrons(self.pseudos, self.structure)

    #@abc.abstractproperty
    #def structure(self):

    #def set_structure(self, structure):
    #    self.structure = structure

    #def change_structure(self, structure):
    #    self.structure = structure

    #def to_abivars(self):
    #def to_dict(self):
    #def from_abivars(cls, d):
    #def from_dict(self, d):

    def copy(self):
        """Shallow copy of self."""
        return copy.copy(self)

    def deepcopy(self):
        """Deep copy of self."""
        return copy.deepcopy(self)


    @abc.abstractmethod
    def make_input(self, *args, **kwargs):
        """Returns an Input instance."""


class StrategyWithInput(object):
    # TODO: Find a better way to do this. I will likely need to refactor the Strategy object
    def __init__(self, abinit_input, deepcopy=True):
        if deepcopy: abinit_input = abinit_input.deepcopy()
        self.abinit_input = abinit_input

    @property
    def pseudos(self):
        # FIXME: pseudos must be order but I need to define an ABC for the Strategies and Inputs.
        # Order pseudos
        pseudos = self.abinit_input.pseudos
        return order_pseudos(pseudos, self.abinit_input.structure)

    def add_extra_abivars(self, abivars):
        """Add variables (dict) to extra_abivars."""
        self.abinit_input.set_variables(**abivars)

    def remove_extra_abivars(self, keys):
        """Remove variables from extra_abivars."""
        self.abinit_input.remove_variables(keys)

    def make_input(self):
        return str(self.abinit_input)


class HtcStrategy(AbstractStrategy):
    """
    Attributes:

        accuracy: Accuracy of the calculation used to define basic parameters of the run.
            such as tolerances, basis set truncation ...
    """
    __metaclass__ = abc.ABCMeta

    # Mapping runlevel --> optdriver variable
    _runl2optdriver = {
        "scf": 0,
        "nscf": 0,
        "relax": 0,
        "dfpt": 1,
        "screening": 3,
        "sigma": 4,
        "bse": 99,
    }

    # Name of the (default) tolerance used by the runlevels.
    _runl2tolname = {
        "scf": 'tolvrs',
        "nscf": 'tolwfr',
        "dfpt": 'toldfe',        # ?
        "screening": 'toldfe',   # dummy
        "sigma": 'toldfe',       # dummy
        "bse": 'toldfe',         # ?
        "relax": 'tolrff',
    }

    # Tolerances for the different levels of accuracy.
    T = collections.namedtuple('Tolerance', "low normal high")
    _tolerances = {
        "toldfe": T(1.e-7,  1.e-8,  1.e-9),
        "tolvrs": T(1.e-7,  1.e-8,  1.e-9),
        "tolwfr": T(1.e-15, 1.e-17, 1.e-19),
        "tolrff": T(0.04,   0.02,   0.01)}
    del T

    def __repr__(self):
        return "<%s at %s, accuracy = %s>" % (self.__class__.__name__, id(self), self.accuracy)

    @abc.abstractproperty
    def runlevel(self):
        """String defining the Runlevel. See _runl2optdriver."""

    @property
    def optdriver(self):
        """The optdriver associated to the calculation."""
        return self._runl2optdriver[self.runlevel]

    def learn(self, **data):
        """Update the data stored in self."""
        if not hasattr(self, "_data"):
            self._data = dict(data)
        else:
            if [k in self._data for k in data].count(True) != 0:
                raise ValueError("Keys %s are already present in data" % str([k for k in data]))
            self._data.update(data)

    @property
    def accuracy(self):
        """Accuracy used by the strategy."""
        try:
            return self._accuracy
        except AttributeError:
            self.set_accuracy("normal")
            return self._accuracy

    def set_accuracy(self, accuracy):
        """Accuracy setter."""
        if hasattr(self, "_accuracy"):
            raise RuntimeError("object already has accuracy %s " % self._accuracy)

        assert accuracy in ["low", "normal", "high"]
        self._accuracy = accuracy

    @property
    def data(self):
        """Data absorbed by the strategy during the workflow."""
        try:
            return self. _data
        except AttributeError:
            return {}

    @property
    def ecut(self):
        """Cutoff energy in Hartree."""
        try:
            # User option.
            return self.extra_abivars["ecut"]
        except KeyError:
            # Compute ecut from the Pseudo Hints.
            hints = [p.hint_for_accuracy(self.accuracy) for p in self.pseudos]
            return max(hint.ecut for hint in hints)

    @property
    def pawecutdg(self):
        """Cutoff energy in Hartree for the dense grid used in PAW calculations."""
        if not self.ispaw:
            return None

        try:
            # User option.
            return self.extra_abivars["pawecutdg"]
        except KeyError:
            raise NotImplementedError("")
            #ratio = max(p.suggested_augratio(accuracy) for p in self.pseudos])
            #ratio = augration_high if high else augratio_norm
            #pawecutdg = ecut * ratio

    @property
    def tolerance(self):
        """Return a dict {varname: varvalue} with the tolerance used for the calculation."""
        # Check user options first.
        for tolname in self._tolerances:
            try:
                return {tolname: self.extra_abivars[tolname]}
            except KeyError:
                pass

        # Use default values depending on the runlevel and the accuracy.
        tolname = self._runl2tolname[self.runlevel]

        return {tolname: getattr(self._tolerances[tolname], self.accuracy)}

    @property
    def need_forces(self):
        """True if forces are required at each SCF step (like the stresses)."""
        return self.runlevel in ["relax",]

    @property
    def need_stress(self):
        """True if the computation of the stress is required."""
        # TODO: here it's easier to check if optcell != 0
        return self.runlevel in ["relax"]

    def add_extra_abivars(self, abivars):
        """Add variables (dict) to extra_abivars."""
        self.extra_abivars.update(abivars)

    def remove_extra_abivars(self, keys):
        for key in keys:
            self.extra_abivars.pop(key)


class ScfStrategy(HtcStrategy):
    """
    Strategy for ground-state SCF calculations.
    """
    def __init__(self, structure, pseudos, ksampling, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, use_symmetries=True, **extra_abivars):
        """
        Args:
            structure: pymatgen structure
            pseudos: List of pseudopotentials.
            ksampling: :class:`Ksampling` object defining the sampling of the BZ.
            accuracy: Accuracy of the calculation.
            spin_mode: Spin polarization mode.
            smearing: string or :class:`Smearing` instance.
            charge: Total charge of the system. Default is 0.
            scf_algorithm: :class:`ElectronsAlgorithm` instance.
            use_symmetries: False if point group symmetries should not be used.
            extra_abivars: Extra variables that will be directly added to the input file.
        """
        super(ScfStrategy, self).__init__()

        self.set_accuracy(accuracy)
        self.structure = structure
        self.pseudos = select_pseudos(pseudos, structure)
        self.ksampling = ksampling
        self.use_symmetries = use_symmetries

        self.electrons = Electrons(spin_mode=spin_mode,
                                   smearing=smearing, algorithm=scf_algorithm,
                                   nband=None, fband=None, charge=charge)

        self.extra_abivars = extra_abivars

    @property
    def runlevel(self):
        return "scf"

    def make_input(self):
        extra = dict(optdriver=self.optdriver, ecut=self.ecut, pawecutdg=self.pawecutdg)
        extra.update(self.tolerance)
        extra.update({"nsym": 1 if not self.use_symmetries else None})

        extra.update(self.extra_abivars)

        inpw = InputWriter(self.structure, self.electrons, self.ksampling, **extra)
        return inpw.get_string()


class NscfStrategy(HtcStrategy):
    """
    Strategy for non-self-consistent calculations.
    """
    def __init__(self, scf_strategy, ksampling, nscf_nband, nscf_algorithm=None, **extra_abivars):
        """
        Args:
            scf_strategy: :class:`ScfStrategy` used for the GS run.
            ksampling: :class:`Ksampling` object defining the sampling of the BZ.
            nscf_nband: Number of bands to compute.
            nscf_algorithm :class:`ElectronsAlgorithm` instance.
            extra_abivars: Extra ABINIT variables that will be directly added to the input file
        """
        super(NscfStrategy, self).__init__()

        self.set_accuracy(scf_strategy.accuracy)
        self.scf_strategy = scf_strategy

        self.nscf_nband = nscf_nband
        self.pseudos = scf_strategy.pseudos
        self.ksampling = ksampling

        if nscf_algorithm is None:
            nscf_algorithm = {"iscf": -2}

        # Electrons used in the GS run.
        scf_electrons = scf_strategy.electrons

        self.electrons = Electrons(
            spin_mode=scf_electrons.spin_mode, smearing=scf_electrons.smearing,
            algorithm=nscf_algorithm, nband=nscf_nband,
            fband=None, charge=scf_electrons.charge, comment=None)

        self.extra_abivars = extra_abivars

    @property
    def runlevel(self):
        return "nscf"

    def make_input(self):
        # Initialize the system section from structure.
        scf_strategy = self.scf_strategy

        extra = dict(optdriver=self.optdriver, ecut=self.ecut, pawecutdg=self.pawecutdg)
        extra.update(self.tolerance)
        extra.update(self.extra_abivars)

        inp = InputWriter(scf_strategy.structure, self.electrons, self.ksampling, **extra)
        return inp.get_string()


class RelaxStrategy(ScfStrategy):
    """Extends ScfStrategy by adding an algorithm for the structural relaxation."""

    def __init__(self, structure, pseudos, ksampling, relax_algo, accuracy="normal", spin_mode="polarized",
                 smearing="fermi_dirac:0.1 eV", charge=0.0, scf_algorithm=None, **extra_abivars):
        """
        Args:
            structure: pymatgen structure
            pseudos: List of pseudopotentials.
            ksampling: :class:`Ksampling` object defining the sampling of the BZ.
            relax_algo: Object defining the algorithm for the structural relaxation.
            accuracy: Accuracy of the calculation.
            spin_mode: Flag defining the spin polarization. Defaults to "polarized"
            smearing: String or :class:`Smearing` instance.
            charge: Total charge of the system. Default is 0.
            scf_algorithm: :class:`ElectronsAlgorithm` instance.
            extra_abivars: Extra ABINIT variables that will be directly added to the input file
        """
        super(RelaxStrategy, self).__init__(
            structure, pseudos, ksampling,
            accuracy=accuracy, spin_mode=spin_mode, smearing=smearing,
            charge=charge, scf_algorithm=scf_algorithm, **extra_abivars)

        self.relax_algo = relax_algo

    @property
    def runlevel(self):
        return "relax"

    def make_input(self):
        # Input for the GS run
        input_str = super(RelaxStrategy, self).make_input()

        # Add the variables for the structural relaxation.
        input = InputWriter(self.relax_algo)
        input_str += input.get_string()

        return input_str


class ScreeningStrategy(HtcStrategy):
    """Strategy for Screening calculations."""
    def __init__(self, scf_strategy, nscf_strategy, screening, **extra_abivars):
        """
        Args:
            scf_strategy: :class:`ScfStrategy` used for the ground-state calculation
            nscf_strategy: :class:`NscStrategy` used for the non-self consistent calculation
            screening: :class:`Screening` instance
            extra_abivars: Extra ABINIT variables added directly to the input file
        """
        super(ScreeningStrategy, self).__init__()

        self.pseudos = scf_strategy.pseudos

        self.scf_strategy = scf_strategy
        self.nscf_strategy = nscf_strategy

        self.screening = screening

        scr_nband = screening.nband

        scf_electrons = scf_strategy.electrons
        nscf_electrons = nscf_strategy.electrons

        if scr_nband > nscf_electrons.nband:
            raise ValueError("Cannot use more that %d bands for the screening" % nscf_electrons.nband)

        self.ksampling = nscf_strategy.ksampling

        if not self.ksampling.is_homogeneous:
            raise ValueError("The k-sampling used for the NSCF run mush be homogeneous")

        self.electrons = Electrons(spin_mode=scf_electrons.spin_mode,
                                   smearing =scf_electrons.smearing,
                                   nband=scr_nband, charge=scf_electrons.charge, comment=None)

        self.extra_abivars = extra_abivars

    @property
    def runlevel(self):
        return "screening"

    def make_input(self):
        # FIXME
        extra = dict(optdriver=self.optdriver, ecut=self.ecut, ecutwfn=self.ecut, pawecutdg=self.pawecutdg)
        extra.update(self.tolerance)
        extra.update(self.extra_abivars)

        inpw = InputWriter(self.scf_strategy.structure, self.electrons, self.ksampling, self.screening, **extra)
        return inpw.get_string()


class SelfEnergyStrategy(HtcStrategy):
    """Strategy for self-energy calculations."""
    def __init__(self, scf_strategy, nscf_strategy, scr_strategy, sigma, **extra_abivars):
        """
        Args:
            scf_strategy: :class:`ScfStrategy` used for the ground-state calculation
            nscf_strategy: :class:`NscfStrategy` used for the non-self consistent calculation
            scr_strategy: :class:`ScrStrategy` used for the screening calculation
            sigma: :class:`SelfEnergy` instance.
            extra_abivars: Extra ABINIT variables added directly to the input file
        """
        # TODO Add consistency check between SCR and SIGMA strategies
        super(SelfEnergyStrategy, self).__init__()

        self.pseudos = scf_strategy.pseudos

        self.scf_strategy = scf_strategy
        self.nscf_strategy = nscf_strategy
        self.scr_strategy = scr_strategy

        self.sigma = sigma

        self.extra_abivars = extra_abivars

        scf_electrons = scf_strategy.electrons
        nscf_electrons = nscf_strategy.electrons

        if sigma.nband > nscf_electrons.nband:
            raise ValueError("Cannot use more that %d bands for the self-energy" % nscf_electrons.nband)

        self.ksampling = nscf_strategy.ksampling

        if not self.ksampling.is_homogeneous:
            raise ValueError("The k-sampling used for the NSCF run mush be homogeneous")

        self.electrons = Electrons(
            spin_mode=scf_electrons.spin_mode, smearing=scf_electrons.smearing,
            nband=sigma.nband, charge=scf_electrons.charge)

    @property
    def runlevel(self):
        return "sigma"

    def make_input(self):
        # FIXME
        extra = dict(optdriver=self.optdriver, ecut=self.ecut, ecutwfn=self.ecut, pawecutdg=self.pawecutdg)
        extra.update(self.tolerance)
        extra.update(self.extra_abivars)

        inpw = InputWriter(self.scf_strategy.structure, self.electrons, self.ksampling, self.sigma, **extra)
        return inpw.get_string()


class MdfBse_Strategy(HtcStrategy):
    """
    Strategy for Bethe-Salpeter calculation based on the
    model dielectric function and the scissors operator
    """
    def __init__(self, scf_strategy, nscf_strategy, exc_ham, **extra_abivars):
        """
        Args:
            scf_strategy: :class:`Strategy` used for the ground-state calculation.
            nscf_strategy: :class:`NscStrategy` used for the non-self consistent calculation.
            exc_ham: :class:`ExcitonicHamiltonian` instance.
            extra_abivars: Extra ABINIT variables added directly to the input file.
        """
        super(MdfBse_Strategy, self).__init__()

        self.pseudos = scf_strategy.pseudos

        self.scf_strategy = scf_strategy
        self.nscf_strategy = nscf_strategy

        self.exc_ham = exc_ham

        self.extra_abivars = extra_abivars

        scf_electrons = scf_strategy.electrons
        nscf_electrons = nscf_strategy.electrons

        if exc_ham.nband > nscf_electrons.nband:
            raise ValueError("Cannot use more that %d bands for the EXC hamiltonian." % nscf_electrons.nband)

        self.ksampling = nscf_strategy.ksampling

        if not self.ksampling.is_homogeneous:
            raise ValueError("The k-sampling used for the NSCF run mush be homogeneous")

        self.electrons = Electrons(
            spin_mode=scf_electrons.spin_mode, smearing=scf_electrons.smearing,
            nband=exc_ham.nband, charge=scf_electrons.charge)

    @property
    def runlevel(self):
        return "bse"

    def make_input(self):
        # FIXME
        extra = dict(optdriver=self.optdriver, ecut=self.ecut, pawecutdg=self.pawecutdg, ecutwfn=self.ecut)
        #extra.update(self.tolerance)
        extra.update(self.extra_abivars)

        inpw = InputWriter(self.scf_strategy.structure, self.electrons, self.ksampling, self.exc_ham, **extra)
        return inpw.get_string()


class InputWriter(object):
    """
    This object receives a list of `AbivarAble` objects, an optional
    dictionary with extra ABINIT variables and produces a (nicely formatted?) string with the input file.
    """
    MAX_SLEN = 100

    def __init__(self, *args, **kwargs):
        self.abiobj_dict = collections.OrderedDict()
        self.extra_abivars = collections.OrderedDict()

        for arg in args:
            if hasattr(arg, "to_abivars"):
                self.add_abiobj(arg)
            else:
                self.add_extra_abivars(arg)

        for k, v in kwargs.items():
            self.add_extra_abivars({k: v})

    def __str__(self):
        """String representation (the section of the abinit input file)."""
        return self.get_string()

    @property
    def abiobjects(self):
        """List of objects stored in self."""
        return self.abiobj_dict.values()

    def add_abiobj(self, obj):
        """Add the object obj to self."""
        if not hasattr(obj, "to_abivars"):
            raise ValueError("%s does not define the method to_abivars" % str(obj))

        cname = obj.__class__.__name__
        if cname in self.abiobj_dict:
            raise ValueError("%s is already stored" % cname)
        self.abiobj_dict[cname] = obj

    def add_extra_abivars(self, abivars):
        """Add variables (dict) to extra_abivars."""
        self.extra_abivars.update(abivars)

    def to_abivars(self):
        """Returns a dictionary with the abinit variables defined by the Card."""
        abivars = {}
        for obj in self.abiobjects:
            abivars.update(obj.to_abivars())

        abivars.update(self.extra_abivars)
        return abivars

    def print_abiobjects(self, stream=sys.stdout):
        lines = [str(obj) for obj in self.abiobjects]
        stream.write("\n".join(lines))

    @staticmethod
    def _format_kv(key, value):
        """Formatter"""
        if value is None:
            # Use ABINIT default.
            return []

        if isinstance(value, collections.Iterable) and not is_string(value):
            arr = np.array(value)
            if len(arr.shape) in [0,1]:
                # scalar or vector.
                token = [key, " ".join(str(i) for i in arr)]

            else:
                # array --> matrix
                matrix = np.reshape(arr, (-1, arr.shape[-1]))
                lines = []
                for idx, row in enumerate(matrix):
                    lines.append(" ".join(str(i) for i in row))
                token = [key + "\n", "\n".join(lines)]

        else:
            token = [key, str(value)]

        return token

    def _cut_lines(self, lines):
        MAX_SLEN = self.MAX_SLEN

        new_lines = []
        for line in lines:
            if len(line) > MAX_SLEN:
                #start, stop = 0, 0
                #while True:
                #    stop = start + MAX_SLEN
                #    if stop > len(line): break
                #    print(start, stop)
                #    if stop > len(line): stop = len(line)
                #    new_lines.append(line[start:stop])
                #    start = stop

                tokens = lines.split()
                cum_nchars, start = 0, 0
                for stop, tok in enumerate(tokens):
                    cum_nchars += len(tok) + 1

                    if cum_nchars > MAX_SLEN:
                        cum_nchars = 0
                        new_lines.append("".join(tokens[start:stop]))
                    else:
                        start = stop

                if cum_nchars:
                    new_lines.append("".join(tokens[start:stop]))

            else:
                new_lines.append(line)

        return new_lines

    def get_string(self, pretty=False):
        """
        Returns a string representation of self. The reason why this
        method is different from the __str__ method is to provide options for pretty printing.

        Args:
            pretty: Set to True for pretty aligned output.
        """
        lines = []
        app = lines.append

        # Write the Abinit objects first.
        for obj in self.abiobjects:
            #print(obj)
            app([80*"#", ""])
            app(["#", "%s" % obj.__class__.__name__])
            app([80*"#", ""])
            for (k, v) in obj.to_abivars().items():
                app(self._format_kv(k, v))

        # Extra variables.
        if self.extra_abivars:
            app([80*"#", ""])
            app(["#", "Extra_Abivars"])
            app([80*"#", ""])
            for (k, v) in self.extra_abivars.items():
                app(self._format_kv(k, v))

        #lines = self._cut_lines(lines)

        if pretty:
            return str_aligned(lines, header=None)
        else:
            return str_delimited(lines, header=None, delimiter=5*" ")


class OpticVar(collections.namedtuple("OpticVar", "name value help")):
    def __str__(self):
        sval = str(self.value)
        return (4*" ").join(sval, "!" + self.help)


class OpticInput(object):
    """
    abo_1WF7      ! Name of the first d/dk response wavefunction file, produced by abinit
    abo_1WF8      ! Name of the second d/dk response wavefunction file, produced by abinit
    abo_1WF9      ! Name of the third d/dk response wavefunction file, produced by abinit
    abo_WFK       ! Name of the ground-state wavefunction file, produced by abinit
    0.01          ! Value of the *smearing factor*, in Hartree
    0.010   1     ! frequency *step* and *maximum* frequency (Ha)
    0.000         ! *Scissor* shift if needed, in Hartree
    0.001         ! *Tolerance* on closeness of singularities (in Hartree)
    3             ! *Number of components* of linear optic tensor to be computed
    11 33 23      ! Linear *coefficients* to be computed (x=1, y=2, z=3)
    2             ! Number of components of nonlinear optic tensor to be computed
    123 222       ! Non-linear coefficients to be computed
    """

    # variable name --> default value.
    _VARIABLES = [
        OpticVar("ddkfile_x",       None, "Name of the first d/dk response wavefunction file"),
        OpticVar("ddkfile_y",       None, "Name of the second d/dk response wavefunction file"),
        OpticVar("ddkfile_z",       None, "Name of the third d/dk response wavefunction file"),
        OpticVar("wfkfile",         None, "Name of the ground-state wavefunction file"),
        OpticVar("zcut",            0.01, "Value of the *smearing factor*, in Hartree"),
        OpticVar("wmesh",     (0.010, 1), "Frequency *step* and *maximum* frequency (Ha)"),
        OpticVar("scissor",        0.000, "*Scissor* shift if needed, in Hartree"),
        OpticVar("sing_tol",       0.001, "*Tolerance* on closeness of singularities (in Hartree)"),
        OpticVar("num_lin_comp",    None, "*Number of components* of linear optic tensor to be computed"),
        OpticVar("lin_comp",        None, "Linear *coefficients* to be computed (x=1, y=2, z=3)"),
        OpticVar("num_nonlin_comp", None, "Number of components of nonlinear optic tensor to be computed"),
        OpticVar("nonlin_comp",     None, "! Non-linear coefficients to be computed"),
    ]

    _VARNAMES = [v.name for v in _VARIABLES]

    def __init__(self, **kwargs):
        # Default values
        self.vars = collections.OrderedDict((v.name, v.value) for v in _VARIABLES)

        # Update the variables with the values passed by the user
        for k, v in kwargs:
            if k not in self.VARNAMES:
                raise ValueError("varname %s not in %s" % (k, str(self.VARNAMES)))
            self.vars[k] = v

    def __init__(self, zcut, wstep, wmax, scissor, sing_tol, linear_components,
                 nonlinear_components=None, ddk_files=None, wfk=None):

        self.vars = vars = collections.OrderedDict(*self.VAR_NAMES)

        if ddk_files is not None:
            assert len(ddk_files) == 3
            assert wfk is not None
            for dir, ddk in zip(["x", "y", "z"], ddk_files):
                vars["ddkfile_" + dir] = os.path.abspath(ddk)

        if wfk is not None:
            vars["wfkfile"] = os.path.abspath(wfk)

        vars["zcut"] = zcut
        vars["wmesh"] = " ".join(map(str, (wstep, wmax)))
        vars["sing_tol"] = sing_tol

        vars["num_lin_comp"] = len(linear_components)
        vars["lin_comp"] = " ".join(str(c) for c in linear_components)

        vars["num_nonlin_comp"] = len(nonlinear_components)
        vars["nonlin_comp"] = " ".join(str(c) for c in nonlinear_components)

    def __init__(self, string):
        self.string = string

    def __str__(self):
        return self.string

    def to_string(self):
        lines = []
        app = lines.append

        for name in self.VARNAMES:
            var = self.vars[name]
            app(str(var))

        return "\n".join(lines)

    def make_input(self):
        return str(self)

    def add_extra_abivars(self, abivars):
        """
        Connection is explicit via the input file
        since we can pass the paths of the output files
        produced by the previous runs.
        """
