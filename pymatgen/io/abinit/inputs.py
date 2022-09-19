"""
This module defines a simplified interface for generating ABINIT input files.
Note that not all the features of Abinit are supported by BasicAbinitInput.
For a more comprehensive implementation, use the AbinitInput object provided by AbiPy.
"""

import abc
import copy
import json
import logging
import os
from collections import namedtuple
from collections.abc import Mapping, MutableMapping
from enum import Enum

import numpy as np
from monty.collections import AttrDict
from monty.json import MSONable
from monty.string import is_string, list_strings

from pymatgen.core.structure import Structure
from pymatgen.io.abinit import abiobjects as aobj
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from pymatgen.io.abinit.variable import InputVariable
from pymatgen.util.serialization import pmg_serialize

logger = logging.getLogger(__file__)


# List of Abinit variables used to specify the structure.
# This variables should not be passed to set_vars since
# they will be generated with structure.to_abivars()
GEOVARS = {
    "acell",
    "rprim",
    "rprimd",
    "angdeg",
    "xred",
    "xcart",
    "xangst",
    "znucl",
    "typat",
    "ntypat",
    "natom",
}

# Variables defining tolerances (used in pop_tolerances)
_TOLVARS = {
    "toldfe",
    "tolvrs",
    "tolwfr",
    "tolrff",
    "toldff",
    "tolimg",
    "tolmxf",
    "tolrde",
}

# Variables defining tolerances for the SCF cycle that are mutally exclusive
_TOLVARS_SCF = {
    "toldfe",
    "tolvrs",
    "tolwfr",
    "tolrff",
    "toldff",
}

# Variables determining if data files should be read in input
_IRDVARS = {
    "irdbseig",
    "irdbsreso",
    "irdhaydock",
    "irdddk",
    "irdden",
    "ird1den",
    "irdqps",
    "irdkss",
    "irdscr",
    "irdsuscep",
    "irdvdw",
    "irdwfk",
    "irdwfkfine",
    "irdwfq",
    "ird1wf",
}

# Name of the (default) tolerance used by the runlevels.
_runl2tolname = {
    "scf": "tolvrs",
    "nscf": "tolwfr",
    "dfpt": "toldfe",  # ?
    "screening": "toldfe",  # dummy
    "sigma": "toldfe",  # dummy
    "bse": "toldfe",  # ?
    "relax": "tolrff",
}

# Tolerances for the different levels of accuracy.

T = namedtuple("T", "low normal high")
_tolerances = {
    "toldfe": T(1.0e-7, 1.0e-8, 1.0e-9),
    "tolvrs": T(1.0e-7, 1.0e-8, 1.0e-9),
    "tolwfr": T(1.0e-15, 1.0e-17, 1.0e-19),
    "tolrff": T(0.04, 0.02, 0.01),
}
del T


# Default values used if user does not specify them
_DEFAULTS = dict(
    kppa=1000,
)


def as_structure(obj):
    """
    Convert obj into a Structure. Accepts:

        - Structure object.
        - Filename
        - Dictionaries (MSONable format or dictionaries with abinit variables).
    """
    if isinstance(obj, Structure):
        return obj

    if is_string(obj):
        return Structure.from_file(obj)

    if isinstance(obj, Mapping):
        if "@module" in obj:
            return Structure.from_dict(obj)
        return aobj.structure_from_abivars(cls=None, **obj)

    raise TypeError(f"Don't know how to convert {type(obj)} into a structure")


class ShiftMode(Enum):
    """
    Class defining the mode to be used for the shifts.
    G: Gamma centered
    M: Monkhorst-Pack ((0.5, 0.5, 0.5))
    S: Symmetric. Respects the chksymbreak with multiple shifts
    O: OneSymmetric. Respects the chksymbreak with a single shift (as in 'S' if a single shift is given, gamma
        centered otherwise.
    """

    GammaCentered = "G"
    MonkhorstPack = "M"
    Symmetric = "S"
    OneSymmetric = "O"

    @classmethod
    def from_object(cls, obj):
        """
        Returns an instance of ShiftMode based on the type of object passed. Converts strings to ShiftMode depending
        on the iniital letter of the string. G for GammaCenterd, M for MonkhorstPack,
        S for Symmetric, O for OneSymmetric.
        Case insensitive.
        """
        if isinstance(obj, cls):
            return obj
        if is_string(obj):
            return cls(obj[0].upper())
        raise TypeError(f"The object provided is not handled: type {type(obj)}")


def _stopping_criterion(runlevel, accuracy):
    """Return the stopping criterion for this runlevel with the given accuracy."""
    tolname = _runl2tolname[runlevel]
    return {tolname: getattr(_tolerances[tolname], accuracy)}


def _find_ecut_pawecutdg(ecut, pawecutdg, pseudos, accuracy):
    """Return a |AttrDict| with the value of ``ecut`` and ``pawecutdg``."""
    # Get ecut and pawecutdg from the pseudo hints.
    if ecut is None or (pawecutdg is None and any(p.ispaw for p in pseudos)):
        has_hints = all(p.has_hints for p in pseudos)

    if ecut is None:
        if has_hints:
            ecut = max(p.hint_for_accuracy(accuracy).ecut for p in pseudos)
        else:
            raise RuntimeError("ecut is None but pseudos do not provide hints for ecut")

    if pawecutdg is None and any(p.ispaw for p in pseudos):
        if has_hints:
            pawecutdg = max(p.hint_for_accuracy(accuracy).pawecutdg for p in pseudos)
        else:
            raise RuntimeError("pawecutdg is None but pseudos do not provide hints")

    return AttrDict(ecut=ecut, pawecutdg=pawecutdg)


def _find_scf_nband(structure, pseudos, electrons, spinat=None):
    """Find the value of ``nband``."""
    if electrons.nband is not None:
        return electrons.nband

    nsppol, smearing = electrons.nsppol, electrons.smearing

    # Number of valence electrons including possible extra charge
    nval = num_valence_electrons(structure, pseudos)
    nval -= electrons.charge

    # First guess (semiconductors)
    nband = nval // 2

    # TODO: Find better algorithm
    # If nband is too small we may kill the job, increase nband and restart
    # but this change could cause problems in the other steps of the calculation
    # if the change is not propagated e.g. phonons in metals.
    if smearing:
        # metallic occupation
        nband = max(np.ceil(nband * 1.2), nband + 10)
    else:
        nband = max(np.ceil(nband * 1.1), nband + 4)

    # Increase number of bands based on the starting magnetization
    if nsppol == 2 and spinat is not None:
        nband += np.ceil(max(np.sum(spinat, axis=0)) / 2.0)

    # Force even nband (easier to divide among procs, mandatory if nspinor == 2)
    nband += nband % 2
    return int(nband)


def _get_shifts(shift_mode, structure):
    """
    Gives the shifts based on the selected shift mode and on the symmetry of the structure.
    G: Gamma centered
    M: Monkhorst-Pack ((0.5, 0.5, 0.5))
    S: Symmetric. Respects the chksymbreak with multiple shifts
    O: OneSymmetric. Respects the chksymbreak with a single shift (as in 'S' if a single shift is given, gamma
        centered otherwise.

    Note: for some cases (e.g. body centered tetragonal), both the Symmetric and OneSymmetric may fail to satisfy the
        ``chksymbreak`` condition (Abinit input variable).
    """
    if shift_mode == ShiftMode.GammaCentered:
        return ((0, 0, 0),)
    if shift_mode == ShiftMode.MonkhorstPack:
        return ((0.5, 0.5, 0.5),)
    if shift_mode == ShiftMode.Symmetric:
        return calc_shiftk(structure)
    if shift_mode == ShiftMode.OneSymmetric:
        shifts = calc_shiftk(structure)
        if len(shifts) == 1:
            return shifts
        return ((0, 0, 0),)

    raise ValueError(f"invalid shift_mode: `{str(shift_mode)}`")


def gs_input(
    structure,
    pseudos,
    kppa=None,
    ecut=None,
    pawecutdg=None,
    scf_nband=None,
    accuracy="normal",
    spin_mode="polarized",
    smearing="fermi_dirac:0.1 eV",
    charge=0.0,
    scf_algorithm=None,
):
    """
    Returns a |BasicAbinitInput| for ground-state calculation.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the SCF run. Defaults to 1000 if not given.
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
                   according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
                   from the list of pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
    """
    multi = ebands_input(
        structure,
        pseudos,
        kppa=kppa,
        ndivsm=0,
        ecut=ecut,
        pawecutdg=pawecutdg,
        scf_nband=scf_nband,
        accuracy=accuracy,
        spin_mode=spin_mode,
        smearing=smearing,
        charge=charge,
        scf_algorithm=scf_algorithm,
    )

    return multi[0]


def ebands_input(
    structure,
    pseudos,
    kppa=None,
    nscf_nband=None,
    ndivsm=15,
    ecut=None,
    pawecutdg=None,
    scf_nband=None,
    accuracy="normal",
    spin_mode="polarized",
    smearing="fermi_dirac:0.1 eV",
    charge=0.0,
    scf_algorithm=None,
    dos_kppa=None,
):
    """
    Returns a |BasicMultiDataset| object for band structure calculations.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the SCF run. Defaults to 1000 if not given.
        nscf_nband: Number of bands included in the NSCF run. Set to scf_nband + 10 if None.
        ndivsm: Number of divisions used to sample the smallest segment of the k-path.
                if 0, only the GS input is returned in multi[0].
        ecut: cutoff energy in Ha (if None, ecut is initialized from the pseudos according to accuracy)
        pawecutdg: cutoff energy in Ha for PAW double-grid (if None, pawecutdg is initialized from the pseudos
            according to accuracy)
        scf_nband: Number of bands for SCF run. If scf_nband is None, nband is automatically initialized
            from the list of pseudos, the structure and the smearing option.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for solving of the SCF cycle.
        dos_kppa: Scalar or List of integers with the number of k-points per atom
            to be used for the computation of the DOS (None if DOS is not wanted).
    """
    structure = as_structure(structure)

    if dos_kppa is not None and not isinstance(dos_kppa, (list, tuple)):
        dos_kppa = [dos_kppa]

    multi = BasicMultiDataset(structure, pseudos, ndtset=2 if dos_kppa is None else 2 + len(dos_kppa))

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy))

    # SCF calculation.
    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa
    scf_ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0)
    scf_electrons = aobj.Electrons(
        spin_mode=spin_mode,
        smearing=smearing,
        algorithm=scf_algorithm,
        charge=charge,
        nband=scf_nband,
        fband=None,
    )

    if scf_electrons.nband is None:
        scf_electrons.nband = _find_scf_nband(structure, multi.pseudos, scf_electrons, multi[0].get("spinat", None))

    multi[0].set_vars(scf_ksampling.to_abivars())
    multi[0].set_vars(scf_electrons.to_abivars())
    multi[0].set_vars(_stopping_criterion("scf", accuracy))
    if ndivsm == 0:
        return multi

    # Band structure calculation.
    nscf_ksampling = aobj.KSampling.path_from_structure(ndivsm, structure)
    nscf_nband = scf_electrons.nband + 10 if nscf_nband is None else nscf_nband
    nscf_electrons = aobj.Electrons(
        spin_mode=spin_mode,
        smearing=smearing,
        algorithm={"iscf": -2},
        charge=charge,
        nband=nscf_nband,
        fband=None,
    )

    multi[1].set_vars(nscf_ksampling.to_abivars())
    multi[1].set_vars(nscf_electrons.to_abivars())
    multi[1].set_vars(_stopping_criterion("nscf", accuracy))

    # DOS calculation with different values of kppa.
    if dos_kppa is not None:
        for i, kppa_ in enumerate(dos_kppa):
            dos_ksampling = aobj.KSampling.automatic_density(structure, kppa_, chksymbreak=0)
            # dos_ksampling = aobj.KSampling.monkhorst(dos_ngkpt, shiftk=dos_shiftk, chksymbreak=0)
            dos_electrons = aobj.Electrons(
                spin_mode=spin_mode,
                smearing=smearing,
                algorithm={"iscf": -2},
                charge=charge,
                nband=nscf_nband,
            )
            dt = 2 + i
            multi[dt].set_vars(dos_ksampling.to_abivars())
            multi[dt].set_vars(dos_electrons.to_abivars())
            multi[dt].set_vars(_stopping_criterion("nscf", accuracy))

    return multi


def ion_ioncell_relax_input(
    structure,
    pseudos,
    kppa=None,
    nband=None,
    ecut=None,
    pawecutdg=None,
    accuracy="normal",
    spin_mode="polarized",
    smearing="fermi_dirac:0.1 eV",
    charge=0.0,
    scf_algorithm=None,
    shift_mode="Monkhorst-pack",
):
    """
    Returns a |BasicMultiDataset| for a structural relaxation. The first dataset optmizes the
    atomic positions at fixed unit cell. The second datasets optimizes both ions and unit cell parameters.

    Args:
        structure: |Structure| object.
        pseudos: List of filenames or list of |Pseudo| objects or |PseudoTable| object.
        kppa: Defines the sampling used for the Brillouin zone.
        nband: Number of bands included in the SCF run.
        accuracy: Accuracy of the calculation.
        spin_mode: Spin polarization.
        smearing: Smearing technique.
        charge: Electronic charge added to the unit cell.
        scf_algorithm: Algorithm used for the solution of the SCF cycle.
    """
    structure = as_structure(structure)
    multi = BasicMultiDataset(structure, pseudos, ndtset=2)

    # Set the cutoff energies.
    multi.set_vars(_find_ecut_pawecutdg(ecut, pawecutdg, multi.pseudos, accuracy))

    kppa = _DEFAULTS.get("kppa") if kppa is None else kppa

    shift_mode = ShiftMode.from_object(shift_mode)
    shifts = _get_shifts(shift_mode, structure)
    ksampling = aobj.KSampling.automatic_density(structure, kppa, chksymbreak=0, shifts=shifts)
    electrons = aobj.Electrons(
        spin_mode=spin_mode,
        smearing=smearing,
        algorithm=scf_algorithm,
        charge=charge,
        nband=nband,
        fband=None,
    )

    if electrons.nband is None:
        electrons.nband = _find_scf_nband(structure, multi.pseudos, electrons, multi[0].get("spinat", None))

    ion_relax = aobj.RelaxationMethod.atoms_only(atoms_constraints=None)
    ioncell_relax = aobj.RelaxationMethod.atoms_and_cell(atoms_constraints=None)

    multi.set_vars(electrons.to_abivars())
    multi.set_vars(ksampling.to_abivars())

    multi[0].set_vars(ion_relax.to_abivars())
    multi[0].set_vars(_stopping_criterion("relax", accuracy))

    multi[1].set_vars(ioncell_relax.to_abivars())
    multi[1].set_vars(_stopping_criterion("relax", accuracy))

    return multi


def calc_shiftk(structure, symprec: float = 0.01, angle_tolerance=5):
    """
    Find the values of ``shiftk`` and ``nshiftk`` appropriated for the sampling of the Brillouin zone.

    When the primitive vectors of the lattice do NOT form a FCC or a BCC lattice,
    the usual (shifted) Monkhorst-Pack grids are formed by using nshiftk=1 and shiftk 0.5 0.5 0.5 .
    This is often the preferred k point sampling. For a non-shifted Monkhorst-Pack grid,
    use `nshiftk=1` and `shiftk 0.0 0.0 0.0`, but there is little reason to do that.

    When the primitive vectors of the lattice form a FCC lattice, with rprim::

            0.0 0.5 0.5
            0.5 0.0 0.5
            0.5 0.5 0.0

    the (very efficient) usual Monkhorst-Pack sampling will be generated by using nshiftk= 4 and shiftk::

        0.5 0.5 0.5
        0.5 0.0 0.0
        0.0 0.5 0.0
        0.0 0.0 0.5

    When the primitive vectors of the lattice form a BCC lattice, with rprim::

           -0.5  0.5  0.5
            0.5 -0.5  0.5
            0.5  0.5 -0.5

    the usual Monkhorst-Pack sampling will be generated by using nshiftk= 2 and shiftk::

            0.25  0.25  0.25
           -0.25 -0.25 -0.25

    However, the simple sampling nshiftk=1 and shiftk 0.5 0.5 0.5 is excellent.

    For hexagonal lattices with hexagonal axes, e.g. rprim::

            1.0  0.0       0.0
           -0.5  sqrt(3)/2 0.0
            0.0  0.0       1.0

    one can use nshiftk= 1 and shiftk 0.0 0.0 0.5
    In rhombohedral axes, e.g. using angdeg 3*60., this corresponds to shiftk 0.5 0.5 0.5,
    to keep the shift along the symmetry axis.

    Returns:
        Suggested value of shiftk.
    """
    # Find lattice type.
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sym = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=angle_tolerance)
    lattice_type, spg_symbol = sym.get_lattice_type(), sym.get_space_group_symbol()

    # Check if the cell is primitive
    is_primitive = len(sym.find_primitive()) == len(structure)

    # Generate the appropriate set of shifts.
    shiftk = None

    if is_primitive:
        if lattice_type == "cubic":
            if "F" in spg_symbol:
                # FCC
                shiftk = [0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5]

            elif "I" in spg_symbol:
                # BCC
                shiftk = [0.25, 0.25, 0.25, -0.25, -0.25, -0.25]
                # shiftk = [0.5, 0.5, 05])

        elif lattice_type == "hexagonal":
            # Find the hexagonal axis and set the shift along it.
            for i, angle in enumerate(structure.lattice.angles):
                if abs(angle - 120) < 1.0:
                    j = (i + 1) % 3
                    k = (i + 2) % 3
                    hex_ax = [ax for ax in range(3) if ax not in [j, k]][0]
                    break
            else:
                raise ValueError("Cannot find hexagonal axis")

            shiftk = [0.0, 0.0, 0.0]
            shiftk[hex_ax] = 0.5

        elif lattice_type == "tetragonal":
            if "I" in spg_symbol:
                # BCT
                shiftk = [0.25, 0.25, 0.25, -0.25, -0.25, -0.25]

    if shiftk is None:
        # Use default value.
        shiftk = [0.5, 0.5, 0.5]

    return np.reshape(shiftk, (-1, 3))


def num_valence_electrons(structure, pseudos):
    """
    Returns the number of valence electrons.

    Args:
        pseudos: List of |Pseudo| objects or list of filenames.
    """
    nval, table = 0, PseudoTable.as_table(pseudos)
    for site in structure:
        pseudo = table.pseudo_with_symbol(site.specie.symbol)
        nval += pseudo.Z_val

    return int(nval) if int(nval) == nval else nval


class AbstractInput(MutableMapping, metaclass=abc.ABCMeta):
    """
    Abstract class defining the methods that must be implemented by Input objects.
    """

    # ABC protocol: __delitem__, __getitem__, __iter__, __len__, __setitem__
    def __delitem__(self, key):
        return self.vars.__delitem__(key)

    def __getitem__(self, key):
        return self.vars.__getitem__(key)

    def __iter__(self):
        return self.vars.__iter__()

    def __len__(self):
        return len(self.vars)

    def __setitem__(self, key, value):
        self._check_varname(key)
        return self.vars.__setitem__(key, value)

    def __repr__(self):
        return f"<{type(self).__name__} at {id(self)}>"

    def __str__(self):
        return self.to_string()

    def write(self, filepath="run.abi"):
        """
        Write the input file to file to ``filepath``.
        """
        dirname = os.path.dirname(os.path.abspath(filepath))
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        # Write the input file.
        with open(filepath, "wt") as fh:
            fh.write(str(self))

    def deepcopy(self):
        """Deep copy of the input."""
        return copy.deepcopy(self)

    def set_vars(self, *args, **kwargs):
        """
        Set the value of the variables.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        for varname, varvalue in kwargs.items():
            self[varname] = varvalue
        return kwargs

    def set_vars_ifnotin(self, *args, **kwargs):
        """
        Set the value of the variables but only if the variable is not already present.
        Return dict with the variables added to the input.

        Example:

            input.set_vars(ecut=10, ionmov=3)
        """
        kwargs.update(dict(*args))
        added = {}
        for varname, varvalue in kwargs.items():
            if varname not in self:
                self[varname] = varvalue
                added[varname] = varvalue
        return added

    def pop_vars(self, keys):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.
        Unlike remove_vars, no exception is raised if the variables are not in the input.

        Args:
            keys: string or list of strings with variable names.

        Example:
            inp.pop_vars(["ionmov", "optcell", "ntime", "dilatmx"])
        """
        return self.remove_vars(keys, strict=False)

    def remove_vars(self, keys, strict=True):
        """
        Remove the variables listed in keys.
        Return dictionary with the variables that have been removed.

        Args:
            keys: string or list of strings with variable names.
            strict: If True, KeyError is raised if at least one variable is not present.
        """
        removed = {}
        for key in list_strings(keys):
            if strict and key not in self:
                raise KeyError(f"key: {key} not in self:\n {list(self)}")
            if key in self:
                removed[key] = self.pop(key)

        return removed

    @abc.abstractproperty
    def vars(self):
        """Dictionary with the input variables. Used to implement dict-like interface."""

    @abc.abstractmethod
    def _check_varname(self, key):
        """Check if key is a valid name. Raise self.Error if not valid."""

    @abc.abstractmethod
    def to_string(self):
        """Returns a string with the input."""


class BasicAbinitInputError(Exception):
    """Base error class for exceptions raised by ``BasicAbinitInput``."""


class BasicAbinitInput(AbstractInput, MSONable):
    """
    This object stores the ABINIT variables for a single dataset.
    """

    Error = BasicAbinitInputError

    def __init__(
        self,
        structure,
        pseudos,
        pseudo_dir=None,
        comment=None,
        abi_args=None,
        abi_kwargs=None,
    ):
        """
        Args:
            structure: Parameters defining the crystalline structure. Accepts |Structure| object
            file with structure (CIF, netcdf file, ...) or dictionary with ABINIT geo variables.
            pseudos: Pseudopotentials to be used for the calculation. Accepts: string or list of strings
                with the name of the pseudopotential files, list of |Pseudo| objects
                or |PseudoTable| object.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
            comment: Optional string with a comment that will be placed at the beginning of the file.
            abi_args: list of tuples (key, value) with the initial set of variables. Default: Empty
            abi_kwargs: Dictionary with the initial set of variables. Default: Empty
        """
        # Internal dict with variables. we use an ordered dict so that
        # variables will be likely grouped by `topics` when we fill the input.
        abi_args = abi_args or []
        for key, _value in abi_args:
            self._check_varname(key)

        abi_kwargs = {} if abi_kwargs is None else abi_kwargs
        for key in abi_kwargs:
            self._check_varname(key)

        args = list(abi_args)[:]
        args.extend(list(abi_kwargs.items()))

        self._vars = dict(args)
        self.set_structure(structure)

        if pseudo_dir is not None:
            pseudo_dir = os.path.abspath(pseudo_dir)
            if not os.path.exists(pseudo_dir):
                raise self.Error(f"Directory {pseudo_dir} does not exist")
            pseudos = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

        try:
            self._pseudos = PseudoTable.as_table(pseudos).get_pseudos_for_structure(self.structure)
        except ValueError as exc:
            raise self.Error(str(exc))

        if comment is not None:
            self.set_comment(comment)

    @pmg_serialize
    def as_dict(self):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        # Use a list of (key, value) to serialize the dict
        abi_args = []
        for key, value in self.items():
            if isinstance(value, np.ndarray):
                value = value.tolist()
            abi_args.append((key, value))

        return dict(
            structure=self.structure.as_dict(),
            pseudos=[p.as_dict() for p in self.pseudos],
            comment=self.comment,
            abi_args=abi_args,
        )

    @property
    def vars(self):
        """Dictionary with variables."""
        return self._vars

    @classmethod
    def from_dict(cls, d):
        """
        JSON interface used in pymatgen for easier serialization.
        """
        pseudos = [Pseudo.from_file(p["filepath"]) for p in d["pseudos"]]
        return cls(d["structure"], pseudos, comment=d["comment"], abi_args=d["abi_args"])

    def add_abiobjects(self, *abi_objects):
        """
        This function receive a list of ``AbiVarable`` objects and add
        the corresponding variables to the input.
        """
        d = {}
        for obj in abi_objects:
            if not hasattr(obj, "to_abivars"):
                raise TypeError(f"type {type(obj)}: {repr(obj)} does not have `to_abivars` method")
            d.update(self.set_vars(obj.to_abivars()))
        return d

    def __setitem__(self, key, value):
        if key in _TOLVARS_SCF and hasattr(self, "_vars") and any(t in self._vars and t != key for t in _TOLVARS_SCF):
            logger.info(f"Replacing previously set tolerance variable: {self.remove_vars(_TOLVARS_SCF, strict=False)}.")

        return super().__setitem__(key, value)

    def _check_varname(self, key):
        if key in GEOVARS:
            raise self.Error(
                "You cannot set the value of a variable associated to the structure.\n"
                "Use Structure objects to prepare the input file."
            )

    def to_string(self, post=None, with_structure=True, with_pseudos=True, exclude=None):
        """
        String representation.

        Args:
            post: String that will be appended to the name of the variables
                Note that post is usually autodetected when we have multiple datatasets
                It is mainly used when we have an input file with a single dataset
                so that we can prevent the code from adding "1" to the name of the variables
                (In this case, indeed, Abinit complains if ndtset=1 is not specified
                and we don't want ndtset=1 simply because the code will start to add
                _DS1_ to all the input and output files.
            with_structure: False if section with structure variables should not be printed.
            with_pseudos: False if JSON section with pseudo data should not be added.
            exclude: List of variable names that should be ignored.
        """
        lines = []
        app = lines.append

        if self.comment:
            app("# " + self.comment.replace("\n", "\n#"))

        post = post if post is not None else ""
        exclude = set(exclude) if exclude is not None else set()

        # Default is no sorting else alphabetical order.
        keys = sorted(k for k, v in self.items() if k not in exclude and v is not None)

        # Extract the items from the dict and add the geo variables at the end
        items = [(k, self[k]) for k in keys]
        if with_structure:
            items.extend(list(aobj.structure_to_abivars(self.structure).items()))

        for name, value in items:
            # Build variable, convert to string and append it
            vname = name + post
            app(str(InputVariable(vname, value)))

        s = "\n".join(lines)
        if not with_pseudos:
            return s

        # Add JSON section with pseudo potentials.
        ppinfo = ["\n\n\n#<JSON>"]
        d = {"pseudos": [p.as_dict() for p in self.pseudos]}
        ppinfo.extend(json.dumps(d, indent=4).splitlines())
        ppinfo.append("</JSON>")

        s += "\n#".join(ppinfo)
        return s

    @property
    def comment(self):
        """Optional string with comment. None if comment is not set."""
        try:
            return self._comment
        except AttributeError:
            return None

    def set_comment(self, comment):
        """Set a comment to be included at the top of the file."""
        self._comment = comment

    @property
    def structure(self):
        """The |Structure| object associated to this input."""
        return self._structure

    def set_structure(self, structure):
        """Set structure."""
        self._structure = as_structure(structure)

        # Check volume
        m = self.structure.lattice.matrix
        if np.dot(np.cross(m[0], m[1]), m[2]) <= 0:
            raise self.Error("The triple product of the lattice vector is negative. Use structure.abi_sanitize.")

        return self._structure

    # Helper functions to facilitate the specification of several variables.
    def set_kmesh(self, ngkpt, shiftk, kptopt=1):
        """
        Set the variables for the sampling of the BZ.

        Args:
            ngkpt: Monkhorst-Pack divisions
            shiftk: List of shifts.
            kptopt: Option for the generation of the mesh.
        """
        shiftk = np.reshape(shiftk, (-1, 3))
        return self.set_vars(ngkpt=ngkpt, kptopt=kptopt, nshiftk=len(shiftk), shiftk=shiftk)

    def set_gamma_sampling(self):
        """Gamma-only sampling of the BZ."""
        return self.set_kmesh(ngkpt=(1, 1, 1), shiftk=(0, 0, 0))

    def set_kpath(self, ndivsm, kptbounds=None, iscf=-2):
        """
        Set the variables for the computation of the electronic band structure.

        Args:
            ndivsm: Number of divisions for the smallest segment.
            kptbounds: k-points defining the path in k-space.
                If None, we use the default high-symmetry k-path defined in the pymatgen database.
        """

        if kptbounds is None:
            from pymatgen.symmetry.bandstructure import HighSymmKpath

            hsym_kpath = HighSymmKpath(self.structure)

            name2frac_coords = hsym_kpath.kpath["kpoints"]
            kpath = hsym_kpath.kpath["path"]

            frac_coords, names = [], []
            for segment in kpath:
                for name in segment:
                    fc = name2frac_coords[name]
                    frac_coords.append(fc)
                    names.append(name)
            kptbounds = np.array(frac_coords)

        kptbounds = np.reshape(kptbounds, (-1, 3))
        # self.pop_vars(["ngkpt", "shiftk"]) ??

        return self.set_vars(kptbounds=kptbounds, kptopt=-(len(kptbounds) - 1), ndivsm=ndivsm, iscf=iscf)

    def set_spin_mode(self, spin_mode):
        """
        Set the variables used to the treat the spin degree of freedom.
        Return dictionary with the variables that have been removed.

        Args:
            spin_mode: :class:`SpinMode` object or string. Possible values for string are:

            - polarized
            - unpolarized
            - afm (anti-ferromagnetic)
            - spinor (non-collinear magnetism)
            - spinor_nomag (non-collinear, no magnetism)
        """
        # Remove all variables used to treat spin
        old_vars = self.pop_vars(["nsppol", "nspden", "nspinor"])
        self.add_abiobjects(aobj.SpinMode.as_spinmode(spin_mode))
        return old_vars

    @property
    def pseudos(self):
        """List of |Pseudo| objects."""
        return self._pseudos

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    def new_with_vars(self, *args, **kwargs):
        """
        Return a new input with the given variables.

        Example:
            new = input.new_with_vars(ecut=20)
        """
        # Avoid modifications in self.
        new = self.deepcopy()
        new.set_vars(*args, **kwargs)
        return new

    def pop_tolerances(self):
        """
        Remove all the tolerance variables present in self.
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_TOLVARS, strict=False)

    def pop_irdvars(self):
        """
        Remove all the `ird*` variables present in self.
        Return dictionary with the variables that have been removed.
        """
        return self.remove_vars(_IRDVARS, strict=False)


class BasicMultiDataset:
    """
    This object is essentially a list of BasicAbinitInput objects.
    that provides an easy-to-use interface to apply global changes to the
    the inputs stored in the objects.

    Let's assume for example that multi contains two ``BasicAbinitInput`` objects and we
    want to set `ecut` to 1 in both dictionaries. The direct approach would be:

        for inp in multi:
            inp.set_vars(ecut=1)

    or alternatively:

        for i in range(multi.ndtset):
            multi[i].set_vars(ecut=1)

    BasicMultiDataset provides its own implementation of __getattr__ so that one can simply use:

        multi.set_vars(ecut=1)

        multi.get("ecut") returns a list of values. It's equivalent to:

            [inp["ecut"] for inp in multi]

        Note that if "ecut" is not present in one of the input of multi, the corresponding entry is set to None.
        A default value can be specified with:

            multi.get("paral_kgb", 0)

    .. warning::

        BasicMultiDataset does not support calculations done with different sets of pseudopotentials.
        The inputs can have different crystalline structures (as long as the atom types are equal)
        but each input in BasicMultiDataset must have the same set of pseudopotentials.
    """

    Error = BasicAbinitInputError

    @classmethod
    def from_inputs(cls, inputs):
        """Build object from a list of BasicAbinitInput objects."""
        for inp in inputs:
            if any(p1 != p2 for p1, p2 in zip(inputs[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")

        # Build BasicMultiDataset from input structures and pseudos and add inputs.
        multi = cls(
            structure=[inp.structure for inp in inputs],
            pseudos=inputs[0].pseudos,
            ndtset=len(inputs),
        )

        # Add variables
        for inp, new_inp in zip(inputs, multi):
            new_inp.set_vars(**inp)

        return multi

    @classmethod
    def replicate_input(cls, input, ndtset):
        """Construct a multidataset with ndtset from the BasicAbinitInput input."""
        multi = cls(input.structure, input.pseudos, ndtset=ndtset)

        for inp in multi:
            inp.set_vars(**input)

        return multi

    def __init__(self, structure: Structure, pseudos, pseudo_dir="", ndtset=1):
        """
        Args:
            structure: file with the structure, |Structure| object or dictionary with ABINIT geo variable
                Accepts also list of objects that can be converted to Structure object.
                In this case, however, ndtset must be equal to the length of the list.
            pseudos: String or list of string with the name of the pseudopotential files.
            pseudo_dir: Name of the directory where the pseudopotential files are located.
            ndtset: Number of datasets.
        """
        # Setup of the pseudopotential files.
        if isinstance(pseudos, Pseudo):
            pseudos = [pseudos]

        elif isinstance(pseudos, PseudoTable):
            pseudos = pseudos

        elif all(isinstance(p, Pseudo) for p in pseudos):
            pseudos = PseudoTable(pseudos)

        else:
            # String(s)
            pseudo_dir = os.path.abspath(pseudo_dir)
            pseudo_paths = [os.path.join(pseudo_dir, p) for p in list_strings(pseudos)]

            missing = [p for p in pseudo_paths if not os.path.exists(p)]
            if missing:
                raise self.Error(f"Cannot find the following pseudopotential files:\n{str(missing)}")

            pseudos = PseudoTable(pseudo_paths)

        # Build the list of BasicAbinitInput objects.
        if ndtset <= 0:
            raise ValueError(f"ndtset {ndtset} cannot be <=0")

        if not isinstance(structure, (list, tuple)):
            self._inputs = [BasicAbinitInput(structure=structure, pseudos=pseudos) for i in range(ndtset)]
        else:
            assert len(structure) == ndtset
            self._inputs = [BasicAbinitInput(structure=s, pseudos=pseudos) for s in structure]

    @property
    def ndtset(self):
        """Number of inputs in self."""
        return len(self)

    @property
    def pseudos(self):
        """Pseudopotential objects."""
        return self[0].pseudos

    @property
    def ispaw(self):
        """True if PAW calculation."""
        return all(p.ispaw for p in self.pseudos)

    @property
    def isnc(self):
        """True if norm-conserving calculation."""
        return all(p.isnc for p in self.pseudos)

    def __len__(self):
        return len(self._inputs)

    def __getitem__(self, key):
        return self._inputs[key]

    def __iter__(self):
        return self._inputs.__iter__()

    def __getattr__(self, name):
        _inputs = object.__getattribute__(self, "_inputs")
        m = getattr(_inputs[0], name)
        if m is None:
            raise AttributeError(
                f"Cannot find attribute {type(self).__name__}. Tried in {name} and then in BasicAbinitInput object"
            )
        isattr = not callable(m)

        def on_all(*args, **kwargs):
            results = []
            for obj in self._inputs:
                a = getattr(obj, name)
                # print("name", name, ", type:", type(a), "callable: ",callable(a))
                if callable(a):
                    results.append(a(*args, **kwargs))
                else:
                    results.append(a)

            return results

        if isattr:
            on_all = on_all()

        return on_all

    def __add__(self, other):
        """self + other"""
        if isinstance(other, BasicAbinitInput):
            new_mds = BasicMultiDataset.from_inputs(self)
            new_mds.append(other)
            return new_mds
        if isinstance(other, BasicMultiDataset):
            new_mds = BasicMultiDataset.from_inputs(self)
            new_mds.extend(other)
            return new_mds

        raise NotImplementedError("Operation not supported")

    def __radd__(self, other):
        if isinstance(other, BasicAbinitInput):
            new_mds = BasicMultiDataset.from_inputs([other])
            new_mds.extend(self)
        elif isinstance(other, BasicMultiDataset):
            new_mds = BasicMultiDataset.from_inputs(other)
            new_mds.extend(self)
        else:
            raise NotImplementedError("Operation not supported")

    def append(self, abinit_input):
        """Add a |BasicAbinitInput| to the list."""
        assert isinstance(abinit_input, BasicAbinitInput)
        if any(p1 != p2 for p1, p2 in zip(abinit_input.pseudos, abinit_input.pseudos)):
            raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.append(abinit_input)

    def extend(self, abinit_inputs):
        """Extends self with a list of |BasicAbinitInput| objects."""
        assert all(isinstance(inp, BasicAbinitInput) for inp in abinit_inputs)
        for inp in abinit_inputs:
            if any(p1 != p2 for p1, p2 in zip(self[0].pseudos, inp.pseudos)):
                raise ValueError("Pseudos must be consistent when from_inputs is invoked.")
        self._inputs.extend(abinit_inputs)

    def addnew_from(self, dtindex):
        """Add a new entry in the multidataset by copying the input with index ``dtindex``."""
        self.append(self[dtindex].deepcopy())

    def split_datasets(self):
        """Return list of |BasicAbinitInput| objects.."""
        return self._inputs

    def deepcopy(self):
        """Deep copy of the BasicMultiDataset."""
        return copy.deepcopy(self)

    @property
    def has_same_structures(self):
        """True if all inputs in BasicMultiDataset are equal."""
        return all(self[0].structure == inp.structure for inp in self)

    def __str__(self):
        return self.to_string()

    def to_string(self, with_pseudos=True):
        """
        String representation i.e. the input file read by Abinit.

        Args:
            with_pseudos: False if JSON section with pseudo data should not be added.
        """
        if self.ndtset > 1:
            # Multi dataset mode.
            lines = [f"ndtset {int(self.ndtset)}"]

            def has_same_variable(kref, vref, other_inp):
                """True if variable kref is present in other_inp with the same value."""
                if kref not in other_inp:
                    return False
                otherv = other_inp[kref]
                return np.array_equal(vref, otherv)

            # Don't repeat variable that are common to the different datasets.
            # Put them in the `Global Variables` section and exclude these variables in inp.to_string
            global_vars = set()
            for k0, v0 in self[0].items():
                isame = True
                for i in range(1, self.ndtset):
                    isame = has_same_variable(k0, v0, self[i])
                    if not isame:
                        break
                if isame:
                    global_vars.add(k0)
            # print("global_vars vars", global_vars)

            w = 92
            if global_vars:
                lines.append(w * "#")
                lines.append("### Global Variables.")
                lines.append(w * "#")
                for key in global_vars:
                    vname = key
                    lines.append(str(InputVariable(vname, self[0][key])))

            has_same_structures = self.has_same_structures
            if has_same_structures:
                # Write structure here and disable structure output in input.to_string
                lines.append(w * "#")
                lines.append("#" + ("STRUCTURE").center(w - 1))
                lines.append(w * "#")
                for key, value in aobj.structure_to_abivars(self[0].structure).items():
                    vname = key
                    lines.append(str(InputVariable(vname, value)))

            for i, inp in enumerate(self):
                header = f"### DATASET {i + 1} ###"
                is_last = i == self.ndtset - 1
                s = inp.to_string(
                    post=str(i + 1),
                    with_pseudos=is_last and with_pseudos,
                    with_structure=not has_same_structures,
                    exclude=global_vars,
                )
                if s:
                    s = f"\n{len(header) * '#'}\n{header}\n{len(header) * '#'}\n{s}\n"

                lines.append(s)

            return "\n".join(lines)

        # single datasets ==> don't append the dataset index to the variables.
        # this trick is needed because Abinit complains if ndtset is not specified
        # and we have variables that end with the dataset index e.g. acell1
        # We don't want to specify ndtset here since abinit will start to add DS# to
        # the input and output files thus complicating the algorithms we have to use to locate the files.
        return self[0].to_string(with_pseudos=with_pseudos)

    def write(self, filepath="run.abi"):
        """
        Write ``ndset`` input files to disk. The name of the file
        is constructed from the dataset index e.g. run0.abi
        """
        root, ext = os.path.splitext(filepath)
        for i, inp in enumerate(self):
            p = root + f"DS{i}" + ext
            inp.write(filepath=p)
