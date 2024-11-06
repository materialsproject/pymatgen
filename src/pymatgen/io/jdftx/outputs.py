"""JDFTx outputs parsing module.

Module for parsing outputs of JDFTx.

Note: JDFTXOutfile will be moved back to its own module once a more broad outputs
class is written.

@mkhorton - this file is ready to review
"""

from __future__ import annotations

import inspect
import pprint
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

from pymatgen.io.jdftx.jdftxoutfileslice import JDFTXOutfileSlice
from pymatgen.io.jdftx.utils import read_outfile_slices

if TYPE_CHECKING:
    from pathlib import Path

    import numpy as np

    from pymatgen.io.jdftx.jelstep import JElSteps
    from pymatgen.io.jdftx.jminsettings import (
        JMinSettingsElectronic,
        JMinSettingsFluid,
        JMinSettingsIonic,
        JMinSettingsLattice,
    )
    from pymatgen.io.jdftx.joutstructures import JOutStructures

__author__ = "Ben Rich, Jacob Clary"


@dataclass
class JDFTXOutfile:
    """JDFTx out file parsing class.

    A class to read and process a JDFTx out file.

    Methods
    -------
    from_file(file_path: str | Path) -> JDFTXOutfile
        Return JDFTXOutfile object from the path to a JDFTx out file.

    Attributes
    ----------
    slices: list[JDFTXOutfileSlice]
        A list of JDFTXOutfileSlice objects. Each slice corresponds to an individual
        call of the JDFTx executable. Subsequent JDFTx calls within the same directory
        and prefix will append outputs to the same out file. More than one slice
        may correspond to restarted calculations, geom + single point calculations,
        or optimizations done with 3rd-party wrappers like ASE.


    Properties
    ----------
    prefix: str
        The prefix of the most recent JDFTx call.

    jstrucs: JOutStructures
        The JOutStructures object from the most recent JDFTx call. This object contains
        a series of JOutStructure objects in its 'slices' attribute, each corresponding
        to a single structure (multiple iff performing a geometric optimization) as well
        as convergence data for the structures as a series.

    jsettings_fluid: JMinSettingsFluid (JminSettings)
        The JMinSettingsFluid object from the most recent JDFTx call. This object
        contains only a 'params' attribute, which is a dictionary of the input parameters
        for the fluid optimization.

    jsettings_electronic: JMinSettingsElectronic (JminSettings)
        The JMinSettingsElectronic object from the most recent JDFTx call. This object
        contains only a 'params' attribute, which is a dictionary of the input parameters
        for the electronic optimization.

    jsettings_lattice: JMinSettingsLattice (JminSettings)
        The JMinSettingsLattice object from the most recent JDFTx call. This object
        contains only a 'params' attribute, which is a dictionary of the input parameters
        for the lattice optimization.

    jsettings_ionic: JMinSettingsIonic (JminSettings)
        The JMinSettingsIonic object from the most recent JDFTx call. This object
        contains only a 'params' attribute, which is a dictionary of the input parameters
        for the ionic optimization.

    xc_func: str
        The exchange-correlation functional used in the most recent JDFTx call.
        See documentation for JDFTx online for a list of available exchange-correlation
        functionals.

    lattice_initial: np.ndarray
        The initial lattice vectors of the most recent JDFTx call as a 3x3 numpy array.
        In units of Angstroms.

    lattice_final: np.ndarray
        The final lattice vectors of the most recent JDFTx call as a 3x3 numpy array.
        In units of Angstroms.

    lattice: np.ndarray
        The lattice vectors of the most recent JDFTx call as a 3x3 numpy array
        (redundant to lattice_final).

    a: float
        Length of the first lattice vector. In units of Angstroms.

    b: float
        Length of the second lattice vector. In units of Angstroms.

    c: float
        Length of the third lattice vector. In units of Angstroms.

    fftgrid: list[int]
        The FFT grid shape used in the most recent JDFTx call. Can be used to properly
        shape densities dumped as binary files.

    geom_opt: bool
        True if the most recent JDFTx call was a geometry optimization (lattice or ionic).

    geom_opt_type: str
        The type of geometry optimization performed in the most recent JDFTx call.
        Options are 'lattice' or 'ionic' if geom_opt, else "single point".
        ('lattice' optimizations perform ionic optimizations as well unless ion
        positions are given in direct coordinates).

    efermi: float
        The Fermi energy in eV of the most recent JDFTx call. Equivalent to "mu".

    egap: float
        The band gap in eV of the most recent JDFTx call. (Only available if eigstats was dumped).

    emin: float
        The minimum energy in eV (smallest Kohn-Sham eigenvalue) of the most recent JDFTx call.
        (Only available if eigstats was dumped).

    emax: float
        The maximum energy in eV (largest Kohn-Sham eigenvalue) of the most recent JDFTx call.
        (Only available if eigstats was dumped).

    homo: float
        The energy in eV of the band-gap lower bound (Highest Occupied Molecular Orbital)
        (Only available if eigstats was dumped).

    lumo: float
        The energy in eV of the band-gap upper bound (Lowet Unoccupied Molecular Orbital)
        (Only available if eigstats was dumped).

    homo_filling: float
        The electron filling at the homo band-state.
        (Only available if eigstats was dumped).

    lumo_filling: float
        The electron filling at the lumo band-state.
        (Only available if eigstats was dumped).

    is_metal: bool
        True if fillings of homo and lumo band-states are off-set by 1 and 0
        by at least an arbitrary tolerance of 0.01 (ie 1 - 0.015 and 0.012 for
        homo/lumo fillings would be metallic, while 1-0.001 and 0 would not be).
        (Only available if eigstats was dumped).

    etype: str
        String representation of total energy-type of system. Commonly "G"
        (grand-canonical potential) for GC calculations, and "F" for canonical
        (fixed electron count) calculations.

    broadening_type: str
        Type of broadening for electronic filling about Fermi-level requested. Either
        "Fermi", "Cold", "MP1", or "Gauss".

    broadening: float
        Magnitude of broadening for electronic filling.

    kgrid: list[int]
        Shape of k-point grid used in calculation. (equivalent to k-point folding)

    truncation_type: str
        Type of coulomb truncation used to prevent interaction between periodic images
        along certain directions. "periodic" means no coulomb truncation was used.

    truncation_radius: float | None
        If spherical truncation_type, this is the radius of the coulomb truncation sphere.

    pwcut: float
        The plane-wave cutoff energy in Hartrees used in the most recent JDFTx call.

    rhocut: float
        The density cutoff energy in Hartrees used in the most recent JDFTx call.

    pp_type: str
        The pseudopotential library used in the most recent JDFTx call.
        Currently only "GBRV" and "SG15" are supported by this output parser.

    total_electrons: float
        The total number of electrons in the most recent JDFTx call (redundant
        to nelectrons).

    semicore_electrons: int
        The number of semicore electrons in the most recent JDFTx call.

    valence_electrons: float
        The number of valence electrons in the most recent JDFTx call.

    total_electrons_uncharged: int
        The total number of electrons in the most recent JDFTx call, uncorrected for
        charge. (ie total_electrons + charge)

    semicore_electrons_uncharged: int
        The number of semicore electrons in the most recent JDFTx call, uncorrected for
        charge. (ie semicore_electrons + charge)

    valence_electrons_uncharged: int
        The number of valence electrons in the most recent JDFTx call, uncorrected for
        charge. (ie valence_electrons + charge)

    nbands: int
        The number of bands used in the most recent JDFTx call.

    atom_elements: list[str]
        The list of each ion's element symbol in the most recent JDFTx call.

    atom_elements_int: list[int]
        The list of ion's atomic numbers in the most recent JDFTx call.

    atom_types: list[str]
        Non-repeating list of each ion's element symbol in the most recent JDFTx call.

    spintype: str
        The spin type used in the most recent JDFTx call. Options are "none", "collinear",

    nspin: int
        The number of spins used in the most recent JDFTx call.

    nat: int
        The number of atoms in the most recent JDFTx call.

    atom_coords_initial: list[list[float]]
        The initial atomic coordinates of the most recent JDFTx call.

    atom_coords_final: list[list[float]]
        The final atomic coordinates of the most recent JDFTx call.

    atom_coords: list[list[float]]
        The atomic coordinates of the most recent JDFTx call.

    has_solvation: bool
        True if the most recent JDFTx call included a solvation calculation.

    fluid: str
        The fluid used in the most recent JDFTx call.

    is_gc: bool
        True if the most recent slice is a grand canonical calculation.

    eopt_type: str
        The type of energy iteration used in the most recent JDFTx call.

    elecmindata: JElSteps
        The JElSteps object from the most recent JDFTx call. This object contains
        a series of JElStep objects in its 'steps' attribute, each corresponding to a
        single energy iteration.

    stress: np.ndarray
        The stress tensor of the most recent JDFTx call as a 3x3 numpy array.
        In units of eV/Angstrom^3.

    strain: np.ndarray
        The strain tensor of the most recent JDFTx call as a 3x3 numpy array.

    nstep: int
        The number of geometric optiization steps in the most recent JDFTx call.

    e: float
        The final energy in eV of the most recent JDFTx call (equivalent to the call's
        etype).

    grad_k: float
        The final norm of the preconditioned gradient for geometric optimization
        of the most recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient
        and Kg is the preconditioned gradient).
        (written as "|grad|_K" in JDFTx output).

    alpha: float
        The step size of the final geometric step in the most recent JDFTx call.

    linmin: float
        The final normalized projection of the geometric step direction onto the
        gradient for the most recent JDFTx call.

    abs_magneticmoment: float | None
        The absolute magnetic moment of the most recent JDFTx call.

    tot_magneticmoment: float | None
        The total magnetic moment of the most recent JDFTx call.

    mu: float
        The Fermi energy in eV of the most recent JDFTx call.

    elec_e: float
        The final energy in eV of the most recent electronic optimization step.

    elec_nstep: int
        The number of electronic optimization steps in the most recent JDFTx call.

    elec_grad_k: float
        The final norm of the preconditioned gradient for electronic optimization
        of the most recent JDFTx call (evaluated as dot(g, Kg), where g is the gradient
        and Kg is the preconditioned gradient).
        (written as "|grad|_K" in JDFTx output).

    elec_alpha: float
        The step size of the final electronic step in the most recent JDFTx call.

    elec_linmin: float
        The final normalized projection of the electronic step direction onto the
        gradient for the most recent JDFTx call.

    Magic Methods
    -------------
    __getitem__(key: str | int) -> Any
        Decides behavior of how JDFTXOutfile objects are indexed. If the key is a string,
        it will return the value of the property with the same name. If the key is an
        integer, it will return the slice of the JDFTXOutfile object at that index.

    __len__() -> int
        Returns the number of slices in the JDFTXOutfile object.

    __getattr__(name: str) -> Any
        Returns the value of the property with the same name as the input string.

    __dir__() -> list[str]
        Returns a list of all the properties of the JDFTXOutfile object.

    __str__() -> str
        Returns a string representation of the JDFTXOutfile object.
    """

    slices: list[JDFTXOutfileSlice] = field(default_factory=list)

    @classmethod
    def from_file(cls, file_path: str | Path, is_bgw: bool = False) -> JDFTXOutfile:
        """Return JDFTXOutfile object.

        Create a JDFTXOutfile object from a JDFTx out file.

        Parameters
        ----------
        file_path: str | Path
            The path to the JDFTx out file

        is_bgw: bool
            Mark True if data must be usable for BGW calculations. This will change
            the behavior of the parser to be stricter with certain criteria.

        Returns
        -------
        instance: JDFTXOutfile
            The JDFTXOutfile object
        """
        texts = read_outfile_slices(file_path)
        slices = [JDFTXOutfileSlice.from_out_slice(text) for text in texts]
        return cls(slices=slices)

    ###########################################################################
    # Properties inherited from most recent JDFTXOutfileSlice
    ###########################################################################

    @property
    def prefix(self) -> str:
        """
        The prefix of the most recent JDFTx call.

        Return prefix from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].prefix
        raise AttributeError("Property prefix inaccessible due to empty slices class field")

    @property
    def jstrucs(self) -> JOutStructures:
        """

        Return jstrucs from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jstrucs
        raise AttributeError("Property jstrucs inaccessible due to empty slices class field")

    @property
    def jsettings_fluid(
        self,
    ) -> JMinSettingsFluid | JMinSettingsElectronic | JMinSettingsLattice | JMinSettingsIonic:
        """
        Return jsettings_fluid from most recent JOutStructure.

        Return jsettings_fluid from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_fluid
        raise AttributeError("Property jsettings_fluid inaccessible due to empty slices class field")

    @property
    def jsettings_electronic(
        self,
    ) -> JMinSettingsFluid | JMinSettingsElectronic | JMinSettingsLattice | JMinSettingsIonic:
        """
        Return jsettings_electronic from most recent JOutStructure.

        Return jsettings_electronic from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_electronic
        raise AttributeError("Property jsettings_electronic inaccessible due to empty slices class field")

    @property
    def jsettings_lattice(
        self,
    ) -> JMinSettingsFluid | JMinSettingsElectronic | JMinSettingsLattice | JMinSettingsIonic:
        """
        Return jsettings_lattice from most recent JOutStructure.

        Return jsettings_lattice from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_lattice
        raise AttributeError("Property jsettings_lattice inaccessible due to empty slices class field")

    @property
    def jsettings_ionic(
        self,
    ) -> JMinSettingsFluid | JMinSettingsElectronic | JMinSettingsLattice | JMinSettingsIonic:
        """
        Return jsettings_ionic from most recent JOutStructure.

        Return jsettings_ionic from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_ionic
        raise AttributeError("Property jsettings_ionic inaccessible due to empty slices class field")

    @property
    def xc_func(self) -> str:
        """
        Return xc_func from most recent JOutStructure.

        Return xc_func from most recent JOutStructure, where xc_func is the name of the exchange
        correlation functional used for the calculation.
        """
        if len(self.slices):
            return self.slices[-1].xc_func
        raise AttributeError("Property xc_func inaccessible due to empty slices class field")

    @property
    def lattice_initial(self) -> np.ndarray:
        """
        Return lattice_initial from most recent JOutStructure.

        Return lattice_initial from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice_initial
        raise AttributeError("Property lattice_initial inaccessible due to empty slices class field")

    @property
    def lattice_final(self) -> np.ndarray:
        """
        Return lattice_final from most recent JOutStructure.

        Return lattice_final from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice_final
        raise AttributeError("Property lattice_final inaccessible due to empty slices class field")

    @property
    def lattice(self) -> np.ndarray:
        """
        Return lattice from most recent JOutStructure.

        Return lattice from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice
        raise AttributeError("Property lattice inaccessible due to empty slices class field")

    @property
    def a(self) -> float:
        """
        Return a from most recent JOutStructure.

        Return a from most recent JOutStructure, where a is the length of first lattice vector
        in Angstroms.
        """
        if len(self.slices):
            return self.slices[-1].a
        raise AttributeError("Property a inaccessible due to empty slices class field")

    @property
    def b(self) -> float:
        """
        Return b from most recent JOutStructure.

        Return b from most recent JOutStructure, where b is the length of the second lattice vector
        in Angstroms.
        """
        if len(self.slices):
            return self.slices[-1].b
        raise AttributeError("Property b inaccessible due to empty slices class field")

    @property
    def c(self) -> float:
        """
        Return c from most recent JOutStructure.

        Return c from most recent JOutStructure, where c is the length of the second lattice vector
        in Angstroms.
        """
        if len(self.slices):
            return self.slices[-1].c
        raise AttributeError("Property c inaccessible due to empty slices class field")

    @property
    def fftgrid(self) -> list[int]:
        """
        Return fftgrid from most recent JOutStructure.

        Return fftgrid from most recent JOutStructure, where fftgrid is the shape of electronic density
        array (as needed for reshaping binary files of real-space arrays).
        """
        if len(self.slices):
            return self.slices[-1].fftgrid
        raise AttributeError("Property fftgrid inaccessible due to empty slices class field")

    @property
    def geom_opt(self) -> bool:
        """
        Return geom_opt from most recent JOutStructure.

        Return geom_opt from most recent JOutStructure, where geom_opt is True if calculation
        included some type of geometric optimization, and False if not (ie single point calculation).
        """
        if len(self.slices):
            return self.slices[-1].geom_opt
        raise AttributeError("Property geom_opt inaccessible due to empty slices class field")

    @property
    def geom_opt_type(self) -> str:
        """
        Return geom_opt_type from most recent JOutStructure.

        Return geom_opt_type from most recent JOutStructure, where geom_opt_type is the type of geometric
        optimization performed (lattice, ionic, or single point).
        """
        if len(self.slices):
            return self.slices[-1].geom_opt_type
        raise AttributeError("Property geom_opt_type inaccessible due to empty slices class field")

    @property
    def efermi(self) -> float | None:
        """
        Return efermi from most recent JOutStructure.

        Return efermi from most recent JOutStructure, where efermi is the energy of the Fermi level
        in eV.
        """
        if len(self.slices):
            return self.slices[-1].efermi
        raise AttributeError("Property efermi inaccessible due to empty slices class field")

    @property
    def egap(self) -> float | None:
        """
        Return egap from most recent JOutStructure.

        Return egap from most recent JOutStructure, where egap is the size of the band gap in eV.
        """
        if len(self.slices):
            return self.slices[-1].egap
        raise AttributeError("Property egap inaccessible due to empty slices class field")

    @property
    def emin(self) -> float | None:
        """
        Return emin from most recent JOutStructure.

        Return emin from most recent JOutStructure, where emin is the lowest Kohn-Sham eigenvalue
        in eV.
        """
        if len(self.slices):
            return self.slices[-1].emin
        raise AttributeError("Property emin inaccessible due to empty slices class field")

    @property
    def emax(self) -> float | None:
        """
        Return emax from most recent JOutStructure.

        Return emax from most recent JOutStructure, where emax is the highest Kohn-Sham eigenvalue
        in eV.
        """
        if len(self.slices):
            return self.slices[-1].emax
        raise AttributeError("Property emax inaccessible due to empty slices class field")

    @property
    def homo(self) -> float | None:
        """
        Return homo from most recent JOutStructure.

        Return homo from most recent JOutStructure, where homo is the energy of last band-state before Fermi level
        (acronym for Highest Occupied Molecular Orbital, even though these are not molecular orbitals and this
        state may not be entirely occupied).
        (None if eigstats are not dumped)
        """
        if len(self.slices):
            return self.slices[-1].homo
        raise AttributeError("Property homo inaccessible due to empty slices class field")

    @property
    def lumo(self) -> float | None:
        """
        Return lumo from most recent JOutStructure.

        Return lumo from most recent JOutStructure, where lumo is the energy of first band-state
        after Fermi level (acronym for Lowest Unoccupied Molecular Orbital, even though these
        are not molecular orbitals, and this state may not be entirely unoccupied)
        (None if eigstats are not dumped)
        """
        if len(self.slices):
            return self.slices[-1].lumo
        raise AttributeError("Property lumo inaccessible due to empty slices class field")

    @property
    def homo_filling(self) -> float | None:
        """
        Return homo_filling from most recent JOutStructure.

        Return homo_filling from most recent JOutStructure, where homo_filling is the filling at the "homo"
        energy level.
        """
        if len(self.slices):
            return self.slices[-1].homo_filling
        raise AttributeError("Property homo_filling inaccessible due to empty slices class field")

    @property
    def lumo_filling(self) -> float | None:
        """
        Return lumo_filling from most recent JOutStructure.

        Return lumo_filling from most recent JOutStructure, where lumo_filling is the filling at the "lumo"
        energy level.
        """
        if len(self.slices):
            return self.slices[-1].lumo_filling
        raise AttributeError("Property lumo_filling inaccessible due to empty slices class field")

    @property
    def is_metal(self) -> bool | None:
        """
        Return is_metal from most recent JOutStructure.

        Return is_metal from most recent JOutStructure, where is_metal is true if fillings of
        homo and lumo band-states are off-set by 1 and 0 by at least an arbitrary tolerance of
        0.01 (ie 1 - 0.015 and 0.012 for homo/lumo fillings would be metallic, while 1-0.001
        and 0 would not be).
        (Only available if eigstats was dumped).
        """
        if len(self.slices):
            return self.slices[-1].is_metal
        raise AttributeError("Property is_metal inaccessible due to empty slices class field")

    @property
    def etype(self) -> str | None:
        """
        Return etype from most recent JOutStructure.

        Return etype from most recent JOutStructure, where etype is the string representation of the
        energy type by which the electronic ensemble was minimized (G, grand-canonical potential for
        grand-canonical ensemble; F, Helmholtz, for canonical ensemble).
        """
        if len(self.slices):
            return self.slices[-1].etype
        raise AttributeError("Property etype inaccessible due to empty slices class field")

    @property
    def broadening_type(self) -> str:
        """
        Return broadening_type from most recent JOutStructure.

        Return broadening_type from most recent JOutStructure, where broadening_type is the function used for
        smearing electronic filling about the Fermi level.
        """
        if len(self.slices):
            return self.slices[-1].broadening_type
        raise AttributeError("Property broadening_type inaccessible due to empty slices class field")

    @property
    def broadening(self) -> float:
        """
        Return broadening from most recent JOutStructure.

        Return broadening from most recent JOutStructure, where broadening is the parameter controlling
        the magnitude of broadening of electronic filling about the Fermi level.
        """
        if len(self.slices):
            return self.slices[-1].broadening
        raise AttributeError("Property broadening inaccessible due to empty slices class field")

    @property
    def kgrid(self) -> list:
        """
        Return kgrid from most recent JOutStructure.

        Return kgrid from most recent JOutStructure, where kgrid is the shape of the k-point mesh used to
        sample the Brillouin-zone of the unit cell (equivalent to kpoint folding).
        """
        if len(self.slices):
            return self.slices[-1].kgrid
        raise AttributeError("Property kgrid inaccessible due to empty slices class field")

    @property
    def truncation_type(self) -> str:
        """
        Return truncation_type from most recent JOutStructure.

        Return truncation_type from most recent JOutStructure, where truncation type is the type of Coloumb
        truncation used to avoid interaction with neighboring periodic images. ("Periodic" if no truncation)
        """
        if len(self.slices):
            return self.slices[-1].truncation_type
        raise AttributeError("Property truncation_type inaccessible due to empty slices class field")

    @property
    def truncation_radius(self) -> float | None:
        """
        Return truncation_radius from most recent JOutStructure.

        Return truncation_radius from most recent JOutStructure, where truncation_radius is the radius of coloumb
        truncation boundary in Bohr (not None iff truncation_type is spherical).
        """
        if len(self.slices):
            return self.slices[-1].truncation_radius
        raise AttributeError("Property truncation_radius inaccessible due to empty slices class field")

    @property
    def pwcut(self) -> float:
        """
        Return pwcut from most recent JOutStructure.

        Return pwcut from most recent JOutStructure, where pwcut is the energy cutoff for planewaves entering
        the basis set in Hartree.
        """
        if len(self.slices):
            return self.slices[-1].pwcut
        raise AttributeError("Property pwcut inaccessible due to empty slices class field")

    @property
    def rhocut(self) -> float:
        """
        Return rhocut from most recent JOutStructure.

        Return rhocut from most recent JOutStructure, where rhocut is the energy cutoff for the resolution
        of the real-space grid in Hartree.
        """
        if len(self.slices):
            return self.slices[-1].rhocut
        raise AttributeError("Property rhocut inaccessible due to empty slices class field")

    @property
    def pp_type(self) -> str | None:
        """
        Return pp_type from most recent JOutStructure.

        Return pp_type from most recent JOutStructure, where pp_type is the name of the pseudopotential
        library used for the calculation. Only "GBRV" and "SG15" are supported by this output parser,
        otherwise pp_type is None.
        """
        if len(self.slices):
            return self.slices[-1].pp_type
        raise AttributeError("Property pp_type inaccessible due to empty slices class field")

    @property
    def total_electrons(self) -> float:
        """
        Return total_electrons from most recent JOutStructure.

        Return total_electrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].total_electrons
        raise AttributeError("Property total_electrons inaccessible due to empty slices class field")

    @property
    def semicore_electrons(self) -> int:
        """
        Return semicore_electrons from most recent JOutStructure.

        Return semicore_electrons from most recent JOutStructure, where semicore_electrons are electrons
        discluded from pseudopotentials but not part of the atom's valence shell.
        """
        if len(self.slices):
            return self.slices[-1].semicore_electrons
        raise AttributeError("Property semicore_electrons inaccessible due to empty slices class field")

    @property
    def valence_electrons(self) -> float:
        """
        Return valence_electrons from most recent JOutStructure.

        Return valence_electrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].valence_electrons
        raise AttributeError("Property valence_electrons inaccessible due to empty slices class field")

    @property
    def total_electrons_uncharged(self) -> int:
        """
        Return total_electrons_uncharged from most recent JOutStructure.

        Return total_electrons_uncharged from most recent JOutStructure, where total_electrons_uncharged indicates
        the number of electrons required to reach a neutral cell charge.
        """
        if len(self.slices):
            return self.slices[-1].total_electrons_uncharged
        raise AttributeError("Property total_electrons_uncharged inaccessible due to empty slices class field")

    @property
    def semicore_electrons_uncharged(self) -> int:
        """
        Return semicore_electrons_uncharged from most recent JOutStructure.

        Return semicore_electrons_uncharged from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].semicore_electrons_uncharged
        raise AttributeError("Property semicore_electrons_uncharged inaccessible due to empty slices class field")

    @property
    def valence_electrons_uncharged(self) -> int:
        """
        Return valence_electrons_uncharged from most recent JOutStructure.

        Return valence_electrons_uncharged from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].valence_electrons_uncharged
        raise AttributeError("Property valence_electrons_uncharged inaccessible due to empty slices class field")

    @property
    def nbands(self) -> int:
        """
        Return nbands from most recent JOutStructure.

        Return nbands from most recent JOutStructure, where nbands is the number of bands used in the calculation.
        """
        if len(self.slices):
            return self.slices[-1].nbands
        raise AttributeError("Property nbands inaccessible due to empty slices class field")

    @property
    def atom_elements(self) -> list[str]:
        """
        Return atom_elements from most recent JOutStructure.

        Return atom_elements from most recent JOutStructure, where atom_elements is the list of each ion's element
        symbol in the most recent JDFTx call.
        """
        if len(self.slices):
            return self.slices[-1].atom_elements
        raise AttributeError("Property atom_elements inaccessible due to empty slices class field")

    @property
    def atom_elements_int(self) -> list[int]:
        """
        Return atom_elements_int from most recent JOutStructure.

        Return atom_elements_int from most recent JOutStructure, where atom_elements_int is the list of ion's atomic
        numbers in the most recent JDFTx call.
        """
        if len(self.slices):
            return self.slices[-1].atom_elements_int
        raise AttributeError("Property atom_elements_int inaccessible due to empty slices class field")

    @property
    def atom_types(self) -> list:
        """
        Return atom_types from most recent JOutStructure.

        Return atom_types from most recent JOutStructure, where atom_types is the non-repeating list of each ion's
        element symbol in the most recent JDFTx call.
        """
        if len(self.slices):
            return self.slices[-1].atom_types
        raise AttributeError("Property atom_types inaccessible due to empty slices class field")

    @property
    def spintype(self) -> str:
        """
        Return spintype from most recent JOutStructure.

        Return spintype from most recent JOutStructure, where spintype indicates the way spin was incorporated.
        """
        if len(self.slices):
            return self.slices[-1].spintype
        raise AttributeError("Property spintype inaccessible due to empty slices class field")

    @property
    def nspin(self) -> int:
        """
        Return nspin from most recent JOutStructure.

        Return nspin from most recent JOutStructure, where nspin indicates the number of spins used in the calculation.
        """
        if len(self.slices):
            return self.slices[-1].nspin
        raise AttributeError("Property nspin inaccessible due to empty slices class field")

    @property
    def nat(self) -> int:
        """
        Return nat from most recent JOutStructure.

        Return nat from most recent JOutStructure, where nat is the number of atoms in the most recent JDFTx call.
        """
        if len(self.slices):
            return self.slices[-1].nat
        raise AttributeError("Property nat inaccessible due to empty slices class field")

    @property
    def atom_coords_initial(self) -> list[list[float]]:
        """
        Return atom_coords_initial from most recent JOutStructure.

        Return atom_coords_initial from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords_initial
        raise AttributeError("Property atom_coords_initial inaccessible due to empty slices class field")

    @property
    def atom_coords_final(self) -> list[list[float]]:
        """
        Return atom_coords_final from most recent JOutStructure.

        Return atom_coords_final from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords_final
        raise AttributeError("Property atom_coords_final inaccessible due to empty slices class field")

    @property
    def atom_coords(self) -> list[list[float]]:
        """
        Return atom_coords from most recent JOutStructure.

        Return atom_coords from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords
        raise AttributeError("Property atom_coords inaccessible due to empty slices class field")

    @property
    def has_solvation(self) -> bool:
        """
        Return has_solvation from most recent JOutStructure.

        Return has_solvation from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].has_solvation
        raise AttributeError("Property has_solvation inaccessible due to empty slices class field")

    @property
    def fluid(self) -> str:
        """
        Return fluid from most recent JOutStructure.

        Return fluid from most recent JOutStructure, where fluid is the name of the implicit solvent used in the
        calculation.
        """
        if len(self.slices):
            return self.slices[-1].fluid
        raise AttributeError("Property fluid inaccessible due to empty slices class field")

    @property
    def is_gc(self) -> bool:
        """
        Return is_gc from most recent JOutStructure.

        Return is_gc from most recent JOutStructure, where is_gc is True if the most recent slice is a grand
        canonical calculation.

        Returns
        -------
        bool
            True if the most recent slice is a grand canonical calculation
        """
        if len(self.slices):
            return self.slices[-1].is_gc
        raise AttributeError("Property is_gc inaccessible due to empty slices class field")

    ###########################################################################
    # Properties inherited from most recent JDFTXOutfileSlice directly through
    # the JDFTXOutfileSlice object's jstrucs class variable.
    ###########################################################################

    @property
    def eopt_type(self) -> str:
        """
        Return eopt_type from most recent JOutStructure.

        Return eopt_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].eopt_type
        raise AttributeError("Property eopt_type inaccessible due to empty jstrucs class field")

    @property
    def elecmindata(self) -> JElSteps:
        """
        Return elecmindata from most recent JOutStructure.

        Return elecmindata from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].elecmindata
        raise AttributeError("Property elecmindata inaccessible due to empty jstrucs class field")

    @property
    def stress(self) -> np.ndarray:
        """
        Return stress from most recent JOutStructure.

        Return stress from most recent JOutStructure, where stress is the stress tensor of the unit cell
        in units eV/A^3.
        """
        if len(self.slices):
            return self.slices[-1].stress
        raise AttributeError("Property stress inaccessible due to empty jstrucs class field")

    @property
    def strain(self) -> np.ndarray:
        """
        Return strain from most recent JOutStructure.

        Return strain from most recent JOutStructure, where strain is the unitless strain tensor.
        """
        if len(self.slices):
            return self.slices[-1].strain
        raise AttributeError("Property strain inaccessible due to empty jstrucs class field")

    @property
    def nstep(self) -> int:
        """
        Return (geometric) step number from most recent JOutStructure.

        Return (geometric) step number from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nstep
        raise AttributeError("Property nstep inaccessible due to empty jstrucs class field")

    @property
    def e(self) -> float:
        """
        Return e from most recent JOutStructure.

        Return e from most recent JOutStructure, where e is the energy of the system's etype in eV.
        """
        if len(self.slices):
            return self.slices[-1].e
        raise AttributeError("Property e inaccessible due to empty jstrucs class field")

    @property
    def grad_k(self) -> float:
        """
        Return (geometric) grad_k from most recent JOutStructure.

        Return (geometric) grad_k from most recent JOutStructure, where grad_k is the final
        norm of the preconditioned gradient for geometric optimization of the most recent JDFTx
        call (evaluated as dot(g, Kg), where g is the gradient and Kg is the preconditioned gradient).
        (written as "|grad|_K" in JDFTx output)..
        """
        if len(self.slices):
            return self.slices[-1].grad_k
        raise AttributeError("Property grad_k inaccessible due to empty jstrucs class field")

    @property
    def alpha(self) -> float:
        """
        Return (geometric) alpha from most recent JOutStructure.

        Return (geometric) alpha from most recent JOutStructure, where alpha is the geometric step size
        along the line minimization.
        """
        if len(self.slices):
            return self.slices[-1].alpha
        raise AttributeError("Property alpha inaccessible due to empty jstrucs class field")

    @property
    def linmin(self) -> float:
        """
        Return (geometric) linmin from most recent JOutStructure.

        Return (geometric) linmin from most recent JOutStructure, where linmin is the final normalized
        projection of the geometric step direction onto the gradient for the most recent JDFTx call.
        """
        if len(self.slices):
            return self.slices[-1].linmin
        raise AttributeError("Property linmin inaccessible due to empty jstrucs class field")

    @property
    def nelectrons(self) -> float:
        """
        Return nelectrons from most recent JOutStructure.

        Return nelectrons from most recent JOutStructure, where nelectrons is the number of electron
        (equivalent to total_electrons).
        """
        if len(self.slices):
            return self.slices[-1].nelectrons
        raise AttributeError("Property nelectrons inaccessible due to empty jstrucs class field")

    @property
    def abs_magneticmoment(self) -> float | None:
        """
        Return abs_magneticmoment from most recent JOutStructure.

        Return abs_magneticmoment from most recent JOutStructure, where abs_magnetic_moment is the absolute
        magnetic moment of electronic density. (None if restricted spin)
        """
        if len(self.slices):
            return self.slices[-1].abs_magneticmoment
        raise AttributeError("Property abs_magneticmoment inaccessible due to empty jstrucs class field")

    @property
    def tot_magneticmoment(self) -> float | None:
        """
        Return tot_magneticmoment from most recent JOutStructure.

        Return tot_magneticmoment from most recent JOutStructure, where tot_magnetic_moment is the total
        magnetic moment of the electronic density. (None if restricted spin)
        """
        if len(self.slices):
            return self.slices[-1].tot_magneticmoment
        raise AttributeError("Property tot_magneticmoment inaccessible due to empty jstrucs class field")

    @property
    def mu(self) -> float:
        """
        Return mu from most recent JOutStructure.

        Return mu from most recent JOutStructure. (Equivalent to efermi)
        """
        if len(self.slices):
            return self.slices[-1].mu
        raise AttributeError("Property mu inaccessible due to empty jstrucs class field")

    ###########################################################################
    # Electronic properties with symbol disambiguation inherited from most
    # recent JDFTXOutfileSlice directly through the JDFTXOutfileSlice
    # object's jstrucs class variable.
    ###########################################################################

    @property
    def elec_nstep(self) -> int:
        """Return the most recent electronic step number.

        Return the most recent electronic step number.

        Returns
        -------
        elec_nstep: int
        """
        if len(self.slices):
            return self.slices[-1].elec_nstep
        raise AttributeError("Property elec_inter inaccessible due to empty jstrucs class field")

    @property
    def elec_e(self) -> float:
        """Return the most recent electronic energy.

        Return the most recent electronic energy.

        Returns
        -------
        elec_e: float
        """
        if len(self.slices):
            return self.slices[-1].elec_e
        raise AttributeError("Property elec_e inaccessible due to empty jstrucs class field")

    @property
    def elec_grad_k(self) -> int:
        """Return the most recent electronic grad_k.

        Return the most recent electronic grad_k. (Equivalent to grad_k but for electronic line minimization)

        Returns
        -------
        grad_k: float
        """
        if len(self.slices):
            return self.slices[-1].elec_grad_k
        raise AttributeError("Property elec_grad_k inaccessible due to empty jstrucs class field")

    @property
    def elec_alpha(self) -> float:
        """Return the most recent electronic alpha.

        Return the most recent electronic alpha.  (Equivalent to alpha but for electronic line minimization)

        Returns
        -------
        alpha: float
        """
        if len(self.slices):
            return self.slices[-1].elec_linmin
        raise AttributeError("Property elec_alpha inaccessible due to empty jstrucs class field")

    @property
    def elec_linmin(self) -> float:
        """Return the most recent electronic linmin.

        Return the most recent electronic linmin.  (Equivalent to linmin but for electronic line minimization)

        Returns
        -------
        linmin: float
        """
        if len(self.slices):
            return self.slices[-1].elec_linmin
        raise AttributeError("Property elec_linmin inaccessible due to empty jstrucs class field")

    ###########################################################################
    # Magic methods
    ###########################################################################

    def __getitem__(self, key: int | str) -> JDFTXOutfileSlice | Any:
        """Return item.

        Return the value of an item.

        Parameters
        ----------
        key: int | str
            The key of the item

        Returns
        -------
        val
            The value of the item
        """
        val = None
        if type(key) is int:
            val = self.slices[key]
        elif type(key) is str:
            val = getattr(self, key)
        else:
            raise TypeError(f"Invalid key type: {type(key)}")
        return val

    def __len__(self) -> int:
        """Return length of JDFTXOutfile object.

        Returns the number of JDFTx calls in the
        JDFTXOutfile object.

        Returns
        -------
        length: int
            The number of geometric optimization steps in the JDFTXOutfile
            object
        """
        return len(self.slices)

    # This method is likely never going to be called as all (currently existing)
    # attributes of the most recent slice are explicitly defined as a class
    # property. However, it is included to reduce the likelihood of errors
    # upon future changes to downstream code.
    def __getattr__(self, name: str) -> Any:
        """Return attribute.

        Return the value of an attribute.

        Parameters
        ----------
        name: str
            The name of the attribute

        Returns
        -------
        val
            The value of the attribute
        """
        if name in self.__dict__:
            return self.__dict__[name]

        # Check if the attribute is a property of the class
        for cls in inspect.getmro(self.__class__):
            if name in cls.__dict__ and isinstance(cls.__dict__[name], property):
                return cls.__dict__[name].__get__(self)

        # Check if the attribute is in self.jstrucs
        if hasattr(self.slices[-1], name):
            return getattr(self.slices[-1], name)

        # If the attribute is not found in either, raise an AttributeError
        raise AttributeError(f"{self.__class__.__name__} not found: {name}")

    def __dir__(self) -> list:
        """List attributes.

        Returns a list of attributes for the object, including those from
        self.slices[-1].

        Returns
        -------
        list
            A list of attribute names
        """
        # Get the default attributes
        default_attrs = dir(self)
        # Get the attributes from self.slices[-1] if slices is not empty
        slice_attrs = dir(self.slices[-1]) if self.slices else []
        # Combine and return unique attributes
        return list(set(default_attrs + slice_attrs))

    def __str__(self) -> str:
        """Return string representation of JDFTXOutfile object.

        Return a string representation of the JDFTXOutfile object.

        Returns
        -------
        str
            The string representation of the JDFTXOutfile object
        """
        return pprint.pformat(self)
