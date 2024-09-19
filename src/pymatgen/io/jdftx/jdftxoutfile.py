"""JDFTx out file parsing class.

A class to read and process a JDFTx out file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import wraps
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

from monty.io import zopen

from atomate2.jdftx.io.jdftxoutfileslice import JDFTXOutfileSlice
from atomate2.jdftx.io.jdftxoutfileslice_helpers import get_start_lines

if TYPE_CHECKING:
    import numpy as np

    from atomate2.jdftx.io.jminsettings import (
        JMinSettingsElectronic,
        JMinSettingsFluid,
        JMinSettingsIonic,
        JMinSettingsLattice,
    )
    from atomate2.jdftx.io.joutstructures import JOutStructures


def check_file_exists(func: Callable) -> Any:
    """Check if file exists.

    Check if file exists (and continue normally) or raise an exception if
    it does not.
    """

    @wraps(func)
    def wrapper(filename: str) -> Any:
        filepath = Path(filename)
        if not filepath.is_file():
            raise OSError(f"'{filename}' file doesn't exist!")
        return func(filename)

    return wrapper


@check_file_exists
def read_file(file_name: str) -> list[str]:
    """
    Read file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read

    Returns
    -------
    text: list[str]
        list of strings from file
    """
    with zopen(file_name, "r") as f:
        text = f.readlines()
    f.close()
    return text


def read_outfile_slices(file_name: str) -> list[list[str]]:
    """
    Read slice of out file into a list of str.

    Parameters
    ----------
    filename: Path or str
        name of file to read
    out_slice_idx: int
        index of slice to read from file

    Returns
    -------
    texts: list[list[str]]
        list of out file slices (individual calls of JDFTx)
    """
    _text = read_file(file_name)
    start_lines = get_start_lines(_text, add_end=True)
    texts = []
    for i in range(len(start_lines) - 1):
        text = _text[start_lines[i] : start_lines[i + 1]]
        texts.append(text)
    return texts


@dataclass
class JDFTXOutfile:
    """JDFTx out file parsing class.

    A class to read and process a JDFTx out file.
    """

    slices: list[JDFTXOutfileSlice] = field(default_factory=list)

    @classmethod
    def from_file(cls, file_path: str) -> JDFTXOutfile:
        """Return JDFTXOutfile object.

        Create a JDFTXOutfile object from a JDFTx out file.

        Parameters
        ----------
        file_path: str
            The path to the JDFTx out file

        Returns
        -------
        instance: JDFTXOutfile
            The JDFTXOutfile object
        """
        texts = read_outfile_slices(file_path)
        slices = [JDFTXOutfileSlice.from_out_slice(text) for text in texts]
        instance = cls()
        instance.slices = slices
        return instance

    ###########################################################################
    # Properties inherited from most recent JDFTXOutfileSlice
    ###########################################################################
    @property
    def prefix(self) -> str:
        """
        Return prefix from most recent JOutStructure.

        Return prefix from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].prefix
        raise AttributeError(
            "Property prefix inaccessible due to \
                             empty slices class field"
        )

    @property
    def jstrucs(self) -> JOutStructures:
        """
        Return jstrucs from most recent JOutStructure.

        Return jstrucs from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jstrucs
        raise AttributeError(
            "Property jstrucs inaccessible due to \
                             empty slices class field"
        )

    @property
    def jsettings_fluid(
        self,
    ) -> (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ):
        """
        Return jsettings_fluid from most recent JOutStructure.

        Return jsettings_fluid from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_fluid
        raise AttributeError(
            "Property jsettings_fluid inaccessible due to \
                             empty slices class field"
        )

    @property
    def jsettings_electronic(
        self,
    ) -> (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ):
        """
        Return jsettings_electronic from most recent JOutStructure.

        Return jsettings_electronic from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_electronic
        raise AttributeError(
            "Property jsettings_electronic inaccessible \
                             due to empty slices class field"
        )

    @property
    def jsettings_lattice(
        self,
    ) -> (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ):
        """
        Return jsettings_lattice from most recent JOutStructure.

        Return jsettings_lattice from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_lattice
        raise AttributeError(
            "Property jsettings_lattice inaccessible \
                             due to empty slices class field"
        )

    @property
    def jsettings_ionic(
        self,
    ) -> (
        JMinSettingsFluid
        | JMinSettingsElectronic
        | JMinSettingsLattice
        | JMinSettingsIonic
    ):
        """
        Return jsettings_ionic from most recent JOutStructure.

        Return jsettings_ionic from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].jsettings_ionic
        raise AttributeError(
            "Property jsettings_ionic inaccessible due to \
                             empty slices class field"
        )

    @property
    def xc_func(self) -> str:
        """
        Return xc_func from most recent JOutStructure.

        Return xc_func from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].xc_func
        raise AttributeError(
            "Property xc_func inaccessible due to empty \
                             slices class field"
        )

    @property
    def lattice_initial(self) -> np.ndarray:
        """
        Return lattice_initial from most recent JOutStructure.

        Return lattice_initial from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice_initial
        raise AttributeError(
            "Property lattice_initial inaccessible due to \
                             empty slices class field"
        )

    @property
    def lattice_final(self) -> np.ndarray:
        """
        Return lattice_final from most recent JOutStructure.

        Return lattice_final from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice_final
        raise AttributeError(
            "Property lattice_final inaccessible due to \
                             empty slices class field"
        )

    @property
    def lattice(self) -> np.ndarray:
        """
        Return lattice from most recent JOutStructure.

        Return lattice from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lattice
        raise AttributeError(
            "Property lattice inaccessible due to \
                             empty slices class field"
        )

    @property
    def a(self) -> float:
        """
        Return a from most recent JOutStructure.

        Return a from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].a
        raise AttributeError(
            "Property a inaccessible due to empty slices \
                             class field"
        )

    @property
    def b(self) -> float:
        """
        Return b from most recent JOutStructure.

        Return b from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].b
        raise AttributeError(
            "Property b inaccessible due to empty \
                             slices class field"
        )

    @property
    def c(self) -> float:
        """
        Return c from most recent JOutStructure.

        Return c from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].c
        raise AttributeError(
            "Property c inaccessible due to empty \
                             slices class field"
        )

    @property
    def fftgrid(self) -> list[int]:
        """
        Return fftgrid from most recent JOutStructure.

        Return fftgrid from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].fftgrid
        raise AttributeError(
            "Property fftgrid inaccessible due to empty \
                             slices class field"
        )

    @property
    def geom_opt(self) -> bool:
        """
        Return geom_opt from most recent JOutStructure.

        Return geom_opt from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].geom_opt
        raise AttributeError(
            "Property geom_opt inaccessible due to \
                             empty slices class field"
        )

    @property
    def geom_opt_type(self) -> str:
        """
        Return geom_opt_type from most recent JOutStructure.

        Return geom_opt_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].geom_opt_type
        raise AttributeError(
            "Property geom_opt_type inaccessible due to \
                             empty slices class field"
        )

    @property
    def efermi(self) -> float:
        """
        Return efermi from most recent JOutStructure.

        Return efermi from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].efermi
        raise AttributeError(
            "Property efermi inaccessible due to \
                             empty slices class field"
        )

    @property
    def egap(self) -> float:
        """
        Return egap from most recent JOutStructure.

        Return egap from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].egap
        raise AttributeError(
            "Property egap inaccessible due to \
                             empty slices class field"
        )

    @property
    def emin(self) -> float:
        """
        Return emin from most recent JOutStructure.

        Return emin from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].emin
        raise AttributeError(
            "Property emin inaccessible due to \
                             empty slices class field"
        )

    @property
    def emax(self) -> float:
        """
        Return emax from most recent JOutStructure.

        Return emax from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].emax
        raise AttributeError(
            "Property emax inaccessible due to \
                             empty slices class field"
        )

    @property
    def homo(self) -> float:
        """
        Return homo from most recent JOutStructure.

        Return homo from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].homo
        raise AttributeError(
            "Property homo inaccessible due to \
                             empty slices class field"
        )

    @property
    def lumo(self) -> float:
        """
        Return lumo from most recent JOutStructure.

        Return lumo from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lumo
        raise AttributeError(
            "Property lumo inaccessible due to \
                             empty slices class field"
        )

    @property
    def homo_filling(self) -> float:
        """
        Return homo_filling from most recent JOutStructure.

        Return homo_filling from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].homo_filling
        raise AttributeError(
            "Property homo_filling inaccessible due to \
                             empty slices class field"
        )

    @property
    def lumo_filling(self) -> float:
        """
        Return lumo_filling from most recent JOutStructure.

        Return lumo_filling from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].lumo_filling
        raise AttributeError(
            "Property lumo_filling inaccessible due to \
                             empty slices class field"
        )

    @property
    def is_metal(self) -> bool:
        """
        Return is_metal from most recent JOutStructure.

        Return is_metal from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].is_metal
        raise AttributeError(
            "Property is_metal inaccessible due to \
                             empty slices class field"
        )

    @property
    def etype(self) -> str:
        """
        Return etype from most recent JOutStructure.

        Return etype from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].etype
        raise AttributeError(
            "Property etype inaccessible due to \
                             empty slices class field"
        )

    @property
    def broadening_type(self) -> str:
        """
        Return broadening_type from most recent JOutStructure.

        Return broadening_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].broadening_type
        raise AttributeError(
            "Property broadening_type inaccessible due \
                             to empty slices class field"
        )

    @property
    def broadening(self) -> float:
        """
        Return broadening from most recent JOutStructure.

        Return broadening from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].broadening
        raise AttributeError(
            "Property broadening inaccessible \
                             due to empty slices class field"
        )

    @property
    def kgrid(self) -> list:
        """
        Return kgrid from most recent JOutStructure.

        Return kgrid from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].kgrid
        raise AttributeError(
            "Property kgrid inaccessible \
                             due to empty slices class field"
        )

    @property
    def truncation_type(self) -> str:
        """
        Return truncation_type from most recent JOutStructure.

        Return truncation_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].truncation_type
        raise AttributeError(
            "Property truncation_type inaccessible \
                             due to empty slices class field"
        )

    @property
    def truncation_radius(self) -> float:
        """
        Return truncation_radius from most recent JOutStructure.

        Return truncation_radius from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].truncation_radius
        raise AttributeError(
            "Property truncation_radius inaccessible \
                             due to empty slices class field"
        )

    @property
    def pwcut(self) -> float:
        """
        Return pwcut from most recent JOutStructure.

        Return pwcut from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].pwcut
        raise AttributeError(
            "Property pwcut inaccessible due to \
                             empty slices class field"
        )

    @property
    def rhocut(self) -> float:
        """
        Return rhocut from most recent JOutStructure.

        Return rhocut from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].rhocut
        raise AttributeError(
            "Property rhocut inaccessible due to empty \
                             slices class field"
        )

    @property
    def pp_type(self) -> str:
        """
        Return pp_type from most recent JOutStructure.

        Return pp_type from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].pp_type
        raise AttributeError(
            "Property pp_type inaccessible due to empty \
                             slices class field"
        )

    @property
    def total_electrons(self) -> float:
        """
        Return total_electrons from most recent JOutStructure.

        Return total_electrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].total_electrons
        raise AttributeError(
            "Property total_electrons inaccessible due to \
                             empty slices class field"
        )

    @property
    def semicore_electrons(self) -> int:
        """
        Return semicore_electrons from most recent JOutStructure.

        Return semicore_electrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].semicore_electrons
        raise AttributeError(
            "Property semicore_electrons inaccessible due \
                             to empty slices class field"
        )

    @property
    def valence_electrons(self) -> float:
        """
        Return valence_electrons from most recent JOutStructure.

        Return valence_electrons from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].valence_electrons
        raise AttributeError(
            "Property valence_electrons inaccessible due \
                             to empty slices class field"
        )

    @property
    def total_electrons_uncharged(self) -> int:
        """
        Return total_electrons_uncharged from most recent JOutStructure.

        Return total_electrons_uncharged from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].total_electrons_uncharged
        raise AttributeError(
            "Property total_electrons_uncharged \
                             inaccessible due to empty slices class field"
        )

    @property
    def semicore_electrons_uncharged(self) -> int:
        """
        Return semicore_electrons_uncharged from most recent JOutStructure.

        Return semicore_electrons_uncharged from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].semicore_electrons_uncharged
        raise AttributeError(
            "Property semicore_electrons_uncharged \
                             inaccessible due to empty slices class field"
        )

    @property
    def valence_electrons_uncharged(self) -> int:
        """
        Return valence_electrons_uncharged from most recent JOutStructure.

        Return valence_electrons_uncharged from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].valence_electrons_uncharged
        raise AttributeError(
            "Property valence_electrons_uncharged \
                             inaccessible due to empty slices class field"
        )

    @property
    def nbands(self) -> int:
        """
        Return Nbands from most recent JOutStructure.

        Return Nbands from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nbands
        raise AttributeError(
            "Property nbands inaccessible due to empty \
                             slices class field"
        )

    @property
    def atom_elements(self) -> list:
        """
        Return atom_elements from most recent JOutStructure.

        Return atom_elements from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_elements
        raise AttributeError(
            "Property atom_elements inaccessible due to \
                             empty slices class field"
        )

    @property
    def atom_elements_int(self) -> list:
        """
        Return atom_elements_int from most recent JOutStructure.

        Return atom_elements_int from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_elements_int
        raise AttributeError(
            "Property atom_elements_int inaccessible due \
                             to empty slices class field"
        )

    @property
    def atom_types(self) -> list:
        """
        Return atom_types from most recent JOutStructure.

        Return atom_types from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_types
        raise AttributeError(
            "Property atom_types inaccessible due to \
                             empty slices class field"
        )

    @property
    def spintype(self) -> str:
        """
        Return spintype from most recent JOutStructure.

        Return spintype from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].spintype
        raise AttributeError(
            "Property spintype inaccessible due to \
                             empty slices class field"
        )

    @property
    def nspin(self) -> int:
        """
        Return nspin from most recent JOutStructure.

        Return nspin from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nspin
        raise AttributeError(
            "Property nspin inaccessible due to \
                             empty slices class field"
        )

    @property
    def nat(self) -> int:
        """
        Return nat from most recent JOutStructure.

        Return nat from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].nat
        raise AttributeError(
            "Property nat inaccessible due to \
                             empty slices class field"
        )

    @property
    def atom_coords_initial(self) -> list[list[float]]:
        """
        Return atom_coords_initial from most recent JOutStructure.

        Return atom_coords_initial from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords_initial
        raise AttributeError(
            "Property atom_coords_initial inaccessible \
                             due to empty slices class field"
        )

    @property
    def atom_coords_final(self) -> list[list[float]]:
        """
        Return atom_coords_final from most recent JOutStructure.

        Return atom_coords_final from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords_final
        raise AttributeError(
            "Property atom_coords_final inaccessible \
                             due to empty slices class field"
        )

    @property
    def atom_coords(self) -> list[list[float]]:
        """
        Return atom_coords from most recent JOutStructure.

        Return atom_coords from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].atom_coords
        raise AttributeError(
            "Property atom_coords inaccessible \
                             due to empty slices class field"
        )

    @property
    def has_solvation(self) -> bool:
        """
        Return has_solvation from most recent JOutStructure.

        Return has_solvation from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].has_solvation
        raise AttributeError(
            "Property has_solvation inaccessible \
                             due to empty slices class field"
        )

    @property
    def fluid(self) -> str:
        """
        Return fluid from most recent JOutStructure.

        Return fluid from most recent JOutStructure.
        """
        if len(self.slices):
            return self.slices[-1].fluid
        raise AttributeError(
            "Property fluid inaccessible \
                             due to empty slices class field"
        )

    @property
    def is_gc(self) -> bool:
        """
        Return is_gc from most recent JOutStructure.

        Return is_gc from most recent JOutStructure.

        Returns
        -------
        bool
            True if the most recent slice is a grand canonical calculation
        """
        if len(self.slices):
            return self.slices[-1].is_gc
        raise AttributeError(
            "Property is_gc inaccessible \
                             due to empty slices class field"
        )

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
        if name not in self.__dict__:
            if not hasattr(self.slices[-1], name):
                raise AttributeError(f"{self.__class__.__name__} not found: {name}")
            return getattr(self.slices[-1], name)
        return self.__dict__[name]

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
