# ruff: noqa: PGH003

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets.base import UserPotcarFunctional, VaspInputSet
from pymatgen.io.vasp.sets.mp import MPRelaxSet

if TYPE_CHECKING:
    from collections.abc import Sequence

    from pymatgen.io.vasp.inputs import Kpoints


@dataclass
class LobsterSet(VaspInputSet):
    """Input set to prepare VASP runs that can be digested by Lobster (See cohp.de).

    Args:
        structure (Structure): input structure.
        isym (int): ISYM entry for INCAR, only isym=-1 and isym=0 are allowed
        ismear (int): ISMEAR entry for INCAR, only ismear=-5 and ismear=0 are allowed
        reciprocal_density (int): density of k-mesh by reciprocal volume
        user_supplied_basis (dict): dict including basis functions for all elements in
            structure, e.g. {"Fe": "3d 3p 4s", "O": "2s 2p"}; if not supplied, a
            standard basis is used
        address_basis_file (str): address to a file similar to
            "BASIS_PBE_54_standard.yaml" in pymatgen.io.lobster.lobster_basis
        user_potcar_settings (dict): dict including potcar settings for all elements in
            structure, e.g. {"Fe": "Fe_pv", "O": "O"}; if not supplied, a standard basis
            is used.
        **kwargs: Other kwargs supported by VaspInputSet.
    """

    isym: int = 0
    ismear: int = -5
    reciprocal_density: int | None = None
    address_basis_file: str | None = None
    user_supplied_basis: dict | None = None

    # newest potcars are preferred
    # Choose PBE_54 unless the user specifies a different potcar_functional
    user_potcar_functional: UserPotcarFunctional = "PBE_54"

    CONFIG = MPRelaxSet.CONFIG
    _valid_potcars: Sequence[str] | None = ("PBE_52", "PBE_54")

    def __post_init__(self):
        super().__post_init__()
        warnings.warn("Make sure that all parameters are okay! This is a brand new implementation.")

        if self.isym not in (-1, 0):
            raise ValueError("Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0")
        if self.ismear not in (-5, 0):
            raise ValueError("Lobster usually works with ismear=-5 or ismear=0")

        self._config_dict["POTCAR"]["W"] = "W_sv"

    @property
    def kpoints_updates(self) -> dict | Kpoints:
        """Updates to the kpoints configuration for this calculation type."""
        # test, if this is okay
        return {"reciprocal_density": self.reciprocal_density or 310}

    @property
    def incar_updates(self) -> dict:
        """Updates to the INCAR config for this calculation type."""
        from pymatgen.io.lobster import Lobsterin

        potcar_symbols = self.potcar_symbols

        # predefined basis! Check if the basis is okay! (charge spilling and bandoverlaps!)
        if self.user_supplied_basis is None and self.address_basis_file is None:
            basis = Lobsterin.get_basis(structure=self.structure, potcar_symbols=potcar_symbols)  # type: ignore
        elif self.address_basis_file is not None:
            basis = Lobsterin.get_basis(
                structure=self.structure,  # type: ignore
                potcar_symbols=potcar_symbols,
                address_basis_file=self.address_basis_file,
            )
        elif self.user_supplied_basis is not None:
            # test if all elements from structure are in user_supplied_basis
            for atom_type in self.structure.symbol_set:  # type: ignore
                if atom_type not in self.user_supplied_basis:
                    raise ValueError(f"There are no basis functions for the atom type {atom_type}")
            basis = [f"{key} {value}" for key, value in self.user_supplied_basis.items()]
        else:
            basis = None

        lobsterin = Lobsterin(settingsdict={"basisfunctions": basis})
        nbands = lobsterin._get_nbands(structure=self.structure)  # type: ignore

        return {
            "EDIFF": 1e-6,
            "NSW": 0,
            "LWAVE": True,
            "ISYM": self.isym,
            "NBANDS": nbands,
            "IBRION": -1,
            "ISMEAR": self.ismear,
            "LORBIT": 11,
            "ICHARG": 0,
            "ALGO": "Normal",
        }
