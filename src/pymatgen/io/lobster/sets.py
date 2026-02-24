"""
This module defines the VaspInputSet abstract base class and a concrete implementation for the parameters developed
and tested by the core team of pymatgen, including the Materials Virtual Lab, Materials Project and the MIT high
throughput project. The basic concept behind an input set is to specify a scheme to generate a consistent set of VASP
inputs from a structure without further user intervention. This ensures comparability across runs.

Read the following carefully before implementing new input sets:

1. 99% of what needs to be done can be done by specifying user_incar_settings to override some of the defaults of
   various input sets. Unless there is an extremely good reason to add a new set, **do not** add one. e.g. if you want
   to turn the Hubbard U off, just set "LDAU": False as a user_incar_setting.
2. All derivative input sets should inherit appropriate configurations (e.g., from MPRelaxSet), and more often than
   not, VaspInputSet should be the superclass. Superclass delegation should be used where possible. In particular,
   you are not supposed to implement your own as_dict or from_dict for derivative sets unless you know what you are
   doing. Improper overriding the as_dict and from_dict protocols is the major cause of implementation headaches. If
   you need an example, look at how the MPStaticSet is initialized.

The above are recommendations. The following are **UNBREAKABLE** rules:

1. All input sets must take in a structure, list of structures or None as the first argument. If None, the input set
   should perform a stateless initialization and before any output can be written, a structure must be set.
2. user_incar_settings, user_kpoints_settings and user_<whatever>_settings are ABSOLUTE. Any new sets you implement
   must obey this. If a user wants to override your settings, you assume he knows what he is doing. Do not
   magically override user supplied settings. You can issue a warning if you think the user is wrong.
3. All input sets must save all supplied args and kwargs as instance variables. e.g. self.arg = arg and
   self.kwargs = kwargs in the __init__. This ensures the as_dict and from_dict work correctly.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from pymatgen.io.lobster import Lobsterin
from pymatgen.io.vasp.sets import BadInputSetWarning, MPRelaxSet, VaspInputSet

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Literal

    UserPotcarFunctional = (
        Literal[
            "PBE",
            "PBE_52",
            "PBE_54",
            "PBE_64",
            "LDA",
            "LDA_52",
            "LDA_54",
            "PW91",
            "LDA_US",
            "PW91_US",
        ]
        | None
    )

MODULE_DIR = os.path.dirname(__file__)


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

    # Latest POTCARs are preferred
    # Choose PBE_54 unless the user specifies a different potcar_functional
    user_potcar_functional: UserPotcarFunctional = "PBE_54"

    CONFIG = MPRelaxSet.CONFIG
    _valid_potcars: Sequence[str] | None = ("PBE_52", "PBE_54", "PBE_64")

    def __post_init__(self) -> None:
        super().__post_init__()
        warnings.warn(
            "Make sure that all parameters are okay! This is a brand new implementation.",
            stacklevel=2,
        )

        if self.user_potcar_functional in ["PBE_52", "PBE_64"]:
            warnings.warn(
                f"Using {self.user_potcar_functional} POTCARs with basis functions generated for PBE_54 POTCARs. "
                "Basis functions for elements with obsoleted, updated or newly added POTCARs in "
                f"{self.user_potcar_functional} will not be available and may cause errors or inaccuracies.",
                BadInputSetWarning,
                stacklevel=2,
            )
        if self.isym not in {-1, 0}:
            raise ValueError("Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0")
        if self.ismear not in {-5, 0}:
            raise ValueError("Lobster usually works with ismear=-5 or ismear=0")

        self._config_dict["POTCAR"]["W"] = "W_sv"

    @property
    def kpoints_updates(self) -> dict[str, int]:
        """Updates to the kpoints configuration for this calculation type."""
        # Test if this is okay
        return {"reciprocal_density": self.reciprocal_density or 310}

    @property
    def incar_updates(self) -> dict[str, Any]:
        """Updates to the INCAR config for this calculation type."""

        potcar_symbols = self.potcar_symbols

        # Predefined basis! Check if the basis is okay! (charge spilling and bandoverlaps!)
        if self.user_supplied_basis is None and self.address_basis_file is None:
            basis = Lobsterin.get_basis(structure=self.structure, potcar_symbols=potcar_symbols)  # type: ignore[arg-type]
        elif self.address_basis_file is not None:
            basis = Lobsterin.get_basis(
                structure=self.structure,  # type: ignore[arg-type]
                potcar_symbols=potcar_symbols,
                address_basis_file=self.address_basis_file,
            )
        elif self.user_supplied_basis is not None:
            # Test if all elements from structure are in user_supplied_basis
            for atom_type in self.structure.symbol_set:  # type: ignore[union-attr]
                if atom_type not in self.user_supplied_basis:
                    raise ValueError(f"There are no basis functions for the atom type {atom_type}")
            basis = [f"{key} {value}" for key, value in self.user_supplied_basis.items()]
        else:
            basis = None

        lobsterin = Lobsterin(settingsdict={"basisfunctions": basis})
        nbands = lobsterin._get_nbands(structure=self.structure)  # type: ignore[arg-type]

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
