from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from pymatgen.io.vasp.sets.base import VaspInputSet, _load_yaml_config

if TYPE_CHECKING:
    from typing import Literal


@dataclass
class MatPESStaticSet(VaspInputSet):
    """Create input files for a MatPES static calculation.

    The goal of MatPES is to generate potential energy surface data. This is a distinctly different
    from the objectives of the MP static calculations, which aims to obtain primarily accurate
    energies and also electronic structure (DOS). For PES data, force accuracy (and to some extent,
    stress accuracy) is of paramount importance.

    The default POTCAR versions have been updated to PBE_54 from the old PBE set used in the
    MPStaticSet. However, **U values** are still based on PBE. The implicit assumption here is that
    the PBE_54 and PBE POTCARs are sufficiently similar that the U values fitted to the old PBE
    functional still applies.

    Args:
        structure (Structure): The Structure to create inputs for. If None, the input
            set is initialized without a Structure but one must be set separately before
            the inputs are generated.
        xc_functional ('R2SCAN'|'PBE'): Exchange-correlation functional to use. Defaults to 'PBE'.
        **kwargs: Keywords supported by VaspInputSet.
    """

    xc_functional: Literal["R2SCAN", "PBE", "PBE+U"] = "PBE"
    prev_incar: dict | str | None = None
    # These are parameters that we will inherit from any previous INCAR supplied. They are mostly parameters related
    # to symmetry and convergence set by Custodian when errors are encountered in a previous run. Given that our goal
    # is to have a strictly homogeneous PES data, all other parameters (e.g., ISMEAR, ALGO, etc.) are not inherited.
    inherit_incar: list[str] | bool = (  # type: ignore  # noqa: PGH003
        "LPEAD",
        "NGX",
        "NGY",
        "NGZ",
        "SYMPREC",
        "IMIX",
        "LMAXMIX",
        "KGAMMA",
        "ISYM",
        "NCORE",
        "NPAR",
        "NELMIN",
        "IOPT",
        "NBANDS",
        "KPAR",
        "AMIN",
        "NELMDL",
        "BMIX",
        "AMIX_MAG",
        "BMIX_MAG",
    )
    CONFIG = _load_yaml_config("MatPESStaticSet")

    def __post_init__(self):
        """Validate inputs."""
        super().__post_init__()
        valid_xc_functionals = ("R2SCAN", "PBE", "PBE+U")
        if self.xc_functional.upper() not in valid_xc_functionals:
            raise ValueError(
                f"Unrecognized xc_functional='{self.xc_functional}'. "
                f"Supported exchange-correlation functionals are {valid_xc_functionals}"
            )

        default_potcars = self.CONFIG["PARENT"].replace("PBE", "PBE_").replace("Base", "")  # PBE64Base -> PBE_64
        self.user_potcar_functional = self.user_potcar_functional or default_potcars
        if self.user_potcar_functional.upper() != default_potcars:
            warnings.warn(
                f"{self.user_potcar_functional=} is inconsistent with the recommended {default_potcars}.", UserWarning
            )

        if self.xc_functional.upper() == "R2SCAN":
            self._config_dict["INCAR"].update({"METAGGA": "R2SCAN", "ALGO": "ALL", "GGA": None})
        if self.xc_functional.upper().endswith("+U"):
            self._config_dict["INCAR"]["LDAU"] = True
