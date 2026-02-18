"""This module provides classes to interface with the Crystallography Open
Database. If you use data from the COD, please cite the following works (as
stipulated by the COD developers).

    Merkys, A., Vaitkus, A., Butkus, J., Okulič-Kazarinas, M., Kairys, V. &
    Gražulis, S. (2016) "COD::CIF::Parser: an error-correcting CIF parser for
    the Perl language". Journal of Applied Crystallography 49.

    Gražulis, S., Merkys, A., Vaitkus, A. & Okulič-Kazarinas, M. (2015)
    "Computing stoichiometric molecular composition from crystal structures".
    Journal of Applied Crystallography 48, 85-91.

    Gražulis, S., Daškevič, A., Merkys, A., Chateigner, D., Lutterotti, L.,
    Quirós, M., Serebryanaya, N. R., Moeck, P., Downs, R. T. & LeBail, A.
    (2012) "Crystallography Open Database (COD): an open-access collection of
    crystal structures and platform for world-wide collaboration". Nucleic
    Acids Research 40, D420-D427.

    Grazulis, S., Chateigner, D., Downs, R. T., Yokochi, A. T., Quiros, M.,
    Lutterotti, L., Manakova, E., Butkus, J., Moeck, P. & Le Bail, A. (2009)
    "Crystallography Open Database - an open-access collection of crystal
    structures". J. Appl. Cryst. 42, 726-729.

    Downs, R. T. & Hall-Wallace, M. (2003) "The American Mineralogist Crystal
    Structure Database". American Mineralogist 88, 247-250.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import requests

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure

if TYPE_CHECKING:
    from typing import Literal


class COD:
    """An interface to the Crystallography Open Database.

    Reference:
        https://wiki.crystallography.net/RESTful_API/
    """

    def __init__(self, timeout: int = 60):
        """Initialize the COD class.

        Args:
            timeout (int): request timeout in seconds.
        """
        self.timeout = timeout
        self.url = "https://www.crystallography.net"
        self.api_url = f"{self.url}/cod/result"

    def get_cod_ids(self, formula: str) -> list[int]:
        """Query the COD for all COD IDs associated with a formula.

        Args:
            formula (str): The formula to request
        """
        # Use hill_formula format as per COD request
        cod_formula = Composition(formula).hill_formula

        # Set up query parameters
        params = {"formula": cod_formula, "format": "json"}

        response = requests.get(self.api_url, params=params, timeout=self.timeout)

        # Raise an exception if the request fails
        response.raise_for_status()

        return [int(entry["file"]) for entry in response.json()]

    def get_structure_by_id(self, cod_id: int, timeout: int | None = None, **kwargs) -> Structure:
        """Query the COD for a structure by ID.

        Args:
            cod_id (int): COD ID.
            timeout (int): DEPRECATED. request timeout in seconds.
            kwargs: kwargs passed to Structure.from_str.

        Returns:
            A Structure.
        """
        # TODO: remove timeout arg and use class level timeout after 2025-10-17
        if timeout is not None:
            warnings.warn(
                "separate timeout arg is deprecated, please use class level timeout", DeprecationWarning, stacklevel=2
            )
        timeout = timeout or self.timeout

        response = requests.get(f"{self.url}/cod/{cod_id}.cif", timeout=timeout)
        return Structure.from_str(response.text, fmt="cif", **kwargs)

    def get_structure_by_formula(
        self,
        formula: str,
        **kwargs,
    ) -> list[dict[Literal["structure", "cod_id", "sg"], str | int | Structure]]:
        """Query the COD for structures by formula.

        Args:
            formula (str): Chemical formula.
            kwargs: All kwargs supported by Structure.from_str.

        Returns:
            A list of dict of: {"structure": Structure, "cod_id": int, "sg": "P n m a"}
        """
        # Prepare the query parameters
        params = {
            "formula": Composition(formula).hill_formula,
            "format": "json",
        }

        response = requests.get(self.api_url, params=params, timeout=self.timeout)
        response.raise_for_status()

        structures: list[dict[Literal["structure", "cod_id", "sg"], str | int | Structure]] = []

        # Parse the JSON response
        for entry in response.json():
            cod_id = entry["file"]
            sg = entry.get("sg")

            try:
                struct = self.get_structure_by_id(cod_id, **kwargs)
                structures.append({"structure": struct, "cod_id": int(cod_id), "sg": sg})

            except Exception:
                warnings.warn(
                    f"Structure.from_str failed while parsing CIF file for COD ID {cod_id}",
                    stacklevel=2,
                )
                raise

        return structures
