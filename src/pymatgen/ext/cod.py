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

import re
import subprocess
import warnings
from shutil import which

import requests
from monty.dev import requires

from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure


class COD:
    """An interface to the Crystallography Open Database."""

    url = "www.crystallography.net"

    def query(self, sql: str) -> str:
        """Perform a query.

        Args:
            sql: SQL string

        Returns:
            Response from SQL query.
        """
        response = subprocess.check_output(["mysql", "-u", "cod_reader", "-h", self.url, "-e", sql, "cod"])
        return response.decode("utf-8")

    @requires(which("mysql"), "mysql must be installed to use this query.")
    def get_cod_ids(self, formula) -> list[int]:
        """Query the COD for all cod ids associated with a formula. Requires
        mysql executable to be in the path.

        Args:
            formula (str): Formula.

        Returns:
            List of cod ids.
        """
        # TODO: Remove dependency on external mysql call. MySQL-python package does not support Py3!

        # Standardize formula to the version used by COD
        cod_formula = Composition(formula).hill_formula
        sql = f'select file from data where formula="- {cod_formula} -"'  # noqa: S608
        text = self.query(sql).split("\n")
        cod_ids = []
        for line in text:
            if match := re.search(r"(\d+)", line):
                cod_ids.append(int(match[1]))
        return cod_ids

    def get_structure_by_id(self, cod_id: int, timeout: int = 600, **kwargs) -> Structure:
        """Query the COD for a structure by id.

        Args:
            cod_id (int): COD id.
            timeout (int): Timeout for the request in seconds. Default = 600.
            kwargs: All kwargs supported by Structure.from_str.

        Returns:
            A Structure.
        """
        response = requests.get(f"https://{self.url}/cod/{cod_id}.cif", timeout=timeout)
        return Structure.from_str(response.text, fmt="cif", **kwargs)

    @requires(which("mysql"), "mysql must be installed to use this query.")
    def get_structure_by_formula(self, formula: str, **kwargs) -> list[dict[str, str | int | Structure]]:
        """Query the COD for structures by formula. Requires mysql executable to
        be in the path.

        Args:
            formula (str): Chemical formula.
            kwargs: All kwargs supported by Structure.from_str.

        Returns:
            A list of dict of the format [{"structure": Structure, "cod_id": int, "sg": "P n m a"}]
        """
        structures: list[dict[str, str | int | Structure]] = []
        sql = f'select file, sg from data where formula="- {Composition(formula).hill_formula} -"'  # noqa: S608
        text = self.query(sql).split("\n")
        text.pop(0)
        for line in text:
            if line.strip():
                cod_id, sg = line.split("\t")
                response = requests.get(f"https://{self.url}/cod/{cod_id.strip()}.cif", timeout=60)
                try:
                    struct = Structure.from_str(response.text, fmt="cif", **kwargs)
                    structures.append({"structure": struct, "cod_id": int(cod_id), "sg": sg})
                except Exception:
                    warnings.warn(f"\nStructure.from_str failed while parsing CIF file:\n{response.text}")
                    raise

        return structures
