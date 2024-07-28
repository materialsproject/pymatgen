"""This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
https://materialsproject.org/dashboard.
"""

from __future__ import annotations

import itertools
import json
import logging
import os
import platform
import sys
import warnings
from typing import TYPE_CHECKING, NamedTuple

import requests
from monty.json import MontyDecoder
from pymatgen.core import SETTINGS
from pymatgen.core import __version__ as PMG_VERSION
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from typing import Callable

    from mp_api.client import MPRester as _MPResterNew
    from pymatgen.core.structure import Structure
    from pymatgen.entries.computed_entries import ComputedStructureEntry
    from pymatgen.ext.matproj_legacy import _MPResterLegacy
    from typing_extensions import Self

logger = logging.getLogger(__name__)

MP_LOG_FILE = os.path.join(os.path.expanduser("~"), ".mprester.log.yaml")


class _MPResterBasic:
    """
    This is Pymatgen's own implement of a MPRester that supports the new MP API. If you are getting your API key from
    the new dashboard of MP, you will need to use this instead of the original MPRester because the new API keys do
    not work with the old MP API (???!).

    The reason why this class exists is because it is the belief of the pymatgen maintainers that access to the MP API
    remains a critical part of pymatgen functionality. However, the new mp-api package is written in a way that
    prevents us from importing the mp-api package because of cyclic dependencies. It is also the opinion of Pymatgen
    maintainers that the implementation of mp-api is too heavy duty for most users (few care about document models,
    etc.). Further, this implementation serves as a simple reference for end developers who want to develop their own
    interfaces to the MP API via the REST urls.

    If you are a power user, feel free to install the mp-api package. All issues regarding that implementation should
    be directed to the maintainers of that repository and not pymatgen. We will support only issues with regards to
    our implementation only.
    """

    def __init__(self, api_key: str | None = None, include_user_agent: bool = True) -> None:
        """
        Args:
            api_key (str): A String API key for accessing the MaterialsProject
                REST interface. Please obtain your API key at
                https://www.materialsproject.org/dashboard. If this is None,
                the code will check if there is a "PMG_MAPI_KEY" setting.
                If so, it will use that environment variable. This makes
                easier for heavy users to simply add this environment variable to
                their setups and MPRester can then be called without any arguments.
            include_user_agent (bool): If True, will include a user agent with the
                HTTP request including information on pymatgen and system version
                making the API request. This helps MP support pymatgen users, and
                is similar to what most web browsers send with each page request.
                Set to False to disable the user agent.
        """
        if api_key is not None:
            self.api_key = api_key
        else:
            self.api_key = SETTINGS.get("PMG_MAPI_KEY", "")
        self.preamble = SETTINGS.get("PMG_MAPI_ENDPOINT", "https://api.materialsproject.org/")

        self.session = requests.Session()
        self.session.headers = {"x-api-key": self.api_key}
        if include_user_agent:
            pymatgen_info = f"pymatgen/{PMG_VERSION}"
            python_info = f"Python/{sys.version.split()[0]}"
            platform_info = f"{platform.system()}/{platform.release()}"
            self.session.headers["user-agent"] = f"{pymatgen_info} ({python_info} {platform_info})"

        # This is a hack, but it fakes most of the functionality of mp-api's summary.search.
        class Summary(NamedTuple):
            search: Callable

        self.summary = Summary(self.summary_search)

    def __getattr__(self, item):
        if item in ("summary", "materials", "thermo"):
            raise AttributeError(
                f"{item} is not an attribute of this implementation of MPRester, which only supports functionality"
                "used by 80% of users. If you are looking for the full functionality MPRester, pls install the mp-api ."
            )

    def __enter__(self) -> Self:
        """Support for "with" context."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Support for "with" context."""
        self.session.close()

    def request(self, sub_url, payload=None, method="GET", mp_decode=True):
        """Helper method to make the requests and perform decoding based on MSONable protocol."""
        response = None
        url = self.preamble + sub_url

        per_page = 1000
        page = 1
        all_data = []
        while True:
            actual_url = f"{url}&_per_page={per_page}&_page={page}"
            try:
                if method == "POST":
                    response = self.session.post(actual_url, data=payload, verify=True)
                else:
                    response = self.session.get(actual_url, params=payload, verify=True)
                if response.status_code in [200, 400]:
                    data = json.loads(response.text, cls=MontyDecoder) if mp_decode else json.loads(response.text)
                else:
                    raise MPRestError(f"REST query returned with error status code {response.status_code}")
                all_data.extend(data["data"])
                if len(data["data"]) < per_page:
                    break
                page += 1
            except Exception as ex:
                msg = f"{ex}. Content: {response.content}" if hasattr(response, "content") else str(ex)
                raise MPRestError(msg)
        return all_data

    def summary_search(self, **kwargs) -> list[dict]:
        """This function mirrors the mp-api's summary.search functionality.

        Args:
            **kwargs: This function only takes kwargs. All kwargs that do not start with an underscore are treated as
                search criteria and those with underscores are treated as params. Example usage:
                MPRester().summary.search(material_ids="mp-19770,mp-19017", _fields="formula_pretty,energy_above_hull")
        """
        criteria = {k: v for k, v in kwargs.items() if not k.startswith("_")}
        params = [f"{k}={v}" for k, v in kwargs.items() if k.startswith("_") and k != "_fields"]
        if "_fields" not in kwargs:
            params.append("_all_fields=True")
        else:
            fields = ",".join(kwargs["_fields"]) if isinstance(kwargs["_fields"], list) else kwargs["_fields"]
            params.extend((f"_fields={fields}", "_all_fields=False"))
        get = "&".join(params)
        logger.info(f"query={get}")
        return self.request(f"materials/summary/?{get}", payload=criteria)

    def get_summary(self, criteria: dict, fields: list | None = None) -> list[dict]:
        """Get  a data corresponding to a criteria.

        Args:
            criteria (dict): Materials Project ID (e.g. mp-1234), e.g. {"formula": "Fe2O3,FeO"}
            fields (list): Fields to query for. If None (the default), all fields are returned.

        Returns:
            List of dict of summary docs.
        """
        get = "_all_fields=True" if fields is None else "_fields=" + ",".join(fields)
        return self.request(f"materials/summary/?{get}", payload=criteria)

    def get_summary_by_material_id(self, material_id: str, fields: list | None = None) -> dict:
        """Get  a data corresponding to a material_id.

        Args:
            material_id (str): Materials Project ID (e.g. mp-1234).
            fields (list): Fields to query for. If None (the default), all fields are returned.

        Returns:
            Dict
        """
        get = "_all_fields=True" if fields is None else "_fields=" + ",".join(fields)
        return self.request(f"materials/summary/?{get}", payload={"material_ids": material_id})[0]

    get_doc = get_summary_by_material_id

    def get_material_ids(self, formula):
        """Get  all materials ids for a formula.

        Args:
            formula (str): A formula (e.g., Fe2O3).

        Returns:
            list[str]: all materials ids.
        """
        return [d["material_id"] for d in self.get_summary({"formula": formula}, fields=["material_id"])]

    # For backwards compatibility and poor spelling.
    get_materials_ids = get_material_ids

    def get_structures(self, chemsys_formula: str, final=True) -> list[Structure]:
        """Get a list of Structures corresponding to a chemical system or formula.

        Args:
            chemsys_formula (str): A chemical system, list of chemical systems
                (e.g., Li-Fe-O, Si-*), or single formula (e.g., Fe2O3, Si*).
            final (bool): Whether to get the final structure, or the list of initial
                (pre-relaxation) structures. Defaults to True.

        Returns:
            List of Structure objects. ([Structure])
        """
        query = f"chemsys={chemsys_formula}" if "-" in chemsys_formula else f"formula={chemsys_formula}"
        prop = "structure" if final else "initial_structure"
        response = self.request(f"materials/summary/?{query}&_all_fields=false&_fields={prop}")

        return [dct[prop] for dct in response]

    def get_structure_by_material_id(self, material_id: str, conventional_unit_cell: bool = False) -> Structure:
        """Get  a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project ID (e.g. mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structures. Defaults to True.
            conventional_unit_cell (bool): Whether to get the standard conventional unit cell

        Returns:
            Structure object.
        """
        prop = "structure"
        response = self.request(f"materials/summary?material_ids={material_id}&_fields={prop}")
        structure = response[0][prop]
        if conventional_unit_cell:
            return SpacegroupAnalyzer(structure).get_conventional_standard_structure()
        return structure

    def get_initial_structures_by_material_id(
        self, material_id: str, conventional_unit_cell: bool = False
    ) -> list[Structure]:
        """Get  a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project ID (e.g. mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structures. Defaults to True.
            conventional_unit_cell (bool): Whether to get the standard conventional unit cell

        Returns:
            Structure object.
        """
        prop = "initial_structures"
        response = self.request(f"materials/summary/{material_id}/?_fields={prop}")
        structures = response[0][prop]
        if conventional_unit_cell:
            return [SpacegroupAnalyzer(s).get_conventional_standard_structure() for s in structures]
        return structures

    def get_entries(
        self,
        criteria,
        compatible_only=True,
        inc_structure=None,
        property_data=None,
        conventional_unit_cell=False,
        sort_by_e_above_hull=False,
    ):
        """Get  a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id or full criteria.

        Args:
            criteria: Chemsys, formula, or mp-id.
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProject2020Compatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            inc_structure (str): If None, entries returned are
                ComputedEntries. If inc_structure="initial",
                ComputedStructureEntries with initial structures are returned.
                Otherwise, ComputedStructureEntries with final structures
                are returned.
            property_data (list): Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell
            sort_by_e_above_hull (bool): Whether to sort the list of entries by
                e_above_hull (will query e_above_hull as a property_data if True).

        Returns:
            List of ComputedStructureEntry objects.
        """
        if criteria.startswith("mp-"):
            query = f"material_ids={criteria}"
        elif "-" in criteria:
            query = f"chemsys={criteria}"
        else:
            query = f"formula={criteria}"

        entries = []
        response = self.request(f"materials/thermo/?_fields=entries&{query}")
        for dct in response:
            entries.extend(dct["entries"].values())

        if compatible_only:
            from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

            # suppress the warning about missing oxidation states
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Failed to guess oxidation states.*")
                entries = MaterialsProject2020Compatibility().process_entries(entries, clean=True)
        return list(set(entries))

    def get_entry_by_material_id(self, material_id: str, *args, **kwargs) -> ComputedStructureEntry:
        r"""Get  a ComputedEntry corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g. mp-1234).
            *args: Pass-through to get_entries.
            **kwargs: Pass-through to get_entries.

        Returns:
            ComputedStructureEntry object.
        """
        return self.get_entries(material_id, *args, **kwargs)[0]

    def get_entries_in_chemsys(self, elements, *args, **kwargs):
        """
        Helper method to get a list of ComputedEntries in a chemical system. For example, elements = ["Li", "Fe", "O"]
        will return a list of all entries in the Li-Fe-O chemical system, i.e., all LixOy, FexOy, LixFey, LixFeyOz,
        Li, Fe and O phases. Extremely useful for creating phase diagrams of entire chemical systems.

        Args:
            elements (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g. "Li-Fe-O" or List of element
                symbols, e.g. ["Li", "Fe", "O"].
            *args: Pass-through to get_entries.
            **kwargs: Pass-through to get_entries.

        Returns:
            List of ComputedEntries.
        """
        chemsys = []
        for i in range(1, len(elements) + 1):
            for els in itertools.combinations(elements, i):
                chemsys.append("-".join(sorted(els)))
        criteria = ",".join(chemsys)

        return self.get_entries(criteria, *args, **kwargs)

    def get_phonon_bandstructure_by_material_id(self, material_id: str):
        """Get phonon bandstructure by material_id.

        Args:
            material_id (str): Materials Project material_id

        Returns:
            PhononBandStructureSymmLine: A phonon band structure.
        """
        prop = "ph_bs"
        response = self.request(f"materials/phonon/?material_ids={material_id}&_fields={prop}")
        return response[0][prop]

    def get_phonon_dos_by_material_id(self, material_id: str):
        """Get phonon density of states by material_id.

        Args:
            material_id (str): Materials Project material_id

        Returns:
            CompletePhononDos: A phonon DOS object.
        """
        prop = "ph_dos"
        response = self.request(f"materials/phonon/?material_ids={material_id}&_fields={prop}")
        return response[0][prop]


class MPRester:
    """A class to conveniently interface with the new and legacy Materials Project REST interface.

    The recommended way to use MPRester is as a context manager to ensure
    that sessions are properly closed after usage:

        with MPRester("API_KEY") as mpr:
            docs = mpr.call_some_method()

    MPRester uses the "requests" package, which provides HTTP connection
    pooling. All connections are made via https for security.

    For more advanced uses of the Materials API, please consult the API
    documentation at https://materialsproject.org/api and https://docs.materialsproject.org.

    This class handles the transition between old and new MP API, making it easy to switch between them
    by passing a new (length 32) or old (15 <= length <= 17) API key. See https://docs.materialsproject.org
    for which API to use.
    """

    def __new__(cls, *args, **kwargs) -> _MPResterNew | _MPResterBasic | _MPResterLegacy:  # type: ignore[misc]
        """
        Args:
           *args: Pass through to either legacy or new MPRester.
           **kwargs: Pass through to either legacy or new MPRester.
        """
        api_key = args[0] if len(args) > 0 else None

        if api_key is None:
            api_key = kwargs.get("api_key", SETTINGS.get("PMG_MAPI_KEY"))
            kwargs["api_key"] = api_key

        if not api_key:
            raise ValueError("Please supply an API key. See https://materialsproject.org/api for details.")

        if len(api_key) != 32:
            from pymatgen.ext.matproj_legacy import _MPResterLegacy

            return _MPResterLegacy(*args, **kwargs)

        try:
            from mp_api.client import MPRester as _MPResterNew

            return _MPResterNew(*args, **kwargs)
        except Exception:
            return _MPResterBasic(*args, **kwargs)


class MPRestError(Exception):
    """Exception class for legacy MPRestAdaptor. Raised when query is malformed."""
