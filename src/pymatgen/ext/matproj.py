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
from functools import partial
from typing import TYPE_CHECKING, NamedTuple

import orjson
import requests
from monty.json import MontyDecoder

from pymatgen.core import SETTINGS
from pymatgen.core import __version__ as PMG_VERSION
from pymatgen.core.composition import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from typing_extensions import Self

    from pymatgen.core.structure import Structure
    from pymatgen.entries.computed_entries import ComputedStructureEntry

logger = logging.getLogger(__name__)

MP_LOG_FILE = os.path.join(os.path.expanduser("~"), ".mprester.log.yaml")

CHUNK_SIZE = 200


class MPRester:
    """
    Pymatgen's implementation of MPRester. Unlike mp-api, this implementation mirrors the exact MP-API end points
    without modification. You just need to refer to https://api.materialsproject.org/docs and use the field names
    exactly. No need to deal with strange renames of various fields. Featurity parity is close to 100% with mp-api.

    Furthermore, we support both the mp-api as well as a simplified syntax. E.g., to query for a summary, you can use
    mpr.summary.search(material_ids="mp-1234") or mpr.materials.summary.search(material_ids="mp-1234").

    If you are a power user that requires some esoteric feature not covered, feel free to install the mp-api package.
    All issues regarding that implementation should be directed to the maintainers of that repository and not
    pymatgen. We will support only issues pertaining to our implementation only.

    Attributes:
    :ivar api_key: API key for authenticating requests to the Materials Project API.
    :type api_key: str
    :ivar preamble: Base endpoint URL for the Materials Project API.
    :type preamble: str
    :ivar session: HTTP session object for managing API requests.
    :type session: requests.Session
    :ivar materials: Placeholder object for dynamically adding endpoints related to materials.
    :type materials: Any
    """

    MATERIALS_DOCS = (
        "summary",
        "core",
        "elasticity",
        "phonon",
        "eos",
        "similarity",
        "xas",
        "grain_boundaries",
        "electronic_structure",
        "tasks",
        "substrates",
        "surface_properties",
        "robocrys",
        "synthesis",
        "magnetism",
        "insertion_electrodes",
        "conversion_electrodes",
        "oxidation_states",
        "provenance",
        "alloys",
        "absorption",
        "chemenv",
        "bonds",
        "piezoelectric",
        "dielectric",
    )

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

        if len(self.api_key) != 32:
            raise ValueError(
                "Invalid or old API key. Please obtain an updated API key at https://materialsproject.org/dashboard."
            )

        self.preamble = SETTINGS.get("PMG_MAPI_ENDPOINT", "https://api.materialsproject.org/")

        self.session = requests.Session()
        self.session.headers = {"x-api-key": self.api_key}
        if include_user_agent:
            pymatgen_info = f"pymatgen/{PMG_VERSION}"
            python_info = f"Python/{sys.version.split()[0]}"
            platform_info = f"{platform.system()}/{platform.release()}"
            self.session.headers["user-agent"] = f"{pymatgen_info} ({python_info} {platform_info})"

        # This is a hack, but it fakes most of the functionality of mp-api's summary.search.
        class Search(NamedTuple):
            search: Callable

        class Materials:
            pass

        self.materials = Materials()

        for doc in MPRester.MATERIALS_DOCS:
            setattr(self, doc, Search(partial(self.search, doc)))
            setattr(self.materials, doc, Search(partial(self.search, doc)))

    def __getattr__(self, item):
        raise AttributeError(
            f"{item} is not an attribute of this implementation of MPRester, which only supports functionality "
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
                    data = json.loads(response.text, cls=MontyDecoder) if mp_decode else orjson.loads(response.text)
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

    def search(self, doc, **kwargs) -> list[dict]:
        """
        Queries a Materials PI end point doc. A notable difference with the mp-api's implementation is that this uses
        the web API to do searches. So the keywords follow the actual API spec, which is as it should be. For
        instance, number of sites is `nsites` and number of elements is `nelements`. The mp-api package has this
        weird renaming that maps `num_elements` to `nelements` and `num_sites` to  `nsites`.

        Parameters:
        - **kwargs: keyword arguments for filtering materials. Fields that do not start with underscores are
        filters, while those that start with underscores are fields to retrieve. Possible filters include:
          - _fields (optional): list of fields to retrieve for each material
          - Other parameters: filter criteria, where each parameter key corresponds to the field to filter and the
            parameter value corresponds to the filter value

        Returns:
        - list of dictionaries, each dictionary representing a material retrieved based on the filtering criteria
        """

        def comma_cat(val):
            return ",".join(val) if isinstance(val, list) else val

        criteria = {k: comma_cat(v) for k, v in kwargs.items() if not k.startswith("_")}
        params = [f"{k}={comma_cat(v)}" for k, v in kwargs.items() if k.startswith("_") and k != "_fields"]
        if "_fields" not in kwargs:
            params.append("_all_fields=True")
        else:
            fields = comma_cat(kwargs["_fields"])
            params.extend((f"_fields={fields}", "_all_fields=False"))
        get = "&".join(params)
        logger.info(f"query={get}")
        return self.request(f"materials/{doc}/?{get}", payload=criteria)

    def summary_search(self, **kwargs) -> list[dict]:
        """
        Mirrors mp-api's mpr.materials.summary.search functionality. Searches for materials based on the specified
        criteria. A notable difference with the mp-api's implementation is that this uses the web API to do searches.
        So the keywords follow the actual API spec, which is as it should be. For instance, number of sites is `nsites`
        and number of elements is `nelements`. The mp-api package has this weird renaming that maps `num_elements` to
        `nelements` and `num_sites` to  `nsites`.

        Parameters:
        - **kwargs: keyword arguments for filtering materials. Fields that do not start with underscores are
        filters, while those that start with underscores are fields to retrieve. Possible filters include:
          - _fields (optional): list of fields to retrieve for each material
          - Other parameters: filter criteria, where each parameter key corresponds to the field to filter and the
            parameter value corresponds to the filter value

        Returns:
        - list of dictionaries, each dictionary representing a material retrieved based on the filtering criteria

        """
        return self.search("summary", **kwargs)

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
        criteria: str | Sequence[str],
        *,
        compatible_only: bool = True,
        property_data: Sequence[str] | None = None,
        summary_data: Sequence[str] | None = None,
        **kwargs,
    ):
        """Get a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id or full criteria.

        Args:
            criteria: Chemsys, formula, or mp-id.
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProject2020Compatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            property_data (list): Specify additional properties to include in entry.data from /materials/thermo
                endpoint of the API.
            summary_data (list): Specify additional properties to include in entry.data from /materials/summary
                endpoint of the API. Note that unlike property_data, summary_data refers to the "best" calculation done
                by MP. There are no guarantees that the summary_data will be consistent with the entry or property
                data since the data can come from say a r2SCAN calculation but the entry is from a GGA calculation.
                The data will be reported in the entry.data["summary"] dictionary.
            **kwargs: Used to catch deprecated kwargs.

        Returns:
            List of ComputedStructureEntry objects.
        """
        if set(kwargs.keys()).intersection({"inc_structure", "conventional_unit_cell", "sort_by_e_above_hull"}):
            warnings.warn(
                "The inc_structure, conventional_unit_cell, and sort_by_e_above_hull arguments are deprecated. "
                "These arguments have no effect and will be removed in 2026.1.1.",
                DeprecationWarning,
                stacklevel=2,
            )

        def proc_crit(val):
            if val.startswith("mp-"):
                return val
            if "-" in val:
                return "-".join(sorted(val.split("-")))
            comp = Composition(val)
            if len(comp) == 1:
                return val
            return comp.reduced_formula

        if isinstance(criteria, str):
            criteria = ",".join([proc_crit(c) for c in criteria.split(",")])
        else:
            criteria = ",".join([proc_crit(c) for c in criteria])

        if criteria.startswith("mp-"):
            query = f"material_ids={criteria}"
        elif "-" in criteria:
            query = f"chemsys={criteria}"
        else:
            query = f"formula={criteria}"

        entries = []
        fields = ["entries", *property_data] if property_data is not None else ["entries"]
        response = self.request(f"materials/thermo/?_fields={','.join(fields)}&{query}")
        for dct in response:
            for e in dct["entries"].values():
                if property_data:
                    for prop in property_data:
                        e.data[prop] = dct[prop]
                entries.append(e)

        if compatible_only:
            from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

            # suppress the warning about missing oxidation states
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Failed to guess oxidation states.*")
                entries = MaterialsProject2020Compatibility().process_entries(entries, clean=True)

        if summary_data and len(entries) > 0:
            for i in range(0, len(entries), CHUNK_SIZE):
                chunked_entries = entries[i : i + CHUNK_SIZE]
                mids = [e.data["material_id"] for e in chunked_entries]
                edata = self.search(
                    "summary",
                    material_ids=mids,
                    _fields=[*summary_data, "material_id"],
                )
                mapped_data = {d["material_id"]: {k: v for k, v in d.items() if k != "material_id"} for d in edata}
                for e in chunked_entries:
                    e.data["summary"] = mapped_data[e.data["material_id"]]

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

    def get_entries_in_chemsys(self, elements: str | list[str], *args, **kwargs):
        """
        Helper method to get a list of ComputedEntries in a chemical system. For example, elements = ["Li", "Fe", "O"]
        will return a list of all entries in the Li-Fe-O chemical system, i.e., all LixOy, FexOy, LixFey, LixFeyOz,
        Li, Fe and O phases. Extremely useful for creating phase diagrams of entire chemical systems.

        Args:
            elements (str | list[str]): Chemical system string comprising element
                symbols separated by dashes, e.g. "Li-Fe-O" or List of element
                symbols, e.g. ["Li", "Fe", "O"].
            *args: Pass-through to get_entries.
            **kwargs: Pass-through to get_entries.

        Returns:
            List of ComputedEntries.
        """
        if isinstance(elements, str):
            elements = elements.split("-")
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


class MPRestError(Exception):
    """Exception class for legacy MPRestAdaptor. Raised when query is malformed."""
