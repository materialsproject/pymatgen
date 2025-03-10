"""This module provides classes to interface with the Materials Project REST
API v1 to enable the creation of data structures and pymatgen objects using
Materials Project data.
"""

from __future__ import annotations

import itertools
import json
import logging
import math
import os
import platform
import re
import sys
import warnings
from enum import Enum, unique
from time import sleep
from typing import TYPE_CHECKING

import requests
from monty.json import MontyDecoder, MontyEncoder
from ruamel.yaml import YAML
from tqdm import tqdm

from pymatgen.core import SETTINGS, Composition, Element, Structure
from pymatgen.core import __version__ as PMG_VERSION
from pymatgen.core.surface import get_symmetrically_equivalent_miller_indices
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any, Literal

    from typing_extensions import Self

    from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
    from pymatgen.phonon.dos import CompletePhononDos

logger = logging.getLogger(__name__)
MP_LOG_FILE = os.path.join(os.path.expanduser("~"), ".mprester.log.yaml")


@unique
class TaskType(Enum):
    """task types available in legacy MP data."""

    GGA_OPT = "GGA Structure Optimization"
    GGAU_OPT = "GGA+U Structure Optimization"
    SCAN_OPT = "SCAN Structure Optimization"
    GGA_LINE = "GGA NSCF Line"
    GGAU_LINE = "GGA+U NSCF Line"
    GGA_UNIFORM = "GGA NSCF Uniform"
    GGAU_UNIFORM = "GGA+U NSCF Uniform"
    GGA_STATIC = "GGA Static"
    GGAU_STATIC = "GGA+U Static"
    GGA_STATIC_DIEL = "GGA Static Dielectric"
    GGAU_STATIC_DIEL = "GGA+U Static Dielectric"
    GGA_DEF = "GGA Deformation"
    GGAU_DEF = "GGA+U Deformation"
    LDA_STATIC_DIEL = "LDA Static Dielectric"


class _MPResterLegacy:
    """A class to conveniently interface with the Materials Project REST interface.
    The recommended way to use MPRester is with the "with" context manager to ensure
    sessions are properly closed after usage.

        with MPRester("API_KEY") as mpr:
            mpr.some_method()

    MPRester uses the "requests" package, which provides for HTTP connection
    pooling. All connections are made via https for security.

    For more advanced uses of the legacy Materials API, please consult the API
    documentation at https://github.com/materialsproject/mapidoc.

    Note that this class is for the *legacy* API. Upcoming changes to the
    Materials Project api are described at https://materialsproject.org/api.
    """

    supported_properties = (
        "energy",
        "energy_per_atom",
        "volume",
        "formation_energy_per_atom",
        "nsites",
        "unit_cell_formula",
        "pretty_formula",
        "is_hubbard",
        "elements",
        "nelements",
        "e_above_hull",
        "hubbards",
        "is_compatible",
        "spacegroup",
        "task_ids",
        "band_gap",
        "density",
        "icsd_id",
        "icsd_ids",
        "cif",
        "total_magnetization",
        "material_id",
        "oxide_type",
        "tags",
        "elasticity",
    )

    supported_task_properties = (
        "energy",
        "energy_per_atom",
        "volume",
        "formation_energy_per_atom",
        "nsites",
        "unit_cell_formula",
        "pretty_formula",
        "is_hubbard",
        "elements",
        "nelements",
        "e_above_hull",
        "hubbards",
        "is_compatible",
        "spacegroup",
        "band_gap",
        "density",
        "icsd_id",
        "cif",
    )

    def __init__(
        self,
        api_key: str | None = None,
        endpoint: str | None = None,
        notify_db_version: bool = True,
        include_user_agent: bool = True,
    ) -> None:
        """
        Args:
            api_key (str): A String API key for accessing the MaterialsProject
                REST interface. Please obtain your API key at
                https://materialsproject.org/dashboard. If this is None,
                the code will check if there is a "PMG_MAPI_KEY" setting.
                If so, it will use that environment variable. This makes
                easier for heavy users to simply add this environment variable to
                their setups and MPRester can then be called without any arguments.
            endpoint (str): Url of endpoint to access the MaterialsProject REST
                interface. Defaults to the standard Materials Project REST
                address at "https://legacy.materialsproject.org/rest/v2", but
                can be changed to other urls implementing a similar interface.
            notify_db_version (bool): If True, the current MP database version will
                be retrieved and logged locally in the ~/.pmgrc.yaml. If the database
                version changes, you will be notified. The current database version is
                also printed on instantiation. These local logs are not sent to
                materialsproject.org and are not associated with your API key, so be
                aware that a notification may not be presented if you run MPRester
                from multiple computing environments.
            include_user_agent (bool): If True, will include a user agent with the
                HTTP request including information on pymatgen and system version
                making the API request. This helps MP support pymatgen users, and
                is similar to what most web browsers send with each page request.
                Set to False to disable the user agent.
        """
        warnings.warn(
            "You are using the legacy MPRester. This version of the MPRester will no longer be updated. "
            "To access the latest data with the new MPRester, obtain a new API key from "
            "https://materialsproject.org/api and consult the docs at https://docs.materialsproject.org/ "
            "for more information.",
            FutureWarning,
            stacklevel=2,
        )
        if api_key is not None:
            self.api_key = api_key
        else:
            self.api_key = SETTINGS.get("PMG_MAPI_KEY", "")

        if endpoint is not None:
            self.preamble = endpoint
        else:
            self.preamble = SETTINGS.get("PMG_MAPI_ENDPOINT", "https://legacy.materialsproject.org/rest/v2")

        if self.preamble != "https://legacy.materialsproject.org/rest/v2":
            warnings.warn(f"Non-default endpoint used: {self.preamble}", stacklevel=2)

        self.session = requests.Session()
        self.session.headers = {"x-api-key": self.api_key}
        if include_user_agent:
            pymatgen_info = f"pymatgen/{PMG_VERSION}"
            python_info = f"Python/{sys.version.split()[0]}"
            platform_info = f"{platform.system()}/{platform.release()}"
            self.session.headers["user-agent"] = f"{pymatgen_info} ({python_info} {platform_info})"

        if notify_db_version:
            yaml = YAML()
            db_version = self.get_database_version()
            logger.debug(f"Connection established to Materials Project database, version {db_version}.")

            try:
                with open(MP_LOG_FILE, encoding="utf-8") as file:
                    dct = dict(yaml.load(file)) or {}
            except (OSError, TypeError):
                # TypeError: 'NoneType' object is not iterable occurs if MP_LOG_FILE exists but is empty
                dct = {}

            if "MAPI_DB_VERSION" not in dct:
                dct["MAPI_DB_VERSION"] = {"LOG": {}, "LAST_ACCESSED": None}
            else:
                # ensure data is parsed as dict, rather than ordered dict,
                # due to change in YAML parsing behavior
                dct["MAPI_DB_VERSION"] = dict(dct["MAPI_DB_VERSION"])

            if "LOG" in dct["MAPI_DB_VERSION"]:
                dct["MAPI_DB_VERSION"]["LOG"] = dict(dct["MAPI_DB_VERSION"]["LOG"])

            # store a log of what database versions are being connected to
            if db_version not in dct["MAPI_DB_VERSION"]["LOG"]:
                dct["MAPI_DB_VERSION"]["LOG"][db_version] = 1
            else:
                dct["MAPI_DB_VERSION"]["LOG"][db_version] += 1

            # alert user if DB version changed
            last_accessed = dct["MAPI_DB_VERSION"]["LAST_ACCESSED"]
            if last_accessed and last_accessed != db_version:
                warnings.warn(
                    f"This database version has changed from the database last accessed ({last_accessed}).\n"
                    f"Please see release notes on materialsproject.org for information about what has changed.",
                    stacklevel=2,
                )
            dct["MAPI_DB_VERSION"]["LAST_ACCESSED"] = db_version

            # write out new database log if possible
            # base Exception is not ideal (perhaps a PermissionError, etc.) but this is not critical
            # and should be allowed to fail regardless of reason
            try:
                with open(MP_LOG_FILE, mode="w", encoding="utf-8") as file:
                    yaml.dump(dct, file)
            except Exception:
                pass

    def __enter__(self) -> Self:
        """Support for "with" context."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Support for "with" context."""
        self.session.close()

    def _make_request(
        self,
        sub_url: str,
        payload: Any = None,
        method: Literal["GET", "POST", "PUT", "DELETE"] = "GET",
        mp_decode: bool = True,
    ) -> Any:
        response = None
        url = f"{self.preamble}{sub_url}"
        try:
            if method == "POST":
                response = self.session.post(url, data=payload, verify=True)
            else:
                response = self.session.get(url, params=payload, verify=True)
            if response.status_code in [200, 400]:
                data = json.loads(response.text, cls=MontyDecoder) if mp_decode else json.loads(response.text)
                if data["valid_response"]:
                    if data.get("warning"):
                        warnings.warn(data["warning"], stacklevel=2)
                    return data["response"]
                raise MPRestError(data["error"])

            raise MPRestError(f"REST query returned with error status code {response.status_code}")

        except Exception as exc:
            msg = f"{exc}. Content: {getattr(response, 'content', str(exc))}"
            raise MPRestError(msg)

    def get_database_version(self) -> str:
        """The Materials Project database is periodically updated and has a
        database version associated with it. When the database is updated,
        consolidated data (information about "a material") may and does
        change, while calculation data about a specific calculation task
        remains unchanged and available for querying via its task_id.

        The database version is set as a date in the format YYYY-MM-DD,
        where "-DD" may be optional. An additional numerical suffix
        might be added if multiple releases happen on the same day.

        Returns:
            str: database version
        """
        dct = self._make_request("/api_check")
        return dct["version"]["db"]

    def get_materials_id_from_task_id(self, task_id) -> str:
        """Get a new MP materials id from a task id (which can be
        equivalent to an old materials id).

        Args:
            task_id (str): A task id.

        Returns:
            str: An MP material id.
        """
        return self._make_request(f"/materials/mid_from_tid/{task_id}")

    def get_materials_id_references(self, material_id):
        """Get all references for a materials id.

        Args:
            material_id (str): A material id.

        Returns:
            str: A BibTeX formatted string.
        """
        return self._make_request(f"/materials/{material_id}/refs")

    def get_data(self, chemsys_formula_id, data_type="vasp", prop=""):
        """Flexible method to get any data using the Materials Project REST
        interface. Generally used by other methods for more specific queries.

        Format of REST return is *always* a list of dict (regardless of the
        number of pieces of data returned. The general format is as follows:

        [{"material_id": material_id, "property_name" : value}, ...]

        This is generally a call to
        https://materialsproject.org/rest/v2/materials/vasp/<prop>.
        See https://github.com/materialsproject/mapidoc for details.

        Args:
            chemsys_formula_id (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).
            data_type (str): Type of data to return. Currently can either be
                "vasp" or "exp".
            prop (str): Property to be obtained. Should be one of the
                MPRester.supported_task_properties. Leave as empty string for a
                general list of useful properties.
        """
        sub_url = f"/materials/{chemsys_formula_id}/{data_type}"
        if prop:
            sub_url += "/" + prop
        return self._make_request(sub_url)

    def get_materials_ids(self, chemsys_formula):
        """Get all materials ids for a formula or chemsys.

        Args:
            chemsys_formula (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3).

        Returns:
            list[str]: all materials ids.
        """
        return self._make_request(f"/materials/{chemsys_formula}/mids", mp_decode=False)

    # For backwards compatibility.
    get_material_id = get_materials_ids

    def get_doc(self, materials_id):
        """Get the entire data document for one materials id. Use this judiciously.

        REST Endpoint: https://materialsproject.org/materials/<mp-id>/doc.

        Args:
            materials_id (str): e.g. mp-1143 for Al2O3

        Returns:
            Dict of JSON document of all data that is displayed on a materials
            details page.
        """
        return self._make_request(f"/materials/{materials_id}/doc", mp_decode=False)

    def get_xas_data(self, material_id, absorbing_element):
        """Get X-ray absorption spectroscopy data for absorbing element in the
        structure corresponding to a material_id. Only X-ray Absorption Near Edge
        Structure (XANES) for K-edge is supported.

        REST Endpoint:
        https://materialsproject.org/materials/<mp-id>/xas/<absorbing_element>.

        Args:
            material_id (str): e.g. mp-1143 for Al2O3
            absorbing_element (str): The absorbing element in the corresponding
                structure. e.g. Al in Al2O3
        """
        element_list = self.get_data(material_id, prop="elements")[0]["elements"]
        if absorbing_element not in element_list:
            raise ValueError(
                f"{absorbing_element} element not contained in corresponding structure with mp_id: {material_id}"
            )
        data = self._make_request(
            f"/materials/{material_id}/xas/{absorbing_element}",
            mp_decode=False,
        )
        return data[0]

    def get_task_data(self, chemsys_formula_id, prop=""):
        """Flexible method to get any data using the Materials Project REST
        interface. Generally used by other methods for more specific queries.
        Unlike the :func:`get_data`_, this method queries the task collection
        for specific run information.

        Format of REST return is *always* a list of dict (regardless of the
        number of pieces of data returned. The general format is as follows:

        [{"material_id": material_id, "property_name" : value}, ...]

        Args:
            chemsys_formula_id (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).
            prop (str): Property to be obtained. Should be one of the
                MPRester.supported_properties. Leave as empty string for a
                general list of useful properties.
        """
        sub_url = f"/tasks/{chemsys_formula_id}"
        if prop:
            sub_url += "/" + prop
        return self._make_request(sub_url)

    def get_structures(self, chemsys_formula_id, final=True):
        """Get a list of Structures corresponding to a chemical system, formula,
        or materials_id.

        Args:
            chemsys_formula_id (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.

        Returns:
            List of Structure objects.
        """
        prop = "final_structure" if final else "initial_structure"
        data = self.get_data(chemsys_formula_id, prop=prop)
        return [d[prop] for d in data]

    def find_structure(self, filename_or_structure):
        """Find matching structures on the Materials Project site.

        Args:
            filename_or_structure: filename or Structure object

        Returns:
            A list of matching materials project ids for structure.

        Raises:
            MPRestError
        """
        if isinstance(filename_or_structure, str):
            struct = Structure.from_file(filename_or_structure)
        elif isinstance(filename_or_structure, Structure):
            struct = filename_or_structure
        else:
            raise MPRestError("Provide filename or Structure object.")
        payload = {"structure": json.dumps(struct.as_dict(), cls=MontyEncoder)}
        response = self.session.post(f"{self.preamble}/find_structure", data=payload)
        if response.status_code in [200, 400]:
            response = json.loads(response.text, cls=MontyDecoder)
            if response["valid_response"]:
                return response["response"]
            raise MPRestError(response["error"])
        raise MPRestError(f"REST error with status code {response.status_code} and error {response.text}")

    def get_entries(
        self,
        chemsys_formula_id_criteria: str | dict[str, Any],
        compatible_only: bool = True,
        inc_structure: bool | Literal["initial"] | None = None,
        property_data: list[str] | None = None,
        conventional_unit_cell: bool = False,
        sort_by_e_above_hull: bool = False,
    ) -> list[ComputedEntry]:
        """Get a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id or full criteria.

        Args:
            chemsys_formula_id_criteria (str/dict): A chemical system
                (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id
                (e.g., mp-1234) or full Mongo-style dict criteria.
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
            List of ComputedEntry or ComputedStructureEntry objects.
        """
        # TODO: This is a very hackish way of doing this. It should be fixed
        # on the REST end.
        params = [
            "run_type",
            "is_hubbard",
            "pseudo_potential",
            "hubbards",
            "potcar_symbols",
            "oxide_type",
        ]
        props = ["energy", "unit_cell_formula", "task_id", *params]
        if sort_by_e_above_hull:
            if property_data and "e_above_hull" not in property_data:
                property_data.append("e_above_hull")
            elif not property_data:
                property_data = ["e_above_hull"]
        if property_data:
            props += property_data
        if inc_structure:
            if inc_structure == "initial":
                props.append("initial_structure")
            else:
                props.append("structure")

        if not isinstance(chemsys_formula_id_criteria, dict):
            criteria = _MPResterLegacy.parse_criteria(chemsys_formula_id_criteria)
        else:
            criteria = chemsys_formula_id_criteria
        data = self.query(criteria, props)

        entries: list[ComputedEntry] = []
        for d in data:
            d["potcar_symbols"] = [
                f"{d['pseudo_potential']['functional']} {label}" for label in d["pseudo_potential"]["labels"]
            ]
            data = {"oxide_type": d["oxide_type"]}
            if property_data:
                data |= {k: d[k] for k in property_data}
            if not inc_structure:
                e = ComputedEntry(
                    d["unit_cell_formula"],
                    d["energy"],
                    parameters={k: d[k] for k in params},
                    data=data,
                    entry_id=d["task_id"],
                )

            else:
                prim = d["initial_structure"] if inc_structure == "initial" else d["structure"]
                if conventional_unit_cell:
                    struct = SpacegroupAnalyzer(prim).get_conventional_standard_structure()
                    energy = d["energy"] * (len(struct) / len(prim))
                else:
                    struct = prim.copy()
                    energy = d["energy"]
                e = ComputedStructureEntry(
                    struct,
                    energy,
                    parameters={k: d[k] for k in params},
                    data=data,
                    entry_id=d["task_id"],
                )
            entries.append(e)
        if compatible_only:
            # suppress the warning about missing oxidation states
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Failed to guess oxidation states.*")
                entries = MaterialsProject2020Compatibility().process_entries(entries, clean=True)
        if sort_by_e_above_hull:
            entries = sorted(entries, key=lambda entry: entry.data["e_above_hull"])
        return entries

    def get_pourbaix_entries(self, chemsys, solid_compat="MaterialsProject2020Compatibility"):
        """A helper function to get all entries necessary to generate
        a Pourbaix diagram from the rest interface.

        Args:
            chemsys (str | list[str]): Chemical system string comprising element
                symbols separated by dashes, e.g. "Li-Fe-O" or List of element
                symbols, e.g. ["Li", "Fe", "O"].
            solid_compat: Compatibility scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProject2020Compatibility) or an instance
                (e.g., MaterialsProject2020Compatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProject2020Compatibility
        """
        # imports are not top-level due to expense
        from pymatgen.analysis.phase_diagram import PhaseDiagram
        from pymatgen.analysis.pourbaix_diagram import IonEntry, PourbaixEntry
        from pymatgen.core.ion import Ion
        from pymatgen.entries.compatibility import (
            Compatibility,
            MaterialsProject2020Compatibility,
            MaterialsProjectAqueousCompatibility,
            MaterialsProjectCompatibility,
        )

        if solid_compat == "MaterialsProjectCompatibility":
            self.solid_compat = MaterialsProjectCompatibility()
        elif solid_compat == "MaterialsProject2020Compatibility":
            self.solid_compat = MaterialsProject2020Compatibility()
        elif isinstance(solid_compat, Compatibility):
            self.solid_compat = solid_compat
        else:
            raise ValueError(
                "Solid compatibility can only be 'MaterialsProjectCompatibility', "
                "'MaterialsProject2020Compatibility', or an instance of a Compatibility class"
            )

        pbx_entries = []

        if isinstance(chemsys, str):
            chemsys = chemsys.split("-")

        # Get ion entries first, because certain ions have reference
        # solids that aren't necessarily in the chemsys (Na2SO4)
        url = "/pourbaix_diagram/reference_data/" + "-".join(chemsys)
        ion_data = self._make_request(url)
        ion_ref_comps = [Composition(d["Reference Solid"]) for d in ion_data]
        ion_ref_elts = list(itertools.chain.from_iterable(i.elements for i in ion_ref_comps))
        ion_ref_entries = self.get_entries_in_chemsys(
            list(set([str(e) for e in ion_ref_elts] + ["O", "H"])),
            property_data=["e_above_hull"],
            compatible_only=False,
        )

        # suppress the warning about supplying the required energies; they will be calculated from the
        # entries we get from MPRester
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="You did not provide the required O2 and H2O energies.",
            )
            compat = MaterialsProjectAqueousCompatibility(solid_compat=self.solid_compat)
        # suppress the warning about missing oxidation states
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Failed to guess oxidation states.*")
            ion_ref_entries = compat.process_entries(ion_ref_entries)
        ion_ref_pd = PhaseDiagram(ion_ref_entries)

        # position the ion energies relative to most stable reference state
        for n, i_d in enumerate(ion_data):
            ion = Ion.from_formula(i_d["Name"])
            refs = [e for e in ion_ref_entries if e.reduced_formula == i_d["Reference Solid"]]
            if not refs:
                raise ValueError("Reference solid not contained in entry list")
            stable_ref = min(refs, key=lambda x: x.data["e_above_hull"])
            rf = stable_ref.composition.get_reduced_composition_and_factor()[1]

            solid_diff = ion_ref_pd.get_form_energy(stable_ref) - i_d["Reference solid energy"] * rf
            elt = i_d["Major_Elements"][0]
            correction_factor = ion.composition[elt] / stable_ref.composition[elt]
            energy = i_d["Energy"] + solid_diff * correction_factor
            ion_entry = IonEntry(ion, energy)
            pbx_entries.append(PourbaixEntry(ion_entry, f"ion-{n}"))

        # Construct the solid Pourbaix entries from filtered ion_ref entries

        extra_elts = set(ion_ref_elts) - {Element(s) for s in chemsys} - {Element("H"), Element("O")}
        for entry in ion_ref_entries:
            entry_elts = set(entry.elements)
            # Ensure no OH chemsys or extraneous elements from ion references
            if not (entry_elts <= {Element("H"), Element("O")} or extra_elts.intersection(entry_elts)):
                # Create new computed entry
                form_e = ion_ref_pd.get_form_energy(entry)
                new_entry = ComputedEntry(entry.composition, form_e, entry_id=entry.entry_id)
                pbx_entry = PourbaixEntry(new_entry)
                pbx_entries.append(pbx_entry)

        return pbx_entries

    def get_structure_by_material_id(
        self, material_id: str, final: bool = True, conventional_unit_cell: bool = False
    ) -> Structure:
        """Get a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project ID (e.g. mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.
            conventional_unit_cell (bool): Whether to get the standard conventional unit cell

        Returns:
            Structure object.
        """
        prop = "final_structure" if final else "initial_structure"
        data = self.get_data(material_id, prop=prop)
        if not data:
            try:
                new_material_id = self.get_materials_id_from_task_id(material_id)
                if new_material_id:
                    warnings.warn(
                        f"The calculation task {material_id} is mapped to canonical mp-id {new_material_id}, "
                        f"so structure for {new_material_id} returned. This is not an error, see "
                        f"documentation. If original task data for {material_id} is required, use "
                        "get_task_data(). To find the canonical mp-id from a task id use "
                        "get_materials_id_from_task_id().",
                        stacklevel=2,
                    )
                return self.get_structure_by_material_id(new_material_id)
            except MPRestError:
                raise MPRestError(
                    f"{material_id=} unknown, if this seems like an error "
                    "please let us know at matsci.org/materials-project"
                )

        structure = data[0][prop]
        if conventional_unit_cell:
            structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
        return structure

    def get_entry_by_material_id(
        self,
        material_id: str,
        compatible_only: bool = True,
        inc_structure: bool | Literal["initial"] | None = None,
        property_data: list[str] | None = None,
        conventional_unit_cell: bool = False,
    ) -> ComputedEntry | ComputedStructureEntry:
        """Get a ComputedEntry corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g. mp-1234).
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

        Raises:
            MPRestError if no data for given material_id is found.

        Returns:
            ComputedEntry or ComputedStructureEntry object.
        """
        data = self.get_entries(
            material_id,
            compatible_only=compatible_only,
            inc_structure=inc_structure,
            property_data=property_data,
            conventional_unit_cell=conventional_unit_cell,
        )
        if len(data) == 0:
            raise MPRestError(f"{material_id = } does not exist")
        return data[0]

    def get_dos_by_material_id(self, material_id):
        """Get a Dos corresponding to a material_id.

        REST Endpoint: https://materialsproject.org/rest/v2/materials/<mp-id>/vasp/dos

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g. mp-1234).

        Returns:
            A Dos object.
        """
        data = self.get_data(material_id, prop="dos")
        return data[0]["dos"]

    def get_bandstructure_by_material_id(self, material_id, line_mode=True):
        """Get a BandStructure corresponding to a material_id.

        REST Endpoint: https://materialsproject.org/rest/v2/materials/<mp-id>/vasp/bandstructure or
        https://materialsproject.org/rest/v2/materials/<mp-id>/vasp/bandstructure_uniform

        Args:
            material_id (str): Materials Project material_id.
            line_mode (bool): If True, fetch a BandStructureSymmLine object
                (default). If False, return the uniform band structure.

        Returns:
            BandStructure
        """
        prop = "bandstructure" if line_mode else "bandstructure_uniform"
        data = self.get_data(material_id, prop=prop)
        return data[0][prop]

    def get_phonon_dos_by_material_id(self, material_id: str) -> CompletePhononDos:
        """Get phonon density of states data corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            CompletePhononDos: A phonon DOS object.
        """
        return self._make_request(f"/materials/{material_id}/phonondos")

    def get_phonon_bandstructure_by_material_id(self, material_id: str) -> PhononBandStructureSymmLine:
        """Get phonon dispersion data corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            PhononBandStructureSymmLine: A phonon band structure.
        """
        return self._make_request(f"/materials/{material_id}/phononbs")

    def get_phonon_ddb_by_material_id(self, material_id: str) -> str:
        """Get ABINIT Derivative Data Base (DDB) output for phonon calculations.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            str: ABINIT DDB file as a string.
        """
        return self._make_request(f"/materials/{material_id}/abinit_ddb")

    def get_entries_in_chemsys(
        self,
        elements,
        compatible_only=True,
        inc_structure=None,
        property_data=None,
        conventional_unit_cell=False,
        additional_criteria=None,
    ):
        """Helper method to get a list of ComputedEntries in a chemical system.

        For example, elements = ["Li", "Fe", "O"] will return a list of all entries in the
        Li-Fe-O chemical system, i.e., all LixOy, FexOy, LixFey, LixFeyOz, Li, Fe and O
        phases. Extremely useful for creating phase diagrams of entire chemical systems.

        Args:
            elements (str | list[str]): Chemical system string comprising element
                symbols separated by dashes, e.g. "Li-Fe-O" or List of element
                symbols, e.g. ["Li", "Fe", "O"].
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
            additional_criteria (dict): Any additional criteria to pass. For instance, if you are only interested in
                stable entries, you can pass {"e_above_hull": {"$lte": 0.001}}.

        Returns:
            List of ComputedEntries.
        """
        if isinstance(elements, str):
            elements = elements.split("-")

        all_chemsyses = []
        for i in range(len(elements)):
            for els in itertools.combinations(elements, i + 1):
                all_chemsyses.append("-".join(sorted(els)))

        criteria = {"chemsys": {"$in": all_chemsyses}}
        if additional_criteria:
            criteria.update(additional_criteria)

        return self.get_entries(
            criteria,
            compatible_only=compatible_only,
            inc_structure=inc_structure,
            property_data=property_data,
            conventional_unit_cell=conventional_unit_cell,
        )

    def get_exp_thermo_data(self, formula):
        """Get a list of ThermoData objects associated with a formula using the
        Materials Project REST interface.

        Args:
            formula (str): A formula to search for.

        Returns:
            List of ThermoData objects.
        """
        return self.get_data(formula, data_type="exp")

    def get_exp_entry(self, formula):
        """Get an ExpEntry object, which is the experimental equivalent of a
        ComputedEntry and can be used for analyses using experimental data.

        Args:
            formula (str): A formula to search for.

        Returns:
            An ExpEntry object.
        """
        return ExpEntry(Composition(formula), self.get_exp_thermo_data(formula))

    def query(
        self,
        criteria,
        properties,
        chunk_size: int = 500,
        max_tries_per_chunk: int = 5,
        mp_decode: bool = True,
        show_progress_bar: bool = True,
    ):
        r"""Perform an advanced query using MongoDB-like syntax for directly
        querying the Materials Project database. This allows one to perform
        queries which are otherwise too cumbersome to perform using the standard
        convenience methods.

        Please consult the Materials API documentation at
        https://github.com/materialsproject/mapidoc, which provides a
        comprehensive explanation of the document schema used in the Materials
        Project (supported criteria and properties) and guidance on how best to
        query for the relevant information you need.

        For queries that request data on more than CHUNK_SIZE materials at once,
        this method will chunk a query by first retrieving a list of material
        IDs that satisfy CRITERIA, and then merging the criteria with a
        restriction to one chunk of materials at a time of size CHUNK_SIZE. You
        can opt out of this behavior by setting CHUNK_SIZE=0. To guard against
        intermittent server errors in the case of many chunks per query,
        possibly-transient server errors will result in re-trying a give chunk
        up to MAX_TRIES_PER_CHUNK times.

        Args:
            criteria (str/dict): Criteria of the query as a string or
                mongo-style dict.

                If string, it supports a powerful but simple string criteria.
                e.g. "Fe2O3" means search for materials with reduced_formula
                Fe2O3. Wild cards are also supported. e.g. "\\*2O" means get
                all materials whose formula can be formed as \\*2O, e.g.
                Li2O, K2O, etc.

                Other syntax examples:
                mp-1234: Interpreted as a Materials ID.
                Fe2O3 or *2O3: Interpreted as reduced formulas.
                Li-Fe-O or *-Fe-O: Interpreted as chemical systems.

                You can mix and match with spaces, which are interpreted as
                "OR". E.g. "mp-1234 FeO" means query for all compounds with
                reduced formula FeO or with materials_id mp-1234.

                Using a full dict syntax, even more powerful queries can be
                constructed. For example, {"elements":{"$in":["Li",
                "Na", "K"], "$all": ["O"]}, "nelements":2} selects all Li, Na
                and K oxides. {"band_gap": {"$gt": 1}} selects all materials
                with band gaps greater than 1 eV.
            properties (list): Properties to request for as a list. For
                example, ["formula", "formation_energy_per_atom"] returns
                the formula and formation energy per atom.
            chunk_size (int): Number of materials for which to fetch data at a
                time. More data-intensive properties may require smaller chunk
                sizes. Use chunk_size=0 to force no chunking -- this is useful
                when fetching only properties such as 'material_id'.
            max_tries_per_chunk (int): How many times to re-try fetching a given
                chunk when the server gives a 5xx error (e.g. a timeout error).
            mp_decode (bool): Whether to do a decoding to a Pymatgen object
                where possible. In some cases, it might be useful to just get
                the raw python dict, i.e., set to False.
            show_progress_bar (bool): Whether to show a progress bar for large queries.
                Defaults to True. Set to False to reduce visual noise.

        Returns:
            List of results. e.g.
            [{u'formula': {u'O': 1, u'Li': 2.0}},
            {u'formula': {u'Na': 2.0, u'O': 2.0}},
            {u'formula': {u'K': 1, u'O': 3.0}},
            ...]
        """
        if not isinstance(criteria, dict):
            criteria = self.parse_criteria(criteria)
        payload = {
            "criteria": json.dumps(criteria),
            "properties": json.dumps(properties),
        }
        if chunk_size == 0:
            return self._make_request("/query", payload=payload, method="POST", mp_decode=mp_decode)

        count_payload = payload.copy()
        count_payload["options"] = json.dumps({"count_only": True})
        num_results = self._make_request("/query", payload=count_payload, method="POST")
        if num_results <= chunk_size:
            return self._make_request("/query", payload=payload, method="POST", mp_decode=mp_decode)

        data = []
        mids = [dct["material_id"] for dct in self.query(criteria, ["material_id"], chunk_size=0)]
        chunks = get_chunks(mids, size=chunk_size)
        progress_bar = tqdm(total=len(mids), disable=not show_progress_bar)
        for chunk in chunks:
            chunk_criteria = criteria.copy()
            chunk_criteria |= {"material_id": {"$in": chunk}}
            n_tries = 0
            while n_tries < max_tries_per_chunk:
                try:
                    data += self.query(
                        chunk_criteria,
                        properties,
                        chunk_size=0,
                        mp_decode=mp_decode,
                    )
                    break
                except MPRestError as exc:
                    if match := re.search(r"error status code (\d+)", str(exc)):
                        if not match[1].startswith("5"):
                            raise
                        n_tries += 1
                        print(
                            "Unknown server error. Trying again in five "
                            f"seconds (will try at most {max_tries_per_chunk} times)..."
                        )
                        sleep(5)
            progress_bar.update(len(chunk))
        return data

    def submit_structures(
        self,
        structures,
        authors,
        projects=None,
        references="",
        remarks=None,
        data=None,
        histories=None,
        created_at=None,
    ) -> list[str]:
        """Submits a list of structures to the Materials Project as SNL files.
        The argument list mirrors the arguments for the StructureNL object,
        except that a list of structures with the same metadata is used as an
        input.

        Note:
            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for the future.

        Args:
            structures: A list of Structure objects
            authors (list): List of {"name":'', "email":''} dicts,
                *list* of Strings as 'John Doe <johndoe@gmail.com>',
                or a single String with commas separating authors
            projects ([str]): List of Strings ['Project A', 'Project B'].
                This applies to all structures.
            references (str): A String in BibTeX format. Again, this applies to
                all structures.
            remarks ([str]): List of Strings ['Remark A', 'Remark B']
            data ([dict]): A list of free form dict. Namespaced at the root
                level with an underscore, e.g. {"_materialsproject":<custom
                data>}. The length of data should be the same as the list of
                structures if not None.
            histories: List of list of dicts - [[{'name':'', 'url':'',
                'description':{}}], ...] The length of histories should be the
                same as the list of structures if not None.
            created_at (datetime): A datetime object

        Returns:
            list[str]: Inserted submission ids.
        """
        from pymatgen.util.provenance import StructureNL

        snl_list = StructureNL.from_structures(
            structures,
            authors,
            projects,
            references,
            remarks,
            data,
            histories,
            created_at,
        )
        return self.submit_snl(snl_list)

    def submit_snl(self, snl):
        """Submits a list of StructureNL to the Materials Project site.

        Note:
            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for the future.

        Args:
            snl (StructureNL/[StructureNL]): A single StructureNL, or a list
            of StructureNL objects

        Returns:
            list[str]: Inserted submission ids.

        Raises:
            MPRestError: If submission fails.
        """
        snl = snl if isinstance(snl, list) else [snl]
        json_data = [s.as_dict() for s in snl]
        payload = {"snl": json.dumps(json_data, cls=MontyEncoder)}
        response = self.session.post(f"{self.preamble}/snl/submit", data=payload)
        if response.status_code in [200, 400]:
            response = json.loads(response.text, cls=MontyDecoder)
            if response["valid_response"]:
                if response.get("warning"):
                    warnings.warn(response["warning"], stacklevel=2)
                return response["inserted_ids"]
            raise MPRestError(response["error"])

        raise MPRestError(f"REST error with status code {response.status_code} and error {response.text}")

    def delete_snl(self, snl_ids):
        """Delete earlier submitted SNLs.

        Note:
            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for the future.

        Args:
            snl_ids: List of SNL ids.

        Raises:
            MPRestError
        """
        payload = {"ids": json.dumps(snl_ids)}
        response = self.session.post(f"{self.preamble}/snl/delete", data=payload)

        if response.status_code in [200, 400]:
            response = json.loads(response.text, cls=MontyDecoder)
            if response["valid_response"]:
                if response.get("warning"):
                    warnings.warn(response["warning"], stacklevel=2)
                return response
            raise MPRestError(response["error"])

        raise MPRestError(f"REST error with status code {response.status_code} and error {response.text}")

    def query_snl(self, criteria):
        """Query for submitted SNLs.

        Note:
            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for the future.

        Args:
            criteria (dict): Query criteria.

        Returns:
            A dict, with a list of submitted SNLs in the "response" key.

        Raises:
            MPRestError
        """
        payload = {"criteria": json.dumps(criteria)}
        response = self.session.post(f"{self.preamble}/snl/query", data=payload)
        if response.status_code in [200, 400]:
            response = json.loads(response.text)
            if response["valid_response"]:
                if response.get("warning"):
                    warnings.warn(response["warning"], stacklevel=2)
                return response["response"]
            raise MPRestError(response["error"])

        raise MPRestError(f"REST error with status code {response.status_code} and error {response.text}")

    def submit_vasp_directory(
        self,
        rootdir,
        authors,
        projects=None,
        references="",
        remarks=None,
        master_data=None,
        master_history=None,
        created_at=None,
        ncpus=None,
    ):
        """Assimilate all vasp run directories beneath a particular
        directory using BorgQueen to obtain structures, and then submits thhem
        to the Materials Project as SNL files. VASP related meta data like
        initial structure and final energies are automatically incorporated.

        Note:
            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for the future.

        Args:
            rootdir (str): Rootdir to start assimilating VASP runs from.
            authors: *List* of {"name":'', "email":''} dicts,
                *list* of Strings as 'John Doe <johndoe@gmail.com>',
                or a single String with commas separating authors. The same
                list of authors should apply to all runs.
            projects ([str]): List of Strings ['Project A', 'Project B'].
                This applies to all structures.
            references (str): A String in BibTeX format. Again, this applies to
                all structures.
            remarks ([str]): List of Strings ['Remark A', 'Remark B']
            master_data (dict): A free form dict. Namespaced at the root
                level with an underscore, e.g. {"_materialsproject":<custom
                data>}. This data is added to all structures detected in the
                directory, in addition to other vasp data on a per structure
                basis.
            master_history: A master history to be added to all entries.
            created_at (datetime): A datetime object
            ncpus (int): Number of cpus to use in using BorgQueen to
                assimilate. Defaults to None, which means serial.
        """
        from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
        from pymatgen.apps.borg.queen import BorgQueen

        drone = VaspToComputedEntryDrone(inc_structure=True, data=["filename", "initial_structure"])
        queen = BorgQueen(drone, number_of_drones=ncpus)
        queen.parallel_assimilate(rootdir)

        structures = []
        metadata = []
        histories = []
        for e in queen.get_data():
            structures.append(e.structure)
            meta_dict = {
                "_vasp": {
                    "parameters": e.parameters,
                    "final_energy": e.energy,
                    "final_energy_per_atom": e.energy_per_atom,
                    "initial_structure": e.data["initial_structure"].as_dict(),
                }
            }
            if "history" in e.parameters:
                histories.append(e.parameters["history"])
            if master_data is not None:
                meta_dict.update(master_data)
            metadata.append(meta_dict)
        if master_history is not None:
            histories = master_history * len(structures)

        return self.submit_structures(
            structures,
            authors,
            projects=projects,
            references=references,
            remarks=remarks,
            data=metadata,
            histories=histories,
            created_at=created_at,
        )

    def get_stability(self, entries):
        """Get the stability of all entries."""
        payload = {"entries": json.dumps(entries, cls=MontyEncoder)}
        response = self.session.post(
            f"{self.preamble}/phase_diagram/calculate_stability",
            data=payload,
        )
        if response.status_code in [200, 400]:
            response = json.loads(response.text, cls=MontyDecoder)
            if response["valid_response"]:
                if response.get("warning"):
                    warnings.warn(response["warning"], stacklevel=2)
                return response["response"]
            raise MPRestError(response["error"])
        raise MPRestError(f"REST error with status code {response.status_code} and error {response.text}")

    def get_cohesive_energy(self, material_id, per_atom=False):
        """Get the cohesive for a material (eV per formula unit). Cohesive energy
            is defined as the difference between the bulk energy and the sum of
            total DFT energy of isolated atoms for atom elements in the bulk.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            per_atom (bool): Whether or not to return cohesive energy per atom

        Returns:
            Cohesive energy (eV).
        """
        entry = self.get_entry_by_material_id(material_id)
        e_bulk = entry.energy / entry.composition.get_integer_formula_and_factor()[1]
        comp_dict = entry.composition.reduced_composition.as_dict()

        isolated_atom_e_sum = 0
        for el in comp_dict:
            ent = self._make_request(f"/element/{el}/tasks/isolated_atom", mp_decode=False)[0]
            isolated_atom_e_sum += ent["output"]["final_energy_per_atom"] * comp_dict[el]
        e_coh_per_formula = isolated_atom_e_sum - e_bulk
        n_atoms = entry.composition.num_atoms
        return e_coh_per_formula / n_atoms if per_atom else e_coh_per_formula

    def get_reaction(self, reactants, products):
        """Get a reaction from the Materials Project.

        Args:
            reactants ([str]): List of formulas
            products ([str]): List of formulas

        Returns:
            rxn
        """
        return self._make_request(
            "/reaction",
            payload={"reactants[]": reactants, "products[]": products},
            mp_decode=False,
        )

    def get_substrates(self, material_id, number=50, orient=None):
        """Get a substrate list for a material id. The list is in order of
        increasing elastic energy if a elastic tensor is available for
        the material_id. Otherwise the list is in order of increasing
        matching area.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            orient (list) : substrate orientation to look for
            number (int) : number of substrates to return
                n=0 returns all available matches

        Returns:
            list of dicts with substrate matches
        """
        req = f"/materials/{material_id}/substrates?n={number}"
        if orient:
            req += f"&orient={' '.join(map(str, orient))}"
        return self._make_request(req)

    def get_all_substrates(self):
        """Get the list of all possible substrates considered in the
        Materials Project substrate database.

        Returns:
            list of material_ids corresponding to possible substrates
        """
        return self._make_request("/materials/all_substrate_ids")

    @due.dcite(
        Doi("10.1038/sdata.2016.80"),
        description="Data Descriptor: Surface energies of elemental crystals. Scientific Data",
    )
    def get_surface_data(self, material_id, miller_index=None, inc_structures=False):
        """Get surface data for a material. Useful for Wulff shapes.

        Reference for surface data:

        Tran, R., Xu, Z., Radhakrishnan, B., Winston, D., Sun, W., Persson, K.
        A., & Ong, S. P. (2016). Data Descriptor: Surface energies of elemental
        crystals. Scientific Data, 3(160080), 1-13.
        https://doi.org/10.1038/sdata.2016.80

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            miller_index (list of integer): The miller index of the surface.
            e.g. [3, 2, 1]. If miller_index is provided, only one dictionary
            of this specific plane will be returned.
            inc_structures (bool): Include final surface slab structures.
                These are unnecessary for Wulff shape construction.

        Returns:
            Surface data for material. Energies are given in SI units (J/m^2).
        """
        req = f"/materials/{material_id}/surfaces"
        if inc_structures:
            req += "?include_structures=true"

        if miller_index:
            surf_data_dict = self._make_request(req)
            surf_list = surf_data_dict["surfaces"]
            ucell = self.get_structure_by_material_id(material_id, conventional_unit_cell=True)
            eq_indices = get_symmetrically_equivalent_miller_indices(ucell, miller_index)
            for one_surf in surf_list:
                if tuple(one_surf["miller_index"]) in eq_indices:
                    return one_surf
            raise ValueError("Bad miller index.")
        return self._make_request(req)

    def get_wulff_shape(self, material_id):
        """Construct a Wulff shape for a material.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.

        Returns:
            pymatgen.analysis.wulff.WulffShape
        """
        from pymatgen.analysis.wulff import WulffShape
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        structure = self.get_structure_by_material_id(material_id)
        surfaces = self.get_surface_data(material_id)["surfaces"]
        lattice = SpacegroupAnalyzer(structure).get_conventional_standard_structure().lattice
        miller_energy_map = {}
        for surf in surfaces:
            miller = tuple(surf["miller_index"])
            # Prefer reconstructed surfaces, which have lower surface energies.
            if (miller not in miller_energy_map) or surf["is_reconstructed"]:
                miller_energy_map[miller] = surf["surface_energy"]
        millers, energies = zip(*miller_energy_map.items(), strict=True)
        return WulffShape(lattice, millers, energies)

    def get_gb_data(
        self,
        material_id=None,
        pretty_formula=None,
        chemsys=None,
        sigma=None,
        gb_plane=None,
        rotation_axis=None,
        include_work_of_separation=False,
    ):
        """Get grain boundary data for a material.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-129'.
            pretty_formula (str): The formula of metals. e.g. 'Fe'
            chemsys (str): The chemical system. e.g. 'Fe-O'
            sigma (int): The sigma value of a certain type of grain boundary
            gb_plane (list of integer): The Miller index of grain boundary plane. e.g. [1, 1, 1]
            rotation_axis (list of integer): The Miller index of rotation axis. e.g.
                [1, 0, 0], [1, 1, 0], and [1, 1, 1] Sigma value is determined by the combination of
                rotation axis and rotation angle. The five degrees of freedom (DOF) of one grain boundary
                include: rotation axis (2 DOFs), rotation angle (1 DOF), and grain boundary plane (2 DOFs).
            include_work_of_separation (bool): whether to include the work of separation
                (in unit of (J/m^2)). If you want to query the work of separation, please
                specify the material_id.

        Returns:
            A list of grain boundaries that satisfy the query conditions (sigma, gb_plane).
            Energies are given in SI units (J/m^2).
        """
        if gb_plane:
            gb_plane = ",".join(str(plane) for plane in gb_plane)
        if rotation_axis:
            rotation_axis = ",".join(str(ax) for ax in rotation_axis)

        payload = {
            "material_id": material_id,
            "pretty_formula": pretty_formula,
            "chemsys": chemsys,
            "sigma": sigma,
            "gb_plane": gb_plane,
            "rotation_axis": rotation_axis,
        }

        if include_work_of_separation and material_id:
            list_of_gbs = self._make_request("/grain_boundaries", payload=payload)
            for gb_dict in list_of_gbs:
                gb_energy = gb_dict["gb_energy"]
                gb_plane_int = gb_dict["gb_plane"]
                surface_energy = self.get_surface_data(material_id=material_id, miller_index=gb_plane_int)[
                    "surface_energy"
                ]
                work_of_sep = 2 * surface_energy - gb_energy  # calculate the work of separation
                gb_dict["work_of_separation"] = work_of_sep
            return list_of_gbs

        return self._make_request("/grain_boundaries", payload=payload)

    def get_interface_reactions(
        self,
        reactant1,
        reactant2,
        open_el=None,
        relative_mu=None,
        use_hull_energy=False,
    ):
        """Get critical reactions between two reactants.

        Get critical reactions ("kinks" in the mixing ratio where
        reaction products change) between two reactants. See the
        `pymatgen.analysis.interface_reactions` module for more info.

        Args:
            reactant1 (str): Chemical formula for reactant
            reactant2 (str): Chemical formula for reactant
            open_el (str): Element in reservoir available to system
            relative_mu (float): Relative chemical potential of element in
                reservoir with respect to pure substance. Must be non-positive.
            use_hull_energy (bool): Whether to use the convex hull energy for a
            given composition for the reaction energy calculation. If false,
            the energy of the ground state structure will be preferred; if a
            ground state can not be found for a composition, the convex hull
            energy will be used with a warning message.

        Returns:
            list: list of dicts of form {ratio,energy,rxn} where `ratio` is the
                reactant mixing ratio, `energy` is the reaction energy
                in eV/atom, and `rxn` is a
                `pymatgen.analysis.reaction_calculator.Reaction`.
        """
        payload = {
            "reactants": f"{reactant1} {reactant2}",
            "open_el": open_el,
            "relative_mu": relative_mu,
            "use_hull_energy": use_hull_energy,
        }
        return self._make_request("/interface_reactions", payload=payload, method="POST")

    def get_download_info(self, material_ids, task_types=None, file_patterns=None):
        """Get a list of URLs to retrieve raw VASP output files from the NoMaD repository.

        Args:
            material_ids (list): list of material identifiers (mp-id's)
            task_types (list): list of task types to include in download (see TaskType Enum class)
            file_patterns (list): list of wildcard file names to include for each task

        Returns:
            a tuple of 1) a dictionary mapping material_ids to task_ids and
            task_types, and 2) a list of URLs to download zip archives from
            NoMaD repository. Each zip archive will contain a manifest.json with
            metadata info, e.g. the task/external_ids that belong to a directory
        """
        # task_id's correspond to NoMaD external_id's
        task_types = [typ.value for typ in (task_types or []) if isinstance(typ, TaskType)]

        meta = {}
        for doc in self.query({"task_id": {"$in": material_ids}}, ["task_id", "blessed_tasks"]):
            for task_type, task_id in doc["blessed_tasks"].items():
                if task_types and task_type not in task_types:
                    continue
                mp_id = doc["task_id"]
                if meta.get(mp_id) is None:
                    meta[mp_id] = [{"task_id": task_id, "task_type": task_type}]
                else:
                    meta[mp_id].append({"task_id": task_id, "task_type": task_type})
        if not meta:
            raise ValueError(f"No tasks found for material id {material_ids}.")

        # return a list of URLs for NoMaD Downloads containing the list of files
        # for every external_id in `task_ids`
        # For reference, please visit https://nomad-lab.eu/prod/rae/api/

        # check if these task ids exist on NOMAD
        prefix = "https://nomad-lab.eu/prod/rae/api/repo/?"
        if file_patterns is not None:
            for file_pattern in file_patterns:
                prefix += f"{file_pattern=}&"
        prefix += "external_id="

        task_ids = [task["task_id"] for task_list in meta.values() for task in task_list]
        nomad_exist_task_ids = self._check_get_download_info_url_by_task_id(prefix=prefix, task_ids=task_ids)
        if len(nomad_exist_task_ids) != len(task_ids):
            self._print_help_message(nomad_exist_task_ids, task_ids, file_patterns, task_types)

        # generate download links for those that exist
        prefix = "https://nomad-lab.eu/prod/rae/api/raw/query?"
        if file_patterns is not None:
            for file_pattern in file_patterns:
                prefix += f"file_pattern={file_pattern}&"
        prefix += "external_id="

        urls = [f"{prefix}{task_ids}" for task_ids in nomad_exist_task_ids]
        return meta, urls

    @staticmethod
    def _print_help_message(nomad_exist_task_ids, task_ids, file_patterns, task_types):
        non_exist_ids = set(task_ids) - set(nomad_exist_task_ids)
        warnings.warn(
            f"For {file_patterns=}] and {task_types=}, \n"
            f"the following ids are not found on NOMAD [{list(non_exist_ids)}]. \n"
            f"If you need to upload them, please contact Patrick Huck at phuck@lbl.gov",
            stacklevel=2,
        )

    def _check_get_download_info_url_by_task_id(self, prefix, task_ids) -> list[str]:
        nomad_exist_task_ids: list[str] = []
        prefix = prefix.replace("/raw/query", "/repo/")
        for task_id in task_ids:
            url = prefix + task_id
            if self._check_nomad_exist(url):
                nomad_exist_task_ids.append(task_id)
        return nomad_exist_task_ids

    @staticmethod
    def _check_nomad_exist(url) -> bool:
        response = requests.get(url=url, timeout=60)
        if response.status_code != 200:
            return False
        content = json.loads(response.text)
        return content["pagination"]["total"] != 0

    @staticmethod
    def parse_criteria(criteria_string) -> dict:
        """Parse a powerful and simple string criteria and generates a proper
        mongo syntax criteria.

        Args:
            criteria_string (str): A string representing a search criteria.
                Also supports wild cards. e.g.
                something like "*2O" gets converted to
                {'pretty_formula': {'$in': [u'B2O', u'Xe2O', u"Li2O", ...]}}

                Other syntax examples:
                    mp-1234: Interpreted as a Materials ID.
                    Fe2O3 or *2O3: Interpreted as reduced formulas.
                    Li-Fe-O or *-Fe-O: Interpreted as chemical systems.

                You can mix and match with spaces, which are interpreted as
                "OR". e.g. "mp-1234 FeO" means query for all compounds with
                reduced formula FeO or with materials_id mp-1234.

        Returns:
            A mongo query dict.
        """
        tokens = criteria_string.split()

        def parse_sym(sym):
            if sym == "*":
                return [el.symbol for el in Element]

            if match := re.match(r"\{(.*)\}", sym):
                return [s.strip() for s in match[1].split(",")]
            return [sym]

        def parse_tok(t):
            if re.match(r"\w+-\d+", t):
                return {"task_id": t}
            if "-" in t:
                elements = [parse_sym(sym) for sym in t.split("-")]
                chem_sys_lst = []
                for cs in itertools.product(*elements):
                    if len(set(cs)) == len(cs):
                        # Check for valid symbols
                        cs = [Element(s).symbol for s in cs]
                        chem_sys_lst.append("-".join(sorted(cs)))
                return {"chemsys": {"$in": chem_sys_lst}}
            all_formulas = set()
            explicit_els = []
            wild_card_els = []
            for sym in re.findall(r"(\*[\.\d]*|\{.*\}[\.\d]*|[A-Z][a-z]*)[\.\d]*", t):
                if ("*" in sym) or ("{" in sym):
                    wild_card_els.append(sym)
                else:
                    match = re.match(r"([A-Z][a-z]*)[\.\d]*", sym)
                    explicit_els.append(match[1])
            n_elements = len(wild_card_els) + len(set(explicit_els))
            parts = re.split(r"(\*|\{.*\})", t)
            parts = [parse_sym(s) for s in parts if s != ""]
            for formula in itertools.product(*parts):
                comp = Composition("".join(formula))
                if len(comp) == n_elements:
                    # Check for valid Elements in keys.
                    for elem in comp:
                        Element(elem.symbol)
                    all_formulas.add(comp.reduced_formula)
            return {"pretty_formula": {"$in": list(all_formulas)}}

        if len(tokens) == 1:
            return parse_tok(tokens[0])
        return {"$or": list(map(parse_tok, tokens))}


class MPRestError(Exception):
    """Exception class for legacy MPRestAdaptor. Raised when query is malformed."""


def get_chunks(sequence: Sequence[Any], size=1):
    """
    Args:
        sequence (Sequence[Any]): Any sequence.
        size (int): Chunk length. Defaults to 1.

    Returns:
        list[Sequence[Any]]: input sequence in chunks of length size.
    """
    chunks = math.ceil(len(sequence) / float(size))
    return [sequence[i * size : (i + 1) * size] for i in range(chunks)]
