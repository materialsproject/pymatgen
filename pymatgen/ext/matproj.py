# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
https://www.materialsproject.org/dashboard.
"""

import sys
import itertools
import json
import platform
import re
import warnings
from time import sleep
from enum import Enum, unique
from collections import defaultdict

import requests
import ruamel.yaml as yaml
from monty.json import MontyDecoder, MontyEncoder
from monty.serialization import dumpfn

from pymatgen import SETTINGS, SETTINGS_FILE, __version__ as pmg_version

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.core.surface import get_symmetrically_equivalent_miller_indices
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.entries.compatibility import (
    MaterialsProjectCompatibility,
    MaterialsProjectAqueousCompatibility,
)
from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.sequence import get_chunks, PBar


@unique
class TaskType(Enum):
    """task types available in MP"""

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


class MPRester:
    """
    A class to conveniently interface with the Materials Project REST
    interface. The recommended way to use MPRester is with the "with" context
    manager to ensure that sessions are properly closed after usage::

        with MPRester("API_KEY") as m:
            do_something

    MPRester uses the "requests" package, which provides for HTTP connection
    pooling. All connections are made via https for security.

    For more advanced uses of the Materials API, please consult the API
    documentation at https://github.com/materialsproject/mapidoc.
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
        api_key=None,
        endpoint=None,
        notify_db_version=True,
        include_user_agent=True,
    ):
        """
        Args:
            api_key (str): A String API key for accessing the MaterialsProject
                REST interface. Please obtain your API key at
                https://www.materialsproject.org/dashboard. If this is None,
                the code will check if there is a "PMG_MAPI_KEY" setting.
                If so, it will use that environment variable. This makes
                easier for heavy users to simply add this environment variable to
                their setups and MPRester can then be called without any arguments.
            endpoint (str): Url of endpoint to access the MaterialsProject REST
                interface. Defaults to the standard Materials Project REST
                address at "https://materialsproject.org/rest/v2", but
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
        if api_key is not None:
            self.api_key = api_key
        else:
            self.api_key = SETTINGS.get("PMG_MAPI_KEY", "")
        if endpoint is not None:
            self.preamble = endpoint
        else:
            self.preamble = SETTINGS.get(
                "PMG_MAPI_ENDPOINT", "https://materialsproject.org/rest/v2"
            )

        if self.preamble != "https://materialsproject.org/rest/v2":
            warnings.warn("Non-default endpoint used: {}".format(self.preamble))

        self.session = requests.Session()
        self.session.headers = {"x-api-key": self.api_key}
        if include_user_agent:
            pymatgen_info = "pymatgen/" + pmg_version
            python_info = "Python/{}.{}.{}".format(
                sys.version_info.major, sys.version_info.minor, sys.version_info.micro
            )
            platform_info = "{}/{}".format(platform.system(), platform.release())
            self.session.headers["user-agent"] = "{} ({} {})".format(
                pymatgen_info, python_info, platform_info
            )

        if notify_db_version:
            db_version = self.get_database_version()
            print(f"Connection established to Materials Project database, version {db_version}.")

            try:
                with open(SETTINGS_FILE, "rt") as f:
                    d = yaml.safe_load(f)
            except IOError:
                d = {}

            if "MAPI_DB_VERSION" not in d:
                d["MAPI_DB_VERSION"] = {"LOG": {}, "LAST_ACCESSED": None}

            # store a log of what database versions are being connected to
            if db_version not in d["MAPI_DB_VERSION"]["LOG"]:
                d["MAPI_DB_VERSION"]["LOG"][db_version] = 1
            else:
                d["MAPI_DB_VERSION"]["LOG"][db_version] += 1

            # alert user if db version changed
            last_accessed = d["MAPI_DB_VERSION"]["LAST_ACCESSED"]
            if last_accessed and last_accessed != db_version:
                print(
                    f"This database version has changed from the database last accessed ({last_accessed}).\n"
                    f"Please see release notes on materialsproject.org for information about what has changed."
                )
            d["MAPI_DB_VERSION"]["LAST_ACCESSED"] = db_version

            # write out new database log
            dumpfn(d, SETTINGS_FILE)

    def __enter__(self):
        """
        Support for "with" context.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Support for "with" context.
        """
        self.session.close()

    def _make_request(self, sub_url, payload=None, method="GET", mp_decode=True):
        response = None
        url = self.preamble + sub_url
        try:
            if method == "POST":
                response = self.session.post(url, data=payload, verify=True)
            else:
                response = self.session.get(url, params=payload, verify=True)
            if response.status_code in [200, 400]:
                if mp_decode:
                    data = json.loads(response.text, cls=MontyDecoder)
                else:
                    data = json.loads(response.text)
                if data["valid_response"]:
                    if data.get("warning"):
                        warnings.warn(data["warning"])
                    return data["response"]
                raise MPRestError(data["error"])

            raise MPRestError(
                "REST query returned with error status code {}".format(
                    response.status_code
                )
            )

        except Exception as ex:
            msg = (
                "{}. Content: {}".format(str(ex), response.content)
                if hasattr(response, "content")
                else str(ex)
            )
            raise MPRestError(msg)

    def get_database_version(self):
        """
        The Materials Project database is periodically updated and has a
        database version associated with it. When the database is updated,
        consolidated data (information about "a material") may and does
        change, while calculation data about a specific calculation task
        remains unchanged and available for querying via its task_id.

        The database version is set as a date in the format YYYY-MM-DD,
        where "-DD" may be optional. An additional numerical suffix
        might be added if multiple releases happen on the same day.

        Returns: database version as a string
        """
        d = self._make_request("/api_check")
        return d["version"]["db"]

    def get_materials_id_from_task_id(self, task_id):
        """
        Returns a new MP materials id from a task id (which can be
        equivalent to an old materials id)

        Args:
            task_id (str): A task id.

        Returns:
            materials_id (str)
        """
        return self._make_request("/materials/mid_from_tid/%s" % task_id)

    def get_materials_id_references(self, material_id):
        """
        Returns all references for a materials id.

        Args:
            material_id (str): A material id.

        Returns:
            BibTeX (str)
        """
        return self._make_request("/materials/%s/refs" % material_id)

    def get_data(self, chemsys_formula_id, data_type="vasp", prop=""):
        """
        Flexible method to get any data using the Materials Project REST
        interface. Generally used by other methods for more specific queries.

        Format of REST return is *always* a list of dict (regardless of the
        number of pieces of data returned. The general format is as follows:

        [{"material_id": material_id, "property_name" : value}, ...]

        This is generally a call to
        https://www.materialsproject.org/rest/v2/materials/vasp/<prop>.
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
        sub_url = "/materials/%s/%s" % (chemsys_formula_id, data_type)
        if prop:
            sub_url += "/" + prop
        return self._make_request(sub_url)

    def get_materials_ids(self, chemsys_formula):
        """
        Get all materials ids for a formula or chemsys.

        Args:
            chemsys_formula (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3).

        Returns:
            ([str]) List of all materials ids.
        """
        return self._make_request(
            "/materials/%s/mids" % chemsys_formula, mp_decode=False
        )

    def get_doc(self, materials_id):
        """
        Get the entire data document for one materials id. Use this judiciously.

        REST Endpoint: https://www.materialsproject.org/materials/<mp-id>/doc.

        Args:
            materials_id (str): E.g., mp-1143 for Al2O3

        Returns:
            Dict of json document of all data that is displayed on a materials
            details page.
        """
        return self._make_request("/materials/%s/doc" % materials_id, mp_decode=False)

    def get_xas_data(self, material_id, absorbing_element):
        """
        Get X-ray absorption spectroscopy data for absorbing element in the
        structure corresponding to a material_id. Only X-ray Absorption Near Edge
        Structure (XANES) for K-edge is supported.

        REST Endpoint:
        https://www.materialsproject.org/materials/<mp-id>/xas/<absorbing_element>.

        Args:
            material_id (str): E.g., mp-1143 for Al2O3
            absorbing_element (str): The absorbing element in the corresponding
                structure. E.g., Al in Al2O3
        """
        element_list = self.get_data(material_id, prop="elements")[0]["elements"]
        if absorbing_element not in element_list:
            raise ValueError(
                "{} element not contained in corresponding structure with "
                "mp_id: {}".format(absorbing_element, material_id)
            )
        data = self._make_request(
            "/materials/{}/xas/{}".format(material_id, absorbing_element),
            mp_decode=False,
        )
        return data[0]

    def get_task_data(self, chemsys_formula_id, prop=""):
        """
        Flexible method to get any data using the Materials Project REST
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
        sub_url = "/tasks/%s" % chemsys_formula_id
        if prop:
            sub_url += "/" + prop
        return self._make_request(sub_url)

    def get_structures(self, chemsys_formula_id, final=True):
        """
        Get a list of Structures corresponding to a chemical system, formula,
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
        """
        Finds matching structures on the Materials Project site.

        Args:
            filename_or_structure: filename or Structure object

        Returns:
            A list of matching materials project ids for structure.

        Raises:
            MPRestError
        """
        try:
            if isinstance(filename_or_structure, str):
                s = Structure.from_file(filename_or_structure)
            elif isinstance(filename_or_structure, Structure):
                s = filename_or_structure
            else:
                raise MPRestError("Provide filename or Structure object.")
            payload = {"structure": json.dumps(s.as_dict(), cls=MontyEncoder)}
            response = self.session.post(
                "{}/find_structure".format(self.preamble), data=payload
            )
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    return resp["response"]
                raise MPRestError(resp["error"])
            raise MPRestError(
                "REST error with status code {} and error {}".format(
                    response.status_code, response.text
                )
            )
        except Exception as ex:
            raise MPRestError(str(ex))

    def get_entries(
        self,
        chemsys_formula_id_criteria,
        compatible_only=True,
        inc_structure=None,
        property_data=None,
        conventional_unit_cell=False,
        sort_by_e_above_hull=False,
    ):
        """
        Get a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id or full criteria.

        Args:
            chemsys_formula_id_criteria (str/dict): A chemical system
                (e.g., Li-Fe-O), or formula (e.g., Fe2O3) or materials_id
                (e.g., mp-1234) or full Mongo-style dict criteria.
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProjectCompatibility class,
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
        props = ["energy", "unit_cell_formula", "task_id"] + params
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
            criteria = MPRester.parse_criteria(chemsys_formula_id_criteria)
        else:
            criteria = chemsys_formula_id_criteria
        data = self.query(criteria, props)

        entries = []
        for d in data:
            d["potcar_symbols"] = [
                "%s %s" % (d["pseudo_potential"]["functional"], l)
                for l in d["pseudo_potential"]["labels"]
            ]
            data = {"oxide_type": d["oxide_type"]}
            if property_data:
                data.update({k: d[k] for k in property_data})
            if not inc_structure:
                e = ComputedEntry(
                    d["unit_cell_formula"],
                    d["energy"],
                    parameters={k: d[k] for k in params},
                    data=data,
                    entry_id=d["task_id"],
                )

            else:
                prim = (
                    d["initial_structure"]
                    if inc_structure == "initial"
                    else d["structure"]
                )
                if conventional_unit_cell:
                    s = SpacegroupAnalyzer(prim).get_conventional_standard_structure()
                    energy = d["energy"] * (len(s) / len(prim))
                else:
                    s = prim.copy()
                    energy = d["energy"]
                e = ComputedStructureEntry(
                    s,
                    energy,
                    parameters={k: d[k] for k in params},
                    data=data,
                    entry_id=d["task_id"],
                )
            entries.append(e)
        if compatible_only:
            from pymatgen.entries.compatibility import MaterialsProjectCompatibility

            entries = MaterialsProjectCompatibility().process_entries(entries)
        if sort_by_e_above_hull:
            entries = sorted(entries, key=lambda entry: entry.data["e_above_hull"])
        return entries

    def get_pourbaix_entries(
        self, chemsys, solid_compat=MaterialsProjectCompatibility
    ):
        """
        A helper function to get all entries necessary to generate
        a pourbaix diagram from the rest interface.

        Args:
            chemsys (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].
            solid_compat: Compatiblity scheme used to pre-process solid DFT energies prior to applying aqueous
                energy adjustments. May be passed as a class (e.g. MaterialsProjectCompatibility) or an instance
                (e.g., MaterialsProjectCompatibility()). If None, solid DFT energies are used as-is.
                Default: MaterialsProjectCompatibility
        """
        from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry
        from pymatgen.analysis.phase_diagram import PhaseDiagram
        from pymatgen.core.ion import Ion

        pbx_entries = []

        if isinstance(chemsys, str):
            chemsys = chemsys.split('-')

        # Get ion entries first, because certain ions have reference
        # solids that aren't necessarily in the chemsys (Na2SO4)
        url = "/pourbaix_diagram/reference_data/" + "-".join(chemsys)
        ion_data = self._make_request(url)
        ion_ref_comps = [Composition(d["Reference Solid"]) for d in ion_data]
        ion_ref_elts = list(
            itertools.chain.from_iterable(i.elements for i in ion_ref_comps)
        )
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
            compat = MaterialsProjectAqueousCompatibility(solid_compat=solid_compat)
        ion_ref_entries = compat.process_entries(ion_ref_entries)
        ion_ref_pd = PhaseDiagram(ion_ref_entries)

        # position the ion energies relative to most stable reference state
        for n, i_d in enumerate(ion_data):
            ion = Ion.from_formula(i_d["Name"])
            refs = [
                e
                for e in ion_ref_entries
                if e.composition.reduced_formula == i_d["Reference Solid"]
            ]
            if not refs:
                raise ValueError("Reference solid not contained in entry list")
            stable_ref = sorted(refs, key=lambda x: x.data["e_above_hull"])[0]
            rf = stable_ref.composition.get_reduced_composition_and_factor()[1]
            solid_diff = (
                ion_ref_pd.get_form_energy(stable_ref)
                - i_d["Reference solid energy"] * rf
            )
            elt = i_d["Major_Elements"][0]
            correction_factor = ion.composition[elt] / stable_ref.composition[elt]
            energy = i_d["Energy"] + solid_diff * correction_factor
            ion_entry = IonEntry(ion, energy)
            pbx_entries.append(PourbaixEntry(ion_entry, "ion-{}".format(n)))

        # Construct the solid pourbaix entries from filtered ion_ref entries
        extra_elts = (
            set(ion_ref_elts)
            - {Element(s) for s in chemsys}
            - {Element("H"), Element("O")}
        )
        for entry in ion_ref_entries:
            entry_elts = set(entry.composition.elements)
            # Ensure no OH chemsys or extraneous elements from ion references
            if not (
                entry_elts <= {Element("H"), Element("O")}
                or extra_elts.intersection(entry_elts)
            ):
                # Create new computed entry
                form_e = ion_ref_pd.get_form_energy(entry)
                new_entry = ComputedEntry(
                    entry.composition, form_e, entry_id=entry.entry_id
                )
                pbx_entry = PourbaixEntry(new_entry)
                pbx_entries.append(pbx_entry)

        return pbx_entries

    def get_structure_by_material_id(
        self, material_id, final=True, conventional_unit_cell=False
    ):
        """
        Get a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell

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
                        "The calculation task {} is mapped to canonical mp-id {}, "
                        "so structure for {} returned. "
                        "This is not an error, see documentation. "
                        "If original task data for {} is required, "
                        "use get_task_data(). To find the canonical mp-id from a task id "
                        "use get_materials_id_from_task_id().".format(
                            material_id, new_material_id, new_material_id, material_id
                        )
                    )
                return self.get_structure_by_material_id(new_material_id)
            except MPRestError:
                raise MPRestError(
                    "material_id {} unknown, if this seems like "
                    "an error please let us know at "
                    "matsci.org/materials-project".format(material_id)
                )
        if conventional_unit_cell:
            data[0][prop] = SpacegroupAnalyzer(
                data[0][prop]
            ).get_conventional_standard_structure()
        return data[0][prop]

    def get_entry_by_material_id(
        self,
        material_id,
        compatible_only=True,
        inc_structure=None,
        property_data=None,
        conventional_unit_cell=False,
    ):
        """
        Get a ComputedEntry corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProjectCompatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            inc_structure (str): If None, entries returned are
                ComputedEntries. If inc_structure="final",
                ComputedStructureEntries with final structures are returned.
                Otherwise, ComputedStructureEntries with initial structures
                are returned.
            property_data (list): Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell

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
        return data[0]

    def get_dos_by_material_id(self, material_id):
        """
        Get a Dos corresponding to a material_id.

        REST Endpoint: https://www.materialsproject.org/rest/v2/materials/<mp-id>/vasp/dos

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).

        Returns:
            A Dos object.
        """
        data = self.get_data(material_id, prop="dos")
        return data[0]["dos"]

    def get_bandstructure_by_material_id(self, material_id, line_mode=True):
        """
        Get a BandStructure corresponding to a material_id.

        REST Endpoint: https://www.materialsproject.org/rest/v2/materials/<mp-id>/vasp/bandstructure or
        https://www.materialsproject.org/rest/v2/materials/<mp-id>/vasp/bandstructure_uniform

        Args:
            material_id (str): Materials Project material_id.
            line_mode (bool): If True, fetch a BandStructureSymmLine object
                (default). If False, return the uniform band structure.

        Returns:
            A BandStructure object.
        """
        prop = "bandstructure" if line_mode else "bandstructure_uniform"
        data = self.get_data(material_id, prop=prop)
        return data[0][prop]

    def get_phonon_dos_by_material_id(self, material_id):
        """
        Get phonon density of states data corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            ﻿CompletePhononDos: A phonon DOS object.
        """
        return self._make_request("/materials/{}/phonondos".format(material_id))

    def get_phonon_bandstructure_by_material_id(self, material_id):
        """
        Get phonon dispersion data corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            PhononBandStructureSymmLine: A phonon band structure.
        """
        return self._make_request("/materials/{}/phononbs".format(material_id))

    def get_phonon_ddb_by_material_id(self, material_id):
        """
        Get ABINIT Derivative Data Base (DDB) output for phonon calculations.

        Args:
            material_id (str): Materials Project material_id.

        Returns:
            str: ABINIT DDB file as a string.
        """
        return self._make_request("/materials/{}/abinit_ddb".format(material_id))

    def get_entries_in_chemsys(
        self,
        elements,
        compatible_only=True,
        inc_structure=None,
        property_data=None,
        conventional_unit_cell=False,
    ):
        """
        Helper method to get a list of ComputedEntries in a chemical system.
        For example, elements = ["Li", "Fe", "O"] will return a list of all
        entries in the Li-Fe-O chemical system, i.e., all LixOy,
        FexOy, LixFey, LixFeyOz, Li, Fe and O phases. Extremely useful for
        creating phase diagrams of entire chemical systems.

        Args:
            elements (str or [str]): Chemical system string comprising element
                symbols separated by dashes, e.g., "Li-Fe-O" or List of element
                symbols, e.g., ["Li", "Fe", "O"].
            compatible_only (bool): Whether to return only "compatible"
                entries. Compatible entries are entries that have been
                processed using the MaterialsProjectCompatibility class,
                which performs adjustments to allow mixing of GGA and GGA+U
                calculations for more accurate phase diagrams and reaction
                energies.
            inc_structure (str): If None, entries returned are
                ComputedEntries. If inc_structure="final",
                ComputedStructureEntries with final structures are returned.
                Otherwise, ComputedStructureEntries with initial structures
                are returned.
            property_data (list): Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.
            conventional_unit_cell (bool): Whether to get the standard
                conventional unit cell

        Returns:
            List of ComputedEntries.

        """
        if isinstance(elements, str):
            elements = elements.split("-")

        all_chemsyses = []
        for i in range(len(elements)):
            for els in itertools.combinations(elements, i + 1):
                all_chemsyses.append("-".join(sorted(els)))

        entries = self.get_entries(
            {"chemsys": {"$in": all_chemsyses}},
            compatible_only=compatible_only,
            inc_structure=inc_structure,
            property_data=property_data,
            conventional_unit_cell=conventional_unit_cell,
        )
        return entries

    def get_exp_thermo_data(self, formula):
        """
        Get a list of ThermoData objects associated with a formula using the
        Materials Project REST interface.

        Args:
            formula (str): A formula to search for.

        Returns:
            List of ThermoData objects.
        """
        return self.get_data(formula, data_type="exp")

    def get_exp_entry(self, formula):
        """
        Returns an ExpEntry object, which is the experimental equivalent of a
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
        chunk_size=500,
        max_tries_per_chunk=5,
        mp_decode=True,
    ):
        r"""

        Performs an advanced query using MongoDB-like syntax for directly
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
                E.g., "Fe2O3" means search for materials with reduced_formula
                Fe2O3. Wild cards are also supported. E.g., "\\*2O" means get
                all materials whose formula can be formed as \\*2O, e.g.,
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

        Returns:
            List of results. E.g.,
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
            return self._make_request(
                "/query", payload=payload, method="POST", mp_decode=mp_decode
            )

        count_payload = payload.copy()
        count_payload["options"] = json.dumps({"count_only": True})
        num_results = self._make_request("/query", payload=count_payload, method="POST")
        if num_results <= chunk_size:
            return self._make_request(
                "/query", payload=payload, method="POST", mp_decode=mp_decode
            )

        data = []
        mids = [
            d["material_id"]
            for d in self.query(criteria, ["material_id"], chunk_size=0)
        ]
        chunks = get_chunks(mids, size=chunk_size)
        progress_bar = PBar(total=len(mids))
        for chunk in chunks:
            chunk_criteria = criteria.copy()
            chunk_criteria.update({"material_id": {"$in": chunk}})
            num_tries = 0
            while num_tries < max_tries_per_chunk:
                try:
                    data.extend(
                        self.query(
                            chunk_criteria,
                            properties,
                            chunk_size=0,
                            mp_decode=mp_decode,
                        )
                    )
                    break
                except MPRestError as e:
                    # pylint: disable=E1101
                    match = re.search(r"error status code (\d+)", e.message)
                    if match:
                        if not match.group(1).startswith("5"):
                            raise e
                        num_tries += 1
                        print(
                            "Unknown server error. Trying again in five "
                            "seconds (will try at most {} times)...".format(
                                max_tries_per_chunk
                            )
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
    ):
        """
        Submits a list of structures to the Materials Project as SNL files.
        The argument list mirrors the arguments for the StructureNL object,
        except that a list of structures with the same metadata is used as an
        input.

        .. note::

            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for
            the future.

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
            A list of inserted submission ids.
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
        self.submit_snl(snl_list)

    def submit_snl(self, snl):
        """
        Submits a list of StructureNL to the Materials Project site.

        .. note::

            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for
            the future.

        Args:
            snl (StructureNL/[StructureNL]): A single StructureNL, or a list
            of StructureNL objects

        Returns:
            A list of inserted submission ids.

        Raises:
            MPRestError
        """
        try:
            snl = snl if isinstance(snl, list) else [snl]
            jsondata = [s.as_dict() for s in snl]
            payload = {"snl": json.dumps(jsondata, cls=MontyEncoder)}
            response = self.session.post(
                "{}/snl/submit".format(self.preamble), data=payload
            )
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp["inserted_ids"]
                raise MPRestError(resp["error"])

            raise MPRestError(
                "REST error with status code {} and error {}".format(
                    response.status_code, response.text
                )
            )

        except Exception as ex:
            raise MPRestError(str(ex))

    def delete_snl(self, snl_ids):
        """
        Delete earlier submitted SNLs.

        .. note::

            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for
            the future.

        Args:
            snl_ids: List of SNL ids.

        Raises:
            MPRestError
        """
        try:
            payload = {"ids": json.dumps(snl_ids)}
            response = self.session.post(
                "{}/snl/delete".format(self.preamble), data=payload
            )

            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp
                raise MPRestError(resp["error"])

            raise MPRestError(
                "REST error with status code {} and error {}".format(
                    response.status_code, response.text
                )
            )

        except Exception as ex:
            raise MPRestError(str(ex))

    def query_snl(self, criteria):
        """
        Query for submitted SNLs.

        .. note::

            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for
            the future.

        Args:
            criteria (dict): Query criteria.

        Returns:
            A dict, with a list of submitted SNLs in the "response" key.

        Raises:
            MPRestError
        """
        try:
            payload = {"criteria": json.dumps(criteria)}
            response = self.session.post(
                "{}/snl/query".format(self.preamble), data=payload
            )
            if response.status_code in [200, 400]:
                resp = json.loads(response.text)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp["response"]
                raise MPRestError(resp["error"])

            raise MPRestError(
                "REST error with status code {} and error {}".format(
                    response.status_code, response.text
                )
            )

        except Exception as ex:
            raise MPRestError(str(ex))

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
        """
        Assimilates all vasp run directories beneath a particular
        directory using BorgQueen to obtain structures, and then submits thhem
        to the Materials Project as SNL files. VASP related meta data like
        initial structure and final energies are automatically incorporated.

        .. note::

            As of now, this MP REST feature is open only to a select group of
            users. Opening up submissions to all users is being planned for
            the future.

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

        drone = VaspToComputedEntryDrone(
            inc_structure=True, data=["filename", "initial_structure"]
        )
        queen = BorgQueen(drone, number_of_drones=ncpus)
        queen.parallel_assimilate(rootdir)

        structures = []
        metadata = []
        histories = []
        for e in queen.get_data():
            structures.append(e.structure)
            m = {
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
                m.update(master_data)
            metadata.append(m)
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
        """
        Returns the stability of all entries.
        """
        try:
            payload = {"entries": json.dumps(entries, cls=MontyEncoder)}
            response = self.session.post(
                "{}/phase_diagram/calculate_stability".format(self.preamble),
                data=payload,
            )
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp["response"]
                raise MPRestError(resp["error"])
            raise MPRestError(
                "REST error with status code {} and error {}".format(
                    response.status_code, response.text
                )
            )
        except Exception as ex:
            raise MPRestError(str(ex))

    def get_cohesive_energy(self, material_id, per_atom=False):
        """
        Gets the cohesive for a material (eV per formula unit). Cohesive energy
            is defined as the difference between the bulk energy and the sum of
            total DFT energy of isolated atoms for atom elements in the bulk.
        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            per_atom (bool): Whether or not to return cohesive energy per atom
        Returns:
            Cohesive energy (eV).
        """
        entry = self.get_entry_by_material_id(material_id)
        ebulk = entry.energy / entry.composition.get_integer_formula_and_factor()[1]
        comp_dict = entry.composition.reduced_composition.as_dict()

        isolated_atom_e_sum, n = 0, 0
        for el in comp_dict.keys():
            e = self._make_request(
                "/element/%s/tasks/isolated_atom" % (el), mp_decode=False
            )[0]
            isolated_atom_e_sum += e["output"]["final_energy_per_atom"] * comp_dict[el]
            n += comp_dict[el]
        ecoh_per_formula = isolated_atom_e_sum - ebulk
        return ecoh_per_formula / n if per_atom else ecoh_per_formula

    def get_reaction(self, reactants, products):
        """
        Gets a reaction from the Materials Project.

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
        """
        Get a substrate list for a material id. The list is in order of
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
        req = "/materials/{}/substrates?n={}".format(material_id, number)
        if orient:
            req += "&orient={}".format(" ".join(map(str, orient)))
        return self._make_request(req)

    def get_all_substrates(self):
        """
        Gets the list of all possible substrates considered in the
        Materials Project substrate database

        Returns:
            list of material_ids corresponding to possible substrates
        """

        return self._make_request("/materials/all_substrate_ids")

    def get_surface_data(self, material_id, miller_index=None, inc_structures=False):
        """
        Gets surface data for a material. Useful for Wulff shapes.

        Reference for surface data:

        Tran, R., Xu, Z., Radhakrishnan, B., Winston, D., Sun, W., Persson, K.
        A., & Ong, S. P. (2016). Data Descripter: Surface energies of elemental
        crystals. Scientific Data, 3(160080), 1–13.
        http://dx.doi.org/10.1038/sdata.2016.80

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            miller_index (list of integer): The miller index of the surface.
            e.g., [3, 2, 1]. If miller_index is provided, only one dictionary
            of this specific plane will be returned.
            inc_structures (bool): Include final surface slab structures.
                These are unnecessary for Wulff shape construction.
        Returns:
            Surface data for material. Energies are given in SI units (J/m^2).
        """
        req = "/materials/{}/surfaces".format(material_id)
        if inc_structures:
            req += "?include_structures=true"

        if miller_index:
            surf_data_dict = self._make_request(req)
            surf_list = surf_data_dict["surfaces"]
            ucell = self.get_structure_by_material_id(
                material_id, conventional_unit_cell=True
            )
            eq_indices = get_symmetrically_equivalent_miller_indices(
                ucell, miller_index
            )
            for one_surf in surf_list:
                if tuple(one_surf["miller_index"]) in eq_indices:
                    return one_surf
            raise ValueError("Bad miller index.")
        return self._make_request(req)

    def get_wulff_shape(self, material_id):
        """
        Constructs a Wulff shape for a material.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
        Returns:
            pymatgen.analysis.wulff.WulffShape
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.analysis.wulff import WulffShape

        structure = self.get_structure_by_material_id(material_id)
        surfaces = self.get_surface_data(material_id)["surfaces"]
        lattice = (
            SpacegroupAnalyzer(structure).get_conventional_standard_structure().lattice
        )
        miller_energy_map = {}
        for surf in surfaces:
            miller = tuple(surf["miller_index"])
            # Prefer reconstructed surfaces, which have lower surface energies.
            if (miller not in miller_energy_map) or surf["is_reconstructed"]:
                miller_energy_map[miller] = surf["surface_energy"]
        millers, energies = zip(*miller_energy_map.items())
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
        """
        Gets grain boundary data for a material.

        Args:
            material_id (str): Materials Project material_id, e.g., 'mp-129'.
            pretty_formula (str): The formula of metals. e.g., 'Fe'
            sigma(int): The sigma value of a certain type of grain boundary
            gb_plane(list of integer): The Miller index of grain
            boundary plane. e.g., [1, 1, 1]
            rotation_axis(list of integer): The Miller index of rotation
            axis. e.g., [1, 0, 0], [1, 1, 0], and [1, 1, 1]
            Sigma value is determined by the combination of rotation axis and
            rotation angle. The five degrees of freedom (DOF) of one grain boundary
            include: rotation axis (2 DOFs), rotation angle (1 DOF), and grain
            boundary plane (2 DOFs).
            include_work_of_separation (bool): whether to include the work of separation
            (in unit of (J/m^2)). If you want to query the work of separation, please
            specify the material_id.


        Returns:
            A list of grain boundaries that satisfy the query conditions (sigma, gb_plane).
            Energies are given in SI units (J/m^2).
        """
        if gb_plane:
            gb_plane = ",".join([str(i) for i in gb_plane])
        if rotation_axis:
            rotation_axis = ",".join([str(i) for i in rotation_axis])

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
            for i, gb_dict in enumerate(list_of_gbs):
                gb_energy = gb_dict["gb_energy"]
                gb_plane_int = gb_dict["gb_plane"]
                surface_energy = self.get_surface_data(
                    material_id=material_id, miller_index=gb_plane_int
                )["surface_energy"]
                wsep = (
                    2 * surface_energy - gb_energy
                )  # calculate the work of separation
                gb_dict["work_of_separation"] = wsep
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
        """
        Gets critical reactions between two reactants.

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
            "reactants": " ".join([reactant1, reactant2]),
            "open_el": open_el,
            "relative_mu": relative_mu,
            "use_hull_energy": use_hull_energy,
        }
        return self._make_request(
            "/interface_reactions", payload=payload, method="POST"
        )

    def get_download_info(self, material_ids, task_types=None, file_patterns=None):
        """
        get a list of URLs to retrieve raw VASP output files from the NoMaD repository

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
        task_types = (
            [t.value for t in task_types if isinstance(t, TaskType)]
            if task_types
            else []
        )

        meta = defaultdict(list)
        for doc in self.query(
            {"material_id": {"$in": material_ids}}, ["material_id", "blessed_tasks"]
        ):

            for task_type, task_id in doc["blessed_tasks"].items():
                if task_types and task_type not in task_types:
                    continue
                meta[doc["material_id"]].append(
                    {"task_id": task_id, "task_type": task_type}
                )

        if not meta:
            raise ValueError("No tasks found.")

        # return a list of URLs for NoMaD Downloads containing the list of files
        # for every external_id in `task_ids`
        prefix = "http://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/mp/api/raw/query?"
        if file_patterns is not None:
            for file_pattern in file_patterns:
                prefix += f"file_pattern={file_pattern}&"
        prefix += "external_id="

        # NOTE: IE has 2kb URL char limit
        nmax = int((2000 - len(prefix)) / 11)  # mp-<7-digit> + , = 11
        task_ids = [t["task_id"] for tl in meta.values() for t in tl]
        chunks = get_chunks(task_ids, size=nmax)
        urls = [prefix + ",".join(tids) for tids in chunks]
        return meta, urls

    @staticmethod
    def parse_criteria(criteria_string):
        """
        Parses a powerful and simple string criteria and generates a proper
        mongo syntax criteria.

        Args:
            criteria_string (str): A string representing a search criteria.
                Also supports wild cards. E.g.,
                something like "*2O" gets converted to
                {'pretty_formula': {'$in': [u'B2O', u'Xe2O', u"Li2O", ...]}}

                Other syntax examples:
                    mp-1234: Interpreted as a Materials ID.
                    Fe2O3 or *2O3: Interpreted as reduced formulas.
                    Li-Fe-O or *-Fe-O: Interpreted as chemical systems.

                You can mix and match with spaces, which are interpreted as
                "OR". E.g., "mp-1234 FeO" means query for all compounds with
                reduced formula FeO or with materials_id mp-1234.

        Returns:
            A mongo query dict.
        """
        toks = criteria_string.split()

        def parse_sym(sym):
            if sym == "*":
                return [el.symbol for el in Element]
            m = re.match(r"\{(.*)\}", sym)
            if m:
                return [s.strip() for s in m.group(1).split(",")]
            return [sym]

        def parse_tok(t):
            if re.match(r"\w+-\d+", t):
                return {"task_id": t}
            if "-" in t:
                elements = [parse_sym(sym) for sym in t.split("-")]
                chemsyss = []
                for cs in itertools.product(*elements):
                    if len(set(cs)) == len(cs):
                        # Check for valid symbols
                        cs = [Element(s).symbol for s in cs]
                        chemsyss.append("-".join(sorted(cs)))
                return {"chemsys": {"$in": chemsyss}}
            all_formulas = set()
            explicit_els = []
            wild_card_els = []
            for sym in re.findall(r"(\*[\.\d]*|\{.*\}[\.\d]*|[A-Z][a-z]*)[\.\d]*", t):
                if ("*" in sym) or ("{" in sym):
                    wild_card_els.append(sym)
                else:
                    m = re.match(r"([A-Z][a-z]*)[\.\d]*", sym)
                    explicit_els.append(m.group(1))
            nelements = len(wild_card_els) + len(set(explicit_els))
            parts = re.split(r"(\*|\{.*\})", t)
            parts = [parse_sym(s) for s in parts if s != ""]
            for f in itertools.product(*parts):
                c = Composition("".join(f))
                if len(c) == nelements:
                    # Check for valid Elements in keys.
                    for e in c.keys():
                        Element(e.symbol)
                    all_formulas.add(c.reduced_formula)
            return {"pretty_formula": {"$in": list(all_formulas)}}

        if len(toks) == 1:
            return parse_tok(toks[0])
        return {"$or": list(map(parse_tok, toks))}


class MPRestError(Exception):
    """
    Exception class for MPRestAdaptor.
    Raised when the query has problems, e.g., bad query format.
    """

    pass
