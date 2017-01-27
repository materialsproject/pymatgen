# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import itertools
import json
import re
import warnings

from monty.json import MontyDecoder, MontyEncoder
from six import string_types

from pymatgen import SETTINGS

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure

from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry
from pymatgen.entries.exp_entries import ExpEntry


"""
This module provides classes to interface with the Materials Project REST
API v2 to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your dashboard at
https://www.materialsproject.org/dashboard.
"""

__author__ = "Shyue Ping Ong, Shreyas Cholia"
__credits__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 22, 2013"


class MPRester(object):
    """
    A class to conveniently interface with the Materials Project REST
    interface. The recommended way to use MPRester is with the "with" context
    manager to ensure that sessions are properly closed after usage::

        with MPRester("API_KEY") as m:
            do_something

    MPRester uses the "requests" package, which provides for HTTP connection
    pooling. All connections are made via https for security.

    .. note::

        The Materials Project recently switched to using string ids with a
        "mp-" prefix for greater flexibility going forward. The MPRester
        should still work as intended if you provide the proper string ids.

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
            address, but can be changed to other urls implementing a similar
            interface.
    """

    supported_properties = ("energy", "energy_per_atom", "volume",
                            "formation_energy_per_atom", "nsites",
                            "unit_cell_formula", "pretty_formula",
                            "is_hubbard", "elements", "nelements",
                            "e_above_hull", "hubbards", "is_compatible",
                            "spacegroup", "task_ids", "band_gap", "density",
                            "icsd_id", "icsd_ids", "cif", "total_magnetization",
                            "material_id", "oxide_type", "tags", "elasticity")

    supported_task_properties = ("energy", "energy_per_atom", "volume",
                                 "formation_energy_per_atom", "nsites",
                                 "unit_cell_formula", "pretty_formula",
                                 "is_hubbard",
                                 "elements", "nelements", "e_above_hull",
                                 "hubbards",
                                 "is_compatible", "spacegroup",
                                 "band_gap", "density", "icsd_id", "cif")

    def __init__(self, api_key=None,
                 endpoint="https://www.materialsproject.org/rest/v2"):
        if api_key is not None:
            self.api_key = api_key
        else:
            self.api_key = SETTINGS.get("PMG_MAPI_KEY", "")
        self.preamble = endpoint
        import requests
        self.session = requests.Session()
        self.session.headers = {"x-api-key": self.api_key}

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

    def _make_request(self, sub_url, payload=None, method="GET",
                      mp_decode=True):
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
                else:
                    raise MPRestError(data["error"])

            raise MPRestError("REST query returned with error status code {}"
                              .format(response.status_code))

        except Exception as ex:
            msg = "{}. Content: {}".format(str(ex), response.content)\
                if hasattr(response, "content") else str(ex)
            raise MPRestError(msg)

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
            A list of matching structures.

        Raises:
            MPRestError
        """
        try:
            if isinstance(filename_or_structure, string_types):
                s = Structure.from_file(filename_or_structure)
            elif isinstance(filename_or_structure, Structure):
                s = filename_or_structure
            else:
                raise MPRestError("Provide filename or Structure object.")
            payload = {'structure': json.dumps(s.as_dict(), cls=MontyEncoder)}
            response = self.session.post(
                '{}/find_structure'.format(self.preamble), data=payload
            )
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp['valid_response']:
                    return resp['response']
                else:
                    raise MPRestError(resp["error"])
            raise MPRestError("REST error with status code {} and error {}"
                              .format(response.status_code, response.text))
        except Exception as ex:
            raise MPRestError(str(ex))

    def get_entries(self, chemsys_formula_id_criteria, compatible_only=True,
                    inc_structure=None, property_data=None):
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
                ComputedEntries. If inc_structure="final",
                ComputedStructureEntries with final structures are returned.
                Otherwise, ComputedStructureEntries with initial structures
                are returned.
            property_data (list): Specify additional properties to include in
                entry.data. If None, no data. Should be a subset of
                supported_properties.

        Returns:
            List of ComputedEntry or ComputedStructureEntry objects.
        """
        # TODO: This is a very hackish way of doing this. It should be fixed
        # on the REST end.
        params = ["run_type", "is_hubbard", "pseudo_potential", "hubbards",
                  "potcar_symbols"]
        if compatible_only:
            props = ["energy", "unit_cell_formula", "task_id"] + params
            if property_data:
                props += property_data
            if inc_structure:
                if inc_structure == "final":
                    props.append("structure")
                else:
                    props.append("initial_structure")

            if not isinstance(chemsys_formula_id_criteria, dict):
                criteria = MPRester.parse_criteria(chemsys_formula_id_criteria)
            else:
                criteria = chemsys_formula_id_criteria

            data = self.query(criteria, props)

            entries = []
            for d in data:
                d["potcar_symbols"] = [
                    "%s %s" % (d["pseudo_potential"]["functional"], l)
                    for l in d["pseudo_potential"]["labels"]]
                data = {k: d[k] for k in property_data} if property_data else None
                if not inc_structure:
                    e = ComputedEntry(d["unit_cell_formula"], d["energy"],
                                      parameters={k: d[k] for k in params},
                                      data=data,
                                      entry_id=d["task_id"])

                else:
                    s = d["structure"] if inc_structure == "final" else d[
                        "initial_structure"]
                    e = ComputedStructureEntry(
                        s, d["energy"],
                        parameters={k: d[k] for k in params},
                        data=data,
                        entry_id=d["task_id"])
                entries.append(e)
            from pymatgen.entries.compatibility import \
                MaterialsProjectCompatibility
            entries = MaterialsProjectCompatibility().process_entries(entries)
        else:
            entries = []
            for d in self.get_data(chemsys_formula_id_criteria,
                                   prop="task_ids"):
                for i in d["task_ids"]:
                    e = self.get_task_data(i, prop="entry")
                    e = e[0]["entry"]
                    if inc_structure:
                        s = self.get_task_data(i,
                                               prop="structure")[0]["structure"]
                        e = ComputedStructureEntry(
                            s, e.energy, e.correction, e.parameters, e.data,
                            e.entry_id)
                    entries.append(e)

        return entries

    def get_structure_by_material_id(self, material_id, final=True):
        """
        Get a Structure corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).
            final (bool): Whether to get the final structure, or the initial
                (pre-relaxation) structure. Defaults to True.

        Returns:
            Structure object.
        """
        prop = "final_structure" if final else "initial_structure"
        data = self.get_data(material_id, prop=prop)
        return data[0][prop]

    def get_entry_by_material_id(self, material_id, compatible_only=True,
                                 inc_structure=None, property_data=None):
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

        Returns:
            ComputedEntry or ComputedStructureEntry object.
        """
        data = self.get_entries(material_id, compatible_only=compatible_only,
                                inc_structure=inc_structure, property_data=property_data)
        return data[0]

    def get_dos_by_material_id(self, material_id):
        """
        Get a Dos corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).

        Returns:
            A Dos object.
        """
        data = self.get_data(material_id, prop="dos")
        return data[0]["dos"]

    def get_bandstructure_by_material_id(self, material_id):
        """
        Get a BandStructure corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (an int).

        Returns:
            A BandStructure object.
        """
        data = self.get_data(material_id, prop="bandstructure")
        return data[0]["bandstructure"]

    def get_entries_in_chemsys(self, elements, compatible_only=True,
                               inc_structure=None, property_data=None):
        """
        Helper method to get a list of ComputedEntries in a chemical system.
        For example, elements = ["Li", "Fe", "O"] will return a list of all
        entries in the Li-Fe-O chemical system, i.e., all LixOy,
        FexOy, LixFey, LixFeyOz, Li, Fe and O phases. Extremely useful for
        creating phase diagrams of entire chemical systems.

        Args:
            elements ([str]): List of element symbols, e.g., ["Li", "Fe",
                "O"].
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

        Returns:
            List of ComputedEntries.
        """
        entries = []
        for i in range(len(elements)):
            for els in itertools.combinations(elements, i + 1):
                entries.extend(
                    self.get_entries(
                        "-".join(els), compatible_only=compatible_only,
                        inc_structure=inc_structure,
                        property_data=property_data))
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

        return ExpEntry(Composition(formula),
                        self.get_exp_thermo_data(formula))

    def query(self, criteria, properties, mp_decode=True):
        """
        Performs an advanced query, which is a Mongo-like syntax for directly
        querying the Materials Project database via the query rest interface.
        Please refer to the Materials Project REST wiki
        https://materialsproject.org/wiki/index.php/The_Materials_API#query
        on the query language and supported criteria and properties.
        Essentially, any supported properties within MPRester should be
        supported in query.

        Query allows an advanced developer to perform queries which are
        otherwise too cumbersome to perform using the standard convenience
        methods.

        It is highly recommended that you consult the Materials API
        documentation at http://bit.ly/materialsapi, which provides a
        comprehensive explanation of the document schema used in the
        Materials Project and how best to query for the relevant information
        you need.

        Args:
            criteria (str/dict): Criteria of the query as a string or
                mongo-style dict.

                If string, it supports a powerful but simple string criteria.
                E.g., "Fe2O3" means search for materials with reduced_formula
                Fe2O3. Wild cards are also supported. E.g., "\*2O" means get
                all materials whose formula can be formed as \*2O, e.g.,
                Li2O, K2O, etc.

                Other syntax examples:
                mp-1234: Interpreted as a Materials ID.
                Fe2O3 or \*2O3: Interpreted as reduced formulas.
                Li-Fe-O or \*-Fe-O: Interpreted as chemical systems.

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
            criteria = MPRester.parse_criteria(criteria)
        payload = {"criteria": json.dumps(criteria),
                   "properties": json.dumps(properties)}
        return self._make_request("/query", payload=payload, method="POST",
                                  mp_decode=mp_decode)

    def submit_structures(self, structures, authors, projects=None,
                          references='', remarks=None, data=None,
                          histories=None, created_at=None):
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
        from pymatgen.matproj.snl import StructureNL
        snl_list = StructureNL.from_structures(structures, authors, projects,
                                               references, remarks, data,
                                               histories, created_at)
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
            response = self.session.post("{}/snl/submit".format(self.preamble),
                                         data=payload)
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp['inserted_ids']
                else:
                    raise MPRestError(resp["error"])

            raise MPRestError("REST error with status code {} and error {}"
                              .format(response.status_code, response.text))

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
                "{}/snl/delete".format(self.preamble), data=payload)

            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp
                else:
                    raise MPRestError(resp["error"])

            raise MPRestError("REST error with status code {} and error {}"
                              .format(response.status_code, response.text))

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
            response = self.session.post("{}/snl/query".format(self.preamble),
                                         data=payload)
            if response.status_code in [200, 400]:
                resp = json.loads(response.text)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp["response"]
                else:
                    raise MPRestError(resp["error"])

            raise MPRestError("REST error with status code {} and error {}"
                              .format(response.status_code, response.text))

        except Exception as ex:
            raise MPRestError(str(ex))

    def submit_vasp_directory(self, rootdir, authors, projects=None,
                              references='', remarks=None, master_data=None,
                              master_history=None, created_at=None,
                              ncpus=None):
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
        drone = VaspToComputedEntryDrone(inc_structure=True,
                                         data=["filename",
                                               "initial_structure"])
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
                    "initial_structure": e.data["initial_structure"].as_dict()
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
            structures, authors, projects=projects, references=references,
            remarks=remarks, data=metadata, histories=histories,
            created_at=created_at)

    def get_stability(self, entries):
        """
        Returns the stability of all entries.
        """
        try:
            payload = {"entries": json.dumps(entries, cls=MontyEncoder)}
            response = self.session.post("{}/phase_diagram/calculate_stability"
                                         .format(self.preamble), data=payload)
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=MontyDecoder)
                if resp["valid_response"]:
                    if resp.get("warning"):
                        warnings.warn(resp["warning"])
                    return resp["response"]
                else:
                    raise MPRestError(resp["error"])
            raise MPRestError("REST error with status code {} and error {}"
                              .format(response.status_code, response.text))
        except Exception as ex:
            raise MPRestError(str(ex))

    def get_reaction(self, reactants, products):
        """
        Gets a reaction from the Materials Project.

        Args:
            reactants ([str]): List of formulas
            products ([str]): List of formulas

        Returns:
            rxn
        """
        return self._make_request("/reaction",
                                  payload={"reactants[]": reactants,
                                           "products[]": products}, mp_decode=False)

    def get_substrates(self, material_id, number=50, orient=None):
        """
        Get a substrate list for a material id. The list is in order of
        increasing elastic energy if a elastic tensor is available for
        the material_id. Otherwise the list is in order of increasing
        matching area.

        Args:
            material_id (str): Materials Project material_id, e.g. 'mp-123'.
            orient (list) : substrate orientation to look for
            number (int) : number of substrates to return;
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

    @staticmethod
    def parse_criteria(criteria_string):
        """
        Parses a powerful and simple string criteria and generates a proper
        mongo syntax criteria.

        Args:
            criteria_string (str): A string representing a search criteria.
                Also supports wild cards. E.g.,
                something like "\*2O" gets converted to
                {'pretty_formula': {'$in': [u'B2O', u'Xe2O', u"Li2O", ...]}}

                Other syntax examples:
                    mp-1234: Interpreted as a Materials ID.
                    Fe2O3 or \*2O3: Interpreted as reduced formulas.
                    Li-Fe-O or \*-Fe-O: Interpreted as chemical systems.

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
            else:
                m = re.match("\{(.*)\}", sym)
                if m:
                    return [s.strip() for s in m.group(1).split(",")]
                else:
                    return [sym]

        def parse_tok(t):
            if re.match("\w+-\d+", t):
                return {"task_id": t}
            elif "-" in t:
                elements = [parse_sym(sym) for sym in t.split("-")]
                chemsyss = []
                for cs in itertools.product(*elements):
                    if len(set(cs)) == len(cs):
                        # Check for valid symbols
                        cs = [Element(s).symbol for s in cs]
                        chemsyss.append("-".join(sorted(cs)))
                return {"chemsys": {"$in": chemsyss}}
            else:
                all_formulas = set()
                explicit_els = []
                wild_card_els = []
                for sym in re.findall(
                        r"(\*[\.\d]*|\{.*\}[\.\d]*|[A-Z][a-z]*)[\.\d]*", t):
                    if ("*" in sym) or ("{" in sym):
                        wild_card_els.append(sym)
                    else:
                        m = re.match("([A-Z][a-z]*)[\.\d]*", sym)
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
        else:
            return {"$or": list(map(parse_tok, toks))}


class MPRestError(Exception):
    """
    Exception class for MPRestAdaptor.
    Raised when the query has problems, e.g., bad query format.
    """
    pass

