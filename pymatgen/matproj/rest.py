#!/usr/bin/env python

"""
This module provides classes to interface with the Materials Project REST
API to enable the creation of data structures and pymatgen objects using
Materials Project data.

To make use of the Materials API, you need to be a registered user of the
Materials Project, and obtain an API key by going to your profile at
https://www.materialsproject.org/profile.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Shreyas Cholia"
__credits__ = "Anubhav Jain"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Feb 22, 2013"

import os
import requests
import json
import warnings

from pymatgen import Composition, PMGJSONDecoder
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.entries.exp_entries import ExpEntry
from pymatgen.io.vaspio_set import DictVaspInputSet
from pymatgen.apps.borg.hive import VaspToComputedEntryDrone
from pymatgen.apps.borg.queen import BorgQueen
from pymatgen.matproj.snl import StructureNL
from pymatgen.serializers.json_coders import PMGJSONEncoder


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
            REST interface. Please apply on the Materials Project website for
            one.
            If this is None, the code will check if there is a "MAPI_KEY"
            environment variable set. If so, it will use that environment
            variable. This makes easier for heavy users to simply add
            this environment variable to their setups and MPRester can
            then be called without any arguments.
        host (str): Url of host to access the MaterialsProject REST interface.
            Defaults to the standard Materials Project REST address, but
            can be changed to other urls implementing a similar interface.
    """

    supported_properties = ("energy", "energy_per_atom", "volume",
                            "formation_energy_per_atom", "nsites",
                            "unit_cell_formula", "pretty_formula",
                            "is_hubbard", "elements", "nelements",
                            "e_above_hull", "hubbards", "is_compatible",
                            "spacegroup", "task_ids", "band_gap", "density",
                            "icsd_id", "cif", "total_magnetization",
                            "material_id")

    def __init__(self, api_key=None, host="www.materialsproject.org"):
        if api_key is not None:
            self.api_key = api_key
        else:
            self.api_key = os.environ.get("MAPI_KEY", "")
        self.preamble = "https://{}/rest/v1".format(host)
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

    def get_materials_id_from_task_id(self, task_id):
        """
        Returns a new MP materials id from a task id (which can be
        equivalent to an old materials id)

        Args:
            task_id (str): A task id.
        """
        data = self.mpquery({"task_ids": task_id}, ["task_id"])
        return data[0]["task_id"]

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
                MPRester.supported_properties. Leave as empty string for a
                general list of useful properties.
        """
        if prop:
            url = "{}/materials/{}/{}/{}".format(
                self.preamble, chemsys_formula_id, data_type, prop)
        else:
            url = "{}/materials/{}/{}".format(
                self.preamble, chemsys_formula_id, data_type)

        response = None
        try:
            response = self.session.get(url)
            if response.status_code in [200, 400]:
                data = json.loads(response.text, cls=PMGJSONDecoder)
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
        if prop:
            url = "{}/tasks/{}/{}".format(
                self.preamble, chemsys_formula_id, prop)
        else:
            url = "{}/tasks/{}".format(
                self.preamble, chemsys_formula_id)

        response = None
        try:
            response = self.session.get(url)
            if response.status_code in [200, 400]:
                data = json.loads(response.text, cls=PMGJSONDecoder)
                if data["valid_response"]:
                    if data.get("warning"):
                        warnings.warn(data["warning"])
                    return data["response"]
                else:
                    raise MPRestError(data["error"])

            raise MPRestError("REST query returned with error status code {}"
                              .format(response.status_code))

        except Exception as ex:
            msg = "{}. Content: {}".format(str(ex), response.content) if hasattr(
                response, "content") else str(ex)
            raise MPRestError(msg)


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

    def get_entries(self, chemsys_formula_id, compatible_only=True,
                    inc_structure=None):
        """
        Get a list of ComputedEntries or ComputedStructureEntries corresponding
        to a chemical system, formula, or materials_id.

        Args:
            chemsys_formula_id (str): A chemical system (e.g., Li-Fe-O),
                or formula (e.g., Fe2O3) or materials_id (e.g., mp-1234).
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

        Returns:
            List of ComputedEntry or ComputedStructureEntry objects.
        """
        # TODO: This is a very hackish way of doing this. It should be fixed
        # on the REST end.
        if compatible_only:
            data = self.get_data(chemsys_formula_id, prop="entry")
            entries = [d["entry"] for d in data]
            if inc_structure:
                for i, e in enumerate(entries):
                    s = self.get_structure_by_material_id(
                        e.entry_id, inc_structure == "final")
                    entries[i] = ComputedStructureEntry(
                        s, e.energy, e.correction, e.parameters, e.data,
                        e.entry_id)
            entries = MaterialsProjectCompatibility().process_entries(entries)

        else:
            entries = []
            for d in self.get_data(chemsys_formula_id, prop="task_ids"):
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

    def get_entry_by_material_id(self, material_id):
        """
        Get a ComputedEntry corresponding to a material_id.

        Args:
            material_id (str): Materials Project material_id (a string,
                e.g., mp-1234).

        Returns:
            ComputedEntry object.
        """
        data = self.get_data(material_id, prop="entry")
        return data[0]["entry"]

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
                               inc_structure=None):
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

        Returns:
            List of ComputedEntries.
        """
        return self.get_entries(
            "-".join(elements), compatible_only=compatible_only,
            inc_structure=inc_structure)

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

    def get_vasp_input_set(self, date_string=None):
        """
        Returns the VaspInputSet used by the Materials Project at a
        particular date.

        Args:
            date_string (str): A date string in the format of "YYYY-MM-DD".
                Defaults to None, which means the VaspInputSet today.

        Returns:
            DictVaspInputSet
        """
        url = "{}/parameters/vasp".format(self.preamble)
        payload = {"date": date_string} if date_string else {}
        try:
            response = self.session.get(url, data=payload)
            if response.status_code in [200, 400]:
                data = json.loads(response.text, cls=PMGJSONDecoder)
                if data["valid_response"]:
                    if data.get("warning"):
                        warnings.warn(data["warning"])
                    return DictVaspInputSet("MPVaspInputSet", data["response"])
                else:
                    raise MPRestError(data["error"])

            raise MPRestError("REST query returned with error status code {}"
                              .format(response.status_code))
        except Exception as ex:
            raise MPRestError(str(ex))

    def mpquery(self, criteria, properties):
        """
        Performs an advanced mpquery, which is a Mongo-like syntax for directly
        querying the Materials Project database via the mpquery rest interface.
        Please refer to the Materials Project REST wiki
        https://materialsproject.org/wiki/index.php/The_Materials_API#mpquery
        on the mpquery language and supported criteria and properties.
        Essentially, any supported properties within MPRester should be
        supported in mpquery.

        Mpquery allows an advanced developer to perform queries which are
        otherwise too cumbersome to perform using the standard convenience
        methods.

        Args:
            criteria (dict): Criteria of the query as a mongo-style dict.
                For example, {"elements":{"$in":["Li", "Na", "K"], "$all": [
                "O"]}, "nelements":2} selects all Li, Na and K oxides.
                {"band_gap": {"$gt": 1}} selects all materials with band gaps
                greater than 1 eV.
            properties (list): Properties to request for as a list. For
                example, ["formula", "formation_energy_per_atom"] returns
                the formula and formation energy per atom.

        Returns:
            List of results. E.g.,
            [{u'formula': {u'O': 1, u'Li': 2.0}},
            {u'formula': {u'Na': 2.0, u'O': 2.0}},
            {u'formula': {u'K': 1, u'O': 3.0}},
            ...]
        """
        try:
            payload = {"criteria": json.dumps(criteria),
                       "properties": json.dumps(properties)}
            response = self.session.post(
                "{}/mpquery".format(self.preamble), data=payload)

            if response.status_code in [200, 400]:
                data = json.loads(response.text, cls=PMGJSONDecoder)
                if data["valid_response"]:
                    if data.get("warning"):
                        warnings.warn(data["warning"])
                    return data["response"]["results"]
                else:
                    raise MPRestError(data["error"])

            raise MPRestError("REST query returned with error status code {}"
                              .format(response.status_code))

        except Exception as ex:
            raise MPRestError(str(ex))

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
            references (str): A String in BibTeX format. Again, this applies to all
                structures.
            remarks ([str]): List of Strings ['Remark A', 'Remark B']
            data ([dict]): A list of free form dict. Namespaced at the root
                level with an underscore, e.g. {"_materialsproject":<custom
                data>}. The length of data should be the same as the list of structures
                if not None.
            histories: List of list of dicts - [[{'name':'', 'url':'',
                'description':{}}], ...] The length of histories should be the
                same as the list of structures if not None.
            created_at (datetime): A datetime object

        Returns:
            A list of inserted submission ids.
        """
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
            jsondata = [s.to_dict for s in snl]
            payload = {"snl": json.dumps(jsondata, cls=PMGJSONEncoder)}
            response = self.session.post("{}/snl/submit".format(self.preamble),
                                         data=payload)
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=PMGJSONDecoder)
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
                resp = json.loads(response.text, cls=PMGJSONDecoder)
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
            references (str): A String in BibTeX format. Again, this applies to all
                structures.
            remarks ([str]): List of Strings ['Remark A', 'Remark B']
            master_data (dict): A free form dict. Namespaced at the root
                level with an underscore, e.g. {"_materialsproject":<custom
                data>}. This data is added to all structures detected in the directory,
                in addition to other vasp data on a per structure basis.
            master_history: A master history to be added to all entries.
            created_at (datetime): A datetime object
            ncpus (int): Number of cpus to use in using BorgQueen to
                assimilate. Defaults to None, which means serial.
       """
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
                    "initial_structure": e.data["initial_structure"].to_dict
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
            payload = {"entries": json.dumps(entries, cls=PMGJSONEncoder)}
            response = self.session.post("{}/phase_diagram/calculate_stability"
                                         .format(self.preamble), data=payload)
            if response.status_code in [200, 400]:
                resp = json.loads(response.text, cls=PMGJSONDecoder)
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


class MPRestError(Exception):
    """
    Exception class for MPRestAdaptor.
    Raised when the query has problems, e.g., bad query format.
    """
    pass
