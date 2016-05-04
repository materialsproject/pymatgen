# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
This module define the various drones used to assimilate data.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 18, 2012"

import abc
import os
import re
import glob
import logging
import fnmatch
import json

import six
from six.moves import zip

from monty.io import zopen
from pymatgen.io.vasp.inputs import Incar, Potcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun, Oszicar, Dynmat
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.entries.computed_entries import ComputedEntry, \
    ComputedStructureEntry
from monty.json import MSONable

logger = logging.getLogger(__name__)


class AbstractDrone(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    Abstract drone class that defines the various methods that must be
    implemented by drones. Because of the quirky nature of Python"s
    multiprocessing, the intermediate data representations has to be in the
    form of python primitives. So all objects that drones work with must be
    MSONable. All drones must also implement the standard MSONable as_dict() and
    from_dict API.
    """

    @abc.abstractmethod
    def assimilate(self, path):
        """
        Assimilate data in a directory path into a pymatgen object. Because of
        the quirky nature of Python"s multiprocessing, the object must support
        pymatgen"s as_dict() for parallel processing.

        Args:
            path: directory path

        Returns:
            An assimilated object
        """
        return

    @abc.abstractmethod
    def get_valid_paths(self, path):
        """
        Checks if path contains valid data for assimilation, and then returns
        the valid paths. The paths returned can be a list of directory or file
        paths, depending on what kind of data you are assimilating. For
        example, if you are assimilating VASP runs, you are only interested in
        directories containing vasprun.xml files. On the other hand, if you are
        interested converting all POSCARs in a directory tree to cifs for
        example, you will want the file paths.

        Args:
            path: input path as a tuple generated from os.walk, i.e.,
                (parent, subdirs, files).

        Returns:
            List of valid dir/file paths for assimilation
        """
        return


class VaspToComputedEntryDrone(AbstractDrone):
    """
    VaspToEntryDrone assimilates directories containing vasp output to
    ComputedEntry/ComputedStructureEntry objects. There are some restrictions
    on the valid directory structures:

    1. There can be only one vasp run in each directory.
    2. Directories designated "relax1", "relax2" are considered to be 2 parts
       of an aflow style run, and only "relax2" is parsed.
    3. The drone parses only the vasprun.xml file.


    Args:
        inc_structure (bool): Set to True if you want
            ComputedStructureEntries to be returned instead of
            ComputedEntries.
        parameters (list): Input parameters to include. It has to be one of
            the properties supported by the Vasprun object. See
            :class:`pymatgen.io.vasp.Vasprun`. If parameters == None,
            a default set of parameters that are necessary for typical
            post-processing will be set.
        data (list): Output data to include. Has to be one of the properties
            supported by the Vasprun object.
    """

    def __init__(self, inc_structure=False, parameters=None, data=None):
        self._inc_structure = inc_structure
        self._parameters = {"is_hubbard", "hubbards", "potcar_spec",
                            "potcar_symbols", "run_type"}
        if parameters:
            self._parameters.update(parameters)
        self._data = data if data else []

    def assimilate(self, path):
        files = os.listdir(path)
        if "relax1" in files and "relax2" in files:
            filepath = glob.glob(os.path.join(path, "relax2",
                                              "vasprun.xml*"))[0]
        else:
            vasprun_files = glob.glob(os.path.join(path, "vasprun.xml*"))
            filepath = None
            if len(vasprun_files) == 1:
                filepath = vasprun_files[0]
            elif len(vasprun_files) > 1:
                """
                This is a bit confusing, since there maybe be multi-steps. By
                default, assimilate will try to find a file simply named
                vasprun.xml, vasprun.xml.bz2, or vasprun.xml.gz.  Failing which
                it will try to get a relax2 from an aflow style run if
                possible. Or else, a randomly chosen file containing
                vasprun.xml is chosen.
                """
                for fname in vasprun_files:
                    if os.path.basename(fname) in ["vasprun.xml",
                                                   "vasprun.xml.gz",
                                                   "vasprun.xml.bz2"]:
                        filepath = fname
                        break
                    if re.search("relax2", fname):
                        filepath = fname
                        break
                    filepath = fname

        try:
            vasprun = Vasprun(filepath)
        except Exception as ex:
            logger.debug("error in {}: {}".format(filepath, ex))
            return None

        entry = vasprun.get_computed_entry(self._inc_structure,
                                           parameters=self._parameters,
                                           data=self._data)
        entry.parameters["history"] = _get_transformation_history(path)
        return entry

    def get_valid_paths(self, path):
        (parent, subdirs, files) = path
        if "relax1" in subdirs and "relax2" in subdirs:
            return [parent]
        if (not parent.endswith("/relax1")) and \
           (not parent.endswith("/relax2")) and (
               len(glob.glob(os.path.join(parent, "vasprun.xml*"))) > 0 or (
               len(glob.glob(os.path.join(parent, "POSCAR*"))) > 0 and
               len(glob.glob(os.path.join(parent, "OSZICAR*"))) > 0)
           ):
            return [parent]
        return []

    def __str__(self):
        return " VaspToComputedEntryDrone"

    def as_dict(self):
        return {"init_args": {"inc_structure": self._inc_structure,
                              "parameters": self._parameters,
                              "data": self._data},
                "version": __version__,
                "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])


class SimpleVaspToComputedEntryDrone(VaspToComputedEntryDrone):
    """
    A simpler VaspToComputedEntryDrone. Instead of parsing vasprun.xml, it
    parses only the INCAR, POTCAR, OSZICAR and KPOINTS files, which are much
    smaller and faster to parse. However, much fewer properties are available
    compared to the standard VaspToComputedEntryDrone.

    Args:
        inc_structure (bool): Set to True if you want
            ComputedStructureEntries to be returned instead of
            ComputedEntries. Structure will be parsed from the CONTCAR.
    """

    def __init__(self, inc_structure=False):
        self._inc_structure = inc_structure
        self._parameters = {"is_hubbard", "hubbards", "potcar_spec",
                            "run_type"}

    def assimilate(self, path):
        files = os.listdir(path)
        try:
            files_to_parse = {}
            if "relax1" in files and "relax2" in files:
                for filename in ("INCAR", "POTCAR", "POSCAR"):
                    search_str = os.path.join(path, "relax1", filename + "*")
                    files_to_parse[filename] = glob.glob(search_str)[0]
                for filename in ("CONTCAR", "OSZICAR"):
                    search_str = os.path.join(path, "relax2", filename + "*")
                    files_to_parse[filename] = glob.glob(search_str)[-1]
            else:
                for filename in (
                    "INCAR", "POTCAR", "CONTCAR", "OSZICAR", "POSCAR", "DYNMAT"
                ):
                    files = glob.glob(os.path.join(path, filename + "*"))
                    if len(files) < 1:
                        continue
                    if len(files) == 1 or filename == "INCAR" or \
                       filename == "POTCAR" or filename == "DYNMAT":
                        files_to_parse[filename] = files[-1]\
                            if filename == "POTCAR" else files[0]
                    elif len(files) > 1:
                        """
                        This is a bit confusing, since there maybe be
                        multiple steps. By default, assimilate will try to find
                        a file simply named filename, filename.bz2, or
                        filename.gz.  Failing which it will try to get a relax2
                        from a custodian double relaxation style run if
                        possible. Or else, a random file is chosen.
                        """
                        for fname in files:
                            if fnmatch.fnmatch(os.path.basename(fname),
                                               "{}(\.gz|\.bz2)*"
                                               .format(filename)):
                                files_to_parse[filename] = fname
                                break
                            if fname == "POSCAR" and \
                                    re.search("relax1", fname):
                                files_to_parse[filename] = fname
                                break
                            if (fname in ("CONTCAR", "OSZICAR") and
                                    re.search("relax2", fname)):
                                files_to_parse[filename] = fname
                                break
                            files_to_parse[filename] = fname

            poscar, contcar, incar, potcar, oszicar, dynmat = [None]*6
            if 'POSCAR' in files_to_parse:
                poscar = Poscar.from_file(files_to_parse["POSCAR"])
            if 'CONTCAR' in files_to_parse:
                contcar = Poscar.from_file(files_to_parse["CONTCAR"])
            if 'INCAR' in files_to_parse:
                incar = Incar.from_file(files_to_parse["INCAR"])
            if 'POTCAR' in files_to_parse:
                potcar = Potcar.from_file(files_to_parse["POTCAR"])
            if 'OSZICAR' in files_to_parse:
                oszicar = Oszicar(files_to_parse["OSZICAR"])
            if 'DYNMAT' in files_to_parse:
                dynmat = Dynmat(files_to_parse["DYNMAT"])

            param = {"hubbards":{}}
            if poscar is not None and incar is not None and "LDAUU" in incar:
                param["hubbards"] = dict(zip(poscar.site_symbols,
                                             incar["LDAUU"]))
            param["is_hubbard"] = (
                incar.get("LDAU", False) and sum(param["hubbards"].values()) > 0
            ) if incar is not None else False
            param["run_type"] = None
            if incar is not None:
                param["run_type"] = "GGA+U" if param["is_hubbard"] else "GGA"
            param["history"] = _get_transformation_history(path)
            param["potcar_spec"] = potcar.spec if potcar is not None else None
            energy = oszicar.final_energy if oszicar is not None else 1e10
            structure = contcar.structure if contcar is not None\
                else poscar.structure
            initial_vol = poscar.structure.volume if poscar is not None else \
                None
            final_vol = contcar.structure.volume if contcar is not None else \
                None
            delta_volume = None
            if initial_vol is not None and final_vol is not None:
                delta_volume = (final_vol / initial_vol - 1)
            data = {"filename": path, "delta_volume": delta_volume}
            if dynmat is not None:
                data['phonon_frequencies'] = dynmat.get_phonon_frequencies()
            if self._inc_structure:
                entry = ComputedStructureEntry(
                    structure, energy, parameters=param, data=data
                )
            else:
                entry = ComputedEntry(
                  structure.composition, energy, parameters=param, data=data
                )
            return entry

        except Exception as ex:
            logger.debug("error in {}: {}".format(path, ex))
            return None

    def __str__(self):
        return "SimpleVaspToComputedEntryDrone"

    def as_dict(self):
        return {"init_args": {"inc_structure": self._inc_structure},
                "version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])


class GaussianToComputedEntryDrone(AbstractDrone):
    """
    GaussianToEntryDrone assimilates directories containing Gaussian output to
    ComputedEntry/ComputedStructureEntry objects. By default, it is assumed
    that Gaussian output files have a ".log" extension.

    Args:
        inc_structure (bool): Set to True if you want
            ComputedStructureEntries to be returned instead of
            ComputedEntries.
        parameters (list): Input parameters to include. It has to be one of
            the properties supported by the GaussianOutput object. See
            :class:`pymatgen.io.gaussianio GaussianOutput`. The parameters
            have to be one of python"s primitive types, i.e., list, dict of
            strings and integers. If parameters == None, a default set of
            parameters will be set.
        data (list): Output data to include. Has to be one of the properties
            supported by the GaussianOutput object. The parameters have to
            be one of python"s primitive types, i.e. list, dict of strings
            and integers. If data == None, a default set will be set.
        file_extensions (list):
            File extensions to be considered as Gaussian output files.
            Defaults to just the typical "log" extension.

    .. note::

        Like the GaussianOutput class, this is still in early beta.
    """

    def __init__(self, inc_structure=False, parameters=None, data=None,
                 file_extensions=(".log",)):
        self._inc_structure = inc_structure
        self._parameters = {"functional", "basis_set", "charge", "spin_mult",
                            "route"}

        if parameters:
            self._parameters.update(parameters)

        self._data = {"stationary_type", "properly_terminated"}
        if data:
            self._data.update(data)

        self._file_extensions = file_extensions

    def assimilate(self, path):
        try:
            gaurun = GaussianOutput(path)
        except Exception as ex:
            logger.debug("error in {}: {}".format(path, ex))
            return None
        param = {}
        for p in self._parameters:
            param[p] = getattr(gaurun, p)
        data = {}
        for d in self._data:
            data[d] = getattr(gaurun, d)
        if self._inc_structure:
            entry = ComputedStructureEntry(gaurun.final_structure,
                                           gaurun.final_energy,
                                           parameters=param,
                                           data=data)
        else:
            entry = ComputedEntry(gaurun.final_structure.composition,
                                  gaurun.final_energy, parameters=param,
                                  data=data)
        return entry

    def get_valid_paths(self, path):
        (parent, subdirs, files) = path
        return [os.path.join(parent, f) for f in files
                if os.path.splitext(f)[1] in self._file_extensions]

    def __str__(self):
        return " GaussianToComputedEntryDrone"

    def as_dict(self):
        return {"init_args": {"inc_structure": self._inc_structure,
                              "parameters": self._parameters,
                              "data": self._data,
                              "file_extensions": self._file_extensions},
                "version": __version__, "@module": self.__class__.__module__,
                "@class": self.__class__.__name__}

    @classmethod
    def from_dict(cls, d):
        return cls(**d["init_args"])


def _get_transformation_history(path):
    """
    Checks for a transformations.json* file and returns the history.
    """
    trans_json = glob.glob(os.path.join(path, "transformations.json*"))
    if trans_json:
        try:
            with zopen(trans_json[0]) as f:
                return json.load(f)["history"]
        except:
            return None
    return None
