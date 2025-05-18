"""This module define the various drones used to assimilate data."""

from __future__ import annotations

import abc
import logging
import os
import warnings
from glob import glob
from typing import TYPE_CHECKING

import orjson
from monty.io import zopen
from monty.json import MSONable

from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar
from pymatgen.io.vasp.outputs import Dynmat, Oszicar, Vasprun

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from typing_extensions import Self

    from pymatgen.util.typing import PathLike


logger = logging.getLogger(__name__)


class AbstractDrone(MSONable, abc.ABC):
    """Abstract drone class that defines the various methods that must be
    implemented by drones. Because of the quirky nature of Python"s
    multiprocessing, the intermediate data representations has to be in the
    form of python primitives. So all objects that drones work with must be
    MSONable. All drones must also implement the standard MSONable as_dict() and
    from_dict API.
    """

    @abc.abstractmethod
    def assimilate(self, path):
        """Assimilate data in a directory path into a pymatgen object. Because of
        the quirky nature of Python's multiprocessing, the object must support
        pymatgen's as_dict() for parallel processing.

        Args:
            path: directory path

        Returns:
            An assimilated object
        """
        return

    @abc.abstractmethod
    def get_valid_paths(self, path):
        """Check if path contains valid data for assimilation, and then returns
        the valid paths. The paths returned can be a list of directory or file
        paths, depending on what kind of data you are assimilating. For
        example, if you are assimilating VASP runs, you are only interested in
        directories containing vasprun.xml files. On the other hand, if you are
        interested converting all POSCARs in a directory tree to CIFs for
        example, you will want the file paths.

        Args:
            path: input path as a tuple generated from os.walk, i.e.,
                (parent, subdirs, files).

        Returns:
            List of valid dir/file paths for assimilation
        """
        return


class VaspToComputedEntryDrone(AbstractDrone):
    """VaspToEntryDrone assimilates directories containing VASP output to
    ComputedEntry/ComputedStructureEntry objects.

    There are some restrictions on the valid directory structures:

    1. There can be only one vasp run in each directory.
    2. Directories designated "relax1", "relax2" are considered to be 2 parts
       of an aflow style run, and only "relax2" is parsed.
    3. The drone parses only the vasprun.xml file.
    """

    def __init__(
        self,
        inc_structure: bool = False,
        parameters: list[str] | None = None,
        data: dict | None = None,
    ) -> None:
        """
        Args:
            inc_structure (bool): Set to True if you want ComputedStructureEntries
                to be returned instead of ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the Vasprun object. See
                pymatgen.io.vasp.Vasprun. If parameters is None,
                a default set of parameters that are necessary for typical
                post-processing will be set.
            data (list): Output data to include. Has to be one of the properties
                supported by the Vasprun object.
        """
        self._inc_structure = inc_structure
        self._parameters = {
            "is_hubbard",
            "hubbards",
            "potcar_spec",
            "potcar_symbols",
            "run_type",
        }
        if parameters:
            self._parameters.update(parameters)
        self._data = data or {}

    def __str__(self) -> str:
        return "VaspToComputedEntryDrone"

    def assimilate(self, path: PathLike) -> ComputedStructureEntry | ComputedEntry | None:
        """Assimilate data in a directory path into a ComputedEntry object.

        Args:
            path: directory path

        Returns:
            ComputedEntry
        """
        files = os.listdir(path)
        if "relax1" in files and "relax2" in files:
            filepath = glob(f"{path}/relax2/vasprun.xml*")[0]
        else:
            vasprun_files = glob(f"{path}/vasprun.xml*")
            if len(vasprun_files) == 1:
                filepath = vasprun_files[0]
            else:
                # Since multiple files are ambiguous, we will always read
                # the last one alphabetically.
                filepath = max(vasprun_files)
                warnings.warn(f"{len(vasprun_files)} vasprun.xml.* found. {filepath} is being parsed.", stacklevel=2)

        try:
            vasp_run = Vasprun(filepath)
        except Exception as exc:
            logger.debug(f"error in {filepath}: {exc}")
            return None

        return vasp_run.get_computed_entry(self._inc_structure, parameters=list(self._parameters), data=self._data)

        # entry.parameters["history"] = _get_transformation_history(path)

    def get_valid_paths(self, path: tuple[str, str, str]) -> list[str]:
        """Check if paths contains vasprun.xml or (POSCAR+OSZICAR).

        Args:
            path: input path as a tuple generated from os.walk, i.e.,
                (parent, subdirs, files).

        Returns:
            List of valid dir/file paths for assimilation
        """
        parent, subdirs, _files = path
        if "relax1" in subdirs and "relax2" in subdirs:
            return [parent]
        if (
            (not parent.endswith("/relax1"))
            and (not parent.endswith("/relax2"))
            and (
                len(glob(f"{parent}/vasprun.xml*")) > 0
                or (len(glob(f"{parent}/POSCAR*")) > 0 and len(glob(f"{parent}/OSZICAR*")) > 0)
            )
        ):
            return [parent]
        return []

    def as_dict(self) -> dict[str, Any]:
        """Get MSONable dict."""
        return {
            "init_args": {
                "inc_structure": self._inc_structure,
                "parameters": self._parameters,
                "data": self._data,
            },
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict Representation.

        Returns:
            VaspToComputedEntryDrone
        """
        return cls(**dct["init_args"])


class SimpleVaspToComputedEntryDrone(VaspToComputedEntryDrone):
    """A simpler VaspToComputedEntryDrone. Instead of parsing vasprun.xml, it
    parses only the INCAR, POTCAR, OSZICAR and KPOINTS files, which are much
    smaller and faster to parse. However, much fewer properties are available
    compared to the standard VaspToComputedEntryDrone.
    """

    def __init__(self, inc_structure: bool = False) -> None:
        """
        Args:
            inc_structure (bool): Return ComputedStructureEntries (True) instead of
                ComputedEntries (False). Structure will be parsed from the CONTCAR.
        """
        self._inc_structure = inc_structure
        self._parameters = {"is_hubbard", "hubbards", "potcar_spec", "run_type"}

    def __str__(self) -> str:
        return "SimpleVaspToComputedEntryDrone"

    def assimilate(self, path: PathLike) -> ComputedStructureEntry | ComputedEntry | None:
        """Assimilate data in a directory path into a ComputedEntry object.

        Args:
            path: directory path

        Returns:
            ComputedEntry
        """
        files = os.listdir(path)
        try:
            files_to_parse = {}
            filenames = {"INCAR", "POTCAR", "CONTCAR", "OSZICAR", "POSCAR", "DYNMAT"}
            if "relax1" in files and "relax2" in files:
                for filename in ("INCAR", "POTCAR", "POSCAR"):
                    search_str = f"{path}/relax1/{filename}*"
                    files_to_parse[filename] = glob(search_str)[0]
                for filename in ("CONTCAR", "OSZICAR"):
                    search_str = f"{path}/relax2/{filename}*"
                    files_to_parse[filename] = glob(search_str)[-1]
            else:
                for filename in filenames:
                    files = sorted(glob(os.path.join(path, f"{filename}*")))
                    if len(files) == 0:
                        continue

                    if len(files) == 1 or filename in {"INCAR", "POTCAR"}:
                        files_to_parse[filename] = files[0]
                    elif len(files) > 1:
                        # Since multiple files are ambiguous, we will always
                        # use the first one for POSCAR and the last one
                        # alphabetically for CONTCAR and OSZICAR.

                        files_to_parse[filename] = files[0] if filename == "POSCAR" else files[-1]
                        warnings.warn(
                            f"{len(files)} files found. {files_to_parse[filename]} is being parsed.", stacklevel=2
                        )

            if not set(files_to_parse).issuperset({"INCAR", "POTCAR", "CONTCAR", "OSZICAR", "POSCAR"}):
                raise ValueError(
                    f"Unable to parse {files_to_parse} as not all necessary files are present! "
                    "SimpleVaspToComputedEntryDrone requires INCAR, POTCAR, CONTCAR, OSZICAR, POSCAR "
                    f"to be present. Only {files} detected"
                )

            poscar = Poscar.from_file(files_to_parse["POSCAR"])
            contcar = Poscar.from_file(files_to_parse["CONTCAR"])
            incar = Incar.from_file(files_to_parse["INCAR"])
            potcar = Potcar.from_file(files_to_parse["POTCAR"])
            oszicar = Oszicar(files_to_parse["OSZICAR"])

            param: dict[str, Any] = {"hubbards": {}}
            if "LDAUU" in incar:
                param["hubbards"] = dict(zip(poscar.site_symbols, incar["LDAUU"], strict=True))
            param["is_hubbard"] = incar.get("LDAU", True) and sum(param["hubbards"].values()) > 0
            param["run_type"] = None
            param["potcar_spec"] = potcar.spec
            energy = oszicar.final_energy
            structure = contcar.structure
            initial_vol = poscar.structure.volume
            final_vol = contcar.structure.volume
            delta_volume = final_vol / initial_vol - 1
            data: dict[str, Any] = {"filename": path, "delta_volume": delta_volume}
            if "DYNMAT" in files_to_parse:
                dynmat = Dynmat(files_to_parse["DYNMAT"])
                data["phonon_frequencies"] = dynmat.get_phonon_frequencies()
            if self._inc_structure:
                return ComputedStructureEntry(structure, energy, parameters=param, data=data)
            return ComputedEntry(structure.composition, energy, parameters=param, data=data)

        except Exception as exc:
            if isinstance(exc, ValueError) and "not all necessary files are present" in str(exc):
                raise

            logger.debug(f"error in {path}: {exc}")
            return None

    def as_dict(self) -> dict[str, Any]:
        """Get MSONable dict."""
        return {
            "init_args": {"inc_structure": self._inc_structure},
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict Representation.

        Returns:
            SimpleVaspToComputedEntryDrone
        """
        return cls(**dct["init_args"])


class GaussianToComputedEntryDrone(AbstractDrone):
    """GaussianToEntryDrone assimilates directories containing Gaussian output to
    ComputedEntry/ComputedStructureEntry objects. By default, it is assumed
    that Gaussian output files have a ".log" extension.

    Note:
        Like the GaussianOutput class, this is still in early beta.
    """

    def __init__(
        self,
        inc_structure: bool = False,
        parameters: list[str] | None = None,
        data: list[str] | None = None,
        file_extensions: Sequence[str] = (".log",),
    ) -> None:
        """
        Args:
            inc_structure (bool): Set to True if you want
                ComputedStructureEntries to be returned instead of
                ComputedEntries.
            parameters (list): Input parameters to include. It has to be one of
                the properties supported by the GaussianOutput object. See
                pymatgen.io.gaussian.GaussianOutput. The parameters
                have to be one of python's primitive types, i.e., list, dict of
                strings and integers. If parameters is None, a default set of
                parameters will be set.
            data (list): Output data to include. Has to be one of the properties
                supported by the GaussianOutput object. The parameters have to
                be one of python's primitive types, i.e. list, dict of strings
                and integers. If data is None, a default set will be set.
            file_extensions (list):
                File extensions to be considered as Gaussian output files.
                Defaults to just the typical "log" extension.
        """
        self._inc_structure = inc_structure
        self._parameters = {
            "functional",
            "basis_set",
            "charge",
            "spin_multiplicity",
            "route_parameters",
        }

        if parameters:
            self._parameters.update(parameters)

        self._data = {"stationary_type", "properly_terminated"}
        if data:
            self._data.update(data)

        self._file_extensions = file_extensions

    def __str__(self) -> str:
        return "GaussianToComputedEntryDrone"

    def assimilate(self, path: PathLike) -> ComputedStructureEntry | ComputedEntry | None:
        """Assimilate data in a directory path into a ComputedEntry object.

        Args:
            path: directory path

        Returns:
            ComputedEntry
        """
        try:
            gau_run = GaussianOutput(path)
        except Exception as exc:
            logger.debug(f"error in {path}: {exc}")
            return None

        param = {p: getattr(gau_run, p) for p in self._parameters}
        data = {d: getattr(gau_run, d) for d in self._data}

        if self._inc_structure:
            return ComputedStructureEntry(
                gau_run.final_structure,
                gau_run.final_energy,
                parameters=param,
                data=data,
            )
        return ComputedEntry(
            gau_run.final_structure.composition,
            gau_run.final_energy,
            parameters=param,
            data=data,
        )

    def get_valid_paths(self, path: tuple[str, str, str]) -> list[str]:
        """Check if path contains files with define extensions.

        Args:
            path: input path as a tuple generated from os.walk, i.e.,
                (parent, subdirs, files).

        Returns:
            List of valid dir/file paths for assimilation
        """
        parent, _subdirs, files = path
        return [os.path.join(parent, file) for file in files if os.path.splitext(file)[1] in self._file_extensions]

    def as_dict(self) -> dict[str, Any]:
        """Get MSONable dict."""
        return {
            "init_args": {
                "inc_structure": self._inc_structure,
                "parameters": self._parameters,
                "data": self._data,
                "file_extensions": self._file_extensions,
            },
            "@module": type(self).__module__,
            "@class": type(self).__name__,
        }

    @classmethod
    def from_dict(cls, dct: dict[str, Any]) -> Self:
        """
        Args:
            dct (dict): Dict Representation.

        Returns:
            GaussianToComputedEntryDrone
        """
        return cls(**dct["init_args"])


def _get_transformation_history(path: PathLike):
    """Check for a transformations.json* file and return the history."""
    if trans_json := glob(f"{path!s}/transformations.json*"):
        try:
            with zopen(trans_json[0], "rb") as file:
                return orjson.loads(file.read())["history"]
        except Exception:
            return None
    return None
