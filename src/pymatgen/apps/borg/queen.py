"""
This module defines the BorgQueen class, which manages drones to assimilate
data using Python's multiprocessing.
"""

from __future__ import annotations

import json
import logging
import os
from multiprocessing import Manager, Pool
from typing import TYPE_CHECKING

from monty.io import zopen
from monty.json import MontyDecoder, MontyEncoder

if TYPE_CHECKING:
    from pymatgen.apps.borg.hive import AbstractDrone
    from pymatgen.util.typing import PathLike

logger = logging.getLogger("BorgQueen")


class BorgQueen:
    """The Borg Queen controls the drones to assimilate data in an entire
    directory tree. Uses multiprocessing to speed up things considerably. It
    also contains convenience methods to save and load data between sessions.
    """

    def __init__(
        self,
        drone: AbstractDrone,
        rootpath: PathLike | None = None,
        number_of_drones: int = 1,
    ) -> None:
        """
        Args:
            drone (AbstractDrone): An implementation of
                AbstractDrone to use for assimilation.
            rootpath (PathLike): The root directory to start assimilation. Leave it
                as None if you want to do assimilation later, or is using the
                BorgQueen to load previously assimilated data.
            number_of_drones (int): Number of drones to parallelize over.
                Typical machines today have up to four processors. Note that you
                won't see a 100% improvement with two drones over one, but you
                will definitely see a significant speedup of at least 50% or so.
                If you are running this over a server with far more processors,
                the speedup will be even greater.
        """
        self._drone = drone
        self._num_drones = number_of_drones
        self._data: list = []

        if rootpath:
            if number_of_drones > 1:
                self.parallel_assimilate(rootpath)
            else:
                self.serial_assimilate(rootpath)

    def parallel_assimilate(self, rootpath: PathLike) -> None:
        """Assimilate the entire subdirectory structure in rootpath."""
        logger.info("Scanning for valid paths...")
        valid_paths = []
        for parent, subdirs, files in os.walk(rootpath):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs, files)))
        manager = Manager()
        data = manager.list()
        status = manager.dict()
        status["count"] = 0
        status["total"] = len(valid_paths)
        logger.info(f"{len(valid_paths)} valid paths found.")
        with Pool(self._num_drones) as pool:
            pool.map(
                order_assimilation,
                ((path, self._drone, data, status) for path in valid_paths),
            )
            for string in data:
                self._data.append(json.loads(string, cls=MontyDecoder))

    def serial_assimilate(self, root: PathLike) -> None:
        """Assimilate the entire subdirectory structure in rootpath serially."""
        valid_paths = []
        for parent, subdirs, files in os.walk(root):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs, files)))
        data: list[str] = []
        total = len(valid_paths)
        for idx, path in enumerate(valid_paths, start=1):
            new_data = self._drone.assimilate(path)
            self._data.append(new_data)
            logger.info(f"{idx}/{total} ({idx / total:.1%}) done")
        for json_str in data:
            self._data.append(json.loads(json_str, cls=MontyDecoder))

    def get_data(self) -> list:
        """Get an list of assimilated objects."""
        return self._data

    def save_data(self, filename: PathLike) -> None:
        """Save the assimilated data to a file.

        Args:
            filename (str): filename to save the assimilated data to. Note
                that if the filename ends with gz or bz2, the relevant gzip
                or bz2 compression will be applied.
        """
        with zopen(filename, mode="wt", encoding="utf-8") as file:
            json.dump(list(self._data), file, cls=MontyEncoder)  # type:ignore[arg-type]

    def load_data(self, filename: PathLike) -> None:
        """Load assimilated data from a file."""
        with zopen(filename, mode="rt", encoding="utf-8") as file:
            self._data = json.load(file, cls=MontyDecoder)


def order_assimilation(args: tuple) -> None:
    """Internal helper method for BorgQueen to process assimilation."""
    path, drone, data, status = args
    if new_data := drone.assimilate(path):
        data.append(json.dumps(new_data, cls=MontyEncoder))
    status["count"] += 1
    count = status["count"]
    total = status["total"]
    logger.info(f"{count}/{total} ({count / total:.2%}) done")
