# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module defines the BorgQueen class, which manages drones to assimilate
data using Python's multiprocessing.
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 18, 2012"


import os
import json
import logging

from monty.io import zopen
from monty.json import MontyEncoder, MontyDecoder

from multiprocessing import Manager, Pool

logger = logging.getLogger("BorgQueen")


class BorgQueen:
    """
    The Borg Queen controls the drones to assimilate data in an entire
    directory tree. Uses multiprocessing to speed up things considerably. It
    also contains convenience methods to save and load data between sessions.
    """

    def __init__(self, drone, rootpath=None, number_of_drones=1):
        """
        Args:
            drone (Drone): An implementation of
                :class:`pymatgen.apps.borg.hive.AbstractDrone` to use for
                assimilation.
            rootpath (str): The root directory to start assimilation. Leave it
                as None if you want to do assimilation later, or is using the
                BorgQueen to load previously assimilated data.
            ndrones (int): Number of drones to parallelize over.
                Typical machines today have up to four processors. Note that you
                won't see a 100% improvement with two drones over one, but you
                will definitely see a significant speedup of at least 50% or so.
                If you are running this over a server with far more processors,
                the speedup will be even greater.
        """
        self._drone = drone
        self._num_drones = number_of_drones
        self._data = []

        if rootpath:
            if number_of_drones > 1:
                self.parallel_assimilate(rootpath)
            else:
                self.serial_assimilate(rootpath)

    def parallel_assimilate(self, rootpath):
        """
        Assimilate the entire subdirectory structure in rootpath.
        """
        logger.info('Scanning for valid paths...')
        valid_paths = []
        for (parent, subdirs, files) in os.walk(rootpath):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs,
                                                            files)))
        manager = Manager()
        data = manager.list()
        status = manager.dict()
        status['count'] = 0
        status['total'] = len(valid_paths)
        logger.info('{} valid paths found.'.format(len(valid_paths)))
        p = Pool(self._num_drones)
        p.map(order_assimilation, ((path, self._drone, data, status)
                                   for path in valid_paths))
        for d in data:
            self._data.append(json.loads(d, cls=MontyDecoder))

    def serial_assimilate(self, rootpath):
        """
        Assimilate the entire subdirectory structure in rootpath serially.
        """
        valid_paths = []
        for (parent, subdirs, files) in os.walk(rootpath):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs,
                                                            files)))
        data = []
        count = 0
        total = len(valid_paths)
        for path in valid_paths:
            newdata = self._drone.assimilate(path)
            self._data.append(newdata)
            count += 1
            logger.info('{}/{} ({:.2f}%) done'.format(count, total,
                                                      count / total * 100))
        for d in data:
            self._data.append(json.loads(d, cls=MontyDecoder))

    def get_data(self):
        """
        Returns an list of assimilated objects
        """
        return self._data

    def save_data(self, filename):
        """
        Save the assimilated data to a file.

        Args:
            filename (str): filename to save the assimilated data to. Note
                that if the filename ends with gz or bz2, the relevant gzip
                or bz2 compression will be applied.
        """
        with zopen(filename, "wt") as f:
            json.dump(list(self._data), f, cls=MontyEncoder)

    def load_data(self, filename):
        """
        Load assimilated data from a file
        """
        with zopen(filename, "rt") as f:
            self._data = json.load(f, cls=MontyDecoder)


def order_assimilation(args):
    """
    Internal helper method for BorgQueen to process assimilation
    """
    (path, drone, data, status) = args
    newdata = drone.assimilate(path)
    if newdata:
        data.append(json.dumps(newdata, cls=MontyEncoder))
    status['count'] += 1
    count = status['count']
    total = status['total']
    logger.info('{}/{} ({:.2f}%) done'.format(count, total,
                                              count / total * 100))
