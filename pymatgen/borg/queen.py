#!/usr/bin/env python

'''
This module defines the BorgQueen class, which manages drones to assimilate 
data using Python's multiprocessing. 
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"


import os
import json

from pymatgen.util.io_utils import file_open_zip_aware

from multiprocessing import Manager, Pool


class BorgQueen(object):
    """
    The Borg Queen controls the drones to assimilate data in an entire directory
    substructure. Uses multiprocessing to speed up things considerably. It also
    contains convenience methods to save and load data between sessions.
    """

    def __init__(self, drone, rootpath=None, number_of_drones=1):
        """
        Args:
            drone:
                The drone to use for assimilation
            rootpath:
                The root directory to start assimilation. Leave it as None if 
                you want to do assimilation later, or is using the BorgQueen
                to load previously assimilated data.
            number_of_drones:
                Number of drones to parallelize over. Typical machines today
                have up to four processors. Note that you won't see a 100%
                improvement with two drones over one, but you will definitely
                see a significant speedup of at least 50% or so. If you are
                running this over a server with far more processors, the
                speedup will be even greater.
        """
        self._drone = drone
        self._num_drones = number_of_drones
        self._data = []
        if rootpath:
            self.parallel_assimilate(rootpath)

    def parallel_assimilate(self, rootpath):
        """
        Assimilate the entire subdirectory structure in rootpath.
        """
        valid_paths = []
        for (parent, subdirs, files) in os.walk(rootpath):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs, files)))

        manager = Manager()
        data = manager.list()

        p = Pool(self._num_drones)
        p.map(order_assimilation, ((path, self._drone, data) for path in valid_paths))
        self._data.extend(data)

    def serial_assimilate(self, rootpath):
        """
        Assimilate the entire subdirectory structure in rootpath serially.
        """
        valid_paths = []
        for (parent, subdirs, files) in os.walk(rootpath):
            valid_paths.extend(self._drone.get_valid_paths((parent, subdirs, files)))
        data = []
        for path in valid_paths:
            order_assimilation((path, self._drone, data))
        self._data.extend(data)

    def get_assimilated_data(self):
        """
        Returns an list of assimilated objects
        """
        return [self._drone.convert(d) for d in self._data]

    def save_data(self, filename):
        """
        Save the assimilated data to a file.
        
        Args:
            filename:
                filename to save the assimilated data to. Note that if the
                filename ends with gz or bz2, the relevant gzip or bz2
                compression will be applied.
        """
        with file_open_zip_aware(filename, "w") as f:
            json.dump(list(self._data), f)

    def load_data(self, filename):
        """
        Load assimilated data from a file
        """
        with file_open_zip_aware(filename, "r") as f:
            self._data = json.load(f)


def order_assimilation(args):
    (path, drone, data) = args
    newdata = drone.assimilate(path)
    if newdata:
        data.append(newdata)
