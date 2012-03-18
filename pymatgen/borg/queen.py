#!/usr/bin/env python

'''
Created on Mar 18, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"


import os

from multiprocessing import Manager, Pool


class BorgQueen(object):
    """
    The Borg Queen controls the drones to assimilate data in an entire directory
    substructure. Uses multiprocessing to speed up things considerably.
    """

    def __init__(self, rootpath, drone, number_of_drones = 1):
        """
        Args:
            rootpath:
                The root directory to start assimilation.
            drones:
                The drones to use for assimilation
        """
        valid_paths = []
        for (parent, subdirs, files) in os.walk(rootpath):
            if drone.is_valid_path((parent, subdirs, files)):
                valid_paths.append(parent)
        manager = Manager()
        data = manager.list()

        p = Pool(number_of_drones)
        p.map(order_assimilation, ((path, drone, data) for path in valid_paths))
        self._data = [drone.convert(d) for d in data]

    def get_data(self):
        return self._data

def order_assimilation(args):
    (path, drone, data) = args
    newdata = drone.assimilate(path)
    if newdata:
        data.append(newdata)
