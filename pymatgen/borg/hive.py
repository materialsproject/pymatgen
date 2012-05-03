#!/usr/bin/env python

'''
This module define the various drones used to assimilate data.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Mar 18, 2012"

import abc
import os
import re
import glob
import logging

from pymatgen.io.vaspio import Vasprun
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry

logger = logging.getLogger(__name__)


class AbstractDrone(object):
    """
    Abstract drone class that defines the various methods that must be implemented
    by drones. Because of the quirky nature of Python's multiprocessing, the 
    representations has to be in the form of python primitives. This can then
    be reverted to the original object with drone.convert.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def assimilate(self, path):
        '''
        Assimilate data in a directory path into a pymatgen object. Because of
        the quirky nature of Python's multiprocessing, the object must support
        pymatgen's to_dict for parallel processing.
        
        Args:
            path:
                directory path
                
        Returns:
            An assimilated object
        '''
        return

    @abc.abstractmethod
    def get_valid_paths(self, path):
        """
        Checks if path contains valid data for assimilation, and then returns
        the valid paths. The paths returned can be a list of directory or file
        paths, depending on what kind of data you are assimilating. For example,
        if you are assimilating VASP runs, you are only interested in
        directories containing vasprun.xml files. On the other hand, if you are
        interested converting all POSCARs in a directory tree to cifs for
        example, you will want the file paths.
        
        Args:
            path:
                input path as a tuple generated from os.walk, i.e.,
                (parent, subdirs, files).
                
        Returns:
            List of valid dir/file paths for assimilation
        """
        return

    @abc.abstractproperty
    def to_dict(self):
        """
        All drones must support a json serializable dict representation
        """
        return


class VaspToComputedEntryDrone(AbstractDrone):
    """
    VaspToEntryDrone assimilates directories containing vasp input to 
    ComputedEntry/ComputedStructureEntry objects. There are some restrictions
    on the valid directory structures:
    
    1. There can be only one vasp run in each directory.
    2. Directories designated "relax1", "relax2" are considered to be 2 parts of
       an aflow style run, and only "relax2" is parsed.
       
    """

    def __init__(self, inc_structure=False, parameters=None, data=None):
        """
        Args:
            inc_structure:
                Set to True if you want ComputedStructureEntries to be returned
                instead of ComputedEntries.
            parameters:
                Input parameters to include. It has to be one of the properties
                supported by the Vasprun object. See pymatgen.io.vaspio Vasprun.
                The parameters have to be one of python's primitive types,
                i.e. list, dict of strings and integers. Complex objects such as
                dos are not supported at this point.
                If parameters == None, a default set of parameters that are 
                necessary for typical post-processing will be set.
            data:
                Output data to include. Has to be one of the properties
                supported by the Vasprun object. The parameters have to be one
                of python's primitive types, i.e. list, dict of strings and
                integers. Complex objects such as dos are not supported at this
                point. e.g., ['filename']
        """
        self._inc_structure = inc_structure
        self._parameters = parameters if parameters else ["is_hubbard", "hubbards", "potcar_symbols", "run_type"]
        self._data = data if data else []

    def assimilate(self, path):
        files = os.listdir(path)
        if 'relax1' in files and 'relax2' in files:
            filepath = glob.glob(os.path.join(path, "relax2", "vasprun.xml*"))[0]
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
                it will try to get a relax2 from an aflow style run if possible.
                Or else, a randomly chosen file containing vasprun.xml is chosen.
                """
                for fname in vasprun_files:
                    if os.path.basename(fname) in ["vasprun.xml", "vasprun.xml.gz", "vasprun.xml.bz2"]:
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
        param = {}
        for p in self._parameters:
            param[p] = getattr(vasprun, p)
        data = {}
        for d in self._data:
            data[d] = getattr(vasprun, d)
        if self._inc_structure:
            entry = ComputedStructureEntry(vasprun.final_structure,
                                   vasprun.final_energy, parameters=param, data=data)
        else:
            entry = ComputedEntry(vasprun.final_structure.composition,
                                   vasprun.final_energy, parameters=param, data=data)
        return entry

    def get_valid_paths(self, path):
        (parent, subdirs, files) = path
        if 'relax1' in subdirs and 'relax2' in subdirs:
            return [parent]
        if (not parent.endswith('/relax1')) and (not parent.endswith('/relax2')) and len(glob.glob(os.path.join(parent, "vasprun.xml*"))) > 0:
            return [parent]
        return []

    def __str__(self):
        return " VaspToComputedEntryDrone"

    @property
    def to_dict(self):
        init_args = {'inc_structure' : self._inc_structure,
                     "parameters": self._parameters,
                     "data": self._data}
        d = {'init_args': init_args, 'version': __version__ }
        d['module'] = self.__class__.__module__
        d['class'] = self.__class__.__name__
        return d

    @staticmethod
    def from_dict(d):
        return VaspToComputedEntryDrone(**d['init_args'])
