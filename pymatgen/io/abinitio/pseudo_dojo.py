from __future__ import division, print_function

import sys
import os
import os.path
import collections
import abc
import json
import shutil
import numpy as np

from pprint import pprint

from pymatgen.util.num_utils import sort_dict
from pymatgen.io.abinitio.workflow import PPConvergenceFactory
from pymatgen.io.abinitio.pseudos import Pseudo

##########################################################################################

class AttrDict(dict):
    "Access dict keys as obj.foo instead of obj['foo']"
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

##########################################################################################
class DojoError(Exception):
    pass

class DojoDrone(object):
    """
    Abstract drone class that defines the various methods that must be implemented by drones. 
    """
    __metaclass__ = abc.ABCMeta

    Error = DojoError

    def __init__(self):
        self.results = {}
        self.errors = {}

    @staticmethod
    def from_dojo_level(dojo_level):
        classes = []
        for cls in DojoDrone.__subclasses__():
            if cls.dojo_level == dojo_level:
                classes.append(cls)
        if len(classes) != 1:
            raise DojoError("Found %d drones for dojo_level %d" % (len(classes), dojo_level))

        return classes[0]()

    def walk(self, top="."):
        "Scan dirs/files starting from top"
        r = self.results
        e = self.errors

        for (dirpath, dirnames, filenames) in os.walk(top):
            for fname in filenames:

                if fname == "ppdojo_results.json":
                    path = os.path.join(dirpath, fname)

                    data = self.assimilate_file(path)

                    pp_name = data["pseudo_name"]
                    pp_path = data["pseudo_path"]

                    if pp_path in r:
                        raise ValueError("%s has been already analyzed" % pp_path)

                    try:
                        r[pp_path] = self.analyze(data)
                    except Exception as exc:
                        e[pp_path] = exc

    def assimilate_file(self, filepath):
        """
        Assimilates data in a file 

        Returns:
            Dictionary with the assimilated results.
        """
        with open(filepath, "r") as fh:
            return json.load(fh)

    #@abc.abstractmethod
    #def get_valid_paths(self, path):
    #    """
    #    Checks if path contains valid data for assimilation, and then returns
    #    the valid paths. The paths returned can be a list of directory or file
    #    paths, depending on what kind of data you are assimilating. 

    #    For example, if you are assimilating VASP runs, you are only interested in
    #    directories containing vasprun.xml files. On the other hand, if you are
    #    interested converting all POSCARs in a directory tree to cifs for
    #    example, you will want the file paths.
    #    Args:
    #        path:
    #            input path as a tuple generated from os.walk, i.e., (parent, subdirs, files).
    #    Returns:
    #        List of valid dir/file paths for assimilation
    #    """
    #    return

    @abc.abstractmethod
    def analyze(self, data):
        "Analyze the data stored in a dictionary. Returns a new dict with postprocessed results"

    def write_dojo_data(self, overwrite_data=False, ignore_errors=False):
        "Writes/updates the dojo_data section of the pseudopotentials"

        if drone.errors and not ignore_errors:
            pprint(drone.errors)
            raise self.Error("Cannot update dojo data since self.errors is not empty")

        for (pp_path, res) in self.results.items():
            new_data = self.results[pp_path]

            pseudo = Pseudo.from_filename(pp_path)

            #old_data = pseudo.read_dojo_data()

            #for okey in old_data:
            #    if okey in new_data and not overwrite:
            #        raise RuntimeError("%k already exists in the old pseudo. Cannot overwrite data" % okey)

            #old_data.update(new_data)

            #pseudo.write_dojo_data(old_data)
                                                                                      
##########################################################################################

class HintsDrone(DojoDrone):
    dojo_level = 0

    def analyze(self, data):
        #pprint(data)

        #etotal    = data["etotal"]      
        #ecut_list = data["ecut_list"]   
        #aug_ratio = data["aug_ratio"]   

        #"pseudo_name" = data["pseudo_name"] 
        #"pseudo_path" = data["pseudo_path"]
        #"ecut_low"    = data["ecut_low"]
        #"ecut_normal" = data["ecut_normal"]
        #"ecut_high"   = data["ecut_high"]

        #de_high, de_normal, de_low = 0.05e-3/Ha_eV,  1e-3/Ha_eV, 10e-3/Ha_eV
        #
        #ihigh   = check_conv(etotal, de_high, min_numpts=1)
        #inormal = check_conv(etotal, de_normal)
        #ilow    = check_conv(etotal, de_low)

        #print("ihigh %d, inormal %d, ilow %d" % (ihigh, inormal, ilow))

        #ecut_high, ecut_normal, ecut_low = 3 * (None,)

        #ihigh != -1:
        #    ecut_low    = ecut_list[ilow] 
        #    ecut_normal = ecut_list[inormal] 
        #    ecut_high   = ecut_list[ihigh] 

        data = AttrDict(data)

        od = AttrDict()
        od.hint_low = {"ecut": data.ecut_low}
        od.hint_normal = {"ecut": data.ecut_normal}
        od.hint_high = {"ecut": data.ecut_high}

        with open("test.json", "w") as fh:
            json.dump(od, fh, indent=4, sort_keys=True)

        return od

##########################################################################################

class DeltaFactorDrone(DojoDrone):
    dojo_level = 1

    def analyze(self, data):

        data = AttrDict(data)
        od = AttrDict()

        v0 = data["v0"]
        b0 = data["b0"]
        bp = data["bp"]
        symbol = data["symbol"]

        #ref_data = DeltaData.
        # Compare ref with data and compute the relative error.

        return od

##########################################################################################

class Dojo(object):
    Error = DojoError

    levels = [1,]

    def accept_pseudos(self, pseudos):

        if isinstance(pseudos, str) or not isinstance(pseudos, collections.Iterable):
            pseudos = [pseudos,]

        self.pseudos = []
        for pseudo in pseudos:
            if not isinstance(pseudo, Pseudo):
                pseudo = Pseudo.from_filename(pseudo)
            self.pseudos.append(pseudo)

        self.masters = [DojoMaster() for i in range(len(self.pseudos))]

        for (master, pseudo) in zip(self.masters, self.pseudos):
            master.accept_pseudo(pseudo)

    def start_training(self, *args, **kwargs):
        for master in self.masters:
            master.start_training()

    def rank_pseudos(self):
        pass

##########################################################################################

class DojoMaster(object):
    Error = DojoError

    def accept_pseudo(self, pseudo):

        if not isinstance(pseudo, Pseudo):
            pseudo = Pseudo.from_filename(pseudo)

        print("Welcome %s-san, I will be your master" % pseudo.name)
        self.pseudo = pseudo

    def start_training(self, *args, **kwargs):
        # Compute hints for ecut, if pseudo doesn't have it.
        #if pseudo.dojo_level is None:
        res = trial_0(self.pseudo)
        pprint(res)

        #drone = DojoDrone.from_dojo_level(0)
        #drone.walk()
        for level in Dojo.levels:
            res = self.challenge_pseudo(self.pseudo, level, *args, **kwargs)

    def challenge_pseudo(self, pseudo, level, *args, **kwargs):
        res = {}
        return res

    def rank_pseudos(self):
        pass

##########################################################################################

def trial_0(pseudo, *args, **kwargs):
    factory = PPConvergenceFactory()

    dirpref = "DOJO0"

    estep = 10 
    eslice = slice(5,None,estep)
    print("Converging %s in iterative mode with eslice %s" % (pseudo.name, eslice))

    w = factory.work_for_pseudo(pseudo, eslice, dirpref=dirpref)
    #print(w)

    if os.path.exists(w.workdir):
        shutil.rmtree(w.workdir)

    w.start()
    w.wait()

    wres = w.get_results()
    w.rmtree()

    #estart, estop, estep = max(wres["ecut_low"] - estep, 5), 10, 1
    estart, estop, estep = max(wres["ecut_low"] - estep, 5), wres["ecut_high"] + estep, 5

    erange = list(np.arange(estart, estop, estep))

    print("Finding optimal values for ecut in the interval %.1f %.1f %1.f" % (estart, estop, estep))
    w = factory.work_for_pseudo(pseudo, erange, dirpref=dirpref)
    print(w)

    w.start()
    w.wait()

    wres = w.get_results()

    #pprint(wres)
    #import json
    #with open(w.path_in_workdir("ppdojo_results.json"), "w") as fh:
    #    json.dump(wres, fh, indent=4, sort_keys=True)
    return wres

##########################################################################################
