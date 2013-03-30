from __future__ import division, print_function

import sys
import os
import os.path
import collections
import abc
import shutil
import numpy as np

from pprint import pprint

from pymatgen.util.num_utils import sort_dict
from pymatgen.io.abinitio.pseudos import Pseudo
from pymatgen.io.abinitio.workflow import PPConvergenceFactory
from pymatgen.io.abinitio.deltaworks import DeltaFactory
from pymatgen.serializers.json_coders import MSONable, json_pretty_dump #, PMGJSONDecoder

##########################################################################################

class AttrDict(dict):
    "Access dict keys as obj.foo instead of obj['foo']"
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

##########################################################################################

class DojoError(Exception):
    pass

class Dojo(object):
    Error = DojoError

    def __init__(self):
        # List of master classes that will be instanciated afterwards.
        # They are ordered according to the master level.
        self.master_classes = [m for m in DojoMaster.__subclasses__()]
        self.master_classes.sort(key = lambda cls : cls.dojo_level)

    def accept_pseudos(self, pseudos):

        if isinstance(pseudos, str) or not isinstance(pseudos, collections.Iterable):
            pseudos = [pseudos,]

        self.pseudos = []
        for pseudo in pseudos:
            if not isinstance(pseudo, Pseudo):
                pseudo = Pseudo.from_filename(pseudo)
            self.pseudos.append(pseudo)

    def start_training(self, **kwargs):

        for pseudo in self.pseudos:
            workdir = "DOJO_" + pseudo.name

            for cls in self.master_classes:
                master = cls()
                if master.accept_pseudo(pseudo):
                    master.start_training(workdir, **kwargs)

    #def rank_pseudos(self):
    #    pass

##########################################################################################

class DojoMaster(object):
    "Subclasses must define the class attribute level"
    __metaclass__ = abc.ABCMeta

    Error = DojoError

    def __init__(self):
        self.reports = []
        self.errors = []

    @staticmethod
    def from_dojo_level(level):
        "Static constructor, returns an instance of DojoMaster for the the given level"
        classes = []
        for cls in DojoMaster.__subclasses__():
            if cls.dojo_level == dojo_level:
                classes.append(cls)

        if len(classes) != 1:
            raise self.Error("Found %d drones for dojo_level %d" % (len(classes), dojo_level))
                                                                                              
        return classes[0]()

    def accept_pseudo(self, pseudo):
        """
        This method is called before testing the pseudo
        Return True if self can train the pseudo. 
        A master can train pseudo if his level == pseudo.dojo_level + 1"
        """
        if not isinstance(pseudo, Pseudo):
            pseudo = Pseudo.from_filename(pseudo)

        ready = False
        if pseudo.dojo_level is None:
            # hints are missing
            ready = (self.dojo_level == 0)
        else:
            ready = (pseudo.dojo_level == self.dojo_level - 1)

        if not ready:
            msg = "%s: Sorry, %s-san, I cannot train you" % (self.__class__.__name__, pseudo.name)
            print(msg)
        else:
            print("%s: Welcome %s-san, I'm your level-%d trainer" % (self.__class__.__name__, pseudo.name, self.dojo_level))
            self.pseudo = pseudo

        return ready

    @abc.abstractmethod
    def challenge(self, workdir, **kwargs):
        "Run the calculation"

    @abc.abstractmethod
    def make_report(self, **kwargs):
        "Returns report, isok"

    def write_dojo_report(self, report, overwrite_data=False, ignore_errors=False):
        "Writes/updates the DOJO_REPORT section of the pseudopotential"

        #if self.errors and not ignore_errors:
        #    pprint(self.errors)
        #    raise self.Error("Cannot update dojo data since self.errors is not empty")
        pseudo = self.pseudo

        # Read old_report from pseudo.
        old_report = pseudo.read_dojo_report()

        for okey in old_report:
            if okey in report and not overwrite:
                raise self.Error("%k already exists in the old pseudo. Cannot overwrite data" % okey)

        # Update the report card with the input report
        old_report.update(report)

        # Write new report
        pseudo.write_dojo_report(old_report)

    def start_training(self, workdir, **kwargs):

        results = self.challenge(workdir, **kwargs)

        report, isok = self.make_report(results, **kwargs)

        if isok:
            self.write_dojo_report(report)
        else:
            raise self.Error("isok: %s" % isok)

##########################################################################################

class HintsMaster(DojoMaster):
    dojo_level = 0
    dojo_key = "hints"

    def challenge(self, workdir, **kwargs):
        pseudo = self.pseudo

        #tols_mev TODO
        factory = PPConvergenceFactory()

        workdir = os.path.join(workdir, "LEVEL_" + str(self.dojo_level))

        estep = 10 
        eslice = slice(5,None,estep)
        print("Converging %s in iterative mode with eslice %s" % (pseudo.name, eslice))

        w = factory.work_for_pseudo(workdir, pseudo, eslice)
        #print(w)

        if os.path.exists(w.workdir):
            shutil.rmtree(w.workdir)

        w.start()
        w.wait()

        wres = w.get_results()
        w.rmtree()

        #estart, estop, estep = max(wres["ecut_low"] - estep, 5), 10, 1
        estart, estop, estep = max(wres["low"]["ecut"] - estep, 5), wres["high"]["ecut"] + estep, 5

        erange = list(np.arange(estart, estop, estep))

        print("Finding optimal values for ecut in the interval %.1f %.1f %1.f" % (estart, estop, estep))
        w = factory.work_for_pseudo(workdir, pseudo, erange)
        print(w)

        w.start()
        w.wait()

        wres = w.get_results()

        json_pretty_dump(wres, w.path_in_workdir("results.json"))

        return wres

    def make_report(self, results, **kwargs):
        isok = True

        d = {}
        for key in ["low", "normal", "high"]:
            try:
                d[key] = results[key]
            except KeyError:
                raise KeyError("%s is missing in input results" % key)

        return {self.dojo_key : d}, isok

##########################################################################################

class DeltaFactorMaster(DojoMaster):
    dojo_level = 1
    dojo_key = "delta_factor"

    def challenge(self, workdir, **kwargs):
        pseudo = self.pseudo

        factory = DeltaFactory()

        workdir = os.path.join(workdir, "LEVEL_" + str(self.dojo_level))

        w = factory.work_for_pseudo(workdir, pseudo, kppa=1)

        w.start()
        w.wait()
                                                                        
        wres = w.get_results()
                                                                        
        json_pretty_dump(wres, w.path_in_workdir("results.json"))

        return wres

    def make_report(self, results, **kwargs):
        isok = True

        d = { 
            #"delta_factor" : ,
            #"e0" :
            #"v0" :
            #"b"  :
            #"bp" :
            #"perr_v0" :
            #"perr_b"  :
            #"perr_bp" :
        }

        #for key in keys:
        #    try:
        #        d[key] = results[key]
        #    except KeyError:
        #        raise KeyError("%s is missing in input results" % key)

        return {self.dojo_key : d}, isok

##########################################################################################

_key2lev = {}
for cls in DojoMaster.__subclasses__():
    _key2lev[cls.dojo_key] = cls.dojo_level

def dojo_key2level(key):
    "Return the trial level from the name found in the pseudo"
    return _key2lev[key] 

##########################################################################################

