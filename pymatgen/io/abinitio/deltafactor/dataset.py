#!/usr/bin/env python
from __future__ import division, print_function

import sys
import os
import os.path
import collections

##########################################################################################

class Record(collections.namedtuple("Record", "symbol v0 b0 bp")):
    "Record with the volume, bulk-modulus and Bp for given symbol (a.u.)"

    def __new__(cls, *args):
        "Extends the base class adding type conversion of arguments."
        new_args = len(args) * [None]

        for (i, arg) in enumerate(args):
            converter = float
            if i == 0: converter = str
            new_args[i] = converter(arg)

        return super(cls, Record).__new__(cls, *new_args)

def read_data_from_filename(filename):
    """
    Reads data (v0, b0, bp) from file filename
    Returns a dict of record objects indexed by element symbol.
    """
    data = collections.OrderedDict()

    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#") or not line: 
                continue
            tokens = line.split()
            symbol = tokens[0]
            data[symbol] = Record(*tokens)
    return data

##########################################################################################

def singleton(cls):
    """
    This decorator can be used to create a singleton out of a class.
    """
    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance

@singleton
class DeltaFactorDataset(object):
    #_refcode = "WIEN2K"

    def __init__(self):
        self.dirpath = os.path.abspath(os.path.dirname(__file__))

        self._data = d = {}
        for entry in os.listdir(self.dirpath):
           file_path = os.path.join(self.dirpath, entry)
           if os.path.isfile(file_path) and file_path.endswith(".txt"):
                codename, ext = os.path.splitext(entry)
                if codename == "README": 
                    continue
                #print(codename)
                d[codename] = read_data_from_filename(file_path)

        self.cif_paths = d = {}

        cif_dirpath = os.path.join(self.dirpath, "CIFs")
        for entry in os.listdir(cif_dirpath):
            if entry.endswith(".cif"):
                symbol, ext = os.path.splitext(entry)
                d[symbol] = os.path.join(cif_dirpath, entry)

    #def compare(self, data)

    def plot(self, codename_or_data, ref_code="WIEN2k"):
        import numpy as np
        import matplotlib.pyplot as plt

        ref_data = self._data[ref_code]

        data = codename_or_data
        if isinstance(codename_or_data, str):
            data = self._data[codename_or_data]

        records = ref_data.values()
        print (records)

        attr_names = ["b0",]

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        for aname in attr_names:
            # Sort records according to the value of the attribute aname.
            records.sort(key = lambda t: getattr(t, aname))
            #print(records)

            ord_symbols = [r.symbol for r in records]
            xticks = []

            for (i, osym) in enumerate(ord_symbols):
                ref_value = getattr(ref_data[osym], aname)
                value = getattr(data[osym], aname)

                err = 100 * (value - ref_value) / ref_value

                ax.plot(i, err, "ro")
                xticks.append(i)

            ax.plot(xticks, np.zeros(len(xticks)), "b-")

            ax.set_xticks(xticks)
            ax.set_xticklabels(ord_symbols, rotation="vertical")

        # Set xticks and labels.
        #ax.grid(True)
        #ax.set_xlabel("Ecut [Ha]")
        #ax.set_ylabel("$\Delta$ Etotal [meV]")

        plt.show()

##########################################################################################
