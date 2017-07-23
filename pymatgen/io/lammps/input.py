# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module implements classes for reading and generating Lammps input.

For the ease of management we divide LAMMPS input into 2 files:

    1.Data file: All structure related settings such as the atomic positions,
            bonds, angles, dihedrals, corresponding parametrizations etc are
            set in the data file.

    2. Control file: This is the main input file that should be fed to the
            lammps binary. The main input file consists of the path to the
            afore-mentioned data file and the job control parameters such as
            the ensemble type(NVT, NPT etc), max number of iterations etc.
"""

import json
import os
import string
from collections import OrderedDict, defaultdict

__author__ = "Kiran Mathew, Brandon Wood"
__email__ = "kmathew@lbl.gov, b.wood@berkeley.edu"
__credits__ = "Navnidhi Rajput"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class LammpsInputTemplate(string.Template):
    delimiter = '$$'


class LammpsInput(defaultdict):

    defaults = {}

    def __init__(self, template_file=None, **kwargs):
        self.template_file = template_file
        self.update(kwargs)

    def __str__(self):
        if self.template_file:
            return self.template_file_str()
        else:
            lines = ""
            for k1, v1 in self.items():
                if isinstance(v1, dict):
                    v1 = v1.values()
                if isinstance(v1, list):
                    for x in v1:
                        lines = "".join([lines, "{} ".format(k1)])
                        lines = "".join([lines, str(x), os.linesep])
                else:
                    lines = "".join([lines, "{} ".format(k1)])
                    lines = "".join([lines, " {}{}".format(str(v1), os.linesep)])
            return lines

    def template_file_str(self):
        with open(self.template_file) as f:
            a = LammpsInputTemplate(f.read())

            # set substitution dict for replacements into the template
            subs_dict = {k: v for k, v in self.items()
                         if v is not None}  # clean null values

            for k, v in self.defaults.items():
                subs_dict.setdefault(k, v)

            # might contain unused parameters as leftover $$
            unclean_template = a.safe_substitute(subs_dict)

            clean_template = filter(lambda l: "$$" not in l,
                                    unclean_template.split('\n'))

            return '\n'.join(clean_template)

    @classmethod
    def from_file(cls, filename, user_lammps_settings=None):
        user_lammps_settings = user_lammps_settings or {}
        try:
            with open(filename) as f:
                config_dict = json.load(f, object_pairs_hook=OrderedDict)
        except ValueError:
            with open(filename, 'r') as f:
                data = f.read().splitlines()
            config_dict = OrderedDict()
            for line in data:
                if line and not line.startswith("#"):
                    spt_line = (line.split(None, 1))
                    if spt_line[0] in config_dict:
                        if isinstance(config_dict[spt_line[0]], list):
                            config_dict[spt_line[0]].append(spt_line[1])
                        else:
                            config_dict[spt_line[0]] = [config_dict[spt_line[0]], spt_line[1]]
                    else:
                        config_dict[spt_line[0]] = spt_line[1]

        config_dict.update(user_lammps_settings)
        return cls(**config_dict)

    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.__str__())

    def as_dict(self):
        d = {"template_file": self.template_file}
        d.update(self)
        return d

    @classmethod
    def from_dict(cls, d):
        template_file = d.pop("template_file")
        return cls(template_file, **d)


