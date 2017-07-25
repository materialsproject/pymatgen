# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module implements classes for reading and generating Lammps input.
"""

import os
from string import Template
from collections import defaultdict

from monty.json import MSONable


__author__ = "Kiran Mathew, Brandon Wood"
__email__ = "kmathew@lbl.gov, b.wood@berkeley.edu"
__credits__ = "Navnidhi Rajput"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class LammpsInputFileTemplate(Template):
    delimiter = "$$"


class LammpsInput(defaultdict, MSONable):

    def __init__(self, template_string, **kwargs):
        self.template_string = template_string
        self.update(kwargs)
        # if 'read_data' is configurable, make it
        if "read_data" not in self and self.template_string.find("read_data") >=0 :
            self["read_data"] = self.template_string.split("read_data")[-1].split("\n")[0].expandtabs().strip()
            self.template_string = self.template_string.replace(self["read_data"], "$${read_data}", 1)

    def __str__(self):
        template_string = LammpsInputFileTemplate(self.template_string)

        # set substitution dict for replacements into the template
        subs_dict = {k: v for k, v in self.items()
                     if v is not None}  # clean null values

        # might contain unused parameters as leftover $$
        unclean_template = template_string.safe_substitute(subs_dict)

        clean_template = filter(lambda l: LammpsInputFileTemplate.delimiter not in l,
                                unclean_template.split('\n'))

        return '\n'.join(clean_template)

    @classmethod
    def from_file(cls, filename, user_settings):
        with open(filename) as f:
            return cls(f.read(), **user_settings)

    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.__str__())

    def as_dict(self):
        d = {"template_string": self.template_string}
        d.update(self)
        return d

    @classmethod
    def from_dict(cls, d):
        template_string = d.pop("template_string")
        return cls(template_string, **d)
