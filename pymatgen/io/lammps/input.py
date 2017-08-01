# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module implements classes for reading and generating Lammps input.
"""

import os
from string import Template

from monty.json import MSONable


__author__ = "Kiran Mathew, Brandon Wood"
__email__ = "kmathew@lbl.gov, b.wood@berkeley.edu"
__credits__ = "Navnidhi Rajput"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class LammpsInput(MSONable):

    def __init__(self, contents, settings, delimiter):
        self.contents = contents
        self.settings = settings or {}
        self.delimiter = delimiter
        # make read_data configurable i.e "read_data $${data_file}"
        self._map_param_to_identifier("read_data", "data_file")
        # log $${log_file}
        self._map_param_to_identifier("log", "log_file")

    def _map_param_to_identifier(self, param, identifier):
        delimited_identifier = self.delimiter+"{"+identifier+"}"
        if delimited_identifier not in self.contents:
            if self.contents.find(param) >= 0:
                self.settings[identifier] = \
                    self.contents.split(param)[-1].split("\n")[0].expandtabs().strip()
                self.contents = \
                    self.contents.replace(self.settings[identifier],
                                          delimited_identifier, 1)
            # if log is missing add it to the input
            elif param == "log":
                self.contents = "log {} \n".format(delimited_identifier)+self.contents

    def __str__(self):
        template = self.get_template(self.__class__.__name__,
                                     delimiter=self.delimiter)
        template_string = template(self.contents)

        unclean_template = template_string.safe_substitute(self.settings)

        clean_template = filter(lambda l: self.delimiter not in l,
                                unclean_template.split('\n'))

        return '\n'.join(clean_template)

    @classmethod
    def from_string(cls, input_string, settings=None, delimiter="$$"):
        return cls(input_string, settings, delimiter)

    @classmethod
    def from_file(cls, input_file, settings=None, delimiter="$$"):
        with open(input_file) as f:
            return cls.from_string(f.read(), settings, delimiter)

    def write_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.__str__())

    @staticmethod
    def get_template(name, delimiter):
        return type(name, (Template,), {"delimiter": delimiter})
