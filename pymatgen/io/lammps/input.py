# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals,\
    absolute_import

"""
This module implements methods for writing LAMMPS input files.

"""

import os
import re
import shutil
import warnings
from string import Template

from pymatgen.io.lammps.data import LammpsData


def write_lammps_inputs(output_dir, script_template, settings=None,
                        data=None, script_filename="in.lammps",
                        make_dir_if_not_present=True, **kwargs):
    """
    Writes input files for a LAMMPS run. Input script is constructed
    from a str template with placeholders to be filled by custom
    settings. Optional data file is either written from a LammpsData
    instance or copied from an existing file.

    Args:
        output_dir (str): Directory to output the input files.
        script_template (str): String template for input script with
            placeholders. The format for placeholders has to be
            '$variable_name', e.g., '$temperature'
        settings (dict): Contains values to be written to the
            placeholders, e.g., {'temperature': 1}. Default to None.
        data (LammpsData or str): Data file as a LammpsData instance or
            path to an existing data file. Default to None, i.e., no
            data file supplied.
        script_filename (str): Filename for the input script.
        make_dir_if_not_present (bool): Set to True if you want the
            directory (and the whole path) to be created if it is not
            present.
        **kwargs: kwargs supported by LammpsData.write_file.

    Returns:

    """
    variables = {} if settings is None else settings
    template = Template(script_template)
    input_script = template.safe_substitute(**variables)
    if make_dir_if_not_present and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(os.path.join(output_dir, script_filename), "w") as f:
        f.write(input_script)
    read_data = re.search(r"read_data\s+(.*)\n", input_script)
    if read_data:
        data_filename = read_data.group(1).split()[0]
        if isinstance(data, LammpsData):
            data.write_file(os.path.join(output_dir, data_filename), **kwargs)
        elif isinstance(data, str) and os.path.exists(data):
            shutil.copyfile(data, os.path.join(output_dir, data_filename))
        else:
            warnings.warn("No data file supplied. Skip writing %s."
                          % data_filename)


