# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements methods for writing LAMMPS input files.
"""

from __future__ import annotations

import os
import re
import shutil
import warnings
from pathlib import Path
from string import Template

from monty.dev import deprecated
from monty.json import MSONable

from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.template import TemplateInputGen

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


class LammpsRun(MSONable):
    """
    Examples for various simple LAMMPS runs with given simulation box,
    force field and a few more settings. Experienced LAMMPS users should
    consider using write_lammps_inputs method with more sophisticated
    templates.

    """

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

    def __init__(self, script_template, settings, data, script_filename):
        """
        Base constructor.

        Args:
            script_template (str): String template for input script
                with placeholders. The format for placeholders has to
                be '$variable_name', e.g., '$temperature'
            settings (dict): Contains values to be written to the
                placeholders, e.g., {'temperature': 1}.
            data (LammpsData or str): Data file as a LammpsData
                instance or path to an existing data file. Default to
                None, i.e., no data file supplied. Useful only when
                read_data cmd is in the script.
            script_filename (str): Filename for the input script.
        """
        self.script_template = script_template
        self.settings = settings
        self.data = data
        self.script_filename = script_filename

    def write_inputs(self, output_dir, **kwargs):
        """
        Writes all input files (input script, and data if needed).
        Other supporting files are not handled at this moment.

        Args:
            output_dir (str): Directory to output the input files.
            **kwargs: kwargs supported by LammpsData.write_file.
        """
        write_lammps_inputs(
            output_dir=output_dir,
            script_template=self.script_template,
            settings=self.settings,
            data=self.data,
            script_filename=self.script_filename,
            **kwargs,
        )

    @classmethod
    def md(cls, data, force_field, temperature, nsteps, other_settings=None):
        r"""
        Example for a simple MD run based on template md.txt.

        Args:
            data (LammpsData or str): Data file as a LammpsData
                instance or path to an existing data file.
            force_field (str): Combined force field related cmds. For
                example, 'pair_style eam\npair_coeff * * Cu_u3.eam'.
            temperature (float): Simulation temperature.
            nsteps (int): No. of steps to run.
            other_settings (dict): other settings to be filled into
                placeholders.
        """
        template_path = os.path.join(cls.template_dir, "md.txt")
        with open(template_path) as f:
            script_template = f.read()
        settings = other_settings.copy() if other_settings is not None else {}
        settings.update({"force_field": force_field, "temperature": temperature, "nsteps": nsteps})
        script_filename = "in.md"
        return cls(
            script_template=script_template,
            settings=settings,
            data=data,
            script_filename=script_filename,
        )


class LammpsTemplateGen(TemplateInputGen):
    """
    Creates an InputSet object for a LAMMPS run based on a template file.
    The input script is constructed by substituting variables into placeholders
    in the template file using python's Template.safe_substitute() function.
    The data file containing coordinates and topology information can be provided
    as a LammpsData instance. Alternatively, you can include a read_data command
    in the template file that points to an existing data file.
    Other supporting files are not handled at the moment.

    To write the input files to a directory, call LammpsTemplateSet.write_input()
    See pymatgen.io.template.py for additional documentation of this method.
    """

    def get_input_set(  # type: ignore
        self,
        script_template: str | Path,
        settings: dict | None = None,
        script_filename: str = "in.lammps",
        data: LammpsData | CombinedData | None = None,
        data_filename: str = "system.data",
    ):
        """
        Args:
            script_template: String template for input script with
                placeholders. The format for placeholders has to be
                '$variable_name', e.g., '$temperature'
            settings: Contains values to be written to the
                placeholders, e.g., {'temperature': 1}. Default to None.
            data: Data file as a LammpsData instance. Default to None, i.e., no
                data file supplied. Note that a matching 'read_data' command
                must be provided in the script template in order for the data
                file to actually be read.
            script_filename: Filename for the input file.
            data_filename: Filename for the data file, if provided.
        """
        input_set = super().get_input_set(template=script_template, variables=settings, filename=script_filename)

        if data:
            input_set.update({data_filename: data})
        return input_set


@deprecated(LammpsTemplateGen, "This method will be retired in the future. Consider using LammpsTemplateSet instead.")
def write_lammps_inputs(
    output_dir,
    script_template,
    settings=None,
    data=None,
    script_filename="in.lammps",
    make_dir_if_not_present=True,
    **kwargs,
):
    """
    Writes input files for a LAMMPS run. Input script is constructed
    from a str template with placeholders to be filled by custom
    settings. Data file is either written from a LammpsData
    instance or copied from an existing file if read_data cmd is
    inspected in the input script. Other supporting files are not
    handled at the moment.

    Args:
        output_dir (str): Directory to output the input files.
        script_template (str): String template for input script with
            placeholders. The format for placeholders has to be
            '$variable_name', e.g., '$temperature'
        settings (dict): Contains values to be written to the
            placeholders, e.g., {'temperature': 1}. Default to None.
        data (LammpsData or str): Data file as a LammpsData instance or
            path to an existing data file. Default to None, i.e., no
            data file supplied. Useful only when read_data cmd is in
            the script.
        script_filename (str): Filename for the input script.
        make_dir_if_not_present (bool): Set to True if you want the
            directory (and the whole path) to be created if it is not
            present.
        **kwargs: kwargs supported by LammpsData.write_file.

    Examples:
        >>> eam_template = '''units           metal
        ... atom_style      atomic
        ...
        ... lattice         fcc 3.615
        ... region          box block 0 20 0 20 0 20
        ... create_box      1 box
        ... create_atoms    1 box
        ...
        ... pair_style      eam
        ... pair_coeff      1 1 Cu_u3.eam
        ...
        ... velocity        all create $temperature 376847 loop geom
        ...
        ... neighbor        1.0 bin
        ... neigh_modify    delay 5 every 1
        ...
        ... fix             1 all nvt temp $temperature $temperature 0.1
        ...
        ... timestep        0.005
        ...
        ... run             $nsteps'''
        >>> write_lammps_inputs('.', eam_template, settings={'temperature': 1600.0, 'nsteps': 100})
        >>> with open('in.lammps') as f:
        ...     script = f.read()
        ...
        >>> print(script)
        units           metal
        atom_style      atomic

        lattice         fcc 3.615
        region          box block 0 20 0 20 0 20
        create_box      1 box
        create_atoms    1 box

        pair_style      eam
        pair_coeff      1 1 Cu_u3.eam

        velocity        all create 1600.0 376847 loop geom

        neighbor        1.0 bin
        neigh_modify    delay 5 every 1

        fix             1 all nvt temp 1600.0 1600.0 0.1

        timestep        0.005

        run             100


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
            warnings.warn(f"No data file supplied. Skip writing {data_filename}.")
