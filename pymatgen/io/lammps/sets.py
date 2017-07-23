import six

from monty.json import MSONable, MontyDecoder

from pymatgen.io.lammps.data import LammpsForceFieldData, LammpsData
from pymatgen.io.lammps.input import LammpsInput


class LammpsInputSet(MSONable):
    """
    Implementation of LammpsInputSet that is initialized from a dict
    settings. It is typically used by other LammpsInputSets for
    initialization from json or yaml source files.

    Args:
        name (str): A name for the input set.
        config_dict (dict): The config dictionary to use.
        lammps_data (LammpsData): LammpsData object
        data_filename (str): name of the the lammps data file
        user_lammps_settings (dict): User lammps settings. This allows a user
            to override lammps settings, e.g., setting a different force field
            or bond type.
    """

    def __init__(self, name, lammps_input, lammps_data=None,
                 data_filename="in.data", user_lammps_settings=None):
        self.name = name
        self.lines = []
        self.lammps_input = lammps_input
        self.lammps_data = lammps_data
        self.data_filename = data_filename
        self.lammps_input["read_data"] = data_filename
        self.user_lammps_settings = user_lammps_settings or None
        if self.user_lammps_settings:
            self.lammps_input.update(self.user_lammps_settings)

    def write_input(self, filename, data_filename=None):
        """
        Get the string representation of the main input file and write it.
        Also writes the data file if the lammps_data attribute is set.

        Args:
            filename (string): name of the input file
            data_filename (string): override the data file name with this
        """
        if data_filename:
            self.lammps_input["read_data"] = data_filename
            self.data_filename = data_filename
        self.lammps_input.write_file(filename)
        # write the data file if present
        if self.lammps_data:
            print("Data file: {}".format(self.data_filename))
            self.lammps_data.write_data_file(filename=self.data_filename)

    @classmethod
    def from_file(cls, name, input_template, lammps_data=None, data_filename="in.data",
                  user_lammps_settings=None, is_forcefield=False):
        """
        Reads lammps style and JSON style input files putting the settings in an ordered dict (config_dict).
        Note: with monty.serialization.loadfn the order of paramters in the
        json file is not preserved

        Args:
            input_template (string): name of the file with the lamps control
                paramters
            lammps_data (string/LammpsData/LammpsForceFieldData): path to the
                data file or an appropriate object
            data_filename (string): name of the the lammps data file
            user_lammps_settings (dict): User lammps settings
            is_forcefield (bool): whether the data file has forcefield and
                topology info in it. This is required only if lammps_data is
                a path to the data file instead of a data object

        Returns:
            LammpsInputSet
        """
        user_lammps_settings = user_lammps_settings or {}
        lammps_input = LammpsInput(input_template, **user_lammps_settings)
        if isinstance(lammps_data, six.string_types):
            if is_forcefield:
                lammps_data = LammpsForceFieldData.from_file(lammps_data)
            else:
                lammps_data = LammpsData.from_file(lammps_data)
        return cls(name, lammps_input, lammps_data=lammps_data, data_filename=data_filename,
                   user_lammps_settings=user_lammps_settings)

    def as_dict(self):
        d = MSONable.as_dict(self)
        return d

    @classmethod
    def from_dict(cls, d):
        decoded = {k: MontyDecoder().process_decoded(v) for k, v in d.items()
                   if k not in ["@module", "@class"]}
        return cls(**decoded)
