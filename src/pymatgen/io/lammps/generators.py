"""
Input set generators for LAMMPS.
This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
For instance, pymatgen will not detect whether a given variable should be adapted based on others
(e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
For additional flexibility and automation, use the atomate2-lammps implementation
(https://github.com/Matgenix/atomate2-lammps).
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from string import Template

from monty.io import zopen
from monty.json import MSONable

from pymatgen.core import Structure
from pymatgen.io.core import InputGenerator
from pymatgen.io.lammps.data import CombinedData, LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile
from pymatgen.io.lammps.sets import LammpsInputSet

__author__ = "Ryan Kingsbury, Guillaume Brunin (Matgenix)"
__copyright__ = "Copyright 2021, The Materials Project"
__version__ = "0.2"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_DIR = f"{MODULE_DIR}/templates"


@dataclass
class BaseLammpsGenerator(InputGenerator):
    r"""
    Base class to generate LAMMPS input sets.
    Uses template files for the input. The variables that can be changed
    in the input template file are those starting with a $ sign, e.g. $nsteps.
    This generic class is specialized for each template in subclasses, e.g. LammpsMinimization.
    You can create a template for your own task following those present in pymatgen/io/lammps/templates.
    The parameters are then replaced based on the values found
    in the settings dictionary that you provide, e.g. `{"nsteps": 1000}`.

    Attributes:
        template: Path (string) to the template file used to create the InputFile for LAMMPS.
        calc_type: Human-readable string used to briefly describe the type of computations performed by LAMMPS.
        settings: Dictionary containing the values of the parameters to replace in the template.
        keep_stages: If True, the string is formatted in a block structure with stage names
        and newlines that differentiate commands in the respective stages of the InputFile.
        If False, stage names are not printed and all commands appear in a single block.

    /!\ This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    inputfile: LammpsInputFile | None = field(default=None)
    template: str | Path | LammpsInputFile = field(default_factory=str)
    data: LammpsData | CombinedData | None = field(default=None)
    settings: dict = field(default_factory=dict)
    calc_type: str = field(default="lammps")
    keep_stages: bool = field(default=True)

    def get_input_set(self, structure: Structure | LammpsData | CombinedData) -> LammpsInputSet:
        """Generate a LammpsInputSet from the structure/data, tailored to the template file."""
                
        data = LammpsData.from_structure(structure, atom_style=self.settings.get('atom_style', "full")) if isinstance(structure, Structure) else structure

        # Load the template
        if Path(self.template).is_file():
            with zopen(self.template, mode="r") as file:
                template_str = file.read()
        elif isinstance(self.template, LammpsInputFile):
            template_str = self.template.get_str()
        else:
            template_str = self.template

        # Replace all variables
        input_str = Template(template_str).safe_substitute(**self.settings)
        # Get LammpsInputFile
        input_file = LammpsInputFile.from_str(input_str, keep_stages=self.keep_stages)

        # Get the LammpsInputSet from the InputFile and data
        return LammpsInputSet(
            inputfile=input_file,
            data=data,
            calc_type=self.calc_type,
            template_file=self.template,
        )   


class LammpsMinimization(BaseLammpsGenerator):
    """
    Generator that yields a LammpsInputSet tailored for minimizing the energy of a system by iteratively
    adjusting atom coordinates.
    Example usage:
    ```
    structure = Structure.from_file("mp-149.cif")
    lmp_minimization = LammpsMinimization(units="atomic").get_input_set(structure)
    ```.

    Do not forget to specify the force field, otherwise LAMMPS will not be able to run!

    This InputSet and InputGenerator implementation is based on templates and is not intended to be very flexible.
    For instance, pymatgen will not detect whether a given variable should be adapted based on others
    (e.g., the number of steps from the temperature), it will not check for convergence nor will it actually run LAMMPS.
    For additional flexibility and automation, use the atomate2-lammps implementation
    (https://github.com/Matgenix/atomate2-lammps).
    """

    def __init__(
        self,
        template: str | None = None,
        units: str = "metal",
        atom_style: str = "full",
        dimension: int = 3,
        boundary: str = "p p p",
        read_data: str = "system.data",
        force_field: str = "Unspecified force field!",
        keep_stages: bool = False,
    ) -> None:
        r"""
        Args:
            structure : Structure | LammpsData
                Structure or LammpsData object for the simulation.
            **kwargs : dict
                Additional keyword arguments to pass to the InputSet from pmg.      
        """
        input_settings = self.settings.input_settings
        atom_style = input_settings.get('Initialization').get('atom_style', "full")
        species = ' '.join(set([s.symbol for s in data.species])) if isinstance(data, Structure) else ''

        if isinstance(data, Path):
            data = LammpsData.from_file(data, atom_style=atom_style)
        if isinstance(data, str):
            data = LammpsData.from_str(data, atom_style=atom_style)
        if isinstance(data, Structure):
            data = LammpsData.from_structure(data, atom_style=atom_style)
            warnings.warn("Structure provided, converting to LammpsData object.")
        #if isinstance(data, LammpsInterchange):
        #    warnings.warn("Interchange is experimental and may not work as expected. Use with caution. Ensure FF units are consistent with LAMMPS.")
            #write unit convertor here
        #    data.to_lammps_datafile("interchange_data.lmp")
        #    data = LammpsData.from_file("interchange_data.lmp", atom_style=atom_style)
            #validate data here: ff coeffs style, atom_style, etc. have to be updated into the input_set_generator.settings'''
            
        if not data.force_field and not force_field:
            raise ValueError("Force field not specified!")
        
        if species:
            input_settings['Outputs'].update({"dump_modify": f"d1 sort id element {species}"})
        
        write_data = {}
        if additional_data:
            write_data.update({"extra.data": additional_data})
        else:
            self.inputfile.remove_stage(stage_name="AdditionalData")
        
        if force_field:
            write_data.update({"forcefield.lammps": force_field})
        else:
            self.inputfile.remove_stage(stage_name="ForceField") #implies FF is in the datafile
            
            
        if self.override_updates:
            return LammpsInputSet(
                inputfile=self.inputfile,
                data=data,
                calc_type=self.calc_type,
                additional_data=additional_data,
                **kwargs
            )
        
        if self.inputfile.contains_command(command='fix'):
            self.inputfile.remove_command(command='fix', stage_name="Ensemble", remove_empty_stages=False)
        
        
        for stage, stage_data in input_settings.items():
            for key, value in stage_data.items():
                if self.inputfile.contains_command(command=key, stage_name=stage) and key not in ['fix']:
                    self.inputfile.set_args(stage_name=stage, command=key, argument=str(value))                        
                else:
                    value = [value] if not isinstance(value, list) else value
                    for val in value:
                        self.inputfile.add_commands(stage_name=stage, commands={key: str(val)})
        
        return LammpsInputSet(
            inputfile=self.inputfile,
            data=data,
            calc_type=self.calc_type,
            additional_data=write_data,
            **kwargs
        )
