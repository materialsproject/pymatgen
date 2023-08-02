---
layout: default
title: pymatgen.io.lammps.inputs.md
nav_exclude: true
---

# pymatgen.io.lammps.inputs module

This module implements methods for reading/manupilating/writing LAMMPS input files.
It does not implement methods for automatically creating inputs based on a structure
and computation type. For this, see the InputSet and InputGenerator in sets.py, or
[https://github.com/Matgenix/atomate2-lammps](https://github.com/Matgenix/atomate2-lammps).


### _class_ pymatgen.io.lammps.inputs.LammpsInputFile(stages: list | None = None)
Bases: [`InputFile`](pymatgen.io.core.md#pymatgen.io.core.InputFile)

Class representing a LAMMPS input settings file, e.g. in.lammps.
Allows for LAMMPS input generation in line/stage wise manner. A stage
here is defined as a block of LAMMPS input commands usually performing a
specific task during the simulation such as energy minimization or
NPT/NVT runs. But more broadly, a stage can also be a block of LAMMPS
input where the simulation box is set up, a set of variables are declared or
quantities are computed.

The LammpsInputFile is defined by the attribute stages,
i.e. a list of dicts each with keys stage_name and commands,
defining the stage names and the corresponding LAMMPS input settings (list of tuples of two strings each).
The structure is the following:


```
``
```

\`
stages = [

> {“stage_name”: “Stage 1”, “commands”: [(cmd1, args1), (cmd2, args2)]},
> {“stage_name”: “Stage 2”, “commands”: [(cmd3, args3)]}

## ]

where cmd’s are the LAMMPS command names (e.g., “units”, or “pair_coeff”)
and the args are the corresponding arguments.
“Stage 1” and “Stage 2” are examples of stage names.


* **param stages**

    list of LAMMPS input settings.



#### add_commands(stage_name: str, commands: str | list[str] | dict)
Method to add a LAMMPS commands and their arguments to a stage of
the LammpsInputFile. The stage name should be provided: a default behavior
is avoided here to avoid mistakes (e.g., the commands are added to the wrong stage).

### Example

In order to add the command `pair_coeff 1 1 morse 0.0580 3.987 3.404`
to the stage “Definition of the potential”, simply use


```
``
```

\`
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands=”pair_coeff 1 1 morse 0.0580 3.987 3.404”

### )

To add multiple commands, use a dict or a list, e.g.,


```
``
```

\`
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands=[“pair_coeff 1 1 morse 0.0580 3.987 3.404”, “units atomic”]

)
your_input_file.add_commands(

> stage_name=”Definition of the potential”,
> commands={“pair_coeff”: “1 1 morse 0.0580 3.987 3.404”, “units”: “atomic”}

### )


* **param stage_name**

    name of the stage to which the command should be added.



* **type stage_name**

    str



* **param commands**

    LAMMPS command, with or without the arguments.



* **type commands**

    str or list or dict



#### add_stage(stage: dict | None = None, commands: str | list[str] | dict[str, str | float] | None = None, stage_name: str | None = None, after_stage: str | int | None = None)
Adds a new stage to the LammpsInputFile, either from a whole stage (dict) or
from a stage_name and commands. Both ways are mutually exclusive.

### Examples

1) In order to add a stage defining the force field to be used, you can use:


```
``
```

\`
your_input_file.add_stage(

> commands=[“pair_coeff 1 1 morse 0.0580 3.987 3.404”, “pair_coeff 1 4 morse 0.0408 1.399 3.204”],
> stage_name=”Definition of the force field”

### )

### or

your_input_file.add_stage(

    {

        “stage_name”: “Definition of the force field”,
        “commands”: [

        > (“pair_coeff”, “1 1 morse 0.0580 3.987 3.404”),
        > (“pair_coeff”, “1 4 morse 0.0408 1.399 3.204”)

        ],

    }

### )

2) Another stage could consist in an energy minimization. In that case, the commands could look like


```
``
```

\`
commands = [

> “thermo 1”,
> “thermo_style custom step lx ly lz press pxx pyy pzz pe”,
> “dump dmp all atom 5 run.dump”,
> “min_style cg”,
> “fix 1 all box/relax iso 0.0 vmax 0.001”,
> “minimize 1.0e-16 1.0e-16 5000 10000”,
> “write_data run.data”

### ]

or a dictionary such as {“thermo”: 1, …}, or a string with a single command (e.g., “units atomic”).


* **param stage**

    if provided, this is the stage that will be added to the LammpsInputFile.stages



* **type stage**

    dict



* **param commands**

    if provided, these are the LAMMPS command(s)
    that will be included in the stage to add.
    Can pass a list of LAMMPS commands with their arguments.
    Also accepts a dictionary of LAMMPS commands and
    corresponding arguments as key, value pairs.
    A single string can also be passed (single command together with its arguments).
    Not used in case a whole stage is given.



* **type commands**

    str or list or dict



* **param stage_name**

    If a stage name is mentioned, the commands are added
    under that stage block, else the new stage is named from numbering.
    If given, stage_name cannot be one of the already present stage names.
    Not used in case a whole stage is given.



* **type stage_name**

    str



* **param after_stage**

    Name of the stage after which this stage should be added.
    If None, the stage is added at the end of the LammpsInputFile.



* **type after_stage**

    str



#### append(lmp_input_file: LammpsInputFile)
Appends a LammpsInputFile to another. The stages are merged,
and the numbering of stages/comments is either kept the same or updated.


* **Parameters**

    **lmp_input_file** (*LammpsInputFile*) – LammpsInputFile to append.



#### contains_command(command: str, stage_name: str | None = None)
Returns whether a given command is present in the LammpsInputFile.
A stage name can be given; in this case the search will happen only for this stage.


* **Parameters**


    * **command** (*str*) – String with the command to find in the input file (e.g., “units”).


    * **stage_name** (*str*) – String giving the stage name where the change should take place.



* **Returns**

    True if the command is present, False if not.



#### _classmethod_ from_file(path: str | Path, ignore_comments: bool = False, keep_stages: bool = False)
Creates an InputFile object from a file.


* **Parameters**


    * **path** (*str** or **path*) – Filename to read, including path.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the input file.


    * **keep_stages** (*bool*) – True if the block structure from the input file should be kept.
    If False, a single block is assumed.



* **Returns**

    LammpsInputFile



#### _classmethod_ from_str(contents: str, ignore_comments: bool = False, keep_stages: bool = False)
Helper method to parse string representation of LammpsInputFile.
If you created the input file by hand, there is no guarantee that the representation
will be perfect as it is difficult to account for all the cosmetic changes you
could have done on your input script. Always check that you have what you want !
By default, a single stage containing all the input settings is created.
If the block structure of your input file should be kept and stored as
different stages, set keep_stages to True.


* **Parameters**


    * **contents** (*str*) – String representation of LammpsInputFile.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the input file.


    * **keep_stages** (*bool*) – True if the block structure from the input file should be kept.
    If False, a single block is assumed.



* **Returns**

    LammpsInputFile



#### _classmethod_ from_string(\*args, \*\*kwds)
from_string is deprecated!
Use from_str instead


#### get_args(command: str, stage_name: str | None = None)
Given a command, returns the corresponding arguments (or list of arguments) in the LammpsInputFile.
A stage name can be given; in this case the search will happen only for this stage.
If a stage name is not given, the command will be searched for through all of them.
If the command is not found, an empty list is returned.


* **Parameters**


    * **command** (*str*) – String with the command to find in the input file (e.g., “units”).


    * **stage_name** (*str*) – String giving the stage name where the change should take place.



* **Returns**

    Value of the argument corresponding to the command.
    List if the same command is used multiple times.



#### get_string(ignore_comments: bool = False, keep_stages: bool = True)
Generates and ² the string representation of the LammpsInputFile.
Stages are separated by empty lines.
The headers of the stages will be put in comments preceding each stage.
Other comments will be put inline within stages, where they have been added.


* **Parameters**


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the InputFile.


    * **keep_stages** (*bool*) – If True, the string is formatted in a block structure with stage names
    and newlines that differentiate commands in the respective stages of the InputFile.
    If False, stage names are not printed and all commands appear in a single block.


Returns: String representation of the LammpsInputFile.


#### merge_stages(stage_names: list[str])
Merges multiple stages of a LammpsInputFile together.
The merged stage will be at the same index as the first of the stages to be merged.
The others will appear in the same order as provided in the list. Other non-merged stages will follow.


* **Parameters**

    **stage_names** (*list*) – list of strings giving the names of the stages to be merged.



#### _property_ ncomments(_: in_ )
Returns the number of comments in the current LammpsInputFile. Includes the blocks of comments as well
as inline comments (comment lines within blocks of LAMMPS commands).


#### _property_ nstages(_: in_ )
Returns the number of stages in the current LammpsInputFile.


#### remove_command(command: str, stage_name: str | list[str] | None = None, remove_empty_stages: bool = True)
Removes a given command from a given stage. If no stage is given, removes all occurrences of the command.
In case removing a command completely empties a stage, the choice whether to keep this stage in the
LammpsInputFile is given by remove_empty_stages.


* **Parameters**


    * **command** (*str*) – command to be removed.


    * **stage_name** (*str** or **list*) – names of the stages where the command should be removed.


    * **remove_empty_stages** (*bool*) – whether to remove the stages emptied by removing the command or not.



#### remove_stage(stage_name: str)
Removes a whole stage from the LammpsInputFile.


* **Parameters**

    **stage_name** (*str*) – name of the stage to remove.



#### rename_stage(stage_name: str, new_name: str)
Renames a stage stage_name from LammpsInputFile into new_name.
First checks that the stage to rename is present, and that
the new name is not already a stage name.


* **Parameters**


    * **stage_name** (*str*) – name of the stage to rename.


    * **new_name** (*str*) – new name of the stage.



#### set_args(command: str, argument: str, stage_name: str | None = None, how: str | int | list[int] = 'all')
Sets the arguments for the given command to the given string.
If the command is not found, nothing is done. Use LammpsInputFile.add_commands instead.
If a stage name is specified, it will be replaced or set only for this stage.
If no stage name is given, it will apply the change in all of them that contain the given command.
If the command is set multiple times in the file/stage, it will be replaced based on “how”:
either the first occurrence, all of them, or the index of the occurrence.


* **Parameters**


    * **command** (*str*) – String representing the command to change, e.g., “units”.


    * **argument** (*str*) – String with the new value for the command, e.g., “atomic”.


    * **stage_name** (*str*) – String giving the stage name where the change should take place.


    * **how** (*str** or **int** or **list*) – “all” for changing all occurrences of the command within the stage_name
    or the whole InputFile, “first” for the first occurrence, int i for the i-th time the command
    is present in the stage_name or the whole InputFile, starting at 0. Can be a list of indexes as well.



#### _property_ stages_names(_: lis_ )
List of names for all the stages present in stages.


#### write_file(filename: str | Path, ignore_comments: bool = False, keep_stages: bool = True)
Writes the input file.


* **Parameters**


    * **filename** (*str** or **path*) – The filename to output to, including path.


    * **ignore_comments** (*bool*) – True if only the commands should be kept from the InputFile.


    * **keep_stages** (*bool*) – True if the block structure from the InputFile should be kept.
    If False, a single block is assumed.



### _class_ pymatgen.io.lammps.inputs.LammpsRun(script_template, settings, data, script_filename)
Bases: `MSONable`

Examples for various simple LAMMPS runs with given simulation box,
force field and a few more settings. Experienced LAMMPS users should
consider using write_lammps_inputs method with more sophisticated
templates.

Base constructor.


* **Parameters**


    * **script_template** (*str*) – String template for input script
    with placeholders. The format for placeholders has to
    be ‘$variable_name’, e.g., ‘$temperature’


    * **settings** (*dict*) – Contains values to be written to the
    placeholders, e.g., {‘temperature’: 1}.


    * **data** ([*LammpsData*](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData)* or **str*) – Data file as a LammpsData
    instance or path to an existing data file. Default to
    None, i.e., no data file supplied. Useful only when
    read_data cmd is in the script.


    * **script_filename** (*str*) – Filename for the input script.



#### _classmethod_ md(data, force_field, temperature, nsteps, other_settings=None)
Example for a simple MD run based on template md.template.


* **Parameters**


    * **data** ([*LammpsData*](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData)* or **str*) – Data file as a LammpsData
    instance or path to an existing data file.


    * **force_field** (*str*) – Combined force field related cmds. For
    example, ‘pair_style eamnpair_coeff \* \* Cu_u3.eam’.


    * **temperature** (*float*) – Simulation temperature.


    * **nsteps** (*int*) – No. of steps to run.


    * **other_settings** (*dict*) – other settings to be filled into
    placeholders.



#### template_dir(_ = '/Users/shyue/repos/pymatgen/pymatgen/io/lammps/templates_ )

#### write_inputs(output_dir, \*\*kwargs)
Writes all input files (input script, and data if needed).
Other supporting files are not handled at this moment.


* **Parameters**


    * **output_dir** (*str*) – Directory to output the input files.


    * **\*\*kwargs** – kwargs supported by LammpsData.write_file.



### _class_ pymatgen.io.lammps.inputs.LammpsTemplateGen()
Bases: [`TemplateInputGen`](pymatgen.io.template.md#pymatgen.io.template.TemplateInputGen)

Creates an InputSet object for a LAMMPS run based on a template file.
The input script is constructed by substituting variables into placeholders
in the template file using python’s Template.safe_substitute() function.
The data file containing coordinates and topology information can be provided
as a LammpsData instance. Alternatively, you can include a read_data command
in the template file that points to an existing data file.
Other supporting files are not handled at the moment.

To write the input files to a directory, call LammpsTemplateSet.write_input()
See pymatgen.io.template.py for additional documentation of this method.


#### get_input_set(script_template: str | Path, settings: dict | None = None, script_filename: str = 'in.lammps', data: [LammpsData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData) | [CombinedData](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData) | None = None, data_filename: str = 'system.data')

* **Parameters**


    * **script_template** – String template for input script with
    placeholders. The format for placeholders has to be
    ‘$variable_name’, e.g., ‘$temperature’


    * **settings** – Contains values to be written to the
    placeholders, e.g., {‘temperature’: 1}. Default to None.


    * **data** – Data file as a LammpsData instance. Default to None, i.e., no
    data file supplied. Note that a matching ‘read_data’ command
    must be provided in the script template in order for the data
    file to actually be read.


    * **script_filename** – Filename for the input file.


    * **data_filename** – Filename for the data file, if provided.