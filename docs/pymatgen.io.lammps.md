---
layout: default
title: pymatgen.io.lammps.md
nav_exclude: true
---

# pymatgen.io.lammps package

IO for LAMMPS.



* [pymatgen.io.lammps.data module](pymatgen.io.lammps.data.md)


    * [`CombinedData`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData)


        * [`CombinedData.as_lammpsdata()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.as_lammpsdata)


        * [`CombinedData.disassemble()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.disassemble)


        * [`CombinedData.from_ff_and_topologies()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.from_ff_and_topologies)


        * [`CombinedData.from_files()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.from_files)


        * [`CombinedData.from_lammpsdata()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.from_lammpsdata)


        * [`CombinedData.from_structure()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.from_structure)


        * [`CombinedData.get_string()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.get_string)


        * [`CombinedData.parse_xyz()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.parse_xyz)


        * [`CombinedData.structure`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.CombinedData.structure)


    * [`ForceField`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField)


        * [`ForceField.masses`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.masses)


        * [`ForceField.force_field`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.force_field)


        * [`ForceField.maps`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.maps)


        * [`ForceField.from_dict()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.from_dict)


        * [`ForceField.from_file()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.from_file)


        * [`ForceField.to_file()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.ForceField.to_file)


    * [`LammpsBox`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox)


        * [`LammpsBox.get_box_shift()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox.get_box_shift)


        * [`LammpsBox.get_string()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox.get_string)


        * [`LammpsBox.to_lattice()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox.to_lattice)


        * [`LammpsBox.volume`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsBox.volume)


    * [`LammpsData`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData)


        * [`LammpsData.disassemble()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.disassemble)


        * [`LammpsData.from_ff_and_topologies()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.from_ff_and_topologies)


        * [`LammpsData.from_file()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.from_file)


        * [`LammpsData.from_structure()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.from_structure)


        * [`LammpsData.get_string()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.get_string)


        * [`LammpsData.set_charge_atom()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.set_charge_atom)


        * [`LammpsData.set_charge_atom_type()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.set_charge_atom_type)


        * [`LammpsData.structure`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.structure)


        * [`LammpsData.write_file()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.LammpsData.write_file)


    * [`Topology`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.Topology)


        * [`Topology.from_bonding()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.Topology.from_bonding)


    * [`lattice_2_lmpbox()`](pymatgen.io.lammps.data.md#pymatgen.io.lammps.data.lattice_2_lmpbox)


* [pymatgen.io.lammps.generators module](pymatgen.io.lammps.generators.md)


    * [`BaseLammpsGenerator`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator)


        * [`BaseLammpsGenerator.calc_type`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.calc_type)


        * [`BaseLammpsGenerator.get_input_set()`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.get_input_set)


        * [`BaseLammpsGenerator.keep_stages`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.keep_stages)


        * [`BaseLammpsGenerator.settings`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.settings)


        * [`BaseLammpsGenerator.template`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.BaseLammpsGenerator.template)


    * [`LammpsMinimization`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization)


        * [`LammpsMinimization.atom_style`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.atom_style)


        * [`LammpsMinimization.boundary`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.boundary)


        * [`LammpsMinimization.dimension`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.dimension)


        * [`LammpsMinimization.force_field`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.force_field)


        * [`LammpsMinimization.read_data`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.read_data)


        * [`LammpsMinimization.settings`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.settings)


        * [`LammpsMinimization.template`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.template)


        * [`LammpsMinimization.units`](pymatgen.io.lammps.generators.md#pymatgen.io.lammps.generators.LammpsMinimization.units)


* [pymatgen.io.lammps.inputs module](pymatgen.io.lammps.inputs.md)


    * [`LammpsInputFile`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile)


        * [`LammpsInputFile.add_commands()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.add_commands)


        * [`LammpsInputFile.add_stage()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.add_stage)


        * [`LammpsInputFile.append()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.append)


        * [`LammpsInputFile.contains_command()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.contains_command)


        * [`LammpsInputFile.from_file()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_file)


        * [`LammpsInputFile.from_str()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_str)


        * [`LammpsInputFile.from_string()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.from_string)


        * [`LammpsInputFile.get_args()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.get_args)


        * [`LammpsInputFile.get_string()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.get_string)


        * [`LammpsInputFile.merge_stages()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.merge_stages)


        * [`LammpsInputFile.ncomments`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.ncomments)


        * [`LammpsInputFile.nstages`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.nstages)


        * [`LammpsInputFile.remove_command()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.remove_command)


        * [`LammpsInputFile.remove_stage()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.remove_stage)


        * [`LammpsInputFile.rename_stage()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.rename_stage)


        * [`LammpsInputFile.set_args()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.set_args)


        * [`LammpsInputFile.stages_names`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.stages_names)


        * [`LammpsInputFile.write_file()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsInputFile.write_file)


    * [`LammpsRun`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsRun)


        * [`LammpsRun.md()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsRun.md)


        * [`LammpsRun.template_dir`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsRun.template_dir)


        * [`LammpsRun.write_inputs()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsRun.write_inputs)


    * [`LammpsTemplateGen`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsTemplateGen)


        * [`LammpsTemplateGen.get_input_set()`](pymatgen.io.lammps.inputs.md#pymatgen.io.lammps.inputs.LammpsTemplateGen.get_input_set)


* [pymatgen.io.lammps.outputs module](pymatgen.io.lammps.outputs.md)


    * [`LammpsDump`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.LammpsDump)


        * [`LammpsDump.as_dict()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.LammpsDump.as_dict)


        * [`LammpsDump.from_dict()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.LammpsDump.from_dict)


        * [`LammpsDump.from_str()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.LammpsDump.from_str)


        * [`LammpsDump.from_string()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.LammpsDump.from_string)


    * [`parse_lammps_dumps()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.parse_lammps_dumps)


    * [`parse_lammps_log()`](pymatgen.io.lammps.outputs.md#pymatgen.io.lammps.outputs.parse_lammps_log)


* [pymatgen.io.lammps.sets module](pymatgen.io.lammps.sets.md)


    * [`LammpsInputSet`](pymatgen.io.lammps.sets.md#pymatgen.io.lammps.sets.LammpsInputSet)


        * [`LammpsInputSet.from_directory()`](pymatgen.io.lammps.sets.md#pymatgen.io.lammps.sets.LammpsInputSet.from_directory)


        * [`LammpsInputSet.validate()`](pymatgen.io.lammps.sets.md#pymatgen.io.lammps.sets.LammpsInputSet.validate)


* [pymatgen.io.lammps.utils module](pymatgen.io.lammps.utils.md)


    * [`LammpsRunner`](pymatgen.io.lammps.utils.md#pymatgen.io.lammps.utils.LammpsRunner)


        * [`LammpsRunner.run()`](pymatgen.io.lammps.utils.md#pymatgen.io.lammps.utils.LammpsRunner.run)


    * [`Polymer`](pymatgen.io.lammps.utils.md#pymatgen.io.lammps.utils.Polymer)