---
layout: default
title: pymatgen.io.feff.md
nav_exclude: true
---

# pymatgen.io.feff package

This package provides the modules to perform FEFF IO.

FEFF: [http://feffproject.org/feffproject-feff.html](http://feffproject.org/feffproject-feff.html)



* [pymatgen.io.feff.inputs module](pymatgen.io.feff.inputs.md)


    * [`Atoms`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms)


        * [`Atoms.atoms_string_from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms.atoms_string_from_file)


        * [`Atoms.cluster`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms.cluster)


        * [`Atoms.cluster_from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms.cluster_from_file)


        * [`Atoms.get_lines()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms.get_lines)


        * [`Atoms.write_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Atoms.write_file)


    * [`FeffParseError`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.FeffParseError)


    * [`Header`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header)


        * [`Header.formula`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.formula)


        * [`Header.from_cif_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.from_cif_file)


        * [`Header.from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.from_file)


        * [`Header.from_str()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.from_str)


        * [`Header.from_string()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.from_string)


        * [`Header.header_string_from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.header_string_from_file)


        * [`Header.structure_symmetry`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.structure_symmetry)


        * [`Header.write_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Header.write_file)


    * [`Paths`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Paths)


        * [`Paths.write_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Paths.write_file)


    * [`Potential`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Potential)


        * [`Potential.pot_dict_from_string()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Potential.pot_dict_from_string)


        * [`Potential.pot_string_from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Potential.pot_string_from_file)


        * [`Potential.write_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Potential.write_file)


    * [`Tags`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags)


        * [`Tags.as_dict()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.as_dict)


        * [`Tags.diff()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.diff)


        * [`Tags.from_dict()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.from_dict)


        * [`Tags.from_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.from_file)


        * [`Tags.get_string()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.get_string)


        * [`Tags.proc_val()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.proc_val)


        * [`Tags.write_file()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.Tags.write_file)


    * [`get_absorbing_atom_symbol_index()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.get_absorbing_atom_symbol_index)


    * [`get_atom_map()`](pymatgen.io.feff.inputs.md#pymatgen.io.feff.inputs.get_atom_map)


* [pymatgen.io.feff.outputs module](pymatgen.io.feff.outputs.md)


    * [`Eels`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels)


        * [`Eels.as_dict()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.as_dict)


        * [`Eels.atomic_background`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.atomic_background)


        * [`Eels.energies`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.energies)


        * [`Eels.fine_structure`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.fine_structure)


        * [`Eels.from_file()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.from_file)


        * [`Eels.total_spectrum`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Eels.total_spectrum)


    * [`LDos`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.LDos)


        * [`LDos.charge_transfer_from_file()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.LDos.charge_transfer_from_file)


        * [`LDos.charge_transfer_to_string()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.LDos.charge_transfer_to_string)


        * [`LDos.from_file()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.LDos.from_file)


    * [`Xmu`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu)


        * [`Xmu.as_dict()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.as_dict)


        * [`Xmu.calc`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.calc)


        * [`Xmu.chi`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.chi)


        * [`Xmu.e_fermi`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.e_fermi)


        * [`Xmu.edge`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.edge)


        * [`Xmu.energies`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.energies)


        * [`Xmu.from_file()`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.from_file)


        * [`Xmu.material_formula`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.material_formula)


        * [`Xmu.mu`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.mu)


        * [`Xmu.mu0`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.mu0)


        * [`Xmu.relative_energies`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.relative_energies)


        * [`Xmu.source`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.source)


        * [`Xmu.wavenumber`](pymatgen.io.feff.outputs.md#pymatgen.io.feff.outputs.Xmu.wavenumber)


* [pymatgen.io.feff.sets module](pymatgen.io.feff.sets.md)


    * [`AbstractFeffInputSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet)


        * [`AbstractFeffInputSet.all_input()`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.all_input)


        * [`AbstractFeffInputSet.atoms`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.atoms)


        * [`AbstractFeffInputSet.header()`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.header)


        * [`AbstractFeffInputSet.potential`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.potential)


        * [`AbstractFeffInputSet.tags`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.tags)


        * [`AbstractFeffInputSet.write_input()`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.AbstractFeffInputSet.write_input)


    * [`FEFFDictSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet)


        * [`FEFFDictSet.atoms`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet.atoms)


        * [`FEFFDictSet.from_directory()`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet.from_directory)


        * [`FEFFDictSet.header()`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet.header)


        * [`FEFFDictSet.potential`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet.potential)


        * [`FEFFDictSet.tags`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.FEFFDictSet.tags)


    * [`MPEELSDictSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPEELSDictSet)


    * [`MPELNESSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPELNESSet)


        * [`MPELNESSet.CONFIG`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPELNESSet.CONFIG)


    * [`MPEXAFSSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPEXAFSSet)


        * [`MPEXAFSSet.CONFIG`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPEXAFSSet.CONFIG)


    * [`MPEXELFSSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPEXELFSSet)


        * [`MPEXELFSSet.CONFIG`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPEXELFSSet.CONFIG)


    * [`MPXANESSet`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPXANESSet)


        * [`MPXANESSet.CONFIG`](pymatgen.io.feff.sets.md#pymatgen.io.feff.sets.MPXANESSet.CONFIG)