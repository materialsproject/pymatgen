---
layout: default
title: pymatgen.alchemy.md
nav_exclude: true
---

# pymatgen.alchemy package

This package provides the modules for performing large scale transformations on
a large number of structures.



* [pymatgen.alchemy.filters module](pymatgen.alchemy.filters.md)


    * [`AbstractStructureFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.AbstractStructureFilter)


        * [`AbstractStructureFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.AbstractStructureFilter.test)


    * [`ChargeBalanceFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ChargeBalanceFilter)


        * [`ChargeBalanceFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ChargeBalanceFilter.test)


    * [`ContainsSpecieFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ContainsSpecieFilter)


        * [`ContainsSpecieFilter.as_dict()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ContainsSpecieFilter.as_dict)


        * [`ContainsSpecieFilter.from_dict()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ContainsSpecieFilter.from_dict)


        * [`ContainsSpecieFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.ContainsSpecieFilter.test)


    * [`RemoveDuplicatesFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.RemoveDuplicatesFilter)


        * [`RemoveDuplicatesFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.RemoveDuplicatesFilter.test)


    * [`RemoveExistingFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.RemoveExistingFilter)


        * [`RemoveExistingFilter.as_dict()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.RemoveExistingFilter.as_dict)


        * [`RemoveExistingFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.RemoveExistingFilter.test)


    * [`SpecieProximityFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpecieProximityFilter)


        * [`SpecieProximityFilter.as_dict()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpecieProximityFilter.as_dict)


        * [`SpecieProximityFilter.from_dict()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpecieProximityFilter.from_dict)


        * [`SpecieProximityFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpecieProximityFilter.test)


    * [`SpeciesMaxDistFilter`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpeciesMaxDistFilter)


        * [`SpeciesMaxDistFilter.test()`](pymatgen.alchemy.filters.md#pymatgen.alchemy.filters.SpeciesMaxDistFilter.test)


* [pymatgen.alchemy.materials module](pymatgen.alchemy.materials.md)


    * [`TransformedStructure`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure)


        * [`TransformedStructure.append_filter()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.append_filter)


        * [`TransformedStructure.append_transformation()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.append_transformation)


        * [`TransformedStructure.as_dict()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.as_dict)


        * [`TransformedStructure.extend_transformations()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.extend_transformations)


        * [`TransformedStructure.from_cif_string()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.from_cif_string)


        * [`TransformedStructure.from_dict()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.from_dict)


        * [`TransformedStructure.from_poscar_string()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.from_poscar_string)


        * [`TransformedStructure.from_snl()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.from_snl)


        * [`TransformedStructure.get_vasp_input()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.get_vasp_input)


        * [`TransformedStructure.redo_next_change()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.redo_next_change)


        * [`TransformedStructure.set_parameter()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.set_parameter)


        * [`TransformedStructure.structures`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.structures)


        * [`TransformedStructure.to_snl()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.to_snl)


        * [`TransformedStructure.undo_last_change()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.undo_last_change)


        * [`TransformedStructure.was_modified`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.was_modified)


        * [`TransformedStructure.write_vasp_input()`](pymatgen.alchemy.materials.md#pymatgen.alchemy.materials.TransformedStructure.write_vasp_input)


* [pymatgen.alchemy.transmuters module](pymatgen.alchemy.transmuters.md)


    * [`CifTransmuter`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.CifTransmuter)


        * [`CifTransmuter.from_filenames()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.CifTransmuter.from_filenames)


    * [`PoscarTransmuter`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.PoscarTransmuter)


        * [`PoscarTransmuter.from_filenames()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.PoscarTransmuter.from_filenames)


    * [`StandardTransmuter`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter)


        * [`StandardTransmuter.add_tags()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.add_tags)


        * [`StandardTransmuter.append_transformation()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.append_transformation)


        * [`StandardTransmuter.append_transformed_structures()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.append_transformed_structures)


        * [`StandardTransmuter.apply_filter()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.apply_filter)


        * [`StandardTransmuter.extend_transformations()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.extend_transformations)


        * [`StandardTransmuter.from_structures()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.from_structures)


        * [`StandardTransmuter.redo_next_change()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.redo_next_change)


        * [`StandardTransmuter.set_parameter()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.set_parameter)


        * [`StandardTransmuter.undo_last_change()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.undo_last_change)


        * [`StandardTransmuter.write_vasp_input()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.StandardTransmuter.write_vasp_input)


    * [`batch_write_vasp_input()`](pymatgen.alchemy.transmuters.md#pymatgen.alchemy.transmuters.batch_write_vasp_input)