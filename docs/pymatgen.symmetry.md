---
layout: default
title: pymatgen.symmetry.md
nav_exclude: true
---

# pymatgen.symmetry package

The symmetry package implements symmetry tools like spacegroup determination, etc.



* [pymatgen.symmetry.analyzer module](pymatgen.symmetry.analyzer.md)


    * [`PointGroupAnalyzer`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer)


        * [`PointGroupAnalyzer.get_equivalent_atoms()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_equivalent_atoms)


        * [`PointGroupAnalyzer.get_pointgroup()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_pointgroup)


        * [`PointGroupAnalyzer.get_rotational_symmetry_number()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_rotational_symmetry_number)


        * [`PointGroupAnalyzer.get_symmetry_operations()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.get_symmetry_operations)


        * [`PointGroupAnalyzer.inversion_op`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.inversion_op)


        * [`PointGroupAnalyzer.is_valid_op()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.is_valid_op)


        * [`PointGroupAnalyzer.symmetrize_molecule()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupAnalyzer.symmetrize_molecule)


    * [`PointGroupOperations`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupOperations)


        * [`PointGroupOperations.sch_symbol`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.PointGroupOperations.sch_symbol)


    * [`SpacegroupAnalyzer`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer)


        * [`SpacegroupAnalyzer.find_primitive()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.find_primitive)


        * [`SpacegroupAnalyzer.get_conventional_standard_structure()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_conventional_standard_structure)


        * [`SpacegroupAnalyzer.get_conventional_to_primitive_transformation_matrix()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_conventional_to_primitive_transformation_matrix)


        * [`SpacegroupAnalyzer.get_crystal_system()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_crystal_system)


        * [`SpacegroupAnalyzer.get_hall()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_hall)


        * [`SpacegroupAnalyzer.get_ir_reciprocal_mesh()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_ir_reciprocal_mesh)


        * [`SpacegroupAnalyzer.get_ir_reciprocal_mesh_map()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_ir_reciprocal_mesh_map)


        * [`SpacegroupAnalyzer.get_kpoint_weights()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_kpoint_weights)


        * [`SpacegroupAnalyzer.get_lattice_type()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_lattice_type)


        * [`SpacegroupAnalyzer.get_point_group_operations()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_point_group_operations)


        * [`SpacegroupAnalyzer.get_point_group_symbol()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_point_group_symbol)


        * [`SpacegroupAnalyzer.get_primitive_standard_structure()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_primitive_standard_structure)


        * [`SpacegroupAnalyzer.get_refined_structure()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_refined_structure)


        * [`SpacegroupAnalyzer.get_space_group_number()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_space_group_number)


        * [`SpacegroupAnalyzer.get_space_group_operations()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_space_group_operations)


        * [`SpacegroupAnalyzer.get_space_group_symbol()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_space_group_symbol)


        * [`SpacegroupAnalyzer.get_symmetrized_structure()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetrized_structure)


        * [`SpacegroupAnalyzer.get_symmetry_dataset()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetry_dataset)


        * [`SpacegroupAnalyzer.get_symmetry_operations()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.get_symmetry_operations)


        * [`SpacegroupAnalyzer.is_laue()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupAnalyzer.is_laue)


    * [`SpacegroupOperations`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupOperations)


        * [`SpacegroupOperations.are_symmetrically_equivalent()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.SpacegroupOperations.are_symmetrically_equivalent)


    * [`cluster_sites()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.cluster_sites)


    * [`generate_full_symmops()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.generate_full_symmops)


    * [`iterative_symmetrize()`](pymatgen.symmetry.analyzer.md#pymatgen.symmetry.analyzer.iterative_symmetrize)


* [pymatgen.symmetry.bandstructure module](pymatgen.symmetry.bandstructure.md)


    * [`HighSymmKpath`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath)


        * [`HighSymmKpath.equiv_labels`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath.equiv_labels)


        * [`HighSymmKpath.get_continuous_path()`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath.get_continuous_path)


        * [`HighSymmKpath.label_index`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath.label_index)


        * [`HighSymmKpath.path_lengths`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath.path_lengths)


        * [`HighSymmKpath.path_type`](pymatgen.symmetry.bandstructure.md#pymatgen.symmetry.bandstructure.HighSymmKpath.path_type)


* [pymatgen.symmetry.groups module](pymatgen.symmetry.groups.md)


    * [`PointGroup`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.PointGroup)


        * [`PointGroup.symbol`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.PointGroup.symbol)


        * [`PointGroup.generators`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.PointGroup.generators)


        * [`PointGroup.symmetry_ops`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.PointGroup.symmetry_ops)


    * [`SpaceGroup`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SpaceGroup)


        * [`SpaceGroup.symbol`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SpaceGroup.symbol)


        * [`SpaceGroup.int_number`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SpaceGroup.int_number)


        * [`SpaceGroup.generators`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SpaceGroup.generators)


        * [`SpaceGroup.order`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SpaceGroup.order)


    * [`SymmetryGroup`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SymmetryGroup)


        * [`SymmetryGroup.is_subgroup()`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SymmetryGroup.is_subgroup)


        * [`SymmetryGroup.is_supergroup()`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SymmetryGroup.is_supergroup)


        * [`SymmetryGroup.symmetry_ops`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SymmetryGroup.symmetry_ops)


        * [`SymmetryGroup.to_latex_string()`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.SymmetryGroup.to_latex_string)


    * [`in_array_list()`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.in_array_list)


    * [`sg_symbol_from_int_number()`](pymatgen.symmetry.groups.md#pymatgen.symmetry.groups.sg_symbol_from_int_number)


* [pymatgen.symmetry.kpath module](pymatgen.symmetry.kpath.md)


    * [`KPathBase`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase)


        * [`KPathBase.get_kpoints()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase.get_kpoints)


        * [`KPathBase.kpath`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase.kpath)


        * [`KPathBase.lattice`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase.lattice)


        * [`KPathBase.rec_lattice`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase.rec_lattice)


        * [`KPathBase.structure`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathBase.structure)


    * [`KPathLatimerMunro`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathLatimerMunro)


        * [`KPathLatimerMunro.LabelPoints()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathLatimerMunro.LabelPoints)


        * [`KPathLatimerMunro.LabelSymbol()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathLatimerMunro.LabelSymbol)


        * [`KPathLatimerMunro.mag_type`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathLatimerMunro.mag_type)


    * [`KPathSeek`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSeek)


    * [`KPathSetyawanCurtarolo`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo)


        * [`KPathSetyawanCurtarolo.bcc()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.bcc)


        * [`KPathSetyawanCurtarolo.bctet1()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.bctet1)


        * [`KPathSetyawanCurtarolo.bctet2()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.bctet2)


        * [`KPathSetyawanCurtarolo.conventional`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.conventional)


        * [`KPathSetyawanCurtarolo.cubic()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.cubic)


        * [`KPathSetyawanCurtarolo.fcc()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.fcc)


        * [`KPathSetyawanCurtarolo.hex()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.hex)


        * [`KPathSetyawanCurtarolo.mcl()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mcl)


        * [`KPathSetyawanCurtarolo.mclc1()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mclc1)


        * [`KPathSetyawanCurtarolo.mclc2()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mclc2)


        * [`KPathSetyawanCurtarolo.mclc3()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mclc3)


        * [`KPathSetyawanCurtarolo.mclc4()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mclc4)


        * [`KPathSetyawanCurtarolo.mclc5()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.mclc5)


        * [`KPathSetyawanCurtarolo.orc()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orc)


        * [`KPathSetyawanCurtarolo.orcc()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orcc)


        * [`KPathSetyawanCurtarolo.orcf1()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orcf1)


        * [`KPathSetyawanCurtarolo.orcf2()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orcf2)


        * [`KPathSetyawanCurtarolo.orcf3()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orcf3)


        * [`KPathSetyawanCurtarolo.orci()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.orci)


        * [`KPathSetyawanCurtarolo.prim`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.prim)


        * [`KPathSetyawanCurtarolo.prim_rec`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.prim_rec)


        * [`KPathSetyawanCurtarolo.rhl1()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.rhl1)


        * [`KPathSetyawanCurtarolo.rhl2()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.rhl2)


        * [`KPathSetyawanCurtarolo.tet()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.tet)


        * [`KPathSetyawanCurtarolo.tria()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.tria)


        * [`KPathSetyawanCurtarolo.trib()`](pymatgen.symmetry.kpath.md#pymatgen.symmetry.kpath.KPathSetyawanCurtarolo.trib)


* [pymatgen.symmetry.maggroups module](pymatgen.symmetry.maggroups.md)


    * [`MagneticSpaceGroup`](pymatgen.symmetry.maggroups.md#pymatgen.symmetry.maggroups.MagneticSpaceGroup)


* [pymatgen.symmetry.settings module](pymatgen.symmetry.settings.md)


    * [`JonesFaithfulTransformation`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation)


        * [`JonesFaithfulTransformation.P`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.P)


        * [`JonesFaithfulTransformation.from_origin_shift()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.from_origin_shift)


        * [`JonesFaithfulTransformation.from_transformation_string()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.from_transformation_string)


        * [`JonesFaithfulTransformation.inverse`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.inverse)


        * [`JonesFaithfulTransformation.p`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.p)


        * [`JonesFaithfulTransformation.parse_transformation_string()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.parse_transformation_string)


        * [`JonesFaithfulTransformation.transform_coords()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.transform_coords)


        * [`JonesFaithfulTransformation.transform_lattice()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.transform_lattice)


        * [`JonesFaithfulTransformation.transform_symmop()`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.transform_symmop)


        * [`JonesFaithfulTransformation.transformation_string`](pymatgen.symmetry.settings.md#pymatgen.symmetry.settings.JonesFaithfulTransformation.transformation_string)


* [pymatgen.symmetry.site_symmetries module](pymatgen.symmetry.site_symmetries.md)


    * [`get_shared_symmetry_operations()`](pymatgen.symmetry.site_symmetries.md#pymatgen.symmetry.site_symmetries.get_shared_symmetry_operations)


    * [`get_site_symmetries()`](pymatgen.symmetry.site_symmetries.md#pymatgen.symmetry.site_symmetries.get_site_symmetries)


* [pymatgen.symmetry.structure module](pymatgen.symmetry.structure.md)


    * [`SymmetrizedStructure`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure)


        * [`SymmetrizedStructure.as_dict()`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure.as_dict)


        * [`SymmetrizedStructure.copy()`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure.copy)


        * [`SymmetrizedStructure.find_equivalent_sites()`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure.find_equivalent_sites)


        * [`SymmetrizedStructure.from_dict()`](pymatgen.symmetry.structure.md#pymatgen.symmetry.structure.SymmetrizedStructure.from_dict)