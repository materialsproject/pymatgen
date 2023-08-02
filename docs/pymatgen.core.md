---
layout: default
title: pymatgen.core.md
nav_exclude: true
---

# pymatgen.core package

This package contains core modules and classes for representing structures and operations on them.



* [pymatgen.core.bonds module](pymatgen.core.bonds.md)


    * [`CovalentBond`](pymatgen.core.bonds.md#pymatgen.core.bonds.CovalentBond)


        * [`CovalentBond.get_bond_order()`](pymatgen.core.bonds.md#pymatgen.core.bonds.CovalentBond.get_bond_order)


        * [`CovalentBond.is_bonded()`](pymatgen.core.bonds.md#pymatgen.core.bonds.CovalentBond.is_bonded)


        * [`CovalentBond.length`](pymatgen.core.bonds.md#pymatgen.core.bonds.CovalentBond.length)


    * [`get_bond_length()`](pymatgen.core.bonds.md#pymatgen.core.bonds.get_bond_length)


    * [`get_bond_order()`](pymatgen.core.bonds.md#pymatgen.core.bonds.get_bond_order)


    * [`obtain_all_bond_lengths()`](pymatgen.core.bonds.md#pymatgen.core.bonds.obtain_all_bond_lengths)


* [pymatgen.core.composition module](pymatgen.core.composition.md)


    * [`ChemicalPotential`](pymatgen.core.composition.md#pymatgen.core.composition.ChemicalPotential)


        * [`ChemicalPotential.get_energy()`](pymatgen.core.composition.md#pymatgen.core.composition.ChemicalPotential.get_energy)


    * [`Composition`](pymatgen.core.composition.md#pymatgen.core.composition.Composition)


        * [`Composition.add_charges_from_oxi_state_guesses()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.add_charges_from_oxi_state_guesses)


        * [`Composition.almost_equals()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.almost_equals)


        * [`Composition.alphabetical_formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.alphabetical_formula)


        * [`Composition.amount_tolerance`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.amount_tolerance)


        * [`Composition.anonymized_formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.anonymized_formula)


        * [`Composition.as_dict()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.as_dict)


        * [`Composition.average_electroneg`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.average_electroneg)


        * [`Composition.chemical_system`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.chemical_system)


        * [`Composition.contains_element_type()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.contains_element_type)


        * [`Composition.copy()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.copy)


        * [`Composition.element_composition`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.element_composition)


        * [`Composition.elements`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.elements)


        * [`Composition.formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.formula)


        * [`Composition.fractional_composition`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.fractional_composition)


        * [`Composition.from_dict()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.from_dict)


        * [`Composition.from_weight_dict()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.from_weight_dict)


        * [`Composition.get_atomic_fraction()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_atomic_fraction)


        * [`Composition.get_el_amt_dict()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_el_amt_dict)


        * [`Composition.get_integer_formula_and_factor()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_integer_formula_and_factor)


        * [`Composition.get_reduced_composition_and_factor()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_reduced_composition_and_factor)


        * [`Composition.get_reduced_formula_and_factor()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_reduced_formula_and_factor)


        * [`Composition.get_wt_fraction()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.get_wt_fraction)


        * [`Composition.hill_formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.hill_formula)


        * [`Composition.is_element`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.is_element)


        * [`Composition.iupac_formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.iupac_formula)


        * [`Composition.num_atoms`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.num_atoms)


        * [`Composition.oxi_prob`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.oxi_prob)


        * [`Composition.oxi_state_guesses()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.oxi_state_guesses)


        * [`Composition.ranked_compositions_from_indeterminate_formula()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.ranked_compositions_from_indeterminate_formula)


        * [`Composition.reduced_composition`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.reduced_composition)


        * [`Composition.reduced_formula`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.reduced_formula)


        * [`Composition.remove_charges()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.remove_charges)


        * [`Composition.replace()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.replace)


        * [`Composition.special_formulas`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.special_formulas)


        * [`Composition.to_data_dict`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.to_data_dict)


        * [`Composition.to_pretty_string()`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.to_pretty_string)


        * [`Composition.to_reduced_dict`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.to_reduced_dict)


        * [`Composition.to_weight_dict`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.to_weight_dict)


        * [`Composition.total_electrons`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.total_electrons)


        * [`Composition.valid`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.valid)


        * [`Composition.weight`](pymatgen.core.composition.md#pymatgen.core.composition.Composition.weight)


    * [`CompositionError`](pymatgen.core.composition.md#pymatgen.core.composition.CompositionError)


    * [`reduce_formula()`](pymatgen.core.composition.md#pymatgen.core.composition.reduce_formula)


* [pymatgen.core.interface module](pymatgen.core.interface.md)


    * [`Interface`](pymatgen.core.interface.md#pymatgen.core.interface.Interface)


        * [`Interface.as_dict()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.as_dict)


        * [`Interface.copy()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.copy)


        * [`Interface.film`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.film)


        * [`Interface.film_indices`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.film_indices)


        * [`Interface.film_layers`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.film_layers)


        * [`Interface.film_sites`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.film_sites)


        * [`Interface.film_termination`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.film_termination)


        * [`Interface.from_dict()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.from_dict)


        * [`Interface.from_slabs()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.from_slabs)


        * [`Interface.gap`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.gap)


        * [`Interface.get_shifts_based_on_adsorbate_sites()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.get_shifts_based_on_adsorbate_sites)


        * [`Interface.get_sorted_structure()`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.get_sorted_structure)


        * [`Interface.in_plane_offset`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.in_plane_offset)


        * [`Interface.substrate`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.substrate)


        * [`Interface.substrate_indices`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.substrate_indices)


        * [`Interface.substrate_layers`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.substrate_layers)


        * [`Interface.substrate_sites`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.substrate_sites)


        * [`Interface.substrate_termination`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.substrate_termination)


        * [`Interface.vacuum_over_film`](pymatgen.core.interface.md#pymatgen.core.interface.Interface.vacuum_over_film)


    * [`count_layers()`](pymatgen.core.interface.md#pymatgen.core.interface.count_layers)


    * [`label_termination()`](pymatgen.core.interface.md#pymatgen.core.interface.label_termination)


* [pymatgen.core.ion module](pymatgen.core.ion.md)


    * [`Ion`](pymatgen.core.ion.md#pymatgen.core.ion.Ion)


        * [`Ion.alphabetical_formula`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.alphabetical_formula)


        * [`Ion.anonymized_formula`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.anonymized_formula)


        * [`Ion.as_dict()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.as_dict)


        * [`Ion.charge`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.charge)


        * [`Ion.composition`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.composition)


        * [`Ion.formula`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.formula)


        * [`Ion.from_dict()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.from_dict)


        * [`Ion.from_formula()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.from_formula)


        * [`Ion.get_reduced_formula_and_factor()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.get_reduced_formula_and_factor)


        * [`Ion.oxi_state_guesses()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.oxi_state_guesses)


        * [`Ion.reduced_formula`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.reduced_formula)


        * [`Ion.to_pretty_string()`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.to_pretty_string)


        * [`Ion.to_reduced_dict`](pymatgen.core.ion.md#pymatgen.core.ion.Ion.to_reduced_dict)


* [pymatgen.core.lattice module](pymatgen.core.lattice.md)


    * [`Lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice)


        * [`Lattice.a`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.a)


        * [`Lattice.abc`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.abc)


        * [`Lattice.alpha`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.alpha)


        * [`Lattice.angles`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.angles)


        * [`Lattice.as_dict()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.as_dict)


        * [`Lattice.b`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.b)


        * [`Lattice.beta`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.beta)


        * [`Lattice.c`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.c)


        * [`Lattice.copy()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.copy)


        * [`Lattice.cubic()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.cubic)


        * [`Lattice.d_hkl()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.d_hkl)


        * [`Lattice.dot()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.dot)


        * [`Lattice.find_all_mappings()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.find_all_mappings)


        * [`Lattice.find_mapping()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.find_mapping)


        * [`Lattice.from_dict()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.from_dict)


        * [`Lattice.from_parameters()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.from_parameters)


        * [`Lattice.gamma`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.gamma)


        * [`Lattice.get_all_distances()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_all_distances)


        * [`Lattice.get_brillouin_zone()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_brillouin_zone)


        * [`Lattice.get_cartesian_coords()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_cartesian_coords)


        * [`Lattice.get_distance_and_image()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_distance_and_image)


        * [`Lattice.get_frac_coords_from_lll()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_frac_coords_from_lll)


        * [`Lattice.get_fractional_coords()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_fractional_coords)


        * [`Lattice.get_lll_frac_coords()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_lll_frac_coords)


        * [`Lattice.get_lll_reduced_lattice()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_lll_reduced_lattice)


        * [`Lattice.get_miller_index_from_coords()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_miller_index_from_coords)


        * [`Lattice.get_niggli_reduced_lattice()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_niggli_reduced_lattice)


        * [`Lattice.get_points_in_sphere()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_points_in_sphere)


        * [`Lattice.get_points_in_sphere_old()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_points_in_sphere_old)


        * [`Lattice.get_points_in_sphere_py()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_points_in_sphere_py)


        * [`Lattice.get_recp_symmetry_operation()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_recp_symmetry_operation)


        * [`Lattice.get_vector_along_lattice_directions()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_vector_along_lattice_directions)


        * [`Lattice.get_wigner_seitz_cell()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.get_wigner_seitz_cell)


        * [`Lattice.hexagonal()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.hexagonal)


        * [`Lattice.inv_matrix`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.inv_matrix)


        * [`Lattice.is_3d_periodic`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.is_3d_periodic)


        * [`Lattice.is_hexagonal()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.is_hexagonal)


        * [`Lattice.is_orthogonal`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.is_orthogonal)


        * [`Lattice.lengths`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.lengths)


        * [`Lattice.lll_inverse`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.lll_inverse)


        * [`Lattice.lll_mapping`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.lll_mapping)


        * [`Lattice.lll_matrix`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.lll_matrix)


        * [`Lattice.matrix`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.matrix)


        * [`Lattice.metric_tensor`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.metric_tensor)


        * [`Lattice.monoclinic()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.monoclinic)


        * [`Lattice.norm()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.norm)


        * [`Lattice.orthorhombic()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.orthorhombic)


        * [`Lattice.parameters`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.parameters)


        * [`Lattice.pbc`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.pbc)


        * [`Lattice.reciprocal_lattice`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.reciprocal_lattice)


        * [`Lattice.reciprocal_lattice_crystallographic`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.reciprocal_lattice_crystallographic)


        * [`Lattice.rhombohedral()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.rhombohedral)


        * [`Lattice.scale()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.scale)


        * [`Lattice.selling_dist()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.selling_dist)


        * [`Lattice.selling_vector`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.selling_vector)


        * [`Lattice.tetragonal()`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.tetragonal)


        * [`Lattice.volume`](pymatgen.core.lattice.md#pymatgen.core.lattice.Lattice.volume)


    * [`find_neighbors()`](pymatgen.core.lattice.md#pymatgen.core.lattice.find_neighbors)


    * [`get_integer_index()`](pymatgen.core.lattice.md#pymatgen.core.lattice.get_integer_index)


    * [`get_points_in_spheres()`](pymatgen.core.lattice.md#pymatgen.core.lattice.get_points_in_spheres)


* [pymatgen.core.libxcfunc module](pymatgen.core.libxcfunc.md)


    * [`LibxcFunc`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc)


        * [`LibxcFunc.GGA_C_AM05`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_AM05)


        * [`LibxcFunc.GGA_C_APBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_APBE)


        * [`LibxcFunc.GGA_C_BGCP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_BGCP)


        * [`LibxcFunc.GGA_C_FT97`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_FT97)


        * [`LibxcFunc.GGA_C_GAM`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_GAM)


        * [`LibxcFunc.GGA_C_HCTH_A`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_HCTH_A)


        * [`LibxcFunc.GGA_C_LM`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_LM)


        * [`LibxcFunc.GGA_C_LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_LYP)


        * [`LibxcFunc.GGA_C_N12`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_N12)


        * [`LibxcFunc.GGA_C_N12_SX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_N12_SX)


        * [`LibxcFunc.GGA_C_OPTC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OPTC)


        * [`LibxcFunc.GGA_C_OP_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OP_B88)


        * [`LibxcFunc.GGA_C_OP_G96`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OP_G96)


        * [`LibxcFunc.GGA_C_OP_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OP_PBE)


        * [`LibxcFunc.GGA_C_OP_PW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OP_PW91)


        * [`LibxcFunc.GGA_C_OP_XALPHA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_OP_XALPHA)


        * [`LibxcFunc.GGA_C_P86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_P86)


        * [`LibxcFunc.GGA_C_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBE)


        * [`LibxcFunc.GGA_C_PBEFE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBEFE)


        * [`LibxcFunc.GGA_C_PBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBEINT)


        * [`LibxcFunc.GGA_C_PBELOC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBELOC)


        * [`LibxcFunc.GGA_C_PBE_JRGX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBE_JRGX)


        * [`LibxcFunc.GGA_C_PBE_SOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PBE_SOL)


        * [`LibxcFunc.GGA_C_PW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_PW91)


        * [`LibxcFunc.GGA_C_Q2D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_Q2D)


        * [`LibxcFunc.GGA_C_REGTPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_REGTPSS)


        * [`LibxcFunc.GGA_C_REVTCA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_REVTCA)


        * [`LibxcFunc.GGA_C_RGE2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_RGE2)


        * [`LibxcFunc.GGA_C_SOGGA11`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_SOGGA11)


        * [`LibxcFunc.GGA_C_SOGGA11_X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_SOGGA11_X)


        * [`LibxcFunc.GGA_C_SPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_SPBE)


        * [`LibxcFunc.GGA_C_TCA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_TCA)


        * [`LibxcFunc.GGA_C_WI`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_WI)


        * [`LibxcFunc.GGA_C_WI0`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_WI0)


        * [`LibxcFunc.GGA_C_WL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_WL)


        * [`LibxcFunc.GGA_C_XPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_XPBE)


        * [`LibxcFunc.GGA_C_ZPBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_ZPBEINT)


        * [`LibxcFunc.GGA_C_ZPBESOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_C_ZPBESOL)


        * [`LibxcFunc.GGA_K_ABSP1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_ABSP1)


        * [`LibxcFunc.GGA_K_ABSP2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_ABSP2)


        * [`LibxcFunc.GGA_K_APBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_APBE)


        * [`LibxcFunc.GGA_K_APBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_APBEINT)


        * [`LibxcFunc.GGA_K_BALTIN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_BALTIN)


        * [`LibxcFunc.GGA_K_DK`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_DK)


        * [`LibxcFunc.GGA_K_ERNZERHOF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_ERNZERHOF)


        * [`LibxcFunc.GGA_K_FR_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_FR_B88)


        * [`LibxcFunc.GGA_K_FR_PW86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_FR_PW86)


        * [`LibxcFunc.GGA_K_GE2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_GE2)


        * [`LibxcFunc.GGA_K_GOLDEN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_GOLDEN)


        * [`LibxcFunc.GGA_K_GP85`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_GP85)


        * [`LibxcFunc.GGA_K_GR`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_GR)


        * [`LibxcFunc.GGA_K_LC94`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_LC94)


        * [`LibxcFunc.GGA_K_LIEB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_LIEB)


        * [`LibxcFunc.GGA_K_LLP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_LLP)


        * [`LibxcFunc.GGA_K_LUDENA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_LUDENA)


        * [`LibxcFunc.GGA_K_MEYER`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_MEYER)


        * [`LibxcFunc.GGA_K_OL1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_OL1)


        * [`LibxcFunc.GGA_K_OL2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_OL2)


        * [`LibxcFunc.GGA_K_PEARSON`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_PEARSON)


        * [`LibxcFunc.GGA_K_PERDEW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_PERDEW)


        * [`LibxcFunc.GGA_K_REVAPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_REVAPBE)


        * [`LibxcFunc.GGA_K_REVAPBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_REVAPBEINT)


        * [`LibxcFunc.GGA_K_TFVW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_TFVW)


        * [`LibxcFunc.GGA_K_THAKKAR`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_THAKKAR)


        * [`LibxcFunc.GGA_K_TW1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_TW1)


        * [`LibxcFunc.GGA_K_TW2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_TW2)


        * [`LibxcFunc.GGA_K_TW3`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_TW3)


        * [`LibxcFunc.GGA_K_TW4`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_TW4)


        * [`LibxcFunc.GGA_K_VJKS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_VJKS)


        * [`LibxcFunc.GGA_K_VSK`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_VSK)


        * [`LibxcFunc.GGA_K_VW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_VW)


        * [`LibxcFunc.GGA_K_YT65`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_K_YT65)


        * [`LibxcFunc.GGA_XC_B97_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_B97_D)


        * [`LibxcFunc.GGA_XC_B97_GGA1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_B97_GGA1)


        * [`LibxcFunc.GGA_XC_EDF1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_EDF1)


        * [`LibxcFunc.GGA_XC_HCTH_120`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_120)


        * [`LibxcFunc.GGA_XC_HCTH_147`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_147)


        * [`LibxcFunc.GGA_XC_HCTH_407`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_407)


        * [`LibxcFunc.GGA_XC_HCTH_407P`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_407P)


        * [`LibxcFunc.GGA_XC_HCTH_93`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_93)


        * [`LibxcFunc.GGA_XC_HCTH_P14`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_P14)


        * [`LibxcFunc.GGA_XC_HCTH_P76`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_HCTH_P76)


        * [`LibxcFunc.GGA_XC_KT2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_KT2)


        * [`LibxcFunc.GGA_XC_MOHLYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_MOHLYP)


        * [`LibxcFunc.GGA_XC_MOHLYP2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_MOHLYP2)


        * [`LibxcFunc.GGA_XC_MPWLYP1W`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_MPWLYP1W)


        * [`LibxcFunc.GGA_XC_OBLYP_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_OBLYP_D)


        * [`LibxcFunc.GGA_XC_OPBE_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_OPBE_D)


        * [`LibxcFunc.GGA_XC_OPWLYP_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_OPWLYP_D)


        * [`LibxcFunc.GGA_XC_PBE1W`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_PBE1W)


        * [`LibxcFunc.GGA_XC_PBELYP1W`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_PBELYP1W)


        * [`LibxcFunc.GGA_XC_TH1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH1)


        * [`LibxcFunc.GGA_XC_TH2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH2)


        * [`LibxcFunc.GGA_XC_TH3`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH3)


        * [`LibxcFunc.GGA_XC_TH4`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH4)


        * [`LibxcFunc.GGA_XC_TH_FC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH_FC)


        * [`LibxcFunc.GGA_XC_TH_FCFO`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH_FCFO)


        * [`LibxcFunc.GGA_XC_TH_FCO`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH_FCO)


        * [`LibxcFunc.GGA_XC_TH_FL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_TH_FL)


        * [`LibxcFunc.GGA_XC_VV10`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_VV10)


        * [`LibxcFunc.GGA_XC_XLYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_XC_XLYP)


        * [`LibxcFunc.GGA_X_2D_B86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_2D_B86)


        * [`LibxcFunc.GGA_X_2D_B86_MGC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_2D_B86_MGC)


        * [`LibxcFunc.GGA_X_2D_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_2D_B88)


        * [`LibxcFunc.GGA_X_2D_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_2D_PBE)


        * [`LibxcFunc.GGA_X_AIRY`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_AIRY)


        * [`LibxcFunc.GGA_X_AK13`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_AK13)


        * [`LibxcFunc.GGA_X_AM05`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_AM05)


        * [`LibxcFunc.GGA_X_APBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_APBE)


        * [`LibxcFunc.GGA_X_B86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_B86)


        * [`LibxcFunc.GGA_X_B86_MGC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_B86_MGC)


        * [`LibxcFunc.GGA_X_B86_R`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_B86_R)


        * [`LibxcFunc.GGA_X_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_B88)


        * [`LibxcFunc.GGA_X_BAYESIAN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_BAYESIAN)


        * [`LibxcFunc.GGA_X_BGCP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_BGCP)


        * [`LibxcFunc.GGA_X_BPCCAC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_BPCCAC)


        * [`LibxcFunc.GGA_X_C09X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_C09X)


        * [`LibxcFunc.GGA_X_CAP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_CAP)


        * [`LibxcFunc.GGA_X_DK87_R1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_DK87_R1)


        * [`LibxcFunc.GGA_X_DK87_R2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_DK87_R2)


        * [`LibxcFunc.GGA_X_EV93`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_EV93)


        * [`LibxcFunc.GGA_X_FT97_A`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_FT97_A)


        * [`LibxcFunc.GGA_X_FT97_B`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_FT97_B)


        * [`LibxcFunc.GGA_X_G96`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_G96)


        * [`LibxcFunc.GGA_X_GAM`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_GAM)


        * [`LibxcFunc.GGA_X_HCTH_A`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HCTH_A)


        * [`LibxcFunc.GGA_X_HERMAN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HERMAN)


        * [`LibxcFunc.GGA_X_HJS_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HJS_B88)


        * [`LibxcFunc.GGA_X_HJS_B88_V2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HJS_B88_V2)


        * [`LibxcFunc.GGA_X_HJS_B97X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HJS_B97X)


        * [`LibxcFunc.GGA_X_HJS_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HJS_PBE)


        * [`LibxcFunc.GGA_X_HJS_PBE_SOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HJS_PBE_SOL)


        * [`LibxcFunc.GGA_X_HTBS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_HTBS)


        * [`LibxcFunc.GGA_X_ITYH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_ITYH)


        * [`LibxcFunc.GGA_X_KT1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_KT1)


        * [`LibxcFunc.GGA_X_LAG`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LAG)


        * [`LibxcFunc.GGA_X_LAMBDA_CH_N`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LAMBDA_CH_N)


        * [`LibxcFunc.GGA_X_LAMBDA_LO_N`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LAMBDA_LO_N)


        * [`LibxcFunc.GGA_X_LAMBDA_OC2_N`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LAMBDA_OC2_N)


        * [`LibxcFunc.GGA_X_LB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LB)


        * [`LibxcFunc.GGA_X_LBM`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LBM)


        * [`LibxcFunc.GGA_X_LG93`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LG93)


        * [`LibxcFunc.GGA_X_LV_RPW86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_LV_RPW86)


        * [`LibxcFunc.GGA_X_MB88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_MB88)


        * [`LibxcFunc.GGA_X_MPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_MPBE)


        * [`LibxcFunc.GGA_X_MPW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_MPW91)


        * [`LibxcFunc.GGA_X_N12`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_N12)


        * [`LibxcFunc.GGA_X_OL2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_OL2)


        * [`LibxcFunc.GGA_X_OPTB88_VDW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_OPTB88_VDW)


        * [`LibxcFunc.GGA_X_OPTPBE_VDW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_OPTPBE_VDW)


        * [`LibxcFunc.GGA_X_OPTX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_OPTX)


        * [`LibxcFunc.GGA_X_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE)


        * [`LibxcFunc.GGA_X_PBEA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBEA)


        * [`LibxcFunc.GGA_X_PBEFE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBEFE)


        * [`LibxcFunc.GGA_X_PBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBEINT)


        * [`LibxcFunc.GGA_X_PBEK1_VDW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBEK1_VDW)


        * [`LibxcFunc.GGA_X_PBE_JSJR`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE_JSJR)


        * [`LibxcFunc.GGA_X_PBE_MOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE_MOL)


        * [`LibxcFunc.GGA_X_PBE_R`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE_R)


        * [`LibxcFunc.GGA_X_PBE_SOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE_SOL)


        * [`LibxcFunc.GGA_X_PBE_TCA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PBE_TCA)


        * [`LibxcFunc.GGA_X_PW86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PW86)


        * [`LibxcFunc.GGA_X_PW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_PW91)


        * [`LibxcFunc.GGA_X_Q2D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_Q2D)


        * [`LibxcFunc.GGA_X_RGE2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_RGE2)


        * [`LibxcFunc.GGA_X_RPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_RPBE)


        * [`LibxcFunc.GGA_X_RPW86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_RPW86)


        * [`LibxcFunc.GGA_X_SFAT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SFAT)


        * [`LibxcFunc.GGA_X_SOGGA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SOGGA)


        * [`LibxcFunc.GGA_X_SOGGA11`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SOGGA11)


        * [`LibxcFunc.GGA_X_SSB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SSB)


        * [`LibxcFunc.GGA_X_SSB_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SSB_D)


        * [`LibxcFunc.GGA_X_SSB_SW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_SSB_SW)


        * [`LibxcFunc.GGA_X_VMT84_GE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_VMT84_GE)


        * [`LibxcFunc.GGA_X_VMT84_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_VMT84_PBE)


        * [`LibxcFunc.GGA_X_VMT_GE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_VMT_GE)


        * [`LibxcFunc.GGA_X_VMT_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_VMT_PBE)


        * [`LibxcFunc.GGA_X_WC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_WC)


        * [`LibxcFunc.GGA_X_WPBEH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_WPBEH)


        * [`LibxcFunc.GGA_X_XPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.GGA_X_XPBE)


        * [`LibxcFunc.HYB_GGA_XC_B1LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B1LYP)


        * [`LibxcFunc.HYB_GGA_XC_B1PW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B1PW91)


        * [`LibxcFunc.HYB_GGA_XC_B1WC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B1WC)


        * [`LibxcFunc.HYB_GGA_XC_B3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B3LYP)


        * [`LibxcFunc.HYB_GGA_XC_B3LYP5`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B3LYP5)


        * [`LibxcFunc.HYB_GGA_XC_B3LYPs`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B3LYPs)


        * [`LibxcFunc.HYB_GGA_XC_B3P86`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B3P86)


        * [`LibxcFunc.HYB_GGA_XC_B3PW91`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B3PW91)


        * [`LibxcFunc.HYB_GGA_XC_B97`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97)


        * [`LibxcFunc.HYB_GGA_XC_B97_1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97_1)


        * [`LibxcFunc.HYB_GGA_XC_B97_1p`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97_1p)


        * [`LibxcFunc.HYB_GGA_XC_B97_2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97_2)


        * [`LibxcFunc.HYB_GGA_XC_B97_3`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97_3)


        * [`LibxcFunc.HYB_GGA_XC_B97_K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_B97_K)


        * [`LibxcFunc.HYB_GGA_XC_BHANDH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_BHANDH)


        * [`LibxcFunc.HYB_GGA_XC_BHANDHLYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_BHANDHLYP)


        * [`LibxcFunc.HYB_GGA_XC_CAMY_B3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_CAMY_B3LYP)


        * [`LibxcFunc.HYB_GGA_XC_CAMY_BLYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_CAMY_BLYP)


        * [`LibxcFunc.HYB_GGA_XC_CAM_B3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_CAM_B3LYP)


        * [`LibxcFunc.HYB_GGA_XC_CAP0`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_CAP0)


        * [`LibxcFunc.HYB_GGA_XC_EDF2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_EDF2)


        * [`LibxcFunc.HYB_GGA_XC_HJS_B88`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HJS_B88)


        * [`LibxcFunc.HYB_GGA_XC_HJS_B97X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HJS_B97X)


        * [`LibxcFunc.HYB_GGA_XC_HJS_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HJS_PBE)


        * [`LibxcFunc.HYB_GGA_XC_HJS_PBE_SOL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HJS_PBE_SOL)


        * [`LibxcFunc.HYB_GGA_XC_HPBEINT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HPBEINT)


        * [`LibxcFunc.HYB_GGA_XC_HSE03`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HSE03)


        * [`LibxcFunc.HYB_GGA_XC_HSE06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_HSE06)


        * [`LibxcFunc.HYB_GGA_XC_LCY_BLYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_LCY_BLYP)


        * [`LibxcFunc.HYB_GGA_XC_LCY_PBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_LCY_PBE)


        * [`LibxcFunc.HYB_GGA_XC_LC_VV10`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_LC_VV10)


        * [`LibxcFunc.HYB_GGA_XC_LRC_WPBE`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_LRC_WPBE)


        * [`LibxcFunc.HYB_GGA_XC_LRC_WPBEH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_LRC_WPBEH)


        * [`LibxcFunc.HYB_GGA_XC_MB3LYP_RC04`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_MB3LYP_RC04)


        * [`LibxcFunc.HYB_GGA_XC_MPW3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_MPW3LYP)


        * [`LibxcFunc.HYB_GGA_XC_MPW3PW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_MPW3PW)


        * [`LibxcFunc.HYB_GGA_XC_MPWLYP1M`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_MPWLYP1M)


        * [`LibxcFunc.HYB_GGA_XC_O3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_O3LYP)


        * [`LibxcFunc.HYB_GGA_XC_PBE0_13`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_PBE0_13)


        * [`LibxcFunc.HYB_GGA_XC_PBEH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_PBEH)


        * [`LibxcFunc.HYB_GGA_XC_REVB3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_REVB3LYP)


        * [`LibxcFunc.HYB_GGA_XC_SB98_1a`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_1a)


        * [`LibxcFunc.HYB_GGA_XC_SB98_1b`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_1b)


        * [`LibxcFunc.HYB_GGA_XC_SB98_1c`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_1c)


        * [`LibxcFunc.HYB_GGA_XC_SB98_2a`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_2a)


        * [`LibxcFunc.HYB_GGA_XC_SB98_2b`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_2b)


        * [`LibxcFunc.HYB_GGA_XC_SB98_2c`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_SB98_2c)


        * [`LibxcFunc.HYB_GGA_XC_TUNED_CAM_B3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_TUNED_CAM_B3LYP)


        * [`LibxcFunc.HYB_GGA_XC_WB97`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_WB97)


        * [`LibxcFunc.HYB_GGA_XC_WB97X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_WB97X)


        * [`LibxcFunc.HYB_GGA_XC_WB97X_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_WB97X_D)


        * [`LibxcFunc.HYB_GGA_XC_WB97X_V`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_WB97X_V)


        * [`LibxcFunc.HYB_GGA_XC_X3LYP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_X3LYP)


        * [`LibxcFunc.HYB_GGA_XC_mPW1K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_mPW1K)


        * [`LibxcFunc.HYB_GGA_XC_mPW1PW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_XC_mPW1PW)


        * [`LibxcFunc.HYB_GGA_X_N12_SX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_X_N12_SX)


        * [`LibxcFunc.HYB_GGA_X_SOGGA11_X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_GGA_X_SOGGA11_X)


        * [`LibxcFunc.HYB_MGGA_XC_B86B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_B86B95)


        * [`LibxcFunc.HYB_MGGA_XC_B88B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_B88B95)


        * [`LibxcFunc.HYB_MGGA_XC_BB1K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_BB1K)


        * [`LibxcFunc.HYB_MGGA_XC_M05`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M05)


        * [`LibxcFunc.HYB_MGGA_XC_M05_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M05_2X)


        * [`LibxcFunc.HYB_MGGA_XC_M06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M06)


        * [`LibxcFunc.HYB_MGGA_XC_M06_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M06_2X)


        * [`LibxcFunc.HYB_MGGA_XC_M06_HF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M06_HF)


        * [`LibxcFunc.HYB_MGGA_XC_M08_HX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M08_HX)


        * [`LibxcFunc.HYB_MGGA_XC_M08_SO`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M08_SO)


        * [`LibxcFunc.HYB_MGGA_XC_M11`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_M11)


        * [`LibxcFunc.HYB_MGGA_XC_MPW1B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_MPW1B95)


        * [`LibxcFunc.HYB_MGGA_XC_MPWB1K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_MPWB1K)


        * [`LibxcFunc.HYB_MGGA_XC_PW6B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_PW6B95)


        * [`LibxcFunc.HYB_MGGA_XC_PW86B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_PW86B95)


        * [`LibxcFunc.HYB_MGGA_XC_PWB6K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_PWB6K)


        * [`LibxcFunc.HYB_MGGA_XC_REVTPSSH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_REVTPSSH)


        * [`LibxcFunc.HYB_MGGA_XC_TPSSH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_TPSSH)


        * [`LibxcFunc.HYB_MGGA_XC_WB97M_V`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_WB97M_V)


        * [`LibxcFunc.HYB_MGGA_XC_X1B95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_X1B95)


        * [`LibxcFunc.HYB_MGGA_XC_XB1K`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_XC_XB1K)


        * [`LibxcFunc.HYB_MGGA_X_DLDF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_DLDF)


        * [`LibxcFunc.HYB_MGGA_X_MN12_SX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_MN12_SX)


        * [`LibxcFunc.HYB_MGGA_X_MN15`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_MN15)


        * [`LibxcFunc.HYB_MGGA_X_MS2H`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_MS2H)


        * [`LibxcFunc.HYB_MGGA_X_MVSH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_MVSH)


        * [`LibxcFunc.HYB_MGGA_X_SCAN0`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.HYB_MGGA_X_SCAN0)


        * [`LibxcFunc.LDA_C_1D_CSC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_1D_CSC)


        * [`LibxcFunc.LDA_C_1D_LOOS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_1D_LOOS)


        * [`LibxcFunc.LDA_C_2D_AMGB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_2D_AMGB)


        * [`LibxcFunc.LDA_C_2D_PRM`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_2D_PRM)


        * [`LibxcFunc.LDA_C_GL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_GL)


        * [`LibxcFunc.LDA_C_GOMBAS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_GOMBAS)


        * [`LibxcFunc.LDA_C_HL`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_HL)


        * [`LibxcFunc.LDA_C_ML1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_ML1)


        * [`LibxcFunc.LDA_C_ML2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_ML2)


        * [`LibxcFunc.LDA_C_OB_PW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_OB_PW)


        * [`LibxcFunc.LDA_C_OB_PZ`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_OB_PZ)


        * [`LibxcFunc.LDA_C_PW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_PW)


        * [`LibxcFunc.LDA_C_PW_MOD`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_PW_MOD)


        * [`LibxcFunc.LDA_C_PW_RPA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_PW_RPA)


        * [`LibxcFunc.LDA_C_PZ`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_PZ)


        * [`LibxcFunc.LDA_C_PZ_MOD`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_PZ_MOD)


        * [`LibxcFunc.LDA_C_RC04`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_RC04)


        * [`LibxcFunc.LDA_C_RPA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_RPA)


        * [`LibxcFunc.LDA_C_VWN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN)


        * [`LibxcFunc.LDA_C_VWN_1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN_1)


        * [`LibxcFunc.LDA_C_VWN_2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN_2)


        * [`LibxcFunc.LDA_C_VWN_3`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN_3)


        * [`LibxcFunc.LDA_C_VWN_4`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN_4)


        * [`LibxcFunc.LDA_C_VWN_RPA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_VWN_RPA)


        * [`LibxcFunc.LDA_C_WIGNER`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_WIGNER)


        * [`LibxcFunc.LDA_C_XALPHA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_XALPHA)


        * [`LibxcFunc.LDA_C_vBH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_C_vBH)


        * [`LibxcFunc.LDA_K_LP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_K_LP)


        * [`LibxcFunc.LDA_K_TF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_K_TF)


        * [`LibxcFunc.LDA_X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_X)


        * [`LibxcFunc.LDA_XC_KSDT`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_XC_KSDT)


        * [`LibxcFunc.LDA_XC_TETER93`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_XC_TETER93)


        * [`LibxcFunc.LDA_XC_ZLP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_XC_ZLP)


        * [`LibxcFunc.LDA_X_1D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_X_1D)


        * [`LibxcFunc.LDA_X_2D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.LDA_X_2D)


        * [`LibxcFunc.MGGA_C_BC95`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_BC95)


        * [`LibxcFunc.MGGA_C_CC06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_CC06)


        * [`LibxcFunc.MGGA_C_CS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_CS)


        * [`LibxcFunc.MGGA_C_DLDF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_DLDF)


        * [`LibxcFunc.MGGA_C_M05`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M05)


        * [`LibxcFunc.MGGA_C_M05_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M05_2X)


        * [`LibxcFunc.MGGA_C_M06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M06)


        * [`LibxcFunc.MGGA_C_M06_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M06_2X)


        * [`LibxcFunc.MGGA_C_M06_HF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M06_HF)


        * [`LibxcFunc.MGGA_C_M06_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M06_L)


        * [`LibxcFunc.MGGA_C_M08_HX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M08_HX)


        * [`LibxcFunc.MGGA_C_M08_SO`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M08_SO)


        * [`LibxcFunc.MGGA_C_M11`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M11)


        * [`LibxcFunc.MGGA_C_M11_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_M11_L)


        * [`LibxcFunc.MGGA_C_MN12_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_MN12_L)


        * [`LibxcFunc.MGGA_C_MN12_SX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_MN12_SX)


        * [`LibxcFunc.MGGA_C_MN15`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_MN15)


        * [`LibxcFunc.MGGA_C_MN15_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_MN15_L)


        * [`LibxcFunc.MGGA_C_PKZB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_PKZB)


        * [`LibxcFunc.MGGA_C_REVTPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_REVTPSS)


        * [`LibxcFunc.MGGA_C_SCAN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_SCAN)


        * [`LibxcFunc.MGGA_C_TPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_TPSS)


        * [`LibxcFunc.MGGA_C_TPSSLOC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_TPSSLOC)


        * [`LibxcFunc.MGGA_C_VSXC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_C_VSXC)


        * [`LibxcFunc.MGGA_XC_B97M_V`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_XC_B97M_V)


        * [`LibxcFunc.MGGA_XC_OTPSS_D`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_XC_OTPSS_D)


        * [`LibxcFunc.MGGA_XC_TPSSLYP1W`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_XC_TPSSLYP1W)


        * [`LibxcFunc.MGGA_XC_ZLP`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_XC_ZLP)


        * [`LibxcFunc.MGGA_X_2D_PRHG07`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_2D_PRHG07)


        * [`LibxcFunc.MGGA_X_2D_PRHG07_PRP10`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_2D_PRHG07_PRP10)


        * [`LibxcFunc.MGGA_X_BJ06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_BJ06)


        * [`LibxcFunc.MGGA_X_BLOC`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_BLOC)


        * [`LibxcFunc.MGGA_X_BR89`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_BR89)


        * [`LibxcFunc.MGGA_X_GVT4`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_GVT4)


        * [`LibxcFunc.MGGA_X_LTA`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_LTA)


        * [`LibxcFunc.MGGA_X_M05`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M05)


        * [`LibxcFunc.MGGA_X_M05_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M05_2X)


        * [`LibxcFunc.MGGA_X_M06`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M06)


        * [`LibxcFunc.MGGA_X_M06_2X`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M06_2X)


        * [`LibxcFunc.MGGA_X_M06_HF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M06_HF)


        * [`LibxcFunc.MGGA_X_M06_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M06_L)


        * [`LibxcFunc.MGGA_X_M08_HX`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M08_HX)


        * [`LibxcFunc.MGGA_X_M08_SO`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M08_SO)


        * [`LibxcFunc.MGGA_X_M11`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M11)


        * [`LibxcFunc.MGGA_X_M11_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_M11_L)


        * [`LibxcFunc.MGGA_X_MBEEF`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MBEEF)


        * [`LibxcFunc.MGGA_X_MBEEFVDW`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MBEEFVDW)


        * [`LibxcFunc.MGGA_X_MK00`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MK00)


        * [`LibxcFunc.MGGA_X_MK00B`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MK00B)


        * [`LibxcFunc.MGGA_X_MN12_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MN12_L)


        * [`LibxcFunc.MGGA_X_MN15_L`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MN15_L)


        * [`LibxcFunc.MGGA_X_MODTPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MODTPSS)


        * [`LibxcFunc.MGGA_X_MS0`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MS0)


        * [`LibxcFunc.MGGA_X_MS1`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MS1)


        * [`LibxcFunc.MGGA_X_MS2`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MS2)


        * [`LibxcFunc.MGGA_X_MVS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_MVS)


        * [`LibxcFunc.MGGA_X_PKZB`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_PKZB)


        * [`LibxcFunc.MGGA_X_REVTPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_REVTPSS)


        * [`LibxcFunc.MGGA_X_RPP09`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_RPP09)


        * [`LibxcFunc.MGGA_X_SCAN`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_SCAN)


        * [`LibxcFunc.MGGA_X_TAU_HCTH`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_TAU_HCTH)


        * [`LibxcFunc.MGGA_X_TB09`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_TB09)


        * [`LibxcFunc.MGGA_X_TPSS`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.MGGA_X_TPSS)


        * [`LibxcFunc.all_families()`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.all_families)


        * [`LibxcFunc.all_kinds()`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.all_kinds)


        * [`LibxcFunc.as_dict()`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.as_dict)


        * [`LibxcFunc.from_dict()`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.from_dict)


        * [`LibxcFunc.info_dict`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.info_dict)


        * [`LibxcFunc.is_c_kind`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_c_kind)


        * [`LibxcFunc.is_gga_family`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_gga_family)


        * [`LibxcFunc.is_hyb_gga_family`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_hyb_gga_family)


        * [`LibxcFunc.is_hyb_mgga_family`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_hyb_mgga_family)


        * [`LibxcFunc.is_k_kind`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_k_kind)


        * [`LibxcFunc.is_lda_family`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_lda_family)


        * [`LibxcFunc.is_mgga_family`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_mgga_family)


        * [`LibxcFunc.is_x_kind`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_x_kind)


        * [`LibxcFunc.is_xc_kind`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.is_xc_kind)


        * [`LibxcFunc.to_json()`](pymatgen.core.libxcfunc.md#pymatgen.core.libxcfunc.LibxcFunc.to_json)


* [pymatgen.core.molecular_orbitals module](pymatgen.core.molecular_orbitals.md)


    * [`MolecularOrbitals`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals)


        * [`MolecularOrbitals.composition`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.composition)


        * [`MolecularOrbitals.elements`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.elements)


        * [`MolecularOrbitals.elec_neg`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.elec_neg)


        * [`MolecularOrbitals.aos`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.aos)


        * [`MolecularOrbitals.band_edges`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.band_edges)


        * [`MolecularOrbitals.aos_as_list()`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.aos_as_list)


        * [`MolecularOrbitals.max_electronegativity()`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.max_electronegativity)


        * [`MolecularOrbitals.obtain_band_edges()`](pymatgen.core.molecular_orbitals.md#pymatgen.core.molecular_orbitals.MolecularOrbitals.obtain_band_edges)


* [pymatgen.core.operations module](pymatgen.core.operations.md)


    * [`MagSymmOp`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp)


        * [`MagSymmOp.as_dict()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.as_dict)


        * [`MagSymmOp.as_xyzt_string()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.as_xyzt_string)


        * [`MagSymmOp.from_dict()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.from_dict)


        * [`MagSymmOp.from_rotation_and_translation_and_time_reversal()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.from_rotation_and_translation_and_time_reversal)


        * [`MagSymmOp.from_symmop()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.from_symmop)


        * [`MagSymmOp.from_xyzt_string()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.from_xyzt_string)


        * [`MagSymmOp.operate_magmom()`](pymatgen.core.operations.md#pymatgen.core.operations.MagSymmOp.operate_magmom)


    * [`SymmOp`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp)


        * [`SymmOp.affine_matrix`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.affine_matrix)


        * [`SymmOp.apply_rotation_only()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.apply_rotation_only)


        * [`SymmOp.are_symmetrically_related()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.are_symmetrically_related)


        * [`SymmOp.are_symmetrically_related_vectors()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.are_symmetrically_related_vectors)


        * [`SymmOp.as_dict()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.as_dict)


        * [`SymmOp.as_xyz_string()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.as_xyz_string)


        * [`SymmOp.from_axis_angle_and_translation()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.from_axis_angle_and_translation)


        * [`SymmOp.from_dict()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.from_dict)


        * [`SymmOp.from_origin_axis_angle()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.from_origin_axis_angle)


        * [`SymmOp.from_rotation_and_translation()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.from_rotation_and_translation)


        * [`SymmOp.from_xyz_string()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.from_xyz_string)


        * [`SymmOp.inverse`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.inverse)


        * [`SymmOp.inversion()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.inversion)


        * [`SymmOp.operate()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.operate)


        * [`SymmOp.operate_multi()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.operate_multi)


        * [`SymmOp.reflection()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.reflection)


        * [`SymmOp.rotation_matrix`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.rotation_matrix)


        * [`SymmOp.rotoreflection()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.rotoreflection)


        * [`SymmOp.transform_tensor()`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.transform_tensor)


        * [`SymmOp.translation_vector`](pymatgen.core.operations.md#pymatgen.core.operations.SymmOp.translation_vector)


* [pymatgen.core.periodic_table module](pymatgen.core.periodic_table.md)


    * [`DummySpecie`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecie)


    * [`DummySpecies`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies)


        * [`DummySpecies.oxi_state`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.oxi_state)


        * [`DummySpecies.Z`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.Z)


        * [`DummySpecies.X`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.X)


        * [`DummySpecies.X`](pymatgen.core.periodic_table.md#id0)


        * [`DummySpecies.Z`](pymatgen.core.periodic_table.md#id1)


        * [`DummySpecies.as_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.as_dict)


        * [`DummySpecies.from_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.from_dict)


        * [`DummySpecies.from_str()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.from_str)


        * [`DummySpecies.oxi_state`](pymatgen.core.periodic_table.md#id2)


        * [`DummySpecies.symbol`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.DummySpecies.symbol)


    * [`Element`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element)


        * [`Element.Z`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Z)


        * [`Element.symbol`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.symbol)


        * [`Element.long_name`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.long_name)


        * [`Element.atomic_radius_calculated`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.atomic_radius_calculated)


        * [`Element.van_der_waals_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.van_der_waals_radius)


        * [`Element.mendeleev_no`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.mendeleev_no)


        * [`Element.electrical_resistivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.electrical_resistivity)


        * [`Element.velocity_of_sound`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.velocity_of_sound)


        * [`Element.reflectivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.reflectivity)


        * [`Element.refractive_index`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.refractive_index)


        * [`Element.poissons_ratio`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.poissons_ratio)


        * [`Element.molar_volume`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.molar_volume)


        * [`Element.electronic_structure`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.electronic_structure)


        * [`Element.atomic_orbitals`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.atomic_orbitals)


        * [`Element.thermal_conductivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.thermal_conductivity)


        * [`Element.boiling_point`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.boiling_point)


        * [`Element.melting_point`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.melting_point)


        * [`Element.critical_temperature`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.critical_temperature)


        * [`Element.superconduction_temperature`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.superconduction_temperature)


        * [`Element.liquid_range`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.liquid_range)


        * [`Element.bulk_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.bulk_modulus)


        * [`Element.youngs_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.youngs_modulus)


        * [`Element.brinell_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.brinell_hardness)


        * [`Element.rigidity_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.rigidity_modulus)


        * [`Element.mineral_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.mineral_hardness)


        * [`Element.vickers_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.vickers_hardness)


        * [`Element.density_of_solid`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.density_of_solid)


        * [`Element.coefficient_of_linear_thermal_expansion`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.coefficient_of_linear_thermal_expansion)


        * [`Element.ground_level`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.ground_level)


        * [`Element.ionization_energies`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.ionization_energies)


        * [`Element.Ac`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ac)


        * [`Element.Ag`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ag)


        * [`Element.Al`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Al)


        * [`Element.Am`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Am)


        * [`Element.Ar`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ar)


        * [`Element.As`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.As)


        * [`Element.At`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.At)


        * [`Element.Au`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Au)


        * [`Element.B`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.B)


        * [`Element.Ba`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ba)


        * [`Element.Be`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Be)


        * [`Element.Bh`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Bh)


        * [`Element.Bi`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Bi)


        * [`Element.Bk`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Bk)


        * [`Element.Br`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Br)


        * [`Element.C`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.C)


        * [`Element.Ca`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ca)


        * [`Element.Cd`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cd)


        * [`Element.Ce`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ce)


        * [`Element.Cf`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cf)


        * [`Element.Cl`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cl)


        * [`Element.Cm`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cm)


        * [`Element.Cn`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cn)


        * [`Element.Co`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Co)


        * [`Element.Cr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cr)


        * [`Element.Cs`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cs)


        * [`Element.Cu`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Cu)


        * [`Element.Db`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Db)


        * [`Element.Ds`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ds)


        * [`Element.Dy`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Dy)


        * [`Element.Er`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Er)


        * [`Element.Es`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Es)


        * [`Element.Eu`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Eu)


        * [`Element.F`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.F)


        * [`Element.Fe`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Fe)


        * [`Element.Fl`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Fl)


        * [`Element.Fm`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Fm)


        * [`Element.Fr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Fr)


        * [`Element.Ga`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ga)


        * [`Element.Gd`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Gd)


        * [`Element.Ge`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ge)


        * [`Element.H`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.H)


        * [`Element.He`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.He)


        * [`Element.Hf`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Hf)


        * [`Element.Hg`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Hg)


        * [`Element.Ho`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ho)


        * [`Element.Hs`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Hs)


        * [`Element.I`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.I)


        * [`Element.In`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.In)


        * [`Element.Ir`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ir)


        * [`Element.K`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.K)


        * [`Element.Kr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Kr)


        * [`Element.La`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.La)


        * [`Element.Li`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Li)


        * [`Element.Lr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Lr)


        * [`Element.Lu`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Lu)


        * [`Element.Lv`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Lv)


        * [`Element.Mc`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Mc)


        * [`Element.Md`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Md)


        * [`Element.Mg`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Mg)


        * [`Element.Mn`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Mn)


        * [`Element.Mo`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Mo)


        * [`Element.Mt`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Mt)


        * [`Element.N`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.N)


        * [`Element.Na`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Na)


        * [`Element.Nb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Nb)


        * [`Element.Nd`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Nd)


        * [`Element.Ne`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ne)


        * [`Element.Nh`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Nh)


        * [`Element.Ni`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ni)


        * [`Element.No`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.No)


        * [`Element.Np`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Np)


        * [`Element.O`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.O)


        * [`Element.Og`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Og)


        * [`Element.Os`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Os)


        * [`Element.P`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.P)


        * [`Element.Pa`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pa)


        * [`Element.Pb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pb)


        * [`Element.Pd`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pd)


        * [`Element.Pm`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pm)


        * [`Element.Po`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Po)


        * [`Element.Pr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pr)


        * [`Element.Pt`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pt)


        * [`Element.Pu`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Pu)


        * [`Element.Ra`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ra)


        * [`Element.Rb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Rb)


        * [`Element.Re`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Re)


        * [`Element.Rf`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Rf)


        * [`Element.Rg`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Rg)


        * [`Element.Rh`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Rh)


        * [`Element.Rn`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Rn)


        * [`Element.Ru`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ru)


        * [`Element.S`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.S)


        * [`Element.Sb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sb)


        * [`Element.Sc`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sc)


        * [`Element.Se`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Se)


        * [`Element.Sg`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sg)


        * [`Element.Si`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Si)


        * [`Element.Sm`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sm)


        * [`Element.Sn`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sn)


        * [`Element.Sr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Sr)


        * [`Element.Ta`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ta)


        * [`Element.Tb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Tb)


        * [`Element.Tc`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Tc)


        * [`Element.Te`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Te)


        * [`Element.Th`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Th)


        * [`Element.Ti`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ti)


        * [`Element.Tl`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Tl)


        * [`Element.Tm`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Tm)


        * [`Element.Ts`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Ts)


        * [`Element.U`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.U)


        * [`Element.V`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.V)


        * [`Element.W`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.W)


        * [`Element.Xe`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Xe)


        * [`Element.Y`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Y)


        * [`Element.Yb`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Yb)


        * [`Element.Zn`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Zn)


        * [`Element.Zr`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Element.Zr)


    * [`ElementBase`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase)


        * [`ElementBase.Z`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.Z)


        * [`ElementBase.symbol`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.symbol)


        * [`ElementBase.long_name`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.long_name)


        * [`ElementBase.atomic_radius_calculated`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.atomic_radius_calculated)


        * [`ElementBase.van_der_waals_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.van_der_waals_radius)


        * [`ElementBase.mendeleev_no`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.mendeleev_no)


        * [`ElementBase.electrical_resistivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.electrical_resistivity)


        * [`ElementBase.velocity_of_sound`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.velocity_of_sound)


        * [`ElementBase.reflectivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.reflectivity)


        * [`ElementBase.refractive_index`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.refractive_index)


        * [`ElementBase.poissons_ratio`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.poissons_ratio)


        * [`ElementBase.molar_volume`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.molar_volume)


        * [`ElementBase.electronic_structure`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.electronic_structure)


        * [`ElementBase.atomic_orbitals`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.atomic_orbitals)


        * [`ElementBase.thermal_conductivity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.thermal_conductivity)


        * [`ElementBase.boiling_point`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.boiling_point)


        * [`ElementBase.melting_point`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.melting_point)


        * [`ElementBase.critical_temperature`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.critical_temperature)


        * [`ElementBase.superconduction_temperature`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.superconduction_temperature)


        * [`ElementBase.liquid_range`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.liquid_range)


        * [`ElementBase.bulk_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.bulk_modulus)


        * [`ElementBase.youngs_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.youngs_modulus)


        * [`ElementBase.brinell_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.brinell_hardness)


        * [`ElementBase.rigidity_modulus`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.rigidity_modulus)


        * [`ElementBase.mineral_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.mineral_hardness)


        * [`ElementBase.vickers_hardness`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.vickers_hardness)


        * [`ElementBase.density_of_solid`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.density_of_solid)


        * [`ElementBase.coefficient_of_linear_thermal_expansion`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.coefficient_of_linear_thermal_expansion)


        * [`ElementBase.ground_level`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.ground_level)


        * [`ElementBase.ionization_energies`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.ionization_energies)


        * [`ElementBase.X`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.X)


        * [`ElementBase.as_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.as_dict)


        * [`ElementBase.atomic_mass`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.atomic_mass)


        * [`ElementBase.atomic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.atomic_radius)


        * [`ElementBase.average_anionic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.average_anionic_radius)


        * [`ElementBase.average_cationic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.average_cationic_radius)


        * [`ElementBase.average_ionic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.average_ionic_radius)


        * [`ElementBase.block`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.block)


        * [`ElementBase.common_oxidation_states`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.common_oxidation_states)


        * [`ElementBase.data`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.data)


        * [`ElementBase.electron_affinity`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.electron_affinity)


        * [`ElementBase.electronic_structure`](pymatgen.core.periodic_table.md#id3)


        * [`ElementBase.from_Z()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.from_Z)


        * [`ElementBase.from_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.from_dict)


        * [`ElementBase.from_name()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.from_name)


        * [`ElementBase.from_row_and_group()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.from_row_and_group)


        * [`ElementBase.full_electronic_structure`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.full_electronic_structure)


        * [`ElementBase.ground_state_term_symbol`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.ground_state_term_symbol)


        * [`ElementBase.group`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.group)


        * [`ElementBase.icsd_oxidation_states`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.icsd_oxidation_states)


        * [`ElementBase.ionic_radii`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.ionic_radii)


        * [`ElementBase.ionization_energy`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.ionization_energy)


        * [`ElementBase.is_actinoid`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_actinoid)


        * [`ElementBase.is_alkali`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_alkali)


        * [`ElementBase.is_alkaline`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_alkaline)


        * [`ElementBase.is_chalcogen`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_chalcogen)


        * [`ElementBase.is_halogen`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_halogen)


        * [`ElementBase.is_lanthanoid`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_lanthanoid)


        * [`ElementBase.is_metal`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_metal)


        * [`ElementBase.is_metalloid`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_metalloid)


        * [`ElementBase.is_noble_gas`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_noble_gas)


        * [`ElementBase.is_post_transition_metal`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_post_transition_metal)


        * [`ElementBase.is_quadrupolar`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_quadrupolar)


        * [`ElementBase.is_rare_earth_metal`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_rare_earth_metal)


        * [`ElementBase.is_transition_metal`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_transition_metal)


        * [`ElementBase.is_valid_symbol()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.is_valid_symbol)


        * [`ElementBase.iupac_ordering`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.iupac_ordering)


        * [`ElementBase.max_oxidation_state`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.max_oxidation_state)


        * [`ElementBase.metallic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.metallic_radius)


        * [`ElementBase.min_oxidation_state`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.min_oxidation_state)


        * [`ElementBase.nmr_quadrupole_moment`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.nmr_quadrupole_moment)


        * [`ElementBase.number`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.number)


        * [`ElementBase.oxidation_states`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.oxidation_states)


        * [`ElementBase.print_periodic_table()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.print_periodic_table)


        * [`ElementBase.row`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.row)


        * [`ElementBase.term_symbols`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.term_symbols)


        * [`ElementBase.valence`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.ElementBase.valence)


    * [`Specie`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Specie)


    * [`Species`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species)


        * [`Species.STRING_MODE`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.STRING_MODE)


        * [`Species.as_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.as_dict)


        * [`Species.element`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.element)


        * [`Species.from_dict()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.from_dict)


        * [`Species.from_str()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.from_str)


        * [`Species.from_string()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.from_string)


        * [`Species.get_crystal_field_spin()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.get_crystal_field_spin)


        * [`Species.get_nmr_quadrupole_moment()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.get_nmr_quadrupole_moment)


        * [`Species.get_shannon_radius()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.get_shannon_radius)


        * [`Species.ionic_radius`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.ionic_radius)


        * [`Species.oxi_state`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.oxi_state)


        * [`Species.properties`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.properties)


        * [`Species.spin`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.spin)


        * [`Species.to_pretty_string()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species.to_pretty_string)


    * [`get_el_sp()`](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.get_el_sp)


* [pymatgen.core.sites module](pymatgen.core.sites.md)


    * [`PeriodicSite`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite)


        * [`PeriodicSite.a`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.a)


        * [`PeriodicSite.as_dict()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.as_dict)


        * [`PeriodicSite.b`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.b)


        * [`PeriodicSite.c`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.c)


        * [`PeriodicSite.coords`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.coords)


        * [`PeriodicSite.distance()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.distance)


        * [`PeriodicSite.distance_and_image()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.distance_and_image)


        * [`PeriodicSite.distance_and_image_from_frac_coords()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.distance_and_image_from_frac_coords)


        * [`PeriodicSite.frac_coords`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.frac_coords)


        * [`PeriodicSite.from_dict()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.from_dict)


        * [`PeriodicSite.is_periodic_image()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.is_periodic_image)


        * [`PeriodicSite.lattice`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.lattice)


        * [`PeriodicSite.to_unit_cell()`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.to_unit_cell)


        * [`PeriodicSite.x`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.x)


        * [`PeriodicSite.y`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.y)


        * [`PeriodicSite.z`](pymatgen.core.sites.md#pymatgen.core.sites.PeriodicSite.z)


    * [`Site`](pymatgen.core.sites.md#pymatgen.core.sites.Site)


        * [`Site.as_dict()`](pymatgen.core.sites.md#pymatgen.core.sites.Site.as_dict)


        * [`Site.distance()`](pymatgen.core.sites.md#pymatgen.core.sites.Site.distance)


        * [`Site.distance_from_point()`](pymatgen.core.sites.md#pymatgen.core.sites.Site.distance_from_point)


        * [`Site.from_dict()`](pymatgen.core.sites.md#pymatgen.core.sites.Site.from_dict)


        * [`Site.is_ordered`](pymatgen.core.sites.md#pymatgen.core.sites.Site.is_ordered)


        * [`Site.label`](pymatgen.core.sites.md#pymatgen.core.sites.Site.label)


        * [`Site.position_atol`](pymatgen.core.sites.md#pymatgen.core.sites.Site.position_atol)


        * [`Site.specie`](pymatgen.core.sites.md#pymatgen.core.sites.Site.specie)


        * [`Site.species`](pymatgen.core.sites.md#pymatgen.core.sites.Site.species)


        * [`Site.species_string`](pymatgen.core.sites.md#pymatgen.core.sites.Site.species_string)


        * [`Site.x`](pymatgen.core.sites.md#pymatgen.core.sites.Site.x)


        * [`Site.y`](pymatgen.core.sites.md#pymatgen.core.sites.Site.y)


        * [`Site.z`](pymatgen.core.sites.md#pymatgen.core.sites.Site.z)


* [pymatgen.core.spectrum module](pymatgen.core.spectrum.md)


    * [`Spectrum`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum)


        * [`Spectrum.XLABEL`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.XLABEL)


        * [`Spectrum.YLABEL`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.YLABEL)


        * [`Spectrum.copy()`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.copy)


        * [`Spectrum.get_interpolated_value()`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.get_interpolated_value)


        * [`Spectrum.normalize()`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.normalize)


        * [`Spectrum.smear()`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.Spectrum.smear)


    * [`lorentzian()`](pymatgen.core.spectrum.md#pymatgen.core.spectrum.lorentzian)


* [pymatgen.core.structure module](pymatgen.core.structure.md)


    * [`IMolecule`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule)


        * [`IMolecule.as_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.as_dict)


        * [`IMolecule.break_bond()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.break_bond)


        * [`IMolecule.center_of_mass`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.center_of_mass)


        * [`IMolecule.charge`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.charge)


        * [`IMolecule.from_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.from_dict)


        * [`IMolecule.from_file()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.from_file)


        * [`IMolecule.from_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.from_sites)


        * [`IMolecule.from_str()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.from_str)


        * [`IMolecule.get_boxed_structure()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_boxed_structure)


        * [`IMolecule.get_centered_molecule()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_centered_molecule)


        * [`IMolecule.get_covalent_bonds()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_covalent_bonds)


        * [`IMolecule.get_distance()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_distance)


        * [`IMolecule.get_neighbors()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_neighbors)


        * [`IMolecule.get_neighbors_in_shell()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_neighbors_in_shell)


        * [`IMolecule.get_sites_in_sphere()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_sites_in_sphere)


        * [`IMolecule.get_zmatrix()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.get_zmatrix)


        * [`IMolecule.nelectrons`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.nelectrons)


        * [`IMolecule.spin_multiplicity`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.spin_multiplicity)


        * [`IMolecule.to()`](pymatgen.core.structure.md#pymatgen.core.structure.IMolecule.to)


    * [`IStructure`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure)


        * [`IStructure.as_dataframe()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.as_dataframe)


        * [`IStructure.as_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.as_dict)


        * [`IStructure.charge`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.charge)


        * [`IStructure.copy()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.copy)


        * [`IStructure.density`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.density)


        * [`IStructure.distance_matrix`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.distance_matrix)


        * [`IStructure.frac_coords`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.frac_coords)


        * [`IStructure.from_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_dict)


        * [`IStructure.from_file()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_file)


        * [`IStructure.from_magnetic_spacegroup()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_magnetic_spacegroup)


        * [`IStructure.from_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_sites)


        * [`IStructure.from_spacegroup()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_spacegroup)


        * [`IStructure.from_str()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.from_str)


        * [`IStructure.get_all_neighbors()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_all_neighbors)


        * [`IStructure.get_all_neighbors_old()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_all_neighbors_old)


        * [`IStructure.get_all_neighbors_py()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_all_neighbors_py)


        * [`IStructure.get_distance()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_distance)


        * [`IStructure.get_miller_index_from_site_indexes()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_miller_index_from_site_indexes)


        * [`IStructure.get_neighbor_list()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_neighbor_list)


        * [`IStructure.get_neighbors()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_neighbors)


        * [`IStructure.get_neighbors_in_shell()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_neighbors_in_shell)


        * [`IStructure.get_neighbors_old()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_neighbors_old)


        * [`IStructure.get_orderings()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_orderings)


        * [`IStructure.get_primitive_structure()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_primitive_structure)


        * [`IStructure.get_reduced_structure()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_reduced_structure)


        * [`IStructure.get_sites_in_sphere()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_sites_in_sphere)


        * [`IStructure.get_sorted_structure()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_sorted_structure)


        * [`IStructure.get_space_group_info()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_space_group_info)


        * [`IStructure.get_symmetric_neighbor_list()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.get_symmetric_neighbor_list)


        * [`IStructure.interpolate()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.interpolate)


        * [`IStructure.is_3d_periodic`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.is_3d_periodic)


        * [`IStructure.lattice`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.lattice)


        * [`IStructure.matches()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.matches)


        * [`IStructure.pbc`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.pbc)


        * [`IStructure.to()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.to)


        * [`IStructure.unset_charge()`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.unset_charge)


        * [`IStructure.volume`](pymatgen.core.structure.md#pymatgen.core.structure.IStructure.volume)


    * [`Molecule`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)


        * [`Molecule.append()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.append)


        * [`Molecule.apply_operation()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.apply_operation)


        * [`Molecule.calculate()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.calculate)


        * [`Molecule.copy()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.copy)


        * [`Molecule.insert()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.insert)


        * [`Molecule.perturb()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.perturb)


        * [`Molecule.relax()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.relax)


        * [`Molecule.remove_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.remove_sites)


        * [`Molecule.remove_species()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.remove_species)


        * [`Molecule.rotate_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.rotate_sites)


        * [`Molecule.set_charge_and_spin()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.set_charge_and_spin)


        * [`Molecule.substitute()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.substitute)


        * [`Molecule.translate_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Molecule.translate_sites)


    * [`Neighbor`](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor)


        * [`Neighbor.as_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor.as_dict)


        * [`Neighbor.coords`](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor.coords)


        * [`Neighbor.from_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor.from_dict)


        * [`Neighbor.properties`](pymatgen.core.structure.md#pymatgen.core.structure.Neighbor.properties)


    * [`PeriodicNeighbor`](pymatgen.core.structure.md#pymatgen.core.structure.PeriodicNeighbor)


        * [`PeriodicNeighbor.as_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.PeriodicNeighbor.as_dict)


        * [`PeriodicNeighbor.coords`](pymatgen.core.structure.md#pymatgen.core.structure.PeriodicNeighbor.coords)


        * [`PeriodicNeighbor.from_dict()`](pymatgen.core.structure.md#pymatgen.core.structure.PeriodicNeighbor.from_dict)


        * [`PeriodicNeighbor.properties`](pymatgen.core.structure.md#pymatgen.core.structure.PeriodicNeighbor.properties)


    * [`SiteCollection`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection)


        * [`SiteCollection.DISTANCE_TOLERANCE`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.DISTANCE_TOLERANCE)


        * [`SiteCollection.add_oxidation_state_by_element()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_oxidation_state_by_element)


        * [`SiteCollection.add_oxidation_state_by_guess()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_oxidation_state_by_guess)


        * [`SiteCollection.add_oxidation_state_by_site()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_oxidation_state_by_site)


        * [`SiteCollection.add_site_property()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_site_property)


        * [`SiteCollection.add_spin_by_element()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_spin_by_element)


        * [`SiteCollection.add_spin_by_site()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.add_spin_by_site)


        * [`SiteCollection.atomic_numbers`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.atomic_numbers)


        * [`SiteCollection.cart_coords`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.cart_coords)


        * [`SiteCollection.charge`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.charge)


        * [`SiteCollection.composition`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.composition)


        * [`SiteCollection.distance_matrix`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.distance_matrix)


        * [`SiteCollection.extract_cluster()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.extract_cluster)


        * [`SiteCollection.formula`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.formula)


        * [`SiteCollection.from_file()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.from_file)


        * [`SiteCollection.from_str()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.from_str)


        * [`SiteCollection.get_angle()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.get_angle)


        * [`SiteCollection.get_dihedral()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.get_dihedral)


        * [`SiteCollection.get_distance()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.get_distance)


        * [`SiteCollection.group_by_types()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.group_by_types)


        * [`SiteCollection.indices_from_symbol()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.indices_from_symbol)


        * [`SiteCollection.is_ordered`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.is_ordered)


        * [`SiteCollection.is_valid()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.is_valid)


        * [`SiteCollection.labels`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.labels)


        * [`SiteCollection.ntypesp`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.ntypesp)


        * [`SiteCollection.num_sites`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.num_sites)


        * [`SiteCollection.remove_oxidation_states()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.remove_oxidation_states)


        * [`SiteCollection.remove_site_property()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.remove_site_property)


        * [`SiteCollection.remove_spin()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.remove_spin)


        * [`SiteCollection.replace_species()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.replace_species)


        * [`SiteCollection.site_properties`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.site_properties)


        * [`SiteCollection.sites`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.sites)


        * [`SiteCollection.species`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.species)


        * [`SiteCollection.species_and_occu`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.species_and_occu)


        * [`SiteCollection.symbol_set`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.symbol_set)


        * [`SiteCollection.to()`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.to)


        * [`SiteCollection.types_of_specie`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.types_of_specie)


        * [`SiteCollection.types_of_species`](pymatgen.core.structure.md#pymatgen.core.structure.SiteCollection.types_of_species)


    * [`Structure`](pymatgen.core.structure.md#pymatgen.core.structure.Structure)


        * [`Structure.append()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.append)


        * [`Structure.apply_operation()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.apply_operation)


        * [`Structure.apply_strain()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.apply_strain)


        * [`Structure.calculate()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.calculate)


        * [`Structure.from_prototype()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.from_prototype)


        * [`Structure.insert()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.insert)


        * [`Structure.lattice`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.lattice)


        * [`Structure.make_supercell()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.make_supercell)


        * [`Structure.merge_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.merge_sites)


        * [`Structure.perturb()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.perturb)


        * [`Structure.relax()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.relax)


        * [`Structure.remove_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.remove_sites)


        * [`Structure.remove_species()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.remove_species)


        * [`Structure.replace()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.replace)


        * [`Structure.rotate_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.rotate_sites)


        * [`Structure.scale_lattice()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.scale_lattice)


        * [`Structure.set_charge()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.set_charge)


        * [`Structure.sort()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.sort)


        * [`Structure.substitute()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.substitute)


        * [`Structure.translate_sites()`](pymatgen.core.structure.md#pymatgen.core.structure.Structure.translate_sites)


    * [`StructureError`](pymatgen.core.structure.md#pymatgen.core.structure.StructureError)


* [pymatgen.core.surface module](pymatgen.core.surface.md)


    * [`ReconstructionGenerator`](pymatgen.core.surface.md#pymatgen.core.surface.ReconstructionGenerator)


        * [`ReconstructionGenerator.slabgen_params`](pymatgen.core.surface.md#pymatgen.core.surface.ReconstructionGenerator.slabgen_params)


        * [`ReconstructionGenerator.build_slabs()`](pymatgen.core.surface.md#pymatgen.core.surface.ReconstructionGenerator.build_slabs)


        * [`ReconstructionGenerator.get_unreconstructed_slabs()`](pymatgen.core.surface.md#pymatgen.core.surface.ReconstructionGenerator.get_unreconstructed_slabs)


    * [`Slab`](pymatgen.core.surface.md#pymatgen.core.surface.Slab)


        * [`Slab.miller_index`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.miller_index)


        * [`Slab.scale_factor`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.scale_factor)


        * [`Slab.shift`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.shift)


        * [`Slab.add_adsorbate_atom()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.add_adsorbate_atom)


        * [`Slab.as_dict()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.as_dict)


        * [`Slab.center_of_mass`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.center_of_mass)


        * [`Slab.copy()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.copy)


        * [`Slab.dipole`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.dipole)


        * [`Slab.from_dict()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.from_dict)


        * [`Slab.get_orthogonal_c_slab()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.get_orthogonal_c_slab)


        * [`Slab.get_sorted_structure()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.get_sorted_structure)


        * [`Slab.get_surface_sites()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.get_surface_sites)


        * [`Slab.get_symmetric_site()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.get_symmetric_site)


        * [`Slab.get_tasker2_slabs()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.get_tasker2_slabs)


        * [`Slab.is_polar()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.is_polar)


        * [`Slab.is_symmetric()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.is_symmetric)


        * [`Slab.normal`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.normal)


        * [`Slab.surface_area`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.surface_area)


        * [`Slab.symmetrically_add_atom()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.symmetrically_add_atom)


        * [`Slab.symmetrically_remove_atoms()`](pymatgen.core.surface.md#pymatgen.core.surface.Slab.symmetrically_remove_atoms)


    * [`SlabGenerator`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator)


        * [`SlabGenerator.oriented_unit_cell`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.oriented_unit_cell)


        * [`SlabGenerator.parent`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.parent)


        * [`SlabGenerator.lll_reduce`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.lll_reduce)


        * [`SlabGenerator.center_slab`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.center_slab)


        * [`SlabGenerator.slab_scale_factor`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.slab_scale_factor)


        * [`SlabGenerator.miller_index`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.miller_index)


        * [`SlabGenerator.min_slab_size`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.min_slab_size)


        * [`SlabGenerator.min_vac_size`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.min_vac_size)


        * [`SlabGenerator.get_slab()`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.get_slab)


        * [`SlabGenerator.get_slabs()`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.get_slabs)


        * [`SlabGenerator.move_to_other_side()`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.move_to_other_side)


        * [`SlabGenerator.nonstoichiometric_symmetrized_slab()`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.nonstoichiometric_symmetrized_slab)


        * [`SlabGenerator.repair_broken_bonds()`](pymatgen.core.surface.md#pymatgen.core.surface.SlabGenerator.repair_broken_bonds)


    * [`center_slab()`](pymatgen.core.surface.md#pymatgen.core.surface.center_slab)


    * [`generate_all_slabs()`](pymatgen.core.surface.md#pymatgen.core.surface.generate_all_slabs)


    * [`get_d()`](pymatgen.core.surface.md#pymatgen.core.surface.get_d)


    * [`get_slab_regions()`](pymatgen.core.surface.md#pymatgen.core.surface.get_slab_regions)


    * [`get_symmetrically_distinct_miller_indices()`](pymatgen.core.surface.md#pymatgen.core.surface.get_symmetrically_distinct_miller_indices)


    * [`get_symmetrically_equivalent_miller_indices()`](pymatgen.core.surface.md#pymatgen.core.surface.get_symmetrically_equivalent_miller_indices)


    * [`hkl_transformation()`](pymatgen.core.surface.md#pymatgen.core.surface.hkl_transformation)


    * [`is_already_analyzed()`](pymatgen.core.surface.md#pymatgen.core.surface.is_already_analyzed)


    * [`miller_index_from_sites()`](pymatgen.core.surface.md#pymatgen.core.surface.miller_index_from_sites)


* [pymatgen.core.tensors module](pymatgen.core.tensors.md)


    * [`SquareTensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor)


        * [`SquareTensor.det`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.det)


        * [`SquareTensor.get_scaled()`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.get_scaled)


        * [`SquareTensor.inv`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.inv)


        * [`SquareTensor.is_rotation()`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.is_rotation)


        * [`SquareTensor.polar_decomposition()`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.polar_decomposition)


        * [`SquareTensor.principal_invariants`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.principal_invariants)


        * [`SquareTensor.refine_rotation()`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.refine_rotation)


        * [`SquareTensor.trans`](pymatgen.core.tensors.md#pymatgen.core.tensors.SquareTensor.trans)


    * [`Tensor`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor)


        * [`Tensor.as_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.as_dict)


        * [`Tensor.average_over_unit_sphere()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.average_over_unit_sphere)


        * [`Tensor.convert_to_ieee()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.convert_to_ieee)


        * [`Tensor.einsum_sequence()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.einsum_sequence)


        * [`Tensor.fit_to_structure()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.fit_to_structure)


        * [`Tensor.from_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.from_dict)


        * [`Tensor.from_values_indices()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.from_values_indices)


        * [`Tensor.from_voigt()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.from_voigt)


        * [`Tensor.get_grouped_indices()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.get_grouped_indices)


        * [`Tensor.get_ieee_rotation()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.get_ieee_rotation)


        * [`Tensor.get_symbol_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.get_symbol_dict)


        * [`Tensor.get_voigt_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.get_voigt_dict)


        * [`Tensor.is_fit_to_structure()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.is_fit_to_structure)


        * [`Tensor.is_symmetric()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.is_symmetric)


        * [`Tensor.is_voigt_symmetric()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.is_voigt_symmetric)


        * [`Tensor.populate()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.populate)


        * [`Tensor.project()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.project)


        * [`Tensor.rotate()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.rotate)


        * [`Tensor.round()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.round)


        * [`Tensor.structure_transform()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.structure_transform)


        * [`Tensor.symbol`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.symbol)


        * [`Tensor.symmetrized`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.symmetrized)


        * [`Tensor.transform()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.transform)


        * [`Tensor.voigt`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.voigt)


        * [`Tensor.voigt_symmetrized`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.voigt_symmetrized)


        * [`Tensor.zeroed()`](pymatgen.core.tensors.md#pymatgen.core.tensors.Tensor.zeroed)


    * [`TensorCollection`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection)


        * [`TensorCollection.as_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.as_dict)


        * [`TensorCollection.convert_to_ieee()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.convert_to_ieee)


        * [`TensorCollection.fit_to_structure()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.fit_to_structure)


        * [`TensorCollection.from_dict()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.from_dict)


        * [`TensorCollection.from_voigt()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.from_voigt)


        * [`TensorCollection.is_fit_to_structure()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.is_fit_to_structure)


        * [`TensorCollection.is_symmetric()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.is_symmetric)


        * [`TensorCollection.is_voigt_symmetric()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.is_voigt_symmetric)


        * [`TensorCollection.ranks`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.ranks)


        * [`TensorCollection.rotate()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.rotate)


        * [`TensorCollection.round()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.round)


        * [`TensorCollection.symmetrized`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.symmetrized)


        * [`TensorCollection.transform()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.transform)


        * [`TensorCollection.voigt`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.voigt)


        * [`TensorCollection.voigt_symmetrized`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.voigt_symmetrized)


        * [`TensorCollection.zeroed()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorCollection.zeroed)


    * [`TensorMapping`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorMapping)


        * [`TensorMapping.items()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorMapping.items)


        * [`TensorMapping.values()`](pymatgen.core.tensors.md#pymatgen.core.tensors.TensorMapping.values)


    * [`get_uvec()`](pymatgen.core.tensors.md#pymatgen.core.tensors.get_uvec)


    * [`symmetry_reduce()`](pymatgen.core.tensors.md#pymatgen.core.tensors.symmetry_reduce)


* [pymatgen.core.trajectory module](pymatgen.core.trajectory.md)


    * [`Trajectory`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory)


        * [`Trajectory.as_dict()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.as_dict)


        * [`Trajectory.extend()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.extend)


        * [`Trajectory.from_file()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.from_file)


        * [`Trajectory.from_molecules()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.from_molecules)


        * [`Trajectory.from_structures()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.from_structures)


        * [`Trajectory.get_molecule()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.get_molecule)


        * [`Trajectory.get_structure()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.get_structure)


        * [`Trajectory.to_displacements()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.to_displacements)


        * [`Trajectory.to_positions()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.to_positions)


        * [`Trajectory.write_Xdatcar()`](pymatgen.core.trajectory.md#pymatgen.core.trajectory.Trajectory.write_Xdatcar)


* [pymatgen.core.units module](pymatgen.core.units.md)


    * [`ArrayWithUnit`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit)


        * [`ArrayWithUnit.Error`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.Error)


        * [`ArrayWithUnit.as_base_units`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.as_base_units)


        * [`ArrayWithUnit.conversions()`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.conversions)


        * [`ArrayWithUnit.supported_units`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.supported_units)


        * [`ArrayWithUnit.to()`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.to)


        * [`ArrayWithUnit.unit`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.unit)


        * [`ArrayWithUnit.unit_type`](pymatgen.core.units.md#pymatgen.core.units.ArrayWithUnit.unit_type)


    * [`Charge`](pymatgen.core.units.md#pymatgen.core.units.Charge)


    * [`Energy`](pymatgen.core.units.md#pymatgen.core.units.Energy)


    * [`FloatWithUnit`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit)


        * [`FloatWithUnit.Error`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.Error)


        * [`FloatWithUnit.as_base_units`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.as_base_units)


        * [`FloatWithUnit.from_str()`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.from_str)


        * [`FloatWithUnit.from_string()`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.from_string)


        * [`FloatWithUnit.supported_units`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.supported_units)


        * [`FloatWithUnit.to()`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.to)


        * [`FloatWithUnit.unit`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.unit)


        * [`FloatWithUnit.unit_type`](pymatgen.core.units.md#pymatgen.core.units.FloatWithUnit.unit_type)


    * [`Length`](pymatgen.core.units.md#pymatgen.core.units.Length)


    * [`Mass`](pymatgen.core.units.md#pymatgen.core.units.Mass)


    * [`Memory`](pymatgen.core.units.md#pymatgen.core.units.Memory)


    * [`Temp`](pymatgen.core.units.md#pymatgen.core.units.Temp)


    * [`Time`](pymatgen.core.units.md#pymatgen.core.units.Time)


    * [`Unit`](pymatgen.core.units.md#pymatgen.core.units.Unit)


        * [`Unit.Error`](pymatgen.core.units.md#pymatgen.core.units.Unit.Error)


        * [`Unit.as_base_units`](pymatgen.core.units.md#pymatgen.core.units.Unit.as_base_units)


        * [`Unit.get_conversion_factor()`](pymatgen.core.units.md#pymatgen.core.units.Unit.get_conversion_factor)


    * [`UnitError`](pymatgen.core.units.md#pymatgen.core.units.UnitError)


    * [`kb`](pymatgen.core.units.md#pymatgen.core.units.kb)


    * [`obj_with_unit()`](pymatgen.core.units.md#pymatgen.core.units.obj_with_unit)


    * [`unitized()`](pymatgen.core.units.md#pymatgen.core.units.unitized)


* [pymatgen.core.xcfunc module](pymatgen.core.xcfunc.md)


    * [`XcFunc`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc)


        * [`XcFunc.abinitixc_to_libxc`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.abinitixc_to_libxc)


        * [`XcFunc.aliases()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.aliases)


        * [`XcFunc.as_dict()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.as_dict)


        * [`XcFunc.asxc()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.asxc)


        * [`XcFunc.defined_aliases`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.defined_aliases)


        * [`XcFunc.from_abinit_ixc()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.from_abinit_ixc)


        * [`XcFunc.from_dict()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.from_dict)


        * [`XcFunc.from_name()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.from_name)


        * [`XcFunc.from_type_name()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.from_type_name)


        * [`XcFunc.name()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.name)


        * [`XcFunc.type()`](pymatgen.core.xcfunc.md#pymatgen.core.xcfunc.XcFunc.type)