---
layout: default
title: pymatgen.util.md
nav_exclude: true
---

# pymatgen.util package

The util package implements various utilities that are commonly used by various
packages.



* [pymatgen.util.convergence module](pymatgen.util.convergence.md)


    * [`SplineInputError`](pymatgen.util.convergence.md#pymatgen.util.convergence.SplineInputError)


    * [`determine_convergence()`](pymatgen.util.convergence.md#pymatgen.util.convergence.determine_convergence)


    * [`exponential()`](pymatgen.util.convergence.md#pymatgen.util.convergence.exponential)


    * [`extrapolate_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.extrapolate_reciprocal)


    * [`extrapolate_simple_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.extrapolate_simple_reciprocal)


    * [`get_derivatives()`](pymatgen.util.convergence.md#pymatgen.util.convergence.get_derivatives)


    * [`get_weights()`](pymatgen.util.convergence.md#pymatgen.util.convergence.get_weights)


    * [`id_generator()`](pymatgen.util.convergence.md#pymatgen.util.convergence.id_generator)


    * [`measure()`](pymatgen.util.convergence.md#pymatgen.util.convergence.measure)


    * [`multi_curve_fit()`](pymatgen.util.convergence.md#pymatgen.util.convergence.multi_curve_fit)


    * [`multi_reciprocal_extra()`](pymatgen.util.convergence.md#pymatgen.util.convergence.multi_reciprocal_extra)


    * [`p0_exponential()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_exponential)


    * [`p0_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_reciprocal)


    * [`p0_simple_2reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_simple_2reciprocal)


    * [`p0_simple_4reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_simple_4reciprocal)


    * [`p0_simple_5reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_simple_5reciprocal)


    * [`p0_simple_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_simple_reciprocal)


    * [`p0_single_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.p0_single_reciprocal)


    * [`print_and_raise_error()`](pymatgen.util.convergence.md#pymatgen.util.convergence.print_and_raise_error)


    * [`print_plot_line()`](pymatgen.util.convergence.md#pymatgen.util.convergence.print_plot_line)


    * [`reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.reciprocal)


    * [`simple_2reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.simple_2reciprocal)


    * [`simple_4reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.simple_4reciprocal)


    * [`simple_5reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.simple_5reciprocal)


    * [`simple_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.simple_reciprocal)


    * [`single_reciprocal()`](pymatgen.util.convergence.md#pymatgen.util.convergence.single_reciprocal)


* [pymatgen.util.coord module](pymatgen.util.coord.md)


    * [`Simplex`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex)


        * [`Simplex.bary_coords()`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.bary_coords)


        * [`Simplex.coords`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.coords)


        * [`Simplex.in_simplex()`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.in_simplex)


        * [`Simplex.line_intersection()`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.line_intersection)


        * [`Simplex.point_from_bary_coords()`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.point_from_bary_coords)


        * [`Simplex.volume`](pymatgen.util.coord.md#pymatgen.util.coord.Simplex.volume)


    * [`all_distances()`](pymatgen.util.coord.md#pymatgen.util.coord.all_distances)


    * [`barycentric_coords()`](pymatgen.util.coord.md#pymatgen.util.coord.barycentric_coords)


    * [`coord_list_mapping()`](pymatgen.util.coord.md#pymatgen.util.coord.coord_list_mapping)


    * [`coord_list_mapping_pbc()`](pymatgen.util.coord.md#pymatgen.util.coord.coord_list_mapping_pbc)


    * [`find_in_coord_list()`](pymatgen.util.coord.md#pymatgen.util.coord.find_in_coord_list)


    * [`find_in_coord_list_pbc()`](pymatgen.util.coord.md#pymatgen.util.coord.find_in_coord_list_pbc)


    * [`get_angle()`](pymatgen.util.coord.md#pymatgen.util.coord.get_angle)


    * [`get_linear_interpolated_value()`](pymatgen.util.coord.md#pymatgen.util.coord.get_linear_interpolated_value)


    * [`in_coord_list()`](pymatgen.util.coord.md#pymatgen.util.coord.in_coord_list)


    * [`in_coord_list_pbc()`](pymatgen.util.coord.md#pymatgen.util.coord.in_coord_list_pbc)


    * [`is_coord_subset()`](pymatgen.util.coord.md#pymatgen.util.coord.is_coord_subset)


    * [`is_coord_subset_pbc()`](pymatgen.util.coord.md#pymatgen.util.coord.is_coord_subset_pbc)


    * [`lattice_points_in_supercell()`](pymatgen.util.coord.md#pymatgen.util.coord.lattice_points_in_supercell)


    * [`pbc_diff()`](pymatgen.util.coord.md#pymatgen.util.coord.pbc_diff)


    * [`pbc_shortest_vectors()`](pymatgen.util.coord.md#pymatgen.util.coord.pbc_shortest_vectors)


* [pymatgen.util.coord_cython module](pymatgen.util.coord_cython.md)


    * [`coord_list_mapping_pbc()`](pymatgen.util.coord_cython.md#pymatgen.util.coord_cython.coord_list_mapping_pbc)


    * [`is_coord_subset_pbc()`](pymatgen.util.coord_cython.md#pymatgen.util.coord_cython.is_coord_subset_pbc)


    * [`pbc_shortest_vectors()`](pymatgen.util.coord_cython.md#pymatgen.util.coord_cython.pbc_shortest_vectors)


* [pymatgen.util.due module](pymatgen.util.due.md)


    * [`BibTeX()`](pymatgen.util.due.md#pymatgen.util.due.BibTeX)


    * [`Doi()`](pymatgen.util.due.md#pymatgen.util.due.Doi)


    * [`InactiveDueCreditCollector`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector)


        * [`InactiveDueCreditCollector.activate()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.activate)


        * [`InactiveDueCreditCollector.active`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.active)


        * [`InactiveDueCreditCollector.add()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.add)


        * [`InactiveDueCreditCollector.cite()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.cite)


        * [`InactiveDueCreditCollector.dcite()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.dcite)


        * [`InactiveDueCreditCollector.dump()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.dump)


        * [`InactiveDueCreditCollector.load()`](pymatgen.util.due.md#pymatgen.util.due.InactiveDueCreditCollector.load)


    * [`Text()`](pymatgen.util.due.md#pymatgen.util.due.Text)


    * [`Url()`](pymatgen.util.due.md#pymatgen.util.due.Url)


* [pymatgen.util.graph_hashing module](pymatgen.util.graph_hashing.md)


    * [`weisfeiler_lehman_graph_hash()`](pymatgen.util.graph_hashing.md#pymatgen.util.graph_hashing.weisfeiler_lehman_graph_hash)


    * [`weisfeiler_lehman_subgraph_hashes()`](pymatgen.util.graph_hashing.md#pymatgen.util.graph_hashing.weisfeiler_lehman_subgraph_hashes)


* [pymatgen.util.io_utils module](pymatgen.util.io_utils.md)


    * [`clean_lines()`](pymatgen.util.io_utils.md#pymatgen.util.io_utils.clean_lines)


    * [`micro_pyawk()`](pymatgen.util.io_utils.md#pymatgen.util.io_utils.micro_pyawk)


* [pymatgen.util.num module](pymatgen.util.num.md)


    * [`abs_cap()`](pymatgen.util.num.md#pymatgen.util.num.abs_cap)


    * [`make_symmetric_matrix_from_upper_tri()`](pymatgen.util.num.md#pymatgen.util.num.make_symmetric_matrix_from_upper_tri)


    * [`maxloc()`](pymatgen.util.num.md#pymatgen.util.num.maxloc)


    * [`min_max_indexes()`](pymatgen.util.num.md#pymatgen.util.num.min_max_indexes)


    * [`minloc()`](pymatgen.util.num.md#pymatgen.util.num.minloc)


    * [`monotonic()`](pymatgen.util.num.md#pymatgen.util.num.monotonic)


    * [`non_decreasing()`](pymatgen.util.num.md#pymatgen.util.num.non_decreasing)


    * [`non_increasing()`](pymatgen.util.num.md#pymatgen.util.num.non_increasing)


    * [`round_to_sigfigs()`](pymatgen.util.num.md#pymatgen.util.num.round_to_sigfigs)


    * [`strictly_decreasing()`](pymatgen.util.num.md#pymatgen.util.num.strictly_decreasing)


    * [`strictly_increasing()`](pymatgen.util.num.md#pymatgen.util.num.strictly_increasing)


* [pymatgen.util.numba module](pymatgen.util.numba.md)


    * [`jit()`](pymatgen.util.numba.md#pymatgen.util.numba.jit)


    * [`njit()`](pymatgen.util.numba.md#pymatgen.util.numba.njit)


* [pymatgen.util.plotting module](pymatgen.util.plotting.md)


    * [`add_fig_kwargs()`](pymatgen.util.plotting.md#pymatgen.util.plotting.add_fig_kwargs)


    * [`format_formula()`](pymatgen.util.plotting.md#pymatgen.util.plotting.format_formula)


    * [`get_ax3d_fig_plt()`](pymatgen.util.plotting.md#pymatgen.util.plotting.get_ax3d_fig_plt)


    * [`get_ax_fig_plt()`](pymatgen.util.plotting.md#pymatgen.util.plotting.get_ax_fig_plt)


    * [`get_axarray_fig_plt()`](pymatgen.util.plotting.md#pymatgen.util.plotting.get_axarray_fig_plt)


    * [`periodic_table_heatmap()`](pymatgen.util.plotting.md#pymatgen.util.plotting.periodic_table_heatmap)


    * [`pretty_plot()`](pymatgen.util.plotting.md#pymatgen.util.plotting.pretty_plot)


    * [`pretty_plot_two_axis()`](pymatgen.util.plotting.md#pymatgen.util.plotting.pretty_plot_two_axis)


    * [`pretty_polyfit_plot()`](pymatgen.util.plotting.md#pymatgen.util.plotting.pretty_polyfit_plot)


    * [`van_arkel_triangle()`](pymatgen.util.plotting.md#pymatgen.util.plotting.van_arkel_triangle)


* [pymatgen.util.provenance module](pymatgen.util.provenance.md)


    * [`Author`](pymatgen.util.provenance.md#pymatgen.util.provenance.Author)


        * [`Author.as_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.Author.as_dict)


        * [`Author.from_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.Author.from_dict)


        * [`Author.parse_author()`](pymatgen.util.provenance.md#pymatgen.util.provenance.Author.parse_author)


    * [`HistoryNode`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode)


        * [`HistoryNode.name`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.name)


        * [`HistoryNode.url`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.url)


        * [`HistoryNode.description`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.description)


        * [`HistoryNode.as_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.as_dict)


        * [`HistoryNode.from_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.from_dict)


        * [`HistoryNode.parse_history_node()`](pymatgen.util.provenance.md#pymatgen.util.provenance.HistoryNode.parse_history_node)


    * [`StructureNL`](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL)


        * [`StructureNL.as_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL.as_dict)


        * [`StructureNL.from_dict()`](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL.from_dict)


        * [`StructureNL.from_structures()`](pymatgen.util.provenance.md#pymatgen.util.provenance.StructureNL.from_structures)


    * [`is_valid_bibtex()`](pymatgen.util.provenance.md#pymatgen.util.provenance.is_valid_bibtex)


* [pymatgen.util.string module](pymatgen.util.string.md)


    * [`Stringify`](pymatgen.util.string.md#pymatgen.util.string.Stringify)


        * [`Stringify.STRING_MODE`](pymatgen.util.string.md#pymatgen.util.string.Stringify.STRING_MODE)


        * [`Stringify.to_html_string()`](pymatgen.util.string.md#pymatgen.util.string.Stringify.to_html_string)


        * [`Stringify.to_latex_string()`](pymatgen.util.string.md#pymatgen.util.string.Stringify.to_latex_string)


        * [`Stringify.to_pretty_string()`](pymatgen.util.string.md#pymatgen.util.string.Stringify.to_pretty_string)


        * [`Stringify.to_unicode_string()`](pymatgen.util.string.md#pymatgen.util.string.Stringify.to_unicode_string)


    * [`charge_string()`](pymatgen.util.string.md#pymatgen.util.string.charge_string)


    * [`disordered_formula()`](pymatgen.util.string.md#pymatgen.util.string.disordered_formula)


    * [`formula_double_format()`](pymatgen.util.string.md#pymatgen.util.string.formula_double_format)


    * [`htmlify()`](pymatgen.util.string.md#pymatgen.util.string.htmlify)


    * [`latexify()`](pymatgen.util.string.md#pymatgen.util.string.latexify)


    * [`latexify_spacegroup()`](pymatgen.util.string.md#pymatgen.util.string.latexify_spacegroup)


    * [`str_delimited()`](pymatgen.util.string.md#pymatgen.util.string.str_delimited)


    * [`stream_has_colours()`](pymatgen.util.string.md#pymatgen.util.string.stream_has_colours)


    * [`transformation_to_string()`](pymatgen.util.string.md#pymatgen.util.string.transformation_to_string)


    * [`unicodeify()`](pymatgen.util.string.md#pymatgen.util.string.unicodeify)


    * [`unicodeify_spacegroup()`](pymatgen.util.string.md#pymatgen.util.string.unicodeify_spacegroup)


    * [`unicodeify_species()`](pymatgen.util.string.md#pymatgen.util.string.unicodeify_species)


* [pymatgen.util.testing module](pymatgen.util.testing.md)


    * [`PymatgenTest`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest)


        * [`PymatgenTest.MODULE_DIR`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.MODULE_DIR)


        * [`PymatgenTest.STRUCTURES_DIR`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.STRUCTURES_DIR)


        * [`PymatgenTest.TEST_FILES_DIR`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.TEST_FILES_DIR)


        * [`PymatgenTest.TEST_STRUCTURES`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.TEST_STRUCTURES)


        * [`PymatgenTest.assert_all_close()`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.assert_all_close)


        * [`PymatgenTest.assert_msonable()`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.assert_msonable)


        * [`PymatgenTest.assert_str_content_equal()`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.assert_str_content_equal)


        * [`PymatgenTest.fn`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.fn)


        * [`PymatgenTest.get_structure()`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.get_structure)


        * [`PymatgenTest.serialize_with_pickle()`](pymatgen.util.testing.md#pymatgen.util.testing.PymatgenTest.serialize_with_pickle)


* [pymatgen.util.typing module](pymatgen.util.typing.md)