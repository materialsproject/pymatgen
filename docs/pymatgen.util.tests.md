---
layout: default
title: pymatgen.util.tests.md
nav_exclude: true
---

# pymatgen.util.tests package


## pymatgen.util.tests.test_convergence module


### _class_ pymatgen.util.tests.test_convergence.ConvergenceTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_determine_convergence()
## pymatgen.util.tests.test_coord module


### _class_ pymatgen.util.tests.test_coord.CoordUtilsTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_all_distances()

#### test_barycentric()

#### test_coord_list_mapping()

#### test_coord_list_mapping_pbc()

#### test_find_in_coord_list()

#### test_find_in_coord_list_pbc()

#### test_get_angle()

#### test_get_linear_interpolated_value()

#### test_in_coord_list()

#### test_in_coord_list_pbc()

#### test_is_coord_subset()

#### test_is_coord_subset_pbc()

#### test_lattice_points_in_supercell()

#### test_pbc_diff()

#### test_pbc_shortest_vectors()

### _class_ pymatgen.util.tests.test_coord.SimplexTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_2dtriangle()

#### test_bary_coords()

#### test_equal()

#### test_in_simplex()

#### test_intersection()

#### test_str()

#### test_to_json()

#### test_volume()
## pymatgen.util.tests.test_graph_hashing module

Copyright (C) 2004-2022, NetworkX Developers
Aric Hagberg <[hagberg@lanl.gov](mailto:hagberg@lanl.gov)>
Dan Schult <[dschult@colgate.edu](mailto:dschult@colgate.edu)>
Pieter Swart <[swart@lanl.gov](mailto:swart@lanl.gov)>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

>
> * Redistributions of source code must retain the above copyright
> notice, this list of conditions and the following disclaimer.


> * Redistributions in binary form must reproduce the above
> copyright notice, this list of conditions and the following
> disclaimer in the documentation and/or other materials provided
> with the distribution.


> * Neither the name of the NetworkX Developers nor the names of its
> contributors may be used to endorse or promote products derived
> from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
“AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### pymatgen.util.tests.test_graph_hashing.test_graph_hash()

### pymatgen.util.tests.test_graph_hashing.test_subgraph_hashes()
## pymatgen.util.tests.test_io_utils module


### _class_ pymatgen.util.tests.test_io_utils.FuncTest(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_micro_pyawk()
## pymatgen.util.tests.test_num_utils module


### _class_ pymatgen.util.tests.test_num_utils.FuncTestCase(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_abs_cap()

#### test_min_max_indexes()

#### test_round()
## pymatgen.util.tests.test_plotting module


### _class_ pymatgen.util.tests.test_plotting.FuncTestCase(methodName='runTest')
Bases: [`PymatgenTest`](pymatgen.util.md#pymatgen.util.testing.PymatgenTest)

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_plot_periodic_heatmap()

#### test_van_arkel_triangle()
## pymatgen.util.tests.test_provenance module


### _class_ pymatgen.util.tests.test_provenance.StructureNLCase(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### setUp()
Hook method for setting up the test fixture before exercising it.


#### test_authors()

#### test_data()

#### test_eq()

#### test_from_structures()

#### test_history_nodes()

#### test_references()

#### test_remarks()

#### test_to_from_dict()
## pymatgen.util.tests.test_string module


### _class_ pymatgen.util.tests.test_string.FuncTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_charge_string()

#### test_disordered_formula()

#### test_formula_double_format()

#### test_htmlify()

#### test_latexify()

#### test_latexify_spacegroup()

#### test_transformation_to_string()

#### test_unicodeify()

### _class_ pymatgen.util.tests.test_string.StringifyTest(methodName='runTest')
Bases: `TestCase`

Create an instance of the class that will use the named test
method when executed. Raises a ValueError if the instance does
not have a method with the specified name.


#### test_to_html_string()

#### test_to_latex_string()

#### test_to_unicode_string()

### _class_ pymatgen.util.tests.test_string.SubStr()
Bases: [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)


### _class_ pymatgen.util.tests.test_string.SupStr()
Bases: [`Stringify`](pymatgen.util.md#pymatgen.util.string.Stringify)


#### STRING_MODE(_ = 'SUPERSCRIPT_ )

#### to_pretty_string()

* **Returns**

    A pretty string representation. By default, the __str__ output is used, but this method can be
    overridden if a different representation from default is desired.


## pymatgen.util.tests.test_typing module


### pymatgen.util.tests.test_typing.test_composition_like()

### pymatgen.util.tests.test_typing.test_entry_like()

### pymatgen.util.tests.test_typing.test_matrix_like()

### pymatgen.util.tests.test_typing.test_path_like()

### pymatgen.util.tests.test_typing.test_species_like()