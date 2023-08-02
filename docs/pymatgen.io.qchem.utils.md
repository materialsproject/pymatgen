---
layout: default
title: pymatgen.io.qchem.utils.md
nav_exclude: true
---

# pymatgen.io.qchem.utils module

Utilities for Qchem io.


### pymatgen.io.qchem.utils.lower_and_check_unique(dict_to_check)
Takes a dictionary and makes all the keys lower case. Also converts all numeric
values (floats, ints) to str and replaces “jobtype” with “job_type” just so that
key specifically can be called elsewhere without ambiguity. Finally, ensures that
multiple identical keys, that differed only due to different capitalizations, are not
present. If there are multiple equivalent keys, an Exception is raised.


* **Parameters**

    **dict_to_check** (*dict*) – The dictionary to check and standardize



* **Returns**

    An identical dictionary but with all keys made

        lower case and no identical keys present.




* **Return type**

    to_return (dict)



### pymatgen.io.qchem.utils.process_parsed_HESS(hess_data)
Takes the information contained in a HESS file and converts it into
the format of the machine-readable 132.0 file which can be printed
out to be read into subsequent optimizations.


### pymatgen.io.qchem.utils.process_parsed_coords(coords)
Takes a set of parsed coordinates, which come as an array of strings,
and returns a numpy array of floats.


### pymatgen.io.qchem.utils.process_parsed_fock_matrix(fock_matrix)
The Fock matrix is parsed as a list, while it should actually be
a square matrix, this function takes the list of finds the right dimensions
in order to reshape the matrix.


### pymatgen.io.qchem.utils.read_matrix_pattern(header_pattern, footer_pattern, elements_pattern, text, postprocess=<class 'str'>)
Parse a matrix to get the quantities in a numpy array.


### pymatgen.io.qchem.utils.read_pattern(text_str, patterns, terminate_on_match=False, postprocess=<class 'str'>)
General pattern reading on an input string.


* **Parameters**


    * **text_str** (*str*) – the input string to search for patterns


    * **patterns** (*dict*) – A dict of patterns, e.g.,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”}.


    * **terminate_on_match** (*bool*) – Whether to terminate when there is at
    least one match in each key in pattern.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


Renders accessible:

    Any attribute in patterns. For example,
    {“energy”: r”energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)”} will set the
    value of matches[“energy”] = [[-1234], [-3453], …], to the
    results from regex and postprocess. Note that the returned values
    are lists of lists, because you can grep multiple items on one line.


### pymatgen.io.qchem.utils.read_table_pattern(text_str, header_pattern, row_pattern, footer_pattern, postprocess=<class 'str'>, attribute_name=None, last_one_only=False)
Parse table-like data. A table composes of three parts: header,
main body, footer. All the data matches “row pattern” in the main body
will be returned.


* **Parameters**


    * **text_str** (*str*) – the input string to search for patterns


    * **header_pattern** (*str*) – The regular expression pattern matches the
    table header. This pattern should match all the text
    immediately before the main body of the table. For multiple
    sections table match the text until the section of
    interest. MULTILINE and DOTALL options are enforced, as a
    result, the “.” meta-character will also match “n” in this
    section.


    * **row_pattern** (*str*) – The regular expression matches a single line in
    the table. Capture interested field using regular expression
    groups.


    * **footer_pattern** (*str*) – The regular expression matches the end of the
    table. E.g. a long dash line.


    * **postprocess** (*callable*) – A post processing function to convert all
    matches. Defaults to str, i.e., no change.


    * **attribute_name** (*str*) – Name of this table. If present the parsed data
    will be attached to “data. e.g. self.data[“efg”] = […]


    * **last_one_only** (*bool*) – All the tables will be parsed, if this option
    is set to True, only the last table will be returned. The
    enclosing list will be removed. i.e. Only a single table will
    be returned. Default to be True.



* **Returns**

    List of tables. 1) A table is a list of rows. 2) A row if either a list of
    attribute values in case the capturing group is defined without name in
    row_pattern, or a dict in case that named capturing groups are defined by
    row_pattern.