import re
import numpy as np
from collections import defaultdict


def read_pattern(text_str, patterns, terminate_on_match=False,
                 postprocess=str):
    """
        General pattern reading on an input string

        Args:
            text_str (str): the input string to search for patterns
            patterns (dict): A dict of patterns, e.g.,
                {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"}.
            terminate_on_match (bool): Whether to terminate when there is at
                least one match in each key in pattern.
            postprocess (callable): A post processing function to convert all
                matches. Defaults to str, i.e., no change.

        Renders accessible:
            Any attribute in patterns. For example,
            {"energy": r"energy\\(sigma->0\\)\\s+=\\s+([\\d\\-.]+)"} will set the
            value of matches["energy"] = [[-1234], [-3453], ...], to the
            results from regex and postprocess. Note that the returned values
            are lists of lists, because you can grep multiple items on one line.
    """

    compiled = {
        key: re.compile(pattern, re.MULTILINE | re.DOTALL)
        for key, pattern in patterns.items()
    }
    matches = defaultdict(list)
    for key, pattern in compiled.items():
        for match in pattern.finditer(text_str):
            matches[key].append([postprocess(i) for i in match.groups()])
            if terminate_on_match:
                break
    return matches


def read_table_pattern(text_str,
                       header_pattern,
                       row_pattern,
                       footer_pattern,
                       postprocess=str,
                       attribute_name=None,
                       last_one_only=False):
    """
    Parse table-like data. A table composes of three parts: header,
    main body, footer. All the data matches "row pattern" in the main body
    will be returned.

    Args:
        text_str (str): the input string to search for patterns
        header_pattern (str): The regular expression pattern matches the
            table header. This pattern should match all the text
            immediately before the main body of the table. For multiple
            sections table match the text until the section of
            interest. MULTILINE and DOTALL options are enforced, as a
            result, the "." meta-character will also match "\n" in this
            section.
        row_pattern (str): The regular expression matches a single line in
            the table. Capture interested field using regular expression
            groups.
        footer_pattern (str): The regular expression matches the end of the
            table. E.g. a long dash line.
        postprocess (callable): A post processing function to convert all
            matches. Defaults to str, i.e., no change.
        attribute_name (str): Name of this table. If present the parsed data
            will be attached to "data. e.g. self.data["efg"] = [...]
        last_one_only (bool): All the tables will be parsed, if this option
            is set to True, only the last table will be returned. The
            enclosing list will be removed. i.e. Only a single table will
            be returned. Default to be True.

    Returns:
        List of tables. 1) A table is a list of rows. 2) A row if either a list of
        attribute values in case the the capturing group is defined without name in
        row_pattern, or a dict in case that named capturing groups are defined by
        row_pattern.
    """

    table_pattern_text = header_pattern + \
        r"\s*(?P<table_body>(?:" + row_pattern + r")+)\s*" + footer_pattern
    table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
    rp = re.compile(row_pattern)
    data = {}
    tables = []
    for mt in table_pattern.finditer(text_str):
        table_body_text = mt.group("table_body")
        table_contents = []
        for ml in rp.finditer(table_body_text):
            d = ml.groupdict()
            if len(d) > 0:
                processed_line = {k: postprocess(v) for k, v in d.items()}
            else:
                processed_line = [postprocess(v) for v in ml.groups()]
            table_contents.append(processed_line)
        tables.append(table_contents)
    if last_one_only:
        retained_data = tables[-1]
    else:
        retained_data = tables
    if attribute_name is not None:
        data[attribute_name] = retained_data
        return data
    return retained_data


def lower_and_check_unique(dict_to_check):
    """
    Takes a dictionary and makes all the keys lower case. Also replaces
    "jobtype" with "job_type" just so that key specifically can be called
    elsewhere without ambiguity. Finally, ensures that multiple identical
    keys, that differed only due to different capitalizations, are not
    present. If there are multiple equivalent keys, an Exception is raised.

    Args:
        dict_to_check (dict): The dictionary to check and standardize

    Returns:
        to_return (dict): An identical dictionary but with all keys made
            lower case and no identical keys present.
    """
    if dict_to_check == None:
        return None
    else:
        to_return = {}
        for key in dict_to_check:
            new_key = key.lower()
            if new_key == "jobtype":
                new_key = "job_type"
            if new_key in to_return:
                raise Exception(
                    "Multiple instances of key " + new_key + " found!")
            else:
                try:
                    to_return[new_key] = dict_to_check.get(key).lower()
                except AttributeError:
                    to_return[new_key] = dict_to_check.get(key)
        return to_return


def process_parsed_coords(coords):
    """
    Takes a set of parsed coordinates, which come as an array of strings,
    and returns a numpy array of floats.
    """
    geometry = np.zeros(shape=(len(coords), 3), dtype=float)
    for ii, entry in enumerate(coords):
        for jj in range(3):
            geometry[ii, jj] = float(entry[jj])
    return geometry
