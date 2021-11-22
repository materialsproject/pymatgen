import os
import pprint


# thanks https://stackoverflow.com/a/36693250
def package_files(directory, extensions):
    """
    Walk package directory to make sure we include all relevant files in
    package.
    """
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            if any([filename.endswith(ext) for ext in extensions]):
                paths.append(os.path.join("..", path, filename))
    return paths


json_yaml_csv_files = package_files("pymatgen", ["yaml", "json", "csv", "yaml.gz", "json.gz", "csv.gz"])

pprint.pprint(json_yaml_csv_files)
