"""Tools and helper functions for abinit calculations"""
import os.path

from pymatgen.util.string_utils import list_strings, StringColorizer


class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def basename(self):
        return os.path.basename(self.path)

    @property
    def dirname(self):
        return os.path.dirname(self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        "True if self is a NetCDF file"
        return self.basename.endswith(".nc")

    def read(self):
        with open(self.path, "r") as f:
            return f.read()

    def readlines(self):
        with open(self.path, "r") as f:
            return f.readlines()

    def write(self, string):
        self.make_dir()
        with open(self.path, "w") as f:
            return f.write(string)

    def writelines(self, lines):
        self.make_dir()
        with open(self.path, "w") as f:
            return f.writelines()

    def make_dir(self):
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)

##########################################################################################

def find_file(files, ext, prefix=None, dataset=None, image=None):
    """
    Given a list of file names, return the file with extension "_" + ext, None if not found.

    The prefix, the dataset index and the image index can be specified

    .. warning::

       There are some border cases that will confuse the algorithm
       since the order of dataset and image is not tested.
       Solving this problem requires the knowledge of ndtset and nimages
       This code, however should work in 99.9% of the cases.
    """
    separator = "_"

    for filename in list_strings(files):
        # Remove Netcdf extension (if any)
        f = filename[:-3] if filename.endswith(".nc") else filename
        if separator not in f: continue
        tokens = f.split(separator)
        if tokens[-1] == ext:
            found = True
            if prefix is not None:  found = found and filename.startswith(prefix)
            if dataset is not None: found = found and "DS" + str(dataset) in tokens
            if image is not None:   found = found and "IMG" + str(image) in tokens
            if found: return filename
    else:
        return None

##########################################################################################


def abinit_output_iscomplete(output_file):
    """Return True if the abinit output file is complete."""
    if not os.path.exists(output_file):
        return False

    chunk = 5 * 1024  # Read only the last 5Kb of data.
    nlines = 10       # Check only in the last 10 lines.

    MAGIC = "Calculation completed." 

    with open(output_file, 'r') as f:
        size = f.tell()
        f.seek(max(size - chunk, 0))
        try:
            for line in f.read().splitlines()[-nlines:]:
                if MAGIC in line:
                    return True
        except:
            pass

    return False


class NullFile(object):
    def __init__(self):
        import os
        return open(os.devnull, 'w')


class NullStream(object):
    def write(*args):
        pass
