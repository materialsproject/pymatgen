"""Tools and helper functions for abinit calculations"""
import os.path

from pymatgen.util.string_utils import list_strings, StringColorizer


class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def path(self):
        """Absolute path of the file."""
        return self._path

    @property
    def basename(self):
        """File basename."""
        return os.path.basename(self.path)

    @property
    def dirname(self):
        """Absolute path of the directory where the file is located."""
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
        """Read data from file."""
        with open(self.path, "r") as f:
            return f.read()

    def readlines(self):
        """Read lines from files."""
        with open(self.path, "r") as f:
            return f.readlines()

    def write(self, string):
        """Write string to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            return f.write(string)

    def writelines(self, lines):
        """Write a list of strings to file."""
        self.make_dir()
        with open(self.path, "w") as f:
            return f.writelines()

    def make_dir(self):
        """Make the directory where the file is located."""
        if not os.path.exists(self.dirname):
            os.makedirs(self.dirname)


class Directory(object):
    """
    Very simple class that provides helper functions
    wrapping the most commonly used functions of os.path.
    """
    def __init__(self, path):
        self._path = os.path.abspath(path)

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def path(self):
        """Absolute path of the directory."""
        return self._path

    @property
    def basename(self):
        """Directory basename."""
        return os.path.basename(self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    def makedirs(self):
        """
        Super-mkdir; create a leaf directory and all intermediate ones.
        Works like mkdir, except that any intermediate path segment (not
        just the rightmost) will be created if it does not exist.
        """
        if not self.exists:
            os.makedirs(self.path)

    def path_in_dir(self, filename):
        """Return the absolute path of filename in the directory."""
        return os.path.join(self.path, filename)

    def list_filepaths(self):
        """Return the list of absolute filepaths in the top-level directory."""
        fnames = [f for f in os.listdir(self.path)]
        return [os.path.join(self.path, f) for f in fnames]

    def has_abifile(self, ext):
        """
        Returns the absolute path of the ABINIT file with extension ext.
        Support both Fortran files and netcdf files. In the later case,
        we check whether a file with extension ext + ".nc" is present 
        in the directory.
        Returns empty string is file is not present.

        Raises:
            ValueError if multiple files with the given ext are found.
            This implies that this method is not compatible with multiple datasets.
        """
        files = []
        for f in self.list_filepaths():
            if f.endswith(ext) or f.endswith(ext + ".nc"):
                files.append(f)

        if not files:
            return ""

        if len(files) > 1:
            # ABINIT users must learn that multiple datasets are bad!
            err_msg = "Found multiple files with the same extensions\n Please avoid the use of mutiple datasets!"
            raise ValueError(err_msg)

        return files[0]


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
