"""Tools and helper functions for abinit calculations"""
import os.path
import collections 

from pymatgen.util.string_utils import list_strings, StringColorizer


##########################################################################################
# Helper functions.

def remove_trees(*paths, **kwargs):
    import shutil
    ignore_errors = kwargs.pop("ignore_errors", True)
    onerror       = kwargs.pop("onerror", None)
                                                                           
    for path in paths:
        shutil.rmtree(path, ignore_errors=ignore_errors, onerror=onerror)

##########################################################################################

class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """
    def __init__(self, basename, dirname="."):
        self.basename = basename
        self.dirname = os.path.abspath(dirname)
        self.path = os.path.join(self.dirname, self.basename)

    def __str__(self):
       return self.read()

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        "True if self is a NetCDF file"
        return self.basename.endswith(".nc")

    def read(self):
        with open(self.path, "r") as f: return f.read()

    def readlines(self):
        with open(self.path, "r") as f: return f.readlines()

    def write(self, string):
        with open(self.path, "w") as f: return f.write(string)
                                        
    def writelines(self, lines):
        with open(self.path, "w") as f: return f.writelines()

##########################################################################################

class _ewc_tuple(collections.namedtuple("ewc_tuple", "errors, warnings, comments")):

    @property
    def num_errors(self):
        "Number of error messages"
        return len(self.errors)

    @property
    def num_warnings(self):
        "Number of warning messages"
        return len(self.warnings)

    @property
    def num_comments(self):
        "Number of comments"
        return len(self.comments)

    def tostream(self, stream):
        "Return a string that can be visualized on stream (with colors if stream support them)."
        str_colorizer = StringColorizer(stream)

        red  = lambda num : str_colorizer(str(num), "red")
        blue = lambda num : str_colorizer(str(num), "blue")

        nums = map(len, [self.errors, self.warnings, self.comments])

        colors = (red, blue, str)

        for (i, n) in enumerate(nums): 
            color = colors[i]
            nums[i] = color(n) if n else str(n)

        return "%s errors, %s warnings, %s comments in main output file" % tuple(nums)

##########################################################################################

def parse_ewc(filename, nafter=5):
    """
    Extract errors, warnings and comments from file filename.
                                                                                           
    :arg nafter: Save nafter lines of trailing context after matching lines.  
    :return: namedtuple instance. The lists of strings with the corresponding messages are
             available in tuple.errors, tuple.warnings, tuple.comments. 
    """
    # TODO
    # we have to standardize the abinit WARNING, COMMENT and ERROR  so that we can parse them easily
    # without having to use nafter.

    errors, warnings, comments = [], [], []
    # Note the space after the name.
    exc_cases = ["ERROR ", "BUG ", "WARNING ", "COMMENT "]

    handlers = {
        "ERROR "   : errors.append,
        "BUG "     : errors.append,
        "WARNING " : warnings.append,
        "COMMENT " : comments.append,
    }

    def exc_case(line):
        for e in exc_cases:
            if e in line: return e
        else:
            return None

    with open(filename, "r") as fh:      
        lines = fh.readlines()
        nlines = len(lines)
        for (lineno, line) in enumerate(lines):
            handle = handlers.get(exc_case(line))
            if handle is None: continue
            context = lines[lineno: min(lineno+nafter, nlines)]
            handle( "".join([c for c in context]) )

    return _ewc_tuple(errors, warnings, comments)

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
        if dataset is not None: found = found and "DS" +  str(dataset) in tokens
        if image is not None:   found = found and "IMG" + str(image)   in tokens
        if found: return filename
  else:
      return None

##########################################################################################

def abinit_output_iscomplete(output_file):
    "Return True if the abinit output file is complete."
    if not os.path.exists(output_file):
        return False

    chunk = 5 * 1024  # Read only the last 5Kb of data.
    nlines = 10       # Check only in the last 10 lines.

    with open(output_file, 'r') as f:
        size = f.tell()
        f.seek(max(size - chunk, 0))
        try:
            for line in f.read().splitlines()[-nlines:]:
                if "Calculation completed." in line:
                    return True
        except:
            pass
    return False
