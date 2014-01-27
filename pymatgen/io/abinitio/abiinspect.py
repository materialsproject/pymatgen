"""
This module provides objects to inspect the status of the Abinit tasks at run-time.
by extracting information from the main output file (text format).
"""
from __future__ import division, print_function

import collections
import numpy as np

from StringIO import StringIO
from pymatgen.io.abinitio import myaml
from pymatgen.util.string_utils import pprint_table


def _magic_parser(stream, magic):
    """
    Parse the section with the SCF cycle

    Returns:
        dict where the key are the name of columns and 
        the values are list of numbers. Note if no section was found.

    .. warning:
        
        The parser is very fragile and should be replaced by YAML.
    """
    #Example (SCF cycle, similar format is used for phonons):
    #
    #  iter   Etot(hartree)      deltaE(h)  residm     vres2
    #  ETOT  1  -8.8604027880849    -8.860E+00 2.458E-02 3.748E+00

    #  At SCF step    5       vres2   =  3.53E-08 < tolvrs=  1.00E-06 =>converged.
    in_doc, fields = 0, None

    for line in stream:
        line = line.strip()

        if line.startswith(magic):
            keys = line.split()
            fields = collections.OrderedDict((k, []) for k in keys)

        if fields is not None:
            #print(line)
            in_doc += 1
            if in_doc == 1:
                continue

            # End of the section.
            if not line: break

            tokens = map(float, line.split()[1:])
            assert len(tokens) == len(keys)
            for l, v in zip(fields.values(), tokens):
                l.append(v)

    return fields


def plottable_from_outfile(filepath):
    """
    Factory function that returns a plottable object by inspecting the main output file of abinit
    Returns None if it is not able to detect the class to instantiate.
    """
    # TODO
    # Figure out how to detect the type of calculations
    # without having to parse the input. Possible approach: YAML doc
    #with YamlTokenizer(filepath) as r:
    #    doc = r.next_doc_with_tag("!CalculationType")
    #    d = yaml.load(doc.text_notag)
    #    calc_type = d["calculation_type"]

    #ctype2class = {
    #    "Ground State": GroundStateScfCycle,
    #    "Phonon": PhononScfCycle,
    #    "Relaxation": Relaxation,
    #}
    #obj = ctype2class.get(calc_type, None)

    obj = GroundStateScfCycle

    if obj is not None:
        return obj.from_file(filepath)
    else:
        return None


class ScfCycle(collections.Mapping):
    """
    It essentially consists of a dictionary mapping string
    to list of floats containing the data at the different iterations.
    """
    def __init__(self, fields):
        self.fields = fields
        #print(fields)

        all_lens = [len(lst) for lst in self.values()]
        self.num_iterations = all_lens[0]
        assert all(n == self.num_iterations for n in all_lens)

    def __getitem__(self, slice):
        return self.fields.__getitem__(slice)
        
    def __iter__(self):
        return self.fields.__iter__()
    
    def __len__(self):
        return len(self.fields)

    def __str__(self):
        """String representation."""
        table = [list(self.fields.keys())]
        for it in range(self.num_iterations):
            row = map(str, (self[k][it] for k in self.keys()))
            table.append(row)

        stream = StringIO()
        pprint_table(table, out=stream)
        stream.seek(0)

        return "".join(stream)

    @property
    def last_iteration(self):
        """Returns a dictionary with the values of the last iteration."""
        return {k: v[-1] for k, v in self.items()}

    @classmethod
    def from_file(cls, filepath):
        """Read the first occurrence of ScfCycle from file."""
        with open(filepath, "r") as stream:
            return cls.from_stream(stream)

    @classmethod
    def from_stream(cls, stream):
        """
        Read the first occurrence of ScfCycle from stream.

        Returns:
            None if no ScfCycle entry is found.
        """
        fields = _magic_parser(stream, magic=cls.MAGIC)

        if fields:
            fields.pop("iter")
            return cls(fields)
        else:
            return None
        
    def plot(self, **kwargs):
        """
        Uses matplotlib to plot the evolution of the SCF cycle.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure
        """
        import matplotlib.pyplot as plt

        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        # Build grid of plots.
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        if title:
            fig.suptitle(title)

        iter_num = np.array(range(self.num_iterations))

        for ((key, values), ax) in zip(self.items(), ax_list):
            ax.grid(True)
            ax.set_xlabel('Iteration')       
            ax.set_xticks(iter_num, minor=False)
            ax.set_ylabel(key)

            xx, yy = iter_num, values
            if self.num_iterations > 1:
                # Don't show the first iteration since it's not very useful.
                xx, yy = xx[1:] + 1, values[1:]

            #print("xx ",xx, "yy ",yy)
            ax.plot(xx, yy, "-o", lw=2.0)

        # Get around a bug in matplotlib.
        if (num_plots % ncols) != 0:
            ax_list[-1].plot(xx, yy, lw=0.0)
            ax_list[-1].axis('off')

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)

        return fig


class GroundStateScfCycle(ScfCycle):
    """Result of the Ground State self-consistent cycle."""
    #yaml_tag = '!GroundStateScfCycle'
    MAGIC = "iter   Etot(hartree)"

    @property
    def last_etotal(self):
        """The total energy at the last iteration."""
        return self["Etot(hartree)"][-1]


class PhononScfCycle(ScfCycle):
    """Result of the Phonon self-consistent cycle."""
    #yaml_tag = '!PhononScfCycle'
    MAGIC = "iter   2DEtotal(Ha)"

    @property
    def last_etotal(self):
        """The 2-nd order derivative of the energy at the last iteration."""
        return self["2DEtotal(Ha)"][-1]


class Relaxation(collections.Iterable):
    """
    A list of `GroundStateScfCycle` objects.

    .. note:

        Forces, stresses  and crystal structures are missing.
        Solving this problem would require the standardization
        of the Abinit output file (YAML).
    """
    def __init__(self, cycles):
        self.cycles = cycles

    def __iter__(self):
        return self.cycles.__iter__()

    def __len__(self):
        return self.cycles.__len__()

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append

        for i, cycle in enumerate(self):
            app("")
            app("RELAXATION STEP: %d" % i)
            app(str(cycle))
            app("")

        return "\n".join(lines)

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from the Abinit main output file."""
        with open(filepath, "r") as stream:
            return cls.from_stream(stream)

    @classmethod
    def from_stream(cls, stream):
        """
        Extract data from stream. Returns None if some error occurred.
        """
        cycles = []
        while True:
            scf_cycle = GroundStateScfCycle.from_stream(stream)
            if scf_cycle is None: break
            cycles.append(scf_cycle)

        return cls(cycles) if cycles else None

    @property
    def history(self):
        """
        Dictionary of lists with the evolution of the data 
        as function of the relaxation step.
        """
        try:
            return self._history
        except AttributeError:
            self._history = history = collections.defaultdict(list)
            for cycle in self:
                d = cycle.last_iteration
                for k, v in d.items():
                    history[k].append(v)

            return self._history

    def plot(self, **kwargs):
        """
        Uses matplotlib to plot the evolution of the structural relaxation.

        ==============  ==============================================================
        kwargs          Meaning
        ==============  ==============================================================
        title           Title of the plot (Default: None).
        show            True to show the figure (Default).
        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
        ==============  ==============================================================

        Returns:
            `matplotlib` figure
        """
        import matplotlib.pyplot as plt

        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        history = self.history
        #print(history)

        relax_step = list(range(len(self)))

        # Build grid of plots.
        num_plots, ncols, nrows = len(list(history.keys())), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
        ax_list = ax_list.ravel()

        if (num_plots % ncols) != 0:
            ax_list[-1].axis('off')

        if title:
            fig.suptitle(title)
        
        for ((key, values), ax) in zip(history.items(), ax_list):
            ax.grid(True)
            ax.set_xlabel('Relaxation Step')       
            ax.set_xticks(relax_step, minor=False)
            ax.set_ylabel(key)

            ax.plot(relax_step, values, "-o", lw=2.0)
        
        if show:
            plt.show()
        
        if savefig is not None:
            fig.savefig(savefig)

        return fig


class YamlTokenizer(collections.Iterator):
    """
    Provides context-manager support so you can use it in a with statement. 
    """
    def __init__(self, filename):
        self.stream = open(filename, "r")
        self.linepos = 0 # The position inside the file.

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        self.stream.close()

    def seek(self, offset, whence=0):
        """
        seek(offset[, whence]) -> None.  Move to new file position.

        Argument offset is a byte count.  Optional argument whence defaults to
        0 (offset from start of file, offset should be >= 0); other values are 1
        (move relative to current position, positive or negative), and 2 (move
        relative to end of file, usually negative, although many platforms allow
        seeking beyond the end of a file).  If the file is opened in text mode,
        only offsets returned by tell() are legal.  Use of other offsets causes
        undefined behavior.
        Note that not all file objects are seekable.
        """
        assert offset == 0
        self.linepos = 0
        return self.stream.seek(offset, whence)

    # Python 3 compatibility
    def __next__(self):
        return self.next()

    def next(self):
        """
        Returns the first YAML document in stream.

        .. warning:

            Assume that the YAML document are closed explicitely with the sentinel '...'
        """
        in_doc, lines, doc_tag = None, [], None

        for line in self.stream:
            self.linepos += 1

            if line.startswith("---"):
                # Include only lines in the form:
                #  "--- !tag"
                #  "---"
                # Other lines are spurious.
                in_doc = False
                l = line[3:].strip().lstrip()

                if l.startswith("!"):
                    # "--- !tag"
                    doc_tag = l
                    in_doc = True
                elif not l:
                    #  "---"
                    in_doc = True
                    doc_tag = None

                if in_doc:
                    lineno = self.linepos

            if in_doc:
                lines.append(line)

            if in_doc and line.startswith("..."):
                return YamlDoc(text="".join(lines), lineno=lineno, tag=doc_tag)

        raise StopIteration("Cannot find next YAML document")

    def all_yaml_docs(self):
        """
        Returns a list with all the YAML docs found in stream. 
        Seek the stream before returning.

        .. warning:

            Assume that all the YAML docs (with the exception of the last one) 
            are closed explicitely with the sentinel '...'
        """
        docs = [doc for doc in self]
        self.seek(0)
        return docs

    def next_doc_with_tag(self, doc_tag):
        """
        Returns the next document with the specified tag.
        Empty string is no doc is found.
        """
        while True:
            try:
                doc = self.next()
                if doc.tag == doc_tag:
                    return doc

            except StopIteration:
                raise

    def all_docs_with_tag(self, doc_tag):
        """
        Returns all the documents with the specified tag.
        """
        docs = []

        while True:
            try:
                doc = self.next_doc_with(doc_tag)
                docs.append(doc)

            except StopIteration:
                break

        self.seek(0)
        return docs


def yaml_read_kpoints(filename, doc_tag="!Kpoints"):
    """Read the K-points from file."""
    with YamlTokenizer(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        d = myaml.load(doc.text_notag)

        return np.array(d["reduced_coordinates_of_qpoints"])


def yaml_read_irred_perts(filename, doc_tag="!IrredPerts"):
    """Read the lisr of irreducible perturbations from file."""
    with YamlTokenizer(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        d = myaml.load(doc.text_notag)

        return d["irred_perts"]


class YamlDoc(object):
    """
    Handy object that stores that YAML document, its main tag and the 
    position inside the file.
    """
    __slots__ = [
        "text",
        "lineno",
        "tag",
    ]

    def __init__(self, text, lineno, tag=None):
        """
        Args:
            text:
                String with the YAML document.
            lineno:
                The line number where the document is located.
            tag:
                The YAML tag associate to the document.
        """
        # Sanitize strings: use "ignore" to skip invalid characters in .encode/.decode like
        text = text.decode("utf-8", "ignore")
        text = text.rstrip().lstrip()
        self.text = text
        self.lineno = lineno
        self.tag = tag.decode("utf-8", "ignore") if tag is not None else tag

    def __str__(self):
        return self.text

    def __eq__(self, other):
        if other is None: return False
        return (self.text == other.text and 
                self.lineno == other.lineno and 
                self.tag == other.tag) 

    def __ne__(self, other):
        return not self == other

    @property
    def text_notag(self):
        """
        Returns the YAML text without the tag.
        Useful if we don't have any constructor registered for the tag
        (we used the tag just to locate the document).
        """
        if self.tag is not None:
           return self.text.replace(self.tag, "")
        else:
           return self.text
