# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides objects to inspect the status of the Abinit tasks at run-time.
by extracting information from the main output file (text format).
"""
from __future__ import unicode_literals, division, print_function

import collections
import numpy as np
import yaml
import six

from six.moves import cStringIO, map, zip
from prettytable import PrettyTable
from monty.collections import AttrDict
from pymatgen.util.plotting_utils import add_fig_kwargs


def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


def _magic_parser(stream, magic):
    """
    Parse the section with the SCF cycle

    Returns:
        dict where the key are the name of columns and
        the values are list of numbers. Note if no section was found.

    .. warning::

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

            tokens = list(map(float, line.split()[1:]))
            assert len(tokens) == len(keys)
            for l, v in zip(fields.values(), tokens):
                l.append(v)

    # Convert values to numpy arrays.
    if fields:
        return collections.OrderedDict([(k, np.array(v)) for k, v in fields.items()])
    else:
        return None


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

    .. attributes::

        num_iterations: Number of iterations performed.
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
        table = PrettyTable(list(self.keys()))
        for it in range(self.num_iterations):
            row = list(map(str, (self[k][it] for k in self.keys())))
            table.add_row(row)

        stream = cStringIO()
        print(table, file=stream)
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
            None if no `ScfCycle` entry is found.
        """
        fields = _magic_parser(stream, magic=cls.MAGIC)

        if fields:
            fields.pop("iter")
            return cls(fields)
        else:
            return None

    @add_fig_kwargs
    def plot(self, fig=None, **kwargs):
        """
        Uses matplotlib to plot the evolution of the SCF cycle. Return `matplotlib` figure

        Args:
            fig: matplotlib figure. If None a new figure is produced.
        """
        import matplotlib.pyplot as plt

        # Build grid of plots.
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = (num_plots//ncols) + (num_plots % ncols)

        if fig is None:
            fig, ax_list = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, squeeze=False)
            ax_list = ax_list.ravel()
        else:
            ax_list = list(fig.axes)

        # Use log scale for these variables.
        use_logscale = set(["residm", "vres2"])

        # Hard-coded y-range for selected variables.
        has_yrange = {
            "deltaE(h)": (-1e-3, +1e-3),
            "deltaE(Ha)": (-1e-3, +1e-3),
            }

        iter_num = np.array(list(range(self.num_iterations)))
        for ((key, values), ax) in zip(self.items(), ax_list):
            ax.grid(True)
            ax.set_xlabel('Iteration')
            ax.set_xticks(iter_num, minor=False)
            ax.set_ylabel(key)

            xx, yy = iter_num, values
            if self.num_iterations > 1:
                # Don't show the first iteration since it's not very useful.
                xx, yy = xx[1:] + 1, values[1:]

            ax.plot(xx, yy, "-o", lw=2.0)

            if key in use_logscale and np.all(yy > 1e-22):
                ax.set_yscale("log")

            if key in has_yrange:
                ymin, ymax = has_yrange[key]
                val_min, val_max = np.min(yy), np.max(yy)
                if abs(val_max - val_min) > abs(ymax - ymin):
                    ax.set_ylim(ymin, ymax)

        # Get around a bug in matplotlib.
        if (num_plots % ncols) != 0:
            ax_list[-1].plot(xx, yy, lw=0.0)
            ax_list[-1].axis('off')

        #plt.legend(loc="best")
        return fig


class GroundStateScfCycle(ScfCycle):
    """Result of the Ground State self-consistent cycle."""
    MAGIC = "iter   Etot(hartree)"

    @property
    def last_etotal(self):
        """The total energy at the last iteration."""
        return self["Etot(hartree)"][-1]


class D2DEScfCycle(ScfCycle):
    """Result of the Phonon self-consistent cycle."""
    MAGIC = "iter   2DEtotal(Ha)"

    @property
    def last_etotal(self):
        """The 2-nd order derivative of the energy at the last iteration."""
        return self["2DEtotal(Ha)"][-1]


class PhononScfCycle(D2DEScfCycle):
    """Iterations of the DFPT SCF cycle for phonons."""


class Relaxation(collections.Iterable):
    """
    A list of :class:`GroundStateScfCycle` objects.

    .. note::

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

    def __getitem__(self, slice):
        return self.cycles[slice]

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
        Dictionary of lists with the evolution of the data as function of the relaxation step.
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

    @add_fig_kwargs
    def plot(self, **kwargs):
        """
        Uses matplotlib to plot the evolution of the structural relaxation.

        Returns:
            `matplotlib` figure
        """
        import matplotlib.pyplot as plt

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

        for (key, values), ax in zip(history.items(), ax_list):
            ax.grid(True)
            ax.set_xlabel('Relaxation Step')
            ax.set_xticks(relax_step, minor=False)
            ax.set_ylabel(key)

            ax.plot(relax_step, values, "-o", lw=2.0)

        return fig


# TODO
#class DfptScfCycle(collections.Iterable):

# TODO
#class HaydockIterations(collections.Iterable):
#    """This object collects info on the different steps of the Haydock technique used in the Bethe-Salpeter code"""
#    @classmethod
#    def from_file(cls, filepath):
#        """Initialize the object from file."""
#        with open(filepath, "r") as stream:
#            return cls.from_stream(stream)
#
#    @classmethod
#    def from_stream(cls, stream):
#        """Extract data from stream. Returns None if some error occurred."""
#        cycles = []
#        while True:
#            scf_cycle = GroundStateScfCycle.from_stream(stream)
#            if scf_cycle is None: break
#            cycles.append(scf_cycle)
#
#        return cls(cycles) if cycles else None
#
#    #def __init__(self):
#
#    def plot(self, **kwargs):
#        """
#        Uses matplotlib to plot the evolution of the structural relaxation.
#        ==============  ==============================================================
#        kwargs          Meaning
#        ==============  ==============================================================
#        title           Title of the plot (Default: None).
#        how            True to show the figure (Default).
#        savefig         'abc.png' or 'abc.eps'* to save the figure to a file.
#        ==============  ==============================================================
#        Returns:
#            `matplotlib` figure
#        """
#        import matplotlib.pyplot as plt
#        title = kwargs.pop("title", None)
#        show = kwargs.pop("show", True)
#        savefig = kwargs.pop("savefig", None)
#        if title: fig.suptitle(title)
#        if savefig is not None: fig.savefig(savefig)
#        if show: plt.show()
#        return fig



##################
## Yaml parsers.
##################

class YamlTokenizerError(Exception):
    """Exceptions raised by :class:`YamlTokenizer`."""


class YamlTokenizer(collections.Iterator):
    """
    Provides context-manager support so you can use it in a with statement.
    """
    Error = YamlTokenizerError

    def __init__(self, filename):
        # The position inside the file.
        self.linepos = 0 
        self.stream = open(filename, "r")

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        try:
            self.stream.close()
        except:
            print("Exception in YAMLTokenizer.close()")
            print(straceback())

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

        .. warning::

            Assume that the YAML document are closed explicitely with the sentinel '...'
        """
        in_doc, lines, doc_tag = None, [], None

        for i, line in enumerate(self.stream):
            self.linepos += 1
            #print(i, line)

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

        .. warning::

            Assume that all the YAML docs (with the exception of the last one)
            are closed explicitely with the sentinel '...'
        """
        docs = [doc for doc in self]
        self.seek(0)
        return docs

    def next_doc_with_tag(self, doc_tag):
        """
        Returns the next document with the specified tag. Empty string is no doc is found.
        """
        while True:
            try:
                doc = six.advance_iterator(self)
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
        d = yaml.load(doc.text_notag)

        return np.array(d["reduced_coordinates_of_qpoints"])


def yaml_read_irred_perts(filename, doc_tag="!IrredPerts"):
    """Read the list of irreducible perturbations from file."""
    with YamlTokenizer(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        d = yaml.load(doc.text_notag)

        return [AttrDict(**pert) for pert in d["irred_perts"]]
        #return d["irred_perts"]


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
            text: String with the YAML document.
            lineno: The line number where the document is located.
            tag: The YAML tag associate to the document.
        """
        # Sanitize strings: use "ignore" to skip invalid characters in .encode/.decode like
        if isinstance(text, bytes):
            text = text.decode("utf-8", "ignore")
        text = text.rstrip().lstrip()
        self.text = text
        self.lineno = lineno
        if isinstance(tag, bytes):
            tag = tag.decode("utf-8", "ignore")
        self.tag = tag

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

    def as_dict(self):
        """Use Yaml to parse the text (without the tag) and returns a dictionary."""
        return yaml.load(self.text_notag)
