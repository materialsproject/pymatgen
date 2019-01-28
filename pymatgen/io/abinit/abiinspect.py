# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides objects to inspect the status of the Abinit tasks at run-time.
by extracting information from the main output file (text format).
"""

import os
import numpy as np
import ruamel.yaml as yaml

try:
    from collections.abc import Mapping, Iterable, Iterator
except ImportError:
    from collections import Mapping, Iterable, Iterator
from collections import OrderedDict
from tabulate import tabulate
from monty.functools import lazy_property
from monty.collections import AttrDict
from pymatgen.util.plotting import add_fig_kwargs, get_axarray_fig_plt


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
            fields = OrderedDict((k, []) for k in keys)

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
        return OrderedDict([(k, np.array(v)) for k, v in fields.items()])
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
    #    d = yaml.safe_load(doc.text_notag)
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

# Use log scale for these variables.
_VARS_SUPPORTING_LOGSCALE = set(["residm", "vres2", "nres2"])

# Hard-coded y-range for selected variables.
_VARS_WITH_YRANGE = {
    "deltaE(h)": (-1e-3, +1e-3),
    "deltaE(Ha)": (-1e-3, +1e-3),
}

class ScfCycle(Mapping):
    """
    It essentially consists of a dictionary mapping string
    to list of floats containing the data at the different iterations.

    .. attributes::

        num_iterations: Number of iterations performed.
    """
    def __init__(self, fields):
        self.fields = fields
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
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        rows = [[it + 1] + list(map(str, (self[k][it] for k in self.keys())))
                for it in range(self.num_iterations)]

        return tabulate(rows, headers=["Iter"] + list(self.keys()))

    @property
    def last_iteration(self):
        """Returns a dictionary with the values of the last iteration."""
        return {k: v[-1] for k, v in self.items()}

    @classmethod
    def from_file(cls, filepath):
        """Read the first occurrence of ScfCycle from file."""
        with open(filepath, "rt") as stream:
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
    def plot(self, ax_list=None, fontsize=12, **kwargs):
        """
        Uses matplotlib to plot the evolution of the SCF cycle.

        Args:
            ax_list: List of axes. If None a new figure is produced.
            fontsize: legend fontsize.
            kwargs: keyword arguments are passed to ax.plot

        Returns: matplotlib figure
        """
        # Build grid of plots.
        num_plots, ncols, nrows = len(self), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_list, fig, plot = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                 sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        iter_num = np.array(list(range(self.num_iterations))) + 1
        label = kwargs.pop("label", None)

        for i, ((key, values), ax) in enumerate(zip(self.items(), ax_list)):
            ax.grid(True)
            ax.set_xlabel('Iteration Step')
            ax.set_xticks(iter_num, minor=False)
            ax.set_ylabel(key)

            xx, yy = iter_num, values
            if self.num_iterations > 1:
                # Don't show the first iteration since it's not very useful.
                xx, yy = xx[1:], values[1:]

            if not kwargs and label is None:
                ax.plot(xx, yy, "-o", lw=2.0)
            else:
                ax.plot(xx, yy, label=label if i == 0 else None, **kwargs)

            if key in _VARS_SUPPORTING_LOGSCALE and np.all(yy > 1e-22):
                ax.set_yscale("log")

            if key in _VARS_WITH_YRANGE:
                ymin, ymax = _VARS_WITH_YRANGE[key]
                val_min, val_max = np.min(yy), np.max(yy)
                if abs(val_max - val_min) > abs(ymax - ymin):
                    ax.set_ylim(ymin, ymax)

            if label is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        # Get around a bug in matplotlib.
        if num_plots % ncols != 0:
            ax_list[-1].plot(xx, yy, lw=0.0)
            ax_list[-1].axis('off')

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


class CyclesPlotter:
    """Relies on the plot method of cycle objects to build multiple subfigures."""

    def __init__(self):
        self.labels = []
        self.cycles = []

    def items(self):
        """To iterate over (label, cycle)."""
        return zip(self.labels, self.cycles)

    def add_label_cycle(self, label, cycle):
        """Add new cycle to the plotter with label `label`."""
        self.labels.append(label)
        self.cycles.append(cycle)

    @add_fig_kwargs
    def combiplot(self, fontsize=8, **kwargs):
        """
        Compare multiple cycels on a grid: one subplot per quantity,
        all cycles on the same subplot.

        Args:
            fontsize: Legend fontsize.
        """
        ax_list = None
        for i, (label, cycle) in enumerate(self.items()):
            fig = cycle.plot(ax_list=ax_list, label=label, fontsize=fontsize,
                             lw=2.0, marker="o", linestyle="-", show=False)
            ax_list = fig.axes

        return fig

    def slideshow(self, **kwargs):
        """
        Produce slides show of the different cycles. One plot per cycle.
        """
        for label, cycle in self.items():
            cycle.plot(title=label, tight_layout=True)


class Relaxation(Iterable):
    """
    A list of :class:`GroundStateScfCycle` objects.

    .. attributes::

        num_iterations: Number of iterations performed.

    .. note::

        Forces, stresses  and crystal structures are missing.
        This object is mainly used to analyze the behavior of the Scf cycles
        during the structural relaxation. A more powerful and detailed analysis
        can be obtained by using the HIST.nc file.
    """
    def __init__(self, cycles):
        self.cycles = cycles
        self.num_iterations = len(self.cycles)

    def __iter__(self):
        return self.cycles.__iter__()

    def __len__(self):
        return self.cycles.__len__()

    def __getitem__(self, slice):
        return self.cycles[slice]

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose=0):
        """String representation."""
        lines = []
        app = lines.append
        for i, cycle in enumerate(self):
            app("")
            app("RELAXATION STEP: %d" % (i + 1))
            app(cycle.to_string(verbose=verbose))

        return "\n".join(lines)

    @classmethod
    def from_file(cls, filepath):
        """Initialize the object from the Abinit main output file."""
        with open(filepath, "rt") as stream:
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

    @lazy_property
    def history(self):
        """
        Ordered Dictionary of lists with the evolution of
        the data as function of the relaxation step.
        """
        history = OrderedDict()
        for cycle in self:
            d = cycle.last_iteration
            for k, v in d.items():
                if k in history:
                    history[k].append(v)
                else:
                    history[k] = [v]

        # Convert to numpy arrays.
        for k, v in history.items():
            history[k] = np.array(v)

        return history

    def slideshow(self, **kwargs):
        """
        Uses matplotlib to plot the evolution of the structural relaxation.

        Args:
            ax_list: List of axes. If None a new figure is produced.

        Returns:
            `matplotlib` figure
        """
        for i, cycle in enumerate(self.cycles):
            cycle.plot(title="Relaxation step %s" % (i + 1),
                       tight_layout=kwargs.pop("tight_layout", True),
                       show=kwargs.pop("show", True))

    @add_fig_kwargs
    def plot(self, ax_list=None, fontsize=12, **kwargs):
        """
        Plot relaxation history i.e. the results of the last iteration of each SCF cycle.

        Args:
            ax_list: List of axes. If None a new figure is produced.
            fontsize: legend fontsize.
            kwargs: keyword arguments are passed to ax.plot

        Returns: matplotlib figure
        """
        history = self.history

        # Build grid of plots.
        num_plots, ncols, nrows = len(history), 1, 1
        if num_plots > 1:
            ncols = 2
            nrows = num_plots // ncols + num_plots % ncols

        ax_list, fig, plot = get_axarray_fig_plt(ax_list, nrows=nrows, ncols=ncols,
                                                 sharex=True, sharey=False, squeeze=False)
        ax_list = np.array(ax_list).ravel()

        iter_num = np.array(list(range(self.num_iterations))) + 1
        label = kwargs.pop("label", None)

        for i, ((key, values), ax) in enumerate(zip(history.items(), ax_list)):
            ax.grid(True)
            ax.set_xlabel('Relaxation Step')
            ax.set_xticks(iter_num, minor=False)
            ax.set_ylabel(key)

            xx, yy = iter_num, values
            if not kwargs and label is None:
                ax.plot(xx, yy, "-o", lw=2.0)
            else:
                ax.plot(xx, yy, label=label if i == 0 else None, **kwargs)

            if key in _VARS_SUPPORTING_LOGSCALE and np.all(yy > 1e-22):
                ax.set_yscale("log")

            if key in _VARS_WITH_YRANGE:
                ymin, ymax = _VARS_WITH_YRANGE[key]
                val_min, val_max = np.min(yy), np.max(yy)
                if abs(val_max - val_min) > abs(ymax - ymin):
                    ax.set_ylim(ymin, ymax)

            if label is not None:
                ax.legend(loc="best", fontsize=fontsize, shadow=True)

        # Get around a bug in matplotlib.
        if num_plots % ncols != 0:
            ax_list[-1].plot(xx, yy, lw=0.0)
            ax_list[-1].axis('off')

        return fig

# TODO
#class HaydockIterations(Iterable):
#    """This object collects info on the different steps of the Haydock technique used in the Bethe-Salpeter code"""
#    @classmethod
#    def from_file(cls, filepath):
#        """Initialize the object from file."""
#        with open(filepath, "rt") as stream:
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


class YamlTokenizer(Iterator):
    """
    Provides context-manager support so you can use it in a with statement.
    """
    Error = YamlTokenizerError

    def __init__(self, filename):
        # The position inside the file.
        self.linepos = 0
        self.filename = filename

        try:
            self.stream = open(filename, "rt")
        except IOError as exc:
            # Look for associated error file.
            root, ext = os.path.splitext(self.filename)
            errfile = root + ".err"
            if os.path.exists(errfile) and errfile != self.filename:
                print("Found error file: %s" % errfile)
                with open(errfile, "rt") as fh:
                    print(fh.read())
            raise exc

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
            print("Python traceback:")
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

        raise StopIteration("Cannot find next YAML document in %s" % self.filename)

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
                doc = next(self)
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
        d = yaml.safe_load(doc.text_notag)

        return np.array(d["reduced_coordinates_of_qpoints"])


def yaml_read_irred_perts(filename, doc_tag="!IrredPerts"):
    """Read the list of irreducible perturbations from file."""
    with YamlTokenizer(filename) as r:
        doc = r.next_doc_with_tag(doc_tag)
        d = yaml.safe_load(doc.text_notag)

        return [AttrDict(**pert) for pert in d["irred_perts"]]
        #return d["irred_perts"]


class YamlDoc:
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
        return yaml.safe_load(self.text_notag)
