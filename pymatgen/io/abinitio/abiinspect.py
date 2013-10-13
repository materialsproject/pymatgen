#!/usr/bin/env python
from __future__ import division, print_function

import collections
import numpy as np


def _magic_parser(stream, magic):
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

            if not line:
                break

            tokens = map(float, line.split()[1:])
            assert len(tokens) == len(keys)
            for l, v in zip(fields.values(), tokens):
                l.append(v)

    return fields


class ScfCycle(collections.Mapping):

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

    @classmethod
    def from_file(cls, filepath):
        """Read the first occurrece of ScfCycle from file."""
        with open(filepath, "r") as stream:
            return cls.from_stream(stream)

    @classmethod
    def from_stream(cls, stream):
        """
        Read the first occurrece of ScfCycle from stream.
        Returns None if no ScfCycle entry is found.
        """
        # TODO: The parser is very fragile and should be replaced by YAML.
        #Example:
        #  iter   Etot(hartree)      deltaE(h)  residm     vres2
        #  ETOT  1  -8.8604027880849    -8.860E+00 2.458E-02 3.748E+00

        #  At SCF step    5       vres2   =  3.53E-08 < tolvrs=  1.00E-06 =>converged.
        fields = _magic_parser(stream, magic=cls.MAGIC)

        if fields:
            fields.pop("iter")
            return cls(fields)
        else:
            return None
        
    def plot(self, **kwargs):
        """
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

        if (num_plots % ncols) != 0:
            ax_list[-1].axis('off')

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

            ax.plot(xx, yy, "-o", lw=2.0)

        if show:
            plt.show()

        if savefig is not None:
            fig.savefig(savefig)


class GroundStateScfCycle(ScfCycle):
    #yaml_tag = '!GroundStateScfCycle'
    MAGIC = "iter   Etot(hartree)"


class PhononScfCycle(ScfCycle):
    #yaml_tag = '!PhononScfCycle'
    MAGIC = "iter   2DEtotal(Ha)"


if __name__ == "__main__":
    import sys
    #scf_cycle = GroundStateScfCycle.from_file(sys.argv[1])
    scf_cycle = PhononScfCycle.from_file(sys.argv[1])

    for key, values in scf_cycle.items():
        print(key, values)

    scf_cycle.plot()
