# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module provides objects for extracting timing data from the ABINIT output files 
It also provides tools to analye and to visualize the parallel efficiency.
"""
from __future__ import unicode_literals, division

import sys
import collections
import numpy as np

from six.moves import zip
from monty.string import is_string, list_strings
from pymatgen.util.num_utils import minloc
from pymatgen.util.plotting_utils import add_fig_kwargs, get_ax_fig_plt

import logging
logger = logging.getLogger(__name__)


def alternate(*iterables):
    """
    [a[0], b[0], ... , a[1], b[1], ..., a[n], b[n] ...]
    >>> alternate([1,4], [2,5], [3,6])
    [1, 2, 3, 4, 5, 6]
    """
    items = []
    for tup in zip(*iterables):
        items.extend([item for item in tup])
    return items



class AbinitTimerParserError(Exception):
    """Errors raised by AbinitTimerParser"""


class AbinitTimerParser(collections.Iterable):
    """
    Responsible for parsing a list of output files, and managing the parsed database.
    """
    # The markers enclosing the data.
    BEGIN_TAG = "-<BEGIN_TIMER"
    END_TAG = "-<END_TIMER>"

    Error = AbinitTimerParserError

    #DEFAULT_MPI_RANK = "0"

    def __init__(self):
        # List of files that have been parsed.
        self._filenames = []

        # timers[filename][mpi_rank]
        # contains the timer extracted from the file filename 
        # associated to the MPI rank mpi_rank.
        self._timers = collections.OrderedDict()

    def __iter__(self):
        return self._timers.__iter__()

    def __len__(self):
        return len(self._timers)

    def parse(self, filenames):
        """
        Read and parse a filename or a list of filenames.

        Files that cannot be opened are ignored. A single filename may also be given.
        Return list of successfully read files.
        """
        filenames = list_strings(filenames)

        read_ok = []
        for fname in filenames:
            try:
                fh = open(fname)
            except IOError:
                logger.warning("Cannot open file %s" % fname)
                continue

            try:
                self._read(fh, fname)
                read_ok.append(fname)

            except self.Error as e:
                logger.warning("exception while parsing file %s:\n%s" % (fname, str(e)))
                continue

            finally:
                fh.close()

        # Add read_ok to the list of files that have been parsed.
        self._filenames.extend(read_ok)
        return read_ok

    def _read(self, fh, fname):
        """Parse the TIMER section"""
        if fname in self._timers:
            raise self.Error("Cannot overwrite timer associated to: %s " % fname)

        data = {}

        def parse_line(line):
            name, vals = line[:25], line[25:].split()
            ctime, cfract, wtime, wfract, ncalls, gflops = vals
            return AbinitTimerSection(name, ctime, cfract, wtime, wfract, ncalls, gflops)

        inside, has_timer = 0, False
        for line in fh:
            #print(line.strip())
            if line.startswith(self.BEGIN_TAG):
                has_timer = True
                sections = []
                info = {}
                inside = 1
                line = line[len(self.BEGIN_TAG):].strip()[:-1]

                info["fname"] = fname
                for tok in line.split(","):
                    (key, val) = [s.strip() for s in tok.split("=")]
                    info[key] = val

            elif line.startswith(self.END_TAG):
                inside = 0
                timer = AbinitTimer(sections, info, cpu_time, wall_time)
                mpi_rank = info["mpi_rank"]

                data[mpi_rank] = timer

            elif inside:
                inside += 1
                line = line[1:].strip()

                if inside == 2:
                    d = dict()
                    for tok in line.split(","):
                        (key, val) = [s.strip() for s in tok.split("=")]
                        d[key] = float(val)
                    cpu_time, wall_time = d["cpu_time"], d["wall_time"]

                elif inside > 5:
                    sections.append(parse_line(line))

                else:
                    try:
                        parse_line(line)
                    except:
                        parser_failed = True

                    if not parser_failed:
                        raise self.Error("line should be empty: " + str(inside) + line)

        if not has_timer:
            raise self.Error("%s: No timer section found" % fname)

        # Add it to the dict
        self._timers[fname] = data

    #def set_default_mpi_rank(mpi_rank): self._default_mpi_rank = mpi_rank
    #def get_default_mpi_rank(mpi_rank): return self._default_mpi_rank

    def timers(self, filename=None, mpi_rank="0"):
        """Return the list of timers associated to the given filename and MPI rank mpi_rank."""
        if filename is not None:
            timers = [self._timers[filename][mpi_rank]]
        else:
            timers = [self._timers[filename][mpi_rank] for filename in self._filenames]

        return timers

    def section_names(self, ordkey="wall_time"):
        """Return the names of sections ordered by ordkey."""
        section_names = [] # Avoid UnboundLocalError

        # FIXME this is not trivial
        for (idx, timer) in enumerate(self.timers()):
            if idx == 0:
                section_names = [s.name for s in timer.order_sections(ordkey)]
                #check = section_names
                #else:
                #  new_set = set( [s.name for s in timer.order_sections(ordkey)])
                #  section_names.intersection_update(new_set)
                #  check = check.union(new_set)

        #if check != section_names:
        #  print("sections", section_names)
        #  print("check",check)

        return section_names

    def get_sections(self, section_name):
        """
        Return the list of sections stored in self.timers() whose name is section_name
        A fake section is returned if the timer does not have sectio_name.
        """
        sections = []
        for timer in self.timers():
            for sect in timer.sections:
                if sect.name == section_name:
                    sections.append(sect)
                    break
            else:
                sections.append(AbinitTimerSection.fake())

        return sections

    def pefficiency(self):
        """
        Analyze the parallel efficiency.
        """
        timers = self.timers()
        #
        # Number of CPUs employed in each calculation.
        ncpus = [timer.ncpus for timer in timers]

        # Find the minimum number of cpus used and its index in timers.
        min_idx = minloc(ncpus)
        min_ncpus = ncpus[min_idx]

        # Reference timer
        ref_t = timers[min_idx]

        # Compute the parallel efficiency (total efficieny and the efficiency of each section)
        peff = {}
        ctime_peff = [(min_ncpus * ref_t.wall_time) / (t.wall_time * ncp) for (t, ncp) in zip(timers, ncpus)]
        wtime_peff = [(min_ncpus * ref_t.cpu_time) / (t.cpu_time * ncp) for (t, ncp) in zip(timers, ncpus)]
        n = len(timers)

        peff["total"] = {}
        peff["total"]["cpu_time"] = ctime_peff
        peff["total"]["wall_time"] = wtime_peff
        peff["total"]["cpu_fract"] = n * [100]
        peff["total"]["wall_fract"] = n * [100]

        for sect_name in self.section_names():
            #print(sect_name)
            ref_sect = ref_t.get_section(sect_name)
            sects = [t.get_section(sect_name) for t in timers]
            try:
                ctime_peff = [(min_ncpus * ref_sect.cpu_time) / (s.cpu_time * ncp) for (s, ncp) in zip(sects, ncpus)]
                wtime_peff = [(min_ncpus * ref_sect.wall_time) / (s.wall_time * ncp) for (s, ncp) in zip(sects, ncpus)]
            except ZeroDivisionError:
                ctime_peff = n * [-1]
                wtime_peff = n * [-1]

            assert sect_name not in peff
            peff[sect_name] = {}
            peff[sect_name]["cpu_time"] = ctime_peff
            peff[sect_name]["wall_time"] = wtime_peff

            peff[sect_name]["cpu_fract"] = [s.cpu_fract for s in sects]
            peff[sect_name]["wall_fract"] = [s.wall_fract for s in sects]

        return ParallelEfficiency(self._filenames, min_idx, peff)

    @add_fig_kwargs
    def plot_efficiency(self, key="wall_time", what="gb", nmax=5, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        timers = self.timers()
        peff = self.pefficiency()

        # Table with the parallel efficiency for all the sections.
        #pprint_table(peff.totable())

        n = len(timers)
        xx = np.arange(n)

        ax.set_color_cycle(['g', 'b', 'c', 'm', 'y', 'k'])

        legend_entries = []

        # Plot sections with good efficiency.
        lines = []
        if "g" in what:
            good = peff.good_sections(key=key, nmax=nmax)
            for g in good:
                #print(g, peff[g])
                yy = peff[g][key]
                line, = ax.plot(xx, yy, "-->", linewidth=3.0, markersize=10)
                lines.append(line)
                legend_entries.append(g)

        # Plot sections with bad efficiency.
        if "b" in what:
            bad = peff.bad_sections(key=key, nmax=nmax)
            for b in bad:
                #print(b, peff[b])
                yy = peff[b][key]
                line, = ax.plot(xx, yy, "-.<", linewidth=3.0, markersize=10)
                lines.append(line)
                legend_entries.append(b)

        if "total" not in legend_entries:
            yy = peff["total"][key]
            total_line, = ax.plot(xx, yy, "r", linewidth=3.0, markersize=10)
            lines.append(total_line)
            legend_entries.append("total")

        ax.legend(lines, legend_entries, loc="best", shadow=True)

        #ax.set_title(title)
        ax.set_xlabel('Total_NCPUs')
        ax.set_ylabel('Efficiency')
        ax.grid(True)

        # Set xticks and labels.
        labels = ["MPI = %d, OMP = %d" % (t.mpi_nprocs, t.omp_nthreads) for t in timers]
        ax.set_xticks(xx)
        ax.set_xticklabels(labels, fontdict=None, minor=False, rotation=15)

        return fig

    @add_fig_kwargs
    def plot_pie(self, key="wall_time", minfract=0.05, ax=None, **kwargs):
        """Pie charts of the different timers."""
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        timers = self.timers()
        n = len(timers)

        # Make square figures and axes
        the_grid = plt.GridSpec(n, 1)

        fig = plt.figure(1, figsize=(6, 6))

        for idx, timer in enumerate(timers):
            plt.subplot(the_grid[idx, 0])
            plt.title(str(timer))
            timer.pie(key=key, minfract=minfract)

        return fig

    @add_fig_kwargs
    def plot_stacked_hist(self, key="wall_time", nmax=5, ax=None, **kwargs):
        """Stacked histogram of the different timers."""
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        mpi_rank = "0"
        timers = self.timers(mpi_rank=mpi_rank)
        n = len(timers)

        names, values = [], []
        rest = np.zeros(n)

        for idx, sname in enumerate(self.section_names(ordkey=key)):
            sections = self.get_sections(sname)
            svals = np.asarray([s.__dict__[key] for s in sections])

            if idx < nmax:
                names.append(sname)
                values.append(svals)
            else:
                rest += svals

        names.append("others (nmax = %d)" % nmax)
        values.append(rest)
        #for (n, vals) in zip(names, values): print(n, vals)

        # The dataset is stored in values.
        # Now create the stacked histogram.

        ind = np.arange(n) # the locations for the groups
        width = 0.35       # the width of the bars

        # this does not work with matplotlib < 1.0
        #plt.rcParams['axes.color_cycle'] = ['r', 'g', 'b', 'c']
        colors = nmax * ['r', 'g', 'b', 'c', 'k', 'y', 'm']

        bars = []
        bottom = np.zeros(n)

        for idx, vals in enumerate(values):
            color = colors[idx]

            bar = plt.bar(ind, vals, width, color=color, bottom=bottom)
            bars.append(bar)

            bottom += vals

        ax.set_ylabel(key)
        #ax.title("Stacked histogram for the %d most important sections" % nmax)

        labels = ["MPI = %d, OMP = %d" % (t.mpi_nprocs, t.omp_nthreads) for t in timers]
        plt.xticks(ind + width / 2.0, labels, rotation=15)
        #plt.yticks(np.arange(0,81,10))

        ax.legend([bar[0] for bar in bars], names, loc="best")

        return fig


class ParallelEfficiency(dict):

    def __init__(self, filenames, ref_idx, *args, **kwargs):
        self.update(*args, **kwargs)
        self.filenames = filenames
        self._ref_idx = ref_idx

    def _order_by_peff(self, key, criterion, reverse=True):

        estimators = {
            "min": min,
            "max": max,
            "mean": lambda items: sum(items) / len(items)
        }

        self.estimator = estimators[criterion]

        data = []
        for (sect_name, peff) in self.items():
            #
            # Ignore values where we had a division by zero.
            if all([v != -1 for v in peff[key]]):
                values = peff[key][:]
                #print(sect_name, values)
                if len(values) > 1:
                    ref_value = values.pop(self._ref_idx)
                    assert ref_value == 1.0

                data.append((sect_name, self.estimator(values)))

        fsort = lambda t: t[1]
        data.sort(key=fsort, reverse=reverse)
        return tuple([sect_name for (sect_name, e) in data])

    def totable(self, stop=None, reverse=True):
        osects = self._order_by_peff("wall_time", criterion="mean", reverse=reverse)

        n = len(self.filenames)
        table = [["AbinitTimerSection"] + alternate(self.filenames, n * ["%"])]
        for sect_name in osects:
            peff = self[sect_name]["wall_time"]
            fract = self[sect_name]["wall_fract"]
            vals = alternate(peff, fract)

            table.append([sect_name] + ["%.2f" % val for val in vals])

        return table

    def good_sections(self, key="wall_time", criterion="mean", nmax=5):
        good_sections = self._order_by_peff(key, criterion=criterion)
        return good_sections[:nmax]

    def bad_sections(self, key="wall_time", criterion="mean", nmax=5):
        bad_sections = self._order_by_peff(key, criterion=criterion, reverse=False)
        return bad_sections[:nmax]


class AbinitTimerSection(object):
    """Record with the timing results associated to a section of code."""
    STR_FIELDS = [
        "name"
    ]

    NUMERIC_FIELDS = [
        "cpu_time",
        "cpu_fract",
        "wall_time",
        "wall_fract",
        "ncalls",
        "gflops",
    ]

    FIELDS = tuple(STR_FIELDS + NUMERIC_FIELDS)

    @classmethod
    def fake(cls):
        return AbinitTimerSection("fake", 0.0, 0.0, 0.0, 0.0, -1, 0.0)

    def __init__(self, name, cpu_time, cpu_fract, wall_time, wall_fract, ncalls, gflops):
        self.name = name.strip()
        self.cpu_time = float(cpu_time)
        self.cpu_fract = float(cpu_fract)
        self.wall_time = float(wall_time)
        self.wall_fract = float(wall_fract)
        self.ncalls = int(ncalls)
        self.gflops = float(gflops)

    def to_tuple(self):
        return tuple([self.__dict__[at] for at in AbinitTimerSection.FIELDS])

    def to_csvline(self, with_header=False):
        """Return a string with data in CSV format"""
        string = ""

        if with_header:
            string += "# " + " ".join(at for at in AbinitTimerSection.FIELDS) + "\n"

        string += ", ".join(str(v) for v in self.to_tuple()) + "\n"
        return string

    def __str__(self):
        string = ""
        for a in AbinitTimerSection.FIELDS: string += a + " = " + self.__dict__[a] + ","
        return string[:-1]


class AbinitTimer(object):
    """Container class used to store the timing results."""

    def __init__(self, sections, info, cpu_time, wall_time):

        self.sections = tuple(sections)
        self.section_names = tuple([s.name for s in self.sections])
        self.info = info
        self.cpu_time = float(cpu_time)
        self.wall_time = float(wall_time)

        self.mpi_nprocs = int(info["mpi_nprocs"])
        self.omp_nthreads = int(info["omp_nthreads"])
        self.mpi_rank = info["mpi_rank"].strip()
        self.fname = info["fname"].strip()

    def __str__(self):
        string = "file = %s, wall_time = %.1f, mpi_nprocs = %d, omp_nthreads = %d" % (
            self.fname, self.wall_time, self.mpi_nprocs, self.omp_nthreads )
        #string += ", rank = " + self.mpi_rank
        return string

    def __cmp__(self, other):
        return cmp(self.wall_time, other.wall_time)

    @property
    def ncpus(self):
        """Total number of CPUs employed."""
        return self.mpi_nprocs * self.omp_nthreads

    def get_section(self, section_name):
        try:
            idx = self.section_names.index(section_name)
        except:
            raise
        sect = self.sections[idx]
        assert sect.name == section_name
        return sect

    def to_csv(self, fileobj=sys.stdout):
        """Write data on file fileobj using CSV format."""
        openclose = is_string(fileobj)

        if openclose:
            fileobj = open(fileobj, "w")

        for (idx, section) in enumerate(self.sections):
            fileobj.write(section.to_csvline(with_header=(idx == 0)))
        fileobj.flush()

        if openclose:
            fileobj.close()

    def totable(self, sort_key="wall_time", stop=None):
        """Return a table (list of lists) with timer data"""
        table = [list(AbinitTimerSection.FIELDS), ]
        ord_sections = self.order_sections(sort_key)

        if stop is not None:
            ord_sections = ord_sections[:stop]

        for osect in ord_sections:
            row = [str(item) for item in osect.to_tuple()]
            table.append(row)

        return table

    def get_values(self, keys):
        """Return a list of values associated to a particular list of keys"""
        if is_string(keys):
            return [s.__dict__[keys] for s in self.sections]
        else:
            values = []
            for k in keys:
                values.append([s.__dict__[k] for s in self.sections])
            return values

    def names_and_values(self, key, minval=None, minfract=None, sorted=True):
        """
        Select the entries whose value[key] is >= minval or whose fraction[key] is >= minfract
        Return the names of the sections and the correspoding value
        """
        values = self.get_values(key)
        names = self.get_values("name")

        new_names, new_values = [], []
        other_val = 0.0

        if minval is not None:
            assert minfract is None

            for (n, v) in zip(names, values):
                if v >= minval:
                    new_names.append(n)
                    new_values.append(v)
                else:
                    other_val += v

            new_names.append("below minval " + str(minval))
            new_values.append(other_val)

        elif minfract is not None:
            assert minval is None

            total = self.sum_sections(key)

            for (n, v) in zip(names, values):
                if v / total >= minfract:
                    new_names.append(n)
                    new_values.append(v)
                else:
                    other_val += v

            new_names.append("below minfract " + str(minfract))
            new_values.append(other_val)

        else:
            # all values
            (new_names, new_values) = (names, values)

        if sorted:
            # Sort new_values and rearrange new_names.
            fsort = lambda t: t[1]
            nandv = [nv for nv in zip(new_names, new_values)]
            nandv.sort(key=fsort)
            new_names, new_values = [n[0] for n in nandv], [n[1] for n in nandv]

        return new_names, new_values

    def _reduce_sections(self, keys, operator):
        return operator(self.get_values(keys))

    def sum_sections(self, keys):
        return self._reduce_sections(keys, sum)

    def order_sections(self, key, reverse=True):
        """Sort sections according to the value of key."""
        fsort = lambda s: s.__dict__[key]
        return sorted(self.sections, key=fsort, reverse=reverse)

    @add_fig_kwargs
    def cpuwall_histogram(self, ax=None, **kwargs):
        ax, fig, plt = get_ax_fig_plt(ax=ax)

        nk = len(self.sections)
        ind = np.arange(nk)  # the x locations for the groups
        width = 0.35         # the width of the bars

        cpu_times = self.get_values("cpu_time")
        rects1 = plt.bar(ind, cpu_times, width, color='r')

        wall_times = self.get_values("wall_time")
        rects2 = plt.bar(ind + width, wall_times, width, color='y')

        # Add ylable and title
        ax.set_ylabel('Time (s)')

        #if title:
        #    plt.title(title)
        #else:
        #    plt.title('CPU-time and Wall-time for the different sections of the code')

        ticks = self.get_values("name")
        ax.set_xticks(ind + width, ticks)

        ax.legend((rects1[0], rects2[0]), ('CPU', 'Wall'), loc="best")

        return fig

    #def hist2(self, key1="wall_time", key2="cpu_time"):

    #    labels = self.get_values("name")
    #    vals1, vals2 = self.get_values([key1, key2])

    #    N = len(vals1)
    #    assert N == len(vals2)

    #    plt.figure(1)
    #    plt.subplot(2, 1, 1) # 2 rows, 1 column, figure 1

    #    n1, bins1, patches1 = plt.hist(vals1, N, facecolor="m")
    #    plt.xlabel(labels)
    #    plt.ylabel(key1)

    #    plt.subplot(2, 1, 2)
    #    n2, bins2, patches2 = plt.hist(vals2, N, facecolor="y")
    #    plt.xlabel(labels)
    #    plt.ylabel(key2)

    #    plt.show()

    def pie(self, key="wall_time", minfract=0.05, title=None):
        import matplotlib.pyplot as plt

        # Don't show section whose value is less that minfract
        labels, vals = self.names_and_values(key, minfract=minfract)

        return plt.pie(vals, explode=None, labels=labels, autopct='%1.1f%%', shadow=True)

    def scatter_hist(self, ax=None, **kwargs):
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        ax, fig, plt = get_ax_fig_plt(ax=ax)

        #title = kwargs.pop("title", None)
        #show = kwargs.pop("show", True)
        #savefig = kwargs.pop("savefig", None)
        #fig = plt.figure(1, figsize=(5.5, 5.5))

        x = np.asarray(self.get_values("cpu_time"))
        y = np.asarray(self.get_values("wall_time"))

        # the scatter plot:
        axScatter = plt.subplot(1, 1, 1)
        axScatter.scatter(x, y)
        axScatter.set_aspect("auto")

        # create new axes on the right and on the top of the current axes
        # The first argument of the new_vertical(new_horizontal) method is
        # the height (width) of the axes to be created in inches.
        divider = make_axes_locatable(axScatter)
        axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
        axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=axScatter)

        # make some labels invisible
        plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)

        # now determine nice limits by hand:
        binwidth = 0.25
        xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
        lim = (int(xymax / binwidth) + 1) * binwidth

        bins = np.arange(-lim, lim + binwidth, binwidth)
        axHistx.hist(x, bins=bins)
        axHisty.hist(y, bins=bins, orientation='horizontal')

        # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
        # thus there is no need to manually adjust the xlim and ylim of these axis.

        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)
        for tl in axHistx.get_xticklabels():
            tl.set_visible(False)
            axHistx.set_yticks([0, 50, 100])

            #axHisty.axis["left"].major_ticklabels.set_visible(False)
            for tl in axHisty.get_yticklabels():
                tl.set_visible(False)
                axHisty.set_xticks([0, 50, 100])

        plt.draw()

        return fig
