---
layout: default
title: pymatgen.io.abinit.abitimer.md
nav_exclude: true
---

# pymatgen.io.abinit.abitimer module

This module provides objects for extracting timing data from the ABINIT output files
It also provides tools to analyze and to visualize the parallel efficiency.


### _class_ pymatgen.io.abinit.abitimer.AbinitTimer(sections, info, cpu_time, wall_time)
Bases: `object`

Container class storing the timing results.


* **Parameters**


    * **sections** – List of sections


    * **info** – Dictionary with extra info.


    * **cpu_time** – Cpu-time in seconds.


    * **wall_time** – Wall-time in seconds.



#### cpuwall_histogram(ax=None, \*\*kwargs)
Plot histogram with cpu- and wall-time on axis ax.


* **Parameters**

    **ax** – matplotlib `Axes` or None if a new figure should be created.


Returns: matplotlib figure

Keyword arguments controlling the display of the figure:

| kwargs

 | Meaning

 |
| ------------ | ------------------------------------------------------------------------------------------------- |  |  |  |  |  |  |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### get_dataframe(sort_key='wall_time', \*\*kwargs)
Return a pandas DataFrame with entries sorted according to sort_key.


#### get_section(section_name)
Return section associated to section_name.


#### get_values(keys)
Return a list of values associated to a particular list of keys.


#### names_and_values(key, minval=None, minfract=None, sorted=True)
Select the entries whose value[key] is >= minval or whose fraction[key] is >= minfract
Return the names of the sections and the corresponding values.


#### _property_ ncpus()
Total number of CPUs employed.


#### order_sections(key, reverse=True)
Sort sections according to the value of key.


#### pie(key='wall_time', minfract=0.05, ax=None, \*\*kwargs)
Plot pie chart for this timer.


* **Parameters**


    * **key** – Keyword used to extract data from the timer.


    * **minfract** – Don’t show sections whose relative weight is less that minfract.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


Returns: matplotlib figure

Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### scatter_hist(ax=None, \*\*kwargs)
Scatter plot + histogram.


* **Parameters**

    **ax** – matplotlib `Axes` or None if a new figure should be created.


Returns: matplotlib figure

Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### sum_sections(keys)
Sum value of keys.


#### to_csv(fileobj=<_io.TextIOWrapper name='<stdout>' mode='w' encoding='utf-8'>)
Write data on file fileobj using CSV format.


#### to_table(sort_key='wall_time', stop=None)
Return a table (list of lists) with timer data.


#### totable(sort_key='wall_time', stop=None)
Return a table (list of lists) with timer data.


### _exception_ pymatgen.io.abinit.abitimer.AbinitTimerParseError()
Bases: [`ParseError`](pymatgen.io.core.md#pymatgen.io.core.ParseError)

Errors raised by AbinitTimerParser.


### _class_ pymatgen.io.abinit.abitimer.AbinitTimerParser()
Bases: `Iterable`

Responsible for parsing a list of output files, extracting the timing results
and analyzing the results.
Assume the Abinit output files have been produced with timopt -1.

### Example

parser = AbinitTimerParser()
parser.parse(list_of_files)

To analyze all

```
*
```

.abo files within top, use:

> parser, paths, okfiles = AbinitTimerParser.walk(top=”.”, ext=”.abo”)

Initialize object.


#### BEGIN_TAG(_ = '-<BEGIN_TIMER_ )

#### END_TAG(_ = '-<END_TIMER>_ )

#### Error()
alias of `AbinitTimerParseError`


#### _property_ filenames()
List of files that have been parsed successfully.


#### get_sections(section_name)
Return the list of sections stored in self.timers() given section_name
A fake section is returned if the timer does not have section_name.


#### parse(filenames)
Read and parse a filename or a list of filenames.
Files that cannot be opened are ignored. A single filename may also be given.

Return: list of successfully read files.


#### pefficiency()
Analyze the parallel efficiency.

Return: `ParallelEfficiency` object.


#### plot_all(show=True, \*\*kwargs)
Call all plot methods provided by the parser.


#### plot_efficiency(key='wall_time', what='good+bad', nmax=5, ax=None, \*\*kwargs)
Plot the parallel efficiency.


* **Parameters**


    * **key** – Parallel efficiency is computed using the wall_time.


    * **what** – Specifies what to plot: good for sections with good parallel efficiency.
    bad for sections with bad efficiency. Options can be concatenated with +.


    * **nmax** – Maximum number of entries in plot


    * **ax** – matplotlib `Axes` or None if a new figure should be created.


| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| linewidth

    | matplotlib linewidth. Default: 2.0

                                                                |
| markersize

   | matplotlib markersize. Default: 10

                                                                |

* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_pie(key='wall_time', minfract=0.05, \*\*kwargs)
Plot pie charts of the different timers.


* **Parameters**


    * **key** – Keyword used to extract data from timers.


    * **minfract** – Don’t show sections whose relative weight is less that minfract.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### plot_stacked_hist(key='wall_time', nmax=5, ax=None, \*\*kwargs)
Plot stacked histogram of the different timers.


* **Parameters**


    * **key** – Keyword used to extract data from the timers. Only the first nmax
    sections with largest value are show.


    * **mmax** – Maximum number of sections to show. Other entries are grouped together
    in the others section.


    * **ax** – matplotlib `Axes` or None if a new figure should be created.



* **Returns**

    matplotlib figure


Keyword arguments controlling the display of the figure:

| kwargs

       | Meaning

                                                                                           |
| ------------ | ------------------------------------------------------------------------------------------------- |
| title

        | Title of the plot (Default: None).

                                                                |
| show

         | True to show the figure (default: True).

                                                          |
| savefig

      | “abc.png” or “abc.eps” to save the figure to a file.

                                              |
| size_kwargs

  | Dictionary with options passed to fig.set_size_inches
e.g. size_kwargs=dict(w=3, h=4)

             |
| tight_layout

 | True to call fig.tight_layout (default: False)

                                                    |
| ax_grid

      | True (False) to add (remove) grid from all axes in fig.
Default: None i.e. fig is left unchanged.

 |
| ax_annotate

  | Add labels to  subplots e.g. (a), (b).
Default: False

                                             |
| fig_close

    | Close figure. Default: False.

                                                                     |

#### section_names(ordkey='wall_time')
Return the names of sections ordered by ordkey.
For the time being, the values are taken from the first timer.


#### summarize(\*\*kwargs)
Return pandas DataFrame with the most important results stored in the timers.


#### timers(filename=None, mpi_rank='0')
Return the list of timers associated to the given filename and MPI rank mpi_rank.


#### _classmethod_ walk(top='.', ext='.abo')
Scan directory tree starting from top, look for files with extension ext and
parse timing data.

Return: (parser, paths, okfiles)

    where parser is the new object, paths is the list of files found and okfiles
    is the list of files that have been parsed successfully.
    (okfiles == paths) if all files have been parsed.


### _class_ pymatgen.io.abinit.abitimer.AbinitTimerSection(name, cpu_time, cpu_fract, wall_time, wall_fract, ncalls, gflops)
Bases: `object`

Record with the timing results associated to a section of code.


* **Parameters**


    * **name** – Name of the sections.


    * **cpu_time** – CPU time in seconds.


    * **cpu_fract** – Percentage of CPU time.


    * **wall_time** – Wall-time in seconds.


    * **wall_fract** – Percentage of wall-time.


    * **ncalls** – Number of calls


    * **gflops** – Gigaflops.



#### FIELDS(_ = ('name', 'wall_time', 'wall_fract', 'cpu_time', 'cpu_fract', 'ncalls', 'gflops'_ )

#### NUMERIC_FIELDS(_ = ('wall_time', 'wall_fract', 'cpu_time', 'cpu_fract', 'ncalls', 'gflops'_ )

#### STR_FIELDS(_ = ('name',_ )

#### _classmethod_ fake()
Return a fake section. Mainly used to fill missing entries if needed.


#### to_csvline(with_header=False)
Return a string with data in CSV format. Add header if with_header.


#### to_dict()
Convert object to dictionary.


#### to_tuple()
Convert object to tuple.


### _class_ pymatgen.io.abinit.abitimer.ParallelEfficiency(filenames, ref_idx, \*args, \*\*kwargs)
Bases: `dict`

Store results concerning the parallel efficiency of the job.


* **Parameters**


    * **filennames** – List of filenames


    * **ref_idx** – Index of the Reference time (calculation done with the smallest number of cpus).



#### bad_sections(key='wall_time', criterion='mean', nmax=5)
Return first nmax sections with worst value of key key using criterion criterion.


#### good_sections(key='wall_time', criterion='mean', nmax=5)
Return first nmax sections with best value of key key using criterion criterion.


#### totable(stop=None, reverse=True)
Return table (list of lists) with timing results.


* **Parameters**


    * **stop** – Include results up to stop. None for all


    * **reverse** – Put items with highest wall_time in first positions if True.



### pymatgen.io.abinit.abitimer.alternate(\*iterables)
[a[0], b[0], … , a[1], b[1], …, a[n], b[n] …]
>>> alternate([1,4], [2,5], [3,6])
[1, 2, 3, 4, 5, 6].