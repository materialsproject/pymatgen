Compatibility
=============

Pymatgen is a tool used for academic research and is actively developed by
a large community of people. As such, releases are frequent, and new features
and capabilities are constantly being added.

However, pymatgen is also used as a library by other tools, and as such breaking
changes such as the removal or renaming of existing functionality, or substantive
changes in the output of existing code, are tried to be kept to a minimum. This is
especially true of all classes contained in the `pymatgen.core` module.

Despite this, it is sometimes necessary to make breaking changes to enable
future development, or because other libraries we depend upon might change
their own requirements. If a breaking change is causing significant issues,
please post on the GitHub Issues page to see if it can be resolved.

Depending on Pymatgen
---------------------

Pymatgen uses `calendar versioning <http://calver.org/>`_ based on a YYYY-MM-DD format.
This has generally worked well since changes to core pymatgen functionality that most
other codes depend on are rare. There have been only two instances in recent memory that
breaking changes have been made.

* v2021.3.4 - Reorganization of pymatgen into namespace packages, which required the removal
  of root-level imports.
* v2019.3.13 - Renaming of `Site.species_and_occu` to `Site.species`.

Where at all possible, the pymatgen maintainers try to allow for a reasonable deprecation
schedule. For example, the `Site.species` change had a deprecation schedule of about 9 months.
However, some changes such as the recent reorganization of pymatgen into namespace packages
cannot be easily done via a deprecation schedule.

As a compromise solution, pymatgen has adopted **temporary** semantic versioning. A v2021.3.5
was released after v2021.3.4 to reverse the changes made, and new versions v2022.0.x were
released that contains the breaking change (removal of root imports). We will continue to release
critical updates to 2021.x.x versions. This allows end users to continue using 2021.x.x versions
without having to deal with the breaking changes. However, it is recommended that users make the
move to be compatible with 2022.0.x before Jan 1 2022, during which pymatgen will only support the
new namespace architecture and the versioning scheme will go back to calendar versioning.

As the developer of a tool that depends on Pymatgen, you can prevent upgrades of the major
Pymatgen version by specifying a version range like `pymatgen>=2021.1,<2022` or, more
succinctly, using the
`compatible release operator <https://www.python.org/dev/peps/pep-0440/#compatible-release>`_
`pymatgen~=2021.1`. This will prevent `pip` (and other package managers) from
pulling in Pymatgen versions with breaking changes that may end up breaking
your tool.

An even more conservative approach is to pin the Pymatgen dependency to a fixed version, for
example `pymatgen==2021.3.3`. While this will always install the the same version of pymatgen,
it can lead to unnecessary dependency conflicts with other tools that depend on (a different
version of) Pymatgen.

Minimum Python Version
----------------------

As a rule of thumb, pymatgen will support whatever versions of Python the latest
version of numpy supports (at the time of writing, this is Python 3.8+). You can
also check what versions of Python are being tested automatically as part of our
continuous integration testing on GitHub. We currently test pymatgen on Mac,
Windows and Linux.

Recent Breaking Changes
-----------------------

v2022.2.1
~~~~~~~~~

Moved defect-specific code under `defects` module.
`#2582 <https://github.com/materialsproject/pymatgen/pull/2582>`_
#. :code:`pymatgen.transformations.defect_transformations`
#. :code:`pymatgen.analysis.structure_matcher.PointDefectComparator`

Removal of deprecated functions:

`#2405 <https://github.com/materialsproject/pymatgen/pull/2405>`_

#. :code:`Plane.orthonormal_vectors_old()`
#. :code:`ElasticTensor.debye_temperature_from_sound_velocities()`
#. :code:`VaspInputSet.all_input()`

`#2397 <https://github.com/materialsproject/pymatgen/pull/2397>`_

#. :code:`get_dimensionality() in pymatgen/analysis/structure_analyzer.py`
#. :code:`ConversionElectrode.as_dict_summary()`
#. :code:`InsertionElectrode.as_dict_summary()`
#. :code:`Lattice.from_lengths_and_angles()`
#. :code:`Lattice.lengths_and_angles()`
#. :code:`Site.species_and_occu()`

v2022.01.08
~~~~~~~~~~~

Dropped Python 3.7 support for compatibility with the latest numpy. `d00945 <https://github.com/materialsproject/pymatgen/commit/d00945491e9b53548ea8a6755a002c2066ad0ac9>`_ `61ec51c <https://github.com/materialsproject/pymatgen/commit/61ec51cc9751d65df0783af3713e2425d733191e>`_

v2022.0.0
~~~~~~~~~

Pymatgen root imports have been removed from v2022.0.0 in preparation for a change to a more modular, extensible
architecture that will allow more developers to contribute.

Specifically, the following "convenience imports" have been removed in favor of
their canonical import::

    from pymatgen import Composition  # now "from pymatgen.core.composition import Composition"
    from pymatgen import Lattice  # now "from pymatgen.core.lattice import Lattice"
    from pymatgen import SymmOp  # now "from pymatgen.core.operations import SymmOp"
    from pymatgen import DummySpecie, DummySpecies, Element, Specie, Species  # now "from pymatgen.core.periodic_table ..."
    from pymatgen import PeriodicSite, Site  # now "from pymatgen.core.sites ..."
    from pymatgen import IMolecule, IStructure, Molecule, Structure  # now "from pymatgen.core.structure ..."
    from pymatgen import ArrayWithUnit, FloatWithUnit, Unit  # now "from pymatgen.core.units ..."
    from pymatgen import Orbital, Spin  # now "from pymatgen.electronic_structure.core ..."
    from pymatgen import MPRester  # now "from pymatgen.ext.matproj ..."

If your existing code uses `from pymatgen import <something>`, you will need to make
modifications.

The easiest way is to use an IDE to run a Search and Replace.
First, replace any `from pymatgen import MPRester` with
`from pymatgen.ext.matproj import MPRester`. Then, replace
`from pymatgen import` with `from pymatgen.core import`. Alternatively, if you
are using a Mac command line, you can try::

    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import/from pymatgen.core import/g'

From a Linux command line, you can try::

    find . -name '*.py' | xargs sed -i 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i 's/from pymatgen import/from pymatgen.core import/g'

This should resolve most import errors and only a few more modifications may
need to be done by hand.

v2021.3.3
~~~~~~~~~

The variable `pymatgen.SETTINGS` has been moved to `pymatgen.settings.SETTINGS`. Since this is
mostly used internally within pymatgen, it is not expected to lead to significant external issues.

v2021.2.8.1
~~~~~~~~~~~

The minimum version of Python was increased from 3.6 to 3.7 following the lead of numpy. However,
at this point there are no exclusively Python 3.7+ features used in pymatgen so pymatgen may still
be able to be installed manually on Python 3.6 systems, although this usage is not supported.

Support for `aconvasp` has been removed since the corresponding tests were failing and this module
was not being maintained.

v2020.10.20
~~~~~~~~~~~

The band structure plotting functionality, `BSPlotter`, has been overhauled to allow plotting of
multiple band structures. This might cause issues for tools relying on the internal structure
of BSPlotter's plot data.
