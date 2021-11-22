.. image:: https://github.com/materialsproject/pymatgen/actions/workflows/test.yml/badge.svg
      :alt: CI Status
      :target: https://github.com/materialsproject/pymatgen/actions/workflows/test.yml

.. image:: https://anaconda.org/conda-forge/pymatgen/badges/downloads.svg
      :alt: Conda Downloads

.. image:: https://coveralls.io/repos/github/materialsproject/pymatgen/badge.svg?branch=master
      :alt: Coveralls Coverage Report
      :target: https://coveralls.io/github/materialsproject/pymatgen?branch=master

Pymatgen (Python Materials Genomics) is a robust, open-source Python library
for materials analysis. These are some of the main features:

1. Highly flexible classes for the representation of Element, Site, Molecule,
   Structure objects.
2. Extensive input/output support, including support for `VASP
   <http://cms.mpi.univie.ac.at/vasp>`_, `ABINIT <http://www.abinit.org>`_, CIF,
   Gaussian, XYZ, and many other file formats.
3. Powerful analysis tools, including generation of phase diagrams, Pourbaix
   diagrams, diffusion analyses, reactions, etc.
4. Electronic structure analyses, such as density of states and band structure.
5. Integration with the Materials Project REST API.

Pymatgen is free to use. However, we also welcome your help to improve this
library by making your own contributions.  These contributions can be in the
form of additional tools or modules you develop, or feature requests and bug
reports. The following are resources for pymatgen:

* `Official documentation <http://pymatgen.org>`_
* Offline docs: HTML are in the pymatgen Github repo's `docs` folder. `Dash <http://kapeli.com/dash>`_ or
  `Zeal <http://zealdocs.org/>`_ docs can be searched and downloaded from "User Contributed Docsets".
* Bug reports or feature requests: Please submit a `GitHub Issue <http://github.com/materialsproject/pymatgen/issues>`_.
* Code contributions via `pull requests <https://github.com/materialsproject/pymatgen/pulls>`_ welcome.
* For help with usage that are unrelated to bugs or feature requests, please use the pymatgen `Discourse page
  <https://discuss.matsci.org/c/pymatgen>`_.
* `matgenb <http://matgenb.materialsvirtuallab.org>`_ provides some Jupyter notebooks demonstrating functionality.
* Follow us on `Twitter <http://twitter.com/pymatgen>`_ to get news and tips.

Major Announcement (v2022.0.*)
==============================

A **backwards incompatible** change has been introduced in v2022.0.*. Pymatgen root-level convenience imports have been
removed from in preparation for a change to a more modular, extensible namespace package architecture that will allow
more developers to contribute. If your existing code uses ``from pymatgen import <something>``, you will need to make
modifications. ``MPRester`` should now be imported from ``pymatgen.ext.matproj``. All other convenience objects such as
``Element``, ``Species``, ``Lattice``, ``Structure``, etc. should be imported from ``pymatgen.core``. There are a few simple ways
you can respond to this change:

* To migrate your code to be compatible with v2022.0.* (it will still be compatible with pymatgen<=2022.0.0 since all
  the imports were already available in previous versions), you need to replace all instances of
  ``from pymatgen import MPRester`` with ``from pymatgen.ext.matproj import MPRester``, followed by replacing all instances
  of ``from pymatgen import`` with ``from pymatgen.core import``. These two steps have to be done in that sequence, since
  MPRester and the other core imports exist in different subpackages. The easiest way is to use an IDE such
  as Pycharm to run a **Replace in Files** on the root directory of your code.
* The pymatgen maintainers have also come up with the following terminal commands you can use to perform the migration.
  On a Mac::

    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i "" 's/from pymatgen import/from pymatgen.core import/g'

  On Linux::

    find . -name '*.py' | xargs sed -i 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '*.py' | xargs sed -i 's/from pymatgen import/from pymatgen.core import/g'

  This should resolve most import errors, though you may have to fix a few issues manually, e.g., if your code contains
  something like ``from pymatgen import Element, MPRester``, which will now need to be split into two lines.

Last but not least, one option is to pin to ``pymatgen==2021.*.*``, which is the last version to contain the root-level
convenience imports, if you are not planning to use future new pymatgen functionality. The new breaking change will
become default from year 2022. Backports to 2021.*.* will still occur for critical bug fixes.

Why use pymatgen?
=================

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers, and is the analysis code powering the
   `Materials Project`_. The analysis it produces survives rigorous scrutiny every single day. Bugs tend to be
   found and corrected quickly. Pymatgen also uses Github Actions for continuous integration, which ensures that every
   new code passes a comprehensive suite of unit tests.
2. **It is well documented.** A fairly comprehensive documentation has been written to help you get to grips with it
   quickly.
3. **It is open.** You are free to use and contribute to pymatgen. It also means that pymatgen is continuously being
   improved. We will attribute any code you contribute to any publication you specify. Contributing to pymatgen means
   your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in pymatgen have been optimized by vectorizing in numpy/scipy.
   This means that coordinate manipulations are extremely fast and are in fact comparable to codes written in other
   languages. Pymatgen also comes with a complete system for handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in the well-established Materials Project.
   It is also actively being developed and maintained by the `Materials Virtual Lab`_, the ABINIT group and many
   other research groups.
6. **A growing ecosystem of developers and add-ons**. Pymatgen has contributions from materials scientists all over the
   world. We also now have an architecture to support add-ons that expand pymatgen's functionality even further. Check
   out the `contributing page <http://pymatgen.org/contributing>`_ and `add-ons page <http://pymatgen.org/addons>`_ for
   details and examples.

Getting pymatgen
================

Before installing pymatgen, you may need to first install a few critical dependencies manually. Please refer to the
official `pymatgen page`_ for installation details and requirements, including instructions for the bleeding edge
developmental version. For people who are absolutely new to Python packages, it is highly recommended you do the
installation using conda, which will make things a lot easier, especially on Windows::

    conda install --channel conda-forge pymatgen

In line with the Scientific Python stack, pymatgen will now support a minimum python version off 3.7 from v2021.1.1.

The version at the `Python Package Index (PyPI) <https://pypi.org/project/pymatgen>`_ is always the latest stable
release that is relatively bug-free and can be installed via pip::

    pip install pymatgen

Some extra functionality (e.g., generation of POTCARs) do require additional setup (see the `pymatgen page`_).

Change Log
==========

The latest change log is available `here <http://pymatgen.org/change_log>`_.

Using pymatgen
==============

Please refer to the official `pymatgen page`_ for tutorials and examples.

How to cite pymatgen
====================

If you use pymatgen in your research, please consider citing the following
work:

    Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier,
    Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier, Kristin A.
    Persson, Gerbrand Ceder. *Python Materials Genomics (pymatgen) : A Robust,
    Open-Source Python Library for Materials Analysis.* Computational
    Materials Science, 2013, 68, 314-319. `doi:10.1016/j.commatsci.2012.10.028
    <https://doi.org/10.1016/j.commatsci.2012.10.028>`_

In addition, some of pymatgen's functionality is based on scientific advances
/ principles developed by the computational materials scientists in our team.
Please refer to `pymatgen's documentation <http://pymatgen.org/>`_ on how to
cite them.

License
=======

Pymatgen is released under the MIT License. The terms of the license are as
follows::

    The MIT License (MIT)
    Copyright (c) 2011-2012 MIT & LBNL

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

About the Pymatgen Development Team
===================================

Shyue Ping Ong of the `Materials Virtual Lab`_ started Pymatgen in 2011, and is
still the project lead.

The Pymatgen Development Team is the set of all contributors to the
pymatgen project, including all subprojects.

Our Copyright Policy
====================

Pymatgen uses a shared copyright model. Each contributor maintains copyright
over their contributions to pymatgen. But, it is important to note that these
contributions are typically only changes to the repositories. Thus, the
pymatgen source code, in its entirety is not the copyright of any
single person or institution. Instead, it is the collective copyright of the
entire pymatgen Development Team. If individual contributors want to maintain a
record of what changes/contributions they have specific copyright on, they
should indicate their copyright in the commit message of the change, when
they commit the change to one of the pymatgen repositories.

With this in mind, the following banner should be used in any source code file
to indicate the copyright and license terms::

    # Copyright (c) Pymatgen Development Team.
    # Distributed under the terms of the MIT License.

.. _`pymatgen page` : http://www.pymatgen.org
.. _`Materials Project` : https://www.materialsproject.org
.. _`Materials Virtual Lab`: http://www.materialsvirtuallab.org
