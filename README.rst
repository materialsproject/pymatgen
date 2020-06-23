**Official docs:** http://www.pymatgen.org

.. image:: https://circleci.com/gh/materialsproject/pymatgen.svg?style=shield&circle-token=:circle-token
        :alt: CircleCI Status
        :target: https://circleci.com/gh/materialsproject/pymatgen

.. image:: https://ci.appveyor.com/api/projects/status/akdyke5jxg6gps45?svg=true
        :alt: AppVeyor Status
        :target: https://ci.appveyor.com/api/projects/status/akdyke5jxg6gps45?svg=true

.. image:: https://anaconda.org/conda-forge/pymatgen/badges/downloads.svg
        :alt: Conda Downloads
        :target: https://anaconda.org/conda-forge/pymatgen/badges/downloads.svg

.. image:: https://coveralls.io/repos/github/materialsproject/pymatgen/badge.svg?branch=master
        :alt: Coveralls Coverage Report
        :target: https://coveralls.io/repos/github/materialsproject/pymatgen/badge.svg?branch=master

.. image:: https://img.shields.io/github/commit-activity/y/materialsproject/pymatgen
        :target: https://img.shields.io/github/commit-activity/y/materialsproject/pymatgen
        :alt: GitHub commit activity

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

* Please report any bugs and issues at pymatgen's `Github Issues
  page <https://github.com/materialsproject/pymatgen/issues>`_.
* For help with any pymatgen issue, please use the pymatgen `Discourse page
  <https://discuss.matsci.org/c/pymatgen>`_. Please note that the pymatgen Google
  group has been deprecated in favor of Discourse.
* `Twitter <http://twitter.com/pymatgen>`_. Follow to get news and tips.
* `matgenb <http://matgenb.materialsvirtuallab.org>`_. For example notebooks.

Why use pymatgen?
=================

There are many materials analysis codes out there, both commerical and free,
but pymatgen offer several advantages:

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers,
   and is the analysis code powering the `Materials Project`_. The analysis it
   produces survives rigorous scrutiny every single day. Bugs tend to be
   found and corrected quickly. Pymatgen also uses
   `CircleCI <https://circleci.com>`_ and `Appveyor <https://www.appveyor.com/>`_
   for continuous integration on the Linux and Windows platforms,
   respectively, which ensures that every commit passes a comprehensive suite
   of unittests. The coverage of the unittests can be seen on
   `coveralls.io <https://coveralls.io/github/materialsproject/pymatgen>`_.
2. **It is well documented.** A fairly comprehensive documentation has been
   written to help you get to grips with it quickly.
3. **It is open.** You are free to use and contribute to pymatgen. It also means
   that pymatgen is continuously being improved. We will attribute any code you
   contribute to any publication you specify. Contributing to pymatgen means
   your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in pymatgen have been
   optimized by vectorizing in numpy/scipy. This means that coordinate
   manipulations are extremely fast and are in fact comparable to codes
   written in other languages. Pymatgen also comes with a complete system for
   handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in
   the well-established Materials Project. It is also actively being developed
   and maintained by the `Materials Virtual Lab`_, the ABINIT group and many
   other research groups.

Getting pymatgen
================

Before installing pymatgen, you may need to first install a few critical
dependencies manually. Please refer to the official `pymatgen page`_ for
installation details and requirements, including instructions for the
bleeding edge developmental version. For people who are absolutely new to
Python packages, it is highly recommended you do the installation using
conda, which will make things a lot easier, especially on Windows:

    conda install --channel conda-forge pymatgen

In line with the Scientific Python stack, pymatgen will now support only
Py3.x from v2019.1.1. Specifically, we now only run testing on Py3.6+ so
this is our officially  supported minimum Python version.

Users who need Python 2.7 support should install v2018.x,
you may also need to enforce an older version of numpy (`pip install numpy==1.16.4 `).

The version at the `Python Package Index (PyPI) <https://pypi.org/project/pymatgen>`_
is always the latest stable release that is relatively bug-free. The easiest
way to install pymatgen on any system is via pip::

    pip install pymatgen

Wheels for Mac and Windows have been built for convenience.

Some extra functionality (e.g., generation of POTCARs) do require additional
setup (please see the `pymatgen page`_).

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
    <http://dx.doi.org/10.1016/j.commatsci.2012.10.028>`_

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
