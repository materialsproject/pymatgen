<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/materialsproject/pymatgen/master/docs/_images/pymatgen-white.svg">
    <img alt="Logo" src="https://raw.githubusercontent.com/materialsproject/pymatgen/master/docs/_images/pymatgen.svg" height="70">
  </picture>
</h1>

<h4 align="center">

[![CI Status](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml/badge.svg)](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml)
[![Coveralls](https://img.shields.io/coveralls/github/materialsproject/pymatgen?logo=coveralls&label=Coverage)](https://coveralls.io/github/materialsproject/pymatgen?branch=master)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pymatgen?logo=pypi&logoColor=white&color=blue&label=PyPI)](https://pypi.org/project/pymatgen)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pymatgen?logo=condaforge&color=blue&label=Conda)](https://anaconda.org/conda-forge/pymatgen)
[![Requires Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)

</h4>

Pymatgen (Python Materials Genomics) is a robust, open-source Python
library for materials analysis. These are some of the main features:

1. Highly flexible classes for the representation of `Element`, `Site`, `Molecule` and `Structure` objects.
2. Extensive input/output support, including support for [VASP](https://cms.mpi.univie.ac.at/vasp), [ABINIT](https://abinit.org), CIF, Gaussian, XYZ, and many other file formats.
3. Powerful analysis tools, including generation of phase diagrams, Pourbaix diagrams, diffusion analyses, reactions, etc.
4. Electronic structure analyses, such as density of states and band structure.
5. Integration with the [Materials Project] REST API.

Pymatgen is free to use. However, we also welcome your help to improve this library by making your contributions. These contributions can be in the form of additional tools or modules you develop, or feature requests and bug reports. The following are resources for `pymatgen`:

- [Official documentation](https://pymatgen.org)
- Bug reports or feature requests: Please submit a [GitHub Issue](https://github.com/materialsproject/pymatgen/issues).
- Code contributions via [pull requests](https://github.com/materialsproject/pymatgen/pulls) are welcome.
- For help with usage that is unrelated to bugs or feature requests, please use the `pymatgen` [MatSci page](https://discuss.matsci.org/c/pymatgen).
- [`matgenb`](https://github.com/materialsvirtuallab/matgenb#introduction) provides some Jupyter notebooks demonstrating functionality.
- Follow us on [Twitter](https://twitter.com/pymatgen) to get news and tips.

## Why use pymatgen?

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers and is the analysis code powering the [Materials Project]. The analysis it produces survives rigorous scrutiny every single day. Bugs tend to be found and corrected quickly. Pymatgen also uses Github Actions for continuous integration, which ensures that every new code passes a comprehensive suite of unit tests.
2. **It is well documented.** A fairly comprehensive documentation has been written to help you get to grips with it quickly.
3. **It is open.** You are free to use and contribute to `pymatgen`. It also means that `pymatgen` is continuously being improved. We will attribute any code you contribute to any publication you specify. Contributing to `pymatgen` means your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in `pymatgen` have been optimized by vectorizing in `numpy`/`scipy`. This means that coordinate manipulations are extremely fast and are in fact comparable to codes written in other languages. Pymatgen also comes with a complete system for handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in the well-established Materials Project. It is also actively being developed and maintained by the [Materials Virtual Lab], the ABINIT group and many other research groups.
6. **A growing ecosystem of developers and add-ons**. Pymatgen has contributions from materials scientists all over the world. We also now have an architecture to support add-ons that expand pymatgen's functionality even further. Check out the [contributing page](https://pymatgen.org/contributing) and [add-ons page](https://pymatgen.org/addons) for details and examples.

## Installation

The version at the [Python Package Index (PyPI)](https://pypi.org/project/pymatgen) is always the latest stable release that is relatively bug-free and can be installed via `pip`:

```sh
pip install pymatgen
```

The minimum Python version is 3.8. Some extra functionality (e.g., generation of POTCARs) does require additional setup (see the [`pymatgen` page]).

## Change Log

Please check [GitHub releases](https://github.com/materialsproject/pymatgen/releases) and [commit history](https://github.com/materialsproject/pymatgen/commits/master) for the latest changes. A legacy changelog is still up at <https://pymatgen.org/change_log>.

## Using pymatgen

Please refer to the official [`pymatgen` page] for tutorials and examples.

## How to cite pymatgen

If you use `pymatgen` in your research, please consider citing the following work:

> Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy
> Hautier, Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier,
> Kristin A. Persson, Gerbrand Ceder. *Python Materials Genomics
> (pymatgen): A Robust, Open-Source Python Library for Materials
> Analysis.* Computational Materials Science, 2013, 68, 314-319.
> [doi:10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)

In addition, some of `pymatgen`'s functionality is based on scientific advances/principles developed by the computational materials scientists in our team. Please refer to [`pymatgen`'s documentation](https://pymatgen.org) on how to cite them.

## License

Pymatgen is released under the MIT License. The terms of the license are as follows:

> The MIT License (MIT) Copyright (c) 2011-2012 MIT & LBNL
>
> Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## About the Pymatgen Development Team

Shyue Ping Ong of the [Materials Virtual Lab] started Pymatgen in 2011 and is still the project lead.

The [`pymatgen` development team] is the set of all contributors to the `pymatgen` project, including all subprojects.

## Our Copyright Policy

Pymatgen uses a shared copyright model. Each contributor maintains copyright over their contributions to `pymatgen`. But, it is important to note that these contributions are typically only changes to the repositories. Thus, the `pymatgen` source code, in its entirety is not the copyright of any single person or institution. Instead, it is the collective copyright of the entire [`pymatgen` Development Team]. If individual contributors want to maintain a record of what changes/contributions they have specific copyright on, they should indicate their copyright in the commit message of the change, when they commit the change to one of the `pymatgen` repositories.

[`pymatgen` page]: https://pymatgen.org
[materials project]: https://materialsproject.org
[`pymatgen` development team]: https://pymatgen.org/team
[materials virtual lab]: https://materialsvirtuallab.org
