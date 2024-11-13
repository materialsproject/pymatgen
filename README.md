<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/materialsproject/pymatgen/raw/master/docs/assets/pymatgen-white.svg">
    <img alt="Logo" src="https://github.com/materialsproject/pymatgen/raw/master/docs/assets/pymatgen.svg"
height="70">
  </picture>
</h1>

<h4 align="center">

[![CI Status](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml/badge.svg)](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/materialsproject/pymatgen/branch/master/graph/badge.svg?token=XC47Un1LV2)](https://codecov.io/gh/materialsproject/pymatgen)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pymatgen?logo=pypi&logoColor=white&color=blue&label=PyPI)](https://pypi.org/project/pymatgen)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pymatgen?logo=condaforge&color=blue&label=Conda)](https://anaconda.org/conda-forge/pymatgen)
[![Requires Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg?logo=python&logoColor=white)](https://python.org/downloads)
[![Paper](https://img.shields.io/badge/J.ComMatSci-2012.10.028-blue?logo=elsevier&logoColor=white)](https://doi.org/10.1016/j.commatsci.2012.10.028)

</h4>

Pymatgen (Python Materials Genomics) is a robust, open-source Python
library for materials analysis. These are some of the main features:

1. Highly flexible classes for the representation of `Element`, `Site`, `Molecule` and `Structure` objects.
2. Extensive input/output support, including support for [VASP](https://cms.mpi.univie.ac.at/vasp), [ABINIT](https://abinit.org), [CIF](https://wikipedia.org/wiki/Crystallographic_Information_File), [Gaussian](https://gaussian.com), [XYZ](https://wikipedia.org/wiki/XYZ_file_format), and many other file formats.
3. Powerful analysis tools, including generation of phase diagrams, Pourbaix diagrams, diffusion analyses, reactions, etc.
4. Electronic structure analyses, such as density of states and band structure.
5. Integration with the [Materials Project] REST API.

Pymatgen is free to use. However, we also welcome your help to improve this library by making your contributions. These contributions can be in the form of additional tools or modules you develop, or feature requests and bug reports. The following are resources for `pymatgen`:

- [Official documentation][`pymatgen` docs]
- Bug reports or feature requests: Please submit a [GitHub issue].
- Code contributions via [pull request] are welcome.
- For questions that are not bugs or feature requests, please use the `pymatgen` [MatSci forum](https://matsci.org/pymatgen) or open a [GitHub discussion].
- [`matgenb`](https://github.com/materialsvirtuallab/matgenb#introduction) provides some example Jupyter notebooks that demonstrate how to use `pymatgen` functionality.

[pull request]: https://github.com/materialsproject/pymatgen/pulls
[github issue]: https://github.com/materialsproject/pymatgen/issues
[github discussion]: https://github.com/materialsproject/pymatgen/discussions

## Why use `pymatgen`?

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers and is the analysis code powering the [Materials Project]. The analysis it produces survives rigorous scrutiny every single day. Bugs tend to be found and corrected quickly. Pymatgen also uses Github Actions for continuous integration, which ensures that every new code passes a comprehensive suite of unit tests.
2. **It is well documented.** A fairly comprehensive documentation has been written to help you get to grips with it quickly.
3. **It is open.** You are free to use and contribute to `pymatgen`. It also means that `pymatgen` is continuously being improved. We will attribute any code you contribute to any publication you specify. Contributing to `pymatgen` means your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in `pymatgen` have been optimized by vectorizing in `numpy`/`scipy`. This means that coordinate manipulations are fast. Pymatgen also comes with a complete system for handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in the well-established Materials Project. It is also actively being developed and maintained by the [Materials Virtual Lab], the ABINIT group and many other research groups.
6. **A growing ecosystem of developers and add-ons**. Pymatgen has contributions from materials scientists all over the world. We also now have an architecture to support add-ons that expand `pymatgen`'s functionality even further. Check out the [contributing page](https://pymatgen.org/contributing) and [add-ons page](https://pymatgen.org/addons) for details and examples.

## Installation

The version at the Python Package Index [PyPI] is always the latest stable release that is relatively bug-free and can be installed via `pip`:

[pypi]: https://pypi.org/project/pymatgen

```sh
pip install pymatgen
```

If you'd like to use the latest unreleased changes on the main branch, you can install directly from GitHub:

```sh
pip install -U git+https://github.com/materialsproject/pymatgen
```

The minimum Python version is 3.10. Some extra functionality (e.g., generation of POTCARs) does require additional setup (see the [`pymatgen` docs]).

## Change Log

See [GitHub releases](https://github.com/materialsproject/pymatgen/releases), [`docs/CHANGES.md`](docs/CHANGES.md) or [commit history](https://github.com/materialsproject/pymatgen/commits/master) in increasing order of details.

## Using pymatgen

Please refer to the official [`pymatgen` docs] for tutorials and examples.

## How to cite pymatgen

If you use `pymatgen` in your research, please consider citing the following work:

> Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy
> Hautier, Michael Kocher, Shreyas Cholia, Dan Gunter, Vincent Chevrier,
> Kristin A. Persson, Gerbrand Ceder. *Python Materials Genomics
> (pymatgen): A Robust, Open-Source Python Library for Materials
> Analysis.* Computational Materials Science, 2013, 68, 314-319.
> [doi:10.1016/j.commatsci.2012.10.028](https://doi.org/10.1016/j.commatsci.2012.10.028)

In addition, some of `pymatgen`'s functionality is based on scientific advances/principles developed by the computational materials scientists in our team. Please refer to the [`pymatgen` docs] on how to cite them.

### Soliciting contributions to 2nd `pymatgen` paper

If you are a long-standing `pymatgen` contributor and would like to be involved in working on an updated `pymatgen` publication,
please fill out this [co-author registration form](https://docs.google.com/forms/d/e/1FAIpQLSecIhD2YjdPGldrRTM8Go3VxVg_vjKjZAOXtIKDG7qckHLYaQ/viewform) or contact [@shyuep, @mkhorton and @janosh](mailto:ongsp@ucsd.edu,m.k.horton@gmail.com,janosh@lbl.gov?subject=Contributing%20to%20updated%20pymatgen%20paper) with questions.

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

Shyue Ping Ong (@shyuep) of the [Materials Virtual Lab] started Pymatgen in 2011 and is still the project lead.
Janosh Riebesell (@janosh) and Matthew Horton (@mkhorton) are co-maintainers.

The [`pymatgen` development team] is the set of all contributors to the `pymatgen` project, including all subprojects.

## Our Copyright Policy

Pymatgen uses a shared copyright model. Each contributor maintains copyright over their contributions to `pymatgen`. But, it is important to note that these contributions are typically only changes to the repositories. Thus, the `pymatgen` source code, in its entirety is not the copyright of any single person or institution. Instead, it is the collective copyright of the entire [`pymatgen` Development Team]. If individual contributors want to maintain a record of what changes/contributions they have specific copyright on, they should indicate their copyright in the commit message of the change, when they commit the change to one of the `pymatgen` repositories.

[`pymatgen` docs]: https://pymatgen.org
[materials project]: https://materialsproject.org
[`pymatgen` development team]: https://pymatgen.org/team
[materials virtual lab]: https://materialsvirtuallab.org
