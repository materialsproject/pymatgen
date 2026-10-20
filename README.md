<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/materialsproject/pymatgen/raw/main/docs/assets/pymatgen-white.svg">
    <img alt="Logo" src="https://github.com/materialsproject/pymatgen/raw/main/docs/assets/pymatgen.svg"
height="70">
  </picture>
</h1>

<h4 align="center">

[![CI Status](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml/badge.svg)](https://github.com/materialsproject/pymatgen/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/materialsproject/pymatgen/branch/main/graph/badge.svg?token=XC47Un1LV2)](https://codecov.io/gh/materialsproject/pymatgen)
[![PyPI Downloads](https://img.shields.io/pypi/dm/pymatgen?logo=pypi&logoColor=white&color=blue&label=PyPI)](https://pypi.org/project/pymatgen)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/pymatgen?logo=condaforge&color=blue&label=Conda)](https://anaconda.org/conda-forge/pymatgen)
[![Paper](https://img.shields.io/badge/J.ComMatSci-2012.10.028-blue?logo=elsevier&logoColor=white)](https://doi.org/10.1016/j.commatsci.2012.10.028)

</h4>

Pymatgen (Python Materials Genomics) is a robust, open-source Python
library for materials analysis.

> [!IMPORTANT]
> The core data objects (`Element`, `Site`, `Molecule`, `Structure`, `Composition`, `Lattice`) and
> file I/O for major electronic structure codes now live in the lean [`pymatgen-core`] package.
> Installing `pymatgen` automatically pulls in `pymatgen-core` and adds the higher-level
> analyses, applications and external integrations described below. If you only need the core
> objects and I/O, depend on `pymatgen-core` directly for a much lighter install.

[`pymatgen-core`]: https://github.com/materialsproject/pymatgen-core

These are some of the main features provided by this package on top of `pymatgen-core`:

1. The full set of core classes from [`pymatgen-core`] (`Element`, `Site`, `Molecule`, `Structure`,
   `Composition`, `Lattice`) re-exported under the familiar `pymatgen.core` namespace, plus
   extensive input/output for [VASP](https://www.vasp.at/),
   [ABINIT](https://abinit.github.io/abinit_web/),
   [CIF](https://wikipedia.org/wiki/Crystallographic_Information_File),
   [Gaussian](https://gaussian.com),
   [XYZ](https://wikipedia.org/wiki/XYZ_file_format) and many other file formats.
2. Powerful analysis tools (`pymatgen.analysis`), including generation of phase diagrams,
   Pourbaix diagrams, reaction analysis, local environments, surfaces, interfaces, defects,
   magnetism, piezoelectricity, elasticity and much more.
3. Electronic structure analyses, such as density of states and band structure.
4. Entry/thermo data objects (`pymatgen.entries`) for thermodynamic and reaction calculations.
5. External integrations (`pymatgen.ext`), including the [Materials Project] REST API and
   other materials databases.
6. End-user applications and a command-line interface (`pymatgen.apps`, `pymatgen.cli`).

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

1. **It is (fairly) robust.** Pymatgen is used by thousands of researchers and is the analysis code powering the
   [Materials Project]. The analysis it produces survives rigorous scrutiny every single day. Bugs tend to be found
   and corrected quickly. Pymatgen also uses Github Actions for continuous integration, which ensures that every
   new code passes a comprehensive suite of unit tests.
2. **It is well documented.** A fairly comprehensive documentation has been written to help you get to grips with
   it quickly.
3. **It is open.** You are free to use and contribute to `pymatgen`. It also means that `pymatgen` is continuously
   being improved. We will attribute any code you contribute to any publication you specify. Contributing to
   `pymatgen` means your research becomes more visible, which translates to greater impact.
4. **It is fast.** Many of the core numerical methods in `pymatgen` have been optimized by vectorizing in
   `numpy`/`scipy`. This means that coordinate manipulations are fast. Pymatgen also comes with a complete system
   for handling periodic boundary conditions.
5. **It will be around.** Pymatgen is not a pet research project. It is used in the well-established Materials
   Project. It is also actively being developed and maintained by the [Materialyze] Lab, the ABINIT group and
   many other research groups.
6. **A growing ecosystem of developers and add-ons**. Pymatgen has contributions from materials scientists all over
   the world. We also now have an architecture to support add-ons that expand `pymatgen`'s functionality even
   further. Check out the [contributing page](https://pymatgen.org/contributing) and [add-ons page](https://pymatgen.org/addons) for details and examples.

## Installation

The version at the Python Package Index [PyPI] is always the latest stable release that is relatively bug-free and can be installed via `pip`:

[pypi]: https://pypi.org/project/pymatgen

```sh
pip install pymatgen
```

This installs both `pymatgen` (analyses, applications and external integrations) and
[`pymatgen-core`] (core data objects and file I/O). For a minimal install with just the core
objects and I/O, depend on `pymatgen-core` instead:

```sh
pip install pymatgen-core
```

Some extra functionality (e.g., generation of POTCARs) does require additional setup (see the [`pymatgen` docs]).

## Change Log

See [GitHub releases](https://github.com/materialsproject/pymatgen/releases), [`docs/CHANGES.md`](docs/CHANGES.md) or [commit history](https://github.com/materialsproject/pymatgen/commits/main) in increasing order of details.

## Using pymatgen

Please refer to the official [`pymatgen` docs] for tutorials and examples. Dr Anubhav Jain (@computron) has also created
a series of [tutorials](https://github.com/computron/pymatgen_tutorials) and [YouTube videos](https://www.youtube.com/playlist?list=PL7gkuUui8u7_M47KrV4tS4pLwhe7mDAjT), which is a good resource, especially for beginners.

## How to cite pymatgen

If you use `pymatgen` in your research, please cite the following [work](https://arxiv.org/abs/2610.23538):

```txt
Janine George, Matthew K. Horton, Aaron D. Kaplan, Seán R. Kavanagh, Janosh Riebesell, Andrew S. Rosen,
Haoyu Yang, Felix Adams, Mahdi Amachraa, Oxana Andriuc, Omar A. Ashour, Saurabh Bajaj, Ella Banyas, Luis
Barroso-Luque, Karlo Berket, Noam Bernstein, Arpit Bhardwaj, Ali Hussain Umar Bhatti, Samuel M. Blau,
Brandon Bocklund, Alexander Bonkowski, Danny Broberg, Guillaume Brunin, Kyle Bystrom, Charles A. Cardot,
Rees Chang, Chi Chen, Wei Chen, Xin Chen, Yuan Chiang, Kamal Choudhary, Orion A. Cohen, Miguel D. Costa,
Stephen T. Dacek III, Maarten de Jong, Tom Demeyere, Alan K. Dozier, Alexander Dunn, Shyam Dwaraknath,
Maxwell T. Dylla, Kristjan Eimre, Alexander R. Epstein, Christina Ertural, Marco Esters, Matthew L.
Evans, Johannes Falke, Adam Fekete, Bernard Field, Alexandr Fonari, Nathan C. Frey, Dale W. Gaines II,
Alex M. Ganose, Matteo Giantomassi, Rhys E. A. Goodall, Petre-Flaviu Gostin, Aayush Gupta, Viet-Anh Ha,
Geoffroy Hautier, Niclas Heinsdorf, Michael F. Herbst, Tingzheng Hou, Patrick Huck, Yuji Ikeda, Adam J.
Jackson, Saivenkataraman Jayaraman, Finn R. Kearns, Sunghyun Kim, Ryan S. Kingsbury, Adam M. Krajewski,
Yu Kumagai, Matthew C. Kuner, William T. Laderer, Weitang Li, Xiang-Guo Li, Yu Hsuan Liang, Hanyu Liu, Ke
Liu, Miao Liu, Ahmad Lutfi, Olivier Marchand, Matthew J. McDermott, Bharat K. Medasani, Fanchen Meng,
Kyle D. Miller, Henrique Miranda, Rohith Srinivaas Mohanakrishnan, Joseph H. Montoya, Guy C. Moore,
Benjamin J. Morgan, Jason M. Munro, Koki Muraoka, Aakash A. Naik, Willis O'Leary, Anthony Onwuli, Jaesuk
Park, Drew F. Parsons, Kai Pei, Guido Petretto, Thomas A. R. Purcell, Ji Qi, Xiaohui Qu, Ihor Radchenko,
Nav Nidhi Rajput, Francesco Ricci, Benjamin S. Rich, William D. Richards, Gian-Marco Rignanese, Xiaoxu
Ruan, Germain Salvato Vallverdu, Jonathan Schmidt, Jimmy-Xuan Shen, Ryuichi Shimogawa, Kohei Shinohara,
Arunima K. Singh, Eric Sivonxay, Stef Smeets, Tess E. Smidt, Andrei Sobolev, Evan Walter Clark
Spotte-Smith, Christopher Stihl, Wenhao Sun, Jack D. Sundberg, Yuta Suzuki, Leopold Talirz, Hanmei Tang,
Eric Taw, Cooper Tezak, Félix Therrien, Scott D. Thiel, Steven B. Torrisi, Mark E. Turiansky, Katharina
Ueltzen, Michiel J. van Setten, Josua Vieten, Sudarshan Vijay, Nicholas Wagner, Frank Wan, Amanda Wang,
Jingyang Wang, Zhenbin Wang, Logan Ward, Thomas Warford, David Waroquiers, Mingjian Wen, Oskar Weser,
Nicholas Winner, Donald Winston, Michael Wolloch, Brandon M. Wood, Rachel Woods-Robinson, Yihan Xiao,
Yaoshu Xie, Firat Yalcin, Ruo Xi Yang, Weike Ye, Pezhman Zarabadi-Poor, Jinzhe Zeng, Hongsheng Zhao,
Kaifeng Zheng, He Zhu, Zhuoying Zhu, Nils E. R. Zimmermann, Daniel Zügner, Kristin A. Persson, Anubhav
Jain, Shyue Ping Ong.
pymatgen: 15 years of community growth, new functionality, and future prospects.
arXiv, 2026, arXiv:2610.23538. doi:10.48550/arXiv.2610.23538
```

Please also consider citing the [original `pymatgen` paper](https://doi.org/10.1016/j.commatsci.2012.10.028):

```txt
Shyue Ping Ong, William Davidson Richards, Anubhav Jain, Geoffroy Hautier, Michael Kocher, Shreyas Cholia, Dan
Gunter, Vincent Chevrier, Kristin A. Persson, Gerbrand Ceder. Python Materials Genomics (pymatgen): A Robust,
Open-Source Python Library for Materials Analysis. Computational Materials Science, 2013, 68, 314-319.
doi:10.1016/j.commatsci.2012.10.028
```

In addition, some of `pymatgen`'s functionality is based on scientific advances/principles developed by the
computational materials scientists in our team. Please refer to the [`pymatgen` docs] on how to cite them.

## License

Pymatgen is released under the MIT License. The terms of the license are as follows:

```txt
The MIT License (MIT) Copyright (c) 2011-2012 MIT & LBNL

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of
the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

## About the Pymatgen Development Team

Shyue Ping Ong (@shyuep) of the [Materialyze] Lab started Pymatgen in 2011 and is still the project lead.

The [`pymatgen` development team] is the set of all contributors to the `pymatgen` project, including all subprojects.

## Our Copyright Policy

Pymatgen uses a shared copyright model. Each contributor maintains copyright over their contributions to `pymatgen`.
But, it is important to note that these contributions are typically only changes to the repositories. Thus, the
`pymatgen` source code, in its entirety is not the copyright of any single person or institution. Instead, it is the
collective copyright of the entire [`pymatgen` Development Team]. If individual contributors want to maintain a
record of what changes/contributions they have specific copyright on, they should indicate their copyright in the
commit message of the change, when they commit the change to one of the `pymatgen` repositories.

[`pymatgen` docs]: https://pymatgen.org
[materials project]: https://materialsproject.org
[`pymatgen` development team]: https://pymatgen.org/team
[materialyze]: https://materialyze.ai
