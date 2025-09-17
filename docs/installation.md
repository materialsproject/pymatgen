---
layout: default
title: Installation
nav_order: 2
---

## Requirements

All required dependencies should be automatically taken care of if you install pymatgen using pip. Otherwise, these packages should be available on [PyPI](http://pypi.python.org).

## Optional dependencies

Optional libraries that are required if you need certain features.

1. sympy: For defect generation and analysis.
2. VTK with Python bindings 5.8+ (<http://www.vtk.org/>): For visualization of crystal structures using the pymatgen.vis package. Note that the VTK package is incompatible with Python 3.x at the moment.
3. Atomistic Simulation Environment or ASE 3.6+: Required for the usage of the adapters in pymatgen.io.aseio between pymatgen's core Structure object and the Atoms object used by ASE. Get it at <https://wiki.fysik.dtu.dk/ase/>. Note that the ASE package is compatible with Python 3.5+ at the moment.
4. OpenBabel with Python bindings (<http://openbabel.org>): Required for the usage of the adapters in pymatgen.io.babelio between pymatgen's Molecule and OpenBabel's OBMol. Opens up input and output support for the very large number of input and output formats supported by OpenBabel.
5. networkx: For graph analysis associated with critic2 topological analysis of electron charge densities, pygraphviz is also required for visualization.
6. pytest - For unittesting. Not optional for developers.
7. numba: Optionally can be installed for faster evaluation of certain functionality, such as the SubstrateAnalyzer. It incurrs an initial slowdown the first time the relevant function is called, as it is compiled, for dramatically faster subsequent evaluations. Note that numba places additional constraints on the versions of numpy and Python that may be used.

## Optional non-Python programs

Optional non-python libraries (because no good python alternative exists at the moment) required only for certain features:

1. `ffmpeg`: For generation of movies in structure_vtk.py. The executable ffmpeg must be in the path. Get it at <http://www.ffmpeg.org>.
2. `enum`: For the use of `pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation` and `pymatgen.command_line.enumlib_caller` module. This library by Gus Hart provides a robust way to enumerate derivative structures. It can be used to completely enumerate all symmetrically distinct ordered structures of disordered structures via EnumerateStructureTransformation. Many other advanced transformations (e.g., MagOrderingTransformation) use EnumerateStructureTransformation. The enum.x and makestr.x executables must be in the path. Get it at <http://github.com/msg-byu/enumlib> and follow the instructions to compile enum.x and makestr.x.
3. `bader`: For use with :class:`pymatgen.command_line.bader_caller.BaderAnalysis`. This library by Henkelmann et al. provides a robust way to calculate the Bader analysis from a CHGCAR. The bader executable must be in the path. Get it at <http://theory.cm.utexas.edu/bader>.
4. `gulp`: For use with :mod:`pymatgen.command_line.gulp_caller`, which is in turn used extensively by :mod:`pymatgen.analysis.defects` to compute empirical defect energies.
5. `aconvasp`: For use with the :mod:`pymatgen.command_line.aconvasp_caller`.
6. [Zeo++](https://zeoplusplus.org): For defect structure generation. This is required in addition to installing the zeo Python package.
7. [`critic2`](https://github.com/aoterodelaroza/critic2): For topological analysis of critical points from electronic charge density. Provides more detailed information compared to bader. For use with :class:`pymatgen.command_line.critic2_caller.Critic2Caller`.
8. [`graphviz`](https://graphviz.org): For visualization of graphs generated using critic2.

## Conda-based install

For these instructions, we will assume the **64-bit** versions of all OSes. For OSX and Linux, both latest Python 3.x and 2.7 are supported. For Windows, only latest Python 3.x is supported. Most common functionality should work out of the box on Windows, but some specialized analyses relying on external programs may require you to compile those programs from source.

### Step 1: Install conda

Download and install the version of conda for your operating system from <http://conda.pydata.org/miniconda.html>. For Windows, **make sure it is the Miniconda3 installer**, and simply double-click the exe file. For Linux or Mac, run:

```bash
# If Mac
bash Miniconda3-latest-MacOSX-x86_64.sh

# If Linux
bash Miniconda3-latest-Linux-x86_64.sh
```

Note that you may need to create a new terminal after this step for the environmental variables added by `conda` to be loaded.

### Step 2b: (Optional) Create a conda environment

If you are working with many Python packages, it is generally recommended you create a separate environment for each of your packages. For example:

```shell
conda create --name my_pymatgen python
source activate my_pymatgen  # OSX or Linux
activate my_pymatgen  # Windows
```

### Step 3: Install pymatgen

You can install pymatgen via conda as well via the conda-forge channel on Anaconda cloud:

```shell
conda install --channel conda-forge pymatgen
```

If the above fails, try using conda to install some critical dependencies and then do pip install::

```shell
conda install --yes numpy scipy matplotlib
pip install pymatgen
```

### Step 4: (Optional) Install `enumlib` and `bader` (only for OSX and Linux)

If you would like to use the enumeration capabilities powered by Gus Hart's `enumlib` or perform Bader charge analysis powered by the Bader analysis code of the Henkelmann group, please try installing these from source using the pmg command line tool as follows::

```shell
pmg config --install enumlib
pmg config --install bader
```

Then put these in your PATH somewhere. You can also download the source of these from the official repos and follow the compile instructions.

## POTCAR Setup

For the code to generate POTCAR files, it needs to know where the VASP pseudopotential files are. We are not allowed to distribute these under the VASP license. The good news is that the `pmg` command line utility includes a config functionality.

After installation, do:

```bash
pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>
```

In the above, `<EXTRACTED_VASP_POTCAR>` is the path to the extracted VASP pseudopotential files as obtained from VASP, and `<MY_PSP>` is the desired path where you would like to store the reformatted, Pymatgen-compatible pseudopotential files. Typically, the `<EXTRACTED_VASP_POTCAR>` directory has the following format:

```txt
potpaw_PBE.54
├── Ac
│   ├── POTCAR
│   └── PSCTR
├── Ag
│   ├── POTCAR
│   └── PSCTR
...
```

If you have done it correctly, your newly generated directory given by `<MY_PSP>` should have the following directory structure:

```txt
<MY_PSP>
├── POT_GGA_PAW_PBE_54
│   ├── POTCAR.Ac.gz
│   ├── POTCAR.Ag.gz
    ...
```

After the `<MY_PSP>` directory is generated, you should add it to your Pymatgen configuration file as follows:

```bash
pmg config --add PMG_VASP_PSP_DIR <MY_PSP>
```

In practice, this entire process might look something like the following:

```bash
pmg config -p /path/to/pseudos/potcar_PBE.54 /path/to/pseudos/pmg_potcars
pmg config -p /path/to/pseudos/potcar_LDA.54 /path/to/pseudos/pmg_potcars
pmg config --add PMG_VASP_PSP_DIR /path/to/pseudos/pmg_potcars
```

If desired, you may specify a default version and type of pseudopotentials as follows:

```bash
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_52
```

For additional options, run the help command to see the full list of choices.

**Note:** The Materials Project currently uses older versions of the VASP pseudopotentials for maximum compatibility with historical data, rather than the current 52/54 pseudopotentials. This setting can be overridden by the user if desired. As such, current versions of pymatgen will check the hashes of your pseudopotentials when constructing input sets to ensure the correct, compatible pseudopotential sets are used so that total energies can be compared to those in the Materials Project database. If you use any functional other than PBE, note that you should not be combining results from these other functionals with Materials Project data. For up-to-date information on this, please consult the Materials Project documentation.

## Contributing to `pymatgen`

See [pymatgen.org/contributing](https://pymatgen.org/contributing).
