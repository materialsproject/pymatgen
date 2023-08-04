---
layout: default
title: Installation
nav_order: 2
---

## Requirements

All required dependencies should be automatically taken care of if you
install pymatgen using easy_install or pip. Otherwise, these packages should
be available on `PyPI <http://pypi.python.org>`_.

## Optional dependencies

Optional libraries that are required if you need certain features.

1. sympy: For defect generation and analysis.
2. VTK with Python bindings 5.8+ (<http://www.vtk.org/>): For visualization of
   crystal structures using the pymatgen.vis package. Note that the VTK
   package is incompatible with Python 3.x at the moment.
3. Atomistic Simulation Environment or ASE 3.6+: Required for the usage of the
   adapters in pymatgen.io.aseio between pymatgen's core Structure object and
   the Atoms object used by ASE. Get it at <https://wiki.fysik.dtu.dk/ase/>.
   Note that the ASE package is compatible with Python 3.5+ at the moment.
4. OpenBabel with Python bindings (<http://openbabel.org>): Required for the
   usage of the adapters in pymatgen.io.babelio between pymatgen's Molecule
   and OpenBabel's OBMol. Opens up input and output support for the very large
   number of input and output formats supported by OpenBabel.
5. networkx: For graph analysis associated with critic2 topological analysis
   of electron charge densities, pygraphviz is also required for visualization.
6. pytest - For unittesting. Not optional for developers.
7. numba: Optionally can be installed for faster evaluation of certain
   functionality, such as the SubstrateAnalyzer. It incurrs an initial
   slowdown the first time the relevant function is called, as it is
   compiled, for dramatically faster subsequent evaluations. Note that
   numba places additional constraints on the versions of numpy and
   Python that may be used.

## Optional non-Python programs

Optional non-python libraries (because no good python alternative exists at
the moment) required only for certain features:

1. ffmpeg: For generation of movies in structure_vtk.py. The executable ffmpeg
   must be in the path. Get it at <http://www.ffmpeg.org>.
2. enum: For the use of
   `pymatgen.transformations.advanced_transformations.EnumerateStructureTransformation`
   and `pymatgen.command_line.enumlib_caller` module. This library by Gus
   Hart provides a robust way to enumerate derivative structures. It can be
   used to completely enumerate all symmetrically distinct ordered structures
   of disordered structures via EnumerateStructureTransformation. Many other
   advanced transformations (e.g., MagOrderingTransformation) use
   EnumerateStructureTransformation. The enum.x and makestr.x
   executables must be in the path. Get it at <http://github.com/msg-byu/enumlib> and
   follow the instructions to compile enum.x and makestr.x.
3. bader: For use with :class:`pymatgen.command_line.bader_caller.BaderAnalysis`.
   This library by Henkelmann et al. provides a robust way to calculate the
   Bader analysis from a CHGCAR. The bader executable must be in the path.
   Get it at <http://theory.cm.utexas.edu/bader>.
4. gulp: For use with :mod:`pymatgen.command_line.gulp_caller`,
   which is in turn used extensively by :mod:`pymatgen.analysis.defects` to
   compute empirical defect energies.
5. aconvasp: For use with the :mod:`pymatgen.command_line.aconvasp_caller`.
6. Zeo++ (<http://zeoplusplus.org>): For defect structure
   generation. This is required in addition to installing the zeo Python
   package.
7. critic2 (<https://github.com/aoterodelaroza/critic2>): For topological
   analysis of critical points from electronic charge density. Provides
   more detailed information compared to bader. For use with
   :class:`pymatgen.command_line.critic2_caller.Critic2Caller`.
8. graphviz (<http://graphviz.org>): For visualization of graphs generated
   using critic2.

## Conda-based install

For these instructions, we will assume the **64-bit** versions of all OSes.
For OSX and Linux, both latest Python 3.x and 2.7 are supported. For Windows,
only latest Python 3.x is supported. Most common functionality should work
out of the box on Windows, but some specialized analyses relying on external
programs may require you to compile those programs from source.

### Step 1: Install conda

Download and install the version of conda for your operating system from
<http://conda.pydata.org/miniconda.html>. For Windows, **make sure it is the
Miniconda3 installer**, and simply double-click the exe file. For Linux or Mac,
run:

```bash
# If Mac
bash Miniconda3-latest-MacOSX-x86_64.sh

# If Linux
bash Miniconda3-latest-Linux-x86_64.sh
```

Note that you may need to create a new terminal after this step in order for
the environmental variables added by conda to be loaded.

### Step 2b: (Optional) Create a conda environment

If you are working with many python packages, it is generally recommended you
create a separate environment for each of your packages. For example:

```shell
conda create --name my_pymatgen python
source activate my_pymatgen  # OSX or Linux
activate my_pymatgen  # Windows
```

### Step 3: Install pymatgen

You can install pymatgen via conda as well via the conda-forge channel on
Anaconda cloud:

```shell
conda install --channel conda-forge pymatgen
```

If the above fails, try using conda to install some critical dependencies and
then do pip install::

```shell
conda install --yes numpy scipy matplotlib
pip install pymatgen
```

### Step 4: (Optional) Install enumlib and bader (only for OSX and Linux)

If you would like to use the enumeration capabilities powered by Gus Hart's
enumlib or perform Bader charge analysis powered by the Bader analysis code
of the Henkelmann group, please try installing these from source using the pmg
command line tool as follows::

```shell
pmg config --install enumlib
pmg config --install bader
```

Then put these in your PATH somewhere. You can also download the source of
these from the official repos and follow the compile instructions.

## POTCAR Setup

For the code to generate POTCAR files, it needs to know where the VASP
pseudopotential files are. We are not allowed to distribute these under the
VASP license. The good news is that the `pmg` command line utility includes a
config functionality.

After installation, do

```bash
pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>
```

In the above, `<EXTRACTED_VASP_POTCAR>` is the location of the directory that
you extracted the downloaded VASP pseudopotential files. Typically, it has
the following format:

```
 - <EXTRACTED_VASP_POTCAR>
 |- POT_GGA_PAW_PBE
 ||- Ac_s
 |||-POTCAR
 |||-...
```

or:

```
 - <EXTRACTED_VASP_POTCAR>
 |- potpaw_PBE
 ||- Ac_s
 |||-POTCAR
 |||-...
```

and follow the instructions. If you have done it correctly, you should get a
resources directory with the following directory structure::

```
- psp_resources
|- POT_GGA_PAW_PBE
||- POTCAR.Ac_s.gz
||- POTCAR.Ac.gz
||- POTCAR.Ag.gz
...
|- POT_GGA_PAW_PW91
...
```

After generating the resources directory, you should add a VASP_PSP_DIR config
variable pointing to the generated directory and you should then be
able to generate POTCARs:

```bash
pmg config --add PMG_VASP_PSP_DIR <MY_PSP>
```

If you are using newer sets of pseudopotential files from VASP, the directory
names may be different, e.g., POT_GGA_PAW_PBE_52. For such cases, please also
add a default functional specification as follows:

```bash
pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_52
```

You can also use this to specify whatever functional you would like to use by
default in pymatgen, e.g., LDA_52, PW91, etc. Type::

```bash
pmg potcar -h
```

to see full list of choices.

**Note:** The Materials Project currently uses older versions of the VASP pseudopotentials
    for maximum compatibility with historical data, rather than the current 52/54
    pseudopotentials. This setting can be overridden by the user if desired.
    As such, current versions of pymatgen will check the hashes of your pseudopotentials
    when constructing input sets to ensure the correct, compatible pseudopotential sets are
    used, so that total energies can be compared to those in the Materials Project database.
    If you use any functional other than PBE, note that you should not be combining results
    from these other functionals with Materials Project data. For up-to-date information
    on this, please consult the Materials Project documentation.

## PyPy Support

[PyPy](https://www.pypy.org) is an alternative Python interpreter for running Python code
and comes with significant speed improvements for common applications. However, historically,
fewer packages offer PyPy support.

It is possible to install and use pymatgen with the PyPy interpreter
but it comes with some important caveats:

* While it is usable, PyPy is not officially supported by pymatgen. We do not run our
  full test suite on PyPy and it's possible some parts of pymatgen will be broken.
* All of pymatgen's dependencies now support PyPy including numpy, scipy, and pandas,
  however matplotlib is difficult to install. If trying PyPy, the current advice
  is to remove the matplotlib dependency, however this means any modules using matplotlib
  will not be importable. The easiest way to install dependencies is using the
  `PyPy builds on conda-forge <https://conda-forge.org/blog/2020/03/10/pypy>`_. For spglib,
  cloning the repository and running ``python setup.py install`` manually is advised.
* Performance improvements are unpredictable. Since pymatgen makes heavy use of numpy
  and custom extensions where appropriate, many code hot spots have already been optimized.

We welcome any developers interested in expanding our PyPy support.

## Setup for Developers (using GitHub)

### Step 1: Preparing your system

#### Windows

1. Download Microsoft Visual Studio 2015 (the free Community Edition) is fine.
2. Install Visual Studio 2015, but *make sure that you select More Options ->
   Programming Languages -> Visual C++ during the installation process*. By
   default, Visual Studio does not install Visual C++, which is needed.

#### Mac OSX

1. Download and install Xcode. Afterwards, install the XCode command line
   tools by typing the following in a terminal::

      xcode-select --install

2. (Optional) Install gfortran. Get an installer at
   <http://gcc.gnu.org/wiki/GFortranBinaries#MacOS>.

#### Linux

1. Usually no preparation is needed as most of the standard compilers should
   already be available.

### Step 2: Install pymatgen in developmental mode

1. Make sure you have git and `git-lfs <https://git-lfs.github.com/>`_ installed.
   Clone the repo at <https://github.com/materialsproject/pymatgen>.

2. Run `git lfs install` in the cloned repo first.

3. In your root pymatgen repo directory, type (you may need to do this with root
   privileges)::

      pip install -e .

4. Install any missing python libraries that are necessary.

I recommend that you start by reading some of the unittests in the tests
subdirectory for each package. The unittests demonstrate the expected behavior
and functionality of the code.

Please read up on pymatgen's :doc:`coding guidelines </contributing>` before
you start coding. It will make integration much easier.

## Installation tips for optional libraries

This section provides a guide for installing various optional libraries used in
pymatgen. Some of the python libraries are rather tricky to build in certain
operating systems, especially for users unfamiliar with building C/C++ code.
Please feel free to send in suggestions to update the instructions based on
your experiences. In all the instructions, it is assumed that you have standard
gcc and other compilers (e.g., Xcode on Macs) already installed.

### VTK on Mac OS X (tested on v7.0)

The easiest is to install cmake from
<http://cmake.org/cmake/resources/software.html>.

Type the following::

    cd VTK (this is the directory you expanded VTK into)
    mkdir build
    cd build
    ccmake .. (this uses cmake in an interactive manner)

Press "t" to toggle advanced mode. Then press "c" to do an initial
configuration. After the list of parameters come out, ensure that the
PYTHON_VERSION is set to 3, the VTK_WRAP_PYTHON is set to ON, and
BUILD_SHARED_LIBS is set to ON. You may also need to modify the python
paths and library paths if they are in non-standard locations. For example, if
you have installed the official version of Python instead of using the
Mac-provided version, you will probably need to edit the CMakeCache Python
links. Example configuration for Python 3.5 installed using conda is given
below (only variables that need to be modified/checked are shown)::

    PYTHON_EXECUTABLE:FILEPATH=/Users/<username>/miniconda3/bin/python3
    PYTHON_INCLUDE_DIR:PATH=/Users/<username>/miniconda3/include/python3.5m
    PYTHON_LIBRARY:FILEPATH=/Users/<username>/miniconda3/lib/libpython3.5m.dylib
    VTK_INSTALL_PYTHON_MODULE_DIR:PATH=/Users/<username>/miniconda3/lib/python3.5/site-packages
    VTK_PYTHON_VERSION:STRING=3
    VTK_WRAP_PYTHON:BOOL=ON

Then press "c" again to configure and finally "g" to generate the required
make files After the CMakeCache.txt file is generated, type::

    make -j 4
    sudo make install

With any luck, you should have vtk with the necessary python wrappers
installed. You can test this by going into a python terminal and trying::

    import vtk

### OpenBabel Mac OS X (tested on v2.3.2)

**Anaconda install**

If you are using anaconda (and have pymatgen installed in your anaconda environment), you should be
able to install openbabel with a single command::

    conda install -c openbabel openbabel

**Manual install**

Openbabel must be compiled with python bindings for integration with pymatgen.
Here are the steps that I took to make it work:

1. Install cmake from <http://cmake.org/cmake/resources/software.html>.

2. Install pcre-8.33 from
   ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.33.tar.gz.

3. Install pkg-config-0.28 using MacPorts or from
   <http://pkgconfig.freedesktop.org/releases/pkg-config-0.28.tar.gz>.

4. Install SWIG from
   <http://prdownloads.sourceforge.net/swig/swig-2.0.10.tar.gz>.

5. Download openbabel 2.3.2 *source code* from
   <https://sourceforge.net/projects/openbabel/files/openbabel/2.3.2/>.

6. Download Eigen version 3.1.2 from
   <http://bitbucket.org/eigen/eigen/get/3.1.2.tar.gz>.

7. Extract your Eigen and openbabel source distributions::

    tar -zxvf openbabel-2.3.2.tar.gz
    tar -zxvf eigen3.tar.gz

8. Now you should have two directories. Assuming that your openbabel src is in
   a directory called "openbabel-2.3.2" and your eigen source is in a directory
   called "eigen3", do the following steps::

    mv openbabel-2.3.2 ob-src
    cd ob-src/scripts/python; rm openbabel.py openbabel-python.cpp; cd ../../..

9. Edit ob-src/scripts/CMakeLists.txt, jump to line 70, change “eigen2_define”
   to “eigen_define”.

10. Let's create a build directory::

        mkdir ob-build
        cd ob-build
        cmake -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DEIGEN3_INCLUDE_DIR=../eigen3 ../ob-src 2>&1 | tee cmake.out

11. Before proceeding further, similar to the VTK installation process in the
    previous section, you may also need to modify the CMakeCache.txt
    file by hand if your python paths and library paths if they are in
    non-standard locations. For example, if you have installed the official
    version of Python instead of using the Mac-provided version,
    you will probably need to edit the CMakeCache Python links. Example
    configuration for Python 2.7 is given below (only variables that need to
    be modified are shown)::

        //Path to a program.
        PYTHON_EXECUTABLE:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/bin/python

        //Path to a file.
        PYTHON_INCLUDE_DIR:PATH=/Library/Frameworks/Python.framework/Versions/2.7/Headers

        //Path to a library.
        PYTHON_LIBRARY:FILEPATH=/Library/Frameworks/Python.framework/Versions/2.7/lib/libpython2.7.dylib

12. If you are using Mavericks (OSX 10.9) and encounter errors relating to <tr1/memory>, you might also need to include
    the following flag in your CMakeCache.txt::

        CMAKE_CXX_FLAGS:STRING=-stdlib=libstdc++

13. Run make and install as follows::

        make -j2
        sudo make install

14. With any luck, you should have openbabel with python bindings installed.
    You can test your installation by trying to import openbabel from the
    python command line. Please note that despite best efforts,
    openbabel seems to install the python bindings into /usr/local/lib even
    if your Python is not the standard Mac version. In that case,
    you may need to add the following into your .bash_profile::

        export PYTHONPATH=/usr/local/lib:$PYTHONPATH

### Zeo++

If you use the defects analysis package, you will need to install Zeo++.

The download and installation instructions for Zeo++ can be found here: <http://www.zeoplusplus.org/>
