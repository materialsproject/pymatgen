## Introduction ##

This readme provides a guide for installing various optional libraries used in pymatgen.  Some of the python libraries are rather tricky to build in certain operating systems, especially for users unfamiliar with building C/C++ code. Please feel free to update the instructions based on your experiences.

## Scipy (tested on v0.10.1) ##

### Mac OS X 10.7 ###

Typical installation of Xcode with python setup.py install seems to work fine. The pre-compiled binary for OSX 10.6 also seems to work.

### Solaris 10 ###

First install solstudio 12.2. Then put the following code in a shell script and run it.

	#!/bin/bash
	PATH=/opt/solstudio12.2/bin:/usr/ccs/bin:/usr/bin:/usr/sfw/bin:/usr/sbin; export PATH
	ATLAS=None; export ATLAS
	BLAS=/opt/solstudio12.2/lib/libsunperf.so; export BLAS
	LAPACK=/opt/solstudio12.2/lib/libsunmath.so; export LAPACK
	python setup.py build
	python setup.py install
	
## Qhull (tested on 2012.1) ##

### Mac OS X 10.7 ###

Typical installation with make fails with the following error:

	cc1plus: error: unrecognized command line option "-Wno-sign-conversion"

Simply removing "-Wno-sign-conversion" where it appears in the Makefile and then doing make followed by make install works fine.

## VTK (tested on v5.8.0) ##

### Mac OS X 10.7 ###

The easiest is to install cmake from http://cmake.org/cmake/resources/software.html

Type the following:
	cd VTK (this is the directory you expanded VTK into)
	cmake -i (this uses cmake in an interactive manner)

For all options, use the defaults, EXCEPT for BUILD_SHARED_LIBS and VTK_WRAP_PYTHON which must be set to ON. You may also need to modify the python paths and library paths if they are in non-standard locations.  After the CMakeCache.txt file is generated, type:

	make (note that this takes a while)
	sudo make install
	
With any luck, you should have vtk with the necessary python wrappers installed.
