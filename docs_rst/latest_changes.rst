Change log
==========

v2022.0.9
---------
* Significant new functionality for handling interfaces between structures (@shyamd, #2149)
* Add input/output for CREST (@arepstein, #2020)
* Add RadialSiteDistortionTransformation (@nwinner, #2108)
* Add Q-Chem NBO functionality (@samblau, #2174)
* Change hkl annotation format in diffraction plots (@flaviu-gostin, #2143)
* Add space group to print output of `SymmetrizedStructure` (@CompRhys, #2139)
* Better error handling in QCOutput (@rkingsbury, #2147, #2165, #2135)
* Add progress bar for applying compatibility scheme (@CompRhys, #2136)
* Allow combining data with multiple molecule IDs in LAMMPS (@htz1992213, #2157)
* Update EDIFF in DFPT input set to be consistent with atomate (@utf, #2172)

* Change names of high-symmetry paths (@munrojm, #2144)
* Change default for filter_solids argument of PourbaixDiagram (@rkingsbury, #2177)

* Fix to improve precision in `FermiDos`, NOTE: this can result in significant changes in some instances (@nwinner, #2109)
* Fix for handling of Exceptions (@kmu, #2150)
* Fix for PourbaixEntry (@JosephMontoya-TRI, #2148)
* Fix for loading of settings from file when environment variables also set (@ardunn, #2164)
* Fix equation for calculation of k-spacing in SCAN sets, NOTE: this now results in a lower k-point density (@ab5424, #2163)
* Fix for parsing of VASP vasprun.xml when ALGO=CHI (@KazMorita, #2171)

* Documentation update for MP2020 corrections scheme (@rkingsbury, #2141)
* Documentation update for SCAN sets (@janosh, #2140)
* Documentation update for using CifWriter (@755452800, #2156)
