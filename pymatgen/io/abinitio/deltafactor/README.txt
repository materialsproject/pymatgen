README FILE
-----------

This file provides a summary of the necessary information to calculate a code's Delta factor. The full procedure and background are described in the article by K. Lejaeghere, V. Van Speybroeck, G. Van Oost, and S. Cottenier: "Error bars for solid-state density-functional theory predictions: an overview by means of the ground-state elemental crystals" (2012).

In addition to this README.txt, the following files should have been included with the Supplementary Material:
- calcDelta.py
- CIFs.tar.gz
- eosfit.py
- exp.txt
- history.tar.gz
- VASP.txt
- VASP-relaxed.txt
- WIEN2k.txt
In order to be able to run calcDelta.py and eosfit.py, python and numpy are needed.

CONTENTS OF THE SUPPLEMENTARY MATERIAL
--------------------------------------

CIFs.tar.gz contains the CIF files of 71 elemental crystals in their ground state structure (except for sulfur and manganese, where a simpler structure that appears at higher temperatures or pressures is taken). For each of these structures, the CIF file presents the crystal at its calculated equilibrium (PBE functional). Spin-orbit coupling has not been incorporated in these calculations.

calcDelta.py is a python script, that calculates the Delta factor between a given code and WIEN2k. Results are saved in Delta-out.txt. The script is called by "python calcDelta.py filename", where 'filename' refers to a file containing V0, B0, and B1 information from the code under study for each elemental crystal. Other options allow to compare 'filename' to another code, or to print output to the screen. By typing "python calcDelta.py --help" a more extensive help can be displayed.

WIEN2k.txt is an input file for calcDelta.py. It is not called in the program arguments by default, however. It consists of similar information as 'filename', but for WIEN2k. These data are used as a reference. The WIEN2k.txt file should be in the work directory at all times in order for calcDelta.py to work.

exp.txt, VASP.txt and VASP-relaxed.txt are three example input files, containing the equation of state properties of the elemental crystals according to experimental or VASP (fixed- and relaxed-geometry respectively) computed data. More examples, as well as an extended history for all code data, are available in history.tar.gz.

eosfit.py is a python script, that calculates the Birch-Murnaghan equation of state for a single set of E(V) data points. It is called by "python eosfit.py filename2", where 'filename2' refers to a file containing the crystal volumes (in A^3/atom) and the corresponding energies (in eV/atom). By typing "python eosfit.py --help" a more extensive help can be printed on screen.

CALCULATING THE DELTA FACTOR
----------------------------

Unpack the CIFs.tar.gz file and translate every CIF file into a structure input file of the code under investigation. Make similar input files for 94%, 96%, 98%, 102%, 104%, and 106% of the equilibrium volume. For each elemental crystal this yields seven structures.

Calculate the energy for each of the 71x7 input files. Include spin polarization for O, Cr and Mn (antiferromagnetic) and Fe, Co, and Ni (ferromagnetic). The antiferromagnetism in the oxygen crystal should take place between the O2 molecules. Make sure the energy has converged in function of the computational parameters. Force optimizations should not be employed, since the computation of Delta requires fixed cell geometries, both with respect to the cell shape and the internal coordinates. Put the energies next to the volumes in a separate text file for each element (expressed per atom) and use eosfit.py to determine the Birch-Murnaghan equation of state parameters. If the predicted equilibrium volume is smaller than 96% or larger than 104% of the WIEN2k equilibrium volume, calculate an additional number of volumes so again 7 volume points are obtained, approximately symmetrical around V0. This ensures that for the pressure derivative of the bulk modulus an accurate value is obtained.

Assemble all sets of Birch-Murnaghan parameters in one text file and put the corresponding element symbol in front of the data. Use calcDelta.py to determine the Delta factor of the code under consideration, compared to WIEN2k.

More elaborate comparisons can be performed as well. A Linux-friendly output scheme has been provided, which allows to execute such analyses. To scale the elementwise Delta (from filename) to the difference with experiment, for example, type:
"python calcDelta.py VASP-relaxed.txt exp.txt --stdout | grep -E -v "\#|np|-" > test.txt"
"python calcDelta.py filename WIEN2k.txt --stdout | grep -E -v "\#|np|-" | join test.txt - | grep -v "N/A" | awk '{print $1, $3/$2*100}' | sort -n -k 2"

To compare the performance of the code from filename to that of VASP (fixed-geometry), type:
"python calcDelta.py VASP.txt WIEN2k.txt --stdout | grep -E -v "\#|np|-" > test.txt"
"python calcDelta.py filename WIEN2k.txt --stdout | grep -E -v "\#|np|-" | join test.txt - | grep -v "N/A" | awk '{print $1, $3/$2}' | sort -n -k 2"

FURTHER INFORMATION
-------------------

K. Lejaeghere, V. Van Speybroeck, G. Van Oost, and S. Cottenier: "Error bars for solid-state density-functional theory predictions: an overview by means of the ground-state elemental crystals" (2012).

Corresponding author: <Stefaan.Cottenier@UGent.be>

