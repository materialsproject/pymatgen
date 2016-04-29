Change log
==========

v3.5.0
------
* Chemical environment analysis package (David Waroquiers).
* Piezoelectric property analysis (Shayam).
* Cythonize certain expensive core functions. 5-10x speedup in large structure matching (Will Richards).
* New NMR parsing functionality for Outcar (Xiaohui Qu).
* Improved io.lammps (Kiran Mathews).
* Update to spglib 1.9.2.
* Element properties now return unitized float where possible.
* Bug fix for get_primitive_standard affecting rhombohedral cells (important for band structures).
* Vasprun.final_energy now returns corrected energy with warning if it is different from final electronic step.
