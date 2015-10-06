Change log
==========

v3.2.3
------
* Massive update to abinit support. Note that pymatgen.io.abinitio has 
  been refactored to pymatgen.io.abinit. (Matteo, Setten)
* NwOutput now supports parsing of Hessian matrices (contributed by Xin 
  Chen)
* Gaussian support now has the ability to read potential energy surface
  and electronic transitions computed with TD-DFT (Germain Salvato 
  Vallverdu)
* Bug fixes for CifWriter with symmetry.
* Bug fixes for surface generation and reactions.
* Monty requirement increased.
