Change log
==========

v2022.0.8
---------
* PR #2130 @rkingsbury ensures that energy corrections applied to each anion
  have unique names (e.g., N vs. Cl vs. Br).
* PR #2133 @rkingsbury adds support for custom vdW radii to `QCInput` and 
  `QChemDictSet`. These radii are used in the construction of PCM cavities and
  when calculating charges.  
* PR #2123 from @gpetretto fixes bug in `get_conventional_standard_structure` 
  method of the `SpacegroupAnalyzer` for triclinic crystals.
* PR #2134 from @ab5424 supports zopen in parsing lammps logs
* PR #2132 from @htz1992213 speeds up LammpsData.as_string for
  non-hybrid data with large coeff sections and adds as_lammpsdata method to
  CombinedData  
* PR #2129 from @richardtran415 improves analysis of surface symmetry of slabs.
* PR #2117 from @nwinner contains bug fixes for bader caller.    
