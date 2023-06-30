Change log
==========

v2023.6.28
----------
* Use lru_cache to speed up get_el_sp by 400x (@v1kko).
* Related to lru_cache of get_el_sp, Species.properties is now deprecated in favor of setting Species(spin=5). The rationale is
  that spin is the only supported property for Species anyway. Species and DummySpecies is now mostly immutable, i.e., setting specie.spin = 5 have no effect. This is as intended since the first version of pymatgen.
* PR #3111 from @xjf729 fix-MoleculeGraph-draw_graph
* PR #3030 from @lbluque Remove superfluous structure argument docstring from `SQSTransformation` init
* PR #3031 from @kavanase Quick fix to allow direct initialisation of the `DictSet` class.
* PR #3015 from @lbluque Optimized cython code in `find_points_in_spheres`, getting ~5x faster runtime.
