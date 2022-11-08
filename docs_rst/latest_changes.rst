Change log
==========

v2022.11.7
----------
* PR #2724 from @janosh: raise ValueError in ``SpacegroupAnalyzer.get_symmetrized_structure()`` if spglib returns no symmetries
* PR #2720 by @utf: Fix tensor mapping
* PR #2562 from @sudarshanv01: In case the Fock-matrix and eigenvalues are requested by the user (though the flags `scf_final_print` or `scf_print`), outputs.py now allows parsing both these quantities.
