Change log
==========

v2023.1.30
----------
* PR #2806 from @samblau qchem
    - Major changes to Q-Chem IO (inputs.py and outputs.py) to accommodate differences and new features in version 6+
    - Additional parsing capabilities for HOMO/LUMO, dipoles, NBO info (hyperbonds and 3C bonds) in outputs.py
    - Utility for processing a parsed binary Hessian scratch file
    - Overdue updates to default values in sets.py and new defaults associated with differences and new features in Q-Chem 6+
* PR #2814 from @jmmshn patch_dos
    ## Added Convenience to obtain the normalized CompleteDos object
    Added tests to make sure calling it multiple time still only gives one result.
