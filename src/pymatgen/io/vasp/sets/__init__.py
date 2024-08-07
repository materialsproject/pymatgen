"""Re-export all VASP input sets, for more convenient imports and to maintain backwards compatible imports
following the split up of sets.py into submodules in #3865.
"""

from __future__ import annotations

from pymatgen.io.vasp.sets.base import (
    MODULE_DIR,
    BadInputSetWarning,
    DictSet,
    UserPotcarFunctional,
    VaspInputGenerator,
    VaspInputSet,
    _load_yaml_config,
    batch_write_input,
    get_structure_from_prev_run,
    get_valid_magmom_struct,
)
from pymatgen.io.vasp.sets.lobster import LobsterSet
from pymatgen.io.vasp.sets.matpes import MatPESStaticSet
from pymatgen.io.vasp.sets.mit import MITMDSet, MITNEBSet, MITRelaxSet
from pymatgen.io.vasp.sets.mp import (
    MPAbsorptionSet,
    MPHSEBSSet,
    MPHSERelaxSet,
    MPMDSet,
    MPMetalRelaxSet,
    MPNMRSet,
    MPNonSCFSet,
    MPRelaxSet,
    MPScanRelaxSet,
    MPScanStaticSet,
    MPSOCSet,
    MPStaticSet,
)
from pymatgen.io.vasp.sets.mvl import (
    MVLElasticSet,
    MVLGBSet,
    MVLGWSet,
    MVLNPTMDSet,
    MVLRelax52Set,
    MVLScanRelaxSet,
    MVLSlabSet,
)
