import warnings

warnings.warn("pymatgen.structure_prediction and submodules has been moved to "
              "pymatgen.analysis.structure_prediction. This stub will be "
              "removed in pmg 2018.01.01.",
              DeprecationWarning)

from pymatgen.analysis.structure_prediction.substitution_probability import *

