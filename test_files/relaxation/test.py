# coding: utf-8
from pymatgen.io.vasp.sets import MPOpticsNonSCFVaspInputSet
vis = MPOpticsNonSCFVaspInputSet(prev_calc_dir=".")

print vis.get_incar(vis.prev_run["structure"])
