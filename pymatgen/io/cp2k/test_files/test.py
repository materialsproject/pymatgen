from pymatgen.io.cp2k.sets import *
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cp2k.outputs import Cp2kOuput
from pymatgen.core.structure import Structure
from custodian.cp2k.jobs import Cp2kJob

if __name__ == '__main__':
    with MPRester() as mp:
        structure = mp.get_structure_by_material_id('mp-149', conventional_unit_cell=True)

    s = HybridRelaxSet(structure)
    s.print_pdos()
    s.print_e_density()
    s.print_hartree_potential()
    s.print_structures()
    s.print_mo_cubes()

    s.write_file('cp2k.inp')