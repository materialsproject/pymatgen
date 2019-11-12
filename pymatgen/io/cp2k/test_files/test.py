from pymatgen.io.cp2k.sets import RelaxSet
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cp2k.outputs import Cp2kOuput


if __name__ == '__main__':

    with MPRester() as mp:
        structure = mp.get_structure_by_material_id('mp-149')

    s = RelaxSet(structure=structure)
    s.write_input('cp2k.inp')

    out = Cp2kOuput(filename='cp2k.out', verbose=False)
