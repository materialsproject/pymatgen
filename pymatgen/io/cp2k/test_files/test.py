from pymatgen.io.cp2k.sets import *
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cp2k.outputs import Cp2kOuput
from pymatgen.core.structure import Structure
from custodian.cp2k.jobs import Cp2kJob

if __name__ == '__main__':
    structure = Structure(lattice=[[3.349399,0,1.933776],[1.116466, 3.157843, 1.933776],[0,0,3.867552]],
                        species=['S','Si'], coords=[[0,0,0],[1.11646617, 0.7894608, 1.93377613]])

    s = HybridStaticSet(structure, sections={'GLOBAL': {'RUN_TYPE': 'ENERGY'}})
    s.write_file('cp2k.inp')

    out = Cp2kOuput(filename='cp2k.out', verbose=False, auto_load=False)
    out._convergence()
    print(out.data)
