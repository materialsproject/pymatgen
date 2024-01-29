# https://github.com/materialsproject/pymatgen/issues/2591

from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

struct = Structure.from_file('CONTCAR')
struct = SpacegroupAnalyzer(struct).get_conventional_standard_structure()
slab = SlabGenerator(
    struct,
    miller_index=[0, 1, 0],
    min_slab_size=10,
    min_vacuum_size=15
    )

for n, s in enumerate(slab.get_slabs(bonds={('Al', 'O'): 2.1}, repair=True)):
    s = s.get_sorted_structure()
    s.to(filename='POSCAR-'+str(n), fmt='poscar')
