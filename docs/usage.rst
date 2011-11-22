Basic usage
===========

Some example scripts have been provided in the scripts directory. In general, most file format conversions, manipulations and io can be done with a few quick lines of code. For example, to read a POSCAR and write a cif:

::

   from pymatgen.io.vaspio import Poscar
   from pymatgen.io.cifio import CifWriter

   p = Poscar('POSCAR')
   w = CifWriter(p.struct)
   w.write_file('mystructure.cif')

For more examples, please take a look in the scripts directory. More examples will be added soon.