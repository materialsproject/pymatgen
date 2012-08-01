import CifFile
import abc

from   pymatgen.io.feffio_set                import *
from   pymatgen.io.vaspio                    import *
from   pymatgen.io.feffio                    import *
from   pymatgen.io.cifio                     import CifParser, CifWriter
from   pymatgen.core.structure               import Structure, Site, PeriodicSite


cif_file='../test_files/CeO_10688.cif'
central_atom='O'

r=CifParser(cif_file)
structure=r.get_structures()[0]
x=FeffInputSet("MaterialsProject")

header = FeffInputSet.get_header(x,structure, cif_file)
print "\n\nHEADER\n"
print header

tags=FeffInputSet.get_fefftags(x,"XANES")
print "\n\nPARAMETERS\n"
print tags

POT=FeffInputSet.get_feffPot(x,structure, central_atom)
print "\n\nPOTENTIALS\n"
print POT

ATOMS=FeffInputSet.get_feffAtoms(x,structure, central_atom)
print"\n\nATOMS\n"
print ATOMS

FeffInputSet.write_input(x, structure, "XANES", cif_file, "./fefftest", central_atom)
