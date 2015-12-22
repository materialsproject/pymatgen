

from pymatgen.core.structure import Structure

class XSF(object):
    """
    Class for parsing XCrysden files.

    TODO: Unittests. Write XCrysden output.
    """

    def __init__(self, structure):
        self.structure = structure

    @classmethod
    def from_string(self, input_string):
        """
        Initialize a `Structure` object from a string with data in XSF format.
        See http://www.xcrysden.org/doc/XSF.html
        """
        # CRYSTAL                                        see (1)
        # these are primitive lattice vectors (in Angstroms)
        # PRIMVEC
        #    0.0000000    2.7100000    2.7100000         see (2)
        #    2.7100000    0.0000000    2.7100000
        #    2.7100000    2.7100000    0.0000000

        # these are conventional lattice vectors (in Angstroms)
        # CONVVEC
        #    5.4200000    0.0000000    0.0000000         see (3)
        #    0.0000000    5.4200000    0.0000000
        #    0.0000000    0.0000000    5.4200000

        # these are atomic coordinates in a primitive unit cell  (in Angstroms)
        # PRIMCOORD
        # 2 1                                            see (4)
        # 16      0.0000000     0.0000000     0.0000000  see (5)
        # 30      1.3550000    -1.3550000    -1.3550000

        lattice, coords, species = [], [], []
        lines = input_string.splitlines()

        for i in range(len(lines)):
            if "PRIMVEC" in lines[i]:
                for j in range(i+1, i+4):
                    lattice.append([float(c) for c in lines[j].split()])

            if "PRIMCOORD" in lines[i]:
                num_sites = int(lines[i+1].split()[0])

                for j in range(i+2, i+2+num_sites):
                    tokens = lines[j].split()
                    species.append(int(tokens[0]))
                    coords.append([float(j) for j in tokens[1:]])
                break
        else:
            raise ValueError("Invalid XSF data")

        s = Structure(lattice, species, coords, coords_are_cartesian=True)
        return XSF(s)