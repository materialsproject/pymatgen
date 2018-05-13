# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import division, unicode_literals, print_function

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


class XSF(object):
    """
    Class for parsing XCrysden files.
    """

    def __init__(self, structure):
        self.structure = structure

    def to_string(self):
        """
        Returns a string with the structure in XSF format
        See http://www.xcrysden.org/doc/XSF.html
        """
        lines = []
        app = lines.append

        app("CRYSTAL")
        app("# Primitive lattice vectors in Angstrom")
        app("PRIMVEC")
        cell = self.structure.lattice.matrix
        for i in range(3):
            app(' %.14f %.14f %.14f' % tuple(cell[i]))

        cart_coords = self.structure.cart_coords
        app("# Cartesian coordinates in Angstrom.")
        app("PRIMCOORD")
        app(" %d 1" % len(cart_coords))

        for a in range(len(cart_coords)):
            sp = "%d" % self.structure.atomic_numbers[a]
            app(sp + ' %20.14f %20.14f %20.14f' % tuple(cart_coords[a]))

        return "\n".join(lines)

    @classmethod
    def from_string(self, input_string, cls=None):
        """
        Initialize a `Structure` object from a string with data in XSF format.

        Args:
            input_string: String with the structure in XSF format.
                See http://www.xcrysden.org/doc/XSF.html
            cls: Structure class to be created. default: pymatgen structure

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
                for j in range(i + 1, i + 4):
                    lattice.append([float(c) for c in lines[j].split()])

            if "PRIMCOORD" in lines[i]:
                num_sites = int(lines[i + 1].split()[0])

                for j in range(i + 2, i + 2 + num_sites):
                    tokens = lines[j].split()
                    species.append(int(tokens[0]))
                    coords.append([float(j) for j in tokens[1:]])
                break
        else:
            raise ValueError("Invalid XSF data")

        if cls is None:
            from pymatgen.core.structure import Structure
            cls = Structure

        s = cls(lattice, species, coords, coords_are_cartesian=True)
        return XSF(s)


class BXSF(object):
    '''
    Class for parsing XCrysden fermi-surface files.
    inspared by https://github.com/MTD-group/Fermi-Surface
    '''

    def __init__(self, efermi, eigenvalues, kpts, lattice_rec):
        self.efermi = efermi
        self.eigenvalues = eigenvalues
        self.kpts = kpts
        self.lattice_rec = lattice_rec

    def to_string(self):
        n_bands = self.eigenvalues[Spin.up].shape[1]
        lines = []
        app = lines.append
        app('BEGIN_INFO')
        app('   #')
        app('   # Case:  unknown system')
        app('   #')
        app('   # Launch as: xcrysden --bxsf example.bxsf')
        app('   #')
        app('   Fermi Energy: {}'.format(self.efermi))
        app(' END_INFO')
        app('')
        app(' BEGIN_BLOCK_BANDGRID_3D')
        app(' Num_bands_are_sum_of_spin_up/down._Change_reci_dimension_based_on_your_unit_cell'
            )
        app('   BEGIN_BANDGRID_3D')
        app('       {}'.format(n_bands))
        app('     {} {} {}'.format(*self.kpts))
        app('     0.0000 0.0000 0.0000')
        for i in range(3):
            app('     {:.4f} {:.4f} {:.4f}'.format(
                *self.lattice_rec.matrix[i] / (2 * np.pi)))
        app('')
        bands_start = {Spin.up: 0, Spin.down: n_bands}
        for spin, v in self.eigenvalues.items():
            for i in range(n_bands):
                app('   BAND:  {}'.format(1 + i + bands_start[spin]))
                band_eigenvalues = np.reshape(
                    v[:, i, 0], (np.prod(self.kpts[:-1]), self.kpts[-1]))
                for j in range(self.kpts[-1]):
                    app('       {}'.format('  '.join(
                        map(str, band_eigenvalues[:, j]))))
                app('')
        app('   END_BANDGRID_3D')
        app(' END_BLOCK_BANDGRID_3D')
        app('')
        return '\n'.join(lines)

    def write_file(self, filename='Xcrysden.bxsf'):
        with open(filename, 'w') as f:
            f.write(self.to_string())

    def __str__(self):
        return self.to_string()

    @classmethod
    def from_vasprun(cls, vasprun):
        kpts = [
            len(set(np.array(vasprun.kpoints.kpts)[:, i])) for i in range(3)
        ]
        return cls(vasprun.efermi, vasprun.eigenvalues, kpts,
                   vasprun.lattice_rec)

    @classmethod
    def from_path(cls, path):
        vasprun, outcar = get_vasprun_outcar(path)
        return cls.from_vasprun(vasprun)
