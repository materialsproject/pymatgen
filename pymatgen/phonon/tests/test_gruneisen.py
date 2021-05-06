import os
import unittest
import copy
from pymatgen.phonon.gruneisen import GruneisenParameter, GruneisenPhononBandStructure
from pymatgen.phonon.plotter import GruneisenPhononBSPlotter, GruneisenPhononBandStructureSymmLine, GruneisenPlotter
from pymatgen.io.phonopy import get_gruneisen_ph_bs_symm_line, get_gruneisenparamter
from pymatgen.util.testing import PymatgenTest
from pymatgen.core import Structure
from pymatgen.io.phonopy import get_gruneisen_ph_bs_symm_line
from pymatgen.phonon.gruneisen import GruneisenPhononBandStructure
from pymatgen.symmetry.bandstructure import HighSymmKpath

def get_kpath(structure: Structure):
    """
    get high-symmetry points in k-space
    Args:
        structure: Structure Object
    Returns:
    """
    kpath = HighSymmKpath(structure, symprec=0.01)
    kpath_save = kpath.kpath
    labels = copy.deepcopy(kpath_save["path"])
    path = copy.deepcopy(kpath_save["path"])

    for ilabelset, labelset in enumerate(labels):
        for ilabel, label in enumerate(labelset):
            path[ilabelset][ilabel] = kpath_save["kpoints"][label]
    return kpath_save, path

class GruneisenParameterTest(PymatgenTest):
    def setUp(self) -> None:
        #bs = GruneisenPhononBandStructure()
        kpath = HighSymmKpath(Structure.from_file(os.path.join(PymatgenTest.TEST_FILES_DIR,"gruneisen/eq/POSCAR")),
                              symprec=0.01)
        structure = kpath.prim
        kpath_dict, kpath_concrete = get_kpath(structure)

        bs_symm_line = get_gruneisen_ph_bs_symm_line(
            gruneisen_path=os.path.join(PymatgenTest.TEST_FILES_DIR,"gruneisen/gruneisen_eq_plus_minus.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR,"gruneisen/eq/POSCAR"),
            labels_dict=kpath_dict["kpoints"],
            fit=True)

        plotter = GruneisenPhononBSPlotter(bs=bs_symm_line)
        plotter.get_plot_gs().show()

class GruneisenMeshTest(PymatgenTest):
    def setUp(self) -> None:
        gruneisenobject = get_gruneisenparamter(
            os.path.join(PymatgenTest.TEST_FILES_DIR,"gruneisen/gruneisen_mesh.yaml"),
            structure_path=os.path.join(PymatgenTest.TEST_FILES_DIR,"gruneisen/eq/POSCAR"))

        plotter = GruneisenPlotter(gruneisenobject)
        plt = plotter.get_plot(units="mev")
        plt.show()
        

    def test_test(self):

        pass

if __name__ == "__main__":
    unittest.main()

