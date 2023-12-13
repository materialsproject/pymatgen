import unittest
from monty.io import zopen
from pymatgen.io.pwmat.inputs import ACExtractor, ACstrExtractor, AtomConfig
from pymatgen.io.pwmat.outputs import Movement


atom_config_path = "/data/home/liuhanyu/hyliu/pwmat_demo/CrI3/scf/atom.config"
#atom_config_path = "/data/home/liuhanyu/hyliu/pwmat_demo/MoS2/scf/atom.config"
ac_extractor = ACExtractor(file_path=atom_config_path)
with zopen(atom_config_path, "r") as f:
    ac_str_extractor = ACstrExtractor(atom_config_str="".join(f.readlines()))
    
    
movement_path = "/data/home/liuhanyu/hyliu/code/mlff/PWmatMLFF_dev/test/Cu/0_300_MOVEMENT"
movement_path = "/data/home/liuhanyu/hyliu/code/mlff/test/demo1/PWdata/data1/MOVEMENT"
movement = Movement(movement_path)


class ACExtractorTest(unittest.TestCase):    
    def test_all(self):
        self.assertEqual(ac_extractor.num_atoms, ac_extractor.get_num_atoms())
        for ii in range(9):
            self.assertEqual(ac_extractor.lattice[ii], ac_extractor.get_lattice()[ii])
        for ii in range(ac_extractor.num_atoms):
            self.assertEqual(ac_extractor.types[ii], ac_extractor.get_types()[ii])
            self.assertEqual(ac_extractor.coords[ii*3+0], ac_extractor.get_coords()[ii*3+0])
            self.assertEqual(ac_extractor.coords[ii*3+1], ac_extractor.get_coords()[ii*3+1])
            self.assertEqual(ac_extractor.coords[ii*3+2], ac_extractor.get_coords()[ii*3+2])
            self.assertEqual(ac_extractor.magmoms[ii], ac_extractor.get_magmoms()[ii])


class ACstrExtractorTest(unittest.TestCase):
    def test_all(self):
        self.assertEqual(ac_extractor.num_atoms, ac_str_extractor.get_num_atoms())
        for ii in range(9):
            self.assertEqual(ac_extractor.lattice[ii], ac_str_extractor.get_lattice()[ii])
        for ii in range(ac_extractor.num_atoms):
            self.assertEqual(ac_extractor.types[ii], ac_str_extractor.get_types()[ii])
            self.assertEqual(ac_extractor.coords[ii*3 + 0], ac_str_extractor.get_coords()[ii*3 + 0])
            self.assertEqual(ac_extractor.coords[ii*3 + 1], ac_str_extractor.get_coords()[ii*3 + 1])
            self.assertEqual(ac_extractor.coords[ii*3 + 2], ac_str_extractor.get_coords()[ii*3 + 2])
            self.assertEqual(ac_extractor.magmoms[ii], ac_str_extractor.get_magmoms()[ii])


class AtomConfigTest(unittest.TestCase):
    def test_all(self):
        atom_config = AtomConfig.from_file(filename=atom_config_path, mag=True)
        assert("magmom" in atom_config.structure.sites[0].properties)
        self.assertEqual("Cr2I6", atom_config.true_names)


class MovementTest(unittest.TestCase):
    def test_all(self):
        self.assertEqual(movement.nionic_steps, len(movement.chunksizes))
        self.assertEqual(movement.nionic_steps, len(movement.chunkstarts))
        self.assertEqual(movement.nionic_steps, len(movement.atom_configs))
        self.assertEqual(movement.atom_configs[0].structure.num_sites, movement.atom_configs[32].structure.num_sites)
        print(movement.ionic_steps[0].keys())
        #print(movement.virials[0])
        

if __name__ == "__main__":
    unittest.main()
