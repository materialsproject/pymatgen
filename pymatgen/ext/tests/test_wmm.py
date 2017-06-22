import unittest

class WmmTest(PymatgenTest):
    def test_get_wmm_kpoints(self):
        si = PymatgenTest.get_structure("Si")
        file_name = os.path.join(test_dir, 'INCAR')
        incar = Incar.from_file(file_name)
        kpoints = Kpoints.get_from_wmm(si, incar=incar)


if __name__ == '__main__':
    unittest.main()
