from pymatgen.io.vasp.outputs import Vasprun
vasprun=Vasprun(filename="/hpc-user/jgeorge/PycharmProjects/2024_03_13_fatband/2024_03_13_fatband/pymatgen/tests/files/cohp/Fatband_SiO2/Test_p_x/vasprun.xml",
            ionic_step_skip=None,
            ionic_step_offset=0,
            parse_dos=True,
            parse_eigen=False,
            parse_projected_eigen=False,
            parse_potcar_file=False,
            occu_tol=1e-8,
            exception_on_bad_xml=True,
        )

efermi=vasprun.efermi
print(efermi)
#print(structure.to("/hpc-user/jgeorge/PycharmProjects/2024_03_13_fatband/2024_03_13_fatband/pymatgen/tests/files/cohp/Fatband_SiO2/Test_p_x/POSCAR"))