import os, sys, time, shutil

base_path = os.getcwd()

for i in range(0, 36):

	mkfol = 'F' + str(i)
	os.mkdir(mkfol)
	shutil.copy('/global/u2/m/maartend/materials_project/strain/test4/F' + str(i) + '/vasprun.xml', mkfol) 
	
#	shutil.copy('POSCAR_' + str(i), mkfol + '/POSCAR')
#	shutil.copy('POTCAR', mkfol)
#	shutil.copy('INCAR', mkfol)
#	shutil.copy('KPOINTS', mkfol)
#	shutil.copy('submit_Franklin.pbs', mkfol)

#for i in range(0, 36):

#	mkfol = 'F' + str(i)
#	os.chdir(mkfol)	
#	os.system('vasp-openmpi.ib-intel09.x86_64 >> vasp.out &')
#	time.sleep(1000)
#	os.system('qsub submit_Franklin.pbs')
#	os.chdir(base_path)
