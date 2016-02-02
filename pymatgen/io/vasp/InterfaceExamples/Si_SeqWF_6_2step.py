#!/usr/bin/env python

import os
import sys
import string
import json

from pymatgen import Structure
from fireworks import Firework, Workflow, LaunchPad
from pymatgen.io.vasp.interfaces import VaspInput, VaspFirework, VaspWorkflow
from fireworks.user_objects.firetasks.script_task import ScriptTask


# get structure for use in VASP input
material='Si_example.cif'
s=Structure.from_file(material)
pf=s.formula

material_file_name="{}_tddft.yaml".format(pf)
mname="{}_tddft".format(pf)
print "creating VASP jobs specifications in ",material_file_name


# and make it dependent on the number of cores
ncores=512
kpar_dft=32
kpar_chi=32
ncore_per_kpoint_dft=ncores/kpar_dft
ncore_per_kpoint_chi=ncores/kpar_chi
remainder=ncore_per_kpoint_dft%ncore_per_kpoint_chi

fw1_task1=ScriptTask.from_str("pwd")
fw1=VaspFirework(fw1_task1,name=mname)

# standard DFT first
vasp1=VaspInput(s)
nelect=vasp1.get_nelect()
nbands=nelect*3/2
remainder=nbands%ncore_per_kpoint_dft
add=ncore_per_kpoint_dft-remainder
nbands=int(nbands+add)
vasp1.ALGO='Normal'
vasp1.KPAR=kpar_dft
vasp1.NBANDS=nbands
vasp1.NEDOS=4096
vasp1.kpts_style='Gamma'
vasp1.kpts=[6,6,6]
vasp1.kpts_shift=[0,0,0]

fw1.add_task(vasp1)
fw1.add_handler('FrozenJobErrorHandler')

fw1_task2=ScriptTask.from_str("pwd; mkdir OrbDir")
fw1_task3=ScriptTask.from_str("mv CHG* CONTCAR DOSCAR EIGENVAL I* K* OSZICAR OUTCAR P* W* X* vasp* OrbDir")
fw1_task4=ScriptTask.from_str("cp OrbDir/WAVECAR . ; cp OrbDir/CONTCAR .")
fw1.add_task(fw1_task2)
fw1.add_task(fw1_task3)
fw1.add_task(fw1_task4)

# create VASP IPA LOPTICS input job yaml file
vasp2=VaspInput('CONTCAR')
vasp2.ISTART=1
vasp2.ALGO='Normal'
vasp2.LHFCALC='.TRUE.'
vasp2.LOPTICS='True'
vasp2.CSHIFT=0.1
vasp2.ISTART=1
vasp2.HFSCREEN=0.2
vasp2.NEDOS=4096
vasp2.KPAR=kpar_dft
vasp2.NBANDS=nbands
vasp2.kpts=[6,6,6]
vasp2.kpts_shift=[0,0,0]
vasp2.kpts_style='Gamma'
vasp2.TIME="0.4"
vasp2.PRECFOCK="Fast"
vasp2.NKRED=3

fw1.add_task(vasp2)

# Create VASP LFE Optics input.  This is for a new FireWorker
vasp3=VaspInput('CONTCAR')
vasp3.ALGO='Chi'
vasp3.LHFCALC='.TRUE.'
vasp3.HFSCREEN=0.2
vasp3.LSPECTRAL='False'
vasp3.NEDOS=4096
vasp3.NOMEGA=1024
vasp3.NBANDS=nbands
vasp3.LOPTICS='False'
vasp3.LRPA='.FALSE.'
vasp3.LFXC='True'
vasp3.NKREDLFX=3
vasp3.NKREDLFY=3
vasp3.NKREDLFZ=3
vasp3.KPAR=kpar_chi
vasp3.kpts_style='Gamma'
vasp3.kpts=[6,6,6]
vasp3.kpts_shift=[0,0,0]

fw2=VaspFirework(vasp3, name=mname)
fw2.copy_files_from_previous('WAVECAR', 'WAVEDER', mode='copy', ignore_errors=True)


# Make Sequential Workflow
#wf1 = VASPWorkflow(fw1,fw2,fw3, name=material_file_name.split('.')[0])
wf=VaspWorkflow(fw1,fw2,name=material_file_name.split('.')[0])
wf.to_file(material_file_name)

# save specification to yaml file for later inspection
# or manual add to lauchpad with lpad script
#wf.add_wf_to_launchpad()

