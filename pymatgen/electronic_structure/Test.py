#!/usr/bin/env python

# In[1]:


import os

import seaborn as sns

from pymatgen.io.lobster import Doscar
from pymatgen.io.vasp.outputs import Vasprun

# from sklearn import preprocessing
# from sklearn.linear_model import LinearRegression


# In[2]:


sns.set_style("white")
sns.set_context("talk")
sns.set_palette(["#0CB1F3", "#F34E0C"])

# In[3]:


os.chdir("/hpc-user/AG-JGeorge/anaik/")

# In[4]:


parent = os.getcwd()
os.chdir("Phonon_Dataset/Results/")

# In[5]:


mpids_lob = [
    f
    for f in os.listdir()
    if not f.startswith("t") and not f.startswith(".") and not f.startswith("__") and os.path.isdir(f)
]
mats = list({ids.split("_")[0] for ids in mpids_lob})
mats.sort()

# In[6]:


mpid = mats[100]
doscar_lobster = Doscar(doscar=f"{mpid}/DOSCAR.lobster.gz", structure_file=f"{mpid}/POSCAR.gz", dftprogram="Vasp")

dos_lobster = doscar_lobster.completedos

vasprun = Vasprun(f"{mpid}/vasprun.xml.gz")
# Sys_elec= round(Outcar('{}/OUTCAR.gz'.format(mpid)).nelect,2)
# Lobout = Lobsterout('{}/lobsterout.gz'.format(mpid))

dos_vasp = vasprun.complete_dos


fp_lobster = dos_lobster.get_dos_fp(type="tdos", n_bins=100, normalize=True, binning=True)


fp_vasp = dos_vasp.get_dos_fp(type="s", max_e=0, min_e=-15, n_bins=100, normalize=True)


dos_vasp.get_dos_fp_similarity(fp_vasp, fp_lobster, tanimoto=False, normalize=True)
