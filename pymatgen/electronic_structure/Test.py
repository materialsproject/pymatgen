#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pymatgen.io.lobster import Doscar
from pymatgen.io.lobster import Lobsterout, Lobsterin
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.electronic_structure.core import Orbital, OrbitalType, Spin
from collections import namedtuple
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.electronic_structure.dos import Dos
from pymatgen.core.periodic_table import Element
import numpy as np
import pandas as pd
import os
import re
import plotly.graph_objects as go
import ast
import warnings
# from sklearn import preprocessing
from scipy.integrate import trapezoid
# from sklearn.linear_model import LinearRegression


# In[2]:


import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpl_patches

sns.set_style("white")
sns.set_context("talk")
sns.set_palette(["#0CB1F3", "#F34E0C"])

# In[3]:


os.chdir('/hpc-user/AG-JGeorge/anaik/')

# In[4]:


parent = os.getcwd()
os.chdir('Phonon_Dataset/Results/')

# In[5]:


mpids_lob = [f for f in os.listdir() if not f.startswith('t') and not f.startswith('.') and not f.startswith('__')
             and os.path.isdir(f)]
mats = list(set([ids.split('_')[0] for ids in mpids_lob]))
mats.sort()

# In[6]:


mpid = mats[100]
doscar_lobster = Doscar(doscar="{}/DOSCAR.lobster.gz".format(mpid),
                        structure_file="{}/POSCAR.gz".format(mpid),
                        dftprogram="Vasp")

dos_lobster = doscar_lobster.completedos

vasprun = Vasprun("{}/vasprun.xml.gz".format(mpid))
# Sys_elec= round(Outcar('{}/OUTCAR.gz'.format(mpid)).nelect,2)
# Lobout = Lobsterout('{}/lobsterout.gz'.format(mpid))

dos_vasp = vasprun.complete_dos



fp_lobster = dos_lobster.get_dos_fp(type='tdos', n_bins=100, normalize=True, binning=True)


fp_vasp = dos_vasp.get_dos_fp(type='s', max_e=0, min_e=-15, n_bins=100, normalize=True)


dos_vasp.get_dos_fp_similarity(fp_vasp, fp_lobster, tanimoto=False, normalize=True)