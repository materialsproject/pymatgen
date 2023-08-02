---
layout: default
title: pymatgen.vis.structure_chemview.md
nav_exclude: true
---

# pymatgen.vis.structure_chemview module

Visualization for structures using chemview.


### pymatgen.vis.structure_chemview.quick_view(structure, bonds=True, conventional=False, transform=None, show_box=True, bond_tol=0.2, stick_radius=0.1)
A function to visualize pymatgen Structure objects in jupyter notebook using chemview package.


* **Parameters**


    * **structure** – pymatgen Structure


    * **bonds** – (bool) visualize bonds. Bonds are found by comparing distances
    to added covalent radii of pairs. Defaults to True.


    * **conventional** – (bool) use conventional cell. Defaults to False.


    * **transform** – (list) can be used to make supercells with pymatgen.Structure.make_supercell method


    * **show_box** – (bool) unit cell is shown. Defaults to True.


    * **bond_tol** – (float) used if bonds=True. Sets the extra distance tolerance when finding bonds.


    * **stick_radius** – (float) radius of bonds.



* **Returns**

    A chemview.MolecularViewer object