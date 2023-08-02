---
layout: default
title: pymatgen.analysis.fragmenter.md
nav_exclude: true
---

# pymatgen.analysis.fragmenter module

Perform fragmentation of molecules.


### _class_ pymatgen.analysis.fragmenter.Fragmenter(molecule, edges=None, depth=1, open_rings=False, use_metal_edge_extender=False, opt_steps=10000, prev_unique_frag_dict=None, assume_previous_thoroughness=True)
Bases: `MSONable`

Molecule fragmenter class.

Standard constructor for molecule fragmentation.


* **Parameters**


    * **molecule** ([*Molecule*](pymatgen.core.structure.md#pymatgen.core.structure.Molecule)) – The molecule to fragment.


    * **edges** (*list*) – List of index pairs that define graph edges, aka molecule bonds. If not set,
    edges will be determined with OpenBabel. Defaults to None.


    * **depth** (*int*) – The number of levels of iterative fragmentation to perform, where each level
    will include fragments obtained by breaking one bond of a fragment one level up.
    Defaults to 1. However, if set to 0, instead all possible fragments are generated
    using an alternative, non-iterative scheme.


    * **open_rings** (*bool*) – Whether or not to open any rings encountered during fragmentation.
    Defaults to False. If true, any bond that fails to yield disconnected graphs when
    broken is instead removed and the entire structure is optimized with OpenBabel in
    order to obtain a good initial guess for an opened geometry that can then be put
    back into QChem to be optimized without the ring just reforming.


    * **use_metal_edge_extender** (*bool*) – Whether or not to attempt to add additional edges from
    O, N, F, or Cl to any Li or Mg atoms present that OpenBabel may have missed. Defaults
    to False. Most important for ionic bonding. Note that additional metal edges may yield
    new “rings” (e.g. -C-O-Li-O- in LiEC) that will not play nicely with ring opening.


    * **opt_steps** (*int*) – Number of optimization steps when opening rings. Defaults to 10000.


    * **prev_unique_frag_dict** (*dict*) – A dictionary of previously identified unique fragments.
    Defaults to None. Typically only used when trying to find the set of unique fragments
    that come from multiple molecules.


    * **assume_previous_thoroughness** (*bool*) – Whether or not to assume that a molecule / fragment
    provided in prev_unique_frag_dict has all of its unique subfragments also provided in
    prev_unique_frag_dict. Defaults to True. This is an essential optimization when trying
    to find the set of unique fragments that come from multiple molecules if all of those
    molecules are being fully iteratively fragmented. However, if you’re passing a
    prev_unique_frag_dict which includes a molecule and its fragments that were generated
    at insufficient depth to find all possible subfragments to a fragmentation calculation
    of a different molecule that you aim to find all possible subfragments of and which has
    common subfragments with the previous molecule, this optimization will cause you to
    miss some unique subfragments.



### pymatgen.analysis.fragmenter.open_ring(mol_graph, bond, opt_steps)
Function to actually open a ring using OpenBabel’s local opt. Given a molecule
graph and a bond, convert the molecule graph into an OpenBabel molecule, remove
the given bond, perform the local opt with the number of steps determined by
self.steps, and then convert the resulting structure back into a molecule graph
to be returned.