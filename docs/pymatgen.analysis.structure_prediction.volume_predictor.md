---
layout: default
title: pymatgen.analysis.structure_prediction.volume_predictor.md
nav_exclude: true
---

# pymatgen.analysis.structure_prediction.volume_predictor module

Predict volumes of crystal structures.


### _class_ pymatgen.analysis.structure_prediction.volume_predictor.DLSVolumePredictor(cutoff=4.0, min_scaling=0.5, max_scaling=1.5)
Bases: `object`

Data-mined lattice scaling (DLS) scheme that relies on data-mined bond
lengths to predict the crystal volume of a given structure.

As of 2/12/19, we suggest this method be used in conjunction with
min_scaling and max_scaling to prevent instances of very large, unphysical
predicted volumes found in a small subset of structures.


* **Parameters**


    * **cutoff** (*float*) – cutoff radius added to site radius for finding
    site pairs. Necessary to increase only if your initial
    structure guess is extremely bad (atoms way too far apart). In
    all other instances, increasing cutoff gives same answer
    but takes more time.


    * **min_scaling** (*float*) – if not None, this will ensure that the new
    volume is at least this fraction of the original (preventing
    too-small volumes)


    * **max_scaling** (*float*) – if not None, this will ensure that the new
    volume is at most this fraction of the original (preventing
    too-large volumes).



#### get_predicted_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), icsd_vol=False)
Given a structure, returns back the structure scaled to predicted
volume.


* **Parameters**

    **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure w/unknown volume



* **Returns**

    a Structure object with predicted volume



#### predict(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), icsd_vol=False)
Given a structure, returns the predicted volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – a crystal structure with an unknown volume.


    * **icsd_vol** (*bool*) – True if the input structure’s volume comes from
    ICSD.



* **Returns**

    a float value of the predicted volume.



### _class_ pymatgen.analysis.structure_prediction.volume_predictor.RLSVolumePredictor(check_isostructural=True, radii_type='ionic-atomic', use_bv=True)
Bases: `object`

Reference lattice scaling (RLS) scheme that predicts the volume of a
structure based on a known crystal structure.


* **Parameters**


    * **check_isostructural** – Whether to test that the two structures are
    isostructural. This algo works best for isostructural compounds.
    Defaults to True.


    * **radii_type** (*str*) – Types of radii to use. You can specify “ionic”
    (only uses ionic radii), “atomic” (only uses atomic radii) or
    “ionic-atomic” (uses either ionic or atomic radii, with a
    preference for ionic where possible).


    * **use_bv** (*bool*) – Whether to use BVAnalyzer to determine oxidation
    states if not present.



#### get_predicted_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), ref_structure)
Given a structure, returns back the structure scaled to predicted
volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure w/unknown volume


    * **ref_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – A reference structure with a similar
    structure but different species.



* **Returns**

    a Structure object with predicted volume



#### predict(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), ref_structure)
Given a structure, returns the predicted volume.


* **Parameters**


    * **structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – structure w/unknown volume


    * **ref_structure** ([*Structure*](pymatgen.core.structure.md#pymatgen.core.structure.Structure)) – A reference structure with a similar
    structure but different species.



* **Returns**

    a float value of the predicted volume