---
layout: default
title: pymatgen.analysis.magnetism.jahnteller.md
nav_exclude: true
---

# pymatgen.analysis.magnetism.jahnteller module

JahnTeller distortion analysis.


### _class_ pymatgen.analysis.magnetism.jahnteller.JahnTellerAnalyzer()
Bases: `object`

Will attempt to classify if structure *may* be Jahn-Teller active.
Class currently uses datafile of hard-coded common Jahn-Teller
active ions.
If structure is annotated with magnetic moments, will estimate
if structure may be high-spin or low-spin.
Class aims for more false-positives than false-negatives.

Init for JahnTellerAnalyzer.


#### get_analysis(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.

Obtain an analysis of a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    analysis of structure, with key ‘strength’ which may be ‘none’, ‘strong’,
    ‘weak’, or ‘unknown’ (Default value = 0.1)



#### get_analysis_and_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Obtain an analysis of a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    analysis of structure, with key ‘strength’ which may be ‘none’, ‘strong’,
    ‘weak’, or ‘unknown’ (Default value = 0.1) and decorated structure



#### get_magnitude_of_effect_from_species(species: str | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species), spin_state: str, motif: str)
Get magnitude of Jahn-Teller effect from provided species, spin state and motif.


* **Parameters**


    * **species** – e.g. Fe2+


    * **spin_state** – “high” or “low”


    * **motif** – “oct” or “tet”


Returns: “none”, “weak” or “strong


#### _static_ get_magnitude_of_effect_from_spin_config(motif: str, spin_config: dict[str, float])
Roughly, the magnitude of Jahn-Teller distortion will be:
\* in octahedral environments, strong if e_g orbitals
unevenly occupied but weak if t_2g orbitals unevenly
occupied
\* in tetrahedral environments always weaker.


* **Parameters**


    * **motif** – “oct” or “tet”


    * **spin_config** – dict of ‘e’ (e_g) and ‘t’ (t2_g)
    with number of electrons in each state


Returns:  “none”, “weak” or “strong”


#### is_jahn_teller_active(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.
Check if a given structure and if it may be Jahn-Teller
active or not. This is a heuristic, and may give false positives and
false negatives (false positives are preferred).


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    boolean, True if might be Jahn-Teller active, False if not



#### _static_ mu_so(species: str | [Species](pymatgen.core.periodic_table.md#pymatgen.core.periodic_table.Species), motif: Literal['oct', 'tet'], spin_state: Literal['high', 'low'])
Calculates the spin-only magnetic moment for a
given species. Only supports transition metals.


* **Parameters**


    * **species** – Species


    * **motif** (*"oct"** | **"tet"*) – Tetrahedron or octahedron crystal site coordination


    * **spin_state** (*"low"** | **"high"*) – Whether the species is in a high or low spin state



* **Returns**

    Spin-only magnetic moment in Bohr magnetons or None if

        species crystal field not defined




* **Return type**

    float



#### tag_structure(structure: [Structure](pymatgen.core.structure.md#pymatgen.core.structure.Structure), calculate_valences: bool = True, guesstimate_spin: bool = False, op_threshold: float = 0.1)
Convenience method, uses get_analysis_and_structure method.
Add a “possible_jt_active” site property on Structure.


* **Parameters**


    * **structure** – input structure


    * **calculate_valences** – whether to attempt to calculate valences or not, structure
    should have oxidation states to perform analysis (Default value = True)


    * **guesstimate_spin** – whether to guesstimate spin state from magnetic moments
    or not, use with caution (Default value = False)


    * **op_threshold** – threshold for order parameter above which to consider site
    to match an octahedral or tetrahedral motif, since Jahn-Teller structures
    can often be
    quite distorted, this threshold is smaller than one might expect



* **Returns**

    Decorated Structure, will be in primitive setting.