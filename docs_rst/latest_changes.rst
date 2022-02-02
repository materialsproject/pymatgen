Change log
==========

v2022.2.1
---------
* Chargemol caller for partial atomic charge analysis (@arosen93)
* ASEAtomAdaptor: (1) Updates to magmom support, (2) Oxidation states support, (3) Charges are now passed (@arosen93)
* Cleanup of deprecated methods. (@janosh) github-actions-setup-python-now-supports-dependency-caching
* Bigfix for gzipped DOSCAR (@JaGeo)
* Updates for QChem Support (@samblau)
* QuantumEspresso k-grid fix input fix. (@vorwerkc)
* `Entry.__repr__()` now ouputs name where available. (@janosh)
* Fixes to Vasprun.final_energy to report `e_0_energy` (the desired energy quantity) for VASP 6+. (@arosen93) 
* `Outcar().final_energy` now prints out `e_0_energy` (also called "energy(sigma->0)" in the OUTCAR) rather than `energy_fr_energy` (also called "free  energy   TOTEN" in the OUTCAR). This is to be consistent with `Vasprun().final_energy` and because it is generally the desired quantity. `Outcar` now has two new attributes: `.final_energy_wo_entrp` and `final_fr_energy`, which correspond to `e_wo_entrp` and `e_fr_energy`, respectively. (@arosen93)
* Improved parsing of coupled cluster calculations in QChem (@espottesmith).
