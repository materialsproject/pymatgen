Change log
==========

v2022.0.15
----------
Welcome to new contributors @blokhin, @pzarabadip, @ml-evs, @wuxiaohua1011, @janssenhenning and @penicillin0. A reminder to all new contributors to 
ensure your information is accurate at https://pymatgen.org/team.html so that 
you are acknowledged appropriately by filling out the linked form.

* Breaking change in PhaseDiagram serialization which will affect any users of BasePhaseDiagram which has now been removed (@shyuep, 2b9911d)

* Speed up nearest-neighbor routines & structure graph generation (@ltalirz, #2239)
* Add two more pre-defined OPTIMADE aliases (@blokhin, #2242)
* Refactor `interface_reactions` module, adding support for Plotly (@mattmcdermott, #2233)

* Update NOMAD access in MPRester (@wuxiaohua1011, #1958)
* General improvements to Phase Diagram code (@CompyRhys, #2263, #2264, #2268)
* Improve appearance of periodic table heatmap (@penicillin0, #2272)
* Small improvements to battery classes (@jmmshn, #2262)
* Fix for Composition.chemical_system to match expected behaviour for compositions with oxidation states (@CompRhys, #2249)
* Fix for bad param in OPTIMADE reponse fields (@ml-evs, #2244)
* Fix for issue in parsing `bandOverlaps.lobster` file (@pzarabadip, #2237)
* Fix for Moladaptor (@orioncohen, #2269)
* Fix for incorrect Potcar hash warnings (@mkhorton, #2273)

* Type hint and correct documentation of Structure.remove_site_properties (@kmu, #2256)
* Type hint improvements across pymatgen (@janosh, #2241, #2247, #2261)
* Add `pymatgen-io-fleur` addon to addons page (@janssenhenning, #2232)

