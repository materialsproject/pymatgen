Change log
==========

v2022.1.20
----------
* Unicode fixes (@janosh)
* YAML deprecation fixes. (@janosh)
* ASE adaptor support for charge, spin multiiciplity and site properties of molecules. (@arosen93).
* New keyword option (`keep_site_properties`) in various `structure.symmetry.analyzer` functions to keep the site properties on the sites after a transformation. (@arosen93)
* Bug fixes for Lobster module (@JaGeo). 
* SCAN / GGA(+U) mixing scheme (@rkingsbury). Mixing scheme code lives in the new file `mixing_scheme.py` and is implemented as a `Compatibility` class.
* Fix for parsing of QuantumExpresso files due to new format (@vorwerkc)
