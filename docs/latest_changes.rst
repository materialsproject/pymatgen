Change log
==========

v3.0.13
-------
* Bug fix for parsing certain types of CIFs.
* MPRester now has get_materials_id_references helper method.
* Minor fix for Vasprun.final_energy.
* Added mp_decode option to MPRester.query to allow option to not decode into
  pymatgen objects.
* New POTCAR hash scheme to more robustly identify unique POTCARs.
* Link to http://bit.ly/materialsapi for information on Materials API
  document schema for use with MPRester.query method.
