Change log
==========

v2023.1.20
----------
* Passthrough kwargs support for Structure.from_file and Structure.from_str
* Allow the `frac_tolerance` to be specified for rounding coordinates in CifParser.
* PR #2803 from @amkrajewski add_weightbasedfunctions
    When working with metallic alloys, weight-fraction-based notations such as Ti64 / Ti-6V-4Al or NiTiNOL60 / Ni-40Ti are commonly employed in both industrial specifications and scientific literature. Regardless of the numerous downsides of this situation, including errors in scientific experiments or NLP-parsing when they are mistaken for atomic fractions or chemical formulas, being able to create a Composition object from them (under correct interpretation) would be a useful pymatgen feature.
    - Composition class method to initialize it from a dictionary of weight fractions
    - Composition property giving a dictionary of weight fractions
    - concise tests for the two above were added
QChem: translate DMSO name in smd_solvent
