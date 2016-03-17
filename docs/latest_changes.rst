Change log
==========

v3.3.5
------
* StructureMatcher can now work with ignored species.
* Added interpolation failure warnings and smooth tolerance for
  scipy.interpolate.splrep in bandstructures (Tess).
* Added DiffusionAnalyzer.get_framework_rms_plot.
* Complete rewrite of Procar class to use ND array access and zero-based
  indexing.
* OrderParameters class for analysis of local structural features
  (Nils Zimmermann).
* Bug fixes for Procar, MPRester and SpaceGroup 64.
* Added Github templates for contributing to pymatgen.
