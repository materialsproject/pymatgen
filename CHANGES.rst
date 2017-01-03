Change log
==========

v4.5.7
------
* PMG settings are now prefixed with PMG_ to ensure proper namespacing.
* Improve error output in command line bader caller.
* Add Py3.6 classifier.
* Misc bug fixes.

v4.5.6
------
* Minor bug fix.
* Fixed elastic energy density

v4.5.5
------
* Fix bad reading of pmgrc.
* Gaussian opt section added allowing for torsion constraints
* Update spglib.

v4.5.4
------
* BSDOSPlotter (Anubhav Jain)
* Fixes to defect analysis (Bharat)
* intrans as an input to BoltztrapAnalyzer. Allows for scissor operation.
* Pmg is now continuously tested on win-64/py35 using Appveyor!

v4.5.3
------
* Added an alternative interstitial finder that works with a grid-based structure-motif search. (Nils Zimmermann)
* Optionnal possibility to specify that the saddle_point in the NEB should have a zero slope. (David Waroquiers)
* Read intensity and normal modes for Gaussian. (Germain Salvato Vallverdu)
* Minor bug fixes.

v4.5.2
------
* Minor bug fix for POTCAR settings.

v4.5.1
------
* You can now specify a different default functional choice for pymatgen by
  setting PMG_DEFAULT_FUNCTIONAL in .pmgrc.yaml. For use with newer
  functional sets, you need to specify PBE_52 or PBE_54 for example.
* Swtich to ISYM 3 by default for HSE.
* Updates to FEFF>
* Misc bug fixes and startup speed improvements.

v4.5.0
------
* Major speed up of initial load.
* Collection of misc changes.


v4.4.12
-------
* Fix for dynamic numpy import.

v4.4.11
-------
* Update to new version of spglib.

v4.4.10
-------
* Minor fixes for proper gzipped structure file support and MVLSlabSet.

v4.4.9
------
* Dependency cleanup. Now, basic pymatgen requires on much fewer
  packages.
* Fixed reading of POSCAR files when more than 20 types of atoms.
* Misc bug fixes.

v4.4.8
------
* Cleanup of entry points and dependencies.

v4.4.7
------
* Update to spglib 1.9.7.1
* Proper use of dependency markers for enum34.

v4.4.6
------
* Update to spglib 1.9.6, which fixes some bugs and is Windows compatible.

v4.4.5
------
* Bug fix for SubstitutionProb.

v4.4.4
------
* Bug fix for electronic structure plotter.

v4.4.3
------
* Bug fix for Diffusion Analyzer.

v4.4.2
------
* Bug fix for BS serialization.
* Cleanup dependencies.

v4.4.1
------
* Massive updates to FEFF support (Kiran Mathews).
* Bug fixes in band structure plotting.

v4.4.0
------
* Much more Pythonic API for modifying Structure/Molecule species. Now,
  strings, slices, and sequences should magically work, in addition to the
  previous API of simple int indices. Examples::

    s[0] = "Fe"
    s[0] = "Fe", [0.5, 0.5, 0.5]  # Replaces site and fractional coordinates.
    s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}  # Replaces site and fractional coordinates and properties.
    s[(0, 2, 3)] = "Fe"  # Replaces sites 0, 2 and 3 with Fe.
    s[0::2] = "Fe"  # Replaces all even index sites with Fe.
    s["Mn"] = "Fe"  # Replaces all Mn in the structure with Fe.
    s["Mn"] = "Fe0.5Co0.5"  # Replaces all Mn in the structure with Fe: 0.5, Co: 0.5, i.e.,creates a disordered structure!

* Massive update to internal representation of Bandstructure objects for
  memory and computational efficiency.
* Bug fixes to CIF parsing in some edge cases. (Will Richards).

v4.3.2
------
* Massive speedup of Bandstructure, especially projected band structures,
  parsing.
* Massive update to pmg cli script, with new query functions as well as a
  more rational command structure.
* Updates to ChemEnv.
* Misc bug fixes.

v4.3.1
------
* Upgrade monty and spglib requirements for bug fixes.
* Updates to feff support (Kiran).

v4.3.0
------
* Massive update to elastic module. (Joey Montaya)
* Pathfinder algorithm for NEB calculations. (Ziqing Rong)
* Wheels for Windows and Mac Py27 and Py35.

v4.2.5
------
* Bug fix for BSPlotter.

v4.2.4
------
* Bug fix for kpoint weight calculation for Monkhorst meshes.

v4.2.3
------
* Minor cleanup.
* Simplified installation. enumlib and bader can now be installed using pmg setup --install.

v4.2.2
------
* Global configuration variables such as VASP\_PSP\_DIR and MAPI\_KEY are now
  stored in "~/.pmgrc.yaml". If you are setting these as environmental
  variables right now, you can easily transition to the new system using::

      pmg config --add VASP_PSP_DIR $VASP_PSP_DIR MAPI_KEY $MAPI_KEY

  This new scheme will provide greater flexibility for user-defined
  global behavior in pymatgen, e.g., tolerances, default input sets for
  transmuters, etc., in future.
* Beta of k-point weight calculator.
* Use default MSONable as and from_dict for all transformations.

v4.2.1
------
* New DopingTransformation that implements an automated doping strategy.
* Updated MIC algorithm that is a lot more robust (Will Richards).
* Major update to chemenv package (David Waroquiers)

v4.2.0
------
* Fix important bug in minimum image distance computation for very skewed cells.
* Major refactoring of WulffShape code.
* Misc bug fixes for elastic tensor and other codes.

v4.1.1
------
* Major refactoring of WulffShape and lammps support.

v4.1.0
------
* Wulff shape generator and analysis.
* Minor bug fixes.

v4.0.2
--------
* Fix kpoint reciprocal density.

v4.0.1
------
* Minor bug fix release.

v4.0.0
------
* Massive update with many deprecated methods removed. Note that this
  may break backwards incompatibility!
* Support for ABINIT 8.
* Improved sulfide compatibility.

v3.7.1
------
* Fix deprecation bug.

v3.7.0
------
* Last version before pymatgen 4.0, where deprecated modules will be removed!
* Massive update to LAMMPS (Kiran Matthews).
* New input sets with a different interface that replaces old input sets.
* Massive update to elastic properties.

v3.6.1
------
* Massive cleanup to Boltztrap interface (Anubhav Jain)
* Refactor of piezoelectric analysis to use tensor base class (Joey)
* More robust CIF parsing.

v3.6.0
------
* Pymatgen now uses spglib directly from Togo's website. Spglib is no longer
  bundled as a dependency.
* Improved support for velocities in Poscar (Germaine Vallverdu)
* Backwards incompatible change in Born charge format in Outcar.
* Fixes for Lammps input serialization

v3.5.3
------
* Misc refactorings and bug fixes, especially for Outcar and Boltztrap classes.

v3.5.2
------
* Minor update to DerivedInputSet interface.

v3.5.1
------
* New derived input sets for generating inputs that depende on previuos
  calculations. Old input sets deprecated.

v3.5.0
------
* Chemical environment analysis package (David Waroquiers).
* Piezoelectric property analysis (Shayam).
* Cythonize certain expensive core functions. 5-10x speedup in large structure matching (Will Richards).
* New NMR parsing functionality for Outcar (Xiaohui Qu).
* Improved io.lammps (Kiran Mathews).
* Update to spglib 1.9.2.
* Element properties now return unitized float where possible.
* Bug fix for get_primitive_standard affecting rhombohedral cells (important for band structures).
* Vasprun.final_energy now returns corrected energy with warning if it is different from final electronic step.

v3.4.0
------
* 10-100x speed up to Structure copying and Site init, which means many
  functionality has seen signifcant speed improvement (e.g., structure
  matching).
* Convenience method Structure.matches now perform similarity matching
  for Structures.
* Bugfix for band gap determination.

v3.3.6
------
* Update to use enum.x instead of multienum.x.
* Minor robustness fixes to VaspInputSet serialization.
* Add a reciprocal density parameter to vasp sets.
* Minor bug fixes to Vasprun parsing. 

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

v3.3.4
------
* Procar now supports parsing of phase factors.
* Miscellaneous bug fixes.

v3.3.3
------
* Bug fixes for Poscar.
* Fix Kpoints pickling.

v3.3.2
------
* Bug fixes for pymatgen.io.abinit
* Other minor big fixes.

v3.3.1
------
* Minor bug fix release for pickle and elastic constants.

v3.3.0
------
* Updated and checked for Python 3.5.* compatibility.
* Element, Spin, Orbital and various other Enum-like classes are now actually
  implemented using Enum (with enum34 dependency for Python < 3.4).
* Speed up Site creation by 20% for ordered sites, with cost in terms of
  slightly slower non-ordered Sites. Since ordered Sites is the far more common
  case, this gives significant boost for large scale manipulations of
  structures.
* Alternative, more pythonic syntax for creating supercells via simply
  Structure * 3 or Structure * (3, 1, 1).
* zeo++ fixes.
* More stable incar settings for MITMDVaspInputSet.

v3.2.10
-------
* Fix missing scripts
* Improvements to units module.
* Speed up EwaldSummation.

v3.2.9
------
* Major PD stability improvements, especially for very high dim hulls with lots
  of entries.
* Improvements to Ewald summation to be close to GULP implementation.
* Deprecate physical constants module in favor of scipy's version.
* Remove many pyhull references to use scipy's ConvexHull implementation.
* Bug fix for sulfide correction.

v3.2.8
------

* Make pyhull optional.
* Sulfur correction added to MaterialsProjectCompatibility for more accurate
  sulfide formation energies.
* ADF io support. (Xin Chen)
* Bug fixes for spacegroup subgroup testing.

v3.2.7
------
* Add warning for limited subgroup testing functionality in Spacegroup.

v3.2.6
------
* Extensive support for elasticity tensor analysis (Joseph Montoya).
* Misc bug fixes and performance improvements.
* Add support for QChem4.3 new format of Batch jobs

v3.2.5
------
* Improved potcar setup via "pmg setup", with MAPI setup.
* Support for new POTCARs issued by VASP.
* Improvements to ABINIT support.
* Improvement to Boltztrap support, e.g., scissor band gap, etc.
* Vasprun now issues warning when unconverged run is detected.

v3.2.4
------

* GaussianOutput can now parse frequencies, normal modes and cartesian forces
  (Xin Chen).
* Support for Aiida<->pymatgen conversion by the Aiida development team (Andrius
  Merkys).
* Specialized BSVasprun parser that is ~2-3x faster than Vasprun.
* Refactor the boltztrap package (merge a few methods together) and add several
  new methods (power factor, seebeck...)
* Support of the new PCM format in QChem 4.3
* Local environment analysis to pmg script.
* Deprecate prettytable in favor of tabulate package.
* Improvements to MITNEBVaspInputSet.
* Misc bug fixes.

v3.2.3
------
* Massive update to abinit support. Note that pymatgen.io.abinitio has 
  been refactored to pymatgen.io.abinit. (Matteo, Setten)
* NwOutput now supports parsing of Hessian matrices (contributed by Xin 
  Chen)
* Gaussian support now has the ability to read potential energy surface
  and electronic transitions computed with TD-DFT (Germain Salvato 
  Vallverdu)
* Bug fixes for CifWriter with symmetry.
* Bug fixes for surface generation and reactions.
* Monty requirement increased.

v3.2.1
------
* Fix wrong U value for Ce and Eu.
* Properly handle empty multiline strings in Cif
* Add ability to get specific data in MPRester.get_entries. Make all get_entry 
  methods consistent  in kwargs.

v3.2.0
------
* Force conversion to an actual list in selective dynamics and velocities in
  Poscar.
* fix small bug in BSPlotter (wrong ylim)
* Elastic tensor parsing in Outcar

v3.1.9
------
* Fix scripts.

v3.1.7
------
* Bug fixes for MPRester.
* Ensure correct monty version requirement in setup.py.

v3.1.6
------
* Rudimentary PWSCF output reading.
* Fix ASE support.
* Support for WAVEDERF and reading multiple dielectricfunctions in vasprun.xml. 
  (Miguel Dias Costa)

v3.1.5
------
* Move vasp.vasp*put to vasp.*puts. Also, maintain backwards compatibility with 
  vaspio.vasp_*put

v3.1.4
------
* Fix missing yaml files that have been moved.

v3.1.3
------
* Major refactoring of pymatgen.io. Now, the io suffix is dropped from all io
  classes. i.e., it is just pymatgen.io.vasp, not pymatgen.io.vaspio. Also, all
  input sets have been moved within the relevant package, e.g.,
  pymatgen.io.vasp.sets. All changes are backwards compatible for now. But
  deprecation messages have been included which states that the stubs will be
  removed in pymatgen 4.0. Pls migrate code when you see the deprecation
  messages.
* Make Composition.anonymized_formula truly chemistry independent (No A2B2
  for peroxides or A2 for diatomic gasses) 
* Allowing CIF data_* header to be prefixed with spaces and tabulations.

v3.1.2
------
* HHI Resource Analysis (by Anubhav Jain).
* Bug fixes for surfaces normalizatino.
* Bug fix for Vasprun parsing of response function keys.
* Dockerfile for generation of an image for pymatgen.
* Updated requirements.txt for latest requests, scipy, numpy.

v3.1.1
------
* Bug fixes for SpacegroupAnalyzer and SlabGenerator.
* Much faster normal vec search.

v3.1.0
------
* Much improved surface generation algorithm that provides for
  orthogonality constraints.
* Transition state analysis tools! (beta)
* Massive improvements in Outcar parsing which provides a powerful grepping
  syntax.
* PWSCFInput generation (beta).
* Reduce default SIGMA to 0.05 for MP input sets.
* Update spglib to 1.7.3 as per recommendation of Togo.
* Many bug fixes and efficiency improvements.

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

v3.0.11
-------
* Lots of abinitio improvements (Matteo).
* Added mp_decode option to MPRester.query to allow option to not decode into pymatgen objects.

v3.0.10
------

* Fix cartesian coord parsing in Poscar class.
* Vasprun now works with non-GGA PBE runs
* Misc bug fixes

v3.0.9
------
* Major bug fixes for CIF parsing (Will Richards).
* Support for {Li,Na} syntax in parse_criteria for MPRester.
* Additional example notebook for ordering and enumeration. 
* More robust checking for oxidation states in EnumerateStructureTRansformation.
* Improvements to Slab polarity checking.

v3.0.8
------
* Massive update to abinitio (Matteo).
* Improvements to OUTCAR parsing (Ioannis Petousis).

v3.0.7
------
* Powerful Slab generation algorithms (beta!).
* Improvements to DiffusionAnalyzer with constant smoothing option.
* Significantly improve look of DOS plots using prettyplotlib.

v3.0.6
------
* Cost analysis module (Anubhav Jain)
* More Py3k fixes.
* Extensive abinitio updates (Matteo).

v3.0.5
------
* Completely revamped symmetry package. The finder.SymmetryFinder and
  pointgroup and spacegroup modules are now deprecated. Instead,
  all symmetry analysis is in the :module:`pymatgen.symmetry.analyzer`_
  module. There is also a completely rewritten support for symmetry groups in
  :module:`pymatgen.symmetry.groups`_. Structure now supports a static
  constructor to generate a structure from a spacegroup (see examples).
* BatteryAnalyzer class (Anubhav Jain) to provide for some common analysis of
  intercalation electrodes.
* Minor bug fixes for structure_matcher, lattice, abinitio.
* MOAB qadapter for abinit. (Liam Damewood)

v3.0.4
------
* Fix missing structures json data.

v3.0.3
------
* Updates to DiffusionAnalyzer for more fine-grained statistics.
* Bug fixes and tweaks to linear assignment
* Improved PymatgenTest class which provides a suite of test structures.
* Speedups to Phase Diagram
* Lots of improvements to Gaussian support (Nicolas Dardenne) and Abinit IO
  (Matteo).
* Lots of Py3k minor updates.
* Updated doc for Diffusion anaylzer. Invert sq_disp_ions for more intuitive handling.

v3.0.2
------
1. Consistent use of unicode throughout pymatgen.
2. Minor bug fixes.

v3.0.1
------
1. Minor bug fixes for cifio.
2. Py3k updates for abinitio.

v3.0.0
------
* Pymatgen is now completely Python 2.7 and Python 3.x compatible!
* Spglib and pyhull have been updated to support Python 3.x.
* Completely rewritten pure python cifio module (courtesy of William Davidson
  Richards) removed dependency on PyCIFRW, which has been causing many issues
  with installation.
* Structure and Molecule now supports a very convenient to() and from_str and
  from_file functionality. Instead of trying to load the appropriate parser,
  you can output and read from the appropriate formats directly. See example
  usage.
* ~50% speedup to LinearAssignment code.
* Continuous integration and contribution guidelines now include Python 3.
* **Backwards incompatible changes**
* matgenie.py has now been renamed simply "pmg" for brevity.
* All deprecated methods in pymatgen 2.x have been removed. E.g.,
  pymatgen.core.structure_modifier is no longer available.
* Pymatgen classes now uses the as_dict() method protocol implemented in the
  Monty package. The to_dict property is now deprecated and will be removed
  in pymatgen v3.1.
* Update main docs page examples with the new Structure to, from formats.

v2.10.6
-------
* Bug fix for np1.9 incompatibility. Now works.
* Use wheel for pymatgen deployments.
* matgenie.py is now renamed to pmg for faster CLI usage.
* Improvements to KPOINTS automatic generation.
* Simpler and faster Structure.get_all_neighbors

v2.10.5
-------
* DiffusionAnalyzer now has non-smoothed option.
* Kpoints generation algorithm now guarantees minimum # of points.
* Compatibility now has a proper explanation dict.
* Vaspruns with NSW == 1 now checked properly for electronic conv.
* make_movie now supports kwargs.

v2.10.3
-------
* MPRester.query now supports a simple but powerful string criteria syntax
  with support for wild cards.
* Improvements to Composition - support for negative compositions, sorting etc.
* Speed ups to StructureMatcher.

v2.10.2
-------
* Bug fix for Projected DOS parsing in new Vasprun.
* Compatibility now has a *explain* method which provides a detailed outline
  of the changes that a Compatibility makes to an Entry.

v2.10.1
-------
* Minor fix for monty requirements in setup.py.

v2.10.0
-------
* Major update: MPRester now uses Materials API v2! Also major refactoring
  of MPRester.
* Vastly improved Vasprun parser using cElementTree. Twice as fast,
  half as much code and easier to maintain.
* Vast improvements to Qchem functionality (Xiaohui Qu).
* Improved handling of Structure manipulations for extremely large
  structures (particularly in terms of memory consumption).
* Bug fix for XYZ parsing for scientific notation.
* Improve monty.serialization for transparent handling of JSON vs YAML.
  Requirements updated to monty>=0.3.3.
* Update numpy requirements to 1.8+. Fixes memory leak.
* Other minor bug fixes.

v2.9.14
-------
* Implements Structure.sort method. Both Structure.sort and the
  get_sorted_structure methods now supports all arguments supported by list
  .sort().
* VaspInputSets configs, as well as several other configs now uses yaml. Note
  the new dependency on pyyaml. It is highly recommended that you install
  pyyaml with the libyaml C bindings.
* Fix missing spglib dependency.
* Use monty.serialization for transparent handling of JSON vs YAML.
  Requirements updated to monty>=0.3.1.

v2.9.13
-------
* Urgent bug fix for missing compatibility yamls.

v2.9.12
-------
* Defect transformations (Bharat).
* Support for optical properties (Geoffroy Hautier and David Waroquiers).
* Improved support for some VASP output files (XDATCAR and OSZICAR).
* Refactored compatibilities now uses YAML for ease of reading.

v2.9.11
-------
* Bug fix for get_xrd_plot.
* Speed up XRD calculator by allowing specification of two theta ranges.
* Minor improvements to Gulp caller.

v2.9.10
-------
* Bug fix for unequal coefficients sizes in XRD.
* Support for Ag radiation in XRD calculator.
* Improved Procar class for extraction of information. (Germain Salvato
  Vallverdu)
* Bug fix for extraction of GGA data from Materials API.

v2.9.9
------
* XRDCalculator now supports disordered structures.
* Minor speed ups and improvements.

v2.9.8
------
* Initial beta version of XRD pattern calculator.
* Pymatgen now uses spglib 1.6.0.
* Update to Vasprun to compute static deilectric constants with DFPT in VASP.
  (Geoffroy Hautier)

v2.9.7
------
* Quick bug-fix release that provides a better solution to Structure handling
  of properties instead of sanitizing MPRester structures.

v2.9.6
------
* Patch to allow 1D phase diagrams (essentially finding the lowest energy
  phase).
* Better error checking for Bandstructure KPOINTs.
* Patch to sanitize structures obtained from MPRester.

v2.9.5
------
* Bug fix for linear assignment, which may sometimes affect Structure
  Matcher results.
* Minor improvement to the way grand canonical PDs work.

v2.9.4
------
* Bug fix for Pourbaix Maker (Sai).
* Streamline use of scratch directories for various calls. Require monty >=
  0.1.2.
* High accuracy mode for Zeo++ (Bharat Medasani).

v2.9.3
------
* Bug fix release for printing TransformedStructures from Substitutor (Will
  Richards).
* Misc improvements in BVAnalyzer, coord_utils and defects (Will Richards,
  David Waroquiers and Bharat Medasani).

v2.9.2
------
* Bug fix release for DummySpecie, which failed when deserializing from
  json and had bad hash function.

v2.9.1
------
* Structure/Molecule now supports Pythonic list-like API for replacing and
  removing sites. See :ref:`quick_start` for examples.

v2.9.0
------
* Updates to support ABINIT 7.6.1 (by Matteo Giantomassi).
* Vastly improved docs.
* Refactoring to move commonly used Python utility functions to `Monty
  package <https://pypi.python.org/pypi/monty>`_, which is now a dependency
  for pymatgen.
* Minor fixes and improvements to DiffusionAnalyzer.
* Many bug fixes and improvements.

v2.8.10
-------
* Refactoring of qchemio module (by Xiaohui Qu).

v2.8.9
------
* qchemio module (by Xiaohui Qu).

v2.8.8
------
* Minor bug fix release for Structure species substitution methods.

v2.8.7
------
* Massive update to pymatgen.io.abinitio package (by Matteo Giantomassi).
* Bug fixes for StructureMatcher's group_structure.
* Misc bug fixes and cleanup.

v2.8.6
------
* Bug fix for VASP io set introduced by the default to sorting of structure
  sites when generating VASP input.

v2.8.4
------
* Completely revamped Compatibility/Correction system which improves
  readability (Shyue Ping Ong/Anubhav Jain/Sai Jayaraman). This change is
  backwards compatible for the most part.

v2.8.3
------
* Big fix release for json dumping for unitized floats.

v2.8.2
------
* Bug fix release to improve CIF parsing for more non-standard CIF files.
  In particular, non-ascii characters are removed and _cgraph* fields are
  removed prior to parsing for better support in PyCiFRW.

v2.8.1
------
* Bug fix release. Incorrect units assigned for ionic radii.
* Improved nwchemio supports COSMO and ESP calculations (Nav Rajput).

v2.8.0
------
* **Units**. Pymatgen now has a new system of managing units,
  defined in pymatgen.core.units. Typical energy, length, time,
  temperature and charge units are supported. Units subclass float,
  which makes the usage transparent in all functions. The value that they
  being are in terms of conversions between different units and addition and
  subtraction of different units of the same type. Some basic quantities
  like ionic radii and atomic masses are now returned in unitized forms for
  easy conversion. Please see :mod:`pymatgen.core.units` and the
  :doc:`examples </examples>` for a demonstration of house to use units in
  pymatgen.
* **Minor backwards-incompatible change**. Structures are now sorted by
  default when generating VASP input files using vaspio_set. Old behavior can
  be obtained by setting sort_structure=False in the constructor. This is
  typically the desired behavior and prevents the generation of large
  POTCARs when atomic species are not grouped together.
* Bug fix for Molecule.substitute. Earlier algorithm was not detecting
  terminal atoms properly.
* Additional conversion tools for ABINIT (by Matteo Giantomassi).

v2.7.9
------
* Minor bug fix release to fix pyhull dependencies to be more friendly.
* Improved structure matcher that allows for more flexible matching. New
  matching between ordered and disordered comparator.

v2.7.7
-------
* Beta new Gulp Caller and Zeo++ interface classes (Bharat . Zeo++ is an open
  source software for performing high-throughput geometry-based analysis of
  porous materials and their voids. Please see
  http://www.maciejharanczyk.info/Zeopp/about.html.
* Specify version of distribute to 0.6.34 for better compatibility.

v2.7.6
------
* Support for VTK 6.x in structure visualization.
* Updated install instructions for openbabel.
* Preliminary pourbaix analysis (Sai Jayaratnam).

v2.7.5
------
* Vastly improved Nwchem IO (by Shyue Ping Ong).
* Much improved ABINIT support (by Matteo Giantomassi).

v2.7.4
------
* Added basic Nwchem (http://www.nwchem-sw.org/) IO support. (by: Shyue Ping
  Ong).
* New MoleculeMatcher class for comparing molecules by RMS. Requires
  openbabel with python bindings. (by: Xiaohui Qu)
* New functional group substitution capability for molecules (by: Lei Cheng
  and Shyue Ping Ong).

v2.7.2
------
* Minor bug fix release to fix some rare errors in very high dimensional
  phase diagrams. **Requires new pyhull version (1.3.8).**

v2.7.1
------
* **Major backwards-incompatible change.** With effect from v2.7.1,
  the default Structure and Molecule classes are now *mutable* objects. All
  functionality in the :mod:`pymatgen.core.structure_modifier` has been
  ported over to the new mutable classes. This change was implemented
  because the immutability of Structure and Molecule has resulted in very
  awkward code to make changes to them. The main cost of this change is that
  Structure and Molecule can no longer be used as dict keys (__hash__ has
  been set to None). However, we believe this is a minor cost given that we
  have rarely seen the use of Structure or Molecule as dict keys in any case.
  For the rare instances where such functionality is needed,
  we have provided the IStructure and IMolecule classes (where I indicates
  immutability) which will perform exactly the same way as the previous
  classes. With this change, the :mod:`pymatgen.core.structure_modifier`
  module is now deprecated and will be removed in a future version.
* read_structure and write_structure now supports pymatgen's json serialized
  structures.
* read_mol and write_mol functions now available (analogues of
  read_structure and write_structure for molecules)

v2.7.0
------
* Beta support for ABINIT input and output via pymatgen.io.abinitio
  (courtesy of the excellent work of Matteo Giantomassi).
* Properties are now checked when comparing two Species for equality.
* MaterialsProjectVaspInputSet is now renamed to MPVaspInputSet for easier
  typing. The old input sets have been deprecated.
* New VaspInputSets for MPStatic, MPNonSCF, MITMD which supports uniform
  grid, bandstructure and molecular dynamics calculations. The MD input set
  uses MIT parameters for speed.
* A beta DiffusionAnalysis class in the apps package.
* A revised KPOINT grid algorithm that generates more reasonable meshes.
* A guided install script is now provided for Mac and Linux users.

v2.6.6
------
* Updates to feffio (credit: Alan Dozier)
* Added detailed installation instructions for various platforms.
* Support for charge and spin multiplicity in Molecule. Expanded methods
  available in Molecule.
* Added supercell matching capabilities to StructureMatcher.
* More robust creation of PhaseDiagrams to take into account potential qhull
  precision errors.

v2.6.5
------
* Added a command_line caller to do Bader charge analysis using Henkelmann
  et al.'s algorithm.
* Bug fix for POSCAR parsing when title line is an empty string.
* Added __rmul__ operator for Composition.
* Vastly expanded available aliases.

v2.6.4
------
* Bug fixes for selective dynamics in Poscar.
* Improved Procar parsing to support both simple and detailed PROCARs.

v2.6.3
------
* Added new MaterialsProject REST interfaces for submit/query/delete_snl
  (currently open in beta for collaborators only).
* Added support for new MaterialsProject REST method get_stability.
* Added aliases for PhaseDiagram, GrandPotentialPhaseDiagram,
  PDAnalyzer and PDPlotter in pymatgen.phasediagrams.
* Improvements to StructureMatcher: stol (site - tolerance) redefined as
  a fraction of the average length per atom. Structures matched in fractional
  space are now also matched in cartesian space and a rms displacement
  normalized by length per atom can be returned using the rms_dist method.

v2.6.2
------

* Site and PeriodicSite now uses a Composition mapping type to represent
  the species and occupancy, instead of a standard dict.
* Bug fix for reading and re-writing out of Potcars.
* VaspInputSet now supports MSONable framework.
* Strain cell option in StructureEditor.
* Miscellaneous bug fixes and speedups.

v2.6.1
------
* Use requests.Session in MPRester for connection pooling and code simplicity.
* Support for "with" context manager in MPRester.
* Updated periodic table data to correct errors in Ru, Tc and other elements.
* New methods in Lattice to obtain Wigner-Seitz cell and Brillouin Zone.
* Miscellaneous bug fixes and speedups.

v2.5.5
------

* Bug fix release for cifio for rhombohedral structures.
* Miscellaneous bug fixes and speedups.

v2.5.4
------
* Vastly improved Gaussian input file parsing that supports more varieties
  of input specifications.
* StructureNL now supports molecules as well as structures.
* Updated atomic and vdw radius for Elements.
* Miscellaneous bug fixes and speedups.

v2.5.3
------
* Bug fix for StructureNotationalLanguage.
* Support for LDA US potential. matgenie.py script option to generate POTCARs.
* Beta version of StructureNotationLanguage, a markup format for Structure
  data with metadata such as authors and references. (Anubhav Jain)
* Vasprun parsing now parses dielectric constant where available. (Geoffroy
  Hautier)
* New custom ipython shell script for pymatgen.
* Miscellaneous bug fixes and speedups.

v2.5.1
------
* Bug fixes for primitive cell finder.
* Remove deprecated use_external_qhull option in PhaseDiagram classes.
* Miscellaneous bug fixes and speedups.

v2.5.0
------
* Added optimization package with linear assignment class.
* Improved robustness of StructureMatcher using linear assignment.
* Improved primitive cell search (faster and more robust).
* Cleanup of deprecated methods, including
  pymatgen.alchemy.materials.TransformedMaterial.undo/redo_last_transformation,
  pymatgen.core.site.Site.distance_and_image_old, Poscar.struct,
  StructureFitter and tests.
* Miscellaneous bug fixes and speedups.

v2.4.3
------
* Bug fix for StructureMatcher.
* Miscellaneous speedups.

v2.4.0
------
* New StructureMatcher that effectively replaces StructureFitter. Orders of
  magnitude faster and more robust. StructureFitter is now deprecated.
* Vastly improved PrimitiveCellTransformation.
* A lot of core methods have been rewritten to take advantage of vectorization
  in numpy, resulting in orders of magnitude improvement in speed.
* Miscellaneous bug fixes and speedups.

v2.3.2
------
* More utilities for working with Periodic Boundary Conditions.
* Improved MPRester that supports more data and a new method of specifying
  the API key for heavy users via a MAPI_KEY environment variable. Please
  refer to the :doc:`usage pages </usage>` for more information.
* Vastly improved POTCAR setup script in scripts directly that is now
  installed as part of a default pymatgen installation.
* Miscellaneous bug fixes and speedups.

v2.3.1
------
* Significant improvements to the high-level interface to the Materials API.
  New interface provides more options to make it easier to get structures and
  entries, better warnings and error handling. It uses the *requests*
  library for a cleaner API.
* Bug fix for VolumetricData parsing and methods such as CHGCAR and LOCPOT.
  Previously, the parsing was done incorrectly because VASP actually provides
  data by running through the x-axis first, followed by y, then z.
* Bug fix for reverse_readline so that it works for gzipped and bzipped
  strucutures (courtesy of Anubhav Jain).
* Fix "lossy" composition to_dict method.  Now composition.to_dict properly
  returns a correct species string as a key for compositions using species,
  instead of just the element symbols.
* Miscellaneous bug fixes.

v2.3.0
------
* Remove usage of scipy and external qhull callers. Now uses pyhull package.
  Please note that this change implies that the pyhull package is now a
  required dependency. If you install pymatgen through the usual
  easy_install or pip install methods, this should be taken care of
  automatically for you. Otherwise, please look for the pyhull package on
  PyPI to download and install it.
* Miscellaneous bug fixes.

v2.2.6
------
* Brand new *beta* bond valence analyzer based on a Maximum A Posteriori
  algo using data-mined ICSD data.
* Speed up and improvements to core classes.
* Improved structure fitter (credits to Geoffroy Hautier).
* Brand new entry_tools module (pymatgen.entries.entry_tools).
* Vastly improved Outcar parser based on reverse parsing that speeds up
  reading of OUTCAR files by orders of magnitude.
* Miscellaneous bug fixes.

v2.2.4
------
* Fixed bug in hexagonal cell KPOINTS file generation.
* New RelaxationAnalyzer to compare structures.
* New *beta* bond valence analyzer.
* Miscellaneous bug fixes.

v2.2.3
------
* New filter framework for filtering structures in pymatgen.alchemy.
* Updated feff io classes to support FEFF 9.6 and other code improvements.
* Miscellaneous bug fixes.

v2.2.2
------
* Bug fix release for REST interface.
* Improvements to unittests.

v2.2.1
------
* Improvements to feffio.
* Master matgenie.py script which replaces many analysis scripts.
* More memory efficient parsing of VolumetricData.
* Beta version of structure prediction classes.
* Changes to MPRester to work with v1 release of the Materials API.
* Miscellaneous bug fixes and speed improvements.

v2.2.0
------
* Beta modules (pymatgen.io.feffio) for io for FEFF, courtesy of Alan Dozier.
* New smartio module that intelligently reads structure input files based on
  file extension.
* Spglib_adaptor module has been renamed to finder for brevity.
* Upgraded spglib to version 1.2.2. Improved handling of spglib install on
  Mac OS X and Solaris.
* Major cleanup of code for PEP8 compliance.
* Cssr module now supports reading of input files.
* Miscellaneous bug fixes and speed improvements.

v2.1.2
------
* Brand new CompoundPD class that allows the plotting of phase diagrams that
  do not have elements as their terminal points.
* Spglib is now completely integrated as part of the setup.py installation.
* Major (but completely backwards compatible) refactoring of sites and vaspio.
* Added a EnumerateStructureTransformation with optional dependency on the enum
  library by Gus Hart. This provides a robust way to enumerate derivative
  structures,
* Implemented LLL lattice reduction algorithm. Also added option to sanitize
  a Structure on copy.
* Bug fix for missing Compatibility file in release distribution.
* Vastly improved StructureFitter which performs cell reduction where necessary
  to speed up fitting.
* Miscellaneous bug fixes and speed improvements.

v2.0.0
------
* Brand new module (pymatgen.matproj.rest) for interfacing with the
  MaterialsProject REST interface.
* Useful aliases for commonly used Objects, similar in style to numpy.
  Supported objects include Element, Composition, Structure, Molecule, Spin
  and Orbital. For example, the following will now work::

      import pymatgen as mg
      # Elemental Si
      fe = mg.Element("Si")
      # Composition of Fe2O3
      comp = mg.Composition("Fe2O3")
      # CsCl structure
      structure = mg.Structure(mg.Lattice.cubic(4.2), ["Cs", "Cl"],
                               [[0, 0, 0], [0.5, 0.5, 0.5]])

* New PDAnalyzer method to generate chemical potential maps.
* Enhanced POSCAR class to support parsing of velocities and more formatting
  options.
* Reorganization of Bandstructure module. Beta support for projected
  bandstructure and eigenvalues in vaspio and electronic_structure.
* Miscellaneous bug fixes and speed improvements.

v1.9.0
------
* Completely new json encoder and decoder that support serialization of almost
  all pymatgen objects.
* Simplification to Borg API utilizing the new json API.
* Bandstructure classes now support spin-polarized runs.
* Beta classes for battery (insertion and conversion) analysis.

v1.8.3
------
* spglib_adaptor now supports disordered structures.
* Update to support new spglib with angle_tolerance.
* Changes to Borg API to support both file and directory style paths.
* Speed up for COMPLETE_ORDERING algo for PartialRemoveSpecieTransformation.


v1.8.1
------
* Revamped transmuter classes for better readability and long term support.
* Much improved speed for PartialRemoveSpecieTransformations.
* Misc bug fixes.

v1.8.0
------
* Support for additional properties on Specie (Spin) and Site (magmom, charge).
* Molecule class to support molecules without periodicity.
* Beta io class for XYZ and GaussianInput.

v1.7.2
------
* Bug fixes for vaspio_set and compatibility classes.

v1.7.0
------
* Complete reorganization of modules for electronic structure.
* Beta of band structure classes.
* Misc improvements to vaspio classes.
* Bug fixes.

v1.6.0
------
* Beta of pymatgen.borg package implemented for high-throughput data assimilation.
* Added ComputedEntry classes for handling calculated data.
* New method of specifying VASP pseudopotential location using a VASP_PSP_DIR
  environment variable.
* Bug fix for pymatgen.symmetry
* Ewald sum speed up by factor of 2 or more.
