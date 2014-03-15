Change log
==========

v2.8.10
-------
1. Refactoring of qchemio module (by Xiaohui Qu).

v2.8.9
------
1. qchemio module (by Xiaohui Qu).

v2.8.8
------
1. Minor bug fix release for Structure species substitution methods.

v2.8.7
------
1. Massive update to pymatgen.io.abinitio package (by Matteo Giantomassi).
2. Bug fixes for StructureMatcher's group_structure.
3. Misc bug fixes and cleanup.

v2.8.6
------
1. Bug fix for VASP io set introduced by the default to sorting of structure
   sites when generating VASP input.

v2.8.4
------
1. Completely revamped Compatibility/Correction system which improves
   readability (Shyue Ping Ong/Anubhav Jain/Sai Jayaraman). This change is
   backwards compatible for the most part.

v2.8.3
------
1. Big fix release for json dumping for unitized floats.

v2.8.2
------
1. Bug fix release to improve CIF parsing for more non-standard CIF files.
   In particular, non-ascii characters are removed and _cgraph* fields are
   removed prior to parsing for better support in PyCiFRW.

v2.8.1
------
1. Bug fix release. Incorrect units assigned for ionic radii.
2. Improved nwchemio supports COSMO and ESP calculations (Nav Rajput).

v2.8.0
------
1. **Units**. Pymatgen now has a new system of managing units,
   defined in pymatgen.core.units. Typical energy, length, time,
   temperature and charge units are supported. Units subclass float,
   which makes the usage transparent in all functions. The value that they
   being are in terms of conversions between different units and addition and
   subtraction of different units of the same type. Some basic quantities
   like ionic radii and atomic masses are now returned in unitized forms for
   easy conversion. Please see :mod:`pymatgen.core.units` and the
   :doc:`examples </examples>` for a demonstration of house to use units in
   pymatgen.
2. **Minor backwards-incompatible change**. Structures are now sorted by
   default when generating VASP input files using vaspio_set. Old behavior can
   be obtained by setting sort_structure=False in the constructor. This is
   typically the desired behavior and prevents the generation of large
   POTCARs when atomic species are not grouped together.
3. Bug fix for Molecule.substitute. Earlier algorithm was not detecting
   terminal atoms properly.
4. Additional conversion tools for ABINIT (by Matteo Giantomassi).

v2.7.9
------
1. Minor bug fix release to fix pyhull dependencies to be more friendly.
2. Improved structure matcher that allows for more flexible matching. New
   matching between ordered and disordered comparator.

v2.7.7
-------
1. Beta new Gulp Caller and Zeo++ interface classes (Bharat . Zeo++ is an open
   source software for performing high-throughput geometry-based analysis of
   porous materials and their voids. Please see
   http://www.maciejharanczyk.info/Zeopp/about.html.
2. Specify version of distribute to 0.6.34 for better compatibility.

v2.7.6
------
1. Support for VTK 6.x in structure visualization.
2. Updated install instructions for openbabel.
3. Preliminary pourbaix analysis (Sai Jayaratnam).

v2.7.5
------
1. Vastly improved Nwchem IO (by Shyue Ping Ong).
2. Much improved ABINIT support (by Matteo Giantomassi).

v2.7.4
------
1. Added basic Nwchem (http://www.nwchem-sw.org/) IO support. (by: Shyue Ping
   Ong).
2. New MoleculeMatcher class for comparing molecules by RMS. Requires
   openbabel with python bindings. (by: Xiaohui Qu)
3. New functional group substitution capability for molecules (by: Lei Cheng
   and Shyue Ping Ong).

v2.7.2
------
1. Minor bug fix release to fix some rare errors in very high dimensional
   phase diagrams. **Requires new pyhull version (1.3.8).**

v2.7.1
------
1. **Major backwards-incompatible change.** With effect from v2.7.1,
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
2. read_structure and write_structure now supports pymatgen's json serialized
   structures.
3. read_mol and write_mol functions now available (analogues of
   read_structure and write_structure for molecules)

v2.7.0
------
1. Beta support for ABINIT input and output via pymatgen.io.abinitio
   (courtesy of the excellent work of Matteo Giantomassi).
2. Properties are now checked when comparing two Species for equality.
3. MaterialsProjectVaspInputSet is now renamed to MPVaspInputSet for easier
   typing. The old input sets have been deprecated.
4. New VaspInputSets for MPStatic, MPNonSCF, MITMD which supports uniform
   grid, bandstructure and molecular dynamics calculations. The MD input set
   uses MIT parameters for speed.
5. A beta DiffusionAnalysis class in the apps package.
6. A revised KPOINT grid algorithm that generates more reasonable meshes.
7. A guided install script is now provided for Mac and Linux users.

v2.6.6
------
1. Updates to feffio (credit: Alan Dozier)
2. Added detailed installation instructions for various platforms.
3. Support for charge and spin multiplicity in Molecule. Expanded methods
   available in Molecule.
4. Added supercell matching capabilities to StructureMatcher.
5. More robust creation of PhaseDiagrams to take into account potential qhull
   precision errors.

v2.6.5
------
1. Added a command_line caller to do Bader charge analysis using Henkelmann
   et al.'s algorithm.
2. Bug fix for POSCAR parsing when title line is an empty string.
3. Added __rmul__ operator for Composition.
4. Vastly expanded available aliases.

v2.6.4
------
1. Bug fixes for selective dynamics in Poscar.
2. Improved Procar parsing to support both simple and detailed PROCARs.

v2.6.3
------
1. Added new MaterialsProject REST interfaces for submit/query/delete_snl
   (currently open in beta for collaborators only).
2. Added support for new MaterialsProject REST method get_stability.
3. Added aliases for PhaseDiagram, GrandPotentialPhaseDiagram,
   PDAnalyzer and PDPlotter in pymatgen.phasediagrams.
4. Improvements to StructureMatcher: stol (site - tolerance) redefined as
   a fraction of the average length per atom. Structures matched in fractional
   space are now also matched in cartesian space and a rms displacement
   normalized by length per atom can be returned using the rms_dist method.

v2.6.2
------

1. Site and PeriodicSite now uses a Composition mapping type to represent
   the species and occupancy, instead of a standard dict.
2. Bug fix for reading and re-writing out of Potcars.
3. VaspInputSet now supports MSONable framework.
4. Strain cell option in StructureEditor.
5. Miscellaneous bug fixes and speedups.

v2.6.1
------
1. Use requests.Session in MPRester for connection pooling and code simplicity.
2. Support for "with" context manager in MPRester.
3. Updated periodic table data to correct errors in Ru, Tc and other elements.
4. New methods in Lattice to obtain Wigner-Seitz cell and Brillouin Zone.
5. Miscellaneous bug fixes and speedups.

v2.5.5
------

1. Bug fix release for cifio for rhombohedral structures.
2. Miscellaneous bug fixes and speedups.

v2.5.4
------
1. Vastly improved Gaussian input file parsing that supports more varieties
   of input specifications.
2. StructureNL now supports molecules as well as structures.
3. Updated atomic and vdw radius for Elements.
4. Miscellaneous bug fixes and speedups.

v2.5.3
------
1. Bug fix for StructureNotationalLanguage.
2. Support for LDA US potential. matgenie.py script option to generate POTCARs.
3. Beta version of StructureNotationLanguage, a markup format for Structure
   data with metadata such as authors and references. (Anubhav Jain)
4. Vasprun parsing now parses dielectric constant where available. (Geoffroy
   Hautier)
5. New custom ipython shell script for pymatgen.
6. Miscellaneous bug fixes and speedups.

v2.5.1
------
1. Bug fixes for primitive cell finder.
2. Remove deprecated use_external_qhull option in PhaseDiagram classes.
3. Miscellaneous bug fixes and speedups.

v2.5.0
------
1. Added optimization package with linear assignment class.
2. Improved robustness of StructureMatcher using linear assignment.
3. Improved primitive cell search (faster and more robust).
4. Cleanup of deprecated methods, including
   pymatgen.alchemy.materials.TransformedMaterial.undo/redo_last_transformation,
   pymatgen.core.site.Site.distance_and_image_old, Poscar.struct,
   StructureFitter and tests.
5. Miscellaneous bug fixes and speedups.

v2.4.3
------
1. Bug fix for StructureMatcher.
2. Miscellaneous speedups.

v2.4.0
------
1. New StructureMatcher that effectively replaces StructureFitter. Orders of
   magnitude faster and more robust. StructureFitter is now deprecated.
2. Vastly improved PrimitiveCellTransformation.
3. A lot of core methods have been rewritten to take advantage of vectorization
   in numpy, resulting in orders of magnitude improvement in speed.
4. Miscellaneous bug fixes and speedups.

v2.3.2
------
1. More utilities for working with Periodic Boundary Conditions.
2. Improved MPRester that supports more data and a new method of specifying
   the API key for heavy users via a MAPI_KEY environment variable. Please
   refer to the :doc:`usage pages </usage>` for more information.
3. Vastly improved POTCAR setup script in scripts directly that is now
   installed as part of a default pymatgen installation.
4. Miscellaneous bug fixes and speedups.

v2.3.1
------
1. Significant improvements to the high-level interface to the Materials API.
   New interface provides more options to make it easier to get structures and
   entries, better warnings and error handling. It uses the *requests*
   library for a cleaner API.
2. Bug fix for VolumetricData parsing and methods such as CHGCAR and LOCPOT.
   Previously, the parsing was done incorrectly because VASP actually provides
   data by running through the x-axis first, followed by y, then z.
3. Bug fix for reverse_readline so that it works for gzipped and bzipped
   strucutures (courtesy of Anubhav Jain).
4. Fix "lossy" composition to_dict method.  Now composition.to_dict properly
   returns a correct species string as a key for compositions using species,
   instead of just the element symbols.
5. Miscellaneous bug fixes.

v2.3.0
------
1. Remove usage of scipy and external qhull callers. Now uses pyhull package.
   Please note that this change implies that the pyhull package is now a
   required dependency. If you install pymatgen through the usual
   easy_install or pip install methods, this should be taken care of
   automatically for you. Otherwise, please look for the pyhull package on
   PyPI to download and install it.
2. Miscellaneous bug fixes.

v2.2.6
------
1. Brand new *beta* bond valence analyzer based on a Maximum A Posteriori
   algo using data-mined ICSD data.
2. Speed up and improvements to core classes.
3. Improved structure fitter (credits to Geoffroy Hautier).
4. Brand new entry_tools module (pymatgen.entries.entry_tools).
5. Vastly improved Outcar parser based on reverse parsing that speeds up
   reading of OUTCAR files by orders of magnitude.
6. Miscellaneous bug fixes.

v2.2.4
------

1. Fixed bug in hexagonal cell KPOINTS file generation.
2. New RelaxationAnalyzer to compare structures.
3. New *beta* bond valence analyzer.
4. Miscellaneous bug fixes.

v2.2.3
------

1. New filter framework for filtering structures in pymatgen.alchemy.
2. Updated feff io classes to support FEFF 9.6 and other code improvements.
3. Miscellaneous bug fixes.

v2.2.2
------

1. Bug fix release for REST interface.
2. Improvements to unittests.

v2.2.1
------

1. Improvements to feffio.
2. Master matgenie.py script which replaces many analysis scripts.
3. More memory efficient parsing of VolumetricData.
4. Beta version of structure prediction classes.
5. Changes to MPRester to work with v1 release of the Materials API.
6. Miscellaneous bug fixes and speed improvements.

v2.2.0
------

1. Beta modules (pymatgen.io.feffio) for io for FEFF, courtesy of Alan Dozier.
2. New smartio module that intelligently reads structure input files based on
   file extension.
3. Spglib_adaptor module has been renamed to finder for brevity.
4. Upgraded spglib to version 1.2.2. Improved handling of spglib install on
   Mac OS X and Solaris.
5. Major cleanup of code for PEP8 compliance.
6. Cssr module now supports reading of input files.
7. Miscellaneous bug fixes and speed improvements.

v2.1.2
------

1. Brand new CompoundPD class that allows the plotting of phase diagrams that
   do not have elements as their terminal points.
2. Spglib is now completely integrated as part of the setup.py installation.
3. Major (but completely backwards compatible) refactoring of sites and vaspio.
4. Added a EnumerateStructureTransformation with optional dependency on the enum
   library by Gus Hart. This provides a robust way to enumerate derivative
   structures,
5. Implemented LLL lattice reduction algorithm. Also added option to sanitize
   a Structure on copy.
6. Bug fix for missing Compatibility file in release distribution.
7. Vastly improved StructureFitter which performs cell reduction where necessary
   to speed up fitting.
8. Miscellaneous bug fixes and speed improvements.

v2.0.0
------

1. Brand new module (pymatgen.matproj.rest) for interfacing with the
   MaterialsProject REST interface.
2. Useful aliases for commonly used Objects, similar in style to numpy.
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

3. New PDAnalyzer method to generate chemical potential maps.
4. Enhanced POSCAR class to support parsing of velocities and more formatting
   options.
5. Reorganization of Bandstructure module. Beta support for projected
   bandstructure and eigenvalues in vaspio and electronic_structure.
6. Miscellaneous bug fixes and speed improvements.

v1.9.0
------

1. Completely new json encoder and decoder that support serialization of almost
   all pymatgen objects.
2. Simplification to Borg API utilizing the new json API.
3. Bandstructure classes now support spin-polarized runs.
4. Beta classes for battery (insertion and conversion) analysis.

v1.8.3
------

1. spglib_adaptor now supports disordered structures.
2. Update to support new spglib with angle_tolerance.
3. Changes to Borg API to support both file and directory style paths.
4. Speed up for COMPLETE_ORDERING algo for PartialRemoveSpecieTransformation.


v1.8.1
------

1. Revamped transmuter classes for better readability and long term support.
2. Much improved speed for PartialRemoveSpecieTransformations.
3. Misc bug fixes.

v1.8.0
------

1. Support for additional properties on Specie (Spin) and Site (magmom, charge).
2. Molecule class to support molecules without periodicity.
3. Beta io class for XYZ and GaussianInput.

v1.7.2
------

1. Bug fixes for vaspio_set and compatibility classes.

v1.7.0
------

1. Complete reorganization of modules for electronic structure.
2. Beta of band structure classes.
3. Misc improvements to vaspio classes.
4. Bug fixes.

v1.6.0
------

1. Beta of pymatgen.borg package implemented for high-throughput data assimilation.
2. Added ComputedEntry classes for handling calculated data.
3. New method of specifying VASP pseudopotential location using a VASP_PSP_DIR
   environment variable.
4. Bug fix for pymatgen.symmetry
5. Ewald sum speed up by factor of 2 or more.
