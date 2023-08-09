---
layout: default
title: Change Log
nav_order: 4
---

# Changelog

## v2023.7.20

- Unreadable string concat ops to f-string by @janosh in <https://github.com/materialsproject/pymatgen/pull/3162>
- Revert `mp-api<0.34.0` pin by @janosh in <https://github.com/materialsproject/pymatgen/pull/3165>
- Fix CI error `"pdentries_test.csv"` not found by @janosh in <https://github.com/materialsproject/pymatgen/pull/3168>
- Fix issues with labels by @stefsmeets in <https://github.com/materialsproject/pymatgen/pull/3169>

## v2023.7.17

- Cython 3.0 support.
- PR #3157 from @mattmcdermott magnetic-analyzer-fix. Fixes bug briefly mentioned in #3070, where recent
  spin property changes resulted in the `MagneticStructureEnumerator` failing. This is apparently due to
  creating structures where only some `Species.spin` properties are defined, causing
  CollinearMagneticStructureEnumerator` to fail.
- PR #3070 from @mattmcdermott magnetic-enumerator-fix. To summarize: changes to default magnetic moments
  introduced in #2727 now mean that structures with only partially defined magnetic moments (e.g., on
  half the sites) cannot be successfully analyzed by `SpaceGroupAnalyzer`. This was encountered when
  performing magnetic ordering enumeration, as the previous default behavior for `
  MagOrderingTransformation` does not implicitly yield spins of 0 on the nonmagnetic sites. This has now
  been fixed.

## v2023.7.14

- Emergency bug fix release to remove use of sys.path in pymatgen.io.ase package.
- Fix "Incompatible POTCAR" error on ComputedEntries with oxidation states.
- New global config variable `PMG_POTCAR_CHECKS` provides means to disable all POTCAR checking.

## v2023.7.11

- Use joblib to speed up expensive enumeration energy computations.
- Minor cleanups.

## v2023.6.28

- Use lru_cache to speed up get_el_sp by 400x (@v1kko).
- Related to lru_cache of get_el_sp, Species.properties is now deprecated in favor of setting Species(spin=5). The rationale is
  that spin is the only supported property for Species anyway. Species and DummySpecies is now mostly immutable, i.e., setting specie.spin = 5 have no effect. This is as intended since the first version of pymatgen.
- PR #3111 from @xjf729 fix-MoleculeGraph-draw_graph
- PR #3030 from @lbluque Remove superfluous structure argument docstring from `SQSTransformation` init
- PR #3031 from @kavanase Quick fix to allow direct initialisation of the `DictSet` class.
- PR #3015 from @lbluque Optimized cython code in `find_points_in_spheres`, getting ~5x faster runtime.

## v2023.6.23

- PR #3062 from @arosen93 asefix
    Closes #3061. @JaGeo
- PR #3030 from @lbluque master
    Remove superfluous structure argument docstring from `SQSTransformation` init
- PR #3031 from @kavanase master
    This is a quick fix to allow direct initialisation of the `DictSet` class, which was possible before but broke in <https://github.com/materialsproject/pymatgen/pull/2972> due to the `Yb_2` check querying `self.CONFIG`, which is only defined if `DictSet` was being initialised from a subclass and not directly.
- PR #3015 from @lbluque neighbors
    Optimized cython code in `find_points_in_spheres`, getting ~5x faster runtime.

## v2023.5.31

- Breaking: Default user_potcar_settings to {"W": "W_sv"} in all input sets if user_potcar_functional == "PBE_54" `#3022 <https://github.com/materialsproject/pymatgen/pull/3022>`_
- Unignore ruff PD011 `#3020 <https://github.com/materialsproject/pymatgen/pull/3020>`_
- Tweak variable names `#3019 <https://github.com/materialsproject/pymatgen/pull/3019>`_
- MaterialsProjectCompatibility issue silencable deprecation warning `#3017 <https://github.com/materialsproject/pymatgen/pull/3017>`_
- Optimize cython find_points_in _spheres `#3015 <https://github.com/materialsproject/pymatgen/pull/3015>`_
- Cp2k 2.0 `#2672 <https://github.com/materialsproject/pymatgen/pull/2672>`_
- Added methods to compute and compare DOS fingerprints `#2772 <https://github.com/materialsproject/pymatgen/pull/2772>`_
- Breaking: Overhaul class PymatgenTest `#3014 <https://github.com/materialsproject/pymatgen/pull/3014>`_
- Fix ValueError when structure.selective_dynamics has type np.array `#3012 <https://github.com/materialsproject/pymatgen/pull/3012>`_
- Clean up `#3010 <https://github.com/materialsproject/pymatgen/pull/3010>`_
- Update .pytest-split-durations `#3005 <https://github.com/materialsproject/pymatgen/pull/3005>`_
- Lookup MPRester API key in settings if None provided as arg `#3004 <https://github.com/materialsproject/pymatgen/pull/3004>`_
- Support writing structures to compressed JSON (.json.gz .json.bz2 .json.xz .json.lzma) `#3003 <https://github.com/materialsproject/pymatgen/pull/3003>`_
- Add LightStructureEnvironments.from_structure_environments() fallback value if ce_and_neighbors is None `#3002 <https://github.com/materialsproject/pymatgen/pull/3002>`_
- Species parse oxi state from symbol str `#2998 <https://github.com/materialsproject/pymatgen/pull/2998>`_
- Re-export SiteCollection + DummySpecies from pymatgen.core `#2995 <https://github.com/materialsproject/pymatgen/pull/2995>`_
- Orbital-resolved icohplist `#2993 <https://github.com/materialsproject/pymatgen/pull/2993>`_
- Hide all type-hint-only imports behind if TYPE_CHECKING `#2992 <https://github.com/materialsproject/pymatgen/pull/2992>`_
- Add type hints for pymatgen.io.ase module `#2991 <https://github.com/materialsproject/pymatgen/pull/2991>`_
- Enable ruff doc rules in CI `#2990 <https://github.com/materialsproject/pymatgen/pull/2990>`_
- Suspected Typo Fix in pymatgen.io.vasp.optics `#2989 <https://github.com/materialsproject/pymatgen/pull/2989>`_
- Doc strings `#2987 <https://github.com/materialsproject/pymatgen/pull/2987>`_
- Fix average error `#2986 <https://github.com/materialsproject/pymatgen/pull/2986>`_
- Drop deprecated SubstrateAnalyzer + ZSLGenerator reexports `#2981 <https://github.com/materialsproject/pymatgen/pull/2981>`_
- Breaking: Default user_potcar_settings to {"W": "W_sv"} in all input sets if user_potcar_functional == "PBE_54" (#3022) `#3022 <https://github.com/materialsproject/pymatgen/pull/3022>`_
- fix unwanted x margins in get_elt_projected_plots_color (closes #562) `#562 <https://github.com/materialsproject/pymatgen/issues/562>`_
- Add LightStructureEnvironments.from_structure_environments() fallback value if ce_and_neighbors is None (#3002) `#2756 <https://github.com/materialsproject/pymatgen/issues/2756>`_
- add doc str explaining need for class ElementBase (closes #2999) `#2999 <https://github.com/materialsproject/pymatgen/issues/2999>`_
- Update docs. `3e3c31c <https://github.com/materialsproject/pymatgen/commit/3e3c31c8d342c84f2c6bbb961c321e458b9accb9>`_
- ruff set isort.split-on-trailing-comma = false `c0ec534 <https://github.com/materialsproject/pymatgen/commit/c0ec53452c3dc87c6cca5edc1c6b2b6218f15569>`_

## v2023.5.10

- Fix mem leak in pbc_shortest_vector cython code. (@stichri)
- Set all cython code to language level 3.

## v2023.5.8

❗ The Yb_2 deprecation release ❗

This release changes the Ytterbium (Yb) pseudo-potential (PSP) from Yb_2 to Yb_3 for all PBE_54 VASP input sets.

Background: The `A-lab <https://newscenter.lbl.gov/2023/04/17/meet-the-autonomous-lab-of-the-future>`_ revealed that as a result of using Yb_2 the energy on Yb compounds is off by a lot, resulting in supposedly stable things being unsynthesizable. While an unfortunate mistake, it's also great to see how experiment can help surface simulation errors.

On pre-PBE_54 input sets, we now issue a warning that Yb_2 will give bad results for most systems since Yb is most often in oxidation state Yb3+.

Reason: The better fix Yb_3 only became available in the PBE_54 PSP set. Requiring it on pre-PBE_54 input sets would mean you can't run Yb compounds.

For more details see `#2968 <https://github.com/materialsproject/pymatgen/pull/2968>`_and `#2969 <https://github.com/materialsproject/pymatgen/pull/2969>`_.

What's Changed

- Fix TypeError: a bytes-like object is required, not 'list' when passing triplet of bools to find_points_in_spheres() pbc kwarg  by @janosh in `#2907 <https://github.com/materialsproject/pymatgen/pull/2907>`_
- Fix ValueError: not enough values to unpack in PDPlotter if no unstable entries in PD by @janosh in `#2908 <https://github.com/materialsproject/pymatgen/pull/2908>`_
- Fix VolumetricData.to_cube() not preserving structure dimensions by @janosh in `#2909 <https://github.com/materialsproject/pymatgen/pull/2909>`_
- Update team.rst by @jmmshn in `#2912 <https://github.com/materialsproject/pymatgen/pull/2912>`_
- Faff by @janosh in `#2915 <https://github.com/materialsproject/pymatgen/pull/2915>`_
- Add formal_chempots option to ChemicalPotentialDiagram to plot the formal chemical potentials rather than the DFT energies by @kavanase in `#2916 <https://github.com/materialsproject/pymatgen/pull/2916>`_
- Modified dosplotter by @kaueltzen in `#2844 <https://github.com/materialsproject/pymatgen/pull/2844>`_
- auto version by @jmmshn in `#2925 <https://github.com/materialsproject/pymatgen/pull/2925>`_
- bug fix for potcar parsing by @jmmshn in `#2910 <https://github.com/materialsproject/pymatgen/pull/2910>`_
- Fix breaking changes from pandas v2 by @janosh in `#2935 <https://github.com/materialsproject/pymatgen/pull/2935>`_
- add kwarg to MoleculeGraph method and fix PackmolSet bug by @orionarcher in `#2927 <https://github.com/materialsproject/pymatgen/pull/2927>`_
- fix on reading multiple route in Gaussian input file by @Ameyanagi in `#2939 <https://github.com/materialsproject/pymatgen/pull/2939>`_
- Fix CI errors by @janosh in `#2940 <https://github.com/materialsproject/pymatgen/pull/2940>`_
- Add ResParser support for reading files with spin values by @ScottNotFound in `#2941 <https://github.com/materialsproject/pymatgen/pull/2941>`_
- Ignore bad unicode characters in Structure.from_file() by @janosh in `#2948 <https://github.com/materialsproject/pymatgen/pull/2948>`_
- Minor modification for symmetrically distinct Miller index generation by @fyalcin in `#2949 <https://github.com/materialsproject/pymatgen/pull/2949>`_
- Fixed Wulff shape for new versions of matplotlib by @CifLord in `#2950 <https://github.com/materialsproject/pymatgen/pull/2950>`_
- Test figure returned by WulffShape.get_plot() contains single Axes3D by @janosh in `#2953 <https://github.com/materialsproject/pymatgen/pull/2953>`_
- Fix Cp2kOutput.spin_polarized() likely not doing what author intended by @janosh in `#2954 <https://github.com/materialsproject/pymatgen/pull/2954>`_
- For MPcules: Molecule Trajectory and graph hashes by @espottesmith in `#2945 <https://github.com/materialsproject/pymatgen/pull/2945>`_
- self.assertArrayEqual->assert by @janosh in `#2955 <https://github.com/materialsproject/pymatgen/pull/2955>`_
- fix GaussianOutput bug with multiple route lines by @xjf729 in `#2937 <https://github.com/materialsproject/pymatgen/pull/2937>`_
- Fix ValueError when passing selective_dynamics to Poscar by @chiang-yuan in `#2951 <https://github.com/materialsproject/pymatgen/pull/2951>`_
- Bump beautifulsoup4 from 4.11.1 to 4.12.2 by @dependabot in `#2962 <https://github.com/materialsproject/pymatgen/pull/2962>`_
- Bump pypa/cibuildwheel from 2.11.4 to 2.12.3 by @dependabot in `#2959 <https://github.com/materialsproject/pymatgen/pull/2959>`_
- Bump uncertainties from 3.1.6 to 3.1.7 by @dependabot in `#2960 <https://github.com/materialsproject/pymatgen/pull/2960>`_
- Bump numpy from 1.23.2 to 1.24.3 by @dependabot in `#2963 <https://github.com/materialsproject/pymatgen/pull/2963>`_
- Bump tabulate from 0.8.10 to 0.9.0 by @dependabot in `#2961 <https://github.com/materialsproject/pymatgen/pull/2961>`_
- Bump pandas from 1.4.4 to 2.0.1 by @dependabot in `#2964 <https://github.com/materialsproject/pymatgen/pull/2964>`_
- Link /addons from new subsection on /contributing page by @janosh in `#2967 <https://github.com/materialsproject/pymatgen/pull/2967>`_
- Breaking: change Yb pseudo-potential on all VASP input sets from Yb_2 to Yb_3 by @janosh in `#2969 <https://github.com/materialsproject/pymatgen/pull/2969>`_
- fix recursion error by adding copy and deepcopy dunder methods by @orionarcher in `#2973 <https://github.com/materialsproject/pymatgen/pull/2973>`_
- Revert to Yb_2 on pre-PBE_54 input sets by @janosh in `#2972

## v2023.3.23

- Misc bug fixes.
- Enable Structure relaxations with TrajectoryObserver (@janosh)
- Breaking: Rename `gen_sl_transform_matricies->gen_sl_transform_matrices` (#2894)

## v2023.3.10

- PR #2882 substrate-optimizations for speed up (@mkhorton)
- Fix very fragile POTCAR parsing.

## v2023.2.28

- 69a548e210 revert adding ubuntu-latest to release job matrix
- e63cab3620 add 3.11 to release job python-versions
- c03dacb94d use `cibuildwheel` to build linux wheels (#2800)
- fe2597d92e Merge setup.cfg into pyproject.toml (#2858)
- 40cbf1d7c4 del class AtomicFile, _maketemp(), ask_yesno() from pymatgen/util/io_utils.py (#2860)
- 0b16987f2c fix reduced formula in Ion (#2864)

## v2023.2.22

- PR #2848 from @ml-evs ml-evs/update_optimade_aliases
    Currently `OptimadeRester` defaults to an outdated list of OPTIMADE database URLs (several of which fail) and the design of the class is such that refreshing these aliases can only be done post-init which means they will not be used if the user provides their own filtered list of aliases, without doing some extra work.
    This PR refreshes the vendored list of aliases (which should be much more stable now since their initial addition 2 years ago), and also adds the option to refresh the aliases on initialization of the class.
    This currently affects the pymatgen OPTIMADE tutorials at <https://github.com/Materials-Consortia/optimade-tutorial-exercises>.

## v2023.1.30

- PR #2806 from @samblau qchem
  - Major changes to Q-Chem IO (inputs.py and outputs.py) to accommodate differences and new features in version 6+
  - Additional parsing capabilities for HOMO/LUMO, dipoles, NBO info (hyperbonds and 3C bonds) in outputs.py
  - Utility for processing a parsed binary Hessian scratch file
  - Overdue updates to default values in sets.py and new defaults associated with differences and new features in Q-Chem 6+
- PR #2814 from @jmmshn patch_dos

## Added Convenience to obtain the normalized CompleteDos object

    Added tests to make sure calling it multiple time still only gives one result.

## v2023.1.20

- Passthrough kwargs support for Structure.from_file and Structure.from_str
- Allow the `frac_tolerance` to be specified for rounding coordinates in CifParser.
- PR #2803 from @amkrajewski add_weightbasedfunctions
    When working with metallic alloys, weight-fraction-based notations such as Ti64 / Ti-6V-4Al or NiTiNOL60 / Ni-40Ti are commonly employed in both industrial specifications and scientific literature. Regardless of the numerous downsides of this situation, including errors in scientific experiments or NLP-parsing when they are mistaken for atomic fractions or chemical formulas, being able to create a Composition object from them (under correct interpretation) would be a useful pymatgen feature.
  - Composition class method to initialize it from a dictionary of weight fractions
  - Composition property giving a dictionary of weight fractions
  - concise tests for the two above were added
QChem: translate DMSO name in smd_solvent

## v2023.1.9

- PR #2792 from @JaGeo bug_fix
- PR #2773 from @ab5424 cbar
- PR #2776 from @MichaelWolloch master
- PR #2762 from @MichaelWolloch master
- PR #2774 from @dgaines2 fix-poscar
- PR #2667 from @nwinner volumetric-data-patch
- PR #2764 from @naik-aakash lobster_lsodos
- PR #2215 from @rkingsbury cmirs
- PR #2742 from @materialsproject pip-dependabot
- PR #2741 from @materialsproject resurrect-req-txt
- PR #2735 from @njzjz patch-1

## v2022.11.7

- PR #2724 from @janosh: raise ValueError in SpacegroupAnalyzer.get_symmetrized_structure() if spglib returns no symmetries
- PR #2720 by @utf: Fix tensor mapping
- PR #2562 from @sudarshanv01: In case the Fock-matrix and eigenvalues are requested by the user (though the flags `scf_final_print` or `scf_print`), outputs.py now allows parsing both these quantities.

## v2022.11.1

- Order of kwargs `fmt` and `filename` in `Structure.to()` swapped for ease of use (note: this can break codes that do not use these options as kwargs).
- @yuzie007 Parse "Atomic configuration" in POTCAR (52 and 54). Useful for estimating a reasonable NBANDS value.
- EnumerateStructureTransformation now supports `m3gnet_relax` or `m3gnet_static` options.

## v2022.10.22

- Allow env settings to override .pmgrc.yaml (@janosh)
- Add EntryLike type (@janosh)
- Update spglib to 2.0+.
- @cnncnnzh Method to plot the atom-resolved phonon band structures.
- @jmmshn More Flexible reproduction of VASP's optical code
- @Ameyanagi Fix the sorting of the FEFF IO module to create ATOMS input.
- @JaGeo Extend the ThermalDisplacementMatrices class to read cif files in P1 format.
- @rkingsbury Changes to FEFF I/O to support the use of non-periodic input structures.
- @jmmshn Merge Waverder and Wavederf
- @jmmshn  Set the structure_charge while parsing Potcar

## v2022.9.21

- @chunweizhu fix the bugs when running `TEMCalculator`
- @munrojm Support for new MPRester.

## v2022.9.8

- @janosh Add AirssProvider.as_dict
- @gpetretto Outcar parsing optimization.
- @ScottNotFound Adds res file io to handle results from airss searches
- @janosh Fixes the `AttributeError` currently raised when passing disordered structures to methods like `get_cn()` and `get_bonded_structure()` of `CrystalNN` and other `NearNeighbors` subclasses.
- @naik-aakash Added new option `standard_with_comp_range` for generating lobsterin files using vasp

## v2022.8.23

- Structure Graphs from Lobster Data (@JaGeo)
- Added 'get_orbit_and_generators'-method to SpaceGroup class (@nheinsdorf)
- Class to handle Thermal displacements matrices (@JaGeo)
- Change default number of significant digits to write VASP POSCAR (@henriquemiranda)
- Misc bug fixes.

## v2022.7.25

- Implemented sufficient methods for new MPRester to cover about 70-80% of common use cases.

## v2022.7.24.1

- Implementation changed to allow seamless use of MPRester regardless of whether new or old API key is used.

## v2022.7.24

- Initial implementation of MPRester2 with new API support. Basic functionality for now.

## v2022.7.19

This will be the final release with the pymatgen.analysis.defects
module included in the standard pymatgen package. This release will
include the older defects code by default, but can also be replaced with
the newer defects code through installation of pymatgen-analysis-defects.

Subsequent versions of pymatgen will require
the additional installation of `pymatgen-analysis-defects <https://github.com/materialsproject/pymatgen-analysis-defects>`_ for all defect-related
functionality via pip install pymatgen-analysis-defects.

Relevant imports will still be from the pymatgen.analysis.defects namespace but the code will now be maintained and developed in this separate repository.

There will be significant changes to the defects code to support new functionality.Existing PyCDT users should use this version of pymatgen or older. Any questions
about this change should be directed to Jimmy-Xuan Shen, @jmmshn.

For more information about other pymatgen "add-on" packages, please see
`this page in our documentation <https://pymatgen.org/addons.html>`_.

- Preparation for the removal of the defects module, PR #2582 by @jmmshn

## v2022.7.8

Welcome to new contributors @naveensrinivasan, @xivh, @dgaines2, @yang-ruoxi, @cajfisher and @mjwen!

- New: Partial periodic boundary conditions, PR #2429 by @gpetretto
- New: Element.from_name(), PR #2567 by @rkingsbury
- New: Materials Project input set for absorption calculations, PR #2320 by @yang-ruoxi
- Enhancement: compressed LAMMPS and XYZ files in pymatgen.io.lammps, PR #2538 by @ab5424
- Enhancement: remove vertical lines from VoltageProfilePlotter.get_plotly_figure(), PR #2552 by @acrutt
- Enhancement: chemical potential plot background color changed, PR #2559 @jmmshn
- Enhancement: ability to change voronoi_distance_cutoff in ChemEnv, PR #2568 by @JaGeo
- Enhancement: Ion.oxi_state_guesses will use correct charge by default, PR #2566 by @rkingsbury
- Enhancement: Remove not converged warning for VASP AIMD runs, PR #2571 by @mjwen
- Fix: generation of continuous line-mode band structures, PR #2533 by @munrojm
- Fix: duplicate site properties for magnetic moments hwen using `AseAtomsAdaptor`, PR #2545 by @arosen93
- Fix: bug in Grüneisen parameter calculation, PR #2543 by @ab5424
- Fix: allow a comment on final line of KPOINTS file, PR #2549 by @xivh
- Fix: for `Composition.replace` with complex mappings, PR #2555 by @jacksund
- Fix: Implement equality method and fix __iter__ for InputSet, PR #2575 by @rkingsbury
- Fix: use negative charge convention for electron in "update_charge_from_potcar", PR #2577 by @jmmshn
- Fix: ensure charge is applied to initial and final structures parsed from vasprun.xml, PR #2579 by @jmmshn
- Chore: Set permissions for GitHub actions, PR #2547 by @naveensrinivasan
- Chore: Included GitHub actions in the Dependabot config, PR #2548 by @naveensrinivasan
- Documentation: fix typos in pymatgen.symmetry.analyzer docstrings, PR #2561 by @dgaines2
- Documentation: clarification about usage of InputFile, PR #2570 by @orionarcher
- Documentation: Improve messages and warnings, PR #2572 and PR #2573 by @cajfisher
- Documentation: fix typo, PR #2580 by @janosh

Notice: functionality from pymatgen.analysis.defects will be incorporated into a separate add-on package in the future,
see deprecation notice.

## v2022.5.26

- Q-Chem updates to NBO and new geometry optimizer, PR #2521 by @samblau
- Bug fix for VolumetricData, PR #2525 by @jmmshn
- Bug fix for MPRester, PR #2531 by @janosh

## v2022.5.19

- Added option for additional criteria to be passed to MPRester.get_entries_in_chemsys (@shyuep).

## v2022.5.18.1

- Initial support for parsing ML MD runs from vasprun.xml (@shyuep).

## v2022.5.18

- Bug fix for sulfide_type. Sometimes symmetry analysis fails because of tolerance issues. A fallback to analyze all sites.

## v2022.5.17

- PR #2518 from @JaGeo. Fixed wrong line in ICOHPLIST.lobster being read to assess whether orbitalwise interactions are included in these files.
- PR #2520 from @arosen93. Adds a new property to the `PointGroupAnalyzer`: the rotational symmetry number.
- PR #2522 from @jmmshn. Fixes PD JSON serialization.
- PR #2514 from @qianchenqc. Replaced the IALGO tag with ALGO as recommended in the vasp documentation <https://www.vasp.at/wiki/index.php/IALGO>.
- PR #2404 from @nheinsdorf. Added a method that gets all the neighbors up a maximum distance for a Structure, and groups these 'bonds' according to their symmetry.
- PR #2509 from @jacksund Fix NMR Set.

## v2022.4.26

- Fix dipole units in recent vasp versions (at least 6.3, maybe even before) (@fraricci)
- Removed complex numbers from the definition of WSWQ (@jmmshn)
- MP database version logging is now no longer logged in the .pmgrc.yaml but rather in the .mprester.log.yaml.
  This avoids the MPRester constantly rewriting a config file and causing users' pymatgen to completely fail.

## v2022.4.19

- Fix for discharged and charged entries in conversion battery. (@peikai)`pylint` in `.pre-commit-config.yaml`.
- Allow skipping of structure reduction in StructureMatcher.group_structures (@lan496)
- Return NotImplemented for composition comparison methods. (@janosh)
- BSPlotter bug fixes (@fraricci)
- Misc bug fixes and deprecation fixes.

## v2022.3.29

- Major update to CP2K module, PR #2475 from @nwinner
- Bug fix to remove problematic import, PR #2477 from @mkhorton

## v2022.3.24

- Emergency bugfix release to fix circular import (@janosh)

## v2022.3.22

- Support kwargs for ASE adaptor. (@arosen93)
- Fix for cation error in Lobster analysis. (@JaGeo)
- Major revampt of Abstract interface for Input classes in IO. (@rkingsbury)
- Orbital-projected band center, band filling, band center, skewness, kurtosis, etc. (@arosen93)
- Misc cleanups. (@janosh)

## v2022.3.7

- Add VASP WSWQ file parsing, PR #2439 from @jmmshn
- Improve chemical potential diagram plotting, PR #2447 from @mattmcdermott
- Update to Lobster calculation settings, PR #2434 from @JaGeo master
- Allow non-integer G-vector cut-off values when parsing WAVECAR, PR #2410 from @arosen93
- Fix for Structure.from_file when file is in YAML format from @janosh fix-structure-from-yml
- Update of linter configuration, PR #2440 from @janosh
- Update to ChemEnv citation, PR #2448 from @JaGeo
- Type annotation fix, PR #2432 from @janosh
- Documentation fix for Structure.apply_operation, PR #2433 from @janosh
- Add caching to compatibility classes as speed optimization, PR #2450 from @munrojm

This release was previously intended for v2022.2.25.

Important note: an update to a library that pymatgen depends upon has led to the
~/.pmgrc.yml configuration file being corrupted for many users. If you are affected,
you may need to re-generate this file. This issue should now be fixed and not re-occur.

## v2022.2.10

- Require Cython during setup. (@jonringer)

## v2022.2.7

- Critical bug fix for pmgrc.yaml being overwritten in MPRester in a non-standard way.
- Change in config file for Lobster basis. Removed the 2p orbitals for Be as they led to problems in our computations and probably should be optional during the projection. (@JaGeo)
- Return None for ISPIN=1 for `Vasprun('vasprun.xml').complete_dos.spin_polarization`.

## v2022.2.1

- Chargemol caller for partial atomic charge analysis (@arosen93)
- ASEAtomAdaptor: (1) Updates to magmom support, (2) Oxidation states support, (3) Charges are now passed (@arosen93)
- Cleanup of deprecated methods. (@janosh)
- Bigfix for gzipped DOSCAR (@JaGeo)
- Updates for QChem Support (@samblau)
- QuantumEspresso k-grid fix input fix. (@vorwerkc)
- `Entry.__repr__()` now outputs name where available. (@janosh)
- Fixes to Vasprun.final_energy to report `e_0_energy` (the desired energy quantity) for VASP 6+. (@arosen93)
- `Outcar().final_energy` now prints out `e_0_energy` (also called "energy(sigma->0)" in the OUTCAR) rather than `energy_fr_energy` (also called "free  energy   TOTEN" in the OUTCAR). This is to be consistent with `Vasprun().final_energy` and because it is generally the desired quantity. `Outcar` now has two new attributes: `.final_energy_wo_entrp` and `final_fr_energy`, which correspond to `e_wo_entrp` and `e_fr_energy`, respectively. (@arosen93)
- Improved parsing of coupled cluster calculations in QChem (@espottesmith).

## v2022.1.24

- Misc bug fixes, e.g., handling of yaml files and type check for MAGMOM flag.

## v2022.1.20

- Unicode fixes (@janosh)
- YAML deprecation fixes. (@janosh)
- ASE adaptor support for charge, spin multiiciplity and site properties of molecules. (@arosen93).
- New keyword option (`keep_site_properties`) in various `structure.symmetry.analyzer` functions to keep the site properties on the sites after a transformation. (@arosen93)
- Bug fixes for Lobster module (@JaGeo).
- SCAN / GGA(+U) mixing scheme (@rkingsbury). Mixing scheme code lives in the new file `mixing_scheme.py` and is implemented as a `Compatibility` class.
- Fix for parsing of QuantumExpresso files due to new format (@vorwerkc)

## v2022.1.9

- Formal support for Python 3.10.
- Misc refactoring and bug fixes. No new functionality.

## v2022.1.8

- First proper new release of 2022 formalizes the switch back to date-based versioning introduced as a temporary measure last year.
- Numpy version pinned to 1.22.0. This is necessary to avoid binary incompatibility.
- With the numpy version, py37 support is dropped.
- ASE io improvements (e.g., magnetic moments and selective dynamics transfer). @arosen93
- New automatic k-point generation scheme, `automatic_density_by_lengths`, which allows the user to specify a density of k-points in each dimension (rather than just for the entire volume). @arosen93
- Build improvements to dynamically generate C code by running Cython on pyx files rather than having hard-generated .c files.

## v2022.0.17

Welcome to new contributor @e-kwsm!

- More robust smart fermi method by @utf in <https://github.com/materialsproject/pymatgen/pull/2303>
- Replace-species by @janosh in <https://github.com/materialsproject/pymatgen/pull/2291>
- Add warning if improper ALGO is used for hybrid calculations by @arosen93 in <https://github.com/materialsproject/pymatgen/pull/2298>
- Wrap supercell to unit cell when performing change of setting by @jmmshn in <https://github.com/materialsproject/pymatgen/pull/2300>
- Clearer handling of the MAGMOM flag in pymatgen.io.vasp.sets by @arosen93 in <https://github.com/materialsproject/pymatgen/pull/2301>
- Add warning if LASPH != True for meta-GGA/hybrid/vdW/+U by @arosen93 in <https://github.com/materialsproject/pymatgen/pull/2297>
- Add ability to request additional OPTIMADE fields by @ml-evs in <https://github.com/materialsproject/pymatgen/pull/2315>
- Add missing elements to MPScanRelaxSet PBE .54 potentials by @arosen93 in <https://github.com/materialsproject/pymatgen/pull/2316>

- Fix write Trajectory XDATACAR with variable lattice by @gpetretto in <https://github.com/materialsproject/pymatgen/pull/2310>
- Fix small cutoff neighbor by @chc273 in <https://github.com/materialsproject/pymatgen/pull/2277>
- Add Composition.replace() by @janosh in <https://github.com/materialsproject/pymatgen/pull/2284>
- Ion bugfixes and enhancements by @rkingsbury in <https://github.com/materialsproject/pymatgen/pull/2287>
- Fix oddly split strings and a few typos by @janosh in <https://github.com/materialsproject/pymatgen/pull/2285>
- InsertionElectrode bug fix and documentation update by @acrutt in <https://github.com/materialsproject/pymatgen/pull/2257>
- Remove accidentally tracked files and unset executable flag by @e-kwsm in <https://github.com/materialsproject/pymatgen/pull/2296>

- Update DOI URLs by @e-kwsm in <https://github.com/materialsproject/pymatgen/pull/2295>
- Documentation update: Fix missing Outcar attributes and update elemental_dos_dos string by @arosen93 in <https://github.com/materialsproject/pymatgen/pull/2293>
- Documentation update for CutOffDictNN by @ltalirz in <https://github.com/materialsproject/pymatgen/pull/2278>

## v2022.0.16

- Fix to allow PhaseDiagram to be JSON-serializable with computed data cached (@mkhorton, #2276)
- Temporarily revert #2239 pending investigation into slow-down in some nearest neighbor finding routines. This does not affect the behavior of any of these classes.

## v2022.0.15

Welcome to new contributors @blokhin, @pzarabadip, @ml-evs, @wuxiaohua1011, @janssenhenning and @penicillin0. A reminder to all new contributors to
ensure your information is accurate at <https://pymatgen.org/team.html> so that
you are acknowledged appropriately by filling out the linked form.

- Breaking change in PhaseDiagram serialization which will affect any users of BasePhaseDiagram which has now been removed (@shyuep, 2b9911d)

- Speed up nearest-neighbor routines & structure graph generation (@ltalirz, #2239)
- Add two more pre-defined OPTIMADE aliases (@blokhin, #2242)
- Refactor `interface_reactions` module, adding support for Plotly (@mattmcdermott, #2233)

- Update NOMAD access in MPRester (@wuxiaohua1011, #1958)
- General improvements to Phase Diagram code (@CompyRhys, #2263, #2264, #2268)
- Improve appearance of periodic table heatmap (@penicillin0, #2272)
- Small improvements to battery classes (@jmmshn, #2262)
- Fix for Composition.chemical_system to match expected behavior for compositions with oxidation states (@CompRhys, #2249)
- Fix for bad param in OPTIMADE response fields (@ml-evs, #2244)
- Fix for issue in parsing `bandOverlaps.lobster` file (@pzarabadip, #2237)
- Fix for Moladaptor (@orioncohen, #2269)
- Fix for incorrect Potcar hash warnings (@mkhorton, #2273)

- Type hint and correct documentation of Structure.remove_site_properties (@kmu, #2256)
- Type hint improvements across pymatgen (@janosh, #2241, #2247, #2261)
- Add `pymatgen-io-fleur` addon to addons page (@janssenhenning, #2232)

## v2022.0.14

- Update OPTIMADE interface to allow querying multiple providers, this changes the
  method signature of OptimadeRester and so is considered a backwards incompatible change (@mkhorton, #2238)

## v2022.0.13

- New feature to plot chemical potential diagrams (@mattmcdermott, #2218), see ArXiv:2104.05986 for example
- Numerous updates to LOBSTER support for new version and including handling COBICAR, SitePotentials and MadelungEnergies (@JaGeo, #2228)
- Updates and fixes for LAMMPS CombinedData (@htz1992213, #2191)
- Bug fix for Bader caller (@nwinner, #2230)
- Documentation fix for Composition (@CompRhys, #2231)

## v2022.0.12

- @chc273 Major bugfix for cython handling of fractional coordinates wrapping.
- @mattmcdermott Bug fix for entry_ID phase diagram plotting bug described in this Issue: #2219
- @FCMeng Fix for PWSCF to distinguish same element with different oxidation state, which might have different pseudopotentials.
- @gmatteo fix minor bug when reading Structure from a netcdf4 file with hdf5 groups

## v2022.0.11

- New features to handle Grüneisen parameters (@JaGeo, @ab5424, @gpetretto, #2190)
- New option to return SymmetrizedStructure in CifParser (@mkhorton, 0d9a455)
- Fix for SubstrateAnalyzer (@shyamd, #2198)
- Fix for BandFillingCorrection (@kavanase, #2193)

## v2022.0.10

- Add spin-dependent eigenvalue band properties (@arosen93, #2187)
- Bug fix for settings loading (@ardunn, #2186)

## v2022.0.9

- Significant new functionality for handling interfaces between structures (@shyamd, #2149)
- Add input/output for CREST (@arepstein, #2020)
- Add RadialSiteDistortionTransformation (@nwinner, #2108)
- Add Q-Chem NBO functionality (@samblau, #2174)
- Change hkl annotation format in diffraction plots (@flaviu-gostin, #2143)
- Add space group to print output of `SymmetrizedStructure` (@CompRhys, #2139)
- Better error handling in QCOutput (@rkingsbury, #2147, #2165, #2135)
- Add progress bar for applying compatibility scheme (@CompRhys, #2136)
- Allow combining data with multiple molecule IDs in LAMMPS (@htz1992213, #2157)
- Update EDIFF in DFPT input set to be consistent with atomate (@utf, #2172)

- Change names of high-symmetry paths (@munrojm, #2144)
- Change default for filter_solids argument of PourbaixDiagram (@rkingsbury, #2177)

- Fix to improve precision in `FermiDos`, NOTE: this can result in significant changes in some instances (@nwinner, #2109)
- Fix for handling of Exceptions (@kmu, #2150)
- Fix for PourbaixEntry (@JosephMontoya-TRI, #2148)
- Fix for loading of settings from file when environment variables also set (@ardunn, #2164)
- Fix equation for calculation of k-spacing in SCAN sets, NOTE: this now results in a lower k-point density (@ab5424, #2163)
- Fix for parsing of VASP vasprun.xml when ALGO=CHI (@KazMorita, #2171)

- Documentation update for MP2020 corrections scheme (@rkingsbury, #2141)
- Documentation update for SCAN sets (@janosh, #2140)
- Documentation update for using CifWriter (@755452800, #2156)

## v2022.0.8

- PR #2130 @rkingsbury ensures that energy corrections applied to each anion
  have unique names (e.g., N vs. Cl vs. Br).
- PR #2133 @rkingsbury adds support for custom vdW radii to `QCInput` and
  `QChemDictSet`. These radii are used in the construction of PCM cavities and
  when calculating charges.
- PR #2123 from @gpetretto fixes bug in `get_conventional_standard_structure`
  method of the `SpacegroupAnalyzer` for triclinic crystals.
- PR #2134 from @ab5424 supports zopen in parsing lammps logs
- PR #2132 from @htz1992213 speeds up LammpsData.as_string for
  non-hybrid data with large coeff sections and adds as_lammpsdata method to
  CombinedData
- PR #2129 from @CifLord improves analysis of surface symmetry of slabs.
- PR #2117 from @nwinner contains bug fixes for bader caller.

## v2022.0.7

- Improved Gaussian Cube I/O (@nwinner, #2121)
- Updated van der Waals radii (@rkingsbury, #2122)
- Update `MaterialsProject2020Compatibility` for multi-anion systems (@rkingsbury, #2128)
- Fixes and improvements to Q-Chem parsing (@samblau, #2125)
- Bug fix for isseus with hard-coded path in `MaterialsProject2020Compatibility` (@CompRhys, #2124)
- Bug fix for DOS serialization (@zooks97, #2119)
- Bug fix for XDATCAR lattice parsing (@nkeilbart, #2115)
- Documentation link fix (@adam-kerrigan, #2127)

## v2022.0.6

- Feature to calculate Selling vectors and distances between Lattices (@bwjustus, #1888)
- XPS Spectrum class added (@shyuep, #2110, see `galore <https://github.com/SMTG-UCL/galore>`_)
- Updated `MaterialsProject2020Compatibility` for formation energy correction (@rkingsbury, #2106)
- Bug fix for detecting broken bonds in slab generation (@fyalcin, #2015)
- Bug fix for electrodes (@jmmshn, #2101)
- Documentation improvement for get_conventional_standard_structure (@tom-wood, #2100)

## v2022.0.5

- Bug fix to remove possibility of duplicate edges in `StructureGraph` (@mkhorton, #2095)

## v2022.0.4 / v2021.3.9

- Element now has `ionization_energies`, `ionization_energy` and
  `electron_affinity` properties.
- Extensive documentation has been added on pymatgen compatibility and the
  new namespace architecture! We have also released a
  `template repo <https://github.com/materialsproject/pymatgen-addon-template>`_
  to help new developers write add-ons for pymatgen! Check out our
  :doc:`contributing page</contributing>` for details.

## v2022.0.3

- Another bug fix release! Now SETTINGS have been moved to pymatgen.core.

## v2022.0.2 (Yanked)

- Bug fix release for missing package data files in v2022.0.1

## v2022.0.1 (Yanked)

- `pymatgen`, `pymatgen.ext`, `pymatgen.io` and `pymatgen.analysis` are now
  namespace packages. Note that this does not affect normal usage of pymatgen
  from v2022.0.0. All imports remain the same. However, it does allow developers
  to write "add-ons" to these subpackages. A full documentation with examples
  and templates is in the works to guide developers on how to write these
  packages.

## v2022.0.0 (Yanked)

- This is the new version of pymatgen going forward. Root-level imports have been removed. Please see
  <https://pymatgen.org/compatibility.html> on how to update your code for compatibility with v2022.

## v2021.3.5

- Backwards incompatible changes in v2021.3.4 have been removed. Instead another semantic version v2022.0.0 has been
  released. Future critical bug fixes will be backported to v2021.x.x, but the main line of development will occur on
  v2022.0.0 onwards.

## v2021.3.4 (Yanked)

- __Backwards incompatible__: Pymatgen root imports have been removed from
  v2021.3.4 in preparation for a change to a more modular, extensible
  architecture that will allow more developers to contribute.

  If your existing code uses `from pymatgen import <something>`, you will need to make
  modifications. The easiest way is to use an IDE to run a Search and Replace.
  First, replace any `from pymatgen import MPRester` with
  `from pymatgen.ext.matproj import MPRester`. Then, replace
  `from pymatgen import` with `from pymatgen.core import`. Alternatively, if you
  are using a Mac command line, you can do::

    find . -name '_.py' | xargs sed -i "" 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '_.py' | xargs sed -i "" 's/from pymatgen import/from pymatgen.core import/g'

  From a Linux command line, you can do::

    find . -name '_.py' | xargs sed -i 's/from pymatgen import MPRester/from pymatgen.ext.matproj import MPRester/g'
    find . -name '_.py' | xargs sed -i 's/from pymatgen import/from pymatgen.core import/g'

  This should resolve most import errors and only a few more modifications may
  need to be done by hand.

  Specifically, the following "convenience imports" have been removed in favor of
  their canonical import::

    from pymatgen import Composition  # now "from pymatgen.core.composition import Composition"
    from pymatgen import Lattice  # now "from pymatgen.core.lattice import Lattice"
    from pymatgen import SymmOp  # now "from pymatgen.core.operations import SymmOp"
    from pymatgen import DummySpecie, DummySpecies, Element, Specie, Species  # now "from pymatgen.core.periodic_table ..."
    from pymatgen import PeriodicSite, Site  # now "from pymatgen.core.sites ..."
    from pymatgen import IMolecule, IStructure, Molecule, Structure  # now "from pymatgen.core.structure ..."
    from pymatgen import ArrayWithUnit, FloatWithUnit, Unit  # now "from pymatgen.core.units ..."
    from pymatgen import Orbital, Spin  # now "from pymatgen.electronic_structure.core ..."
    from pymatgen import MPRester  # now "from pymatgen.ext.matproj ..."

## v2021.3.3

- __Backwards incompatible__: pymatgen.SETTINGS have been moved to
  pymatgen.settings.SETTINGS. In general, this should not lead to many breakages
  since most of these settings are used within pymatgen itself.
- __Backwards incompatible__: pymatgen.loadfn and get_structure_from_mp have been
  removed since no one was using them.
- critic2_caller has been refactored. (@samblau)
- Improved hash for Composition (@CompRhys)
- Fixes Outcar parsing for VASP 6.2.0. (@MichaelWolloch)
- Allow None for Gaussian functional, bset, charge and multiplicity (@eimrek)

## v2021.2.16

- Add a new interface to OPTIMADE-compliant APIs in pymatgen.ext.optimade (@mkhorton, #2066)
- Addresses missing text file, all_cg.txt, in package
- Note that a previous released increased the suggested minimum numpy version and suggested minimum Python version
- Previous release also dropped support for aconvasp since this the interface has not been maintained

## v2021.2.14

- Misc bug fixes.

## v2021.2.12

- Misc bug fixes.

## v2021.2.8.1

- Patch release to restore `CompositionError` to preserve backwards compatibility.

## v2021.2.8

- Addition of new job types to Q-Chem IO (@espottesmith, #2055),
  note `metal_edge_extender` has been moved into `local_env` for this change
- Improvements to string utils, new Stringify mixin with
  to_pretty_string(), to_latex_string(), to_unicode_string(), to_html_string() (@shyuep)
- Improvements to build system (@shyuep, @ltalirz, see #2046)
- Entry is now immutable, removing "in_place" option for normalize (@mkhorton, @mattmcdermott, #2060)
- Bug fix for co-ordination geometry finder (@davidwaroquiers, #2035)
- Bug fix for GibbsComputedStructureEntry (@mattmcdermott)

## v2021.1.28

- Ability to read Lobster wavefunctions (@JaGeo, #2034)
- Method to estimate number of bands for VASP calculation (@rwoodsrobinson, #2044)
- Q-Chem cube file plotting and improvements to output parsring (@samblau, #2032)
- Improvements to PhaseDiagram hashing and equality checking (@CompRhys, #2014)
- Improvements to pymatgen import speed (@mkhorton, #2031)
- Bug fix for k-path generation (@munrojm, #2037)
- Bug fix for parsing of core potentials from VASP (@utf, #2033)

## v2020.12.31

- End of 2020 release with minor bug fixes for cli scripts.

## v2020.12.18

- New IsayevNN nearest-neighbor algorithm (@utf, #2011)
- Improvements to electrode objects (@jmmshn, #2016)
- Improvements to Element and PhaseDiagram (@jmmshn, #2005)
- Bug fix to increase minimum version of setuptools which was causing incompatible versions of numpy to be installed for some users (@shyuep, see issue #2010)
- Bug fix to VASP run type detection (@rkingsbury, #2007)

## v2020.12.3

- Site insertion algorithm based on charge density (@jmmshn, #1997)
- Allow calculation of Fermi level from occupancies in VASP calculation (@rkingsbury, #2000)
- Improvement to legibility of 3D phase diagram plots (@bayesfactor, #1999)
- Improvement to allow general input for exciting (@vorwerkc, #1975)
- Improvements to code formatting (@mkhorton, #2008)
- Bug fix for VASP run type detection (@rkingsbury, #1996)

## v2020.11.11

- Bug fix for PhononBandStructureSymmLine. (@gpetretto)
- Improved robustness in ABINIT input generation. (@gpetretto)
- Other minor bug fixes.

## v2020.10.20

1. Cp2K support (@nwinner)
2. Better BSPlotter (@fraricci)
3. Better deprecation warnings.
4. Bug fix for Py3.9 support.
5. Bug fix for neutron diffraction get_plot.

## v2020.10.9

- Cube parsing and Cube integration to Bader (@nwinner, #1967)
- Improvements to PhaseDiagram (@CompRhys, #1899)
- Improvements to VASP sets to calculate NGX/Y/Z, NGX/Y/ZF (@jmmshn, #1959)
- Changes to MPRelaxSet, default to low spin for Co (@shyuep, #1976)
- Changes to MPScanSet (@rkingsbury, #1952)
- Rename of `Specie` to `Species`, `Specie` will be retained for backwards compatibility (@shyuep, #1963)
- Bug fix for VASP sets (@utf, #1979)
- Bug fix for PDPlotter (@mattmcdermott, #1973)
- Bug fix for EnergyAdjustment (@rkingsbury, #1960)

## v2020.9.14

- New Plotly backend for PhaseDiagram plotting (@mattmcdermott, #1936)
- New reporting and logging of Materials Project database version in MPRester (@mkhorton, #1945)
- Improvements and bug fixes with mcsqs integration (@rwoodsrobinson, #1942)
- Improvements to PackmolRunner (@rkingsbury, #1947)
- Improvements to ComputerEntry (@rkingsbury, #1948)
- Improvements for MPScanSet (@rkingsbury, #1940)
- Bug fix for Surface and Composition (@gpetretto, #1937)
- Bug fix for EwaldSummation serialization (@lbluque, #1932)
- Bug fix for SeeK k-path (@Ian496, #1930)
- Fix for deprecation warning in MPRester (@rkingsbury, #1951)

## v2020.8.13

- New GibbsComputedStructureEntry (@mattmcdermott, #1921)
- Changes to MPScanRelaxSet and new MPScanStaticSet (@rkingsbury, #1917)
- Changes to LobsterSet (@JaGeo, #1928)
- Bug fix and change for MPRelaxSet (@mkhorton, 9eb3ac2)
- Bug fix for JMolNN (@utf, #1920)
- Bug fix for Element valences (@rkurchin, #1926)
- Bug fix for BabelMolAdaptor (@smheidrich, #1924)
- Bug fix for Gaussion IO (@eimrek, #1918)

## v2020.8.3

- Change neighbor-finding algorithm extension to C instead of C++ for better cross-platform robustness (@chc273)
- Add I/O for JARVIS Atoms (@knc6)

## v2020.7.18

- Add validation and extrapolation for stitching XAS (@yimingcheng)
- Better error handling and possibly verbose warning to get_structure_by_material_id

## v2020.7.16

- Bug fix for boltztrap2 spin support. (@fraricci)

## v2020.7.14

- EwaldSummation is now MSONable (@lbluque).
- Fix for QChem freq parsing (@samblau)
- Much improved linting and workflows.

## v2020.7.10

- Bug fix: serialization of slabs (@utf)
- Bug fix: enumlib url (@wsyxbcl)
- Bug fix: change in tolerance for Lattice comparison (@mbjumar)
- Bug fix: k-path division by zero (@mfherbst)
- New: support for openbabel 3.0 (@orioncohen)

## v2020.7.3

- Make Slabs properly serializable in as_dict. Fixes #1892.
- Fixes for Critic2Caller (@yuuukuma)
- Add cost data for He, H, Ar, Ne, Kr, Tc (@computron)
- Parse scientific notation in OUTCAR (possibly without spaces in between)
- Spin support for boltztrap2 (@fraricci)
- New static method to generate basis functions Lobster (@JaGeo)
- SLME and spillage analysis (@knc6)

## v2020.6.8

- New: Support for parsing WAVECARS with spin-orbit coupling (@mturiansky, #1861)
- New: Support to convert WAVECAR to wannier90 UNK files (@mturiansky, #1861)
- New: Site-weighted XAS spectrum (@yimingchen95, #1837)
- Fixed: Elfcar serialization (@ayushgupta, #1859)
- Fixed: Units in label for phonon plot (@ab5424, #1857)
- Fixed: StructureMatcher serialization (@lbluque, #1850)
- Fixed: Comment string in KPOINTS file (@arosen93, #1842)
- Fixed: parsing of dielectric function in VASP output (@computron, #1836)

## v2020.4.29

- Improved SQS caller. (@rwoodsrobinson)
- VolumetricData speedup (@mturiansk)
- Misc bug fixes

## v2020.4.2

- New high-symmetry k-path algorithm (@munrojm, @kt-latimer)
- New TEM diffraction calculator (@welltemperedpaprika, @thefrankwan, @shyamd)
- New plotly plotting option for Wulff shapes (@CifLord)
- Improvements to SQS caller (@rwoodsrobinson)
- Various bug fixes and improvements (@mfherbst, @chc273,
  @jacksund, @espottesmith, @hongyi-zhao, @montoyjh,
  @dongsenfo, @dynikon) including significant BrunnerNN, EconNN fixes (@utf),
  see individual pull requests for details.

## v2020.3.13

- Added angle_tolerance to CifWriter.
- Change default float precision in CifWriter to 8. Adds float_prec kwarg to
  allow setting of arbitrary precision.
- Rudimentary pymatgen.io.vasp.help.VaspDoc class for obtaining help from VASP wiki.
- Massive documentation cleanup.
- Reorganization of Entry, ComputedEntry (@ayushsgupta).
- Bug fix for PourbaixDiagram (@montoyjh).
- Read WAVECAR from gamma-point only VASP executable. (@bernstei)

## v2020.3.2

- New MonteCarloRattleTransformation and phonopy integration (@utf)
- New structure connectivity features in Chemenv analysis (@davidwaroquiers)
- Bug fixes (@CifLord, @chc273, @JaGeo, @dskoda, @rkingsbury,
  @jmmshn, @espottesmith, @gVallverdu, @yimingchen95, @fraricci)

## v2020.1.28

- Plugin architecture for pymatgen.
- Improvements to pymatgen.analysis.xas.spectrum.XAS class. (@yiming)
- Fixes for ISYM uniform bug and auto-NEDSO (@fraricci)
- Improvements to ReactionDiagram.
- Chemenv improvements (@davidwaroquiers)
- Misc bug fixes.

## v2020.1.10

- New connectivity analysis in Chemenv (@davidwaroquiers)
- Improvements to DOSPlotter (@uthpalah)
- Improvements to writing VASP input sets (@rkingsbury)
- Bug fix for PhaseDiagram (@montoyjh)

## v2019.12.22

- Improvements to reaction calculator (@mattmcdermott)
- VASP input set for SCAN from Materials Project, MPScanSet (@rkingsbury)
- Bug fixes and documentation improvements (@LindaHung-TRI, @rkingsbury, @kwaters4, @rwoodsrobinson, @JaGeo, @nishiyamat, @smheidrich)

## v2019.12.3

- Respect KSPACING in INCAR.
- Bug fixes.

## v2019.11.11

- Extend grosspop class (@Jageo)
- Add option to VaspInputSet to write output with POTCAR.spec
- Add sort_structure option to Poscar.
- Added ability to make gaussian input file without a geometry (@WardLT)
- Misc big fixes.

## v2019.10.16

1. Major refactoring of ABINIT IO to remove workflow-based packages (@gmatteo)
2. Use caching in MinimumVIRENN class. (Alex Ganose)
3. Changes to Lobster module and lobsterset (@jageo)
4. Eigenval object for EIGENVAL output file (@mturiansky)

## v2019.10.4

1. Fix compile args.

## v2019.10.3

- Faster get_all_neighbors based on @chc273's improvements. get_all_neighbors
  now returns a Site-like object with nn_distance, image and index attributes.
  Much easier to use.
- Bug fix for XCrySDen parser (@stevetorr)
- Added optional mid_struct to direct interpolation (@jmmshn)

## v2019.10.2

- IRSpectra class (@henriquemiranda)
- Much faster get_neighbors written in Cython (@chc273).
- VolumetricData allows for sum or subtraction of data with different
  structures, with warnings.

## v2019.9.16

- Updates to annotation, docstrings, etc. Linting service now provided on Github
  Actions as well as CircleCI.

## v2019.9.12

- Massive updates to type annotations, especially for core classes.
- pycodestyle, pydocstyle and mypy will henchforth be enforced for all new PRs.

## v2019.9.8

- Supplemental release to address missing incar_parameters.json

## v2019.9.7

- New fast Pourbaix algorithm (@montoyjh)
- VASP Incar parameter checking (@CifLord)
- New VASP input set for Lobster, read support for GROSSPOP file (@JaGeo)
- New CombinedData class  for LAMMPS (@htz1992213)
- Improvements to molecule fragmenter (@samblau)
- Various bug fixes and improvements (@dongsenfo, @shyuep, @ardunn, @nathan-diodan, @rkingsbury, @kmu)

## v2019.8.23

- pycodestyle now enforced, except on tests. Developers should install
  pycodestyle and the pre-commit hook (copy pre-commit to .git/hooks)
  provided in the repo to check before commits. CI now checks for code style
  and PRs must pass pycodestyle.
- chemsys str input now allowed in get_entries_in_chemsys (@rkingsbury)
- ComputedEntry and subclasses now support a normalize().
- Speed improvements in fragmeter using igraph. (@samblau)

## v2019.8.14

- Update DOSCAR from lobster (@JaGEO)
- PerturbStructureTransformation (@rees-c)
- Misc bug fixes.

## v2019.7.30

- Bug fixes (@shyuep, @mfherbst)
- More type hint annotations (@shyuep)
- Improvements to BabelMolAdaptor (@smheidrich)
- Convenience Transformations for AdsorbateSiteFinder (@mkhorton)

## v2019.7.21

- Add CubicSupercellTransformation and PerturbedSupercellsTransformation (@rees-c, @utf)
- Add interface for ShengBTE (@rees-c, @utf)
- Add interface for Vampire (@ncfrey)
- Improved Lobster interface (@JaGeo)
- Bug fixes (@sthartman, @dwinston, @utf)
- New functionality for calculation of Heisenberg exchange parameters (@ncfrey)
- Improvements to Miller indices handling and Lattice (@CifLord)

## v2019.7.2

- Improvements to grain boundary transformations and Rester (@Tinaatucsd)
- Improvements to AdsorbateSiteFinder (@oxana-a)
- Improvements to Waveder support (@JRSuckert)
- Improvements to run type detection (@darnoceloc)
- Add XAS data to Rester (@yimingchen95)
- Fix to ATAT input/output (@dongsenfo)
- Initial support for Prismatic input (@mkhorton)

## v2019.6.20

- New interface class (@sivonxay, @kylebystrom, @shyamd)
- Updates to SlabGenerator (@CifLord)
- Updates to PiezoTensor (@dongsenfo)
- Add support for parsing on-site density matrix to Outcar (@mkhorton, @mhsiron, @clegaspi)
- Fixes for magnetic space groups (@simonward86)
- Fixes for Lobster class (@JaGeo)
- Fix for FEFF (@stevetorr)
- Fix for Waveder (@JRSuckert)

## v2019.6.5

- Linear scaling get_all_neighbors. Tested to be faster for > 100 atoms (@chc273).
- Lobsterin class to handle input for Lobster (@JaGeo).
- Strict options for composition parsing (@mkhorton).
- Bug fix for CovalentBondNN.get_bonded_structure (@lan496).

## v2019.5.28

- New VASP Input Set "from previous" interface (@utf)
- ELFCAR support (@mkhorton)
- Improvements to plotting of band structures and densities of states (@ShuaishuaiYuan)
- Convenience functions added to Composition including chemical system convention (@mkhorton)
- Various bug fixes (@mkhorton, @utf)
- Improvements to MEGNET API (@shyuep)
- Improvements to Structure interpolation (@mturiansky)

## v2019.5.8

- Numerous updates and improvements to defect classes (@dbroberg)
- New API for MEGNET models, see <http://megnet.crystals.ai> (@shyuep)
- Update to NMR symmeterization (@dongsenfo)
- Change CIF indexing (@kmu)
- Add BoltzTraP mode to NonSCF input sets (@utf)

## v2019.5.1

- Small speeds to Structure.get_all_neighbors.
- Big fixes for gulp_caller. (@kmu)
- Plot fatbands from Lobster. (@jageo)
- Speed up get_ir_mesh (@utf)
- Parsing of plasma frequencies from Outcar.
- Miscellaneous bug fixes.

## v2019.4.11

- Improvements to MimimumDistanceNN (@jmmshn)
- Improvements to Lobster. (@JaGeo)
- Implement a metal warning for ISMEAR < 1 and NSW > 0.
- Misc bug fixes to input sets, including detection of metal systems and
  checking for standardization.

## v2019.3.27

- Bug fixes for OrderDisorderComparator (@utf), custom k-points
in MPNonSCFSet (@dyllamt), battery app (@jmmshn), MPSOCSet (@mkhorton),
more
- Improvements to COHP (@JaGeo)
- Support to read WAVEDER files (@knc6)
- Addition of van Arkel-Ketelaar triangle plots (@CifLord)
- Addition of optional user agent to MPRester API calls, see documentation
for more information (@dwinston)

## v2019.3.13

- Streamlined Site, PeriodicSite, Molecule and Structure code by abandoning
  immutability for Site and PeriodicSite.
- VaspInput class now supports a run_vasp method, which can be used to code
  runnable python scripts for running simple calculations (custodian still
  recommended for more complex calculations.). For example, the following is a
  kpoint convergence script that can be submitted in a queue

.. code-block:: pycon

    from pymatgen import MPRester
    from pymatgen.io.vasp.sets import MPRelaxSet


    VASP_CMD = ["mpirun", "-machinefile", "$PBS_NODEFILE", "-np", "16", "vasp"]


    def main():
        mpr = MPRester()
        structure = mpr.get_structures("Li2O")[0]
        for k_dens in [100, 200, 400, 800]:
            vis = MPRelaxSet(structure,
                user_kpoints_settings={"reciprocal_density": k_dens})
            vi = vis.get_vasp_input()
            kpoints = vi["KPOINTS"].kpts[0][0]
            d = f"Li2O_kpoints_{kpoints}"

            # Directly run vasp.
            vi.run_vasp(d, vasp_cmd=VASP_CMD)
            # Use the final structure as the new initial structure to speed up calculations.
            structure = Vasprun(f"{d}/vasprun.xml").final_structure


    if __name__ == "__main__":
        main()

- Many pymatgen from_file methods now support pathlib.Path as well as strings.
- Misc bug fixes.

## v2019.2.28

- Type hints now available for core classes.
- New pymatgen.util.typing module for useful types.
- Misc bug fixes.

## v2019.2.24

- New EntrySet class for easy manipulation of entries to grab subsets,
  remove non-ground-states, etc. Makes it easier to grab a large set of entries and work with sub chemical systems. Also MSONable for caching.
- Performance improvements in core classes and Poscar (@ExpHP).
- New/changed methods for IcohpCollection and Completecohp

## v2019.2.4

- New Trajectory class for MD simulations (@sivonxay)
- Lattice.get_vector_along_lattice_directions (@blondgeek)
- Misc bug fixes.

## v2019.1.24

- Python 3 only!
- Improvements to local environment code including VESTA bond emulation (@utf)
- Update Cohp analysis (@JaGEO)
- Updates to Boltztrap2 (@fraricci)

## v2019.1.13

- Pymatgen is now Py3 ONLY. If you need Py27 support, please use versions
  < 2019.1.1.
- PARCHG parsing from WAVECAR (@mturiansky)
- Improvements to defect generation algorithms (@dbroberg)
- Simplifications to COHP plotting (@JaGeo)

## v2018.12.12

- Support for IUPAC ordering of elements in Composition formulae (@utf)
- Various bug fixes including returning integer miller indices, catching negative values in Composition and fixes to graph analysis (@utf), fix to Composition serialization (@jmmshen), defect analysis (@HanmeiTang), removing sites in surfaces (@yiming-xu), and fix to support the new PROCAR format in VASP (@dkorotin)
- `PMG_MAPI_ENDPOINT` environment variable added to support different endpoints for the Materials Project REST interface (@mkhorton)

## v2018.11.30

- MPRester.query now supports bulk queries for large scale requests.
  (@dwinston)
- MVLRelax52Set which uses VASP 52 pseudopotentials. (@HanmeiTang)
- EPH calculations in ABINIT (@gmatteo)
- New ScaleToRelaxedTransformation (@CifLord)
- New dimensionality finder, and consolidation of existing algorithms (@utf)
- New dopant predictor built on structure predictor (@utf)
- Misc bug fixes (@HanmeiTang, @utf, @tamuhey, @mkhorton, @yiming-xu, @CifLord)

## v2018.11.6

- Ionic radius based CrystalNN (@computron)
- InterfacialReactivity (@dbroberg)
- Misc bug fixes

## v2018.10.18

- New bond fragmenter and bond dissociation analysis modules (@samblau)
- Improvements to MoleculeGraph (@espottesmith)
- Fix: bug in triclinic tensor conversion to IEEE standard (@montoyjh)
- Fix: insertion battery summary dictionary format (@jmmshn)
- Speed improvements to certain tests (@shyuep, @samblau)

## v2018.9.30

- Fix: increased cut-off to VoronoiNN to avoid scipy crash (@utf)
- Fix: Outcar parsing issues with certain values of electrostatic potential (@sivonxay)
- Fix: bug in EnumlibAdaptor/EnumerateStructureTransformation involving incorrect
  stoichiometries in some instances (#1286) (@shyuep)
- Fix: fractional coordinate finite precision errors in CifParser, now
  also includes additional warnings for implicit hydrogens (@mkhorton)
- New features and improvements to GBGenerator (@ucsdlxg, @shyuep)
- New analysis options in StructureGraph, speed up tests (@mkhorton)
- New utility function to pretty print disordered formulae, along with a
  ordered-to-disordered structure transformation (@mkhorton)
- Ability to use pymatgen's StructureMatcher against AFLOW's library of
  crystallographic prototypes (@mkhorton)
- Make Chgcar serializable to/from dict for database insertion (@jmmshn)

## v2018.9.19

- Fix to composition handling in `MolecularOrbitals` (@dyllamt)
- Fix to allow mixed compressed/uncompressed loading of VASP band structures (@ajjackson)
- New features and fixes to `chemenv` analysis module (@davidwaroquiers)
- Fix to include structure predictor data with pip/conda-installed pymatgen (@shyamd)
- Fixes to `Defect` objects, including allowing rotational supercell transformations (@dbroberg)
- Fix to `BSDOSPlotter` to correctly fill in parts of DOS (@fraricci)
- Added '@' notation parsing in `Composition` (@tamuhey)
- BibTex reference extraction updated in `CifParser` to support ICSD CIFs (@shyamd)
- Various updates to speed up and fix test suite (@shyuep, @fraricci)
- Improvements to BoltzTraP 2 support (@shyuep, @fraricci)

## v2018.9.12

- Use boltztrap2 (@fraricci)
- Refactoring of tensor code to core (@montoyjh)
- Support for new Lobster version (@JaGeo)
- Misc bug fixes

## v2018.8.10

- Bug fix for pymatgen.analysis.gb and pymatgen.io.lammps.

## v2018.8.7

- Massive refactoring of LAMMPS support. (@adengz)
- Allow kwargs passthrough for Structure.to.
- Updates to ABINIT support (@gmatteo)
- GrainBoundaryTransformation class. (@Tinaatucsd)

## v2018.7.15

- Grain boundary generator (Xiangguo Li @ucsdlxg)
- Massive updates to defect code and new DefectTransformation
  (@shyamd)
- Bug fix for OUTCAR parsing with more than one space in
  electrostatic potential.
- get_fermi_interextrapolated to support wider range of
  input doping (@albalu)
- Update to cython compile to support Py3.7.
- Update VoronoiNN cutoff dynamically (@computron)

## v2018.6.27

- Improved local_env and MoleculeGraph (@WardLT, @espottesmith)
- Improve BabelMolAdaptor with conformer search and other functions (@Qi-Max)
- Improved surface analysis (@CifLord)

## v2018.6.11

- Updates to ABINIT support for 8.1.3
- Updates to Interface analyzer.
- Fix bug in deserialization of ComputedStructureEntry.
- Misc bug fixes.

## v2018.5.22

- Misc bug fixes.

## v2018.5.21

- Bug-fix for missing HHI data file.
- Misc bug fixes.

## v2018.5.14

- Dash docs now available for pymatgen. See pymatgen.org "Offline docs" section
  for details.
- Better CrystalNN. (Anubhav Jain)
- Fixes for elastic module. (Joseph Montoya)

## v2018.5.3

- Improvements to qchem (@samblau).
- Improvements to nwchem to support tddft input and parsing (@shyuep).
- Improvements to CrystalNN (@computron).
- Add methods for getting phonon BS, DOS, and DDB output (@dwinston).

## v2018.4.20

- Neutron diffraciton calculator (Yuta)
- Non-existent electronegativity (e.g., He and Ne) are now returned as NaN
  instead of infinity.
- CifParser now handles Elements that are in all caps, which is found in some
  databases. (Gpretto)
- Improvements to local_env (Anubhav Jain)
- Improvements to Qchem ()
- Inputs sets for NMR (Shyam)
- New ChargeDensityAnalyzer class to find interstitial sites from charge density (Hanmei)

## v2018.4.6

- Updated debye temperature formulation (Joey Montoya)
- Add bandgap option for FermiDos for scissoring (Alireza Faghaninia)
- Improved Pourbaix code (Joey Montoya)
- Local env code improvements (Nils)

## v2018.3.22

- Bug fixes to structure, phase diagram module, enumlib adaptor, local env analysis.

## v2018.3.14

- ReactionDiagram for calculating possible reactions between two compositions.
- Misc bug fixes for EnumlibAdaptor and MagOrderingTransformation

## v2018.3.13

- Support for VESTA lattice vector definitions.
- GaussianOutput read now bond_orders of a NBO calculations (@gVallverdu)
- Bug fixes to phonons, abinit support.

## v2018.3.2

- Various algorithms for nearest neighbor analysis (Hillary Pan)
- Cleanup of local_env modules (Nils)
- Enhancements to surface packages (Richard)
- Misc bud fixes

## v2018.2.13

- Improved chemenv parameters and bug fixes (David Waroquiers).
- Improved Qchem IO (Shyam).
- Improved interfacial reactions.
- local_env update (Nils).
- Improved ABINIT support (@gmatteo).
- Misc bug fixes.

## v2018.1.29

- Improvements to local_env (Nils)
- Term symbols for Element (Weike Ye).
- Timeout for enumlib (Horton).

## v2018.1.19

- Phonon plotting and analysis improvements (Guido Petretto).
- Voronoi site finder (Hanmei Tang)
- Some bug fixes for Gaussian (Marco Esters)
- Misc improvements.

## v2017.12.30

- Added detailed Shannon radii information and method.
- Magoms for lanthanides (Weike Ye)
- Chemenv improvements (David Waroquiers)
- Ewald summation improvements (Logan Ward)
- Update to ABINIT support (G Matteo)

## v2017.12.16

- Add optical absorption coefficient method
- Improve plot_element_profile

## v2017.12.15

- Deprecated methods cleanup for 2018. Note that this may break some legacy
  code. Please make sure you update your code!
- Better dielectric parsing for VASP 5.4.4 to include both density-density and
  velocity-velocity functions.
- Orbital-resolved COHPs support (Macro Esters)
- Convenient plot_element_profile method in PDPlotter.
- Input set for SCAN functional calculations.
- Misc bug fixes and code improvements.

## v2017.12.8

- Pymatgen no longer depends on pyhull.
- MPRester method to get interface reaction kinks between two reactants.
- Misc improvements.

## v2017.12.6

- Support for HDF5 output for VolumetricData (CHGCAR, LOCPOT, etc.).
- Support for Crystal Orbital Hamilton Populations (COHPs) (@marcoesters)
- REST interface for Pourbaix data
- Support for optical property parsing in Vasprun.
- Improvements to LammpsData
- Misc bug fixes.

## v2017.11.30

- Fix for severe enumlib_caller bug. This causes enumerations not to be carried
  out properly due to bad accounting of symmetry of ordered sites. It results
  in too few orderings.
- New method to extract clusters of atoms from a Molecule based on bonds.

## v2017.11.27

- Improvements to FEFF
- MPRester now supports surface data.
- Improvement to DiscretizeOccupanciesTransformation.

## v2017.11.9

- Massive rewrite of LAMMPSData to support more functionality (Zhi Deng)
- Misc bug fixes.

## v2017.11.6

- Better exception handling in EnumlibAdaptor and
  EnumerateStructureTransformation.
- Allow bypassing of ewald calculation in EnumerateStructureTransformation.
- get_symmetry_operations API convenience method for PointGroupAnalyzer.
- New DiscretizeOccupanciesTransformation to help automate ordering of
  disordered structures.
- Fix POTCAR check for POSCAR.
- Minor updates to periodic table data.
- Misc bug fixes.

## v2017.10.16

- Added many more OPs and made normalization procedure more robust (Nils Zimmermann)
- Molecular orbitals functionality in Element (Maxwell Dylla)
- Improvements in chemenv (David Waroquiers)
- Add I/O for ATAT’s mcsqs lattice format (Matthew Horton)

## v2017.9.29

- critic2 command line caller for topological analysis (M. Horton)
- Refactor coord_util -> coord.

## v2017.9.23

- Gibbs free energy of a material with respect to Pourbaix stable domains.
- Phonopy io now supports structure conversions.
- EnumerateStructureTransformation now implements a useful occupancy rounding.
- MVLNPTMDSet
- Improved PDPlotter options.
- Misc bug fixes.

## v2017.9.3

- VDW support (Marco Esters)
- Bug fix release.

## v2017.9.1

- Massive refactoring of PhaseDiagram. Now, PDAnalyzer is completely defunct
  and all analysis is carried out within PhaseDiagram itself, e.g.,
  pd.get_e_above_hull as opposed to PDAnalyzer(pd).get_e_above_hull.
- Refactoring of structure prediction. Now in
  pymatgen.analysis.structure_prediction.
- New core Spectrum object and associated pymatgen.vis.plotters.SpectrumPlotter.
- Parsing energies from gen_scfman module in Qchem 5 (Brandon Wood)
- Improvements to LAMMPSData, vasp IO.

## v2017.8.21

- Minor bug fixes.

## v2017.8.20

- Input sets for GW and BSE calculations (Zhenbin Wang) and grain boundary
  calculations (Hui Zheng). Input sets now support overriding of POTCAR
  settings.
- Haven ratio calculation (Iek-Heng Chu).
- LAMMPS io updates (Kiran Matthews).
- Oxidation state guessing algorithms based on ICSD data (Anubhav Jain).
- New local_env module for local environment analysis. (Nils Zimmerman).
- pymatgen.util.plotting.periodic_table_heatmap (Iek-Heng Chu).
- Improvements to surface code for tasker 3 to 2 reconstructions.
- pymatgen.analysis.interface_reactions.py for analyzing interfacial reactions
  (Yihan Xiao).

## v2017.8.16

- PointGroupAnalyzer now allows for symmetrization of molecules. (@mcocdawc)
- QuasiharmonicDebyeApprox with anharmonic contribution. (Brandon)
- Improvements to LAMMPS io. (Kiran)
- Misc bug fixes.

## v2017.8.14

- Fixes and minor improvements to elastic, bader and defect analyses.

## v2017.8.4

- Major refactoring and improvements to lammps io. (Kiran)
- Major improvements to BaderAnalysis. (Joey and Zhi)
- Major improvements to Magmom support in cifs, SOC calculations, etc.
  (Matthew Horton)
- Add remove_site_property function. Add magmom for Eu3+ and Eu2+.
- BoltztrapAnalyzer/Plotter support for seebeck effective mass and complexity
  factor (fraricci)

## v2017.7.21

- Misc bug fixes to elastic (J. Montaya),
- Decrease default symprec in SpacegroupAnalyzer to 0.01, which should be
  sufficiently flexible for a lot of non-DFT applications.

## v2017.7.4

- Bug fixes for oxide corrections for MP queried entries, and pickling of Potcars.
- Default to LPEAD=T for LEPSILON=T.

## v2017.6.24

- New package pymatgen.ext supporting external interfaces. Materials Project
  REST interface has been moved to pymatgen.ext.matproj. Backwards compatibility
  will be maintained until 2018.
- Two new interfaces have been added: i) Support for John Hopkin's Mueller
  group's efficient k-point servelet (J Montaya). ii) Support for
  Crystallography Open Database structure queries and downloads. (S. P. Ong).
  See the examples page for usage in getting structures from online sources.

## v2017.6.22

- Speed up pmg load times.
- Selective dynamics parsing for Vasprun (Joseph Montaya)
- Allow element radius updates in get_dimensionality (Viet-Anh Ha).
- Dielectric function parse for vasp 5.4.4 (Zhenbin Wang).
- Parsing for CIF implicit hydrogens (Xiaohui Qu).

## v2017.6.8

- Switch to date-based version for pymatgen.
- Electronegativities now available for all elements except for He, Ne and
  Ar, which are set to infinity with a warning.
- Bond lengths are now set to sum of atomic radii with warning if not available.
- Bug fixes to boltztrap, symmetry for trigonal-hex systems, etc.

## v4.7.7

- Magnetic symmetry and CIF support. (Horton)
- Improved PWSCF Input file generation.
- Misc bug fixes.

## v4.7.6

- Fix serious bug in PointGroupAnalyzer that resulted in wrong point groups assigned to non-centered molecules.
- Useful get_structure_from_mp at the root level for quick retrieval of common structures for analysis.
- More efficient kpoint grids.
- Misc bug fixes.

## v4.7.5

- MultiXYZ support (Xiaohui Xu)
- Misc bug fixes and cleanup.

## v4.7.4

- New ferroelectric analysis module (Tess).
- Magmom support and MagSymmOp (Matthew Horton).
- Improved CIF Parsing.

## v4.7.3

- Sympy now a dependency.
- Massive improvements to elastic package. (Joseph Montoya)
- Symmetrized structures now contain Wyckoff symbols.
- More robust CIF parsing and MITRelaxSet parameters (Will).

## v4.7.2

- Support for Abinit 8.2.2, including support for DFPT calculations. (Matteo)

## v4.7.1

- Pathfinder speedup
- Minor bug fix for plots.

## v4.7.0

- Improvements to BSDOSPlotter.
- Enhancements to Phase diagram analysis and reaction calculator.
- Enhancements to surface slab and adsorption. (Richard and Joey)
- Support NpT ensemble in diffusion analysis.

## v4.6.2

- Improve Spacegroup class support for alternative settings. Add a get_settings class method.
- Improvements to FEFF support.
- Improvements to EOS class.

## v4.6.1

- Phonon bandstructure plotting and analysis. (Guido Petretto)
- New capabilities for performing adsorption on slabs. (Joey Montoya)
- Remove pathlib dependency.

## v4.6.0

- Improve support for alternative settings in SpaceGroup.
- Fix respect for user_incar_settings in MPNonSCFSet and MPSOCSet
- Support for argcomplete in pmg script.
- Speed ups to Ewald summation.
- Add functionality to parse frequency dependent dielectric function.
- Improvements to Bolztrap support.

## v4.5.7

- PMG settings are now prefixed with PMG_ to ensure proper namespacing.
- Improve error output in command line bader caller.
- Add Py3.6 classifier.
- Misc bug fixes.

## v4.5.6

- Minor bug fix.
- Fixed elastic energy density

## v4.5.5

- Fix bad reading of pmgrc.
- Gaussian opt section added allowing for torsion constraints
- Update spglib.

## v4.5.4

- BSDOSPlotter (Anubhav Jain)
- Fixes to defect analysis (Bharat)
- intrans as an input to BoltztrapAnalyzer. Allows for scissor operation.
- Pmg is now continuously tested on win-64/py35 using Appveyor!

## v4.5.3

- Added an alternative interstitial finder that works with a grid-based structure-motif search. (Nils Zimmermann)
- Optional possibility to specify that the saddle_point in the NEB should have a zero slope. (David Waroquiers)
- Read intensity and normal modes for Gaussian. (Germain Salvato Vallverdu)
- Minor bug fixes.

## v4.5.2

- Minor bug fix for POTCAR settings.

## v4.5.1

- You can now specify a different default functional choice for pymatgen by
  setting PMG_DEFAULT_FUNCTIONAL in .pmgrc.yaml. For use with newer
  functional sets, you need to specify PBE_52 or PBE_54 for example.
- Switch to ISYM 3 by default for HSE.
- Updates to FEFF>
- Misc bug fixes and startup speed improvements.

## v4.5.0

- Major speed up of initial load.
- Collection of misc changes.

## v4.4.12

- Fix for dynamic numpy import.

## v4.4.11

- Update to new version of spglib.

## v4.4.10

- Minor fixes for proper gzipped structure file support and MVLSlabSet.

## v4.4.9

- Dependency cleanup. Now, basic pymatgen requires on much fewer
  packages.
- Fixed reading of POSCAR files when more than 20 types of atoms.
- Misc bug fixes.

## v4.4.8

- Cleanup of entry points and dependencies.

## v4.4.7

- Update to spglib 1.9.7.1
- Proper use of dependency markers for enum34.

## v4.4.6

- Update to spglib 1.9.6, which fixes some bugs and is Windows compatible.

## v4.4.5

- Bug fix for SubstitutionProb.

## v4.4.4

- Bug fix for electronic structure plotter.

## v4.4.3

- Bug fix for Diffusion Analyzer.

## v4.4.2

- Bug fix for BS serialization.
- Cleanup dependencies.

## v4.4.1

- Massive updates to FEFF support (Kiran Mathews).
- Bug fixes in band structure plotting.

## v4.4.0

- Much more Pythonic API for modifying Structure/Molecule species. Now,
  strings, slices, and sequences should magically work, in addition to the
  previous API of simple int indices. Examples::

    s[0] = "Fe"
    s[0] = "Fe", [0.5, 0.5, 0.5]  # Replaces site and fractional coordinates.
    s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}  # Replaces site and fractional coordinates and properties.
    s[(0, 2, 3)] = "Fe"  # Replaces sites 0, 2 and 3 with Fe.
    s[0::2] = "Fe"  # Replaces all even index sites with Fe.
    s["Mn"] = "Fe"  # Replaces all Mn in the structure with Fe.
    s["Mn"] = "Fe0.5Co0.5"  # Replaces all Mn in the structure with Fe: 0.5, Co: 0.5, i.e.,creates a disordered structure!

- Massive update to internal representation of Bandstructure objects for
  memory and computational efficiency.
- Bug fixes to CIF parsing in some edge cases. (Will Richards).

## v4.3.2

- Massive speedup of Bandstructure, especially projected band structures,
  parsing.
- Massive update to pmg cli script, with new query functions as well as a
  more rational command structure.
- Updates to ChemEnv.
- Misc bug fixes.

## v4.3.1

- Upgrade monty and spglib requirements for bug fixes.
- Updates to feff support (Kiran).

## v4.3.0

- Massive update to elastic module. (Joey Montaya)
- Pathfinder algorithm for NEB calculations. (Ziqing Rong)
- Wheels for Windows and Mac Py27 and Py35.

## v4.2.5

- Bug fix for BSPlotter.

## v4.2.4

- Bug fix for kpoint weight calculation for Monkhorst meshes.

## v4.2.3

- Minor cleanup.
- Simplified installation. enumlib and bader can now be installed using pmg setup --install.

## v4.2.2

- Global configuration variables such as VASP\_PSP\_DIR and MAPI\_KEY are now
  stored in "~/.pmgrc.yaml". If you are setting these as environmental
  variables right now, you can easily transition to the new system using::

      pmg config --add VASP_PSP_DIR $VASP_PSP_DIR MAPI_KEY $MAPI_KEY

  This new scheme will provide greater flexibility for user-defined
  global behavior in pymatgen, e.g., tolerances, default input sets for
  transmuters, etc., in future.
- Beta of k-point weight calculator.
- Use default MSONable as and from_dict for all transformations.

## v4.2.1

- New DopingTransformation that implements an automated doping strategy.
- Updated MIC algorithm that is a lot more robust (Will Richards).
- Major update to chemenv package (David Waroquiers)

## v4.2.0

- Fix important bug in minimum image distance computation for very skewed cells.
- Major refactoring of WulffShape code.
- Misc bug fixes for elastic tensor and other codes.

## v4.1.1

- Major refactoring of WulffShape and lammps support.

## v4.1.0

- Wulff shape generator and analysis.
- Minor bug fixes.

## v4.0.2

- Fix kpoint reciprocal density.

## v4.0.1

- Minor bug fix release.

## v4.0.0

- Massive update with many deprecated methods removed. Note that this
  may break backwards incompatibility!
- Support for ABINIT 8.
- Improved sulfide compatibility.

## v3.7.1

- Fix deprecation bug.

## v3.7.0

- Last version before pymatgen 4.0, where deprecated modules will be removed!
- Massive update to LAMMPS (Kiran Matthews).
- New input sets with a different interface that replaces old input sets.
- Massive update to elastic properties.

## v3.6.1

- Massive cleanup to Boltztrap interface (Anubhav Jain)
- Refactor of piezoelectric analysis to use tensor base class (Joey)
- More robust CIF parsing.

## v3.6.0

- Pymatgen now uses spglib directly from Togo's website. Spglib is no longer
  bundled as a dependency.
- Improved support for velocities in Poscar (Germaine Vallverdu)
- Backwards incompatible change in Born charge format in Outcar.
- Fixes for Lammps input serialization

## v3.5.3

- Misc refactorings and bug fixes, especially for Outcar and Boltztrap classes.

## v3.5.2

- Minor update to DerivedInputSet interface.

## v3.5.1

- New derived input sets for generating inputs that depende on previuos
  calculations. Old input sets deprecated.

## v3.5.0

- Chemical environment analysis package (David Waroquiers).
- Piezoelectric property analysis (Shayam).
- Cythonize certain expensive core functions. 5-10x speedup in large structure matching (Will Richards).
- New NMR parsing functionality for Outcar (Xiaohui Qu).
- Improved io.lammps (Kiran Mathews).
- Update to spglib 1.9.2.
- Element properties now return unitized float where possible.
- Bug fix for get_primitive_standard affecting rhombohedral cells (important for band structures).
- Vasprun.final_energy now returns corrected energy with warning if it is different from final electronic step.

## v3.4.0

- 10-100x speed up to Structure copying and Site init, which means many
  functionality has seen significant speed improvement (e.g., structure
  matching).
- Convenience method Structure.matches now perform similarity matching
  for Structures.
- Bugfix for band gap determination.

## v3.3.6

- Update to use enum.x instead of multienum.x.
- Minor robustness fixes to VaspInputSet serialization.
- Add a reciprocal density parameter to vasp sets.
- Minor bug fixes to Vasprun parsing.

## v3.3.5

- StructureMatcher can now work with ignored species.
- Added interpolation failure warnings and smooth tolerance for
  scipy.interpolate.splrep in bandstructures (Tess).
- Added DiffusionAnalyzer.get_framework_rms_plot.
- Complete rewrite of Procar class to use NDarray access and zero-based
  indexing.
- OrderParameters class for analysis of local structural features
  (Nils Zimmermann).
- Bug fixes for Procar, MPRester and SpaceGroup 64.
- Added Github templates for contributing to pymatgen.

## v3.3.4

- Procar now supports parsing of phase factors.
- Miscellaneous bug fixes.

## v3.3.3

- Bug fixes for Poscar.
- Fix Kpoints pickling.

## v3.3.2

- Bug fixes for pymatgen.io.abinit
- Other minor big fixes.

## v3.3.1

- Minor bug fix release for pickle and elastic constants.

## v3.3.0

- Updated and checked for Python 3.5.* compatibility.
- Element, Spin, Orbital and various other Enum-like classes are now actually
  implemented using Enum (with enum34 dependency for Python < 3.4).
- Speed up Site creation by 20% for ordered sites, with cost in terms of
  slightly slower non-ordered Sites. Since ordered Sites is the far more common
  case, this gives significant boost for large scale manipulations of
  structures.
- Alternative, more pythonic syntax for creating supercells via simply
  Structure _3 or Structure_ (3, 1, 1).
- zeo++ fixes.
- More stable incar settings for MITMDVaspInputSet.

## v3.2.10

- Fix missing scripts
- Improvements to units module.
- Speed up EwaldSummation.

## v3.2.9

- Major PD stability improvements, especially for very high dim hulls with lots
  of entries.
- Improvements to Ewald summation to be close to GULP implementation.
- Deprecate physical constants module in favor of scipy's version.
- Remove many pyhull references to use scipy's ConvexHull implementation.
- Bug fix for sulfide correction.

## v3.2.8

- Make pyhull optional.
- Sulfur correction added to MaterialsProjectCompatibility for more accurate
  sulfide formation energies.
- ADF io support. (Xin Chen)
- Bug fixes for spacegroup subgroup testing.

## v3.2.7

- Add warning for limited subgroup testing functionality in Spacegroup.

## v3.2.6

- Extensive support for elasticity tensor analysis (Joseph Montoya).
- Misc bug fixes and performance improvements.
- Add support for QChem4.3 new format of Batch jobs

## v3.2.5

- Improved potcar setup via "pmg setup", with MAPI setup.
- Support for new POTCARs issued by VASP.
- Improvements to ABINIT support.
- Improvement to Boltztrap support, e.g., scissor band gap, etc.
- Vasprun now issues warning when unconverged run is detected.

## v3.2.4

- GaussianOutput can now parse frequencies, normal modes and Cartesian forces
  (Xin Chen).
- Support for Aiida<->pymatgen conversion by the Aiida development team (Andrius
  Merkys).
- Specialized BSVasprun parser that is ~2-3x faster than Vasprun.
- Refactor the boltztrap package (merge a few methods together) and add several
  new methods (power factor, seebeck...)
- Support of the new PCM format in QChem 4.3
- Local environment analysis to pmg script.
- Deprecate prettytable in favor of tabulate package.
- Improvements to MITNEBVaspInputSet.
- Misc bug fixes.

## v3.2.3

- Massive update to abinit support. Note that pymatgen.io.abinitio has
  been refactored to pymatgen.io.abinit. (Matteo, Setten)
- NwOutput now supports parsing of Hessian matrices (contributed by Xin
  Chen)
- Gaussian support now has the ability to read potential energy surface
  and electronic transitions computed with TD-DFT (Germain Salvato
  Vallverdu)
- Bug fixes for CifWriter with symmetry.
- Bug fixes for surface generation and reactions.
- Monty requirement increased.

## v3.2.1

- Fix wrong U value for Ce and Eu.
- Properly handle empty multiline strings in Cif
- Add ability to get specific data in MPRester.get_entries. Make all get_entry
  methods consistent  in kwargs.

## v3.2.0

- Force conversion to an actual list in selective dynamics and velocities in
  Poscar.
- fix small bug in BSPlotter (wrong ylim)
- Elastic tensor parsing in Outcar

## v3.1.9

- Fix scripts.

## v3.1.7

- Bug fixes for MPRester.
- Ensure correct monty version requirement in setup.py.

## v3.1.6

- Rudimentary PWSCF output reading.
- Fix ASE support.
- Support for WAVEDERF and reading multiple dielectricfunctions in vasprun.xml.
  (Miguel Dias Costa)

## v3.1.5

- Move vasp.vasp*put to vasp.*puts. Also, maintain backwards compatibility with
  vaspio.vasp_*put

## v3.1.4

- Fix missing yaml files that have been moved.

## v3.1.3

- Major refactoring of pymatgen.io. Now, the io suffix is dropped from all io
  classes. i.e., it is just pymatgen.io.vasp, not pymatgen.io.vaspio. Also, all
  input sets have been moved within the relevant package, e.g.,
  pymatgen.io.vasp.sets. All changes are backwards compatible for now. But
  deprecation messages have been included which states that the stubs will be
  removed in pymatgen 4.0. Pls migrate code when you see the deprecation
  messages.
- Make Composition.anonymized_formula truly chemistry independent (No A2B2
  for peroxides or A2 for diatomic gasses)
- Allowing CIF data_* header to be prefixed with spaces and tabulations.

## v3.1.2

- HHI Resource Analysis (by Anubhav Jain).
- Bug fixes for surfaces normalizatino.
- Bug fix for Vasprun parsing of response function keys.
- Dockerfile for generation of an image for pymatgen.
- Updated requirements.txt for latest requests, scipy, numpy.

## v3.1.1

- Bug fixes for SpacegroupAnalyzer and SlabGenerator.
- Much faster normal vec search.

## v3.1.0

- Much improved surface generation algorithm that provides for
  orthogonality constraints.
- Transition state analysis tools! (beta)
- Massive improvements in Outcar parsing which provides a powerful grepping
  syntax.
- PWSCFInput generation (beta).
- Reduce default SIGMA to 0.05 for MP input sets.
- Update spglib to 1.7.3 as per recommendation of Togo.
- Many bug fixes and efficiency improvements.

## v3.0.13

- Bug fix for parsing certain types of CIFs.
- MPRester now has get_materials_id_references helper method.
- Minor fix for Vasprun.final_energy.
- Added mp_decode option to MPRester.query to allow option to not decode into
  pymatgen objects.
- New POTCAR hash scheme to more robustly identify unique POTCARs.
- Link to <http://bit.ly/materialsapi> for information on Materials API
  document schema for use with MPRester.query method.

## v3.0.11

- Lots of abinitio improvements (Matteo).
- Added mp_decode option to MPRester.query to allow option to not decode into pymatgen objects.

## v3.0.10

- Fix Cartesian coord parsing in Poscar class.
- Vasprun now works with non-GGA PBE runs
- Misc bug fixes

## v3.0.9

- Major bug fixes for CIF parsing (Will Richards).
- Support for {Li,Na} syntax in parse_criteria for MPRester.
- Additional example notebook for ordering and enumeration.
- More robust checking for oxidation states in EnumerateStructureTRansformation.
- Improvements to Slab polarity checking.

## v3.0.8

- Massive update to abinitio (Matteo).
- Improvements to OUTCAR parsing (Ioannis Petousis).

## v3.0.7

- Powerful Slab generation algorithms (beta!).
- Improvements to DiffusionAnalyzer with constant smoothing option.
- Significantly improve look of DOS plots using prettyplotlib.

## v3.0.6

- Cost analysis module (Anubhav Jain)
- More Py3k fixes.
- Extensive abinitio updates (Matteo).

## v3.0.5

- Completely revamped symmetry package. The finder.SymmetryFinder and
  pointgroup and spacegroup modules are now deprecated. Instead,
  all symmetry analysis is in the :module:`pymatgen.symmetry.analyzer`_
  module. There is also a completely rewritten support for symmetry groups in
  :module:`pymatgen.symmetry.groups`_. Structure now supports a static
  constructor to generate a structure from a spacegroup (see examples).
- BatteryAnalyzer class (Anubhav Jain) to provide for some common analysis of
  intercalation electrodes.
- Minor bug fixes for structure_matcher, lattice, abinitio.
- MOAB qadapter for abinit. (Liam Damewood)

## v3.0.4

- Fix missing structures json data.

## v3.0.3

- Updates to DiffusionAnalyzer for more fine-grained statistics.
- Bug fixes and tweaks to linear assignment
- Improved PymatgenTest class which provides a suite of test structures.
- Speedups to Phase Diagram
- Lots of improvements to Gaussian support (Nicolas Dardenne) and Abinit IO
  (Matteo).
- Lots of Py3k minor updates.
- Updated doc for Diffusion analyzer. Invert sq_disp_ions for more intuitive handling.

## v3.0.2

1. Consistent use of unicode throughout pymatgen.
2. Minor bug fixes.

## v3.0.1

1. Minor bug fixes for cifio.
2. Py3k updates for abinitio.

## v3.0.0

- Pymatgen is now completely Python 2.7 and Python 3.x compatible!
- Spglib and pyhull have been updated to support Python 3.x.
- Completely rewritten pure python cifio module (courtesy of William Davidson
  Richards) removed dependency on PyCIFRW, which has been causing many issues
  with installation.
- Structure and Molecule now supports a very convenient to() and from_str and
  from_file functionality. Instead of trying to load the appropriate parser,
  you can output and read from the appropriate formats directly. See example
  usage.
- ~50% speedup to LinearAssignment code.
- Continuous integration and contribution guidelines now include Python 3.
- **Backwards incompatible changes**
- matgenie.py has now been renamed simply "pmg" for brevity.
- All deprecated methods in pymatgen 2.x have been removed. E.g.,
  pymatgen.core.structure_modifier is no longer available.
- Pymatgen classes now uses the as_dict() method protocol implemented in the
  Monty package. The to_dict property is now deprecated and will be removed
  in pymatgen v3.1.
- Update main docs page examples with the new Structure to, from formats.

## v2.10.6

- Bug fix for np1.9 incompatibility. Now works.
- Use wheel for pymatgen deployments.
- matgenie.py is now renamed to pmg for faster CLI usage.
- Improvements to KPOINTS automatic generation.
- Simpler and faster Structure.get_all_neighbors

## v2.10.5

- DiffusionAnalyzer now has non-smoothed option.
- Kpoints generation algorithm now guarantees minimum # of points.
- Compatibility now has a proper explanation dict.
- Vaspruns with NSW == 1 now checked properly for electronic conv.
- make_movie now supports kwargs.

## v2.10.3

- MPRester.query now supports a simple but powerful string criteria syntax
  with support for wild cards.
- Improvements to Composition - support for negative compositions, sorting etc.
- Speed ups to StructureMatcher.

## v2.10.2

- Bug fix for Projected DOS parsing in new Vasprun.
- Compatibility now has a _explain_ method which provides a detailed outline
  of the changes that a Compatibility makes to an Entry.

## v2.10.1

- Minor fix for monty requirements in setup.py.

## v2.10.0

- Major update: MPRester now uses Materials API v2! Also major refactoring
  of MPRester.
- Vastly improved Vasprun parser using cElementTree. Twice as fast,
  half as much code and easier to maintain.
- Vast improvements to Qchem functionality (Xiaohui Qu).
- Improved handling of Structure manipulations for extremely large
  structures (particularly in terms of memory consumption).
- Bug fix for XYZ parsing for scientific notation.
- Improve monty.serialization for transparent handling of JSON vs YAML.
  Requirements updated to monty>=0.3.3.
- Update numpy requirements to 1.8+. Fixes memory leak.
- Other minor bug fixes.

## v2.9.14

- Implements Structure.sort method. Both Structure.sort and the
  get_sorted_structure methods now supports all arguments supported by list
  .sort().
- VaspInputSets configs, as well as several other configs now uses yaml. Note
  the new dependency on pyyaml. It is highly recommended that you install
  pyyaml with the libyaml C bindings.
- Fix missing spglib dependency.
- Use monty.serialization for transparent handling of JSON vs YAML.
  Requirements updated to monty>=0.3.1.

## v2.9.13

- Urgent bug fix for missing compatibility yamls.

## v2.9.12

- Defect transformations (Bharat).
- Support for optical properties (Geoffroy Hautier and David Waroquiers).
- Improved support for some VASP output files (XDATCAR and OSZICAR).
- Refactored compatibilities now uses YAML for ease of reading.

## v2.9.11

- Bug fix for get_xrd_plot.
- Speed up XRD calculator by allowing specification of two theta ranges.
- Minor improvements to Gulp caller.

## v2.9.10

- Bug fix for unequal coefficients sizes in XRD.
- Support for Ag radiation in XRD calculator.
- Improved Procar class for extraction of information. (Germain Salvato
  Vallverdu)
- Bug fix for extraction of GGA data from Materials API.

## v2.9.9

- XRDCalculator now supports disordered structures.
- Minor speed ups and improvements.

## v2.9.8

- Initial beta version of XRD pattern calculator.
- Pymatgen now uses spglib 1.6.0.
- Update to Vasprun to compute static deilectric constants with DFPT in VASP.
  (Geoffroy Hautier)

## v2.9.7

- Quick bug-fix release that provides a better solution to Structure handling
  of properties instead of sanitizing MPRester structures.

## v2.9.6

- Patch to allow 1D phase diagrams (essentially finding the lowest energy
  phase).
- Better error checking for Bandstructure KPOINTs.
- Patch to sanitize structures obtained from MPRester.

## v2.9.5

- Bug fix for linear assignment, which may sometimes affect Structure
  Matcher results.
- Minor improvement to the way grand canonical PDs work.

## v2.9.4

- Bug fix for Pourbaix Maker (Sai).
- Streamline use of scratch directories for various calls. Require monty >=
  0.1.2.
- High accuracy mode for Zeo++ (Bharat Medasani).

## v2.9.3

- Bug fix release for printing TransformedStructures from Substitutor (Will
  Richards).
- Misc improvements in BVAnalyzer, coord_utils and defects (Will Richards,
  David Waroquiers and Bharat Medasani).

## v2.9.2

- Bug fix release for DummySpecie, which failed when deserializing from
  json and had bad hash function.

## v2.9.1

- Structure/Molecule now supports Pythonic list-like API for replacing and
  removing sites. See :ref:`quick_start` for examples.

## v2.9.0

- Updates to support ABINIT 7.6.1 (by Matteo Giantomassi).
- Vastly improved docs.
- Refactoring to move commonly used Python utility functions to `Monty
  package <https://pypi.python.org/pypi/monty>`_, which is now a dependency
  for pymatgen.
- Minor fixes and improvements to DiffusionAnalyzer.
- Many bug fixes and improvements.

## v2.8.10

- Refactoring of qchemio module (by Xiaohui Qu).

## v2.8.9

- qchemio module (by Xiaohui Qu).

## v2.8.8

- Minor bug fix release for Structure species substitution methods.

## v2.8.7

- Massive update to pymatgen.io.abinitio package (by Matteo Giantomassi).
- Bug fixes for StructureMatcher's group_structure.
- Misc bug fixes and cleanup.

## v2.8.6

- Bug fix for VASP io set introduced by the default to sorting of structure
  sites when generating VASP input.

## v2.8.4

- Completely revamped Compatibility/Correction system which improves
  readability (Shyue Ping Ong/Anubhav Jain/Sai Jayaraman). This change is
  backwards compatible for the most part.

## v2.8.3

- Big fix release for json dumping for unitized floats.

## v2.8.2

- Bug fix release to improve CIF parsing for more non-standard CIF files.
  In particular, non-ascii characters are removed and _cgraph* fields are
  removed prior to parsing for better support in PyCiFRW.

## v2.8.1

- Bug fix release. Incorrect units assigned for ionic radii.
- Improved nwchemio supports COSMO and ESP calculations (Nav Rajput).

## v2.8.0

- __Units__. Pymatgen now has a new system of managing units,
  defined in pymatgen.core.units. Typical energy, length, time,
  temperature and charge units are supported. Units subclass float,
  which makes the usage transparent in all functions. The value that they
  being are in terms of conversions between different units and addition and
  subtraction of different units of the same type. Some basic quantities
  like ionic radii and atomic masses are now returned in unitized forms for
  easy conversion. Please see :mod:`pymatgen.core.units` and the
  :doc:`examples </examples>` for a demonstration of house to use units in
  pymatgen.
- __Minor backwards-incompatible change__. Structures are now sorted by
  default when generating VASP input files using vaspio_set. Old behavior can
  be obtained by setting sort_structure=False in the constructor. This is
  typically the desired behavior and prevents the generation of large
  POTCARs when atomic species are not grouped together.
- Bug fix for Molecule.substitute. Earlier algorithm was not detecting
  terminal atoms properly.
- Additional conversion tools for ABINIT (by Matteo Giantomassi).

## v2.7.9

- Minor bug fix release to fix pyhull dependencies to be more friendly.
- Improved structure matcher that allows for more flexible matching. New
  matching between ordered and disordered comparator.

## v2.7.7

- Beta new Gulp Caller and Zeo++ interface classes (Bharat . Zeo++ is an open
  source software for performing high-throughput geometry-based analysis of
  porous materials and their voids. Please see
  <http://www.maciejharanczyk.info/Zeopp/about.html>.
- Specify version of distribute to 0.6.34 for better compatibility.

## v2.7.6

- Support for VTK 6.x in structure visualization.
- Updated install instructions for openbabel.
- Preliminary pourbaix analysis (Sai Jayaratnam).

## v2.7.5

- Vastly improved Nwchem IO (by Shyue Ping Ong).
- Much improved ABINIT support (by Matteo Giantomassi).

## v2.7.4

- Added basic Nwchem (<http://www.nwchem-sw.org/>) IO support. (by: Shyue Ping
  Ong).
- New MoleculeMatcher class for comparing molecules by RMS. Requires
  openbabel with python bindings. (by: Xiaohui Qu)
- New functional group substitution capability for molecules (by: Lei Cheng
  and Shyue Ping Ong).

## v2.7.2

- Minor bug fix release to fix some rare errors in very high dimensional
  phase diagrams. **Requires new pyhull version (1.3.8).**

## v2.7.1

- __Major backwards-incompatible change.__ With effect from v2.7.1,
  the default Structure and Molecule classes are now _mutable_ objects. All
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
- read_structure and write_structure now supports pymatgen's JSON-serialized
  structures.
- read_mol and write_mol functions now available (analogues of
  read_structure and write_structure for molecules)

## v2.7.0

- Beta support for ABINIT input and output via pymatgen.io.abinitio
  (courtesy of the excellent work of Matteo Giantomassi).
- Properties are now checked when comparing two Species for equality.
- MaterialsProjectVaspInputSet is now renamed to MPVaspInputSet for easier
  typing. The old input sets have been deprecated.
- New VaspInputSets for MPStatic, MPNonSCF, MITMD which supports uniform
  grid, bandstructure and molecular dynamics calculations. The MD input set
  uses MIT parameters for speed.
- A beta DiffusionAnalysis class in the apps package.
- A revised KPOINT grid algorithm that generates more reasonable meshes.
- A guided install script is now provided for Mac and Linux users.

## v2.6.6

- Updates to feffio (credit: Alan Dozier)
- Added detailed installation instructions for various platforms.
- Support for charge and spin multiplicity in Molecule. Expanded methods
  available in Molecule.
- Added supercell matching capabilities to StructureMatcher.
- More robust creation of PhaseDiagrams to take into account potential qhull
  precision errors.

## v2.6.5

- Added a command_line caller to do Bader charge analysis using Henkelmann
  et al.'s algorithm.
- Bug fix for POSCAR parsing when title line is an empty string.
- Added __rmul__ operator for Composition.
- Vastly expanded available aliases.

## v2.6.4

- Bug fixes for selective dynamics in Poscar.
- Improved Procar parsing to support both simple and detailed PROCARs.

## v2.6.3

- Added new MaterialsProject REST interfaces for submit/query/delete_snl
  (currently open in beta for collaborators only).
- Added support for new MaterialsProject REST method get_stability.
- Added aliases for PhaseDiagram, GrandPotentialPhaseDiagram,
  PDAnalyzer and PDPlotter in pymatgen.phasediagrams.
- Improvements to StructureMatcher: stol (site - tolerance) redefined as
  a fraction of the average length per atom. Structures matched in fractional
  space are now also matched in Cartesian space and a rms displacement
  normalized by length per atom can be returned using the rms_dist method.

## v2.6.2

- Site and PeriodicSite now uses a Composition mapping type to represent
  the species and occupancy, instead of a standard dict.
- Bug fix for reading and re-writing out of Potcars.
- VaspInputSet now supports MSONable framework.
- Strain cell option in StructureEditor.
- Miscellaneous bug fixes and speedups.

## v2.6.1

- Use requests.Session in MPRester for connection pooling and code simplicity.
- Support for "with" context manager in MPRester.
- Updated periodic table data to correct errors in Ru, Tc and other elements.
- New methods in Lattice to obtain Wigner-Seitz cell and Brillouin Zone.
- Miscellaneous bug fixes and speedups.

## v2.5.5

- Bug fix release for cifio for rhombohedral structures.
- Miscellaneous bug fixes and speedups.

## v2.5.4

- Vastly improved Gaussian input file parsing that supports more varieties
  of input specifications.
- StructureNL now supports molecules as well as structures.
- Updated atomic and vdw radius for Elements.
- Miscellaneous bug fixes and speedups.

## v2.5.3

- Bug fix for StructureNotationalLanguage.
- Support for LDA US potential. matgenie.py script option to generate POTCARs.
- Beta version of StructureNotationLanguage, a markup format for Structure
  data with metadata such as authors and references. (Anubhav Jain)
- Vasprun parsing now parses dielectric constant where available. (Geoffroy
  Hautier)
- New custom ipython shell script for pymatgen.
- Miscellaneous bug fixes and speedups.

## v2.5.1

- Bug fixes for primitive cell finder.
- Remove deprecated use_external_qhull option in PhaseDiagram classes.
- Miscellaneous bug fixes and speedups.

## v2.5.0

- Added optimization package with linear assignment class.
- Improved robustness of StructureMatcher using linear assignment.
- Improved primitive cell search (faster and more robust).
- Cleanup of deprecated methods, including
  pymatgen.alchemy.materials.TransformedMaterial.undo/redo_last_transformation,
  pymatgen.core.site.Site.distance_and_image_old, Poscar.struct,
  StructureFitter and tests.
- Miscellaneous bug fixes and speedups.

## v2.4.3

- Bug fix for StructureMatcher.
- Miscellaneous speedups.

## v2.4.0

- New StructureMatcher that effectively replaces StructureFitter. Orders of
  magnitude faster and more robust. StructureFitter is now deprecated.
- Vastly improved PrimitiveCellTransformation.
- A lot of core methods have been rewritten to take advantage of vectorization
  in numpy, resulting in orders of magnitude improvement in speed.
- Miscellaneous bug fixes and speedups.

## v2.3.2

- More utilities for working with Periodic Boundary Conditions.
- Improved MPRester that supports more data and a new method of specifying
  the API key for heavy users via a MAPI_KEY environment variable. Please
  refer to the :doc:`usage pages </usage>` for more information.
- Vastly improved POTCAR setup script in scripts directly that is now
  installed as part of a default pymatgen installation.
- Miscellaneous bug fixes and speedups.

## v2.3.1

- Significant improvements to the high-level interface to the Materials API.
  New interface provides more options to make it easier to get structures and
  entries, better warnings and error handling. It uses the *requests*
  library for a cleaner API.
- Bug fix for VolumetricData parsing and methods such as CHGCAR and LOCPOT.
  Previously, the parsing was done incorrectly because VASP actually provides
  data by running through the x-axis first, followed by y, then z.
- Bug fix for reverse_readline so that it works for gzipped and bzipped
  strucutures (courtesy of Anubhav Jain).
- Fix "lossy" composition to_dict method. Now composition.to_dict properly
  returns a correct species string as a key for compositions using species,
  instead of just the element symbols.
- Miscellaneous bug fixes.

## v2.3.0

- Remove usage of scipy and external qhull callers. Now uses pyhull package.
  Please note that this change implies that the pyhull package is now a
  required dependency. If you install pymatgen through the usual
  easy_install or pip install methods, this should be taken care of
  automatically for you. Otherwise, please look for the pyhull package on
  PyPI to download and install it.
- Miscellaneous bug fixes.

## v2.2.6

- Brand new _beta_ bond valence analyzer based on a Maximum A Posteriori
  algo using data-mined ICSD data.
- Speed up and improvements to core classes.
- Improved structure fitter (credits to Geoffroy Hautier).
- Brand new entry_tools module (pymatgen.entries.entry_tools).
- Vastly improved Outcar parser based on reverse parsing that speeds up
  reading of OUTCAR files by orders of magnitude.
- Miscellaneous bug fixes.

## v2.2.4

- Fixed bug in hexagonal cell KPOINTS file generation.
- New RelaxationAnalyzer to compare structures.
- New _beta_ bond valence analyzer.
- Miscellaneous bug fixes.

## v2.2.3

- New filter framework for filtering structures in pymatgen.alchemy.
- Updated feff io classes to support FEFF 9.6 and other code improvements.
- Miscellaneous bug fixes.

## v2.2.2

- Bug fix release for REST interface.
- Improvements to unittests.

## v2.2.1

- Improvements to feffio.
- Master matgenie.py script which replaces many analysis scripts.
- More memory efficient parsing of VolumetricData.
- Beta version of structure prediction classes.
- Changes to MPRester to work with v1 release of the Materials API.
- Miscellaneous bug fixes and speed improvements.

## v2.2.0

- Beta modules (pymatgen.io.feffio) for io for FEFF, courtesy of Alan Dozier.
- New smartio module that intelligently reads structure input files based on
  file extension.
- Spglib_adaptor module has been renamed to finder for brevity.
- Upgraded spglib to version 1.2.2. Improved handling of spglib install on
  Mac OS X and Solaris.
- Major cleanup of code for PEP8 compliance.
- Cssr module now supports reading of input files.
- Miscellaneous bug fixes and speed improvements.

## v2.1.2

- Brand new CompoundPD class that allows the plotting of phase diagrams that
  do not have elements as their terminal points.
- Spglib is now completely integrated as part of the setup.py installation.
- Major (but completely backwards compatible) refactoring of sites and vaspio.
- Added a EnumerateStructureTransformation with optional dependency on the enum
  library by Gus Hart. This provides a robust way to enumerate derivative
  structures,
- Implemented LLL lattice reduction algorithm. Also added option to sanitize
  a Structure on copy.
- Bug fix for missing Compatibility file in release distribution.
- Vastly improved StructureFitter which performs cell reduction where necessary
  to speed up fitting.
- Miscellaneous bug fixes and speed improvements.

## v2.0.0

- Brand new module (pymatgen.matproj.rest) for interfacing with the
  MaterialsProject REST interface.
- Useful aliases for commonly used Objects, similar in style to numpy.
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

- New PDAnalyzer method to generate chemical potential maps.
- Enhanced POSCAR class to support parsing of velocities and more formatting
  options.
- Reorganization of Bandstructure module. Beta support for projected
  bandstructure and eigenvalues in vaspio and electronic_structure.
- Miscellaneous bug fixes and speed improvements.

## v1.9.0

- Completely new json encoder and decoder that support serialization of almost
  all pymatgen objects.
- Simplification to Borg API utilizing the new json API.
- Bandstructure classes now support spin-polarized runs.
- Beta classes for battery (insertion and conversion) analysis.

## v1.8.3

- spglib_adaptor now supports disordered structures.
- Update to support new spglib with angle_tolerance.
- Changes to Borg API to support both file and directory style paths.
- Speed up for COMPLETE_ORDERING algo for PartialRemoveSpecieTransformation.

## v1.8.1

- Revamped transmuter classes for better readability and long term support.
- Much improved speed for PartialRemoveSpecieTransformations.
- Misc bug fixes.

## v1.8.0

- Support for additional properties on Specie (Spin) and Site (magmom, charge).
- Molecule class to support molecules without periodicity.
- Beta io class for XYZ and GaussianInput.

## v1.7.2

- Bug fixes for vaspio_set and compatibility classes.

## v1.7.0

- Complete reorganization of modules for electronic structure.
- Beta of band structure classes.
- Misc improvements to vaspio classes.
- Bug fixes.

## v1.6.0

- Beta of pymatgen.borg package implemented for high-throughput data assimilation.
- Added ComputedEntry classes for handling calculated data.
- New method of specifying VASP pseudopotential location using a VASP_PSP_DIR
  environment variable.
- Bug fix for pymatgen.symmetry
- Ewald sum speed up by factor of 2 or more.
