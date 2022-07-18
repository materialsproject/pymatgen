Change log
==========

v2022.7.8
---------
Welcome to new contributors @naveensrinivasan, @xivh, @dgaines2, @yang-ruoxi, @cajfisher and @mjwen!

* New: Partial periodic boundary conditions, PR #2429 by @gpetretto
* New: Element.from_name(), PR #2567 by @rkingsbury
* New: Materials Project input set for absorption calculations, PR #2320 by @yang-ruoxi
* Enhancement: compressed LAMMPS and XYZ files in pymatgen.io.lammps, PR #2538 by @ab5424
* Enhancement: remove vertical lines from VoltageProfilePlotter.get_plotly_figure(), PR #2552 by @acrutt
* Enhancement: chemical potential plot background color changed, PR #2559 @jmmshn
* Enhancement: ability to change voronoi_distance_cutoff in ChemEnv, PR #2568 by @JaGeo
* Enhancement: Ion.oxi_state_guesses will use correct charge by default, PR #2566 by @rkingsbury
* Enhancement: Remove not converged warning for VASP AIMD runs, PR #2571 by @mjwen
* Fix: generation of continuous line-mode band structures, PR #2533 by @munrojm
* Fix: duplicate site properties for magnetic moments hwen using `AseAtomsAdaptor`, PR #2545 by @arosen93
* Fix: bug in Grüneisen parameter calculation, PR #2543 by by @ab5424
* Fix: allow a comment on final line of KPOINTS file, PR #2549 by @xivh
* Fix: for `Composition.replace` with complex mappings, PR #2555 by @jacksund
* Fix: Implement equality method and fix __iter__ for InputSet, PR #2575 by @rkingsbury
* Fix: use negative charge convention for electron in "update_charge_from_potcar", PR #2577 by @jmmshn
* Fix: ensure charge is applied to initial and final structures parsed from vasprun.xml, PR #2579 by @jmmshn
* Chore: Set permissions for GitHub actions, PR #2547 by @naveensrinivasan
* Chore: Included GitHub actions in the Dependabot config, PR #2548 by by @naveensrinivasan
* Documentation: fix typos in pymatgen.symmetry.analyzer docstrings, PR #2561 by @dgaines2
* Documentation: clarification about usage of InputFile, PR #2570 by by @orionarcher
* Documentation: Improve messages and warnings, PR #2572 and PR #2573 by @cajfisher
* Documentation: fix typo, PR #2580 by @janosh

Notice: functionality from pymatgen.analysis.defects will be incorporated into a separate add-on package in the future,
see deprecation notice.
