Change log
==========

v2023.7.17
----------
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
