Change log
==========

v2023.5.31
----------
- Breaking: Default ``user_potcar_settings`` to ``{"W": "W_sv"}`` in all input sets if ``user_potcar_functional == "PBE_54"`` `#3022 <https://github.com/materialsproject/pymatgen/pull/3022>`_
- Unignore ``ruff`` ``PD011`` `#3020 <https://github.com/materialsproject/pymatgen/pull/3020>`_
- Tweak variable names `#3019 <https://github.com/materialsproject/pymatgen/pull/3019>`_
- MaterialsProjectCompatibility issue silencable deprecation warning `#3017 <https://github.com/materialsproject/pymatgen/pull/3017>`_
- Optimize cython find_points_in _spheres `#3015 <https://github.com/materialsproject/pymatgen/pull/3015>`_
- Cp2k 2.0 `#2672 <https://github.com/materialsproject/pymatgen/pull/2672>`_
- Added methods to compute and compare DOS fingerprints `#2772 <https://github.com/materialsproject/pymatgen/pull/2772>`_
- Breaking: Overhaul ``class PymatgenTest`` `#3014 <https://github.com/materialsproject/pymatgen/pull/3014>`_
- Fix ``ValueError`` when ``structure.selective_dynamics`` has type ``np.array`` `#3012 <https://github.com/materialsproject/pymatgen/pull/3012>`_
- Clean up `#3010 <https://github.com/materialsproject/pymatgen/pull/3010>`_
- Update ``.pytest-split-durations`` `#3005 <https://github.com/materialsproject/pymatgen/pull/3005>`_
- Lookup ``MPRester`` API key in settings if ``None`` provided as arg `#3004 <https://github.com/materialsproject/pymatgen/pull/3004>`_
- Support writing structures to compressed JSON (.json.gz .json.bz2 .json.xz .json.lzma) `#3003 <https://github.com/materialsproject/pymatgen/pull/3003>`_
- Add LightStructureEnvironments.from_structure_environments() fallback value if ``ce_and_neighbors`` is None `#3002 <https://github.com/materialsproject/pymatgen/pull/3002>`_
- ``Species`` parse oxi state from symbol str `#2998 <https://github.com/materialsproject/pymatgen/pull/2998>`_
- Re-export ``SiteCollection`` + ``DummySpecies`` from ``pymatgen.core`` `#2995 <https://github.com/materialsproject/pymatgen/pull/2995>`_
- Orbital-resolved icohplist `#2993 <https://github.com/materialsproject/pymatgen/pull/2993>`_
- Hide all type-hint-only imports behind ``if TYPE_CHECKING`` `#2992 <https://github.com/materialsproject/pymatgen/pull/2992>`_
- Add type hints for ``pymatgen.io.ase`` module `#2991 <https://github.com/materialsproject/pymatgen/pull/2991>`_
- Enable ruff doc rules in CI `#2990 <https://github.com/materialsproject/pymatgen/pull/2990>`_
- Suspected Typo Fix in ``pymatgen.io.vasp.optics`` `#2989 <https://github.com/materialsproject/pymatgen/pull/2989>`_
- Doc strings `#2987 <https://github.com/materialsproject/pymatgen/pull/2987>`_
- Fix average error `#2986 <https://github.com/materialsproject/pymatgen/pull/2986>`_
- Drop deprecated SubstrateAnalyzer + ZSLGenerator reexports `#2981 <https://github.com/materialsproject/pymatgen/pull/2981>`_
- Breaking: Default ``user_potcar_settings`` to ``{"W": "W_sv"}`` in all input sets if ``user_potcar_functional == "PBE_54"`` (#3022) `#3022 <https://github.com/materialsproject/pymatgen/pull/3022>`_
- fix unwanted x margins in get_elt_projected_plots_color (closes #562) `#562 <https://github.com/materialsproject/pymatgen/issues/562>`_
- Add LightStructureEnvironments.from_structure_environments() fallback value if ``ce_and_neighbors`` is None (#3002) `#2756 <https://github.com/materialsproject/pymatgen/issues/2756>`_
- add doc str explaining need for class ElementBase (closes #2999) `#2999 <https://github.com/materialsproject/pymatgen/issues/2999>`_
- Update docs. `3e3c31c <https://github.com/materialsproject/pymatgen/commit/3e3c31c8d342c84f2c6bbb961c321e458b9accb9>`_
- ruff set isort.split-on-trailing-comma = false `c0ec534 <https://github.com/materialsproject/pymatgen/commit/c0ec53452c3dc87c6cca5edc1c6b2b6218f15569>`_
