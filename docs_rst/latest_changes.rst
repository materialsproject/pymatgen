Change log
==========

v2023.2.22
----------
* PR #2848 from @ml-evs ml-evs/update_optimade_aliases
    Currently `OptimadeRester` defaults to an outdated list of OPTIMADE database URLs (several of which fail) and the design of the class is such that refreshing these aliases can only be done post-init which means they will not be used if the user provides their own filtered list of aliases, without doing some extra work.
    This PR refreshes the vendored list of aliases (which should be much more stable now since their initial addition 2 years ago), and also adds the option to refresh the aliases on initialization of the class.
    This currently affects the pymatgen OPTIMADE tutorials at https://github.com/Materials-Consortia/optimade-tutorial-exercises.
