Change log
==========

v2019.8.23
----------
* pycodestyle now enforced, except on tests. Developers should install
  pycodestyle and the pre-commit hook (copy pre-commit to .git/hooks) 
  provided in the repo to check before commits. CI now checks for code style
  and PRs must pass pycodestyle.
* chemsys str input now allowed in get_entries_in_chemsys (@rkingsbury)
* ComputedEntry and subclasses now support a normalize().
* Speed improvements in fragmeter using igraph. (@samblau)
