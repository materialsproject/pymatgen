# Temporary PR Notes: CLI Startup Perf

## Scope

- Optimize `pmg` startup by lazily importing subcommand handlers.
- Reduce import-time overhead in `pmg analyze`, `pmg plot`, `pmg structure`, and `pmg potcar`.
- Keep behavior unchanged while shifting heavy scientific imports onto the command paths that actually use them.

## Main Code Changes

- `src/pymatgen/cli/pmg.py`
  - Add lazy handler dispatch via `_lazy_func(...)`.
  - Defer `Structure` and `Incar` imports for `view` and `diff`.
  - Remove eager `SETTINGS` / `Potcar` dependency from parser construction.

- `src/pymatgen/cli/pmg_plot.py`
  - Move `XRDCalculator`, `Structure`, `DosPlotter`, `Vasprun`, `Chgcar`, `SpacegroupAnalyzer`, and `pretty_plot` into per-command functions.

- `src/pymatgen/cli/pmg_analyze.py`
  - Move `BorgQueen`, VASP drones, and `Outcar` imports into runtime paths.

- `src/pymatgen/cli/pmg_structure.py`
  - Move `Structure`, `SpacegroupAnalyzer`, and `StructureMatcher` imports into runtime paths.

- `src/pymatgen/cli/pmg_potcar.py`
  - Move `SETTINGS` and `Potcar` imports into runtime paths.
  - Resolve default functional at execution time instead of parser construction time.

## Benchmarks

### `pmg` cold import

- Before optimization: about `0.97 s`
- After optimization: about `0.055 s` mean over 8 cold runs

### `pymatgen.cli.pmg_plot` import

- Before optimization: about `2.13 s`
- After optimization: about `0.56 s`

### Cold-run subcommand samples

- `pmg plot --chgint ...`
  - before: about `2.214 s`
  - after: about `1.824 s`

- `pmg plot --dos ...`
  - after lazy-import refactor: about `1.280 s`

- `pmg plot --xrd ...`
  - after lazy-import refactor: about `1.562 s`
  - observed some run-to-run variance

- `pmg structure --symmetry ...`
  - before: about `0.408 s`
  - after: about `0.314 s`

## Validation

- `uv run pytest tests/cli/test_pmg.py tests/cli/test_pmg_plot.py tests/cli/test_pmg_structure.py tests/cli/test_pmg_analyze.py`
- `uv run pytest tests/cli/test_pmg_plot.py`
- `uv run pytest tests/cli/test_pmg_analyze.py tests/cli/test_pmg_structure.py`

All relevant CLI tests passed during local verification.

## PR Framing Ideas

- Improve `pmg` startup performance by lazily importing CLI handlers and heavy scientific dependencies.
- Reduce import-time overhead for `pmg analyze`, `pmg plot`, `pmg structure`, and `pmg potcar` without changing command behavior.
