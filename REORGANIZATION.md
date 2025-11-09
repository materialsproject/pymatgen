# Pymatgen Reorganization Guide

This document describes the reorganization of pymatgen into a mono-repo structure with four packages:

1. **pymatgen-core**: Core functionality (structures, I/O, etc.) - no dependency on pymatgen-analysis or pymatgen-apps
2. **pymatgen-analysis**: Analysis tools (phase diagrams, structure matching, etc.) - depends on pymatgen-core
3. **pymatgen-apps**: Applications and CLI tools (battery analysis, borg workflows, command-line interfaces) - depends on pymatgen-core and pymatgen-analysis
4. **pymatgen** (meta-package): Installs pymatgen-core, pymatgen-analysis, and pymatgen-apps

## Directory Structure

```
pymatgen/
├── pymatgen-core/
│   ├── src/
│   │   └── pymatgen/
│   │       ├── core/
│   │       ├── io/
│   │       ├── entries/
│   │       └── ... (everything except analysis, apps, cli)
│   ├── pyproject.toml
│   └── setup.py
├── pymatgen-analysis/
│   ├── src/
│   │   └── pymatgen/
│   │       └── analysis/
│   ├── pyproject.toml
│   └── setup.py
├── pymatgen-apps/
│   ├── src/
│   │   └── pymatgen/
│   │       ├── apps/
│   │       └── cli/
│   ├── pyproject.toml
│   └── setup.py
├── pymatgen-meta/
│   ├── src/
│   │   └── pymatgen/
│   │       └── __init__.py
│   ├── pyproject.toml
│   └── setup.py
└── src/  (original source, can be removed after reorganization)
```

## Steps to Complete Reorganization

### 1. Copy Files

Run the reorganization script to copy files:

```bash
python3 reorganize.py
```

This will:
- Copy `src/pymatgen/analysis` to `pymatgen-analysis/src/pymatgen/analysis` (excluding `structure_matcher.py`)
- Move `src/pymatgen/analysis/structure_matcher.py` to `pymatgen-core/src/pymatgen/core/structure_matcher.py`
- Copy `src/pymatgen/apps` and `src/pymatgen/cli` to `pymatgen-apps/src/pymatgen/`
- Copy everything else to `pymatgen-core/src/pymatgen/`

**Note**: `structure_matcher` has been moved from `pymatgen.analysis.structure_matcher` to `pymatgen.core.structure_matcher`. A backward-compatibility stub is created in the analysis package to maintain existing imports.

### 2. Create Structure Matcher Stub

Create a backward-compatibility stub for `structure_matcher`:

```bash
python3 create_structure_matcher_stub.py
```

This creates a stub file in `pymatgen-analysis/src/pymatgen/analysis/structure_matcher.py` that imports from `pymatgen.core.structure_matcher` to maintain backward compatibility. This allows existing code using `from pymatgen.analysis.structure_matcher import StructureMatcher` to continue working.

### 3. Refactor StructureMatcher Imports

Update imports in pymatgen-core files to use the new location:

```bash
python3 refactor_structure_matcher_imports.py
```

This updates all imports from `pymatgen.analysis.structure_matcher` to `pymatgen.core.structure_matcher` in files that will be part of pymatgen-core.

**Note**: This step has already been completed in the source files. The script is provided for reference.

### 4. Make Imports Optional

Some files in pymatgen-core import from pymatgen.analysis (other than structure_matcher). These need to be made optional:

```bash
python3 make_imports_optional.py
```

This script updates module-level imports to use try/except blocks. You may also need to add runtime checks for `None` values where these imports are used.

### 5. Files That Need Manual Updates

The following files import from pymatgen.analysis and may need additional runtime checks:

- `pymatgen-core/src/pymatgen/core/interface.py` - Uses AdsorbateSiteFinder
- `pymatgen-core/src/pymatgen/transformations/standard_transformations.py` - Uses multiple analysis modules (StructureMatcher already refactored to core)
- `pymatgen-core/src/pymatgen/entries/compatibility.py` - Uses structure_analyzer functions
- `pymatgen-core/src/pymatgen/entries/correction_calculator.py` - Uses ComputedReaction
- `pymatgen-core/src/pymatgen/entries/mixing_scheme.py` - Uses PhaseDiagram (StructureMatcher already refactored to core)
- `pymatgen-core/src/pymatgen/io/vasp/sets.py` - (StructureMatcher already refactored to core)
- `pymatgen-core/src/pymatgen/command_line/vampire_caller.py` - Uses HeisenbergMapper
- `pymatgen-core/src/pymatgen/command_line/gulp_caller.py` - Uses BVAnalyzer
- `pymatgen-core/src/pymatgen/cli/pmg_structure.py` - (Note: cli is in pymatgen-apps, not core)
- `pymatgen-core/src/pymatgen/cli/pmg_plot.py` - Uses XRDCalculator
- `pymatgen-core/src/pymatgen/transformations/site_transformations.py` - Uses EwaldMinimizer, EwaldSummation
- `pymatgen-core/src/pymatgen/phonon/thermal_displacements.py` - (StructureMatcher already refactored to core)
- `pymatgen-core/src/pymatgen/io/nwchem.py` - Uses ExcitationSpectrum
- `pymatgen-core/src/pymatgen/entries/entry_tools.py` - Uses PDEntry (StructureMatcher already refactored to core)
- `pymatgen-core/src/pymatgen/apps/battery/conversion_battery.py` - Uses PhaseDiagram, BalancedReaction

For each of these, you should:
1. Make the import optional (try/except)
2. Add runtime checks to raise helpful error messages if the analysis package is not installed

Example:

```python
try:
    from pymatgen.analysis.structure_matcher import StructureMatcher
except ImportError:
    StructureMatcher = None


def some_function():
    if StructureMatcher is None:
        raise ImportError(
            "This functionality requires pymatgen-analysis. "
            "Install it with: pip install pymatgen-analysis"
        )
    # ... use StructureMatcher
```

### 4. Testing

After reorganization, test each package:

```bash
# Test pymatgen-core
cd pymatgen-core
pip install -e .
python -c "import pymatgen.core; print('Core works!')"

# Test pymatgen-analysis
cd ../pymatgen-analysis
pip install -e .
python -c "import pymatgen.analysis; print('Analysis works!')"

# Test meta-package
cd ../pymatgen-meta
pip install -e .
python -c "import pymatgen; print('Meta-package works!')"
```

### 5. Update Tests

Tests should be reorganized similarly:
- `tests/core/`, `tests/io/`, etc. → `pymatgen-core/tests/`
- `tests/analysis/` → `pymatgen-analysis/tests/`

### 6. Update Documentation

Update documentation to reflect the new package structure and installation instructions.

## Notes

- The meta-package (`pymatgen`) maintains backward compatibility by re-exporting from all sub-packages
- `pymatgen-core` has NO dependency on `pymatgen-analysis` or `pymatgen-apps` - this is critical
- `structure_matcher` has been moved from `pymatgen.analysis.structure_matcher` to `pymatgen.core.structure_matcher`
  - A backward-compatibility stub is provided in `pymatgen-analysis` to maintain existing imports
  - New code should use `from pymatgen.core.structure_matcher import StructureMatcher`
- Users can install just `pymatgen-core` if they don't need analysis tools or CLI
- The version numbers should be kept in sync across all four packages
