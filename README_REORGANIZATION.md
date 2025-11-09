# Pymatgen Reorganization - Quick Start

The codebase has been reorganized into a mono-repo structure. Here's what you need to do:

## Quick Start

1. **Copy the files** (run this first):
   ```bash
   python3 reorganize.py
   ```

2. **Create structure_matcher stub** (run this after copying):
   ```bash
   python3 create_structure_matcher_stub.py
   ```

3. **Refactor StructureMatcher imports** (already done in source):
   ```bash
   python3 refactor_structure_matcher_imports.py
   ```
   Note: This has already been completed in the source files.

4. **Make imports optional** (run this after copying):
   ```bash
   python3 make_imports_optional.py
   ```

5. **Review and fix imports manually** - See `REORGANIZATION.md` for a list of files that need manual updates to handle optional imports gracefully.

## What's Been Created

### Package Structure
- ✅ `pymatgen-core/` - Core package (no dependency on analysis or apps)
- ✅ `pymatgen-analysis/` - Analysis package (depends on core)
- ✅ `pymatgen-apps/` - Applications and CLI package (depends on core and analysis)
- ✅ `pymatgen-meta/` - Meta-package that installs all three

### Configuration Files
- ✅ `pymatgen-core/pyproject.toml` - Core package configuration
- ✅ `pymatgen-analysis/pyproject.toml` - Analysis package configuration
- ✅ `pymatgen-apps/pyproject.toml` - Apps package configuration
- ✅ `pymatgen-meta/pyproject.toml` - Meta-package configuration
- ✅ `pymatgen-core/setup.py` - Cython extensions for core
- ✅ `pymatgen-analysis/setup.py` - Setup for analysis
- ✅ `pymatgen-apps/setup.py` - Setup for apps
- ✅ `pymatgen-meta/setup.py` - Setup for meta-package
- ✅ `pymatgen-meta/src/pymatgen/__init__.py` - Meta-package imports

### Scripts
- ✅ `reorganize.py` - Copies files to new structure
- ✅ `make_imports_optional.py` - Makes analysis imports optional in core

### Documentation
- ✅ `REORGANIZATION.md` - Detailed reorganization guide

## Next Steps

After running the scripts:

1. **Test each package** to ensure they work independently
2. **Add runtime checks** for optional imports (see REORGANIZATION.md)
3. **Update tests** to match the new structure
4. **Update CI/CD** to build and test all three packages

## Important Notes

- `pymatgen-core` MUST NOT depend on `pymatgen-analysis` or `pymatgen-apps`
- `pymatgen-apps` depends on both `pymatgen-core` and `pymatgen-analysis`
- **`structure_matcher` has been moved**: from `pymatgen.analysis.structure_matcher` to `pymatgen.core.structure_matcher`
  - A backward-compatibility stub allows existing imports to continue working
  - New code should use `from pymatgen.core.structure_matcher import StructureMatcher`
- The meta-package maintains backward compatibility
- All four packages should share the same version number
- Users can install just `pymatgen-core` if they don't need analysis tools or CLI
