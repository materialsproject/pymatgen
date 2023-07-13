## Summary

Major changes:

- feature 1: ...
- fix 1: ...

## Todos

If this is work in progress, what else needs to be done?

- feature 2: ...
- fix 2:

## Checklist

- [ ] Google format doc strings added. Check with `ruff`.
- [ ] Type annotations included. Check with `mypy`.
- [ ] Tests added for new features/fixes.
- [ ] If applicable, new classes/functions/modules have [`duecredit`](https://github.com/duecredit/duecredit) `@due.dcite` decorators to reference relevant papers by DOI ([example](https://github.com/materialsproject/pymatgen/blob/91dbe6ee9ed01d781a9388bf147648e20c6d58e0/pymatgen/core/lattice.py#L1168-L1172))

Tip: Install `pre-commit` hooks to auto-check types and linting before every commit:

```sh
pip install -U pre-commit
pre-commit install
```
