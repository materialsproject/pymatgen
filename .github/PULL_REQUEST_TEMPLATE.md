## Summary

Include a summary of major changes in bullet points:

- Feature 1
- Feature 2
- Fix 1
- Fix 2

## Todo (if any)

If this is a work-in-progress, write something about what else needs to be done

- Feature 1 supports A, but not B.

## Checklist

Work-in-progress pull requests are encouraged, but please put \[WIP\] in the pull request title.

Before a pull request can be merged, the following items must be checked:

- [ ] Doc strings have been added in the [Google docstring format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html). Run [pydocstyle](http://www.pydocstyle.org/en/2.1.1/index.html) on your code.
- [ ] Type annotations are *highly* encouraged. Run [`mypy path/to/file.py`](https://github.com/python/mypy) to type check your code.
- [ ] Tests have been added for any new functionality or bug fixes.
- [ ] All linting and tests pass.

Note that the CI system will run all the above checks. But it will be much more efficient if you already fix most errors prior to submitting the PR. We highly recommended installing `pre-commit` hooks. Simply Run

```sh
pip install -U pre-commit
pre-commit install
```

in the repo's root directory. Afterwards linters will run before every commit and abort if any issues pop up.
