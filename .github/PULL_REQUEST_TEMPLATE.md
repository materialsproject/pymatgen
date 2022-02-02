## Summary

Include a summary of major changes in bullet points:

* Feature 1
* Feature 2
* Fix 1
* Fix 2

## Additional dependencies introduced (if any)

* List all new dependencies needed and justify why. While adding dependencies that bring
significantly useful functionality is perfectly fine, adding ones that
add trivial functionality, e.g., to use one single easily implementable
function, is frowned upon. Provide a justification why that dependency is needed.
Especially frowned upon are circular dependencies, e.g., depending on derivative
modules of pymatgen such as custodian or Fireworks.

## TODO (if any)

If this is a work-in-progress, write something about what else needs
to be done

* Feature 1 supports A, but not B.

## Checklist

Work-in-progress pull requests are encouraged, but please put [WIP]
in the pull request title.

Before a pull request can be merged, the following items must be checked:

- [ ] Code is in the [standard Python style](https://www.python.org/dev/peps/pep-0008/). The easiest way to handle this
      is to run the following in the **correct sequence** on your local machine. Start with running
      [black](https://black.readthedocs.io/en/stable/index.html) on your new code. This will automatically reformat
      your code to PEP8 conventions and removes most issues. Then run
      [pycodestyle](https://pycodestyle.readthedocs.io/en/latest/), followed by
      [flake8](http://flake8.pycqa.org/en/latest/).
- [ ] Docstrings have been added in the [Google docstring format](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).
      Run [pydocstyle](http://www.pydocstyle.org/en/2.1.1/index.html) on your code.
- [ ] Type annotations are **highly** encouraged. Run [mypy](http://mypy-lang.org/) to type check your code.
- [ ] Tests have been added for any new functionality or bug fixes.
- [ ] All linting and tests pass.

Note that the CI system will run all the above checks. But it will be much more efficient if you already fix most
errors prior to submitting the PR. It is highly recommended that you use the pre-commit hook provided in the pymatgen
repository. Simply `cp pre-commit .git/hooks` and a check will be run prior to allowing commits.
