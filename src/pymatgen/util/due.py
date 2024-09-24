"""Stub file for a guaranteed safe import of duecredit constructs: if duecredit
is not available.

Then use in your code as

    from .due import due, Doi, BibTeX, Text

See https://github.com/duecredit/duecredit/blob/master/README.md for examples.

Origin:     Originally a part of the duecredit
Copyright:  2015-2021 DueCredit developers
License:    BSD-2
"""

from __future__ import annotations

import logging

__version__ = "0.0.9"


class InactiveDueCreditCollector:
    """Just a stub at the Collector which would not do anything."""

    def _donothing(self, *args, **kwargs):
        """Perform no good and no bad."""

    def dcite(self, *args, **kwargs):
        """If I could cite I would."""

        def non_decorating_decorator(func):
            return func

        return non_decorating_decorator

    active = False
    activate = add = cite = dump = load = _donothing

    def __repr__(self):
        return f"{type(self).__name__}()"


def _donothing_func(*args, **kwargs):
    """Perform no good and no bad."""


try:
    from duecredit import BibTeX, Doi, Text, Url, due

    if "due" in locals() and not hasattr(due, "cite"):
        raise RuntimeError("Imported due lacks .cite. DueCredit is now disabled")
except Exception as exc:
    if not isinstance(exc, ImportError):
        logging.getLogger("duecredit").exception("Failed to import duecredit")
    # Initiate due stub
    due = InactiveDueCreditCollector()
    BibTeX = Doi = Url = Text = _donothing_func
