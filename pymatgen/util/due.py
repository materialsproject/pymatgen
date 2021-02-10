"""

Stub file for a guaranteed safe import of duecredit constructs:  if duecredit
is not available.

See  https://github.com/duecredit/duecredit/blob/master/README.md

Origin:     Originally a part of the duecredit
Copyright:  2015-2019  DueCredit developers
License:    BSD-2
"""

__version__ = '0.0.8'


class InactiveDueCreditCollector(object):
    """Just a stub at the Collector which would not do anything"""
    def _donothing(self, *args, **kwargs):
        """Perform no good and no bad"""
        pass

    def dcite(self, *args, **kwargs):
        """If I could cite I would"""
        def nondecorating_decorator(func):
            return func
        return nondecorating_decorator

    active = False
    activate = add = cite = dump = load = _donothing

    def __repr__(self):
        return self.__class__.__name__ + '()'


def _donothing_func(*args, **kwargs):
    """Perform no good and no bad"""
    pass


try:
    from duecredit import due, BibTeX, Doi, Url, Text
    if 'due' in locals() and not hasattr(due, 'cite'):
        raise RuntimeError(
            "Imported due lacks .cite. DueCredit is now disabled")
except Exception as e:
    if not isinstance(e, ImportError):
        import logging
        logging.getLogger("duecredit").error(
            "Failed to import duecredit due to %s" % str(e))
    # Initiate due stub
    due = InactiveDueCreditCollector()
    BibTeX = Doi = Url = Text = _donothing_func

# Emacs mode definitions
# Local Variables:
# mode: python
# py-indent-offset: 4
# tab-width: 4
# indent-tabs-mode: nil
# End:
