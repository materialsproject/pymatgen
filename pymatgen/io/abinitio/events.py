"""
This module defines the events signaled by abinit during the execution. It also
provides a parser to extract these events form the main output file and the log file.
"""
from __future__ import division, print_function

import os.path
import collections
import yaml

from pymatgen.util.string_utils import WildCard
from pymatgen.io.abinitio.abiinspect import YamlTokenizer, YamlDoc

__all__ = [
    "EventParser",
]

class AbinitEvent(yaml.YAMLObject):
    """
    Example (YAML syntax)::

        Normal warning without any handler:

        --- !Warning
        message: | 
            This is a normal warning that won't 
            trigger any handler in the python code!
        src_file: routine_name
        src_line:  112
        ...

        Critical warning that will trigger some action in the python code.

        --- !ScfConvergeWarning
        message: |
            The human-readable message goes here!
        src_file: foo.F90
        src_line: 112
        tolname: tolwfr
        actual_tol: 1.0e-8
        required_tol: 1.0e-10
        nstep: 50
        ...

    The algorithm to extract the YAML sections is very simple.

    1) Find --- at the beginning of line.
    2) If we have a tag that ends with "Warning", "Error", "Bug", "Comment
       we know we have encountered a new ABINIT event 
    3) The key defines the event class to instantiate.

    Note that now --- and ... become reserved words (whey they are placed at
    the begining of a line) since they are used to mark the beginning and 
    the end of YAML documents.
    """
    def __init__(self, message, scr_file, src_line):
        self.message = message
        self.src_file = src_file
        self.src_line = src_line

    def __str__(self):
        header = "%s:%s" % (self.src_file, self.src_line)
        return "\n".join((header, self.message))

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def baseclass(self):
        """The baseclass of self."""
        for cls in _BASE_CLASSES:
            if isinstance(self, cls):
                return cls

        err_msg = "Cannot determine the base class of %s" % (self.__class__.__name__)
        raise ValueError(err_msg)

    def iscritical(self):
        """
        True if event is critical namely that if this event should be analyzed in
        more detail to understand what action should be performed
        """
        return False

    def action(self):
        """
        Returns a dictionary whose values that can be used to decide
        which actions should be performed e.g the SCF data at the last
        iteration can be used to decide whether the calculations should
        be restarted or not.
        """
        return {}


class AbinitComment(AbinitEvent):
    """Base class for Comment events"""
    yaml_tag = u'!COMMENT'


class AbinitError(AbinitEvent):
    """Base class for Error events"""
    yaml_tag = u'!ERROR'

    @property
    def iscritical(self):
        return True


class AbinitBug(AbinitEvent):
    """Base class for Bug events"""
    yaml_tag = u'!BUG'

    @property
    def iscritical(self):
        return True


class AbinitWarning(AbinitEvent):
    yaml_tag = u'!WARNING'
    """
    Base class for Warning events (the most important class).
    Developers should subclass this class to define the different exceptions
    raised by the code and the possible actions that can be performed.
    """
    # FIXME: for the moment we tag a warning as critical, then, once we migrate to the
    # JSON-like format, only CriticalWarnings will trigger some kind of action.
    @property
    def iscritical(self):
        return True

#class ScfConvergenceWarning(AbinitWarning):
#yaml_tag = u'ScfConvergenceWarning'

#class NscfConvergenceWarning(AbinitWarning):
#yaml_tag = u'NscfConvergenceWarning'

#class RelaxConvergenceWarning(AbinitWarning):
#yaml_tag = u'RelaxConvergenceWarning'

#class PhononConvergenceWarning(AbinitWarning):
#yaml_tag = u'PhononConvergenceWarning'

#class QPSConvergenceWarning(AbinitWarning):
#yaml_tag = u'QPSConvergenceWarning'

#class HaydockConvergenceWarning(AbinitWarning):
#yaml_tag = u'HaydockConvergenceWarning'

# Register the concrete base classes.
_BASE_CLASSES = [
    AbinitComment,
    AbinitError,
    AbinitBug,
    AbinitWarning,
]


class EventReport(collections.Iterable):
    """Iterable storing the events raised by an ABINIT calculation."""

    def __init__(self, filename, events=None):
        """
        Args:
            filename:
                Name of the file
            events:
                List of Event objects
        """
        self.filename = os.path.abspath(filename)
        self._events = []
        self._events_by_baseclass = collections.defaultdict(list)

        if events is not None:
            for ev in events:
                self.append(ev)

    def __len__(self):
        return len(self._events)

    def __iter__(self):
        return self._events.__iter__()

    def __str__(self):
        lines = [self.filename + ":",]
        for event in self:
            lines.append(str(event))

        return "\n".join(l for l in lines)

    def append(self, event):
        """Add an event to the list."""
        self._events.append(event)
        self._events_by_baseclass[event.baseclass].append(event)

    def set_run_completed(self, bool_value):
        """Set the value of _run_completed."""
        self._run_completed = bool_value

    @property
    def run_completed(self):
        """
        Returns True if the calculation terminated.
        """
        try:
            return self._run_completed
        except AttributeError:
            return False

    @property
    def critical_events(self):
        """List of critical events."""
        return [e for e in self if e.iscritical]

    @property
    def comments(self):
        """List of comments found."""
        return self.select(AbinitComment)

    @property
    def errors(self):
        """List of errors found."""
        return self.select(AbinitError)

    @property
    def bugs(self):
        """List of bugs found."""
        return self.select(AbinitBug)

    @property
    def warnings(self):
        """List of warnings found."""
        return self.select(AbinitWarning)

    @property
    def num_warnings(self):
        """Number of warnings reported."""
        return len(self.warnings)

    @property
    def num_errors(self):
        """Number of errors reported."""
        return len(self.errors)

    @property
    def num_comments(self):
        """Number of comments reported."""
        return len(self.comments)

    def select(self, base_class, only_critical=False):
        """
        Return list of events that inherits from class base_class

        Args:
            only_critical:
                if True, only critical events are returned.
        """
        if only_critical:
            return [e for e in self._events_by_baseclass[base_class] if e.iscritical]
        else:
            return self._events_by_baseclass[base_class][:]



class EventParserError(Exception):
    """Base class for the exceptions raised by EventParser."""


class EventParser(object):
    """
    Parses the output or the log file produced by abinit and extract the list of events.
    """
    Error = EventParserError

    @staticmethod
    def parse(filename):
        """
        This is the new parser, it will be used when we implement
        the new format in abinit.
        """
        run_completed = False
        filename = os.path.abspath(filename)
        report = EventReport(filename)

        #w = WildCard("*Error|*Warning|*Comment")
        w = WildCard("*ERROR|*WARNING|*COMMENT")

        with YamlTokenizer(filename) as tokens:

            for doc in tokens:
                #print(80*"*")
                #print("doc.tag", doc.tag)
                #print("doc",doc)
                #print(80*"*")
                if w.match(doc.tag):
                    #print("got doc.tag", doc.tag)
                    event = yaml.load(doc.text)
                    event.lineno = doc.lineno
                    report.append(event)

                #if doc.tag == "!FinalSummary":
                #    run_completed = True

        # TODO: Add YAML doc with FinalSummary.
        # Check whether the calculation completed.
        MAGIC = "Calculation completed."
        with open(filename) as fh:
            for line in fh:
                if MAGIC in line:
                    run_completed = True
                    break

        report.set_run_completed(run_completed)
        return report

    def report_exception(self, filename, exc):
        """
        This method is used when self.parser raises an Exception so that
        we can report a customized `EventReport` object with info the exception.
        """
        return EventReport(filename, events=[Error(str(exc))])

