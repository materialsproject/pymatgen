"""
This module defines the events signaled by abinit during the execution. It also
provides a parser to extract these events form the main output file and the log file.
"""
from __future__ import division, print_function

import os.path
import collections

__all__ = [
    "EventParser",
]

class AbinitEvent(object):
    """
    Example (YAML syntax)::

        Normal warning without any handler:

        --- !Warning
        message: "This is a normal warning that won't trigger any handler in the python code!"
        src_file: routine_name
        src_line:  112
        ...

        Critical warning that will trigger some action in the python code.

        --- !ScfConvergeWarning
        message: "The human-readable message goes here!"
        src_file: foo.F90
        src_line: 112
        tolname: tolwfr
        actual_tol: 1.0e-8
        required_tol: 1.0e-10
        nstep: 50
        ...

    The algorithm to extract the YAML sections is very simple.

    1) Find --- at the beginning of line.
    2) If the next lines ends with "Warning:", "Error:", "Bug:", "Comment
       we know we have encountered a new ABINIT event 
    3) The key defines the event class to instantiate.

    Note that now --- and ... become reserved words (whey they are placed at
    the begining of a line) since they are used to mark the beginning and 
    the end of YAML documents.
    """
    def __init__(self, **kwargs):
        self._kwargs = kwargs.copy()
        self.message = kwargs.pop("message")
        self.lineno  = kwargs.pop("lineno")
        self.data = kwargs

    def __str__(self):
        return "%s:\n%s" % (self.lineno, self.message)

    @property
    def name(self):
        return self.__class__.__name__

    @property
    def baseclass(self):
        for cls in _BASE_CLASSES:
            if isinstance(self, cls):
                return cls
        raise ValueError("Cannot determine the base class of %s" %
                         self.__class__.__name__)

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


class AbinitError(AbinitEvent):
    """Base class for Error events"""
    @property
    def iscritical(self):
        return True


class AbinitBug(AbinitEvent):
    """Base class for Bug events"""
    @property
    def iscritical(self):
        return True


class AbinitWarning(AbinitEvent):
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
#class NscfConvergenceWarning(AbinitWarning):
#class RelaxConvergenceWarning(AbinitWarning):
#class PhononConvergenceWarning(AbinitWarning):
#class QPSConvergenceWarning(AbinitWarning):
#class HaydockConvergenceWarning(AbinitWarning):

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
        lines = [self.filename+":",]
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
    def parse(filename, nafter=5):
        """
        Read and parse the main output file or the file produced for abinit

        Args:
            filename:
                path to the file.
            nafter:
                Save nafter lines of trailing context after matching lines.

        Returns:
            `EventReport` instance.
        """
        filename = os.path.abspath(filename)

        # TODO
        # we have to standardize the abinit WARNING, COMMENT and ERROR  
        # so that we can parse them easily without having to use nafter.

        # Note the space after the name.
        exc_cases = ["ERROR ", "BUG ", "WARNING ", "COMMENT "]

        errors, bugs, warnings, comments = [], [], [], []

        handlers = {
            "ERROR "  : errors.append,
            "BUG "    : bugs.append,
            "WARNING ": warnings.append,
            "COMMENT ": comments.append,
        }

        def exc_case(line):
            for e in exc_cases:
                if e in line: return e
            else:
                return None

        MAGIC = "Calculation completed."
        run_completed = False

        with open(filename, "r") as fh:
            lines = fh.readlines()
            nlines = len(lines)
            for (lineno, line) in enumerate(lines):
                if MAGIC in line:
                    run_completed = True
                handle = handlers.get(exc_case(line))
                if handle is None: continue
                context = lines[lineno: min(lineno+nafter, nlines)]
                handle((lineno, "".join(c for c in context)))

        events = EventReport(filename)

        for lineno, s in errors:
            events.append(AbinitError(lineno=lineno, message=s))

        for lineno, s in bugs:
            events.append(AbinitBug(lineno=lineno, message=s))

        for lineno, s in warnings:
            events.append(AbinitWarning(lineno=lineno, message=s))

        for lineno, s in comments:
            events.append(AbinitComment(lineno=lineno, message=s))

        events.set_run_completed(run_completed)
        return events

    @staticmethod
    def new_event_parser(filename):
        """
        This is the new parser, it will be used when we implement
        the new format in abinit.
        """
        filename = os.path.abspath(filename)

        START_TAG = "<Event>"
        STOP_TAG = "<\Event>"

        exc_cases = ["ERROR ", "BUG ", "WARNING ", "COMMENT "]

        # Find YAML documents with events
        def is_event_start(line):
            if not line.startswith("..."):
                return False

            i = line.index("!")
            if i == -1: return False

            class_name = line[i:]
            if filter(class_name.endswith, exc_cases):
                return True

            return False

        doc_list = []
        with open(filename, "r") as fh:
            in_event = False
            for line in fh:

                if is_event_start(line):
                    in_event = False
                    doc_list.append(doc)

                if is_event_start(line):
                    in_event = True
                    doc = []

                if in_event:
                    doc.append(line)

        if doc:
            doc_list.append(doc)

        events = EventReport(filename)

        MAGIC = "Calculation completed."
        run_completed = False

        with open(filename) as fh:
            in_event, s = False, None
            for l, line in enumerate(fh):
                if MAGIC in line:
                    run_completed = True

                if not in_event:
                    if line.startswith(START_TAG):
                        in_event = True
                        if s is None:
                            # First event found.
                            s = ""
                        else:
                            # Parse the previous string to generate the appropriate event.
                            events.append(Event.from_string(s, lineno))
                        lineno = l
                else:
                    if line.startswith(STOP_TAG):
                        in_event = False
                    else:
                        # Add new line.
                        s += line

            if s:
                events.append(Event.from_string(s, lineno))

            events.set_run_completed(run_completed)
            return events

    def report_exception(self, filename, exc):
        """
        This method is used when self.parser raises an Exception so that
        we can report a customized `EventReport` object with info the exception.
        """
        return EventReport(filename, events=[Error(str(exc))])
