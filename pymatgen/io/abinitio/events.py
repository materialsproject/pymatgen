"""
This module defines the events signaled by abinit during the execution. It also
provides a parser to extract these events form the main output file and the log file.
"""
from __future__ import division, print_function

import os.path
import collections
import yaml
from pymatgen.io.abinitio import myaml

from pymatgen.util.string_utils import WildCard
from pymatgen.io.abinitio.abiinspect import YamlTokenizer, YamlDoc

__all__ = [
    "EventsParser",
]

def straceback():
    """Returns a string with the traceback."""
    import traceback
    return traceback.format_exc()


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

    1) We use YamlTokenizer to extract the documents from the output file 
    2) If we have a tag that ends with "Warning", "Error", "Bug", "Comment
       we know we have encountered a new ABINIT event 
    3) We parse the document with yaml.load(doc.text) and we get the object

    Note that:
        # --- and ... become reserved words (whey they are placed at
          the begining of a line) since they are used to mark the beginning and 
          the end of YAML documents.

        # All the possible events should subclass `AbinitEvent` and define 
          the class attribute yaml_tag so that yaml.load will know how to 
          build the instance.
    """
    def __init__(self, message, src_file, src_line):
        """
        Basic constructor for `AbinitEvent`. 

        Args:
            message:
                String with human-readable message providing info on the event.
            src_file:
                String with the name of the Fortran file where the event is raised.
            src_line
                Integer giving the line number in src_file.
        """
        self.message = message
        self.src_file = src_file
        self.src_line = src_line

    def __str__(self):
        header = "%s:%s" % (self.src_file, self.src_line)
        return "\n".join((header, self.message))

    @property
    def name(self):
        """Name of the event (class name)"""
        return self.__class__.__name__

    @property
    def baseclass(self):
        """The baseclass of self."""
        for cls in _BASE_CLASSES:
            if isinstance(self, cls):
                return cls

        err_msg = "Cannot determine the base class of %s" % (self.__class__.__name__)
        raise ValueError(err_msg)

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


class AbinitYamlError(AbinitError):
    """Raised if the YAML parser cannot parse the document and the doc tag is an Error."""


class AbinitBug(AbinitEvent):
    """Base class for Bug events"""
    yaml_tag = u'!BUG'


class AbinitWarning(AbinitEvent):
    """
    Base class for Warning events (the most important class).
    Developers should subclass this class to define the different exceptions
    raised by the code and the possible actions that can be performed.
    """
    yaml_tag = u'!WARNING'


class AbinitYamlWarning(AbinitWarning):
    """
    Raised if the YAML parser cannot parse the document and the doc tas is a Warning.
    """


class ScfConvergenceWarning(AbinitWarning):
    """Warning raised when the GS SCF cycle did not converge."""
    yaml_tag = u'!ScfConvergenceWarning'


class NscfConvergenceWarning(AbinitWarning):
    """Warning raised when the GS NSCF cycle did not converge."""
    yaml_tag = u'!NscfConvergenceWarning'


class RelaxConvergenceWarning(AbinitWarning):
    """Warning raised when the structural relaxation did not converge."""
    yaml_tag = u'!RelaxConvergenceWarning'


# TODO: for the time being we don't discern between GS and PhononCalculations.
#class PhononConvergenceWarning(AbinitWarning):
#    """Warning raised when the phonon calculation did not converge."""
#    yaml_tag = u'!PhononConvergenceWarning'


class QPSConvergenceWarning(AbinitWarning):
    """Warning raised when the QPS iteration (GW) did not converge."""
    yaml_tag = u'!QPSConvergenceWarning'


class HaydockConvergenceWarning(AbinitWarning):
    """Warning raised when the Haydock method (BSE) did not converge."""
    yaml_tag = u'!HaydockConvergenceWarning'


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
        lines = []
        app = lines.append

        app("Event Report for file: %s" % self.filename)
        for i, event in enumerate(self):
            app("%d) %s" % (i+1, str(event)))

        app("num_errors: %s, num_warnings: %s, num_comments: %s" % (
            self.num_errors, self.num_warnings, self.num_comments))
        app("run_completed: %s" % self.run_completed)

        return "\n".join(lines)

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

    def select(self, base_class):
        """
        Return the list of events that inherits from class base_class

        Args:
            only_critical:
                if True, only critical events are returned.
        """
        return self._events_by_baseclass[base_class][:]

    def filter_types(self, event_types):
        evts = []
        for event in self:
            if type(event) in event_types:
                evts.append(event)
        return evts


class EventsParserError(Exception):
    """Base class for the exceptions raised by EventsParser."""


class EventsParser(object):
    """
    Parses the output or the log file produced by abinit and extract the list of events.
    """
    Error = EventsParserError

    @staticmethod
    def parse(filename):
        """
        This is the new parser, it will be used when we implement
        the new format in abinit.
        """
        run_completed = False
        filename = os.path.abspath(filename)
        report = EventReport(filename)

        # TODO Use CamelCase for the Fortran messages.
        w = WildCard("*Error|*Warning|*Comment|*ERROR|*WARNING|*COMMENT")

        with YamlTokenizer(filename) as tokens:

            for doc in tokens:
                #print(80*"*")
                #print("doc.tag", doc.tag)
                #print("doc",doc)
                #print(80*"*")
                if w.match(doc.tag):
                    #print("got doc.tag", doc.tag,"--")
                    try:
                        event = myaml.load(doc.text)
                    except:
                        # Wrong YAML doc. Check tha doc tag and instantiate the proper event.
                        message = "Malformatted YAML document at line: %d\n" % doc.lineno
                        message += doc.text
                        # This call is very expensive when we have many exceptions due to malformatted YAML docs. 
                        #message += "Traceback:\n %s" % straceback()

                        if "error" in doc.tag.lower():
                            print("It seems an error",doc.tag)
                            event = AbinitYamlError(message=message, src_file=__file__, src_line=0)
                        else:
                            event = AbinitYamlWarning(message=message, src_file=__file__, src_line=0)

                    event.lineno = doc.lineno
                    report.append(event)

                # Check whether the calculation completed.
                if doc.tag == "!FinalSummary":
                    run_completed = True

        # TODO: Add YAML doc with FinalSummary.
        #MAGIC = "Calculation completed."
        #with open(filename) as fh:
        #    for line in fh:
        #        if MAGIC in line:
        #            run_completed = True
        #            break

        report.set_run_completed(run_completed)
        return report

    def report_exception(self, filename, exc):
        """
        This method is used when self.parser raises an Exception so that
        we can report a customized `EventReport` object with info the exception.
        """
        return EventReport(filename, events=[Error(str(exc))])

