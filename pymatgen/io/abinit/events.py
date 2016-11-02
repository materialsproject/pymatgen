# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
"""
This module defines the events signaled by abinit during the execution. It also
provides a parser to extract these events form the main output file and the log file.
"""
from __future__ import unicode_literals, division, print_function, absolute_import

import sys
import os.path
import datetime
import collections
import yaml
import six
import abc
import logging
import inspect
import numpy as np

from monty.string import indent, is_string, list_strings
from monty.fnmatch import WildCard
from monty.termcolor import colored
from monty.inspect import all_subclasses
from monty.json import MontyDecoder
from pymatgen.core.structure import Structure
from monty.json import MSONable
from pymatgen.serializers.json_coders import pmg_serialize
from .abiinspect import YamlTokenizer

logger = logging.getLogger(__name__)

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
    color = None

    def __init__(self, src_file, src_line, message):
        """
        Basic constructor for :class:`AbinitEvent`.

        Args:
            message: String with human-readable message providing info on the event.
            src_file: String with the name of the Fortran file where the event is raised.
            src_line Integer giving the line number in src_file.
        """
        #print("src_file", src_file, "src_line", src_line)
        self.message = message
        self.src_file = src_file
        self.src_line = src_line

    @pmg_serialize
    def as_dict(self):
        # This is needed because the events printed in the main output file do not define scr_file and src_line
        src_file = getattr(self, "src_file", "Unknown")
        src_line = getattr(self, "src_line", 0)
        return dict(message=self.message, src_file=src_file, src_line=src_line, yaml_tag=self.yaml_tag)

    @classmethod
    def from_dict(cls, d):
        cls = as_event_class(d.get("yaml_tag"))
        return cls(**{k: v for k, v in d.items() if k != "yaml_tag" and not k.startswith("@")})

    @property
    def header(self):
        try:
            return "<%s at %s:%s>" % (self.name, self.src_file, self.src_line)
        except AttributeError:
            # This is needed because the events printed in the main output file do not define scr_file and src_line
            return "<%s at %s:%s>" % (self.name, "Unknown", 0)

    def __repr__(self):
        return self.header

    def __str__(self):
        return "\n".join((self.header, self.message))

    def __eq__(self, other):
        if other is None: return False
        return self.message == other.message

    def __ne__(self, other):
        return not self.__eq__(other)

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

        raise ValueError("Cannot determine the base class of %s" % self.__class__.__name__)

    def correct(self, task):
        """
        This method is called when an error is detected in a :class:`Task`
        It should perform any corrective measures relating to the detected error.
        The idea is similar to the one used in custodian but the handler receives
        a :class:`Task` object so that we have access to its methods.

        Returns:
        (dict) JSON serializable dict that describes the errors and actions taken. E.g.
        {"errors": list_of_errors, "actions": list_of_actions_taken}.
        If this is an unfixable error, actions should be set to None.
        """
        return 0


class AbinitComment(AbinitEvent):
    """Base class for Comment events"""
    yaml_tag = '!COMMENT'
    color = "blue"


class AbinitError(AbinitEvent):
    """Base class for Error events"""
    yaml_tag = '!ERROR'
    color = "red"


class AbinitYamlError(AbinitError):
    """
    Raised if the YAML parser cannot parse the document and the doc tag is an Error.
    It's an AbinitError because the msg produced by the code is not valid YAML!
    """


class AbinitBug(AbinitEvent):
    """Base class for Bug events"""
    yaml_tag = '!BUG'
    color = "red"


class AbinitWarning(AbinitEvent):
    """
    Base class for Warning events (the most important class).
    Developers should subclass this class to define the different exceptions
    raised by the code and the possible actions that can be performed.
    """
    yaml_tag = '!WARNING'
    color = "magenta"


class AbinitCriticalWarning(AbinitWarning):
    color = "red"


class AbinitYamlWarning(AbinitCriticalWarning):
    """
    Raised if the YAML parser cannot parse the document and the doc tas is a Warning.
    """

###############################
# Warnings triggering restart #
###############################

class ScfConvergenceWarning(AbinitCriticalWarning):
    """Warning raised when the GS SCF cycle did not converge."""
    yaml_tag = '!ScfConvergenceWarning'


class NscfConvergenceWarning(AbinitCriticalWarning):
    """Warning raised when the GS NSCF cycle did not converge."""
    yaml_tag = '!NscfConvergenceWarning'


class RelaxConvergenceWarning(AbinitCriticalWarning):
    """Warning raised when the structural relaxation did not converge."""
    yaml_tag = '!RelaxConvergenceWarning'


# TODO: for the time being we don't discern between GS and PhononCalculations.
#class PhononConvergenceWarning(AbinitCriticalWarning):
#    """Warning raised when the phonon calculation did not converge."""
#    yaml_tag = u'!PhononConvergenceWarning'


class QPSConvergenceWarning(AbinitCriticalWarning):
    """Warning raised when the QPS iteration (GW) did not converge."""
    yaml_tag = '!QPSConvergenceWarning'


class HaydockConvergenceWarning(AbinitCriticalWarning):
    """Warning raised when the Haydock method (BSE) did not converge."""
    yaml_tag = '!HaydockConvergenceWarning'


# Error classes providing a correct method.

# Register the concrete base classes.
_BASE_CLASSES = [
    AbinitComment,
    AbinitError,
    AbinitBug,
    AbinitWarning,
]


class EventReport(collections.Iterable, MSONable):
    """
    Iterable storing the events raised by an ABINIT calculation.

    Attributes::

        stat: information about a file as returned by os.stat
    """
    def __init__(self, filename, events=None):
        """
        List of ABINIT events.

        Args:
            filename: Name of the file
            events: List of Event objects
        """
        self.filename = os.path.abspath(filename)
        self.stat = os.stat(self.filename)
        self.start_datetime, self.end_datetime = None, None

        self._events = []
        self._events_by_baseclass = collections.defaultdict(list)

        if events is not None:
            for ev in events:
                self.append(ev)

    def __len__(self):
        return len(self._events)

    def __iter__(self):
        return self._events.__iter__()

    def __getitem__(self, slice):
        return self._events[slice]

    def __str__(self):
        #has_colours = stream_has_colours(stream)
        has_colours = True

        lines = []
        app = lines.append

        app("Events found in %s\n" % self.filename)
        for i, event in enumerate(self):
            if has_colours:
                app("[%d] %s" % (i+1, colored(event.header, color=event.color)))
                app(indent(event.message, 4))
            else:
                app("[%d] %s" % (i+1, str(event)))

        app("num_errors: %s, num_warnings: %s, num_comments: %s, completed: %s\n" % (
            self.num_errors, self.num_warnings, self.num_comments, self.run_completed))

        return "\n".join(lines)

    def append(self, event):
        """Add an event to the list."""
        self._events.append(event)
        self._events_by_baseclass[event.baseclass].append(event)

    def set_run_completed(self, boolean, start_datetime, end_datetime):
        """Set the value of _run_completed."""
        self._run_completed = boolean

        if (start_datetime, end_datetime) != (None, None):
            # start_datetime: Sat Feb 28 23:54:27 2015
            # end_datetime: Sat Feb 28 23:54:30 2015
            try:
                fmt = "%a %b %d %H:%M:%S %Y"
                self.start_datetime = datetime.datetime.strptime(start_datetime, fmt)
                self.end_datetime = datetime.datetime.strptime(end_datetime, fmt)
            except Exception as exc:
                # Maybe LOCALE != en_US
                logger.warning(str(exc))

    @property
    def run_etime(self):
        """Wall-time of the run as `timedelta` object."""
        if self.start_datetime is None or self.end_datetime is None:
            return None

        return self.end_datetime - self.start_datetime

    @property
    def run_completed(self):
        """True if the calculation terminated."""
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
        """List of errors + bugs found."""
        return self.select(AbinitError) + self.select(AbinitBug)

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
        """
        return self._events_by_baseclass[base_class]

    def filter_types(self, event_types):
        events = []
        for ev in self:
            if type(ev) in event_types: events.append(ev)
        return self.__class__(filename=self.filename, events=events)

    def get_events_of_type(self, event_class):
        """Return a list of events of the given class."""
        return [ev for ev in self if type(ev) == event_class]

    @pmg_serialize
    def as_dict(self):
        return dict(filename=self.filename, events=[e.as_dict() for e in self._events])

    @classmethod
    def from_dict(cls, d):
        return cls(filename=d["filename"], events=[AbinitEvent.from_dict(e) for e in d["events"]])


class EventsParserError(Exception):
    """Base class for the exceptions raised by :class:`EventsParser`."""


class EventsParser(object):
    """
    Parses the output or the log file produced by ABINIT and extract the list of events.
    """
    Error = EventsParserError

    def parse(self, filename, verbose=0):
        """
        Parse the given file. Return :class:`EventReport`.
        """
        run_completed, start_datetime, end_datetime = False, None, None
        filename = os.path.abspath(filename)
        report = EventReport(filename)

        w = WildCard("*Error|*Warning|*Comment|*Bug|*ERROR|*WARNING|*COMMENT|*BUG")

        with YamlTokenizer(filename) as tokens:
            for doc in tokens:
                if w.match(doc.tag):
                    #print("got doc.tag", doc.tag,"--")
                    try:
                        #print(doc.text)
                        event = yaml.load(doc.text)
                        #print(event.yaml_tag, type(event))
                    except:
                        #raise
                        # Wrong YAML doc. Check tha doc tag and instantiate the proper event.
                        message = "Malformatted YAML document at line: %d\n" % doc.lineno
                        message += doc.text

                        # This call is very expensive when we have many exceptions due to malformatted YAML docs.
                        if verbose:
                            message += "Traceback:\n %s" % straceback()

                        if "error" in doc.tag.lower():
                            print("It seems an error", doc.tag)
                            event = AbinitYamlError(message=message, src_file=__file__, src_line=0)
                        else:
                            event = AbinitYamlWarning(message=message, src_file=__file__, src_line=0)

                    event.lineno = doc.lineno
                    report.append(event)

                # Check whether the calculation completed.
                if doc.tag == "!FinalSummary":
                    #print(doc)
                    run_completed = True
                    d = doc.as_dict()
                    #print(d)
                    start_datetime, end_datetime = d["start_datetime"], d["end_datetime"]

        report.set_run_completed(run_completed, start_datetime, end_datetime)
        return report

    def report_exception(self, filename, exc):
        """
        This method is used when self.parser raises an Exception so that
        we can report a customized :class:`EventReport` object with info the exception.
        """
        # Build fake event.
        event = AbinitError(src_file="Unknown", src_line=0, message=str(exc))
        return EventReport(filename, events=[event])


class EventHandler(six.with_metaclass(abc.ABCMeta, MSONable, object)):
    """
    Abstract base class defining the interface for an EventHandler.

    The__init__ should always provide default values for its arguments so that we can
    easily instantiate the handlers with:

        handlers = [cls() for cls in get_event_handler_classes()]

    The defaul values should be chosen so to cover the most typical cases.

    Each EventHandler should define the class attribute `can_change_physics`
    that is true if the handler changes `important` parameters of the
    run that are tightly connected to the physics of the system.

    For example, an `EventHandler` that changes the value of `dilatmx` and
    prepare the restart is not changing the physics. Similarly a handler
    that changes the mixing algorithm. On the contrary, a handler that
    changes the value of the smearing is modifying an important physical
    parameter, and the user should be made aware of this so that
    there's an explicit agreement between the user and the code.

    The default handlers are those that do not change the physics,
    other handlers can be installed by the user when constructing with the flow with

        TODO

    .. warning::

        The EventHandler should perform any action at the level of the input files
        needed to solve the problem and then prepare the task for a new submission
        The handler should never try to resubmit the task. The submission must be
        delegated to the scheduler or Fireworks.
    """

    event_class = AbinitEvent
    """AbinitEvent subclass associated to this handler."""

    #can_change_physics

    FIXED = 1
    NOT_FIXED = 0

    def __init__(self):
        """Simple init for compatibility with introspection in as_dict/from_dict"""
        return super(EventHandler,self).__init__()

    @classmethod
    def cls2str(cls):
        lines = []
        app = lines.append

        ecls = cls.event_class
        app("event name = %s" % ecls.yaml_tag)
        app("event documentation: ")
        lines.extend(ecls.__doc__.split("\n"))
        app("handler documentation: ")
        lines.extend(cls.__doc__.split("\n"))

        return "\n".join(lines)

    def __str__(self):
        return "<%s>" % self.__class__.__name__

    def can_handle(self, event):
        """True if this handler is associated to the given :class:`AbinitEvent`"""
        return self.event_class == event.__class__

    # TODO: defined CorrectionRecord object and provide helper functions to build it

    def count(self, task):
        """
        Return the number of times the event associated to this handler
        has been already fixed in the :class:`Task`.
        """
        return len([c for c in task.corrections if c["event"]["@class"] == self.event_class])

    @abc.abstractmethod
    def handle_task_event(self, task, event):
        """
        Method to handle Abinit events.

        Args:
            task: :class:`Task` object.
            event: :class:`AbinitEvent` found in the log file.

        Return:
            0 if no action has been applied, 1 if the problem has been fixed.
        """

    @pmg_serialize
    def as_dict(self):
        #@Guido this introspection is nice but it's not safe
        d = {}
        if hasattr(self, "__init__"):
            for c in inspect.getargspec(self.__init__).args:
                if c != "self":
                    d[c] = self.__getattribute__(c)
        return d

    @classmethod
    def from_dict(cls, d):
        kwargs = {k: v for k, v in d.items() if k in inspect.getargspec(cls.__init__).args}
        return cls(**kwargs)

    @classmethod
    def compare_inputs(cls, new_input, old_input):

        def vars_dict(d):
            """
            make a simple dictionary and convert numpy arrays to lists
            """
            new_d = {}
            for key, value in d.items():
                if isinstance(value, np.ndarray): value = value.tolist()
                new_d[key] = value

            return new_d

        new_vars = vars_dict(new_input)
        old_vars = vars_dict(old_input)

        new_keys = set(new_vars.keys())
        old_keys = set(old_vars.keys())
        intersect = new_keys.intersection(old_keys)

        added_keys = new_keys - intersect
        removed_keys = old_keys - intersect
        changed_keys = set(v for v in intersect if new_vars[v] != old_vars[v])

        log_diff = {}
        if added_keys:
            log_diff['_set'] = {k: new_vars[k] for k in added_keys}

        if changed_keys:
            log_diff['_update'] = ({k: {'new': new_vars[k], 'old': old_vars[k]} for k in changed_keys})

        if new_input.structure != old_input.structure:
            log_diff['_change_structure'] = new_input.structure.as_dict()

        if removed_keys:
            log_diff['_pop'] = {k: old_vars[k] for k in removed_keys}

        return log_diff


class Correction(MSONable):

    def __init__(self, handler, actions, event, reset=False):
        self.handler = handler
        self.actions = actions
        self.event = event
        self.reset = reset

    @pmg_serialize
    def as_dict(self):
        return dict(handler=self.handler.as_dict(), actions=self.actions, event=self.event.as_dict(), reset=self.reset)

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        return cls(handler=dec.process_decoded(d['handler']), actions=d['actions'],
                   event=dec.process_decoded(d['event']), reset=d['reset'])


#class WarningHandler(EventHandler):
#    """Base class for handlers associated to ABINIT warnings."""
#    event_class = AbinitWarning
#
#class BugHandler(EventHandler):
#    """Base class for handlers associated to ABINIT bugs."""
#    event_class = AbinitBug


class ErrorHandler(EventHandler):
    """Base class for handlers associated to ABINIT errors."""
    event_class = AbinitError

_ABC_EVHANDLER_CLASSES = set([ErrorHandler,])


# Public API
def autodoc_event_handlers(stream=sys.stdout):
    """
    Print to the given string, the documentation for the events
    and the associated handlers.
    """
    lines = []
    for cls in all_subclasses(EventHandler):
        if cls in _ABC_EVHANDLER_CLASSES: continue
        event_class = cls.event_class
        lines.extend(cls.cls2str().split("\n"))

        # Here we enforce the abstract protocol of the class
        # The unit test in tests_events will detect the problem.
        if not hasattr(cls, "can_change_physics"):
            raise RuntimeError("%s: can_change_physics must be defined" % cls)

    stream.write("\n".join(lines) + "\n")


def get_event_handler_classes(categories=None):
    """Return the list of handler classes."""
    classes = [c for c in all_subclasses(EventHandler) if c not in _ABC_EVHANDLER_CLASSES]
    return classes


def as_event_class(obj):
    """
    Convert obj into a subclass of AbinitEvent.
    obj can be either a class or a string with the class name or the YAML tag
    """
    if is_string(obj):
        for c in all_subclasses(AbinitEvent):
            if c.__name__ == obj or c.yaml_tag == obj: return c
        raise ValueError("Cannot find event class associated to %s" % obj)

    # Assume class.
    assert obj in all_subclasses(AbinitEvent)
    return obj


############################################
########## Concrete classes ################
############################################

class DilatmxError(AbinitError):
    """
    This Error occurs in variable cell calculations when the increase in the
    unit cell volume is too large.
    """
    yaml_tag = '!DilatmxError'


class DilatmxErrorHandler(ErrorHandler):
    """
    Handle DilatmxError. Abinit produces a netcdf file with the last structure before aborting
    The handler changes the structure in the input with the last configuration and modify the value of dilatmx.
    """
    event_class = DilatmxError

    can_change_physics = False

    def __init__(self, max_dilatmx=1.3):
        self.max_dilatmx = max_dilatmx

    def handle_task_event(self, task, event):
        # Read the last structure dumped by ABINIT before aborting.
        filepath = task.outdir.has_abiext("DILATMX_STRUCT.nc")
        last_structure = Structure.from_file(filepath)

        task._change_structure(last_structure)

        #read the suggested dilatmx
        # new_dilatmx = 1.05
        # if new_dilatmx > self.max_dilatmx:
        #     msg = "Suggested dilatmx ({}) exceeds maximux configured value ({}).".format(new_dilatmx, self.max_dilatmx)
        #     return self.NOT_FIXED
        # task.strategy.abinit_input.set_vars(dilatmx=new_dilatmx)
        msg = "Take last structure from DILATMX_STRUCT.nc, will try to restart with dilatmx %s" % task.get_inpvar("dilatmx")
        task.log_correction(event, msg)
        # Note that we change the structure but we don't try restart from the previous WFK|DEN file
        # because Abinit called mpi_abort and therefore no final WFK|DEN file has been produced.

        return self.FIXED

    def handle_input_event(self, abiinput, outdir, event):
        try:
            old_abiinput = abiinput.deepcopy()
            # Read the last structure dumped by ABINIT before aborting.
            filepath = outdir.has_abiext("DILATMX_STRUCT.nc")
            last_structure = Structure.from_file(filepath)
            abiinput.set_structure(last_structure)
            #FIXME restart from DEN files not always working with interpolation
            return Correction(self, self.compare_inputs(abiinput, old_abiinput), event, reset=True)
            # return Correction(self, self.compare_inputs(abiinput, old_abiinput), event, event=False)
        except Exception as exc:
            logger.warning('Error while trying to apply the handler {}.'.format(str(self)), exc)
            return None


class TolSymError(AbinitError):
    """
    Class of errors raised by Abinit when it cannot detect the symmetries of the system.
    The handler assumes the structure makes sense and the error is just due to numerical inaccuracies.
    We increase the value of tolsym in the input file (default 1-8) so that Abinit can find the space group
    and re-symmetrize the input structure.
    """
    yaml_tag = '!TolSymError'


class TolSymErrorHandler(ErrorHandler):
    """
    Increase the value of tolsym in the input file.
    """
    event_class = TolSymError

    can_change_physics = False

    def __init__(self, max_nfixes=3):
        self.max_nfixes = max_nfixes

    def handle_task_event(self, task, event):
        # TODO: Add limit on the number of fixes one can do for the same error
        # For example in this case, the scheduler will stop after 20 submissions
        if self.count(task) > self.max_nfixes:
            return self.NOT_FIXED

        old_tolsym = task.get_inpvar("tolsym")
        new_tolsym = 1e-6 if old_tolsym is None else old_tolsym * 10
        task.set_vars(tolsym=new_tolsym)

        task.log_correction(event, "Increasing tolsym from %s to %s" % (old_tolsym, new_tolsym))
        return self.FIXED

    def handle_input_event(self, abiinput, outdir, event):
        try:
            old_abiinput = abiinput.deepcopy()
            old_tolsym = abiinput["tolsym"]
            new_tolsym = 1e-6 if old_tolsym is None else old_tolsym * 10
            abiinput.set_vars(tolsym=new_tolsym)
            return Correction(self, self.compare_inputs(abiinput, old_abiinput), event, reset=False)
        except Exception as exc:
            logger.warning('Error while trying to apply the handler {}.'.format(str(self)), exc)
            return None


class MemanaError(AbinitError):
    """
    Class of errors raised by the memory analyzer.
    (the section that estimates the memory requirements from the input parameters).
    """
    yaml_tag = '!MemanaError'


class MemanaErrorHandler(ErrorHandler):
    """
    Set mem_test to 0 to bypass the memory check.
    """
    event_class = MemanaError

    can_change_physics = False

    def handle_task_event(self, task, event):
        task.set_vars(mem_test=0)
        task.log_correction(event, "Find MemanaError. Setting mem_test to 0 in input file.")
        return self.FIXED

    def handle_input_event(self, abiinput, outdir, event):
        try:
            old_abiinput = abiinput.deepcopy()
            abiinput.set_vars(mem_test=0)
            return Correction(self, self.compare_inputs(abiinput, old_abiinput), event, reset=False)
        except Exception as exc:
            logger.warning('Error while trying to apply the handler {}.'.format(str(self)), exc)
            return None


class MemoryError(AbinitError):
    """
    This error occurs when a checked allocation fails in Abinit
    The only way to go is to increase memory
    """
    yaml_tag = '!MemoryError'


class MemoryErrorHandler(ErrorHandler):
    """
    Handle MemoryError. Increase the resources requirements
    """
    event_class = MemoryError

    can_change_physics = False

    def handle_task_event(self, task, event):

      task.manager.increase_resources()
      return self.FIXED

    def handle_input_event(self, abiinput, outdir, event):
      """
      Shouldn't do anything on the input
      """
      return None
