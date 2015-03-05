# coding: utf-8
"""
This module defines the Node class that is inherited by Task, Work and Flow objects.
"""
from __future__ import division, print_function, unicode_literals

import sys
import os
import time
import collections
import abc
import six

from pprint import pprint
from atomicfile import AtomicFile
from monty.termcolor import colored
from monty.serialization import loadfn
from monty.string import is_string
from monty.io import FileLock
from monty.collections import AttrDict, Namespace
from monty.functools import lazy_property
from pymatgen.serializers.json_coders import PMGSONable, json_pretty_dump, pmg_serialize
from .utils import File, Directory, irdvars_for_ext, abi_extensions


import logging
logger = logging.getLogger(__name__)


__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


def _2attrs(item):
    return item if item is None or isinstance(list, tuple) else (item,)


class Status(int):
    """This object is an integer representing the status of the `Node`."""

    # Possible status of the node. See monty.termocolor for the meaning of color, on_color and attrs.
    _STATUS_INFO = [
        #(value, name, color, on_color, attrs)
        (1,  "Initialized",   None     , None, None),         # Node has been initialized
        (2,  "Locked",        None     , None, None),         # Task is locked an must be explicitly unlocked by an external subject (Work).
        (3,  "Ready",         None     , None, None),         # Node is ready i.e. all the depencies of the node have status S_OK
        (4,  "Submitted",     "blue"   , None, None),         # Node has been submitted (The `Task` is running or we have started to finalize the Work)
        (5,  "Running",       "magenta", None, None),         # Node is running.
        (6,  "Done",          None     , None, None),         # Node done, This does not imply that results are ok or that the calculation completed successfully
        (7,  "AbiCritical",   "red"    , None, None),         # Node raised an Error by ABINIT.
        (8,  "QCritical",     "red"    , "on_white", None),   # Node raised an Error by submitting submission script, or by executing it
        (9,  "Unconverged",   "red"    , "on_yellow", None),  # This usually means that an iterative algorithm didn't converge.
        (10, "Error",         "red"    , None, None),         # Node raised an unrecoverable error, usually raised when an attempt to fix one of other types failed.
        (11, "Completed",     "green"  , None, None),         # Execution completed successfully.
        #(11, "Completed",     "green"  , None, "underline"),   
    ]
    _STATUS2STR = collections.OrderedDict([(t[0], t[1]) for t in _STATUS_INFO])
    _STATUS2COLOR_OPTS = collections.OrderedDict([(t[0], {"color": t[2], "on_color": t[3], "attrs": _2attrs(t[4])}) for t in _STATUS_INFO])

    def __repr__(self):
        return "<%s: %s, at %s>" % (self.__class__.__name__, str(self), id(self))

    def __str__(self):
        """String representation."""
        return self._STATUS2STR[self]

    @classmethod
    def as_status(cls, obj):
        """Convert obj into Status."""
        return obj if isinstance(obj, cls) else cls.from_string(obj)

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in cls._STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)

    @property
    def is_critical(self):
        """True if status is critical."""
        return str(self) in ("AbiCritical", "QCritical", "Uncoverged", "Error") 

    @property
    def color_opts(self):
        return self._STATUS2COLOR_OPTS[self]

    @property
    def colored(self):
        """Return colorized text used to print the status if the stream supports it."""
        return colored(str(self), **self.color_opts) 


class Dependency(object):
    """
    This object describes the dependencies among the nodes of a calculation.

    A `Dependency` consists of a `Node` that produces a list of products (files) 
    that are used by the other nodes (`Task` or `Work`) to start the calculation.
    One usually creates the object by calling work.register 

    Example:

        # Register the SCF task in work.
        scf_task = work.register(scf_strategy)

        # Register the NSCF calculation and its dependency on the SCF run via deps.
        nscf_task = work.register(nscf_strategy, deps={scf_task: "DEN"})
    """
    def __init__(self, node, exts=None):
        """
        Args:
            node: The task or the worfklow associated to the dependency or string with a filepath.
            exts: Extensions of the output files that are needed for running the other tasks.
        """
        self._node = Node.as_node(node)

        if exts and is_string(exts):
            exts = exts.split()

        # Extract extensions.
        self.exts = [e for e in exts if not e.startswith("@")]

        # Save getters
        self.getters = [e for e in exts if e.startswith("@")]
        #if self.getters: print(self.getters)

    def __hash__(self):
        return hash(self._node)

    def __repr__(self):
        return "node %s will produce: %s " % (repr(self.node), repr(self.exts))

    def __str__(self):
        return "node %s will produce: %s " % (str(self.node), str(self.exts))

    @property
    def info(self):
        return str(self.node)

    @property
    def node(self):
        """The :class:`Node` associated to the dependency."""
        return self._node

    @property
    def status(self):
        """The status of the dependency, i.e. the status of the :class:`Node`."""
        return self.node.status

    @lazy_property
    def products(self):
        """List of output files produces by self."""
        _products = []
        for ext in self.exts:
            prod = Product(ext, self.node.opath_from_ext(ext))
            _products.append(prod)

        return _products

    def apply_getters(self, task):
        if not self.getters: return

        for getter in self.getters:
            if getter == "@structure":
                new_structure = self.node.get_final_structure()
                task._change_structure(new_structure)
            else:
                raise ValueError("Wrong getter %s" % getter)

    def connecting_vars(self):
        """
        Returns a dictionary with the variables that must be added to the 
        input file in order to connect this :class:`Node` to its dependencies.
        """
        vars = {}
        for prod in self.products:
            vars.update(prod.connecting_vars())

        return vars

    def get_filepaths_and_exts(self):
        """Returns the paths of the output files produced by self and its extensions"""
        filepaths = [prod.filepath for prod in self.products]
        exts = [prod.ext for prod in self.products]

        return filepaths, exts


class Product(object):
    """
    A product represents an output file produced by ABINIT instance.
    This file is needed to start another `Task` or another `Work`.
    """
    def __init__(self, ext, path):
        """
        Args:
            ext: ABINIT file extension
            path: (asbolute) filepath
        """
        if ext not in abi_extensions():
            raise ValueError("Extension %s has not been registered in the internal database" % str(ext))

        self.ext = ext
        self.file = File(path)

    @classmethod
    def from_file(cls, filepath):
        """Build a :class:`Product` instance from a filepath."""
        # Find the abinit extension.
        for i in range(len(filepath)):
            if filepath[i:] in abi_extensions():
                ext = filepath[i:]
                break
        else:
            raise ValueError("Cannot detect abinit extension in %s" % filepath)
        
        return cls(ext, filepath)

    def __str__(self):
        return "File=%s, Extension=%s, " % (self.file.path, self.ext)

    @property
    def filepath(self):
        """Absolute path of the file."""
        return self.file.path

    def connecting_vars(self):
        """
        Returns a dictionary with the ABINIT variables that 
        must be used to make the code use this file.
        """
        return irdvars_for_ext(self.ext)


class GridFsFile(AttrDict):
    """Information on a file that will stored in the MongoDb gridfs collection."""
    def __init__(self, path, fs_id=None, mode="b"):
        super(GridFsFile, self).__init__(path=path, fs_id=fs_id, mode=mode)


class NodeResults(dict, PMGSONable):
    """Dictionary used to store the most important results produced by a :class:`Node`."""
    JSON_SCHEMA = {
        "type": "object",
        "properties": {
            "node_id": {"type": "integer", "required": True},
            "node_finalized": {"type": "boolean", "required": True},
            "node_history": {"type": "array", "required": True},
            "node_class": {"type": "string", "required": True},
            "node_name": {"type": "string", "required": True},
            "node_status": {"type": "string", "required": True},
            "in": {"type": "object", "required": True, "description": "dictionary with input parameters"},
            "out": {"type": "object", "required": True, "description": "dictionary with the output results"},
            "exceptions": {"type": "array", "required": True},
            "files": {"type": "object", "required": True},
        },
    }

    @classmethod
    def from_node(cls, node):
        """Initialize an instance of `NodeResults` from a `Node` subclass."""
        kwargs = dict(
            node_id=node.node_id,
            node_finalized=node.finalized,
            node_history=list(node.history),
            node_name=node.name, 
            node_class=node.__class__.__name__,
            node_status=str(node.status),
        )

        return node.Results(node, **kwargs)

    def __init__(self, node, **kwargs):
        super(NodeResults, self).__init__(**kwargs)
        self.node = node

        if "in" not in self: self["in"] = Namespace()
        if "out" not in self: self["out"] = Namespace()
        if "exceptions" not in self: self["exceptions"] = []
        if "files" not in self: self["files"] = Namespace()

    @property
    def exceptions(self):
        return self["exceptions"]

    @property
    def gridfs_files(self):
        """List with the absolute paths of the files to be put in GridFs."""
        return self["files"]

    def register_gridfs_files(self, **kwargs):
        """
        This function registers the files that will be saved in GridFS.
        kwargs is a dictionary mapping the key associated to the file (usually the extension)
        to the absolute path. By default, files are assumed to be in binary form, for formatted files
        one should pass a tuple ("filepath", "t").

        Example::

            results.register_gridfs(GSR="path/to/GSR.nc", text_file=("/path/to/txt_file", "t"))

        The GSR file is a binary file, whereas text_file is a text file.
        """
        d = {}
        for k, v in kwargs.items():
            mode = "b" 
            if isinstance(v, (list, tuple)): v, mode = v
            d[k] = GridFsFile(path=v, mode=mode)
            
        self["files"].update(d)
        return self

    def push_exceptions(self, *exceptions):
        for exc in exceptions:
            newstr = str(exc)
            if newstr not in self.exceptions:
                self["exceptions"] += [newstr,]

    @pmg_serialize
    def as_dict(self):
        return self.copy()
                                                                                
    @classmethod
    def from_dict(cls, d):
        return cls({k: v for k, v in d.items() if k not in ("@module", "@class")})

    def json_dump(self, filename):
        json_pretty_dump(self.as_dict(), filename)

    @classmethod
    def json_load(cls, filename):
        return cls.from_dict(loadfn(filename))

    def validate_json_schema(self):
        import validictory
        d = self.as_dict()
        try:
            validictory.validate(d, self.JSON_SCHEMA)
            return True
        except ValueError as exc:
            pprint(d)
            print(exc)
            return False

    def update_collection(self, collection):
        """
        Update a mongodb collection.
        """
        node = self.node 
        flow = node if node.is_flow else node.flow

        # Build the key used to store the entry in the document.
        key = node.name
        if node.is_task:
            key = "w" + str(node.pos[0]) + "_t" + str(node.pos[1])
        elif node.is_work:
            key = "w" + str(node.pos)

        db = collection.database

        # Save files with GridFs first in order to get the ID.
        if self.gridfs_files:
            import gridfs
            fs = gridfs.GridFS(db)
            for ext, gridfile in self.gridfs_files.items():
                logger.info("gridfs: about to put file:", str(gridfile))
                # Here we set gridfile.fs_id that will be stored in the mondodb document
                try:
                    with open(gridfile.path, "r" + gridfile.mode) as f:
                        gridfile.fs_id = fs.put(f, filename=gridfile.path)
                except IOError as exc:
                    logger.critical(str(exc))

        if flow.mongo_id is None:
            # Flow does not have a mongo_id, allocate doc for the flow and save its id.
            flow.mongo_id = collection.insert({})
            print("Creating flow.mongo_id", flow.mongo_id, type(flow.mongo_id))

        # Get the document from flow.mongo_id and update it.
        doc = collection.find_one({"_id": flow.mongo_id})
        if key in doc:
            raise ValueError("%s is already in doc!" % key)
        doc[key] = self.as_dict()

        collection.save(doc)
        #collection.update({'_id':mongo_id}, {"$set": doc}, upsert=False)


class Node(six.with_metaclass(abc.ABCMeta, object)):
    """
    Abstract base class defining the interface that must be 
    implemented by the nodes of the calculation.

    Nodes are hashable and can be tested for equality
    """
    Results = NodeResults

    # Possible status of the node.
    S_INIT = Status.from_string("Initialized")
    S_LOCKED = Status.from_string("Locked")
    S_READY = Status.from_string("Ready")
    S_SUB = Status.from_string("Submitted")
    S_RUN = Status.from_string("Running")
    S_DONE = Status.from_string("Done")
    S_ABICRITICAL = Status.from_string("AbiCritical")
    S_QCRITICAL = Status.from_string("QCritical")
    S_UNCONVERGED = Status.from_string("Unconverged")
    S_ERROR = Status.from_string("Error")
    S_OK = Status.from_string("Completed")

    ALL_STATUS = [
        S_INIT,
        S_LOCKED,
        S_READY,
        S_SUB,
        S_RUN,
        S_DONE,
        S_ABICRITICAL,
        S_QCRITICAL,
        S_UNCONVERGED,
        S_ERROR,
        S_OK,
    ]

    def __init__(self):
        # Node identifier.
        self._node_id = get_newnode_id()

        # List of dependencies
        self._deps = []

        # List of files (products) needed by this node.
        self._required_files = []

        # Used to push additional info during the execution. 
        self.history = NodeHistory(maxlen=100)

        # Actions performed to fix abicritical events.
        self._corrections = collections.deque(maxlen=100)

        # Set to true if the node has been finalized.
        self._finalized = False
        self._status = self.S_INIT

    def __eq__(self, other):
        if not isinstance(other, Node): return False
        return self.node_id == other.node_id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.node_id)

    def __repr__(self):
        try:
            return "<%s, node_id=%s, workdir=%s>" % (
                self.__class__.__name__, self.node_id, os.path.relpath(self.workdir))

        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, node_id=%s, workdir=None>" % (self.__class__.__name__, self.node_id)
                                                                                            
    @classmethod
    def as_node(cls, obj):
        """
        Convert obj into a Node instance.

        Return:
            obj if obj is a Node instance,
            cast obj to :class:`FileNode` instance of obj is a string.
            None if obj is None
        """
        if isinstance(obj, cls):
            return obj
        elif is_string(obj):
            # Assume filepath.
            return FileNode(obj)
        elif obj is None:
            return obj
        else:
            raise TypeError("Don't know how to convert %s to Node instance." % obj)

    @property
    def name(self):
        """
        The name of the node 
        (only used for facilitating its identification in the user interface).
        """
        try:
            return self._name
        except AttributeError:
            return os.path.relpath(self.workdir)

    def set_name(self, name):
        """Set the name of the Node."""
        self._name = name

    @property
    def node_id(self):
        """Node identifier."""
        return self._node_id

    def set_node_id(self, node_id):
        """Set the node identifier. Use it carefully!"""
        self._node_id = node_id

    @property
    def finalized(self):
        """True if the `Node` has been finalized."""
        return self._finalized

    @finalized.setter
    def finalized(self, boolean):
        self._finalized = boolean
        self.history.info("Finalized")

    @property
    def corrections(self):
        """
        List of dictionaries with infornation on the actions performed to solve `AbiCritical` Events.
        Each dictionary contains the `AbinitEvent` who triggered the correction and 
        a human-readable message with the description of the operation performed.
        """
        return self._corrections

    @property
    def num_corrections(self):
        return len(self.corrections)

    def log_correction(self, event, action):
        """
        This method should be called once we have fixed the problem associated to this event.
        It adds a new entry in the correction history of the node.

        Args:
            event: `AbinitEvent` that triggered the correction.
            action (str): Human-readable string with info on the action perfomed to solve the problem.
        """
        action = str(action)
        self.history.info(action)

        self._corrections.append(dict(
            event=event.as_dict(), 
            action=action,
        ))

    @property
    def is_file(self):
        """True if this node is a file"""
        return isinstance(self, FileNode)

    @property
    def is_task(self):
        """True if this node is a Task"""
        from .tasks import Task
        return isinstance(self, Task)

    @property
    def is_work(self):
        """True if this node is a Work"""
        from .works import Work
        return isinstance(self, Work)

    @property
    def is_flow(self):
        """True if this node is a Flow"""
        from .flows import Flow
        return isinstance(self, Flow)

    @property
    def deps(self):
        """
        List of :class:`Dependency` objects defining the dependencies 
        of this `Node`. Empty list if this :class:`Node` does not have dependencies.
        """
        return self._deps

    def add_deps(self, deps):
        """
        Add a list of dependencies to the :class:`Node`.

        Args:
            deps: List of :class:`Dependency` objects specifying the dependencies of the node.
                  or dictionary mapping nodes to file extensions e.g. {task: "DEN"}
        """
        if isinstance(deps, collections.Mapping):
            # Convert dictionary into list of dependencies.
            deps = [Dependency(node, exts) for node, exts in deps.items()]

        # We want a list
        if not isinstance(deps, (list, tuple)):
            deps = [deps]

        assert all(isinstance(d, Dependency) for d in deps)

        # Add the dependencies to the node
        self._deps.extend(deps)

        if self.is_work:
            # The task in the work should inherit the same dependency.
            for task in self:
                task.add_deps(deps)

    def remove_deps(self, deps):
        """
        Remove a list of dependencies from the :class:`Node`.

        Args:
            deps: List of :class:`Dependency` objects specifying the  dependencies of the node.
        """
        if not isinstance(deps, (list, tuple)):
            deps = [deps]

        assert all(isinstance(d, Dependency) for d in deps)

        self._deps = [d for d in self._deps if d not in deps]
                                                                                      
        if self.is_work:
            # remove the same list of dependencies from the task in the work
            for task in self:
                task.remove_deps(deps)                                                                                                                                        

    @property
    def deps_status(self):
        """Returns a list with the status of the dependencies."""
        if not self.deps:
            return [self.S_OK]
                                                                  
        return [d.status for d in self.deps]

    def depends_on(self, other):
        """True if this node depends on the other node."""
        return other in [d.node for d in self.deps]

    def get_parents(self):
        """Return the list of nodes in the :class:`Flow` required by this :class:`Node`"""
        parents = []
        for work in self.flow:
            if self.depends_on(work): parents.append(work)
            for task in work:
                if self.depends_on(task): parents.append(task)
        return parents

    def get_children(self):
        """Return the list of nodes in the :class:`Flow` that depends on this :class:`Node`"""
        children = []
        for work in self.flow:
            if work.depends_on(self): children.append(work)
            for task in work:
                if task.depends_on(self): children.append(task)
        return children

    def str_deps(self):
        """Return the string representation of the dependencies of the node."""
        lines = []
        app = lines.append

        app("Dependencies of node %s:" % str(self))
        for i, dep in enumerate(self.deps):
            app("%d) %s, status=%s" % (i, dep.info, str(dep.status)))

        return "\n".join(lines)

    def set_cleanup_exts(self, exts=None):
        """
        Set the list of file extensions that should be removed when the task reaches S_OK.

            Args:
                exts: List of file extensions, if exts is None a default list is provided.
        """
        if exts is None: exts = ["WFK", "SUS", "SCR"]
        self._cleanup_exts = set(exts)

    @property
    def cleanup_exts(self):
        """Set of file extensions to remove."""
        try:
            return self._cleanup_exts
        except AttributeError:
            return set()

    #@abc.abstractmethod
    #def set_status(self, status, msg=None):
    #    """
    #    Set and return the status of the None
    #                                                                                     
    #    Args:
    #        status: Status object or string representation of the status
    #        info_msg: string with human-readable message used in the case of errors (optional)
    #    """

    @abc.abstractproperty
    def status(self):
        """The status of the `Node`."""

    @abc.abstractmethod
    def check_status(self):
        """Check the status of the `Node`."""


class FileNode(Node):
    """
    A Node that consists of a file. May be not yet existing

    Mainly used to connect :class:`Task` objects to external files produced in previous runs.
    """
    def __init__(self, filename):
        super(FileNode, self).__init__()
        self.filepath = os.path.abspath(filename)

        # Directories with input|output|temporary data.
        self.workdir = os.path.dirname(self.filepath)

        self.indir = Directory(self.workdir)
        self.outdir = Directory(self.workdir)
        self.tmpdir = Directory(self.workdir)

    @property
    def products(self):
        return [Product.from_file(self.filepath)]

    def opath_from_ext(self, ext):
        return self.filepath

    @property
    def status(self):
        return self.S_OK if os.path.exists(self.filepath) else self.S_ERROR

    def check_status(self):
        return self.status

    def get_results(self, **kwargs):
        results = super(FileNode, self).get_results(**kwargs)
        #results.register_gridfs_files(filepath=self.filepath)
        return results


class HistoryRecord(object):
    """
    A `HistoryRecord` instance represents an entry in the :class:`NodeHistory`.

    `HistoryRecord` instances are created every time something is logged. 
    They contain all the information pertinent to the event being logged. 
    The main information passed in is in msg and args, which are combined
    using str(msg) % args to create the message field of the record. 
    The record also includes information such as when the record was created,
    the source line where the logging call was made

    .. attribute:: levelno      
    
        Numeric logging level for the message (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    .. attribute:: levelname       

        Text logging level for the message ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")

    .. attribute:: pathname      

        Full pathname of the source file where the logging call was issued (if available)

    .. attribute:: filename 

        Filename portion of pathname

    .. attribute:: module

        Module (name portion of filename)

    .. attribute:: lineno 

        Source line number where the logging call was issued (if available)

    .. attribute:: func_name

        Function name

    .. attribute:: created

        Time when the HistoryRecord was created (time.time() return value)

    .. attribute:: asctime 

        Textual time when the HistoryRecord was created

    .. attribute:: message 
        The result of record.getMessage(), computed just as the record is emitted
    """
    def __init__(self, level, pathname, lineno, msg, args, exc_info, func=None):
        """
        Initialize a logging record with interesting information.
        """
        #
        # The following statement allows passing of a dictionary as a sole
        # argument, so that you can do something like
        #  logging.debug("a %(a)d b %(b)s", {'a':1, 'b':2})
        # Suggested by Stefan Behnel.
        # Note that without the test for args[0], we get a problem because
        # during formatting, we test to see if the arg is present using
        # 'if self.args:'. If the event being logged is e.g. 'Value is %d'
        # and if the passed arg fails 'if self.args:' then no formatting
        # is done. For example, logger.warn('Value is %d', 0) would log
        # 'Value is %d' instead of 'Value is 0'.
        # For the use case of passing a dictionary, this should not be a problem.
        if args and len(args) == 1 and isinstance(args[0], dict) and args[0]:
            args = args[0]
        self.args = args
        self.levelno = level
        self.pathname = pathname
        self.msg = msg

        self.levelname = "FOOBAR" #getLevelName(level)

        try:
            self.filename = os.path.basename(pathname)
            self.module = os.path.splitext(self.filename)[0]
        except (TypeError, ValueError, AttributeError):
            self.filename = pathname
            self.module = "Unknown module"

        self.exc_info = exc_info
        self.exc_text = None      # used to cache the traceback text
        self.lineno = lineno
        self.func_name = func
        self.created = time.time()
        self.asctime = time.asctime()

    def __repr__(self):
        return '<%s, %s, %s, %s,\n"%s">' % (self.__class__.__name__, self.levelno, self.pathname, self.lineno, self.msg)

    def __str__(self):
        return self.get_message(metadata=False)

    def get_message(self,  metadata=False, asctime=True):
        """
        Return the message after merging any user-supplied arguments with the message.

        Args:
            metadata: True if function and module name should be added.
            asctime: True if time string should be added.
        """
        msg = self.msg if is_string(self.msg) else str(self.msg)
        if self.args: msg = msg % self.args

        if asctime:
            msg = "[" + self.asctime + "] " + msg

        # Add metadata
        if metadata:
            msg += "\nCalled in function %s in %s:%s" % (self.func_name, self.pathname, self.lineno)

        return msg


class NodeHistory(collections.deque):
    """Logger-like object"""

    def __str__(self):
        return self.to_string()

    def to_string(self, metadata=False):
        """Returns  a string with the history. Set metadata to True to have info on function and module."""
        return "\n".join(rec.get_message(metadata=metadata) for rec in self)

    def info(self, msg, *args, **kwargs):
        """Log 'msg % args' with the info severity level"""
        self._log("INFO", msg, args, kwargs)

    def warning(self, msg, *args, **kwargs):
        """Log 'msg % args' with the warning severity level"""
        self._log("WARNING", msg, args, kwargs)

    def critical(self, msg, *args, **kwargs):
        """Log 'msg % args' with the critical severity level"""
        self._log("CRITICAL", msg, args, kwargs)

    def find_caller(self):
        """
        find the stack frame of the caller so that we can note the source
        file name, line number and function name.

        See also: 

            http://farmdev.com/src/secrets/framehack/
        """
        # next bit filched from 1.5.2's inspect.py
        def currentframe():
            """Return the frame object for the caller's stack frame."""
            try:
                raise Exception
            except:
                return sys.exc_info()[2].tb_frame.f_back

        # CPython implementation detail: This function should be used for internal and specialized purposes only. 
        # It is not guaranteed to exist in all implementations of Python.
        if hasattr(sys, '_getframe'): currentframe = lambda: sys._getframe(3)
        # done filching

        # _srcfile is used when walking the stack to check when we've got the first caller stack frame.
        if hasattr(sys, 'frozen'): #support for py2exe
            _srcfile = "logging%s__init__%s" % (os.sep, __file__[-4:])
        elif __file__[-4:].lower() in ['.pyc', '.pyo']:
            _srcfile = __file__[:-4] + '.py'
        else:
            _srcfile = __file__

        _srcfile = os.path.normcase(_srcfile)

        f = currentframe()
        # On some versions of IronPython, currentframe() returns None if
        # IronPython isn't run with -X:Frames.
        if f is not None:
            f = f.f_back
        rv = "(unknown file)", 0, "(unknown function)"

        while hasattr(f, "f_code"):
            co = f.f_code
            filename = os.path.normcase(co.co_filename)
            if filename == _srcfile:
                f = f.f_back
                continue
            rv = (co.co_filename, f.f_lineno, co.co_name)
            break

        return rv

    def _log(self, level, msg, args, exc_info=None, extra=None):
        """Low-level logging routine which creates a :class:`HistoryRecord`."""
        #if _srcfile:
        if True:
            # IronPython doesn't track Python frames, so findCaller raises an
            # exception on some versions of IronPython. We trap it here so that
            # IronPython can use logging.
            try:
                fn, lno, func = self.find_caller()
            except ValueError:
                fn, lno, func = "(unknown file)", 0, "(unknown function)"
        else:
            fn, lno, func = "(unknown file)", 0, "(unknown function)"

        if exc_info and not isinstance(exc_info, tuple):
            exc_info = sys.exc_info()

        self.append(HistoryRecord(level, fn, lno, msg, args, exc_info, func=func))


# The code below initializes a counter from a file when the module is imported 
# and save the counter's updated value automatically when the program terminates 
# without relying on the application making an explicit call into this module at termination.

_COUNTER = None
_COUNTER_FILE = os.path.join(os.getenv("HOME"), ".abinit", "abipy", "nodecounter")


def init_counter():
    global _COUNTER

    # Make dir and file if not present.
    if not os.path.exists(os.path.dirname(_COUNTER_FILE)):
        os.makedirs(os.path.dirname(_COUNTER_FILE))

    if not os.path.exists(_COUNTER_FILE):
        with open(_COUNTER_FILE, "wt") as fh:
            fh.write("%d\n" % -1)

    if _COUNTER is None:
        with open(_COUNTER_FILE, "r") as fh:
            _COUNTER = int(fh.read())


def get_newnode_id():
    """
    Returns a new node identifier used for :class:`Task`, :class:`Work` and :class:`Flow` objects.

    .. warning:

        The id is unique inside the same python process so be careful when 
        Works and Tasks are constructed at run-time or when threads are used.
    """
    init_counter()

    global _COUNTER
    _COUNTER += 1
    return _COUNTER


def save_lastnode_id():
    """Save the id of the last node created."""
    init_counter()

    with FileLock(_COUNTER_FILE):
        with AtomicFile(_COUNTER_FILE, mode="w") as fh:
            fh.write("%d\n" % _COUNTER)


# Register function atexit
import atexit
atexit.register(save_lastnode_id)

