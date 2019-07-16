# coding: utf-8
"""
This module defines the Node class that is inherited by Task, Work and Flow objects.
"""

import sys
import os
import time
import collections
import abc
import numpy as np

from pprint import pprint
from pymatgen.util.io_utils import AtomicFile
from pydispatch import dispatcher
from monty.termcolor import colored
from monty.serialization import loadfn
from monty.string import is_string
from monty.io import FileLock
from monty.collections import AttrDict, Namespace
from monty.functools import lazy_property
from monty.json import MSONable
from pymatgen.util.serialization import json_pretty_dump, pmg_serialize
from .utils import File, Directory, Dirviz, irdvars_for_ext, abi_extensions


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
        (2,  "Locked",        "grey"   , None, None),         # Task is locked an must be explicitly unlocked by an external subject (Work).
        (3,  "Ready",         None     , None, None),         # Node is ready i.e. all the depencies of the node have status S_OK
        (4,  "Submitted",     "blue"   , None, None),         # Node has been submitted (The `Task` is running or we have started to finalize the Work)
        (5,  "Running",       "magenta", None, None),         # Node is running.
        (6,  "Done",          None     , None, None),         # Node done, This does not imply that results are ok or that the calculation completed successfully
        (7,  "AbiCritical",   "red"    , None, None),         # Node raised an Error by ABINIT.
        (8,  "QCritical",     "red"    , "on_white", None),   # Node raised an Error by submitting submission script, or by executing it
        (9,  "Unconverged",   "red"    , "on_yellow", None),  # This usually means that an iterative algorithm didn't converge.
        (10, "Error",         "red"    , None, None),         # Node raised an unrecoverable error, usually raised when an attempt to fix one of other types failed.
        (11, "Completed",     "green"  , None, None),         # Execution completed successfully.
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
        if obj is None: return None
        return obj if isinstance(obj, cls) else cls.from_string(obj)

    @classmethod
    def from_string(cls, s):
        """Return a `Status` instance from its string representation."""
        for num, text in cls._STATUS2STR.items():
            if text == s:
                return cls(num)
        else:
            raise ValueError("Wrong string %s" % s)

    @classmethod
    def all_status_strings(cls):
        """List of strings with all possible values status."""
        return [info[1] for info in cls._STATUS_INFO]

    @property
    def is_critical(self):
        """True if status is critical."""
        return str(self) in ("AbiCritical", "QCritical", "Unconverged", "Error")

    @property
    def color_opts(self):
        return self._STATUS2COLOR_OPTS[self]

    @property
    def colored(self):
        """Return colorized text used to print the status if the stream supports it."""
        return colored(str(self), **self.color_opts)


class Dependency:
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

        if exts and is_string(exts): exts = exts.split()

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
        """
        This function is called when we specify the task dependencies with the syntax:

            deps={node: "@property"}

        In this case the task has to the get `property` from `node` before starting the calculation.

        At present, the following properties are supported:

            - @structure
        """
        if not self.getters: return

        for getter in self.getters:
            if getter == "@structure":
                task.history.info("Getting structure from %s" % self.node)
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


class Product:
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
        super().__init__(path=path, fs_id=fs_id, mode=mode)


class NodeResults(dict, MSONable):
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
        super().__init__(**kwargs)
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


def check_spectator(node_method):
    """
    Decorator for :class:`Node` methods. Raise `SpectatorNodeError`.
    """
    from functools import wraps
    @wraps(node_method)
    def wrapper(*args, **kwargs):
        node = args[0]
        if node.in_spectator_mode:
            #raise node.SpectatorError("You should not call this method when the node in spectator_mode")
            #warnings.warn("You should not call %s when the node in spectator_mode" % node_method)
            import warnings

        return node_method(*args, **kwargs)

    return wrapper


class NodeError(Exception):
    """Base Exception raised by :class:`Node` subclasses"""


class SpectatorNodeError(NodeError):
    """
    Exception raised by :class:`Node` methods when the node is in spectator mode
    and we are calling a method with side effects.
    """


class Node(metaclass=abc.ABCMeta):
    """
    Abstract base class defining the interface that must be
    implemented by the nodes of the calculation.

    Nodes are hashable and can be tested for equality
    """
    Results = NodeResults

    Error  = NodeError
    SpectatorError = SpectatorNodeError

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
    #S_CANCELLED = Status.from_string("Cancelled")
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
        #S_CANCELLED,
        S_ERROR,
        S_OK,
    ]

    # Color used to plot the network in networkx
    color_rgb = np.array((105, 105, 105)) / 255

    def __init__(self):
        self._in_spectator_mode = False
        # Node identifier.
        self._node_id = get_newnode_id()

        # List of dependencies
        self._deps = []

        # List of files (products) needed by this node.
        self._required_files = []

        # Used to push additional info during the execution.
        self.history = NodeHistory(maxlen=80)

        # Actions performed to fix abicritical events.
        self._corrections = NodeCorrections()

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
                self.__class__.__name__, self.node_id, self.relworkdir)
        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, node_id=%s, workdir=None>" % (self.__class__.__name__, self.node_id)

    #def __setattr__(self, name, value):
    #    if self.in_spectator_mode:
    #        raise RuntimeError("You should not call __setattr__ in spectator_mode")
    #    return super().__setattr__(name,value)

    @lazy_property
    def color_hex(self):
        """Node color as Hex Triplet https://en.wikipedia.org/wiki/Web_colors#Hex_triplet"""
        def clamp(x):
            return max(0, min(int(x), 255))

        r, g, b = np.trunc(self.color_rgb * 255)
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

    def isinstance(self, class_or_string):
        """
        Check whether the node is a instance of `class_or_string`.
        Unlinke the standard isinstance builtin, the method accepts either a class or a string.
        In the later case, the string is compared with self.__class__.__name__ (case insensitive).
        """
        if class_or_string is None:
            return False
        import inspect
        if inspect.isclass(class_or_string):
            return isinstance(self, class_or_string)
        else:
            return self.__class__.__name__.lower() == class_or_string.lower()

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
            if self.is_task:
                try:
                    return self.pos_str
                except:
                    return os.path.basename(self.workdir)
            else:
                return os.path.basename(self.workdir)

    @property
    def relworkdir(self):
        """Return a relative version of the workdir"""
        if getattr(self, "workdir", None) is None:
            return None
        try:
            return os.path.relpath(self.workdir)
        except OSError:
            # current working directory may not be defined!
            return self.workdir

    def set_name(self, name):
        """Set the name of the Node."""
        self._name = name

    @property
    def node_id(self):
        """Node identifier."""
        return self._node_id

    @check_spectator
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
        self.history.info("Finalized set to %s" % self._finalized)

    @property
    def in_spectator_mode(self):
        return self._in_spectator_mode

    @in_spectator_mode.setter
    def in_spectator_mode(self, mode):
        self._in_spectator_mode = bool(mode)
        #self.history.info("in_spectator_mode set to %s" % mode)

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
            event: :class:`AbinitEvent` that triggered the correction.
            action (str): Human-readable string with info on the action perfomed to solve the problem.
        """
        # TODO: Create CorrectionObject
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

    @check_spectator
    def add_deps(self, deps):
        """
        Add a list of dependencies to the :class:`Node`.

        Args:
            deps: List of :class:`Dependency` objects specifying the dependencies of the node.
                  or dictionary mapping nodes to file extensions e.g. {task: "DEN"}
        """
        if isinstance(deps, collections.abc.Mapping):
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

        # If we have a FileNode as dependency, add self to its children
        # Node.get_parents will use this list if node.is_isfile.
        for dep in (d for d in deps if d.node.is_file):
            dep.node.add_filechild(self)

    @check_spectator
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
        return [d.node for d in self.deps]
        #parents = []
        #for work in self.flow:
        #    if self.depends_on(work): parents.append(work)
        #    for task in work:
        #        if self.depends_on(task): parents.append(task)
        #return parents

    def get_children(self):
        """
        Return the list of nodes in the :class:`Flow` that depends on this :class:`Node`

        .. note::

            This routine assumes the entire flow has been allocated.
        """
        # Specialized branch for FileNode.
        if self.is_file:
            return self.filechildren

        # Inspect the entire flow to get children.
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

    def get_vars_dataframe(self, *varnames):
        """
        Return pandas DataFrame with the value of the variables specified in `varnames`.
        Can be used for task/works/flow. It's recursive!

        .. example:

            flow.get_vars_dataframe("ecut", "ngkpt")
            work.get_vars_dataframe("acell", "usepawu")
        """
        import pandas as pd
        if self.is_task:
            df = pd.DataFrame([{v: self.input.get(v, None) for v in varnames}], index=[self.name], columns=varnames)
            df["class"] = self.__class__.__name__
            return df

        elif self.is_work:
            frames = [task.get_vars_dataframe(*varnames) for task in self]
            return pd.concat(frames)

        elif self.is_flow:
            frames = [work.get_vars_dataframe(*varnames) for work in self]
            return pd.concat(frames)

        else:
            #print("Ignoring node of type: `%s`" % type(self))
            return pd.DataFrame(index=[self.name])

    def get_graphviz_dirtree(self, engine="automatic", **kwargs):
        """
        Generate directory graph in the DOT language. The graph show the files and directories
        in the node workdir.

        Returns: graphviz.Digraph <https://graphviz.readthedocs.io/en/stable/api.html#digraph>
        """
        if engine == "automatic":
            engine = "fdp"

        return Dirviz(self.workdir).get_cluster_graph(engine=engine, **kwargs)

    def set_gc(self, gc):
        """
        Set the garbage collector.
        """
        assert isinstance(gc, GarbageCollector)
        self._gc = gc

    @property
    def gc(self):
        """
        Garbage collector. None if garbage collection is deactivated.
        Use flow.set_garbage_collector to initialize the object.
        """
        try:
            return self._gc
        except AttributeError:
            #if not self.is_flow and self.flow.gc: return self.flow.gc
            return None

    @property
    def event_handlers(self):
        """
        The list of handlers registered for this node.
        If the node is not a `Flow` and does not have its own list of
        `handlers` the handlers registered at the level of the flow are returned.

        This trick allows one to registered different handlers at the level of the Task
        for testing purposes. By default, we have a common list of handlers for all the nodes in the flow.
        This choice facilitates the automatic installation of the handlers when we use callbacks to generate
        new Works and Tasks!
        """
        if self.is_flow:
            return self._event_handlers

        try:
            return self._event_handlers
        except AttributeError:
            return self.flow._event_handlers

    @check_spectator
    def install_event_handlers(self, categories=None, handlers=None):
        """
        Install the `EventHandlers for this `Node`. If no argument is provided
        the default list of handlers is installed.

        Args:
            categories: List of categories to install e.g. base + can_change_physics
            handlers: explicit list of :class:`EventHandler` instances.
                      This is the most flexible way to install handlers.

        .. note::

            categories and handlers are mutually exclusive.
        """
        if categories is not None and handlers is not None:
            raise ValueError("categories and handlers are mutually exclusive!")

        from .events import get_event_handler_classes
        if categories:
            raise NotImplementedError()
            handlers = [cls() for cls in get_event_handler_classes(categories=categories)]
        else:
            handlers = handlers or [cls() for cls in get_event_handler_classes()]

        self._event_handlers = handlers

    def show_event_handlers(self, stream=sys.stdout, verbose=0):
        """Print to `stream` the event handlers installed for this flow."""
        lines = ["List of event handlers installed:"]
        for handler in self.event_handlers:
            if verbose:
                lines.extend(handler.__class__.cls2str().split("\n"))
            else:
                lines.extend(str(handler).split("\n"))

        stream.write("\n".join(lines))
        stream.write("\n")

    def send_signal(self, signal):
        """
        Send signal from this node to all connected receivers unless the node is in spectator mode.

        signal -- (hashable) signal value, see `dispatcher` connect for details

        Return a list of tuple pairs [(receiver, response), ... ]
        or None if the node is in spectator mode.

        if any receiver raises an error, the error propagates back
        through send, terminating the dispatch loop, so it is quite
        possible to not have all receivers called if a raises an error.
        """
        if self.in_spectator_mode: return None
        logger.debug("Node %s broadcasts signal %s" % (self, signal))
        dispatcher.send(signal=signal, sender=self)

   ##########################
   ### Abstract protocol ####
   ##########################

    @property
    @abc.abstractmethod
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
    color_rgb = np.array((102, 51, 255)) / 255

    def __init__(self, filename):
        super().__init__()
        self.filepath = os.path.abspath(filename)

        # Directories with input|output|temporary data.
        self.workdir = os.path.dirname(self.filepath)

        self.indir = Directory(self.workdir)
        self.outdir = Directory(self.workdir)
        self.tmpdir = Directory(self.workdir)

        self._filechildren = []

    def __repr__(self):
        try:
            return "<%s, node_id=%s, rpath=%s>" % (
                self.__class__.__name__, self.node_id, os.path.relpath(self.filepath))
        except AttributeError:
            # this usually happens when workdir has not been initialized
            return "<%s, node_id=%s, path=%s>" % (self.__class__.__name__, self.node_id, self.filepath)

    @lazy_property
    def basename(self):
        """Basename of the file."""
        return os.path.basename(self.filepath)

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
        results = super().get_results(**kwargs)
        #results.register_gridfs_files(filepath=self.filepath)
        return results

    def add_filechild(self, node):
        """Add a node (usually Task) to the children of this FileNode."""
        self._filechildren.append(node)

    @property
    def filechildren(self):
        """List with the children (nodes) of this FileNode."""
        return self._filechildren

    # This part provides IO capabilities to FileNode with API similar to the one implemented in Task.
    # We may need it at runtime to extract information from netcdf files e.g.
    # a NscfTask will change the FFT grid to match the one used in the GsTask.

    def abiopen(self):
        from abipy import abilab
        return abilab.abiopen(self.filepath)

    def open_gsr(self):
        return self._abiopen_abiext("_GSR.nc")

    def _abiopen_abiext(self, abiext):
        import glob
        from abipy import abilab
        if not self.filepath.endswith(abiext):
            msg = """\n
File type does not match the abinit file extension.
Caller asked for abiext: `%s` whereas filepath: `%s`.
Continuing anyway assuming that the netcdf file provides the API/dims/vars neeeded by the caller.
""" % (abiext, self.filepath)
            logger.warning(msg)
            self.history.warning(msg)

        #try to find file in the same path
        filepath = os.path.dirname(self.filepath)
        glob_result = glob.glob(os.path.join(filepath,"*%s"%abiext))
        if len(glob_result): return abilab.abiopen(glob_result[0])
        return self.abiopen()


class HistoryRecord:
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
        # Remove milliseconds
        i = self.asctime.find(".")
        if i != -1: self.asctime = self.asctime[:i]

    def __repr__(self):
        return '<%s, %s, %s, %s,\n"%s">' % (self.__class__.__name__, self.levelno, self.pathname, self.lineno, self.msg)

    def __str__(self):
        return self.get_message(metadata=False)

    def get_message(self, metadata=False, asctime=True):
        """
        Return the message after merging any user-supplied arguments with the message.

        Args:
            metadata: True if function and module name should be added.
            asctime: True if time string should be added.
        """
        msg = self.msg if is_string(self.msg) else str(self.msg)
        if self.args:
            try:
                msg = msg % self.args
            except:
                msg += str(self.args)

        if asctime: msg = "[" + self.asctime + "] " + msg

        # Add metadata
        if metadata:
            msg += "\nCalled by %s at %s:%s\n" % (self.func_name, self.pathname, self.lineno)

        return msg

    @pmg_serialize
    def as_dict(self):
        return {'level': self.levelno, 'pathname': self.pathname, 'lineno': self.lineno, 'msg': self.msg,
                'args': self.args, 'exc_info': self.exc_info, 'func': self.func_name}

    @classmethod
    def from_dict(cls, d):
        return cls(level=d['level'], pathname=d['pathname'], lineno=int(d['lineno']), msg=d['msg'], args=d['args'],
                   exc_info=d['exc_info'], func=d['func'])


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

    def _log(self, level, msg, args, exc_info=None, extra=None):
        """Low-level logging routine which creates a :class:`HistoryRecord`."""
        if exc_info and not isinstance(exc_info, tuple):
            exc_info = sys.exc_info()

        self.append(HistoryRecord(level, "unknown filename", 0, msg, args, exc_info, func="unknown func"))


class NodeCorrections(list):
    """Iterable storing the correctios performed by the :class:`EventHandler`"""
    #TODO
    # Correction should have a human-readable message
    # and a list of operatins in JSON format (Modder?) so that
    # we can read them and re-apply the corrections to another task if needed.

    #def count_event_class(self, event_class):
    #    """
    #    Return the number of times the event class has been already fixed.
    #    """
    #    #return len([c for c in self if c["event"]["@class"] == str(event_class)])

    #def _find(self, event_class)


class GarbageCollector:
    """This object stores information on the """
    def __init__(self, exts, policy):
        self.exts, self.policy = set(exts), policy


# The code below initializes a counter from a file when the module is imported
# and save the counter's updated value automatically when the program terminates
# without relying on the application making an explicit call into this module at termination.

_COUNTER = None
_COUNTER_FILE = os.path.join(os.path.expanduser("~"), ".abinit", "abipy", "nodecounter")


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
            s = fh.read().strip()
            if not s: s = "-1"
            _COUNTER = int(s)


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

