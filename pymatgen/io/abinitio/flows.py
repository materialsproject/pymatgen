# coding: utf-8
"""
Abinit Flows
"""
from __future__ import unicode_literals, division, print_function

import os
import sys
import time
import collections
import warnings
import shutil
import pickle
import copy

from six.moves import map 
from atomicfile import AtomicFile
from pprint import pprint
from prettytable import PrettyTable
from monty.io import FileLock
from monty.termcolor import cprint, colored, stream_has_colours
from pymatgen.serializers.pickle_coders import pmg_pickle_load, pmg_pickle_dump 
from .tasks import Dependency, Status, Node, NodeResults, Task, ScfTask, PhononTask, TaskManager, NscfTask
#from .tasks imort AnaddbTask, QpMergeTask
from .utils import Directory, Editor
from .abiinspect import yaml_read_irred_perts
from .workflows import Workflow, BandStructureWorkflow, PhononWorkflow, G0W0_Workflow, QptdmWorkflow

try:
    from pydispatch import dispatcher
except ImportError:
    pass

import logging
logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"


__all__ = [
    "AbinitFlow",
    "G0W0WithQptdmFlow",
    "bandstructure_flow",
    "g0w0_flow",
    "phonon_flow",
]


class FlowResults(NodeResults):

    JSON_SCHEMA = NodeResults.JSON_SCHEMA.copy() 
    #JSON_SCHEMA["properties"] = {
    #    "queries": {"type": "string", "required": True},
    #}

    @classmethod
    def from_node(cls, flow):
        """Initialize an instance from a WorkFlow instance."""
        new = super(FlowResults, cls).from_node(flow)

        #new.update(
        #    #input=flow.strategy
        #)

        # Will put all files found in outdir in GridFs 
        d = {os.path.basename(f): f for f in flow.outdir.list_filepaths()}

        # Add the pickle file.
        pickle_path = os.path.join(flow.workdir, flow.PICKLE_FNAME)
        d["pickle"] = pickle_path if flow.pickle_protocol != 0 else (pickle_path, "t")
        new.add_gridfs_files(**d)

        return new


class AbinitFlow(Node):
    """
    This object is a container of workflows. Its main task is managing the 
    possible inter-depedencies among the workflows and the creation of
    dynamic worfflows that are generates by callbacks registered by the user.

    .. attributes:

        creation_date:
            String with the creation_date

        pickle_protocol: 
            Protocol for Pickle database (default: -1 i.e. latest protocol)
    """
    VERSION = "0.1"
    PICKLE_FNAME = "__AbinitFlow__.pickle"

    Results = FlowResults

    def __init__(self, workdir, manager=None, pickle_protocol=-1):
        """
        Args:
            workdir:
                String specifying the directory where the workflows will be produced.
            manager:
                `TaskManager` object responsible for the submission of the jobs.
                If manager is None, the object is initialized from the yaml file
                located either in the working directory or in the user configuration dir.
            pickle_procol:
                Pickle protocol version used for saving the status of the object.
                -1 denotes the latest version supported by the python interpreter.
        """
        super(AbinitFlow, self).__init__()

        self.set_workdir(workdir)

        self.creation_date = time.asctime()

        if manager is None: manager = TaskManager.from_user_config()
        self.manager = manager.deepcopy()

        # List of workflows.
        self._works = []

        self._waited = 0

        # List of callbacks that must be executed when the dependencies reach S_OK
        self._callbacks = []

        self.pickle_protocol = int(pickle_protocol)

        # ID used to access mongodb
        self._mongo_id = None

        # TODO
        # Signal slots: a dictionary with the list 
        # of callbacks indexed by node_id and SIGNAL_TYPE.
        # When the node changes its status, it broadcast a signal.
        # The flow is listening to all the nodes of the calculation
        # [node_id][SIGNAL] = list_of_signal_handlers
        #self._sig_slots =  slots = {}
        #for work in self:
        #    slots[work] = {s: [] for s in work.S_ALL}

        #for task in self.iflat_tasks():
        #    slots[task] = {s: [] for s in work.S_ALL}

    # This is needed for fireworks although Node.__str__ and __repr__ are much more readable.
    #def __repr__(self):
    #    return self.workdir

    #def __str__(self):
    #    return repr(self)

    def as_dict(self, **kwargs):
        """
        JSON serialization, note that we only need to save 
        a string with the working directory since the object will be 
        reconstructed from the pickle file located in workdir
        """
        return {"workdir": self.workdir}

    # This is needed for fireworks.
    to_dict = as_dict

    @classmethod
    def from_dict(cls, d, **kwargs):
        """Reconstruct the flow from the pickle file."""
        return cls.pickle_load(d["workdir"], **kwargs)

    def set_workdir(self, workdir, chroot=False):
        """
        Set the working directory. Cannot be set more than once unless chroot is True
        """
        if not chroot and hasattr(self, "workdir") and self.workdir != workdir:
            raise ValueError("self.workdir != workdir: %s, %s" % (self.workdir,  workdir))

        # Directories with (input|output|temporary) data.
        self.workdir = os.path.abspath(workdir)
        self.indir = Directory(os.path.join(self.workdir, "indata"))
        self.outdir = Directory(os.path.join(self.workdir, "outdata"))
        self.tmpdir = Directory(os.path.join(self.workdir, "tmpdata"))

    @classmethod
    def pickle_load(cls, filepath, disable_signals=False):
        """
        Loads the object from a pickle file and performs initial setup.

        Args:
            filepath:
                Filename or directory name. It filepath is a directory, we 
                scan the directory tree starting from filepath and we 
                read the first pickle database. Raise RuntimeError if multiple
                databases are found.
            disable_signals:
                If True, the nodes of the flow are not connected by signals.
                This option is usually used when we want to read a flow 
                in read-only mode and we want to avoid any possible side effect.
        """
        if os.path.isdir(filepath):
            # Walk through each directory inside path and find the pickle database.
            for dirpath, dirnames, filenames in os.walk(filepath):
                fnames = [f for f in filenames if f == cls.PICKLE_FNAME]
                if fnames:
                    if len(fnames) == 1:
                        filepath = os.path.join(dirpath, fnames[0])
                        break  # Exit os.walk
                    else:
                        err_msg = "Found multiple databases:\n %s" % str(fnames)
                        raise RuntimeError(err_msg)
            else:
                err_msg = "Cannot find %s inside directory %s" % (cls.PICKLE_FNAME, filepath)
                raise ValueError(err_msg)

        with FileLock(filepath):
            with open(filepath, "rb") as fh:
                flow = pmg_pickle_load(fh)

        # Check if versions match.
        if flow.VERSION != cls.VERSION:
            msg = ("File flow version %s != latest version %s\n."
                   "Regenerate the flow to solve the problem " % (flow.VERSION, cls.VERSION))
            warnings.warn(msg)

        if not disable_signals:
            flow.connect_signals()

        # Recompute the status of each task since tasks that
        # have been submitted previously might be completed.
        flow.check_status()
        return flow

    def __len__(self):
        return len(self.works)

    def __iter__(self):
        return self.works.__iter__()

    def __getitem__(self, slice):
        return self.works[slice]

    @property
    def mongo_id(self):
        return self._mongo_id

    @mongo_id.setter
    def mongo_id(self, value):
        if self.mongo_id is not None:
            raise RuntimeError("Cannot change mongo_id %s" % self.mongo_id)
        self._mongo_id = value

    def validate_json_schema(self):
        """Validate the JSON schema. Return list of errors."""
        errors = []

        for work in self:
            for task in work:
                if not task.get_results().validate_json_schema(): 
                    errors.append(task)
            if not work.get_results().validate_json_schema(): 
                errors.append(work)
        if not self.get_results().validate_json_schema(): 
            errors.append(self)

        return errors

    @property
    def works(self):
        """List of `Workflow` objects contained in self.."""
        return self._works

    @property
    def all_ok(self):
        """True if all the tasks in workflows have reached S_OK."""
        return all(work.all_ok for work in self)

    @property
    def all_tasks(self):
        return self.iflat_tasks()

    @property
    def num_tasks(self):
        """Total number of tasks"""
        return len(list(self.iflat_tasks()))

    @property
    def errored_tasks(self):
        """List of errored tasks."""
        return list(self.iflat_tasks(status=self.S_ERROR))

    @property
    def num_errored_tasks(self):
        """The number of tasks whose status is `S_ERROR`."""
        return len(self.errored_tasks)

    @property
    def unconverged_tasks(self):
        """List of unconverged tasks."""
        return list(self.iflat_tasks(status=self.S_UNCONVERGED))

    @property
    def num_unconverged_tasks(self):
        """The number of tasks whose status is `S_UNCONVERGED`."""
        return len(self.unconverged_tasks)

    @property
    def status_counter(self):
        """
        Returns a `Counter` object that counts the number of tasks with 
        given status (use the string representation of the status as key).
        """
        # Count the number of tasks with given status in each workflow.
        counter = self[0].status_counter
        for work in self[1:]:
            counter += work.status_counter

        return counter

    @property
    def ncores_reserved(self):
        """
        Returns the number of cores reserved in this moment.
        A core is reserved if the task is not running but
        we have submitted the task to the queue manager.
        """
        return sum(work.ncores_reserved for work in self)

    @property
    def ncores_allocated(self):
        """
        Returns the number of cores allocated in this moment.
        A core is allocated if it's running a task or if we have
        submitted a task to the queue manager but the job is still pending.
        """
        return sum(work.ncores_allocated for work in self)

    @property
    def ncores_inuse(self):
        """
        Returns the number of cores used in this moment.
        A core is used if there's a job that is running on it.
        """
        return sum(work.ncores_inuse for work in self)

    @property
    def has_chrooted(self):
        """
        Returns a string that evaluates to True if we have changed 
        the workdir for visualization purposes e.g. we are using sshfs.
        to mount the remote directory where the `Flow` is located.
        The string gives the previous workdir of the flow.
        """
        try:
            return self._chrooted_from
        except AttributeError:
            return ""

    def chroot(self, new_workdir):
        """
        Change the workir of the `Flow`. Mainly used for
        allowing the user to open the GUI on the local host
        and access the flow from remote via sshfs.

        .. note:
            Calling this method will make the flow go in read-only mode.
        """
        self._chrooted_from = self.workdir
        self.set_workdir(new_workdir, chroot=True)

        for i, work in enumerate(self):
            new_wdir = os.path.join(self.workdir, "w" + str(i))
            work.chroot(new_wdir)

    def groupby_status(self):
        """
        Returns a ordered dictionary mapping the task status to 
        the list of named tuples (task, work_index, task_index).
        """
        Entry = collections.namedtuple("Entry", "task wi ti")
        d = collections.defaultdict(list)

        for task, wi, ti in self.iflat_tasks_wti():
            d[task.status].append(Entry(task, wi, ti))

        # Sort keys according to their status.
        return collections.OrderedDict([(k, d[k]) for k in sorted(list(d.keys()))])

    def iflat_tasks_wti(self, status=None, op="=="):
        """
        Generator to iterate over all the tasks of the `Flow`.
        Yields

            (task, work_index, task_index)

        If status is not None, only the tasks whose status satisfies
        the condition (task.status op status) are selected
        status can be either one of the flags defined in the `Task` class 
        (e.g Task.S_OK) or a string e.g "S_OK" 
        """
        return self._iflat_tasks_wti(status=status, op=op, with_wti=True)

    def iflat_tasks(self, status=None, op="=="):
        """
        Generator to iterate over all the tasks of the `Flow`.

        If status is not None, only the tasks whose status satisfies
        the condition (task.status op status) are selected
        status can be either one of the flags defined in the `Task` class 
        (e.g Task.S_OK) or a string e.g "S_OK" 
        """
        return self._iflat_tasks_wti(status=status, op=op, with_wti=False)

    def _iflat_tasks_wti(self, status=None, op="==", with_wti=True):
        """
        Generators that produces a flat sequence of task.
        if status is not None, only the tasks with the specified status are selected.

        Returns:
            (task, work_index, task_index) if with_wti is True else task
        """
        if status is None:
            for wi, work in enumerate(self):
                for ti, task in enumerate(work):
                    if with_wti:
                        yield task, wi, ti
                    else:
                        yield task

        else:
            # Get the operator from the string.
            import operator
            op = {
                "==": operator.eq,
                "!=": operator.ne,
                ">": operator.gt,
                ">=": operator.ge,
                "<": operator.lt,
                "<=": operator.le,
            }[op]

            # Accept Task.S_FLAG or string.
            status = Status.as_status(status)

            for wi, work in enumerate(self):
                for ti, task in enumerate(work):
                    if op(task.status, status):
                        if with_wti:
                            yield task, wi, ti
                        else:
                            yield task

    def check_dependencies(self):
        """Test the dependencies of the nodes for possible deadlocks."""
        deadlocks = []

        for task in self.all_tasks:
            for dep in task.deps:
                if dep.node.depends_on(task):
                    deadlocks.append((task, dep.node))

        if deadlocks:
            lines = ["Detect wrong list of dependecies that will lead to a deadlock:"]
            lines.extend(["%s <--> %s" % nodes for nodes in deadlocks])
            raise ValueError("\n".join(lines))

    #def detect_deadlock(self):
    #    eu_tasks = list(self.errored_tasks) + list(self.unconverged_tasks)
    #     if not eu_tasks:
    #        return []

    #    deadlocked = []
    #    for task in self.all_tasks:
    #        if any(task.depends_on(eu_task) for eu_task in eu_tasks):
    #            deadlocked.append(task)

    #    return deadlocked

    def check_status(self):
        """Check the status of the workflows in self."""
        for work in self:
            work.check_status()

    #def set_status(self, status):

    @property
    def status(self):
        """The status of the flow i.e. the minimum of the status of its tasks and its works"""
        return min(work.get_all_status(only_min=True) for work in self)

    def fix_critical(self):
        self.fix_queue_critical()
        self.fix_abi_critical()

    def fix_abi_critical(self):
        """
        Fixer for critical events originating form abinit
        """
        for task in self.iflat_tasks(status=Task.S_ABICRITICAL):
            #todo
            if task.fix_abicritical():
                task.reset_from_scratch()
                # task.set_status(Task.S_READY)
            else:
                info_msg = 'We encountered an abi critial envent that could not be fixed'
                logger.warning(info_msg)
                task.set_status(status=task.S_ERROR)

    def fix_queue_critical(self):
        """
        Fixer for errors originating from the scheduler.

        General strategy, first try to increase resources in order to fix the problem,
        if this is not possible, call a task specific method to attempt to decrease the demands.
        """
        from pymatgen.io.gwwrapper.scheduler_error_parsers import NodeFailureError, MemoryCancelError, TimeCancelError

        for task in self.iflat_tasks(status=Task.S_QUEUECRITICAL):
            logger.info("Will try to fix task %s" % str(task))

            if not task.queue_errors:
                # queue error but no errors detected, try to solve by increasing resources
                # if resources are at maximum the tast is definitively turned to errored
                if self.manager.increase_resources():  # acts either on the policy or on the qadapter
                    task.reset_from_scratch()
                    return True
                else:
                    info_msg = 'unknown queue error, could not increase resources any further'
                    task.set_status(task.S_ERROR, info_msg)
                    return False
            else:
                for error in task.queue_errors:
                    logger.info('fixing : %s' % str(error))
                    if isinstance(error, NodeFailureError):
                        # if the problematic node is know exclude it
                        if error.nodes is not None:
                            task.manager.qadapter.exclude_nodes(error.nodes)
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='increased resources')
                        else:
                            info_msg = 'Node error detected but no was node identified. Unrecoverable error.'
                            return task.set_status(task.S_ERROR, info_msg)
                    elif isinstance(error, MemoryCancelError):
                        # ask the qadapter to provide more resources, i.e. more cpu's so more total memory
                        if task.manager.increase_resources():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='increased mem')
                        # if the max is reached, try to increase the memory per cpu:
                        elif task.manager.qadapter.increase_mem():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='increased mem')
                        # if this failed ask the task to provide a method to reduce the memory demand
                        elif task.reduce_memory_demand():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='decreased mem demand')
                        else:
                            info_msg = 'Memory error detected but the memory could not be increased neigther could the ' \
                                       'memory demand be decreased. Unrecoverable error.'
                            return task.set_status(task.S_ERROR, info_msg)
                    elif isinstance(error, TimeCancelError):
                        # ask the qadapter to provide more memory
                        if task.manager.qadapter.increase_time():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='increased wall time')
                        # if this fails ask the qadapter to increase the number of cpus
                        elif task.manager.increase_resources():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='increased number of cpus')
                        # if this failed ask the task to provide a method to speed up the task
                        elif task.speed_up():
                            task.reset_from_scratch()
                            return task.set_status(task.S_READY, info_msg='task speedup')
                        else:
                            info_msg = 'Time cancel error detected but the time could not be increased neigther could ' \
                                       'the time demand be decreased by speedup of increasing the number of cpus. ' \
                                       'Unrecoverable error.'
                            return task.set_status(task.S_ERROR, info_msg)
                    else:
                        info_msg = 'No solution provided for error %s. Unrecoverable error.' % error.name
                        logger.debug(info_msg)
                        return task.set_status(task.S_ERROR, info_msg)

    def show_status(self, stream=sys.stdout, verbose=0):
        """
        Report the status of the workflows and the status 
        of the different tasks on the specified stream.

        if not verbose, no full entry for works that are completed is printed.
        """
        has_colours = stream_has_colours(stream)
        red = "red" if has_colours else None

        for i, work in enumerate(self):
            print(80*"=", file=stream)
            print("Workflow #%d: %s, Finalized=%s\n" % (i, work, work.finalized), file=stream)

            if verbose == 0 and work.finalized:
                continue

            table =PrettyTable([
                "Task", "Status", "Queue-id", "Errors", "Warnings", "Comments", 
                "MPI", "OMP", "Restarts", "Task-Class", "Run-Etime"])

            tot_num_errors = 0
            for task in work:
                task_name = os.path.basename(task.name)

                # Parse the events in the main output.
                report = task.get_event_report()

                events = map(str, 3*["N/A"])
                if report is not None: 
                    events = map(str, [report.num_errors, report.num_warnings, report.num_comments])
                events = list(events)

                cpu_info = list(map(str, [task.mpi_procs, task.omp_threads]))
                task_info = list(map(str, [task.num_restarts, task.__class__.__name__, task.run_etime()]))

                if task.status.is_critical:
                    tot_num_errors += 1
                    task_name = colored(task_name, red)

                if has_colours:
                    table.add_row([task_name, task.status.colored, str(task.queue_id)] +  events + cpu_info + task_info)
                else:
                    table.add_row([task_name, str(task.status), str(task.queue_id)] +  events + cpu_info + task_info)

            # Print table and write colorized line with the total number of errors.
            print(table, file=stream)
            if tot_num_errors:
                cprint("Total number of errors: %d" % tot_num_errors, red, file=stream)

    def get_results(self, **kwargs):
        results = self.Results.from_node(self)
        results.update(self.get_dict_for_mongodb_queries())
        return results

    def get_dict_for_mongodb_queries(self):
        """
        This function returns a dictionary with the attributes that will be 
        put in the mongodb document to facilitate the query. 
        Subclasses may want to replace or extend the default behaviour.
        """
        d = {}
        return d
        # TODO
        all_structures = [task.strategy.structure for task in self.iflat_tasks()]
        all_pseudos = [task.strategy.pseudos for task in self.iflat_tasks()]

    def look_before_you_leap(self):
        """
        This method should be called before running the calculation to make
        sure that the most important requirements are satisfied.
        """
        errors = []
        if self.has_db:
            try:
                self.manager.db_connector.get_collection()
            except Exception as exc:
                errors.append("""
                    ERROR while trying to connect to the MongoDB database:
                        Exception:
                            %s
                        Connector:
                            %s
                    """ % (exc, self.manager.db_connector))

        return "\n".join(errors)

    @property
    def has_db(self):
        """True if flow uses MongoDB to store the results."""
        return self.manager.has_db

    def db_insert(self):
        """
        Insert results in the mongdob database.
        """
        assert self.has_db
        # Connect to MongoDb and get the collection.
        coll = self.manager.db_connector.get_collection()
        print("Mongodb collection %s with count %d", coll, coll.count())

        start = time.time()

        for work in self:
            for task in work:
                results = task.get_results()
                pprint(results)
                results.update_collection(coll)
            results = work.get_results()
            pprint(results)
            results.update_collection(coll)

        msg = "MongoDb update done in %s [s]" % time.time() - start
        print(msg)

        results = self.get_results()
        pprint(results)
        results.update_collection(coll)

        # Update the pickle file to save the mongo ids.
        self.pickle_dump()

        for d in coll.find():
            pprint(d)

    def open_files(self, what="o", wti=None, status=None, op="==", editor=None):
        """
        Open the files of the flow inside an editor (command line interface).

        Args:
            what:
                string with the list of characters selecting the file type
                Possible choices:
                    i ==> input_file,
                    o ==> output_file,
                    f ==> files_file,
                    j ==> job_file,
                    l ==> log_file,
                    e ==> stderr_file,
                    q ==> qerr_file,
            wti:
                tuple with the (work, task_index) to select
                or string in the form w_start:w_stop,task_start:task_stop
            status:
                if not None, only the tasks with this status are select
            op:
                status operator. Requires status. A task is selected 
                if task.status op status evaluates to true.
            editor:
                Select the editor. None to use the default editor ($EDITOR shell env var)
        """
        #TODO: Add support for wti
        if wti is not None:
            raise NotImplementedError("wti option is not available!")

        def get_files(task, wi, ti):
            """Helper function used to select the files of a task."""
            choices = {
                "i": task.input_file,
                "o": task.output_file,
                "f": task.files_file,
                "j": task.job_file,
                "l": task.log_file,
                "e": task.stderr_file,
                "q": task.qerr_file,
                #"q": task.qout_file,
            }

            selected = []
            for c in what:
                try:
                    selected.append(getattr(choices[c], "path"))
                except KeyError:
                    warnings.warn("Wrong keyword %s" % c)

            return selected

        # Build list of files to analyze.
        files = []
        for (task, wi, ti) in self.iflat_tasks_wti(status=status, op=op):
            lst = get_files(task, wi, ti)
            if lst:
                files.extend(lst)

        #logger.info("Will edit %d files: %s" % (len(files), str(files)))
        return Editor(editor=editor).edit_files(files)

    def cancel(self):
        """
        Cancel all the tasks that are in the queue.

        Returns:
            Number of jobs cancelled, negative value if error
        """
        if self.has_chrooted:
            # TODO: Use paramiko to kill the job?
            warnings.warn("Cannot cancel the flow via sshfs!")
            return -1

        # If we are running with the scheduler, we must send a SIGKILL signal.
        pid_file = os.path.join(self.workdir, "_PyFlowScheduler.pid")
        if os.path.exists(pid_file):
            with open(pid_file, "r") as fh:
                pid = int(fh.readline())
                
            retcode = os.system("kill -9 %d" % pid)
            logger.info("Sent SIGKILL to the scheduler, retcode = %s" % retcode)
            try:
                os.remove(pid_file)
            except IOError:
                pass

        num_cancelled = 0
        for task in self.iflat_tasks():
            num_cancelled += task.cancel()

        return num_cancelled

    def get_njobs_in_queue(self, username=None):
        """
        returns the number of jobs in the queue,
        returns None when the number of jobs cannot be determined.

        Args:
            username: (str) the username of the jobs to count (default is to autodetect)
        """
        return self.manager.qadapter.get_njobs_in_queue(username=username)

    def rmtree(self, ignore_errors=False, onerror=None):
        """Remove workdir (same API as shutil.rmtree)."""
        shutil.rmtree(self.workdir, ignore_errors=ignore_errors, onerror=onerror)

    def build(self, *args, **kwargs):
        """Make directories and files of the `Flow`."""
        self.indir.makedirs()
        self.outdir.makedirs()
        self.tmpdir.makedirs()

        for work in self:
            work.build(*args, **kwargs)

    def build_and_pickle_dump(self):
        """
        Build dirs and file of the `Flow` and save the object in pickle format.

        Returns:
            0 if success
        """
        self.build()
        return self.pickle_dump()

    def pickle_dump(self):
        """
        Save the status of the object in pickle format.

        Returns:
            0 if success
        """
        if self.has_chrooted:
            warnings.warn("Cannot pickle_dump since we have chrooted from %s" % self.has_chrooted)
            return -1

        protocol = self.pickle_protocol
        filepath = os.path.join(self.workdir, self.PICKLE_FNAME)

        # Atomic transaction with FileLock.
        with FileLock(filepath):
            #with open(filepath, mode="wb") as fh:
            with AtomicFile(filepath, mode="wb") as fh:
                pmg_pickle_dump(self, fh, protocol=protocol)

        return 0

    def register_task(self, input, deps=None, manager=None, task_class=None):
        """
        Utility function that generates a `Workflow` made of a single task

        Args:
            input:
                Abinit Input file or `Strategy` object of `Task` object.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.
            task_class:
                Task subclass to instantiate. Default: `AbinitTask` 

        Returns:   
            The generated `Task`.
        """
        work = Workflow(manager=manager)
        task = work.register(input, deps=deps, task_class=task_class)
        self.register_work(work)

        return task

    def register_work(self, work, deps=None, manager=None, workdir=None):
        """
        Register a new `Workflow` and add it to the internal list, 
        taking into account possible dependencies.

        Args:
            work:
                `Workflow` object.
            deps:
                List of `Dependency` objects specifying the dependency of this node.
                An empy list of deps implies that this node has no dependencies.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the workflow.
            workdir:
                The name of the directory used for the `Workflow`.

        Returns:   
            The registered `Workflow`.
        """
        # Directory of the workflow.
        if workdir is None:
            work_workdir = os.path.join(self.workdir, "w" + str(len(self)))
        else:
            work_workdir = os.path.join(self.workdir, os.path.basename(workdir))

        work.set_workdir(work_workdir)

        if manager is not None:
            work.set_manager(manager)

        self.works.append(work)

        if deps:
            deps = [Dependency(node, exts) for node, exts in deps.items()]
            work.add_deps(deps)

        return work

    def register_work_from_cbk(self, cbk_name, cbk_data, deps, work_class, manager=None):
        """
        Registers a callback function that will generate the Tasks of the `Workflow`.

        Args:
            cbk_name:
                Name of the callback function (must be a bound method of self)
            cbk_data
                Additional data passed to the callback function.
            deps:
                List of `Dependency` objects specifying the dependency of the workflow.
            work_class:
                `Workflow` class to instantiate.
            manager:
                The `TaskManager` responsible for the submission of the task. 
                If manager is None, we use the `TaskManager` specified during the creation of the `Flow`.
                                                                                                            
        Returns:   
            The `Workflow` that will be finalized by the callback.
        """
        # TODO: pass a workflow factory instead of a class
        # Directory of the workflow.
        work_workdir = os.path.join(self.workdir, "w" + str(len(self)))

        # Create an empty workflow and register the callback
        work = work_class(workdir=work_workdir, manager=manager)
        
        self._works.append(work)
                                                                                                            
        deps = [Dependency(node, exts) for node, exts in deps.items()]
        if not deps:
            raise ValueError("A callback must have deps!")

        work.add_deps(deps)

        # Wrap the callable in a Callback object and save 
        # useful info such as the index of the workflow and the callback data.
        cbk = FlowCallback(cbk_name, self, deps=deps, cbk_data=cbk_data)
        self._callbacks.append(cbk)
                                                                                                            
        return work

    def allocate(self):
        """
        Allocate the `AbinitFlow` i.e. assign the `workdir` and (optionally) 
        the `TaskManager` to the different tasks in the Flow.
        """
        for work in self:
            # Each workflow has a reference to its flow.
            work.allocate(manager=self.manager)
            work.set_flow(self)
            # Each task has a reference to its workflow.
            for task in work:
                task.set_work(work)

        self.check_dependencies()

        return self

    def show_dependencies(self):
        for work in self:
            work.show_intrawork_deps()

    def on_dep_ok(self, signal, sender):
        # TODO
        # Replace this callback with dynamic dispatch
        # on_all_S_OK for workflow
        # on_S_OK for task
        logger.info("on_dep_ok with sender %s, signal %s" % (str(sender), signal))

        for i, cbk in enumerate(self._callbacks):

            if not cbk.handle_sender(sender):
                logger.info("%s does not handle sender %s" % (cbk, sender))
                continue

            if not cbk.can_execute():
                logger.info("Cannot execute %s" % cbk)
                continue 

            # Execute the callback and disable it
            logger.info("about to execute callback %s" % str(cbk))
            cbk()
            cbk.disable()

            # Update the database.
            self.pickle_dump()

    #def finalize(self):
    #    """This method is called when the flow is completed."""

    def connect_signals(self):
        """
        Connect the signals within the workflow.
        self is responsible for catching the important signals raised from 
        its task and raise new signals when some particular condition occurs.
        """
        # Connect the signals inside each Workflow.
        for work in self:
            work.connect_signals()

        # Observe the nodes that must reach S_OK in order to call the callbacks.
        for cbk in self._callbacks:
            for dep in cbk.deps:
                logger.info("connecting %s \nwith sender %s, signal %s" % (str(cbk), dep.node, dep.node.S_OK))
                dispatcher.connect(self.on_dep_ok, signal=dep.node.S_OK, sender=dep.node, weak=False)

        # Associate to each signal the callback _on_signal
        # (bound method of the node that will be called by `AbinitFlow`
        # Each node will set its attribute _done_signal to True to tell
        # the flow that this callback should be disabled.

        # Register the callbacks for the Workflows.
        #for work in self:
        #    slot = self._sig_slots[work]
        #    for signal in S_ALL:
        #        done_signal = getattr(work, "_done_ " + signal, False)
        #        if not done_sig:
        #            cbk_name = "_on_" + str(signal)
        #            cbk = getattr(work, cbk_name, None)
        #            if cbk is None: continue
        #            slot[work][signal].append(cbk)
        #            print("connecting %s\nwith sender %s, signal %s" % (str(cbk), dep.node, dep.node.S_OK))
        #            dispatcher.connect(self.on_dep_ok, signal=signal, sender=dep.node, weak=False)

        # Register the callbacks for the Tasks.

        #self.show_receivers()

    def show_receivers(self, sender=None, signal=None):
        sender = sender if sender is not None else dispatcher.Any
        signal = signal if signal is not None else dispatcher.Any
        print("*** live receivers ***")
        for rec in dispatcher.liveReceivers(dispatcher.getReceivers(sender, signal)):
            print("receiver -->", rec)
        print("*** end live receivers ***")

    #def get_results(self, **kwargs)

    def rapidfire(self, check_status=False, **kwargs):
        """
        Use PyLauncher to submits tasks in rapidfire mode.
        kwargs contains the options passed to the launcher.

        Return the number of tasks submitted.
        """
        from .launcher import PyLauncher

        if check_status:
            self.check_status()

        return PyLauncher(self, **kwargs).rapidfire()

    def make_scheduler(self, **kwargs):
        """
        Build a return a scheduler to run the flow.

        kwargs:
            if empty we use the user configuration file.
            if filepath in kwargs we init the scheduler from file.
            else pass **kwargs to PyFlowScheduler.__init__
        """
        from .launcher import PyFlowScheduler

        if not kwargs:
            # User config if kwargs is empty
            sched = PyFlowScheduler.from_user_config()
        else:
            # Use from_file if filepath if present, else call __init__
            filepath = kwargs.pop("filepath", None)
            if filepath is not None:
                assert not kwargs
                sched = PyFlowScheduler.from_file(filepath)
            else:
                sched = PyFlowScheduler.from_file(**kwargs)

        sched.add_flow(self)
        return sched


class G0W0WithQptdmFlow(AbinitFlow):
    """
    Build an `AbinitFlow` for one-shot G0W0 calculations.
    The computation of the q-points for the screening is parallelized with qptdm
    i.e. we run independent calculation for each q-point and then we merge
    the final results.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).
        scr_input:
            Input for the SCR run.
        sigma_inputs:
            Input(s) for the SIGMA run(s).
    """
    def __init__(self, workdir, manager, scf_input, nscf_input, scr_input, sigma_inputs):
        super(G0W0WithQptdmFlow, self).__init__(workdir, manager)

        # Register the first workflow (GS + NSCF calculation)
        bands_work = self.register_work(BandStructureWorkflow(scf_input, nscf_input))

        # Register the callback that will be executed the workflow for the SCR with qptdm.
        scr_work = self.register_work_from_cbk(cbk_name="cbk_qptdm_workflow", cbk_data={"input": scr_input},
                                               deps={bands_work.nscf_task: "WFK"}, work_class=QptdmWorkflow)

        # The last workflow contains a list of SIGMA tasks
        # that will use the data produced in the previous two workflows.
        if not isinstance(sigma_inputs, (list, tuple)):
            sigma_inputs = [sigma_inputs]

        sigma_work = Workflow()
        for sigma_input in sigma_inputs:
            sigma_work.register(sigma_input, deps={bands_work.nscf_task: "WFK", scr_work: "SCR"})
        self.register_work(sigma_work)

        self.allocate()

    def cbk_qptdm_workflow(self, cbk):
        """
        This callback is executed by the flow when bands_work.nscf_task reaches S_OK.

        It computes the list of q-points for the W(q,G,G'), creates nqpt tasks
        in the second workflow (QptdmWorkflow), and connect the signals.
        """
        scr_input = cbk.data["input"]
        # Use the WFK file produced by the second
        # Task in the first Workflow (NSCF step).
        nscf_task = self[0][1]
        wfk_file = nscf_task.outdir.has_abiext("WFK")

        work = self[1]
        work.set_manager(self.manager)
        work.create_tasks(wfk_file, scr_input)
        work.add_deps(cbk.deps)
        work.connect_signals()
        work.build()

        return work


class FlowCallbackError(Exception):
    """Exceptions raised by FlowCallback."""


class FlowCallback(object):
    """
    This object implements the callbacks executeed by the flow when
    particular conditions are fulfilled. See on_dep_ok method of Flow.

    .. note:
        I decided to implement callbacks via this object instead of a standard
        approach based on bound methods because:

            1) pickle (v<=3) does not support the pickling/unplickling of bound methods

            2) There's some extra logic and extra data needed for the proper functioning
               of a callback at the flow level and this object provides an easy-to-use interface.
    """
    Error = FlowCallbackError

    def __init__(self, func_name, flow, deps, cbk_data):
        """
        Args:
            func_name:
                String with the name of the callback to execute.
                func_name must be a bound method of flow with signature:

                    func_name(self, cbk)

                where self is the Flow instance and cbk is the callback
            flow:
                Reference to the `Flow`
            deps:
                List of dependencies associated to the callback
                The callback is executed when all dependencies reach S_OK.
            cbk_data:
                Dictionary with additional data that will be passed to the callback via self.
        """
        self.func_name = func_name
        self.flow = flow
        self.deps = deps
        self.data = cbk_data or {}
        self._disabled = False

    def __str__(self):
        return "%s: %s bound to %s" % (self.__class__.__name__, self.func_name, self.flow)

    def __call__(self):
        """Execute the callback."""
        if self.can_execute():
            # Get the bound method of the flow from func_name.
            # We use this trick because pickle (format <=3)
            # does not support bound methods.
            try:
                func = getattr(self.flow, self.func_name)
            except AttributeError as exc:
                raise self.Error(str(exc))

            return func(self)

        else:
            raise self.Error("You tried to __call_ a callback that cannot be executed!")

    def can_execute(self):
        """True if we can execute the callback."""
        return not self._disabled and [dep.status == dep.node.S_OK for dep in self.deps]

    def disable(self):
        """
        True if the callback has been disabled.
        This usually happens when the callback has been executed.
        """
        self._disabled = True

    def handle_sender(self, sender):
        """
        True if the callback is associated to the sender
        i.e. if the node who sent the signal appears in the 
        dependencies of the callback.
        """
        return sender in [d.node for d in self.deps]


# Factory functions.
def bandstructure_flow(workdir, manager, scf_input, nscf_input, dos_inputs=None):
    """
    Build an `AbinitFlow` for band structure calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).
        dos_inputs:
            Input(s) for the NSCF run (dos run).

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = BandStructureWorkflow(scf_input, nscf_input, dos_inputs=dos_inputs)
    flow.register_work(work)
    return flow.allocate()


def g0w0_flow(workdir, manager, scf_input, nscf_input, scr_input, sigma_inputs):
    """
    Build an `AbinitFlow` for one-shot $G_0W_0$ calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` object used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        nscf_input:
            Input for the NSCF run (band structure run).
        scr_input:
            Input for the SCR run.
        sigma_inputs:
            List of inputs for the SIGMA run.

    Returns:
        `AbinitFlow`
    """
    flow = AbinitFlow(workdir, manager)
    work = G0W0_Workflow(scf_input, nscf_input, scr_input, sigma_inputs)
    flow.register_work(work)
    return flow.allocate()


def phonon_flow(workdir, manager, scf_input, ph_inputs, with_nscf=False, ana_input=None):
    """
    Build an `AbinitFlow` for phonon calculations.

    Args:
        workdir:
            Working directory.
        manager:
            `TaskManager` used to submit the jobs
        scf_input:
            Input for the GS SCF run.
        ph_inputs:
            List of Inputs for the phonon runs.

    Returns:
        `AbinitFlow`
    """
    natom = len(scf_input.structure)

    # Create the container that will manage the different workflows.
    flow = AbinitFlow(workdir, manager)

    # Register the first workflow (GS calculation)
    scf_task = flow.register_task(scf_input, task_class=ScfTask)

    if with_nscf:
        nscf_input = copy.deepcopy(scf_input)
        nscf_input.set_variable('iscf', -3)

    # Build a temporary workflow with a shell manager just to run
    # ABINIT to get the list of irreducible pertubations for this q-point.
    shell_manager = manager.to_shell_manager(mpi_procs=1)

    if not isinstance(ph_inputs, (list, tuple)):
        ph_inputs = [ph_inputs]

    for i, ph_input in enumerate(ph_inputs):
        fake_input = ph_input.deepcopy()

        # Run abinit on the front-end to get the list of irreducible pertubations.
        tmp_dir = os.path.join(workdir, "__ph_run" + str(i) + "__")
        w = PhononWorkflow(workdir=tmp_dir, manager=shell_manager)
        fake_task = w.register(fake_input)

        # Use the magic value paral_rf = -1 to get the list of irreducible perturbations for this q-point.
        vars = dict(paral_rf=-1,
                    rfatpol=[1, natom],  # Set of atoms to displace.
                    rfdir=[1, 1, 1],     # Along this set of reduced coordinate axis.
                    )

        fake_task.strategy.add_extra_abivars(vars)
        w.allocate()
        w.start(wait=True)

        # Parse the file to get the perturbations.
        irred_perts = yaml_read_irred_perts(fake_task.log_file.path)
        print(irred_perts)

        w.rmtree()

        # Now we can build the final list of workflows:
        # One workflow per q-point, each workflow computes all 
        # the irreducible perturbations for a singe q-point.

        if with_nscf:
            nscf_input.set_variable('nqpt', 1)
            nscf_input.set_variable('qpt', irred_perts[0]['qpt'])
            nscf_task = flow.register_task(nscf_input, deps={scf_task: "DEN"}, task_class=NscfTask)

        work_qpt = PhononWorkflow()

        for irred_pert in irred_perts:
            print(irred_pert)
            new_input = ph_input.deepcopy()

            #rfatpol   1 1   # Only the first atom is displaced
            #rfdir   1 0 0   # Along the first reduced coordinate axis
            qpt = irred_pert["qpt"]
            idir = irred_pert["idir"]
            ipert = irred_pert["ipert"]

            # TODO this will work for phonons, but not for the other types of perturbations.
            rfdir = 3 * [0]
            rfdir[idir -1] = 1
            rfatpol = [ipert, ipert]

            new_input.set_variables(
                #rfpert=1,
                qpt=qpt,
                rfdir=rfdir,
                rfatpol=rfatpol,
            )

            if with_nscf:
                work_qpt.register(new_input, deps={nscf_task: "WFQ", scf_task: "WFK"}, task_class=PhononTask)
            else:
                work_qpt.register(new_input, deps={scf_task: "WFK"}, task_class=PhononTask)

        flow.register_work(work_qpt)

        #if ana_input is not None:
        #    merge_input = {}
        #    qp_merge_task = flow.register_task(merge_input, deps={work_qpt: "DDB"}, task_class=QpMergeTask)
        #    flow.register_task(ana_input, deps={qp_merge_task: "DDB"}, task_class=AnaddbTask)
                                            
    return flow.allocate()
