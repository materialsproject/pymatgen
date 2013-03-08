"""
Classes for defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import collections
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder
from pymatgen.io.smartio import read_structure
from pymatgen.util.string_utils import stream_has_colours

from .netcdf import GSR_Reader
from .pseudos import PseudoDatabase, PseudoTable, PseudoExtraInfo
from .abinit_input import Input, Electrons, System, Control, Kpoints
from .utils import parse_ewc, abinit_output_iscomplete

#from abipy.htc.jobfile import JobFile
#from ..eos import EOS

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

__all__ = [
    "Workflow",
]

##########################################################################################

class AbinitTask(object):

    # Prefixes for Abinit (input, output, temporary) files.
    Prefix = collections.namedtuple("Prefix", "idata odata tdata")
    pj = os.path.join

    prefix = Prefix("in", pj("output","out"), pj("temporary","tmp"))
    del Prefix, pj

    # Basenames for Abinit input and output files.
    Basename = collections.namedtuple("Basename", "input output files_file log_file stderr_file jobfile")
    basename = Basename("run.input", "run.output", "run.files", "log", "stderr", "job.sh")
    del Basename

    def __init__(self, input, workdir, 
                 varpaths    = None, 
                 user_options= None,
                ):

        self.workdir = os.path.abspath(workdir)

        self.input = input.copy()

        self.need_files = []
        if varpaths is not None: 
            self.need_files.extend(varpaths.values())
            self.input = self.input.add_vars_to_control(varpaths)

        # Files required for the execution.
        self.input_file  = File(self.basename.input      , dirname = self.workdir)
        self.output_file = File(self.basename.output     , dirname = self.workdir)
        self.files_file  = File(self.basename.files_file , dirname = self.workdir)
        self.log_file    = File(self.basename.log_file   , dirname = self.workdir)
        self.stderr_file = File(self.basename.stderr_file, dirname = self.workdir)

        # Find number of processors ....

        # Set jobfile variables.
        # TODO Rename JobFile to ShellScript(File)
        #self.jobfile = File(self.workdir, self.basename.jobfile)
        #excecutable = "abinit"
        self.jobfile = JobFile(name   = self.jobfile_path, 
                               input  = self.files_file.path, 
                               log    = self.log_file.path,
                               stderr = self.stderr_file.path,
                              )
    def __str__(self):
        lines = []
        app = lines.append

        #basenames = ["input", "files_file", "jobfile"]
        #for s in basenames:
        #    apath = getattr(self, s + "_path") 
        #    app(s + ": " + apath)

        app("outfiles_dir: " + self.outfiles_dir)
        app("tmpfiles_dir: " + self.tmpfiles_dir)

        return "\n".join(lines)

    def __repr__(self):
        return "<%s at %s, task_workdir = %s>" % (
            self.__class__.__name__, id(self), os.path.basename(self.workdir))

    @property
    def name(self):
        return self.input_file.basename

    @property
    def jobfile_path(self):
        "Absolute path of the job file (shell script)."
        return os.path.join(self.workdir, self.basename.jobfile)    

    @property
    def outfiles_dir(self):
        head, tail = os.path.split(self.prefix.odata)
        return os.path.join(self.workdir, head)

    @property
    def tmpfiles_dir(self):
        head, tail = os.path.split(self.prefix.tdata)
        return os.path.join(self.workdir, head)

    def odata_path_from_ext(self, ext):
        "Return the path of the output file with extension ext"
        if ext[0] != "_": ext = "_" + ext
        return os.path.join(self.workdir, (self.prefix.odata + ext))

    @property
    def pseudos(self):
        "List of pseudos used in the calculation."
        return self.input.pseudos

    @property
    def filesfile_string(self):
        "String with the list of files and prefixes needed to execute abinit."
        lines = []
        app = lines.append
        pj = os.path.join

        app(self.input_file.path)                 # Path to the input file
        app(self.output_file.path)                # Path to the output file
        app(pj(self.workdir, self.prefix.idata))  # Prefix for in data
        app(pj(self.workdir, self.prefix.odata))  # Prefix for out data
        app(pj(self.workdir, self.prefix.tdata))  # Prefix for tmp data

        # Paths to the pseudopotential files.
        for pseudo in self.pseudos:
            app(pseudo.path)

        return "\n".join(lines)

    @property
    def to_dict(self):
        raise NotimplementedError("")
        d = {k: v.to_dict for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["input"] = self.input 
        d["workdir"] = workdir 
        return d
                                                                    
    @staticmethod
    def from_dict(d):
        raise NotimplementedError("")
        return AbinitTask(**d)
        #return AbinitTask(d["input"], d["workdir"],)

    def iscomplete(self):
        "True if the output file is complete."
        return abinit_output_iscomplete(self.output_file.path)

    @property
    def isnc(self):
        return all(p.isnc for p in self.pseudos)

    @property
    def ispaw(self):
        return all(p.ispaw for p in self.pseudos)

    #def in_files(self):
    #    "Return all the input data files."
    #    files = list()
    #    for file in os.listdir(dirname(self.idat_root)):
    #        if file.startswith(basename(self.idat_root)):
    #            files.append(join(dirname(self.idat_root), file))
    #    return files
                                                                  
    def outfiles(self):
        "Return all the output data files produced."
        files = list()
        for file in os.listdir(self.outfiles_dir):
            if file.startswith(os.path.basename(self.prefix.odata)):
                files.append(os.path.join(self.outfiles_dir, file))
        return files
                                                                  
    def tmpfiles(self):
        "Return all the input data files produced."
        files = list()
        for file in os.listdir(self.tmpfiles_dir):
            if file.startswith(os.path.basename(self.prefix.tdata)):
                files.append(os.path.join(self.tmpfiles_dir, file))
        return files

    def path_in_workdir(self, filename):
        "Create the absolute path of filename in the workind directory."
        return os.path.join(self.workdir, filename)

    def make_dirs_and_files(self, *args, **kwargs):
        """
        Write Abinit input files and directory.
        """
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        if not os.path.exists(self.outfiles_dir):
            os.makedirs(self.outfiles_dir)

        if not os.path.exists(self.tmpfiles_dir):
            os.makedirs(self.tmpfiles_dir)

        self.input_file.write(str(self.input))
        self.files_file.write(self.filesfile_string)

        with open(self.jobfile_path, "w") as f:
                f.write(str(self.jobfile))

    def destroy(self, *args, **kwargs):
        """
        Remove all calculation files and directories.
                                                                                   
        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.pop('verbose', 0):
            print('Removing directory tree: %s' % self.workdir)

        _remove_tree(self.workdir, **kwargs)

    def read_mainlog_ewc(self, nafter=5):
       """
       Read errors, warnings and comments from the main output and the log file.
                                                                                        
       :return: Two namedtuple instances: main, log. 
                The lists of strings with the corresponding messages are
                available in main.errors, main.warnings, main.comments, log.errors etc.
       """
       main = parse_ewc(self.output_file.path, nafter=nafter)
       log  = parse_ewc(self.log_file.path, nafter=nafter)
       return main, log

    def get_status(self):
        return TaskStatus(self)

    def show_status(self, stream=sys.stdout, *args, **kwargs):
        self.status().show(stream=stream, *args, **kwargs)

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations. 
        The default implementation creates the workspace directory and the input files.
        """

    def teardown(self, *args, **kwargs):
        """
        Method called once the calculations completed.
        The default implementation does nothing.
        """

    def get_parallel_hints(self, max_numproc):
        """
        Run abinit in sequential to obtain a set of possible
        configuration for the number of processors
        """
        raise NotImplementedError("")

        # number of MPI processors, number of OpenMP processors, memory estimate in Gb.
        ParaHint = collections.namedtuple('BaseParallelHint', "mpi_nproc omp_nproc memory_gb")

        return hint

    def run(self, *args, **kwargs):
        """
        Run the calculation by executing the job file from the shell.

        Keyword arguments:
            verbose: (0)
                Print message if verbose is not zero.

        .. warning::
             This method must be thread safe since we may want to run several indipendent
             calculations with different python threads. 
        """
        if kwargs.get('verbose'):
            print('Running ' + self.input_path)

        self.setup(*args, **kwargs)

        import subprocess
        returncode = subprocess.call((self.jobfile.shell, self.jobfile.path), cwd=self.workdir)

        self.teardown(*args, **kwargs)

        return returncode

##########################################################################################

class TaskStatus(object):
    """
    Object used to inquire the status of the task and to have access 
    to the error and warning messages
    """
    OK   = 1 
    WAIT = 2 
    RUN  = 4
    ERR  = 8

    #: Rank associated to the different status (in order of increasing priority)
    _level2rank = {
      "completed" : OK, 
      "unstarted" : WAIT, 
      "running"   : RUN, 
      "error"     : ERR,
    }

    levels = Enum(s for s in _level2rank)

    def __init__(self, task):

        self._task = task

        if task.iscomplete():
            # The main output file seems completed.
            level_str = 'completed'

            #TODO
            # Read the number of errors, warnings and comments 
            # for the (last) main output and the log file.
            #main, log =  task.read_mainlog_ewc()
                                                                
            #main_info = main.tostream(stream)
            #log_info  = log.tostream(stream)


        # TODO err file!
        elif task.output_file.exists and task.log_file.exists:
            # Inspect the main output and the log file for ERROR messages.
            main, log = task.read_mainlog_ewc()

            level_str = 'running'
            if main.errors or log.errors:
                level_str = 'error'

            with open(task.stderr_file.path, "r") as f:
                lines = f.readlines()
                if lines: level_str = 'error'
        else:
            level_str = 'unstarted'

        self._level_str = level_str

        assert self._level_str in TaskStatus.levels

    @property
    def rank(self):
        return self._level2rank[self._level_str]

    @property
    def is_completed(self):
        return self._level_str == "completed"

    @property
    def is_unstarted(self):
        return self._level_str == "unstarted"

    @property
    def is_running(self):
        return self._level_str == "running"

    @property
    def is_error(self):
        return self._level_str == "error"

    def __repr__(self):
        return "<%s at %s, task = %s, level = %s>" % (
            self.__class__.__name__, id(self), repr(self._task), self._level_str)

    def __str__(self):
        return str(self._task) + self._level_str

    # Rich comparison support. Mainly used for selecting the 
    # most critical status when we have several tasks executed inside a workflow.
    def __lt__(self, other):
        return self.rank < other.rank

    def __le__(self, other):
        return self.rank <= other.rank

    def __eq__(self, other):
        return self.rank == other.rank

    def __ne__(self, other):
        return not self == other

    def __gt__(self, other):
        return self.rank > other.rank

    def __ge__(self, other):
        return self.rank >= other.rank

    def show(self, stream=sys.stdout, *args, **kwargs):

        from ..tools import StringColorizer
        str_colorizer = StringColorizer(stream)
                                                                       
        _status2txtcolor = {
          "completed" : lambda string : str_colorizer(string, "green"),
          "error"     : lambda string : str_colorizer(string, "red"),
          "running"   : lambda string : str_colorizer(string, "blue"),
          "unstarted" : lambda string : str_colorizer(string, "cyan"),
        }

        color = lambda stat : str
        if stream_has_colours(stream):
            color = lambda stat : _status2txtcolor[stat](stat)

        stream.write( self._task.name + ': ' + color(self._level_str) + "\n")

    def print_ewc_messages(self, stream=sys.stdout, **kwargs):

        verbose = kwargs.pop("verbose", 0)

        # Read the number of errors, warnings and comments 
        # for the (last) main output and the log file.
        main, log =  self._task.read_mainlog_ewc()
                                                            
        main_info = main.tostream(stream)
        log_info  = log.tostream(stream)
                                                            
        stream.write("\n".join([main_info, log_info]) + "\n")

        if verbose:
            for w in main.warnings: 
                stream.write(w + "\n")
            if verbose > 1:
                for w in log.warnings: 
                    stream.write(w + "\n")

        if self.is_error:
            for e in log.errors: 
                    stream.write(e + "\n")

            if not log.errors: # Read the standard error.
                with open(self._task.stderr_file.path, "r") as f:
                    lines = f.readlines()
                    stream.write("\n".join(lines))

##########################################################################################

class TaskDependencies(object):
    """
    This object describes the depencies among the AbinitTask instances 
    belonging to a Workflow.
    """

    _ext2getvars = {
        "DEN" : "getden_path",
        "WFK" : "getwfk_path",
        "SCR" : "getscr_path",
        "QPS" : "getqps_path",
    }

    def __init__(self, task, task_id, *odata_required):
        self.task    = task
        self.task_id = task_id

        self.odata_required = []
        if odata_required:
            self.odata_required = odata_required

    def with_odata(self, *odata_required):
        return TaskDependencies(self.task, self.task_id, *odata_required)

    def get_varpaths(self):
        vars = {}
        for ext in self.odata_required:
            varname = self._ext2getvars[ext]
            path = self.task.odata_path_from_ext(ext)
            vars.update({varname : path})
        return vars

##########################################################################################

class Workflow(MSONable):
    """
    A workflow is a list of (possibly connected) tasks.
    """
    def __init__(self, workdir, user_options=None):
        """
        Args:
            workdir:
            user_options:
        """
        self.workdir = os.path.abspath(workdir)

        # TODO Use default configuration.
        #if os.path.exists(self.workdir):
        #    raise ValueError("Workir %s already exists!")

        if user_options is not None:
            self.user_options = user_options
        else:
            # TODO Use default configuration.
            self.user_options = user_options

        self._tasks = []

        self._deps = collections.OrderedDict()

    def __len__(self):
        return len(self._tasks)

    def __iter__(self):
        return self._tasks.__iter__()

    def __getitem__(self, slice):
        return self._tasks[slice]

    def __repr__(self):
        return "<%s at %s, workdir = %s>" % (self.__class__.__name__, id(self), str(self.workdir))

    @property
    def to_dict(self):
        raise NotImplementedError("")
        #d = {k: v for k, v in self.items()}
        #d["@module"] = self.__class__.__module__
        #d["@class"] = self.__class__.__name__
        #return d

    @classmethod
    def from_dict(cls, d):
        raise NotImplementedError("")
        #i = cls()
        #for (k, v) in d.items():
        #    if k not in ("@module", "@class"):
        #        i[k] = v
        #return i

    @property
    def isnc(self):
        "True if norm-conserving calculation"
        return all(task.isnc for task in self)
                                                
    @property
    def ispaw(self):
        "True if PAW calculation"
        return all(task.ispaw for task in self)

    def path_in_workdir(self, filename):
        "Create the absolute path of filename in the workind directory."
        return os.path.join(self.workdir, filename)

    def setup(self, *args, **kwargs):
        """
        Method called before running the calculations. 
        The default implementation creates the workspace directory and the input files.
        """
        self.make_dirs_and_files(*args, **kwargs)

    def teardown(self, *args, **kwargs):
        """
        Method called once the calculations completed.
        The default implementation does nothing.
        """

    def show_inputs(self, stream=sys.stdout):
        lines = []
        app = lines.append

        width = 120
        for task in self:
            app("\n")
            app(repr(task))
            app("\ninput: %s" % task.input_file.path)
            app("\n")
            app(str(task.input))
            app(width*"=" + "\n")

        stream.write("\n".join(lines))

    def register_input(self, input, depends=()):

        task_workdir = os.path.join(self.workdir, "task_" + str(len(self) + 1))
        
        # Handle possible dependencies.
        varpaths = None
        if depends:
            if not isinstance(depends, collections.Iterable): 
                depends = [depends]

            varpaths = {}
            for dep in depends:
                varpaths.update( dep.get_varpaths() )
                #print("varpaths %s" % str(varpaths))

        new_task = AbinitTask(input, task_workdir, 
                              varpaths     = varpaths,
                              user_options = self.user_options, 
                             )

        # Add it to the list and return the ID of the task 
        # so that client code can specify possible dependencies.
        #if new_task in self._deps:
        #    raise ValueError("task is already registered!")

        #if new_task in self:
            #print "task is already in self!"
            #raise ValueError("task is already in self!")

        self._tasks.append(new_task)

        newtask_id = len(self)
        self._deps[newtask_id] = depends

        return TaskDependencies(new_task, newtask_id)

    def make_dirs_and_files(self, *args, **kwargs):
        self._tasks = tuple(self._tasks)

        # Create top level directory.
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # Create task workdirs and files.
        for task in self:
            task.make_dirs_and_files(*args, **kwargs)

    def get_status(self, only_highest_rank=True):
        "Get the status of the tasks in self."

        status_list = [task.get_status() for task in self]

        if only_highest_rank:
            return max(status_list)
        else:
            return status_list

    def destroy(self, *args, **kwargs):
        """
        Remove all calculation files and directories.
                                                                                   
        Keyword arguments:
            force: (False)
                Do not ask confirmation.
            verbose: (0)
                Print message if verbose is not zero.
        """
        if kwargs.pop('verbose', 0):
            print('Removing directory tree: %s' % self.workdir)
                                                                                    
        _remove_tree(self.workdir, **kwargs)

    def start(self, *args, **kwargs):
        """
        Start the workflow. Call setup first, then the 
        the tasks are executed. Finally, the teardown method
        is called.

            Args:
                num_pythreads = Number of python threads to use (defaults to 1)
        """

        num_pythreads = int(kwargs.pop("num_pythreads", 1))

        self.setup(*args, **kwargs)

        # TODO Run only the calculations that are not done 
        # and whose dependencies are satisfied.
        if num_pythreads == 1:
            for task in self:
                task.run(*args, **kwargs)

        else:
            # Threaded version.
            from threading import Thread
            from Queue import Queue
            print("Threaded version with num_pythreads %s" % num_pythreads)

            def worker():
                while True:
                    func, args, kwargs = q.get()
                    func(*args, **kwargs)
                    q.task_done()
                                                         
            q = Queue()
            for i in range(num_pythreads):
                t = Thread(target=worker)
                t.setDaemon(True)
                t.start()

            for task in self:
                q.put((task.run, args, kwargs))
                                                                    
            # Block until all tasks are done. 
            q.join()  

        self.teardown(*args, **kwargs)

##########################################################################################

class PseudoEcutTest_Workflow(Workflow):

    def __init__(self, pseudo, ecut_list, 
                 workdir       = None,
                 spin_mode     = "unpolarized", 
                 smearing      = None,
                 acell         = 8, 
                ):
        
        user_options = {}

        self.pseudo = pseudo
        workdir = workdir if workdir else pseudo.basename

        Workflow.__init__(self, workdir, user_options=user_options)
        
        # Define System: one atom in a box of lenghts acell. 
        box = System.atom_in_periodic_box(pseudo, acell=acell)

        # Gamma-only sampling.
        gamma_only = Kpoints.gamma_only()

        # Default setup for electrons: no smearing.
        electrons = Electrons(spin_mode=spin_mode, smearing=smearing)

        # Generate inputs with different values of ecut and register the task.
        self.ecut_list = np.array(ecut_list)

        for ecut in self.ecut_list:
            control = Control(box, electrons, gamma_only,
                              ecut        = ecut, 
                              prtwf       = 0,
                              want_forces = False,
                              want_stress = False,
                             )

            self.register_input(Input(control, box, gamma_only, electrons))

    def start(self, *args, **kwargs):

        self.setup(*args, **kwargs)

        for task in self:
            pass
            #task.run(*args, **kwargs)

        self.teardown(*args, **kwargs)

    def teardown(self, *args, **kwargs):

        etotal = self.get_etotal()

        if not self.isnc:
            raise NotImplementedError("PAW convergence tests are not supported yet")

        etotal_dict = {1.0 : etotal}

        status = self.get_status()
        print(repr(status))

        #input = ["I really\n", "dislike\n", "XML\n"]

        # TODO handle possible problems with the SCF cycle
        pp_info = PseudoExtraInfo.from_data(self.ecut_list, etotal_dict)

        #pp_info.show_etotal()

        # Append PP_EXTRA_INFO section to the pseudo file
        # provided that the results seem ok. If results seem strange, 
        # the XML section is written on a separate file.
        xml_string = pp_info.toxml()
        #print self.pseudo.basename, xml_string

        #path = self.pseudo.path
        #if pp_info.strange_data is not None:
        #    path = self.pseudo.path + ".strange"

        with open(self.path_in_workdir("pp_info.xml"), "w") as fh:
            fh.write(xml_string)
        #new = pp_info.from_string(xml_string)
        #print new
        #print new.toxml()

    def get_etotal(self):

        etotal = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            gsr_path = task.odata_path_from_ext("GSR") 

            with GSR_Reader(gsr_path) as gsr_data:
                etotal.append(gsr_data.get_value("etotal"))

        return np.array(etotal)

##########################################################################################

#class KPathDescriptor(object):
#    # TODO
#    def __init__(self, kptbounds, ndivsm=20, labels=None, structure=None):
#        self.kptbounds = np.reshape(kptbounds, (-1,3))
#
#        if labels is not None:
#            assert len(labels) == len(kptbounds)
#            self._labels = labels
#        else:
#            if structure is None:
#                raise ValueError("lattice must be specified")
#
#    def get_labels(self):
#        raise NotImplementedError("")
#
#    def toabivars(self):
#        d = {"kptbounds" : kptbounds,
#             "ndivsm"    : num_kpts,
#             "kptopt"    : -len(kptbounds)+1,
#             #"labels"   : labels,  # TODO
#            }
#        return d

##########################################################################################

class BandStructure_Workflow(Workflow):

    def __init__(self, structure, pptable_or_pseudos, scf_ngkpt, nscf_nband, kpath_bounds, 
                 workdir   = None, 
                 spin_mode = "polarized",
                 smearing  = None,
                 ndivsm    = 20, 
                 dos_ngkpt = None
                ):

        user_options = {}

        workdir = workdir if workdir else structure.formula + "_Band_Structure"

        Workflow.__init__(self, workdir, user_options=user_options)

        scf_input = Input.SCF_groundstate(structure, pptable_or_pseudos, scf_ngkpt, 
                                          spin_mode = spin_mode,
                                          smearing  = smearing,
                                         )

        scf_dep = self.register_input(scf_input)

        nscf_input = Input.NSCF_kpath_from_SCF(scf_input, nscf_nband, kpath_bounds, ndivsm=ndivsm)

        self.register_input(nscf_input, depends=scf_dep.with_odata("DEN"))

        # Add DOS computation
        if dos_ngkpt is not None:

            dos_input = Input.NSCF_kmesh_from_SCF(scf_input, nscf_nband, dos_ngkpt)

            self.register_input(dos_input, depends=scf_dep.with_odata("DEN"))

##########################################################################################

class RelaxWorkflow(Workflow):

    def __init__(self, structure, pptable_or_pseudos, ngkpt,
                 workdir   = None, 
                 spin_mode = "polarized",
                 smearing  = None,
                ):
                                                                                                   
        user_options = {}
                                                                                                   
        workdir = workdir if workdir else structure.formula + "_Relax"
                                                                                                   
        Workflow.__init__(self, workdir, user_options=user_options)

        ions_relax = Relax(
                    iomov     = 3, 
                    optcell   = 2, 
                    dilatmx   = 1.1, 
                    ecutsm    = 0.5, 
                    ntime     = 80, 
                    strtarget = None,
                    strfact   = 100,
                    )                                                                                                                  

        ions_dep = self.register_input(ions_relax)

        cell_ions_relax = Relax()

        #self.register_input(cell_ions, depends=ions_dep.with_odata("WFK"))

##########################################################################################

def scale_lattice(lattice, new_volume):
    """
    Return a new lattice with volume new_volume by performing a 
    scaling of the lattice vectors so that length proportions and angles are preserved.
    """
    versors = lattice.matrix / lattice.abc

    geo_factor = abs(np.dot(np.cross(versors[0], versors[1]), versors[2]))

    ratios = lattice.abc / lattice.c

    new_c = (new_volume / ( geo_factor * np.prod(ratios))) ** (1/3.)

    new_lattice = Lattice(versors * (new_c * ratios))
    #assert abs(new_lattice.volume - new_volume) < 1.0e-12
    return new_lattice

##########################################################################################

class DeltaTest(Workflow):

    def __init__(self, structure_or_cif, pptable_or_pseudos, ngkpt,
                 workdir   = None, 
                 spin_mode = "polarized",
                 smearing  = None,
                 ecutsm    = 0.05,
                ):

        if isinstance(structure_or_cif, Structure):
            structure = structure_or_cif
        else:
            # Assume CIF file

            structure = read_structure(structure_or_cif)

        self._input_structure = structure

        user_options = {}
                                                                                                   
        workdir = workdir if workdir else structure.formula + "_EOS"
                                                                                                   
        Workflow.__init__(self, workdir, user_options=user_options)

        v0 = structure.volume

        self.volumes = v0 * np.arange(90, 112, 2) / 100.

        for vol in self.volumes:

            new_lattice = scale_lattice(structure.lattice, vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

            scf_input = Input.SCF_groundstate(new_structure, pptable_or_pseudos, ngkpt, 
                                              spin_mode = spin_mode,
                                              smearing  = smearing,
                                              # **kwargs
                                              ecutsm = ecutsm,
                                             )
            self.register_input(scf_input)

    def teardown(self):
        eos_fit = self.compute_eos_fit()
        #eos_fit.dump(self.path_in_workdir("eos_fit.json"))

    def compute_eos_fit(self):
        etotal = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            gsr_path = task.odata_path_from_ext("GSR") 
                                                            
            with GSR_Reader(gsr_path) as gsr_data:
                etotal.append(gsr_data.get_value("etotal"))
                                                            
        return EOS.Murnaghan.fit(self.volumes, Ha2eV(etotal))

##########################################################################################

class PversusV(Workflow):
    "Calculates P(V)."

    def __init__(self, structure_or_cif, pptable_or_pseudos, ngkpt,
                 workdir   = None, 
                 spin_mode = "polarized",
                 smearing  = None,
                 ecutsm    = 0.05,  # 1.0
                 dilatmx   = 1.1,
                 strtargets = None,
                ):
        raise NotImplementedError()

        user_options = {}
                                                                                                   
        workdir = workdir if workdir else structure.formula + "_PvsV"
                                                                                                   
        Workflow.__init__(self, workdir, user_options=user_options)

        if isinstance(structure_or_cif, Structure):
            structure = cif_or_stucture
        else:
            # Assume CIF file
            from pymatgen import read_structure
            structure = read_structure(structure_or_cif)

        #ndtset 5

        #strtarget1 3*3.4D-5 3*0.0
        #strtarget2 3*1.7D-5 3*0.0
        #strtarget3 3*0.0 3*0.0
        #strtarget4 3*-1.7D-5 3*0.0
        #strtarget5 3*-3.4D-5 3*0.0
        # -1.0 GPa, -0.5 GPa, 0.0 GPa, 0.5 GPa, and 1.0 GPa

        #ionmov 2
        #optcell 2
        #ntime 20
        #ecutsm 1.0
        #dilatmx 1.1
        #tolmxf 1.0D-6
        #toldff 1.0D-7 (or tolvrs 1.0D-12 if all ions are on special positions)

        for target in strtargets:

            input = Input.Relax(new_structure, pptable_or_pseudos, ngkpt, 
                                spin_mode = spin_mode,
                                smearing  = smearing,
                                # **kwargs
                                ecutsm = 0.05,
                                )

            self.register_input(scf_input)

    def teardown(self):

        pressures, volumes = [], []
        for task in self:
            # Open GSR file and read etotal (Hartree)
            gsr_path = task.odata_path_from_ext("GSR") 
                                                            
            with GSR_Reader(gsr_path) as gsr_data:

                pressures.append(gsr_data.get_value("pressure"))
                volumes.append(gsr_data.get_value("volume"))

##########################################################################################

class G0W0_Workflow(Workflow):

    def __init__(self, structure, pptable_or_pseudos, scf_ngkpt, nscf_ngkpt, 
                 ppmodel_or_freqmesh, ecuteps, ecutsigx, nband_screening, nband_sigma,
                 workdir   = None, 
                 spin_mode = "polarized",
                 smearing  = None,
                ):

        user_options = {}
        workdir = workdir if workdir else "G0W0_" + structure.formula

        Workflow.__init__(self, workdir, user_options=user_options)

        scf_input = Input.SCF_groundstate(structure, pptable_or_pseudos, scf_ngkpt, 
                                          spin_mode = spin_mode, 
                                          smearing  = smearing,
                                         )

        scf_dep = self.register_input(scf_input)

        max_nband = max(nband_screening, nband_sigma) 

        nband = int(max_nband + 0.05 * max_nband)

        nscf_input = Input.NSCF_kmesh_from_SCF(scf_input, nband, scf_ngkpt)

        nscf_dep = self.register_input(nscf_input, depends=scf_dep.with_odata("DEN"))

        screen_input = Input.SCR_from_NSCF(nscf_input, ecuteps, ppmodel_or_freqmesh, nband_screening, smearing=smearing)

        screen_dep = self.register_input(screen_input, depends=nscf_dep.with_odata("WFK"))

        kptgw = [0,0,0]
        bdgw = [1,4]
        sigma_input = Input.SIGMA_from_SCR(screen_input, nband_sigma, ecuteps, ecutsigx, kptgw, bdgw, smearing=smearing)

        depends = [nscf_dep.with_odata("WFK"), screen_dep.with_odata("SCR"),]

        self.register_input(sigma_input, depends=depends)

     #@staticmethod
     #def with_GodbyNeedsPPM()

     #@staticmethod
     #def with_contour_deformation()

##########################################################################################
# TODO 
# Should become an API common for the different codes that requires pseudos
def get_abinit_psp_dir(code="ABINIT"):
    import ConfigParser
    import pymatgen
    if code + "_PSP_DIR" in os.environ:
        return os.environ[code + "_PSP_DIR"]
    elif os.path.exists(os.path.join(os.path.dirname(pymatgen.__file__), "pymatgen.cfg")):
        module_dir = os.path.dirname(pymatgen.__file__)
        config = ConfigParser.SafeConfigParser()
        config.readfp(open(os.path.join(module_dir, "pymatgen.cfg")))
        return config.get(code, "pspdir")

    return None

def test_abinitio(structure):
    v = np.array([13.72, 14.83, 16.0, 17.23, 18.52])
    e = np.array([-56.29, -56.41, -56.46, -56.46, -56.42])

    #eos = EOS.Murnaghan()
    #fit = eos.fit(v, e)
    #print fit
    #sys.exit(1)

    psp_dir = get_abinit_psp_dir()

    pp_database = PseudoDatabase(dirpath=psp_dir, force_reload=False)

    GGA_HGHK_PPTABLE = pp_database.GGA_HGHK_PPTABLE

    pptable_or_pseudos = GGA_HGHK_PPTABLE
    pseudo = GGA_HGHK_PPTABLE[14][0]

    user_options = {}

    pptest = PseudoEcutTest_Workflow(pseudo, range(10,40,2))

    #pptest.show_inputs()
    pptest.start()
    sys.exit(1)

    scf_ngkpt = [4,4,4]
    kpath_bounds = [0,0,0, 0.5,0.5,0.5]
    nscf_nband = 333

    bands_wf = BandStructure_Workflow(structure, pptable_or_pseudos, scf_ngkpt, nscf_nband, kpath_bounds)
    #bands_wf.start()

    #bands_wf.show_inputs()
    #sys.exit(1)

    #for task in bands_wf:
    #    print task.input
    #    status = task.get_status()
    #    status.show()
    #    status.print_ewc_messages()
    #sys.exit(1)

    ecuteps   = 10
    ecutsigx  = 30
    nband_screening = 100
    nband_sigma = 50

    nscf_ngkpt = [2,2,2]

    ppmodel_or_freqmesh = PPModel.Godby()
    #ppmodel_or_freqmesh = ScreeningFrequencyMesh(nomega_real=None, maxomega_real=None, nomega_imag=None)

    g0w0_wf = G0W0_Workflow(structure, pptable_or_pseudos, scf_ngkpt, nscf_ngkpt, 
                            ppmodel_or_freqmesh, ecuteps, ecutsigx, nband_screening, nband_sigma)

    g0w0_wf.show_inputs()

    delta_wf = DeltaTest(structure, pseudo, [2,2,2])
    delta_wf.show_inputs()

##########################################################################################
# Helper functions.

def _remove_tree(*paths, **kwargs):
    import shutil
    ignore_errors = kwargs.pop("ignore_errors", True)
    onerror       = kwargs.pop("onerror", None)
                                                                           
    for path in paths:
        shutil.rmtree(path, ignore_errors=ignore_errors, onerror=onerror)

##########################################################################################

class File(object):
    """
    Very simple class used to store file basenames, absolute paths and directory names.

    Provides wrappers for the most commonly used os.path functions.
    """
    def __init__(self, basename, dirname="."):
        self.basename = basename
        self.dirname = os.path.abspath(dirname)
        self.path = os.path.join(self.dirname, self.basename)

    def __str__(self):
       return self.read()

    def __repr__(self):
        return "<%s at %s, %s>" % (self.__class__.__name__, id(self), self.path)

    @property
    def exists(self):
        "True if file exists."
        return os.path.exists(self.path)

    @property
    def isncfile(self):
        "True if self is a NetCDF file"
        return self.basename.endswith(".nc")

    def read(self):
        with open(self.path, "r") as f: return f.read()

    def readlines(self):
        with open(self.path, "r") as f: return f.readlines()

    def write(self, string):
        with open(self.path, "w") as f: return f.write(string)
                                        
    def writelines(self, lines):
        with open(self.path, "w") as f: return f.writelines()

##########################################################################################

if __name__ == "__main__":
    test_abinitio()
