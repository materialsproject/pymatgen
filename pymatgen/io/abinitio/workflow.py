"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import shutil
import abc
import json
import collections
import functools
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr, Ha2eV, Ha_eV, Ha2meV
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder
from pymatgen.io.smartio import read_structure
#from pymatgen.util.filelock import FileLock
from pymatgen.util.num_utils import iterator_from_slice 

from .netcdf import GSR_Reader
from .pseudos import Pseudo, PseudoDatabase, PseudoTable, PseudoExtraInfo, get_abinit_psp_dir
from .abinit_input import Input, Electrons, System, Control, Kpoints, Smearing
from .task import AbinitTask, TaskDependencies
from .utils import parse_ewc, abinit_output_iscomplete

from .jobfile import JobFile

#import logging
#logger = logging.getLogger(__name__)

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__email__ = "gmatteo at gmail.com"
__status__ = "Development"
__date__ = "$Feb 21, 2013M$"

#__all__ = [
#]

##########################################################################################

def map_method(method):
    "Decorator that calls item.method for all items in a iterable object."
    @functools.wraps(method)
    def wrapped(iter_obj, *args, **kwargs):
        return [getattr(item, method.__name__)(*args, **kwargs) for item in iter_obj]
    return wrapped


class WorkError(Exception):
    pass

class BaseWork(object):
    __metaclass__ = abc.ABCMeta

    Error = WorkError

    @abc.abstractproperty
    def processes(self):
        "Return a list of objects that support the subprocess.Popen protocol."

    # interface modeled after subprocess.Popen
    def poll(self):
        "Check if all child processes have terminated. Set and return returncode attribute."
        return [p.poll() for p in self.processes]

    def wait(self):
        "Wait for child processed to terminate. Set and return returncode attributes."
        return [p.wait() for p in self.processes]

    def communicate(self, input=None):
        """
        Interact with processes: Send data to stdin. Read data from stdout and stderr, until end-of-file is reached. 
        Wait for process to terminate. The optional input argument should be a string to be sent to the 
        child processed, or None, if no data should be sent to the children

        communicate() returns a list of tuples (stdoutdata, stderrdata).
        """
        return [p.communicate(input) for p in self.processes]

    def returncodes(self):
        """
        The children return codes, set by poll() and wait() (and indirectly by communicate()). 
        A None value indicates that the process hasn't terminated yet.
        A negative value -N indicates that the child was terminated by signal N (Unix only).
        """
        return [p.returncode for p in self.processes]

    @abc.abstractmethod
    def setup(self, *args, **kwargs):
        """
        Method called before submitting the calculations. 
        """
                                                         
    def _setup(self, *args, **kwargs):
        self.setup(*args, **kwargs)

    def get_results(self, *args, **kwargs):
        """
        Method called once the calculations completes.
        The base version returns a dictionary task_name : TaskResults for each task in self.
        """
        results = collections.OrderedDict()
        for task in self:
            results[task.name] = task.get_results(*args, **kwargs)
        return results
                                                         
##########################################################################################

class Work(BaseWork, MSONable):
    """
    A work is a list of (possibly connected) tasks.
    """
    Error = WorkError

    def __init__(self, workdir, **kwargs):
        """
        Args:
            workdir:
        """
        self.workdir = os.path.abspath(workdir)

        self.kwargs = kwargs

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
        The default implementation is empty.
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

    def _add_task(self, task, depends=()):
        self._tasks.append(task)
                                         
        task_id = len(self)

        if depends:
            self._deps[task_id] = depends

        return task_id

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
                              varpaths = varpaths,
                             )

        # Add it to the list and return the ID of the task 
        # so that client code can specify possible dependencies.
        #if new_task in self._deps:
        #    raise ValueError("task is already registered!")

        newtask_id = self._add_task(new_task, depends=depends)

        #self._tasks.append(new_task)
        #newtask_id = len(self)
        #self._deps[newtask_id] = depends

        return TaskDependencies(new_task, newtask_id)

    def build(self, *args, **kwargs):
        "Creates the top level directory"
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

    def get_status(self, only_highest_rank=False):
        "Get the status of the tasks in self."

        status_list = [task.get_status() for task in self]

        if only_highest_rank:
            return max(status_list)
        else:
            return status_list

    @property
    def processes(self):
        return [task.process for task in self]

    def rmtree(self, *args, **kwargs):
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
                                                                                    
        shutil.rmtree(self.workdir)

    def submit_tasks(self, *args, **kwargs):

        num_pythreads = int(kwargs.pop("num_pythreads", 1))

        # TODO Run only the calculations that are not done 
        # and whose dependencies are satisfied.
        if num_pythreads == 1:
            for task in self:
                task.start(*args, **kwargs)
                #stdoutdata, stderrdata = task.communicate()

        else:
            # Threaded version.
            from threading import Thread
            from Queue import Queue
            print("Threaded version with num_pythreads %s" % num_pythreads)

            def worker():
                while True:
                    task, args, kwargs = qin.get()
                    task.start(*args, **kwargs)
                    #stdoutdata, stderrdata = task.communicate()
                    qin.task_done()
                                                         
            qin = Queue()
            for i in range(num_pythreads):
                t = Thread(target=worker)
                t.setDaemon(True)
                t.start()

            for task in self:
                qin.put((task, args, kwargs))
                                                                    
            # Block until all tasks are done. 
            qin.join()  

    def start(self, *args, **kwargs):
        """
        Start the work. Call setup first, then the 
        the tasks are executed. Finally, the get_results method is called.

            Args:
                num_pythreads = Number of python threads to use (defaults to 1)
        """
        # Build dirs and files.
        self.build(*args, **kwargs)

        # Initial setup
        self._setup(*args, **kwargs)

        # Submit tasks (non-blocking)
        self.submit_tasks(*args, **kwargs)

    def read_etotal(self):
        """
        Reads the total energy from the GSR file produced by the task.

        Return a numpy array with the total energies in Hartree 
        The array element is set to np.inf if an exception is raised while reading the GSR file.
        """
        etotal = []
        for task in self:
            # Open the GSR file and read etotal (Hartree)
            gsr_path = task.odata_path_from_ext("GSR") 
            try:
                with GSR_Reader(gsr_path) as ncdata:
                    etotal.append(ncdata.get_value("etotal"))
            except:
                etotal.append(np.inf)
                                                            
        return np.array(etotal)

##########################################################################################

class IterativeWork(Work):
    """
    TODO
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self, workdir, input_generator, max_niter=25):
        """
        Args:
            workdir:
            input_generator:
            max_niter:
                Maximum number of iterations. A negative value or zero value 
                is equivalent to having an infinite number of iterations.
        """
        super(IterativeWork, self).__init__(workdir)

        self.input_generator = input_generator

        self.max_niter = max_niter

    def next_task(self):
        """
        Generate and register a new task

        Return: task object
        """
        try:
            next_input = next(self.input_generator)
        except StopIteration:
            raise StopIteration

        self.register_input(next_input)
        assert len(self) == self.niter

        return self[-1]

    def submit_tasks(self, *args, **kwargs):
        """
        Run the tasks till self.exit_iteration says to exit or the number of iterations exceeds self.max_niter

        Return dictionary with the final results
        """
        self.niter = 1
        iteration_results = {"converged": False}

        #if kwargs.pop("num_pythreads", 1) != 1

        while True:

            if self.max_niter > 0 and self.niter > self.max_niter: 
                print("niter %d > max_niter %d" % (self.niter, self.max_niter))
                break 

            try:
                task = self.next_task()
            except StopIteration:
                break

            # Start the task and block till completion.
            task.start(*args,**kwargs)

            task.communicate()

            data = self.exit_iteration(*args, **kwargs)

            if data["exit"]:
                break

            self.niter += 1

    @abc.abstractmethod
    def exit_iteration(self, *args, **kwargs):
        """
        Return a dictionary with the results produced at the given iteration.
        The dictionary must contains an entry "converged" that evaluates to 
        True if the iteration should be stopped.
        """

##########################################################################################

def check_conv(values, tol, min_numpts=1, mode="abs", vinf=None):
    """
    Given a list of values and a tolerance tol, returns the leftmost index for which

        abs(value[i] - vinf) < tol if mode == "abs"

        or 

        abs(value[i] - vinf) / vinf < tol if mode == "rel"

    returns -1 if convergence is not achieved. By default, vinf = values[-1]

    Args:
        tol:
            Tolerance
        min_numpts:
            Minimum number of points that must be converged. 
        mode:
            "abs" for absolute convergence, "rel" for relative convergence.
        vinf:
            Used to specify an alternative value instead of values[-1].
    """
    vinf = values[-1] if vinf is None else vinf

    if mode == "abs":
        vdiff = [abs(v - vinf) for v in values]
    elif mode == "rel":
        vdiff = [abs(v - vinf) / vinf for v in values]
    else:
        raise ValueError("Wrong mode %s" % mode)

    numpts = len(vdiff)
    i = -2

    if (numpts > min_numpts) and vdiff[-2] < tol:
        for i in range(numpts-1, -1, -1):
            if vdiff[i] > tol:
                break
        if (numpts - i -1) < min_numpts: i = -2
 
    return i + 1

class PseudoConvergence(Work):

    def __init__(self, workdir, pseudo, ecut_list, 
                 spin_mode     = "polarized", 
                 acell         = 3*(8,), 
                 smearing      = None,
                ):

        super(PseudoConvergence, self).__init__(workdir)

        # Temporary object used to build the input files
        generator = PseudoIterativeConvergence(workdir, pseudo, ecut_list, 
                                               spin_mode     = spin_mode, 
                                               acell         = acell, 
                                               smearing      = smearing,
                                              )
        self.pseudo = pseudo
        if not isinstance(pseudo, Pseudo):
            self.pseudo = Pseudo.from_filename(pseudo)

        self.ecut_list = []
        for ecut in ecut_list:
            input = generator.input_with_ecut(ecut)
            self.ecut_list.append(ecut)
            self.register_input(input)

    def get_results(self, *args, **kwargs):

        # Get the results of the tasks.
        work_results = super(PseudoConvergence, self).get_results(*args, **kwargs)

        etotal = self.read_etotal()

        data = {
            "ecut_list"   : self.ecut_list,
            "etotal"      : list(etotal),
            "pseudo_name" : self.pseudo.name,
            "pseudo_path" : self.pseudo.path,
            #"ecut_low"    : ecut_low,
            #"ecut_normal" : ecut_normal,
            #"ecut_high"   : ecut_high,
        }

        work_results.update(data)

        return work_results

class PseudoIterativeConvergence(IterativeWork):

    def __init__(self, workdir, pseudo, ecut_list_or_slice, 
                 spin_mode     = "polarized", 
                 acell         = 3*(8,), 
                 smearing      = None,
                 max_niter     = 100,
                ):
        """
        Args:
            pseudo:
                string or Pseudo instance
            ecut_list_or_slice:
                List of cutoff energies or slice object (mainly used for infinite iterations). 
            spin_mode:
                Defined how the electronic spin will be treated.
            acell:
                Lengths of the periodic box in Bohr.
            smearing:
                Smearing instance. Default: FemiDirac with T=0.1 eV
        """
        self.pseudo = pseudo
        if not isinstance(pseudo, Pseudo):
            self.pseudo = Pseudo.from_filename(pseudo)

        if smearing is None: 
            smearing = Smearing.FermiDirac(0.1/Ha_eV)

        self._spin_mode = spin_mode
        self._smearing  = smearing
        self._acell     = acell

        if isinstance(ecut_list_or_slice, slice):
            self.ecut_iterator = iterator_from_slice(ecut_list_or_slice)
        else:
            self.ecut_iterator = iter(ecut_list_or_slice)

        # Construct a generator that returns AbinitInput objects.
        def input_generator():
            for ecut in self.ecut_iterator:
                yield self.input_with_ecut(ecut)

        super(PseudoIterativeConvergence, self).__init__(workdir, input_generator(), max_niter=max_niter)

        if not self.isnc:
            raise NotImplementedError("PAW convergence tests are not supported yet")

    def input_with_ecut(self, ecut):
        "Return an Input instance with given cutoff energy ecut"

        # Define System: one atom in a box of lenghts acell. 
        boxed_atom = System.boxed_atom(self.pseudo, acell=self._acell)

        # Gamma-only sampling.
        gamma_only = Kpoints.gamma_only()

        # Default setup for electrons: no smearing.
        electrons = Electrons(spin_mode=self._spin_mode, smearing=self._smearing)

        # Generate inputs with different values of ecut and register the task.
        control = Control(boxed_atom, electrons, gamma_only,
                          ecut        = ecut, 
                          prtwf       = 0,
                          want_forces = False,
                          want_stress = False,
                         )

        return Input(control, boxed_atom, gamma_only, electrons)

    @property
    def ecut_list(self):
        "The list of cutoff energies computed so far"
        return [float(task.input.get_variable("ecut")) for task in self]

    def check_etotal_convergence(self, *args, **kwargs):

        ecut_list = self.ecut_list
        etotal = self.read_etotal()

        num_ene = len(etotal)
        etotal_inf = etotal[-1]

        lines = []
        app = lines.append

        app(" pseudo: %s" % self.pseudo.name)
        app(" idx ecut etotal (et-e_inf) [meV]")
        for idx, (ec, et) in enumerate(zip(ecut_list, etotal)):
            line = "%d %.3f %.3f %.3f" % (idx, ec, et, (et-etotal_inf)* Ha_eV * 1.e+3)
            app(line)

        stream = sys.stdout
        stream.writelines("\n".join(lines)+"\n")

        de_high, de_normal, de_low = 0.05e-3/Ha_eV,  1e-3/Ha_eV, 10e-3/Ha_eV
        
        ihigh   = check_conv(etotal, de_high, min_numpts=1)
        inormal = check_conv(etotal, de_normal)
        ilow    = check_conv(etotal, de_low)

        print("ihigh %d, inormal %d, ilow %d" % (ihigh, inormal, ilow))

        ecut_high, ecut_normal, ecut_low = 3 * (None,)
        exit = (ihigh != -1)
        if exit:
            ecut_low    = self.ecut_list[ilow] 
            ecut_normal = self.ecut_list[inormal] 
            ecut_high   = self.ecut_list[ihigh] 

        data = {
            "exit"        : ihigh != -1,
            "ecut_list"   : ecut_list,
            "etotal"      : list(etotal),
            "ecut_low"    : ecut_low,
            "ecut_normal" : ecut_normal,
            "ecut_high"   : ecut_high,
        }

        return data

    def exit_iteration(self, *args, **kwargs):
        return self.check_etotal_convergence(self, *args, **kwargs)
                                                                                       
    def get_results(self, *args, **kwargs):

        # Get the results of the tasks.
        work_results = super(PseudoIterativeConvergence, self).get_results(*args, **kwargs)

        data = self.check_etotal_convergence()

        work_results.update(data)

        #etotal = self.read_etotal()
        #etotal_dict = {1.0 : etotal}
        #status_list = self.get_status()
        #status = max(status_list)
        #num_errors, num_warnings, num_comments = 0, 0, 0 
        #for task in self:
        #    main, log =  task.read_mainlog_ewc()
        #    num_errors  =  max(num_errors,   main.num_errors)
        #    num_warnings = max(num_warnings, main.num_warnings)
        #    num_comments = max(num_comments, main.num_comments)
        #work_results.update({
        #    "num_errors":   num_errors,
        #    "num_warnings": num_warnings,
        #    "num_comments": num_comments,
        #})

        #with open(self.path_in_workdir("results.json"), "w") as fh:
        #    json.dump(work_results, fh)

        return work_results

##########################################################################################

class BandStructure(Work):

    def __init__(self, workdir, structure, pptable_or_pseudos, scf_ngkpt, nscf_nband, kpath_bounds, 
                 spin_mode = "polarized",
                 smearing  = None,
                 ndivsm    = 20, 
                 dos_ngkpt = None
                ):

        super(BandStructure, self).__init__(workdir)

        scf_input = Input.SCF_groundstate(structure, pptable_or_pseudos, 
                                          ngkpt     = scf_ngkpt, 
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

class Relaxation(Work):

    def __init__(self, workdir, structure, pptable_or_pseudos, ngkpt,
                 spin_mode = "polarized",
                 smearing  = None,
                ):
                                                                                                   
        super(Relaxation, self).__init__(workdir)

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

class DeltaTest(Work):

    def __init__(self, workdir, structure_or_cif, pptable_or_pseudos, kppa,
                 spin_mode = "polarized",
                 smearing  = None,
                 accuracy  = "normal",
                 ecutsm    = 0.05,
                ):

        if isinstance(structure_or_cif, Structure):
            structure = structure_or_cif
        else:
            # Assume CIF file
            structure = read_structure(structure_or_cif)

        self._input_structure = structure

        super(DeltaTest, self).__init__(workdir)

        v0 = structure.volume

        self.volumes = v0 * np.arange(90, 112, 2) / 100.

        for vol in self.volumes:

            new_lattice = structure.lattice.scale(vol)

            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)

            scf_input = Input.SCF_groundstate(new_structure, pptable_or_pseudos, 
                                              kppa      = kppa,
                                              spin_mode = spin_mode,
                                              smearing  = smearing,
                                              # **kwargs
                                              accuracy  = accuracy,
                                              ecutsm    = ecutsm,
                                              prtwf     = 0,
                                             )
            self.register_input(scf_input)

    def get_results(self, *args, **kwargs):
        num_sites = self._input_structure.num_sites

        etotal = Ha2eV(self.read_etotal())

        from .eos import EOS
        eos_fit = EOS.Murnaghan().fit(self.volumes, etotal)

        print(eos_fit)

        eos_fit.plot(show=False, savefig=self.path_in_workdir("eos.pdf"))

        results = {
            "etotal" : list(etotal),
            "volumes": list(self.volumes),
            "natom"  : num_sites,
            "v0"     : eos_fit.v0,
            "b"      : eos_fit.b,
            "bp"     : eos_fit.bp,
        }

        with open(self.path_in_workdir("results.json"), "w") as fh:
            json.dump(results, fh)

        # Write data for the computation of the delta factor
        with open(self.path_in_workdir("deltadata.txt"), "w") as fh:
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")
            for (v, e) in zip(self.volumes, etotal):
                fh.write("%s %s\n" % (v/num_sites, e/num_sites))

##########################################################################################

class PvsV(Work):
    "Calculates P(V)."

    def __init__(self, workdir, structure_or_cif, pptable_or_pseudos, ngkpt,
                 spin_mode = "polarized",
                 smearing  = None,
                 ecutsm    = 0.05,  # 1.0
                 dilatmx   = 1.1,
                 strtargets = None,
                ):
        raise NotImplementedError()

        super(PvsV, self).__init__(workdir)

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

    def get_results(self):
        pressures, volumes = [], []
        for task in self:
            # Open GSR file and read etotal (Hartree)
            gsr_path = task.odata_path_from_ext("GSR") 
                                                            
            with GSR_Reader(gsr_path) as gsr_data:

                pressures.append(gsr_data.get_value("pressure"))
                volumes.append(gsr_data.get_value("volume"))

##########################################################################################

class G0W0(Work):

    def __init__(self, workdir, structure, pptable_or_pseudos, scf_ngkpt, nscf_ngkpt, 
                 ppmodel_or_freqmesh, ecuteps, ecutsigx, nband_screening, nband_sigma,
                 spin_mode = "polarized",
                 smearing  = None,
                ):

        super(G0W0, self).__init__(workdir)

        scf_input = Input.SCF_groundstate(structure, pptable_or_pseudos, 
                                          ngkpt     = scf_ngkpt, 
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
     #def with_contour_deformation()

##########################################################################################

class PPConvergenceFactory(object):
    "Factory object"

    def work_for_pseudo(self, pseudo, ecut_range, dirname=".", spin_mode="polarized", smearing=None):
        """
        Return a Work object from the given pseudopotential.
        """
        workdir = os.path.join(os.path.abspath(dirname), pseudo.name)

        if smearing is None:
            smearing = Smearing.FermiDirac(0.1 / Ha_eV)

        if isinstance(ecut_range, slice):
            work = PseudoIterativeConvergence(workdir, pseudo, ecut_range,
                                              spin_mode  = spin_mode,
                                              smearing   = smearing,
                                              acell      = 3*(8,), 
                                             )

        else:
            work = PseudoConvergence(workdir, pseudo, ecut_range,
                                     spin_mode  = spin_mode,
                                     smearing   = smearing,
                                     acell      = 3*(8,), 
                                    )

        return work

##########################################################################################

class SimpleParallelRunner(object):
    "Execute a list of work objects in parallel with python threads"

    def __init__(self, max_numthreads):
        """
            Args:
                max_numthreads: The maximum number of threads that can be used
        """
        self.max_numthreads = max_numthreads 

    def run_works(self, works, *args, **kwargs):
        "Call the start method of the object contained in the iterable works."

        # Remove num_pythreads from kwargs since it's one of the options 
        # accepted by the start methods.
        kwargs.pop("num_pythreads", None)

        num_works = len(works)

        if num_works == 0:
            # Likely an iterative work, run it with one thread.
            num_works = 1

        num_pythreads = min(num_works, self.max_numthreads)

        print("%s: Executing works with num_pythreads %d" % (self.__class__.__name__, num_pythreads))
        from threading import Thread
        from Queue import Queue

        def worker():
            while True:
                func, args, kwargs = q.get()
                results = func(*args, **kwargs)
                q.task_done()
                                                     
        q = Queue()
        for i in range(num_pythreads):
            t = Thread(target=worker)
            t.setDaemon(True)
            t.start()

        # If works is a Work instance, we have to call setup here
        #if hasattr(works, "setup"):
        #    works.setup(*args, **kwargs)

        for work in works:
            q.put((work.start, args, kwargs))

        # Block until all tasks are done. 
        q.join()  

        work.communicate()

        # If works is a Work instance, we have to call get_results here
        #results = {}
        #if hasattr(works, "get_results"):
        #    results = works.get_results(*args, **kwargs)

        #return results

##########################################################################################
