"""
Classes defining Abinit calculations and workflows
"""
from __future__ import division, print_function

import sys
import os
import os.path
import abc
import collections
import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.design_patterns import Enum
from pymatgen.core.physical_constants import Bohr2Ang, Ang2Bohr, Ha2eV, Ha_eV, Ha2meV
from pymatgen.serializers.json_coders import MSONable, PMGJSONDecoder
from pymatgen.io.smartio import read_structure
from pymatgen.util.filelock import FileLock

from .netcdf import GSR_Reader
from .pseudos import Pseudo, PseudoDatabase, PseudoTable, PseudoExtraInfo, get_abinit_psp_dir
from .abinit_input import Input, Electrons, System, Control, Kpoints
from .task import AbinitTask, TaskDependencies
from .utils import parse_ewc, abinit_output_iscomplete, remove_tree

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
class WorkError(Exception):
    pass

class Work(MSONable):
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
        The default implementation does nothing.
        """

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

        # Create top level directory.
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)

        # Create task workdirs and files.
        for task in self:
            task.build(*args, **kwargs)

    def get_status(self, only_highest_rank=False):
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
                                                                                    
        remove_tree(self.workdir, **kwargs)

    def run(self, *args, **kwargs):

        num_pythreads = int(kwargs.pop("num_pythreads", 1))

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

    def start(self, *args, **kwargs):
        """
        Start the work. Call setup first, then the 
        the tasks are executed. Finally, the teardown method is called.

            Args:
                num_pythreads = Number of python threads to use (defaults to 1)
        """
        # Acquire the lock
        lock = FileLock(self.path_in_workdir("__lock__"))

        try:
            lock.acquire()
        except FileLock.Exception:
            raise

        # Initial setup
        self.setup(*args, **kwargs)

        # Execute 
        results = self.run(*args, **kwargs)

        # Finalize
        self.teardown(*args, results=results)

        # Release the lock
        lock.release()

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

    def __init__(self, workdir, input_generator, max_niter=10):
        """
        Args:
            workdir:
            input_generator:
            max_niter:
                Maximum number of iterations. 
                A negative value is equivalent to having an infinite number of iterations.
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

    def run(self, *args, **kwargs):

        self.niter = 1
        results = {"converged": False}

        while True:
            if self.max_niter > 0 and self.niter > self.max_niter: 
                raise StopIteration("niter %d > max_niter %d" % (self.niter, self.max_niter))

            try:
                task = self.next_task()
            except StopIteration:
                break

            task.build()

            task.run(*args,**kwargs)

            results = self.exit_iteration(*args, **kwargs)

            if results["converged"]:
                break

            self.niter += 1

        return results

    @abc.abstractmethod
    def exit_iteration(self, *args, **kwargs):
        "return exit, data"

##########################################################################################
def iterator_from_slice(s):
    "Constructs an iterator given a slice object"

    start = s.start if s.start is not None else 0
    step  = s.step  if s.step is not None else 1

    if s.stop is None: 
        # Infinite iterators.
        import itertools
        return itertools.count(start=start, step=step)
    else:
        # xrange-like iterator.
        return iter(range(start, s.stop, step))

class PseudoEcutConvergence(IterativeWork):

    def __init__(self, workdir, pseudo, ecut_list_or_slice, 
                 spin_mode     = "polarized", 
                 acell         = 3*(8,), 
                 smearing      = None,
                ):
        """
        Args:
            pseudo:
                string or Pseudo instance
        """
        self.pseudo = pseudo
        if not isinstance(pseudo, Pseudo):
            self.pseudo = Pseudo.from_filename(pseudo)

        self._spin_mode = spin_mode
        self._smearing  = smearing
        self._acell     = acell

        if isinstance(ecut_list_or_slice, slice):
            self.ecut_iterator = iterator_from_slice(ecut_list_or_slice)
        else:
            self.ecut_iterator = iter(ecut_list_or_slice)

        def input_generator():
            for ecut in self.ecut_iterator:
                yield self.input_with_ecut(ecut)

        super(PseudoEcutConvergence, self).__init__(workdir, input_generator(), max_niter=-1)

        if not self.isnc:
            raise NotImplementedError("PAW convergence tests are not supported yet")

    def input_with_ecut(self, ecut):

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
        return [float(task.input.get_variable("ecut")) for task in self]

    def exit_iteration(self, *args, **kwargs):

        ecut_list = self.ecut_list
        etotal = self.read_etotal()

        num_ene = len(etotal)
        etotal_inf = etotal[-1]
                                                                
        #ediff_mev = [(e-etotal_inf)* Ha_eV * 1.e+3 for e in etotal]
        print(" idx ecut, etotal (et-e_inf) [meV]")
        for idx, (ec, et) in enumerate(zip(ecut_list, etotal)):
            print(idx, ec, et, (et-etotal_inf)* Ha_eV * 1.e+3)

        ecut_high, ecut_normal, ecut_low, conv_idx = 4 * (None,)
        strange_data = 0

        de_high, de_normal, de_low = 0.1e-3/Ha_eV,  1e-3/Ha_eV, 1e-3/Ha_eV

        if (num_ene > 3) and (abs(etotal[-2] - etotal_inf) < de_high):

            for i in range(num_ene-2, -1, -1):
                etot  = etotal[i] 
                ediff = etot - etotal_inf
                if ediff < 0.0: strange_data += 1

                if ecut_high is None and ediff > de_high:
                    conv_idx =  i+1
                    ecut_high = ecut_list[i+1]
                                                                                      
                if ecut_normal is None and ediff > de_normal:
                    ecut_normal = ecut_list[i+1]
                                                                                      
                if ecut_low is None and ediff > de_low:
                    ecut_low = ecut_list[i+1]
                                                                                      
            if conv_idx is None or (num_ene - conv_idx) < 2:
                print("Not converged %s " % str(conv_idx))
                strange_data += 1

        results = {
            "converged"  : conv_idx,
            "ecut_list"  : ecut_list,
            "etotal"     : list(etotal),
            "ecut_low"   : ecut_low,
            "ecut_normal": ecut_normal,
            "ecut_high"  : ecut_high,
        }

        return results

    def teardown(self, *args, **kwargs):
        import json

        results = kwargs.pop("results", {})
        print(results)

        if results:
            with open(self.path_in_workdir("results.json"), "w") as fh:
                json.dump(results, fh)

        etotal = self.read_etotal()

        etotal_dict = {1.0 : etotal}

        status_list = self.get_status()

        status = max(status_list)

        #num_warnings = max(s.nwarnings for s in status_list) 
        #num_errors = max(s.nerrors for s in status_list) 
        #num_comments = max(s.comments for s in status_list) 

        num_errors, num_warnings, num_comments = 0, 0, 0 
        for task in self:
            main, log =  task.read_mainlog_ewc()

            num_errors  =  max(num_errors,   main.num_errors)
            num_warnings = max(num_warnings, main.num_warnings)
            num_comments = max(num_comments, main.num_comments)

        with open(self.path_in_workdir("data.txt"), "w") as fh:
            info = "max_errors = %d, max_warnings %d, max_comments = %d\n" % (num_errors, num_warnings, num_comments)
            print(info)
            fh.write(info)
            for (ecut, ene) in zip(self.ecut_list, etotal):
                fh.write(str(ecut) + " " + str(ene) + "\n")

        #print(repr(status))
        #input = ["I really\n", "dislike\n", "XML\n"]

        # TODO handle possible problems with the SCF cycle
        strange_data = 0
        if num_errors != 0 or num_warnings != 0: strange_data = 999

        pp_info = PseudoExtraInfo.from_data(self.ecut_list, etotal_dict, strange_data)
        #pp_info.show_etotal()

        # Append PP_EXTRA_INFO section to the pseudo file
        # provided that the results seem ok. If results seem strange, 
        # the XML section is written on a separate file.
        xml_string = pp_info.toxml()
        #print self.pseudo.basename, xml_string

        with open(self.path_in_workdir("pp_info.xml"), "w") as fh:
            fh.write(xml_string)
        #new = pp_info.from_string(xml_string)
        #print new
        #print new.toxml()

##########################################################################################

class BandStructure(Work):

    def __init__(self, workdir, structure, pptable_or_pseudos, scf_ngkpt, nscf_nband, kpath_bounds, 
                 spin_mode = "polarized",
                 smearing  = None,
                 ndivsm    = 20, 
                 dos_ngkpt = None
                ):

        super(BandStructure, self).__init__(workdir)

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

    def __init__(self, workdir, structure_or_cif, pptable_or_pseudos, ngkpt,
                 spin_mode = "polarized",
                 smearing  = None,
                 accuracy  = "normal",
                 ecutsm    = 0.05,
                 ecut      = None
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

            scf_input = Input.SCF_groundstate(new_structure, pptable_or_pseudos, ngkpt, 
                                              spin_mode = spin_mode,
                                              smearing  = smearing,
                                              # **kwargs
                                              accuracy  = accuracy,
                                              ecutsm    = ecutsm,
                                              ecut      = ecut,
                                             )
            self.register_input(scf_input)

    def run(self, *args, **kwargs):
        pass

    def teardown(self, *args, **kwargs):
        num_sites = self._input_structure.num_sites

        etotal = Ha2eV(self.read_etotal())

        with open(self.path_in_workdir("data.txt"), "w") as fh:
            fh.write("# Volume/natom [Ang^3] Etotal/natom [eV]\n")

            for (v, e) in zip(self.volumes, etotal):
                line = "%s %s\n" % (v/num_sites, e/num_sites)
                fh.write(line)

        from .eos import EOS
        eos_fit = EOS.Murnaghan().fit(self.volumes/num_sites, etotal/num_sites)
        print(eos_fit)
        eos_fit.plot()

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

    def teardown(self):

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
     #def with_contour_deformation()

##########################################################################################

def test_abinitio(structure):
    psp_dir = get_abinit_psp_dir()

    pp_database = PseudoDatabase(dirpath=psp_dir, force_reload=False)

    GGA_HGHK_PPTABLE = pp_database.GGA_HGHK_PPTABLE

    pptable_or_pseudos = GGA_HGHK_PPTABLE
    pseudo = GGA_HGHK_PPTABLE[14][0]

    pptest_wf = PseudoEcutTest("Test_pseudo", pseudo, list(range(10,40,2)))

    #pptest.show_inputs()
    pptest_wf.start()
    sys.exit(1)

    scf_ngkpt = [4,4,4]
    kpath_bounds = [0,0,0, 0.5,0.5,0.5]
    nscf_nband = 333

    bands_wf = BandStructure("Test_bands",structure, pptable_or_pseudos, scf_ngkpt, nscf_nband, kpath_bounds)
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

    g0w0_wf = G0W0("Test_G0W0", structure, pptable_or_pseudos, scf_ngkpt, nscf_ngkpt, 
                            ppmodel_or_freqmesh, ecuteps, ecutsigx, nband_screening, nband_sigma)

    g0w0_wf.show_inputs()

    delta_wf = DeltaTest(structure, pseudo, [2,2,2])
    delta_wf.show_inputs()

##########################################################################################

class EcutTest(object):

    def __init__(self, table, spin_mode):
        works = []

        table_name = table.name

        #table = [p for p in table if p.Z < 5]
        table = [p for p in table if p.Z < 2]

        for pseudo in table:
            ecut_list = list(range(10,31,5))
            #ecut_list = list(range(50,451,50))
            #ecut_list = slice(20, None, 40)

            #element = pseudo.element
            #print(element.Z, element)
            #if element.block not in ["s", "p"]:
            #    ecut_list = range(30,121,1),

            workdir = "PPTEST_" + table_name + "_" + pseudo.basename

            pp_wf = PseudoEcutConvergence(workdir, pseudo, ecut_list,
                        spin_mode  = spin_mode,
                        smearing   = None,
                        acell      = 3*(8,), 
                       )

            #print(pp_wf.show_inputs())
            #sys.exit(1)

            works.append(pp_wf)

        self.works = works

    def build(self, *args, **kwargs):
        for work in self.works:
            work.build(*args, **kwargs)

    def start(self, *args, **kwargs):

        num_pythreads = kwargs.pop("num_pythreads", 1)

        if num_pythreads == 1:
            for work in self.works:
                work.start(*args, **kwargs)
            return

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

        args, kwargs = [], {}
        for work in self.works:
            q.put((work.start, args, kwargs))

        #Block until all tasks are done. 
        q.join()  

##########################################################################################

if __name__ == "__main__":
    test_abinitio()
