__author__ = 'David Cossey, Gary Kedziora'
__copyright__ = 'Approved for public release: AFRL DSRC Case Number 88ABW-2016-0260'
__version__ = '0.5'
__maintainer__ = 'David Cossey'
__email__ = 'dcossey014@gmail.com; gary.kedziora@engilitycorp.com'
__date__ = '9/14/15'

import os
import sys
import six
import glob
import shutil
import string
import logging
import datetime
import traceback

from monty.serialization import loadfn

from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import MinimalVaspInputSet
from fireworks import Firework, FireTaskBase, FWAction, explicit_serialize, Workflow, LaunchPad

from custodian import Custodian
from custodian.vasp.jobs import VaspJob
from custodian_afrl.vasp.validators import VasprunXMLValidator

logger = logging.getLogger(__name__)
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

def load_class(mod, name):
    mod = __import__(mod, globals(), locals(), [name], 0)
    return getattr(mod, name)


class VaspInput:
    """
    This is an interface class to WriteVaspInputTask(FireTaskBase).  

    Required Params:
        -s (obj):               Pymatgen Structure Object used to for
                                creating VASP input files

    Optional Params:
        -vis (Vasp Input Set):  Input Set from pymatgen.io.vasp.sets to use
                                for default input parameters to VASP
        -isp (Input Set param): Custom VASP input parameters to overwrite the
                                default parameters from the VASP Input Set file.
        -custom_params ({dict}): Used to pass custom KPOINT settings for 
                                a VASP calculation.
        -config_file (str):     Location of a configuration file used to overwrite
                                VASP input set parameters.  If none is given, 
                                it will check for vasp_interface_defaults.yaml in 
                                the user's HOME directory.
    """

    def __init__(self,s=None,vis='MinimalVaspInputSet',isp={}, custom_params={}, config_file=None):
        # Set ignore list for ignoring certain __dict__ items in member functions
        self.__dict__['ignore_list'] = ['input','custom_settings',
                'user_incar','vasp_is','ignore_list', 'input_set_params', 
                'kpts_params', 'params']

        if config_file:
            config_dict = loadfn(config_file)
        elif os.path.exists(os.path.join(os.environ['HOME'], 'vasp_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'], 'vasp_interface_defaults.yaml'))
        else:
            config_dict = {}
            self.__dict__['vasp_is'] = vis
            self.__dict__['user_incar'] = isp
            self.__dict__['custom_settings'] = custom_params

            print("\n*********************************************\n"
                    "Configuration File not supplied and\n"
                    "${HOME}/vasp_interface_defaults.yaml does\n"
                    "not exist. You can copy this file from\n"
                    "pymatgen/io/vasp/InterfaceExamples/test_files/\n
                    vasp_interface_defaults.yaml and modify it \n"
                    "for your project's needs.\n"
                    "*********************************************\n\n"
                    "Will continue with MinimalVaspInputSet parameters.\n")
            
        if config_dict:
            self.__dict__['vasp_is'] = config_dict.get('VASP_INPUT_SET', 'MinimalVaspInputSet')
            self.__dict__['user_incar'] = {} if config_dict['INPUT_SET_PARAMS']['user_incar_settings'] == None else config_dict['INPUT_SET_PARAMS']
            self.__dict__['custom_settings'] = config_dict.get('CUSTOM_PARAMS', {}) 
            if isp:
                self.user_incar['user_incar_settings'].update(isp)
            if custom_params:
                self.custom_settings['user_kpts_settings'].update(custom_params)

        self.__dict__['input']=WriteVaspInputTask(structure=s, vasp_input_set=self.vasp_is,
                input_set_params=self.user_incar, custom_params=self.custom_settings)

        self.__dict__['params'] = [u'NELMDL', u'LMUSIC', u'WC', u'LCHIMAG', u'INIWAV', u'APACO', 
                  u'IBRION', u'ICORELEVEL', u'MODEL_GW', u'IRESTART', u'MAGMOM', 
                  u'PARAM2', u'POMASS', u'LRPA', u'ADDGRID', u'OFIELD_KAPPA', 
                  u'MAXMIX', u'EMAX', u'LFXHEG', u'LBERRY', u'ENAUG', u'LHFONE', 
                  u'ISYM', u'NREBOOT', u'MODEL_ALPHA', u'AEXX', u'LTCTE', 
                  u'IDIOT', u'ENINI', u'LELF', u'HFALPHA', u'NBANDS', u'LMAXPAW', 
                  u'LSPIRAL', u'ISMEAR', u'ENCUTGWSOFT', u'EBREAK', u'NBANDSO', 
                  u'IWAVPR', u'LCHARG', u'WEIMIN', u'NEDOS', u'LCOMPAT', 
                  u'EFIELD', u'LVEL', u'LDAU', u'LMAXMIX', u'NOMEGAR', 
                  u'ODDONLYGW', u'ISPIN', u'EVENONLY', u'BMIX_MAG', u'LUSE_VDW',
                  u'NMAXFOCKAE', u'SYMPREC', u'LSYMGRAD', u'LVHAR', u'IBSE', 
                  u'LVTOT', u'LGAUGE', u'TEEND', u'AVECCONST', u'MAGDIPOL', 
                  u'ODDONLY', u'ORBITALMAG', u'MCALPHA', u'LZEROZ', u'SIGMA',
                  u'SAXIS', u'NGZ', u'TURBO', u'NBLOCK', u'TELESCOPE', u'NCORE',
                  u'NRMM', u'IDIPOL', u'HFSCREENC', u'EVENONLYGW', u'LRHFCALC',
                  u'EFERMI', u'NBANDSGWLOW', u'LSCAAWARE', u'OMEGATL', u'EDIFF', 
                  u'LEPSILON', u'FOURORBIT', u'DIPOL', u'LTRIPLET', u'LDIPOL', 
                  u'LSORBIT', u'NUPDOWN', u'PRECFOCK', u'HFSCREEN', u'LTHOMAS', 
                  u'LSCALAPACK', u'NBLK', u'NFREE', u'OFIELD_Q6_FAR', u'ENCUT4O', 
                  u'GGA_COMPAT', u'NELECT', u'LMAXMP2', u'LCORR', u'NBANDSV', 
                  u'EXXOEP', u'AMIX', u'TIME', u'LORBIT', u'RWIGS', u'QSPIRAL', 
                  u'NKREDLFX', u'NKREDLFY', u'NKREDLFZ', u'NKREDLF', u'ALGO', 
                  u'IALGO', u'LPARD', u'PSTRESS', u'PARAM3', u'LMODELHF', 
                  u'PARAM1', u'CSHIFT', u'NGYF', u'LWAVE', u'SELFENERGY', u'NPACO', 
                  u'NELMIN', u'LSPECTRAL', u'ENCUTLF', u'EMIN', u'LMAGBLOCH', 
                  u'LSCALU', u'ICHARG', u'NSIM', u'LHARTREE', u'MODEL_EPS0', 
                  u'I_CONSTRAINED_M', u'GGA', u'OMEGAMIN', u'TEBEG', u'DEPER', 
                  u'NMIN', u'DEG_THRESHOLD', u'NLSPLINE', u'LNICSALL', u'MREMOVE', 
                  u'SCALEE', u'INIMIX', u'AMIN', u'NGZF', u'LMETAGGA', u'NPAR', 
                  u'LUSEW', u'MIXPRE', u'MAXMEM', u'ALDAC', u'L2ORDER', u'NELM', 
                  u'OMEGAGRID', u'PREC', u'TAUPAR', u'NWRITE', u'ALDAX', u'EREF', 
                  u'EDIFFG', u'LHFCALC', u'NBANDSGW', u'SHIFTRED', u'LSUBROT', 
                  u'LFXC', u'KSPACING', u'DIM', u'LFOCKAEDFT', u'NUCIND', u'KPAR',
                  u'OFIELD_Q6_NEAR', u'NKREDX', u'NKREDY', u'NKREDZ', u'KINTER', 
                  u'LASPH', u'VOSKOWN', u'ISIF', u'KGAMMA', u'SMASS', u'ENMAX', 
                  u'LFXCEPS', u'AMIX_MAG', u'OFIELD_K', u'BMIX', u'ISTART', 
                  u'AGGAC', u'NGX', u'LDIAG', u'NGY', u'LOPTICS', u'LASYNC', 
                  u'Zab_VDW', u'AGGAX', u'NSW', u'ENCUTGW', u'ROPT', u'STM', 
                  u'LNABLA', u'NOMEGA', u'SCISSOR', u'POTIM', u'LPLANE', u'LREAL', 
                  u'KBLOCK', u'NGXF', u'DARWINR', u'MAGATOM', u'DARWINV', u'NKRED',
                  u'LORBITALREAL', u'EPSILON', u'LNONCOLLINEAR', u'OMEGAMAX', 
                  u'SYSTEM', u'IMIX', u'MAGPOS', u'LMAXFOCK', u'OFIELD_A', 
                  u'LMONO', u'ANTIRES', u'LTETE', u'LADDER', u'ENCUT', u'ZVAL']

        self.__dict__['kpts_params'] = ['kpts_style', 'kpts', 'kpts_shift']
        self.__dict__['input_set_params'] = ['constrain_total_magmom','sort_structure',
                'ediff_per_atom','potcar_functional','force_gamma','reduce_structure']

        #Set dictionaries in case they do not exsist
        try:
            self.input['input_set_params']['user_incar_settings']
        except:
            self.input['input_set_params']['user_incar_settings'] = {}
        try:
            self.input['custom_params']['user_kpts_settings']
        except:
            self.input['custom_params']['user_kpts_settings'] = {}
        
        
    def __setattr__(self, key, val):
        '''
        Automatically adds Class Object Attributes directly to their
        appropriate dictionary for the WriteVaspInputTask

        Example:
        This will put NEDOS=2048 into the WriteVaspInputTask['input_set_params']\
                    ['user_incar_settings']['NEDOS']=2048
            input = VaspInput(s=s)
            input.NEDOS=2048
        
        This will put KPTS=[12,12,12] into WriteVaspInputTask['custom_params']\
                    ['user_kpts_settings']['kpts']=[12,12,12]
            input.kpts=[12,12,12]
        '''

        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(),val)

    def __setitem__(self, key, val):
        '''
        Automatically adds Class Dictionary Entries directly to their 
        appropriate dictionary for WriteVaspInputTask.  
        
        Example:
        This will put NEDOS=2048 into the WriteVaspInputTask['input_set_params']\
                    ['user_incar_settings']['NEDOS']=2048
            input = VaspInput(s=s)
            input['NEDOS']=2048
        '''
        self.proc_key_val(key.strip(), val.strip()) if isinstance(
                val, six.string_types) else self.proc_key_val(key.strip(),val)

    def proc_key_val(self, key, val):
        '''
        Internal method for sorting Class attributes and Class Dictionary
        entries into their appropriate WriteVaspInputTask Dictionary
        '''

        incar_dict = self.input.get('input_set_params', 
                {}).get('user_incar_settings', {})
        kpts_dict = self.input.get('custom_params', 
                {}).get('user_kpts_settings', {})

        if key not in self.params and key not in self.kpts_params and \
                key not in self.ignore_list and key not in self.input_set_params:
            print "Parameter '{}' not found in any parameter list".format(key)
            print "Please check spelling and try again."
            print("Adding {} to base <VaspInput Obj> Attribute\n\n".format(key))

        if key in self.params:
            incar_dict.update({key: val})
        elif key in self.kpts_params:
            kpts_dict.update({key: val})
        elif key in self.input_set_params:
            self.input['input_set_params'].update({key:val})
        else:
            self.__dict__.update({key:val})


    def write_input(self,dir='.'):
        self.input.write_input_task(dir)

    def get_nelect(self):
        return self.input.get_nelect_task()

    def add_to_allowed_list(self,parameter):
        self.__dict__['params'].append(parameter)


class VaspFirework():
    """
    This is an interface to Fireworks for VaspInput workflows

    Required Params:
        -vasp_task (obj):   FireTask Object used to create a Firework.
                            2 FireTasks are automatically created upon constructing
                            the VaspFirework object from a WriteVaspInputTask: 
                            [<WriteVaspInputTask()>, <RunCustodianTask()>].

    Optional Params:
        -name (str):        Name given to the Firework
    """

    def __init__(self,vasp_task=None,name='vaspfw',handlers=None, handler_params=None, config_file=None):
        self.name = name
        self.handlers=handlers if handlers else []
        self.handler_params=handler_params if handler_params else {}

        if config_file:
            config_dict = loadfn(config_file)
        elif os.path.exists(os.path.join(os.environ['HOME'], 'vasp_interface_defaults.yaml')):
            config_dict = loadfn(os.path.join(os.environ['HOME'], 'vasp_interface_defaults.yaml'))
        else:
            config_dict = {}

        if config_dict:
            self.custodian_opts = config_dict.get('CUSTODIAN_PARAMS', {})
            if self.custodian_opts.get('handlers', []):
                self.handlers.extend(self.custodian_opts.get('handlers', []))
            self.handler_params.update(self.custodian_opts.get('handler_params', {}))

        self.tasks=[vasp_task.input,RunCustodianTask(handlers=self.handlers, 
                handler_params=self.handler_params)] if isinstance(vasp_task, 
                VaspInput) else [vasp_task]
        self.Firework=Firework(self.tasks,name=self.name)
        self.LaunchPad=LaunchPad.from_file(os.path.join(os.environ["HOME"], ".fireworks", "my_launchpad.yaml"))

    def add_task(self, task):
        '''
        Function used to add another FireTask to the Firework.  If
        given task is a VaspInput Object, it will create the WriteInputTask
        as well as add another RunCustodianTask() to the FireTask list.
        '''

        if isinstance(task, VaspInput):
            self.tasks.extend([task.input, 
                    RunCustodianTask(handlers=self.handlers,
                            handler_params=self.handler_params)])
        else:
            self.tasks.append(task)
        self.Firework=Firework(self.tasks,name=self.name)

    def to_file(self,file_name):
        self.Firework.to_file(file_name)

    def add_fw_to_launchpad(self):
        self.Workflow=Workflow([self.Firework])
        self.LaunchPad.add_wf(self.Workflow)

    def copy_files_from_previous(self, *args, **kwargs):
        '''
        Function used to automatically insert VASPTransferTask 
        between WriteVaspInputTask() and RunCustodianTask() in the 
        FireTask list.  This is used to copy files from a previous 
        VASP Firework Calculation or from a known directory.  
        Known directories require the full path of the directory.
        This function should not be used to copy files from a 
        previous Firetask within the same Firework.

        Parameters:
            -mode (str):    Mode to use in transferring files:
                            move, mv, copy, cp, copy2, copytree, copyfile
                            Default:  'copy'
            -ignore_errors (bool):  Whether to quit upon error.
                            Default: True

            -dir (str):     Directory to transfer files from if not copying
                            directly from a previous Firework
                            Default: None

        Example:
            fw2.copy_files_from_previous('file1', 'file2', ... , 'filen',  mode='copy', ignore_errors=False)

        Example to copy files from directory other than previous Firework:
            fw2.copy_files_from_previous('file1', 'file2', ... , 'filen',
                                        dir='/path/to/dir', mode='copy', ignore_errors=False)
        '''

        mode = kwargs.get('mode', 'copy')
        ignore_errors = kwargs.get('ignore_errors', True)
        files = args[0] if isinstance(args[0], list) else [k for k in args]
        self.dir = kwargs.get('dir', None)

        if self.dir:
            self.transfer=VaspTransferTask(files=files, dir=self.dir, 
                                mode=mode, ignore_errors=ignore_errors)
        else:
            self.transfer=VaspTransferTask(files=files, mode=mode, ignore_errors=ignore_errors)

        if isinstance(self.tasks[-1], RunCustodianTask):
            self.tasks.pop(-1)
            self.add_task(self.transfer)
            self.add_task(RunCustodianTask(handlers=self.handlers, 
                                handler_params=self.handler_params))
        else:
            self.add_task(self.transfer)

    def add_handler(self, handler, **kwargs):
        '''
        Member function to add handler and handler options to 
        all Fireworks in the VaspFirework.tasks list.

        Example (assuming fw1 is the firework object):
        fw1.add_handler('WalltimeHandler', wall_time=3600, buffer_time=300)
        fw1.add_handler('FrozenJobErrorHandler')
        '''
        for i,j in enumerate(self.Firework.tasks):
            if isinstance(j, RunCustodianTask):
                cust_handlers = j.get('handlers', [])
                cust_params = j.get('handler_params', {})
                if handler not in cust_handlers:
                    j['handlers'].append(handler)
                    j['handler_params'].update(kwargs)
        self.Firework.spec['_tasks']=[t.to_dict() for t in
                self.Firework.tasks] #re-initialize FW.spec


class VaspWorkflow():
    """
    A VaspWorkflow encapsulates multiple VaspFirework objects into a Single Workflow.  
    If the kwarg "dependency" is not set, it will create a Sequential Workflow where the next 
    Firework in the Workflow will not start before the currently running Firework in the 
    Workflow completes.

    Parameters:
        -args (obj):        List of VaspFirework objects
        -deps_dict {dict}:  Specifies the dependency of the VaspInput objects given.
                            If no dependency is given, Firworks are assumed to be
                            sequentially dependent.
        -name (str):        Name to be given to the Workflow


    Example:
        VaspWorkflow(FW1, FW2, FW3, FW4, deps_dict={FW1: [FW2, FW3], FW2: [FW4], FW3: [FW4]}, name='Example WF')
        
        This will create a Workflow containing the 4 given VaspFirework objects
        with a Workflow name of 'Example WF' and the given dependencies.
        Dependency Dictionary Explanation:
            In the above example, FW2 and FW3 will not start before FW1 is complete.
            Likewise, FW4 depends on the completion of FW2 and FW3 before starting.
    """

    def __init__(self, *args, **kwargs):
        '''
        :param args:       (VaspFirework objects) objects to create Workflow from.  No limit
                           on the amount of VaspInput objects to be given.  Entered as just
                           comma separated objects passed to class.
        :param deps_dict:  (dict) specifies the dependency of the VaspInput objects given.  
                           If no dependency is given, Firworks are assumed to be 
                           sequentially dependent.
        :param name        (str) Name given to Workflow
        '''
        self.fws = []
        self.name = kwargs.get('name', 'Sequential WF')
        self.deps_dict = kwargs.get('deps_dict', {})
        self.dependency = {}
        if self.deps_dict:
            for i in self.deps_dict.keys():
                fw_deps = []
                for j in self.deps_dict[i]:
                    fw_deps.append(j.Firework)                    
                self.dependency[i.Firework]=fw_deps
        self.deps = True if self.dependency else False
        for id, fw_task in enumerate(args):
            self.fws.append(fw_task.Firework)
            if not self.deps and id != 0:
                self.dependency[self.fws[id-1]]=[fw_task.Firework]
        self.wf = Workflow(self.fws, self.dependency, name=self.name)
        self.LaunchPad=LaunchPad.from_file(os.path.join(os.environ["HOME"], ".fireworks", "my_launchpad.yaml"))


    def add_fw(self, fw_task, deps=None):
        self.fws.append(fw_task.Firework)
        if deps:
            for i in deps.keys():
                fw_deps = []
                for j in deps[i]:
                    fw_deps.append(j.Firework)
                self.dependency[i.Firework]=fw_deps
        else:
            id = len(self.fws) - 2
            self.dependency[self.fws[id]]=[fw_task.Firework]
        self.wf=Workflow(self.fws, self.dependency, name=self.name)


    def to_file(self,filename):
        self.wf.to_file(filename)


    def add_wf_to_launchpad(self):
        self.LaunchPad.add_wf(self.Workflow)


@explicit_serialize
class WriteVaspInputTask(FireTaskBase):
    """ 
    Writes VASP Input files.

    Required params:
        structure (dict): An input structure in pymatgen's Structure.to_dict
            format.
        vasp_input_set (str): A string name for the VASP input set. E.g.,
            "MPVaspInputSet" or "MITVaspInputSet".

    Optional params:
        input_set_params (dict): If the input set requires some additional
            parameters, specify them using input_set_params. Note that if you
            want to change the user_incar_settings, you need to provide
            {"input_set_params": {"user_incar_settings": ...}}. This is
            because input_set_params caters to not just user_incar_settings,
            but all other potential args to the input set.
    """
    required_params = ["structure", "vasp_input_set"]
    optional_params = ["input_set_params","custom_params"]


    def run_task(self, fw_spec):

        prev_dir = fw_spec.get('PREV_DIR', None)
        self.custom_params = self.get('custom_params', None)   

        if isinstance(self["structure"], Structure):
            s = self["structure"]
        elif isinstance(self["structure"], dict):
            s = Structure.from_dict(self["structure"])
        else:
            s = Structure.from_file(os.path.join(prev_dir, self["structure"]))


        vis = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])(
                         **self.get("input_set_params", {}))
        vis.write_input(s, ".")


        # Write Custom KPOINTS settings if necessary
        ksettings = self.custom_params.get('user_kpts_settings', None) if isinstance(
                self.custom_params, dict) else None
        if ksettings:
            style = ksettings.get('kpts_style', 'Gamma')
            kpoints = ksettings.get('kpts', [16,16,16])
            shift = ksettings.get('kpts_shift', [0,0,0])
            k = Kpoints(kpts=[kpoints], kpts_shift=shift)
            k.style = style
            k.write_file("KPOINTS")
            

    def write_input_task(self, dir='.'):
        self.custom_params = self.get('custom_params', None)

        if isinstance(self["structure"], Structure):
            s = self["structure"]
        elif isinstance(self["structure"], dict):
            s = Structure.from_dict(self["structure"])
        else:
            s = Structure.from_file(self["structure"])


        if os.environ.get('VASP_PSP_DIR') == None:
            print "VASP_PSP_DIR not set.  Checking User's HOME directory for VASP potentials."
            if os.path.exists(os.path.join(os.environ.get('HOME'), 'Potentials')):
                os.environ['VASP_PSP_DIR'] = os.path.join(os.environ.get('HOME'), 'Potentials')
            else:
                print "VASP Potentials not found!"
                print "Please copy the Potentials Folder"
                print "from VASP into your HOME directory"
                sys.exit()

        vis = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])(
                         **self.get("input_set_params", {}))
        vis.write_input(s, dir)

        # Write Custom KPOINTS settings if necessary
        ksettings = self.custom_params.get('user_kpts_settings', None) if isinstance(
                self.custom_params, dict) else None
        if ksettings:
            style = ksettings.get('kpts_style', 'Gamma')
            kpoints = ksettings.get('kpts', [16,16,16])
            shift = ksettings.get('kpts_shift', [0,0,0])
            k = Kpoints(kpts=[kpoints], kpts_shift=shift)
            k.style = style
            filename = os.path.join(dir, 'KPOINTS')
            k.write_file(filename)
        
        print "Wrote VASP input files to '{}'".format(dir)

    def get_nelect_task(self):
        
        if isinstance(self["structure"], Structure):
            s = self["structure"]
        elif isinstance(self["structure"], dict):
            s = Structure.from_dict(self["structure"])
        else:
            s = Structure.from_file(self["structure"])


        if os.environ.get('VASP_PSP_DIR') == None:
            if os.path.exists(os.path.join(os.environ.get('HOME'), 'Potentials')):
                os.environ['VASP_PSP_DIR'] = os.path.join(os.environ.get('HOME'), 'Potentials')
            else:
                print "VASP_PSP_DIR not set.  Check User's HOME directory for VASP potentials."
                print "Potentials not found!"
                sys.exit()

        vis = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])(
                         **self.get("input_set_params", {}))
        return vis.get_nelect(s)

@explicit_serialize
class VaspTransferTask(FireTaskBase):
    """
    A FireTask to Transfer files from a preivous VASP run in the same Workflow.
    Required params:
        - mode: (str) - move, mv, copy, cp, copy2, copytree, copyfile
        - files: [list] - list of files to pull from previous Firework in the Workflow

    Optional params:
        - ignore_errors (bool) - Whether or not to raise Errors and terminate task.
        - dir (str)  - Directory from which to copy files from, otherwise PREV_DIR
                       is used.
    """
    required_params = ['mode', 'files']
    optional_params = ['ignore_errors', 'dir']
    
    fn_list = {
        "move": shutil.move,
        "mv": shutil.move,
        "copy": shutil.copy,
        "cp": shutil.copy,
        "copy2": shutil.copy2,
        "copytree": shutil.copytree,
        "copyfile": shutil.copyfile,        
    }
    
    def run_task(self, fw_spec):
        prev_dir = fw_spec.get('PREV_DIR', self.get('dir', None))
        if not prev_dir or not os.path.exists(prev_dir):
            raise RuntimeError(
                    'Source Directory not found or was not specified\n'
                    'Looked for: {}'.format(prev_dir))
        ignore_errors = self.get('ignore_errors', False)
        mode = self.get('mode', 'copy')

        for f in self['files']:
            try:
                if '*' in f:
                    file_list = glob.glob(os.path.join(prev_dir,f))
                    for src in file_list:
                        VaspTransferTask.fn_list[mode](src, '.')
                else:
                    src = os.path.join(prev_dir, f)
                    VaspTransferTask.fn_list[mode](src, '.')
            except:
                if not ignore_errors:
                    traceback.print_exc()
                    raise ValueError(
                            "There was an error performing operation {} from {} "
                            "to '.'".format(mode, self["files"]))


@explicit_serialize
class RunCustodianTask(FireTaskBase):
    """
    Runs VASP using Custodian.

    Optional Params:
        handlers ([str]): List of error handler names to use. See custodian
        .vasp.handlers for list of applicable handlers. Note that
        default args are assumed for all handlers for simplicity. A
        special option is "all", which simply uses a set of common
        handlers for relaxation jobs.
    """

    # Insert way to Add Handler Options and 
    def run_task(self, fw_spec):
        FORMAT = '%(asctime)s  %(levelname)s:  %(pathname)s\n\t%(module)s.%(funcName)s:  %(message)s'
        logging.basicConfig(format=FORMAT, level=logging.INFO, filename="run.log")
        job = VaspJob(["vasp"],default_vasp_input_set=MinimalVaspInputSet(),auto_npar=False)

        #Set up Handlers to be used
        hnames = self.get('handlers')
        logging.info("handler names:  {}".format(hnames))
        handlers = [load_class("custodian.vasp.handlers", n)(**self.get('handler_params', \
                                {})) for n in hnames]

        c = Custodian(handlers=handlers, validators=[VasprunXMLValidator()], jobs=[job])
        output = c.run()
        return FWAction(stored_data=output, mod_spec=[{'_set': {'PREV_DIR': os.getcwd()}}])

