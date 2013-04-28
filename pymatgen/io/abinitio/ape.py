#!/usr/bin/env python
from __future__ import division, print_function

import os
import os.path
import sys
import time
import collections
import shutil
import json
import numpy as np

from pprint import pprint
from scipy.interpolate import UnivariateSpline
from subprocess import Popen, PIPE

try:
    #from pymatgen.io.abinitio.configurations import configurations as dft_neutral_confs
    from .configurations import configurations as dft_neutral_confs
except:
    head, x = os.path.split(os.path.abspath(__file__))
    sys.path.insert(0, head)
    from configurations import configurations as dft_neutral_confs

__version__ = "0.1"
__status__ = "Development"
__date__ = "$April 26, 2013M$"

##########################################################################################
# Helper functions

class ApeError(Exception):
    "Base class for APE Exceptions"

def parse_orbital(orbstr):
    import re
    m = re.match("(\d+)([spdfg]+)(\d+)", orbstr)

    if m:
        return int(m.group(1)), m.group(2), float(m.group(3))

    raise ValueError("Don't know how to interpret %s" % orbstr)


def ape_read_waves(dirname):
    "Read the APE radial wavefunctions located in directory dirname"
    waves = {}
    for filename in os.listdir(dirname):
        if not filename.startswith("wf-"): 
            continue

        path = os.path.join(dirname, filename)
        tokens = filename.split("-")
        assert len(tokens) == 2
        state = tokens[1]
        # TODO check spin and spinor
        # Load [r, Radial wavefunctions and first derivative] in data.
        data = np.loadtxt(path)

        waves[state] = RadialWaveFunction(state, state, data[:,0], data[:,1])
    return waves
                                                                          
def ape_read_potentials(dirname):
    "Read the APE radial potentials located in directory dirname."
    pots = {}
    for filename in os.listdir(dirname):
        if not filename.startswith("v_"):
            continue
        path = os.path.join(dirname, filename)
                                                                          
        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
                                                                       
        pots[filename] = RadialFunction(filename, data[:,0], data[:,1])

    return pots
                                                                          
def ape_read_densities(dirname):
    "Read APE AE densities and tau located in directory dirname."
    dens = {}
    for filename in os.listdir(dirname):
        if filename not in ["density", "tau",]:
            continue
        path = os.path.join(dirname, filename)

        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
        dens[filename] = RadialFunction(filename, data[:,0], data[:,1])

    return dens

def ape_read_logders(dirname):
    "Read the APE log derivatives located in directory dirname."
    ae_logders, pp_logders = {}, {}
    for filename in os.listdir(dirname):
        if not filename.startswith("ld-"):
            continue
        path = os.path.join(dirname, filename)

        tokens = filename.split("-")
        assert len(tokens) == 2
        l_name = tokens[1]

        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
        ae_logders[l_name] = RadialFunction(l_name, data[:,0], data[:,1])
        pp_logders[l_name] = RadialFunction(l_name, data[:,0], data[:,2])

    return ae_logders, pp_logders

##########################################################################################

_char2l = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
    "i": 6,
}

def _asl(obj):
    try:
        return _char2l[obj]
    except KeyError:
        return int(obj)

class QState(collections.namedtuple("QState", "n, l, occ, eig, j, s")):
    "Quantum numbers, occupancies and eigenvalue of the atomic orbital."
    # TODO
    #Spin +1, -1 or 1,2 or 0,1?

    def __new__(cls, n, l, occ, eig=None, j=None, s=None):
        "Intercepts super.__new__ adding type conversion and default values"
        return super(QState, cls).__new__(cls, n, _asl(l), occ, eig=eig, j=j, s=s)

    #@classmethod
    #def asqstate(cls, obj):
    #    if isinstance(obj, cls):
    #        return obj
    #    #if isinstance(obj, str):
    #    raise ValueError("Dont' know how to cast %s: %s" % (obj.__class__.__name__, obj))

    #def __str__(self):

    # Rich comparison support. 
    # Note that the ordering is based on the quantum numbers and not on energies!
    #def __gt__(self, other):
    #def __lt__(self, other):

    def to_apeinput(self):
        if self.s is None:
            string = "%s | %s | %s " % (self.n, self.l, self.occ)
        elif self.s == 2:
            raise NotImplementedError("Spin")
        elif self.j is not None:
            raise NotImplementedError("Spinor")
        else:
            raise ValueError("Don't know how to convert %s" % self)

        return [string,]

##########################################################################################

class AtomicConfiguration(object):
    "Atomic configuration defining the AE atom."

    def __init__(self, Z, states):
        self.Z = Z
        self.states = states

    @classmethod
    def from_symbol(cls, symbol):
        """
        symbol: str or int
            Can be a chemical symbol (str) or an atomic number (int).
        """
        # Table
        # 'element': (atomic number, (n, l, occ, energy))
        # Eg 'Be': (4, [(1, 0, 2, -3.856411), (2, 0, 2, -0.20574400000000001)]),

        Z, items = dft_neutral_confs[symbol]
        states = [QState(n=t[0], l=t[1], occ=t[2]) for t in items]
        return cls(Z, states)

    def __str__(self):
        lines = ["%s: " % self.Z]
        lines += [str(state) for state in self]
        return "\n".join(lines)

    def __iter__(self):
        return self.states.__iter__()

    @property
    def spin_mode(self):
        """
        unpolarized: Spin-unpolarized calculation.
        polarized: Spin-polarized calculation.
        """
        for state in self:
            # FIXME
            if state.s is not None and state.s == 2:
                return "polarized"
        return "unpolarized"

    #@property
    #def title(self):
    #    return self.symbol

    #@property
    #def echarge(self):
    #    echarge = 0.0
    #    for state in self:
    #        echarge -= (2 * state.l + 1) * state.occ
    #    return echarge

    @property
    def isneutral(self):
        "True if neutral configuration"
        return (self.echarge - self.Z) < 1.e-10

    def add_state(self, state):
        #self._push(asqstate(state))
        self._push(state)

    #def add_empty_state(self, state):
    #    self._push(asqstate(state))

    def _push(self, state):
        # TODO check that ordering in the input does not matter!
        if state in self.states:
            raise ValueError("state %s is already in self" % state)
        self.states.append(state)

    def _pop(self, state):
        try:
            self.states.remove(state)
        except ValueError:
            raise

    #def add(self, n, l, df=+1, s=None):
    #    "Add (remove) electrons."

    #@property
    #def num_nodes(self):
    #    "Exact number of nodes computed from the quantum numbers."

    def to_apeinput(self):
        lines  = ["NuclearCharge = %s " % self.Z]
        lines += ["SpinMode = %s" % self.spin_mode]
        lines += ["%Orbitals"]
        for state in self:
            lines += state.to_apeinput()
        lines += ["%"]
        return lines

#class SpinPolarizedConfiguration(AeAtom):
#class DiracConfiguration(AeAtom):

##########################################################################################

class RadialFunction(object):
    "A RadialFunction has a radial mesh and values defined on this mesh."

    def __init__(self, name, rmesh, values):
        """
        Args:
            name:
                Name of the function (string).
            rmesh:
                Points of the radial mesh.
            values:
                Values of the function on the radial mesh
        """
        self.name = name 
        self.rmesh = np.ascontiguousarray(rmesh)
        self.values = np.ascontiguousarray(values)

        assert len(self.rmesh) == len(self.values)

    def __len__(self):
        return len(self.values)

    def __iter__(self):
        "Iterate over (rpoint, value)"
        return zip(self.rmesh, self.values)

    def __getitem__(self, ridx):
        return self.rmesh[ridx], self.values[ridx]

    def __str__(self):
        return "<%s, name = %s>" % (self.__class__.__name__, self.name)

    #def __add__(self):
    #def __sub__(self):
    #def __mul__(self):

    @property
    def rmax(self):
        "Outermost point of the radial mesh"
        return self.rmesh[-1]

    @property
    def rsize(self):
        "Size of the radial mesh"
        return len(self.rmesh)

    @property
    def minmax_ridx(self): 
        """
        Returns the indices of the values in a list with the maximum and minimum value.
        """
        minimum = min(enumerate(self.values), key=lambda s: s[1]) 
        maximum = max(enumerate(self.values), key=lambda s: s[1]) 
        return minimum[0], maximum[0]

    @property
    def inodes(self):
        inodes = []
        for i in range(len(self.values)-1):
            if self.values[i] * self.values[i+1] <= 0:
                inodes.append(i)
        return inodes

    @property
    def spline(self):
        "Cubic spline"
        try:
            return self._spline
        except AttributeError:
            self._spline = UnivariateSpline(self.rmesh, self.values, s=0)
            return self._spline

    @property
    def roots(self):
        "Return the zeros of the spline."
        return self.spline.roots()

    def derivatives(self, r):
        "Return all derivatives of the spline at the point r."
        return self.spline.derivatives(r)

    def integral(self, a=None, b=None):
        """
        Return definite integral of the spline of (r**2 values**2) between two given points a and b
            Args:
                a = rmesh[0] if a is None
                b = rmesh[-1] if a is None
        """
        a = self.rmesh[0] if a is None else a
        b = self.rmesh[-1] if b is None else b
        r2v2_spline = UnivariateSpline(self.rmesh, (self.rmesh * self.values) ** 2, s=0)
        return r2v2_spline.integral(a, b)

    #@property
    #def peaks(self):

    def ir_small(self, abs_tol=0.01):
        """
        Returns the rightmost index where the abs value of the wf becomes greater than abs_tol 

        .. warning:
            Assumes that self.values are tending to zero for r --> infinity.
        """
        for i in range(len(self.rmesh)-1,-1,-1):
            if abs(self.values[i]) > abs_tol:
                break
        return i

##########################################################################################

class RadialWaveFunction(RadialFunction):
    TOL_BOUND = 1.e-10

    def __init__(self, state, name, rmesh, values):
        super(RadialWaveFunction, self).__init__(name, rmesh, values)
        self.state = state

    @property
    def isbound(self):
        "True if self is a bound state."
        back = min(10, len(self))
        return all(abs(self.values[-back:]) < self.TOL_BOUND)

#class DiracWaveFunction(RadialFunction):

##########################################################################################

class ApeRadialMesh(dict):
    """
    APE uses a mesh to store functions. In order to change the mesh parameters you can use the following options:
    """
    # Supported variables
    _KEYS = [
        "MeshType",
        "MeshStartingPoint",
        "MeshOutmostPoint",
        "MeshNumberOfPoints", 
        "MeshDerivMethod",
        "MeshFiniteDiffOrder",
    ]

    def __init__(self, **kwargs):
        super(ApeRadialMesh, self).__init__(**kwargs)

        for k in self:
            if k not in self._KEYS:
                raise ValueError("%s is not a registered key" % k)

    def to_apeinput(self):
        return["%s = %s" % kv for kv in self.items()]

##########################################################################################

class AeSolver(object):
    """
    An AeSolver receives an AtomicConfiguration for the AE atom and solves the AE problem
    using the algorithm and the parameters passed to the constructor.
    """
    ape_exe = "ape"

    def __init__(self, ae_conf, workdir=None, xcfunctional='lda_x+lda_c_pw', 
                 wave_equation="scalar_rel", ape_radmesh=None, ape_control=None, ape_verbose=30):
        """
        XCFunctional (integer, lda_x+lda_c_pw)

        WaveEquation (integer, schrodinger)

        When performing atomic calculations APE can solve either the Kohn-Sham equations, 
        the Dirac-Kohn-Sham equations or the scalar-relativistic Khon-Sham equations. Valid options are:

        - schrodinger: Kohn-Sham equations.
        - dirac: Dirac-Kohn-Sham equations.
        - scalar_rel: scalar-relativistic Kohn-Sham equations.
        """
        self.ae_conf = ae_conf
        self.xcfunctional = xcfunctional
        self.wave_equation = wave_equation
        self.ape_radmesh = ape_radmesh if ape_radmesh else ApeRadialMesh()
        self.ape_control = ape_control if ape_control else ApeControl()
        self.ape_verbose = ape_verbose

        self.workdir = os.path.abspath("helloape")
        #if workdir is not None:
        #    self.workdir = os.path.abspath(workdir)

    @property
    def ape_output(self):
        return os.path.join(self.workdir, "ape.out")

    @property
    def ape_input(self):
        return os.path.join(self.workdir, "ape.in")

    @property
    def ape_stderr(self):
        return os.path.join(self.workdir, "ape.stderr")

    @property
    def aedir(self):
        "Absolute path to the ae directory containing the all-electron results." 
        return os.path.join(self.workdir, "ae")

    # Lazy properties.
    @property
    def ae_waves(self):
        try:
            return self._ae_waves
        except AttributeError:
            self._ae_waves = ape_read_waves(self.aedir)
            return self._ae_waves

    @property
    def ae_potentials(self):
        try:
            return self._ae_potentials
        except AttributeError:
            self._ae_potentials = ape_read_potentials(self.aedir)
            return self._ae_potentials

    @property
    def ae_densities(self):
        try:
            return self._ae_densities
        except AttributeError:
            self._ae_densities = ape_read_densities(self.aedir)
            return self._ae_densities

    def to_apeinput(self):
        lines = ["# Generalities"]
        #["Title = %s " % self.ae_conf.title]
        lines += ["CalculationMode = ae"]
        lines += ["XCFunctional = %s" % self.xcfunctional]
        lines += ["Verbose = %s" % self.ape_verbose]

        lines += ["\n# Hamiltonian"]
        lines += ["WaveEquation = %s" % self.wave_equation]

        lines += ["\n# Atomic configuration."]
        lines.extend(self.ae_conf.to_apeinput())

        lines += ["\n# Radial mesh."]
        lines += self.ape_radmesh.to_apeinput()

        lines += ["\n# APE Control"]
        lines += self.ape_control.to_apeinput()
        return lines

    def show_input(self, stream=sys.stdout):
        lines  = ["AE INPUT".center(80,"*")]
        lines += self.to_apeinput()
        lines += ["END AE INPUT".center(80,"*")]
        stream.writelines("\n".join(lines)+"\n")

    def solve(self):
        assert not hasattr(self, "exit_code")

        print("Solving the all-electron problem ...")
        start = time.time()

        # Write the input file.
        os.makedirs(self.workdir)
        with open(self.ape_input, "w") as fh:
            fh.writelines("\n".join(self.to_apeinput()))

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.ape_exe, self.ape_input, self.ape_output, self.ape_stderr)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        self.exit_code = process.wait()

        print("Run completed in %.1f s" % (time.time() - start))

        if self.exit_code != 0:
            with open(self.ape_stderr, "r") as stderr:
                raise ApeError("APE returned exit_code %s\n stderr: %s" % (self.exit_code, stderr.read()))

        with open(self.ape_stderr, "r") as stderr:
            error = stderr.read()
            if error:
                raise ApeError("APE exit_code %s\n stderr: %s" % (self.exit_code, error))

        # TODO Analize the output for possible warnings.
        #self.validate()

    def plot_waves(self, show=True, savefig=None): 
        _plot_waves(self.ae_waves, show=show, savefig=savefig)

##########################################################################################

class ApeControl(dict):
    """
    The self consistent field procedure will stop when one of the convergence criteria is fulfilled. 
    At each iteration the new guess potential is built mixing the input and output potentials.
    """
    # Supported variables
    _KEYS = [
        # SCF
        "MaximumIter", 
        "ConvAbsDens",
        "ConvRelDens",
        "ConvAbsEnergy",
        "ConvRelEnergy",
        "SmearingFunction",   # This one seems a very delicate point
        "MixingScheme",
        "Mixing",
        "MixNumberSteps",
        "MaximumIter",
        # Eigensolver
        "EigenSolverTolerance",
        "ODEIntTolerance",
        "ODESteppingFunction",
        "ODEMaxSteps",
        "EigenSolverTolerance",
    ]

    def __init__(self, **kwargs):
        super(ApeControl, self).__init__(**kwargs)
                                                                   
        for k in self:
            if k not in self._KEYS:
                raise ValueError("%s is not a registered key" % k)
                                                                   
    def to_apeinput(self):
        return["%s = %s" % kv for kv in self.items()]

##########################################################################################

class ApePPComponents(object):

    @classmethod
    def from_strings(cls, *strings):
        """
        Instanciate the object from a list of strings 
        Example: "3s:1.2:tm"
        """
        states, core_radii, schemes = [], [], []

        occ0 = 0.0
        for s in strings:
            try:
                tokens = s.split(":")
                assert len(tokens) == 3
                assert len(tokens[0]) == 2
                n = tokens[0][0]
                l = tokens[0][1]
                rc = float(tokens[1])
                scheme = tokens[2]

                # Add them to the lists.
                states.append(QState(n, l, occ0))
                core_radii.append(rc)
                schemes.append(scheme)

            except:
                raise ValueError("Malformatted string %s" % s)

        return cls(states, core_radii, schemes)

    def __init__(self, states, core_radii, schemes):
        self.states = states
        self.core_radii = core_radii
        self.schemes = schemes

    def to_apeinput(self):
        lines = ["%PPComponents"]
        for (state, rcore, scheme) in zip(self.states, self.core_radii, self.schemes):
            lines += [" %s | %s | %s | %s " % (state.n, state.l, rcore, scheme)]
        lines += ["%"]
        return lines

##########################################################################################

class PseudoGenerator(object):
    """
    A PseudoGenerator uses the data produced by the AeSolver to construct 
    a pseudopotential using the parameters passed to the constructor.
    """
    ape_exe = "ape"

    def __init__(self, workdir, ae_solver, pp_components, core_correction=0, llocal=-1, pptest_orbitals=None, ppformat="abinit6"):

        self.workdir = os.path.abspath(workdir)
        self.ae_solver = ae_solver
        self.pp_components = pp_components

        self.core_correction = core_correction 
        self.llocal = llocal
        self.pptest_orbitals = pptest_orbitals
        self.ppformat = ppformat

    @property
    def ae_conf(self):
        return self.ae_solver.ae_conf

    @property
    def ape_input(self):
        return os.path.join(self.workdir, "ape.in")

    @property
    def ape_output(self):
        return os.path.join(self.workdir, "ape.out")
                                                        
    @property
    def ape_stderr(self):
        return os.path.join(self.workdir, "ape.stderr")

    @property
    def output(self):
        "Output file generate by APE (List of strings)"
        if not hasattr(self, "_output"):
            with open(self.ape_output, "r") as fh:
                self._output = fh.readlines()
        return self._output[:]

    @property
    def ppdir(self):
        "Absolute path to the ae directory containing the all-electron results." 
        return os.path.join(self.workdir, "pp")

    @property
    def testsdir(self):
        "Absolute path to the directory containing test results." 
        return os.path.join(self.workdir, "tests")

    @property
    def ae_waves(self):
        return self.ae_solver.ae_waves

    # Lazy properties.
    @property
    def pp_waves(self):
        try:
            return self._pp_waves
        except AttributeError:
            self._pp_waves = ape_read_waves(self.ppdir)
            return self._pp_waves
                                                       
    @property
    def pp_potentials(self):
        try:
            return self._pp_potentials
        except AttributeError:
            self._pp_potentials = ape_read_potentials(self.ppdir)
            return self._pp_potentials
                                                       
    @property
    def pp_densities(self):
        try:
            return self._pp_densities
        except AttributeError:
            self._pp_densities = ape_read_densities(self.ppdir)
            return self._pp_densities

    @property
    def ae_logders(self):
        try:
            return self._ae_logders
        except AttributeError:
            self._ae_logders, self._pp_logders = ape_read_logders(self.testsdir)
            return self._ae_logders

    @property
    def pp_logders(self):
        try:
            return self._pp_logders
        except AttributeError:
            self._ae_logders, self._pp_logders = ape_read_logders(self.testsdir)
            return self._pp_logders

    @property
    def ghosts(self):
        try:
            return self._ghosts
        except AttributeError:
            self._ghosts = self._check_ghosts()
            return self._ghosts

    @property
    def dipoles(self):
        try:
            return self._dipoles
        except AttributeError:
            self._dipoles = ape_read_dipoles(self.output)
            return self._dipoles

    def path_in_workdir(self, filename):
        return os.path.join(self.workdir, filename)

    def to_apeinput(self):
        # Get the AE input and replace CalculationMode
        lines = self.ae_solver.to_apeinput()
        i = lines.index("CalculationMode = ae")
        lines[i] = "CalculationMode = pp + pp_test"

        # Add the variables defining the pseudization
        lines += ["# PseudoPotentials"]
        lines += self.pp_components.to_apeinput()

        lines += ["CoreCorrection = %s" % self.core_correction]
        lines += ["Llocal = %s" % self.llocal]

        # Add the variables for the PP tests.
        #PPTests = "ld + dm"
        #LogDerivativeRadius = 3.0
        #LogDerivativeEnergyMax = 1.5
        #LogDerivativeEnergyMin = -1.5
        #PPCalculationTolerance = 1e-6

        #if pptest_configurations is not None:
        #PPTestSCF = "yes"
        #PPTestOrbitals
        lines += ["PPOutputFileFormat = %s" % self.ppformat]

        #for l in lines: print(l)
        # FIXME This seems not to work!
        #lines += ["PPTestAEDir = '%s'" % self.ae_solver.aedir]
        return lines 

    def show_input(self, stream=sys.stdout):
        lines  = ["PP INPUT".center(80,"*")]
        lines += self.to_apeinput()
        lines += ["END PP INPUT".center(80,"*")]
        stream.writelines("\n".join(lines)+"\n")

    def pseudize(self):
        assert not hasattr(self, "exit_code")

        print("Starting pseudization process ...")
        start = time.time()

        # Write the input file.
        os.makedirs(self.workdir)

        # FIXME: Have to copy data dir
        #shutil.copy(os.path.join(self.ae_solver.aedir, "data"), self.workdir)
        shutil.copytree(self.ae_solver.aedir, os.path.join(self.workdir, "ae"))

        with open(self.ape_input, "w") as fh:
            fh.writelines("\n".join(self.to_apeinput()))

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.ape_exe, self.ape_input, self.ape_output, self.ape_stderr)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        self.exit_code = process.wait()

        print("Pseudization completed in %.1f s" % (time.time() - start))

        if self.exit_code != 0:
            with open(self.ape_stderr, "r") as stderr:
                raise ApeError("APE exit_code %s\n stderr: %s" % (self.exit_code, stderr.read()))

        with open(self.ape_stderr, "r") as stderr:
            error = stderr.read()
            if error:
                raise ApeError("APE exit_code %s\n stderr: %s" % (self.exit_code, error))

        # TODO Analize the output for possible warnings.

        # Check ghost-states 
        self._check_ghosts()

        # Check PP eigenvalues
        self._check_ppeigen()                                                

        # Check logarithmic derivative.

        # Dump quality control file.
        #with open(self.path_in_workdir("quality.json"), "w") as fh:
        #    json.dump(self.quality_control(), fh, indent=4, sort_keys=4)

    def _check_ppeigen(self):
        return ape_check_ppeigen(self.output)

    def _check_ghosts(self):
        """Check for presence of ghost states. Reads data from ape.out."""
        return ape_check_ghosts(self.output)

    def check_logders(self):
        "Check the quality of the log derivatives."
        merits = {}
        for (state, pp_ld) in self.pp_logders.items():
            ae_ld = self.ae_logders[state]
            rmesh = ae_ld.rmesh
            f = np.abs(np.tan(ae_ld.values) - np.tan(pp_ld.values))
            spline = UnivariateSpline(rmesh, f, s=0)
            merits[state] = spline.integral(rmesh[0], rmesh[-1])  / (rmesh[-1] - rmesh[0])

        return merits

    #def quality_control(self):
    #    return {
    #        "ghosts"        : self.ghosts,
    #        "ppeig_quality" : self._check_ppeigen,
    #        "logder_quality": self.check_logders(),
    #    }

    def plot_waves(self, show=True, savefig=None): 
        _plot_waves(self.ae_waves, pp_waves=self.pp_waves, show=show, savefig=savefig)

    def plot_logders(self, show=True, savefig=None): 
        _plot_logders(self.ae_logders, self.pp_logders, show=show, savefig=savefig)

    #def plot_potentials(self, vname):
    #def plot_densities(self):
    #def plot_core_correction(self):

##########################################################################################

class PseudoCandidate(object):

    def __init__(self, pp_generator):
        self.pp_generator = pp_generator

    # Rich comparison support.
    # This part is not trivial since we want to order pseudos
    # according to their quality factor. Pseudos with ghosts are obviously 
    # very bad but then we have to compare the error in the logder and in 
    # the eigvalues and we have to assign some priority. 
    # For the time-being, we give precedence to logder since a good logder
    # implies good eigens
    def __eq__(self, other):
        return self.quality == other.quality

    def __ne__(self, other):
        return not self == other

    #def __ge__(self, other):
    #    return self.quality >= other.quality
    #def __gt__(self, other):
    #    return self.quality > other.quality
    #def __le__(self, other):
    #    return self.quality <= other.quality
    #def __lt__(self, other):
    #    return self.quality < other.quality

    @property
    def quality(self):
        if self.ghosts:
            return -np.inf

##########################################################################################

def _plot_waves(ae_waves, pp_waves=None, show=True, savefig=None): 
    """
    Uses Matplotlib to plot the radial wavefunction (AE vs PP)

    Args:
        ae_waves:
        pp_waves:
        show:
            True to show the figure
        savefig:
            'abc.png' or 'abc.eps'* to save the figure to a file.
    """
    import matplotlib.pyplot as plt

    fig = plt.figure()

    num_waves = len(ae_waves) if pp_waves is None else len(pp_waves)

    spl_idx = 0
    for (state, ae_phi) in ae_waves.items():

        if pp_waves is not None and state not in pp_waves:
            continue
        
        spl_idx += 1
        ax = fig.add_subplot(num_waves, 1, spl_idx)

        lines, legends = [], []

        ir = len(ae_phi) + 1
        if ae_phi.isbound:
            ir = ae_phi.ir_small() + 1

        line, = ax.plot(ae_phi.rmesh[:ir], ae_phi.values[:ir], "-b", linewidth=2.0, markersize=1)

        lines.append(line)
        legends.append("AE = %s" % state)

        if pp_waves is not None:
            pp_phi = pp_waves[state]
            line, = ax.plot(pp_phi.rmesh[:ir], pp_phi.values[:ir], "^r", linewidth=2.0, markersize=4)
            lines.append(line)
            legends.append("PP = %s" % state)

        ax.legend(lines, legends, 'lower right', shadow=True)

        # Set ticks and labels.
        ax.set_ylabel("$\phi(r)$")

        #ax.set_xticks(ecut_list)
        ax.grid(True)
        if spl_idx == num_waves:
            ax.set_xlabel("r [Bohr]")

    if show:
        plt.show()
                             
    if savefig is not None:
        fig.savefig(os.path.abspath(savefig))

##########################################################################################

def _plot_logders(ae_logders, pp_logders, show=True, savefig=None): 
    """
    Uses Matplotlib to plot the logarithmic derivatives.
                                                                                         
    Args:
        show:
            True to show the figure
        savefig:
            'abc.png' or 'abc.eps'* to save the figure to a file.
    """
    import matplotlib.pyplot as plt
    assert len(ae_logders) == len(pp_logders) 
    fig = plt.figure()
                                                                                         
    num_logds, spl_idx = len(ae_logders), 0

    for (state, pp_logd) in pp_logders.items():
        spl_idx += 1
        ax = fig.add_subplot(num_logds, 1, spl_idx)
                                                                                         
        lines, legends = [], []
                                                                                         
        ae_logd = ae_logders[state]

        line, = ax.plot(ae_logd.rmesh, ae_logd.values, "-b", linewidth=2.0, markersize=1)
        lines.append(line)
        legends.append("AE logder %s" % state)
                                                                                         
        line, = ax.plot(pp_logd.rmesh, pp_logd.values, "^r", linewidth=2.0, markersize=4)
        lines.append(line)
        legends.append("PP logder %s" % state)
                                                                                         
        ax.legend(lines, legends, 'lower left', shadow=True)

        # Set ticks and labels.
        ax.set_ylabel("LOGDER")

        ax.grid(True)
        #ax.set_ylabel("$\Delta$ Etotal [meV]")
        #ax.set_xticks(ecut_list)
        if spl_idx == num_logds:
            ax.set_xlabel("r [Bohr]")

    if show:
        plt.show()
                             
    if savefig is not None:
        fig.savefig(os.path.abspath(savefig))

##########################################################################################

def ape_check_ghosts(out_lines):
    """
    Check for presence of ghost states. Reads data from ape.out.
    Example:

          Ghost state analysis:
            State: 3s
              KB energy < 0; Eref < E0       =>  No ghost states
              Local potential eigenvalues:   -0.1182 (E0)     0.0000 (E1)
              Reference energy:              -0.3981 (Eref)
            State: 3d
              KB energy < 0; Eref = E0 = 0   =>  Unable to determine
              Local potential eigenvalues:    0.0000 (E0)     0.0000 (E1)
              Reference energy:               0.0000 (Eref)
            State: 4f
              KB energy < 0; Eref > E0       =>  Ghost state found
              Local potential eigenvalues:   -3.0076 (E0)     0.0000 (E1)
              Reference energy:               0.0000 (Eref)

            Localization radii [b]:
    """
    out_lines = out_lines[:]

    SENTINEL = "Ghost state analysis:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+1:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    ghosts = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break

        if line.startswith("State:"):
            state = line.split("State:")[1].strip()
            line = out_lines.pop(0).strip()
            result = line.split("=>")[-1].lower().strip()
            result = result.split()[0]
            # See states.f90
            assert result in ["no", "ghost", "unable", "illdefined"]
            if result == "ghost":
                ghosts[state] = result

    return ghosts

##########################################################################################

def ape_check_ppeigen(out_lines):
    """
    Check the quality of the pseudo eigenvalues.

    Pseudopotentials Self-Consistency:
      State  Eigenvalue [H ]    Norm Test   Slope Test
        3s        -0.39812      1.0000003   0.9999983
        3p        -0.15331      0.9999993   0.9999941
        3d         0.00000      1.0033912   1.0017468
        4f         0.00000      0.9912183   0.9956052

    """
    class PseudoScfTest(collections.namedtuple("PseudoScfTest", "state, eig, norm, slope")):
        def __new__(cls, **kwargs):
            # Type conversion
            for k, v in kwargs.items():
                if k == "state":
                    kwargs[k] = str(v)
                else:
                    kwargs[k] = float(v)
            return super(PseudoScfTest, cls).__new__(cls, **kwargs)

    out_lines = out_lines[:]

    SENTINEL = "Pseudopotentials Self-Consistency:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+2:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    scf_tests = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break

        state, eig, norm, slope = line.split()

        scf_tests[state] = PseudoScfTest(state=state, eig=eig, norm=norm, slope=slope)


class Dipole(object):
    TOL_LDIFF = 0.002

    def __init__(self, istate, ostate, aeres, ppres):
        self.istate = istate
        self.ostate = ostate
        self.aeres = float(aeres)
        self.ppres = float(ppres)
        self.aempp = self.aeres - self.ppres

    @property
    def fulfills_lrule(self):
        if self.istate.lselect(self.ostate) and abs(self.aempp) > TOL_LDIFF:
            return False
        return True

def ape_read_dipoles(out_lines):
    #Dipole Matrix Elements:
    #      States           AE         PS
    #     3s -- 3p        2.2998     2.3074
    #     3s -- 3d        0.0003     0.0003
    #     3s -- 4f        0.0000     0.0000
    #     3p -- 3d        0.0008     0.0008
    #     3p -- 4f        0.0000     0.0000
    #     3d -- 4f       98.7760    98.7760

    out_lines = out_lines[:]
                                                                           
    SENTINEL = "Dipole Matrix Elements:"
    for (i, line) in enumerate(out_lines):
        if line.strip() == SENTINEL:
            out_lines = out_lines[i+2:]
            break
    else:
        raise ApeError("Cannot find sentinel %s in APE output" % SENTINEL)

    dipoles = {}
    while True:
        line = out_lines.pop(0).strip()
        if not line: 
            break
        tokens = line.split()

        ae_r, ps_r = map(float, [tokens[-2], tokens[-1]])
        trans = " ".join(tokens[:3])
        dipoles[trans] = {"AE": ae_r, "PS": ps_r, "AE-PS": (ae_r-ps_r)}

    return dipoles

##########################################################################################

if __name__ == "__main__":
    with_plots, show_input = False, False

    # AE configuration: neutral silicon + empty 3d-4f states
    aconf = AtomicConfiguration.from_symbol("Si")

    aconf.add_state(QState(n=3, l="d", occ=0.0))
    aconf.add_state(QState(n=4, l="f", occ=0.0))

    # Solve AE problem.
    ae_solver = AeSolver(aconf)
    ae_solver.solve()

    if show_input:
        ae_solver.show_input()

    #if with_plots:
    #    ae_solver.plot_waves()

    #for (state, wave) in ae_solver.ae_waves.items():
        #print("state %s, nroots %s inodes %s" % (state, len(wave.roots), len(wave.inodes)))
        #print("state %s, integral %s" %(state, wave.integral()))
        #if state == "1s": 
        #   for (r, v) in wave:
        #       print(v, wave.spline(r))

    # Define the parameters for the pseudization.
    pp_components = ApePPComponents.from_strings("3s:1.2:tm", "3p:1.27:tm", "3d:1.5:tm", "4f:1.9:tm")

    # For each possible local L:
    #     1) Pseudize.
    #     2) Check ghosts
    #     3) Plot wfs and logders

    for llocal in [-1,0,1,2,3,4]:
    #for llocal in [4]:
        workdir = "pptest_loc%d" % llocal
        pp_gen = PseudoGenerator(workdir, ae_solver, pp_components, core_correction=0, llocal=llocal)

        pp_gen.pseudize()

        if show_input:
            pp_gen.show_input()

        if pp_gen.ghosts:
            print("Detected ghosts for states %s" % pp_gen.ghosts.keys())

        #pprint(pp_gen.dipoles)
        print("LOGDER merit factors:")
        pprint(pp_gen.check_logders())

        # Plots
        if with_plots:
            pp_gen.plot_waves()
            pp_gen.plot_logders()

    sys.exit(0)
