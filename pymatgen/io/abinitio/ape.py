#!/usr/bin/env python
from __future__ import division, print_function

#import sys
import os
import os.path
import collections
import shutil
import numpy as np

from pprint import pprint
from scipy.interpolate import UnivariateSpline
from subprocess import Popen, PIPE

__author__ = "Matteo Giantomassi"
__copyright__ = "Copyright 2013, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Matteo Giantomassi"
__status__ = "Development"
__date__ = "$April 26, 2013M$"

##########################################################################################
# Helper functions

def ape_read_waves(dirname):
    "Read all APE radial wavefunctions located in directory dirname"
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
    "Read all APE radial potentials located in directory dirname."
    pots = {}
    for filename in os.listdir(dirname):
        if not filename.startswith("v_"):
            continue
        path = os.path.join(dirname, filename)
                                                                          
        tokens = filename.split("_")
        assert len(tokens) == 2
        pot_name = tokens[1]

        #TODO check spin and spinor
        # Load [r, v(r)] in data 
        data = np.loadtxt(path)
                                                                       
        pots[pot_name] = RadialFunction(pot_name, data[:,0], data[:,1])

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
    "Read log derivatives located in directory dirname."
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

#Table
# 'element': (atomic number, (n, l, occ, energy))
#see gpaw/atom/configuration
# Eg 'Be': (4, [(1, 0, 2, -3.856411), (2, 0, 2, -0.20574400000000001)]),
#class Configuration(collections.namedtuple("Configuration", "Z, n, l, occ")):
#    pass

#  type qn_t
#    integer  :: n    
#    integer  :: l    
#    real(R8) :: j
#    real(R8) :: s    
#    real(R8) :: m
#    real(R8) :: sg
#    integer  :: k
#  end type qn_t

# Spinor 
#class State(object):
#    def __init__(self, n, l, j, occ):
#        self.n = n
#        self.l = l
#        self.j = j
#        self.occ = occ

#  type state_t
#    ! General information about the state
#    type(qn_t)       :: qn    ! state quantum numbers
#    real(R8)         :: occ   ! occupation
#    real(R8)         :: ev    ! eigenvalue
#    character(len=6) :: label ! a label to identify the state
#    integer          :: wave_eq ! wave-equation used to obtain the wave-functions
#
#    ! The wavefunctions
#    integer :: np, wf_dim
#    real(R8), pointer :: wf(:,:) !  Schrodinger equation:
#                                 !   wf(:,1) -> wavefunction
#                                 !  Scalar-relativistic equation
#                                 !   wf(:,1) -> wavefunction
#                                 !  Dirac equation:
#                                 !   wf(:,1) -> spinor major component
#                                 !   wf(:,2) -> spinor minor component
#                                 !  Dirac equation with spin-polarization:
#                                 !   wf(:,1) -> spinor major component +
#                                 !   wf(:,2) -> spinor minor component +
#                                 !   wf(:,3) -> spinor major component -
#                                 !   wf(:,4) -> spinor minor component -
#    real(R8), pointer :: wfp(:,:) ! Derivative of the wavefunction
#
#    ! Some information about the wavefunctions
#    real(R8) :: peak ! outermost peak position
#    real(R8) :: node ! outermost node position
#  end type state_t
#
#   

#class Configuration(object):
#
#    def __init__(self):
#        pass
#
#    @classmethod
#    def neutral(self):
#        pass
#
#    def to_apeinput(self):
#        lines = ["%Orbitals",]
#        lines += ["%"]
#        return lines

class RadialFunction(object):
    "A RadialFunction has a radial mesh rmesh and values defined on this mesh."
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

    def __repr__(self):
        return "<%s at %s, name = %s>" % (self.__class__.__name__, id(self), self.name)

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
        Returns the rightmost index where the abs wavefunction value becomes greater than abs_tol 
        .. warning:
            Assumes that self.values are tending to zero for r --> infinity.
        """
        for i in range(len(self.rmesh)-1,-1,-1):
            if abs(self.values[i]) > abs_tol:
                break
        return i

##########################################################################################

class RadialWaveFunction(RadialFunction):

    def __init__(self, state, name, rmesh, values):
        super(RadialWaveFunction, self).__init__(name, rmesh, values)
        self.state = state

    @property
    def isbound(self):
        "True if self is a bound state."
        back = min(10, len(self))
        return all(abs(self.values[-back:]) < 1.0e-12)

#class DiracWaveFunction(RadialFunction):

##########################################################################################

class AeAtom(object):
    """
    Class for representing the solution of the all-electron atom problem

    Parameters:
    
    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).

    orbitals
    """
    def __init__(self, symbol, orbitals=None):
        self.symbol = symbol
        #self.nuclear_charge = z_from_symbol(symbol) 

        #if orbitals is None:
        #    raise NotImplementedError("")
        #else:
        #    self.orbitals = orbitals

    #def __repr__(self):

    @property
    def title(self):
        return self.symbol

    @classmethod
    def neutral(cls, symbol):
        return cls(symbol, orbitals=neutral_configuration)

    #def add(self, n, l, df=+1, s=None):
    #    "Add (remove) electrons."

    def to_apeinput(self):
        lines = []
        lines += ["Title = %s" % self.title]
        #lines += ["NuclearCharge = %s " % self.nuclear_charge]
        #lines += self.orbitals.to_apeinput()
        return lines

##########################################################################################

class ApeRadialMesh(dict):
    """
    APE uses a mesh to store functions. In order to change the mesh parameters you can use the following options:

        MeshType (integer, log1). Possible values are:
            - lin: Linear
            - log1: Logarithmic [ri = b*exp(a*i)]
            - log2: Logarithmic [ri = b*(exp(a*i) - 1)]

        MeshStartingPoint (real, sqrt(NuclearCharge)*1e-5 a.u.):
            Sets the starting point of the mesh.

        MeshOutmostPoint (real, sqrt(NuclearCharge)*30 a.u.):
            Sets the outmost point of the mesh.

        MeshNumberOfPoints (integer, sqrt(NuclearCharge)*200):
            Sets the mesh number of points.

        MeshDerivMethod (int, cubic_spline)
            Numerical method used to compute derivatives. Possible choices are:

            - cubic_spline: Use cubic splines to interpolate the function and then uses this interpolation to compute the derivatives.
            - finite_diff: Finite differences. The number of points used is controled by the MeshFiniteDiffOrder variable.

        MeshFiniteDiffOrder (int, 4)
        This variable controls the numbers of points used for the finite differences discretization. 
        To compute the derivative at a given point, the finite difference formula will use MeshFiniteDiffOrder 
        points to the left and MeshFiniteDiffOrder points to the right. This means that in total MeshFiniteDiffOrder*2 + 1 points will be used.
    """

    _KEYS = [
        "MeshType",
        "MeshStartingPoint",
        "MeshOutmostPoint",
        "MeshNumberOfPoints", 
        "MeshDerivMethod",
        "MeshFiniteDiffOrder",
    ]

    def __init__(self, **kwargs):
        super(RadialMesh, self).__init__(**kwargs)

        for k in self:
            if k not in self._KEYS:
                raise ValueError("%s is not a registered key" % k)

    def to_apeinput(self):
        lines = []
        for k,v in self.items():
            lines += ["%s = %s" % k, v]
        return lines

##########################################################################################

class ApeError(Exception):
    "Base class for APE Exceptions"

class AeSolver(object):
    """
    An AeSolver has a all-electron atom and solves the AE equation with the 
    algorith and the paramenters passed to the constructor.
    """
    ape_exe = "ape"

    def __init__(self, ae_atom, workdir=None, xcfunctional='lda_x+lda_c_pw', spin_mode="unpolarized", 
                 wave_equation="schrodinger", ape_radmesh=None, ape_control=None, verbose=0):
        """
        XCFunctional (integer, lda_x+lda_c_pw)

        WaveEquation (integer, schrodinger)

        spin_mode:
            - unpolarized: Spin-unpolarized calculation.
            - polarized: Spin-polarized calculation.

        When performing atomic calculations APE can solve either the Kohn-Sham equations, 
        the Dirac-Kohn-Sham equations or the scalar-relativistic Khon-Sham equations. Valid options are:

        - schrodinger: Kohn-Sham equations.
        - dirac: Dirac-Kohn-Sham equations.
        - scalar_rel: scalar-relativistic Kohn-Sham equations.

        6.2.2 SpinMode (integer, unpolarized)
            - unpolarized: Spin-unpolarized calculation.
            - polarized: Spin-polarized calculation.
        """
        self.ae_atom = ae_atom
        self.xcfunctional = xcfunctional
        self.spin_mode = spin_mode
        self.wave_equations = wave_equation

        self.ape_radmesh = ape_radmesh
        self.ape_control = ape_control
        self.verbose = verbose

        self.workdir = os.path.abspath("helloape")
        #if workdir is not None:
        #    self.workdir = os.path.abspath(workdir)

    @property
    def ape_output(self):
        return os.path.join(self.workdir, "ape.out")

    @property
    def ape_stderr(self):
        return os.path.join(self.workdir, "ape.stderr")

    @property
    def aedir(self):
        "Absolute path to the ae directory containing the all-electron results." 
        return os.path.join(self.workdir, "ae")

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
        #lines = []
        #lines += self.ae_atom.to_apeinput() 
        ##xc_type='LDA', spin_mode="unpolarized", wave_equation="schrodinger", 
        #lines += self.ape_radmesh.to_apeinput()
        #lines += self.ape_control.to_apeinput()
        #return lines
        # HACK
        with open("ae_si", "r") as fh:
            return fh.readlines()

    def solve(self):
        assert not hasattr(self, "exit_code")
        # Write the input file.
        os.makedirs(self.workdir)
        self.ape_input = os.path.join(self.workdir, "ape.in")

        with open(self.ape_input, "w") as fh:
            fh.writelines(self.to_apeinput())

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.ape_exe, self.ape_input, self.ape_output, self.ape_stderr)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        self.exit_code = process.wait()

        if self.exit_code != 0:
            with open(self.ape_stderr, "r") as stderr:
                raise ApeError("APE returned exit_code %s\n stderr:" % (self.exit_code, stderr.read()))

        # TODO Analize the output for possible warnings.
        #self.validate()

    def plot_waves(self, show=True, savefig=None): 
        ape_plot_waves(self.ae_waves, show=show, savefig=savefig)

##########################################################################################

class ApeControl(dict):
    """
    6.7 SCF

    The self consistent field procedure will stop when one of the convergence criteria is fulfilled. 
    At each iteration the new guess potential is built mixing the input and output potentials.

    MaximumIter (integer, 300):      
        Maximum number of SCF iterations. When that number is reached the SCF procedure will stop, 
        even if none of the criteria are fulfilled, and the calculation will continue as normal. 0 means unlimited.

    ConvAbsDens (real, 0.0):
        Absolute convergence of the density. 0 means do not use this criterion.

    ConvRelDens (real, 1e-8):
        Relative convergence of the density. 0 means do not use this criterion.

    ConvAbsEnergy (real, 0.0):
        Absolute convergence of the total energy. 0 means do not use this criterion.

    ConvRelEnergy (real, 0.0)
        Relative convergence of the total energy. 0 means do not use this criterion.

    SmearingFunction (integer, fixed_occ):
        Select how the states should be occupied. The option are:

            - fixed_occ: the occupancies of the states are fixed.
            - semiconducting: semiconducting occupancies, i.e., the lowest lying states are occupied until no more electrons are left
            - averill_painter: F.W. Averill and G.S. Painter, Phys. Rev. B 46, 2498 (1992)

    MixingScheme (integer, broyden):
        Selects the mixing procedure to be used during the SCF cycle. Possible values are:
            - linear: Linear mixing.
            - broyden: Broyden mixing.

    Mixing (real, 0.3):
        Determines the amount of the new potential which is to be mixed with the old one.

    MixNumberSteps (integer, 3):
        Number of steps used by Broyden mixing to extrapolate the new potential.

    Wave-equations solver
    APE solves the Schrodinger and Dirac equations using the ODE solver from GSL.

    6.8.1 EigenSolverTolerance (real, 1.0e-8 a.u.)
        The eigensolver will improve the eigenvalues untill the error estimate of each eigenvalue is smaller than this tolerance.

    6.8.2 ODEIntTolerance (real, 1.0e-12)

    6.8.3 ODESteppingFunction (integer, rkpd8). Possible values are:
        - rk2: Embedded 2nd order Runge-Kutta method
        - rk4: 4th order (classical) Runge-Kutta method
        - rkf4: Embedded 4th order Runge-Kutta-Fehlberg method
        - rkck4: Embedded 4th order Runge-Kutta Cash-Karp method
        - rkpd8: Embedded 8th order Runge-Kutta Prince-Dormand method
    For more information about the stepping functions please refer to the GSL documentation.

    6.8.4 ODEMaxSteps (integer, 500000):
        Sometimes it may happen that the ODE solver takes too many steps, because of some problem. 
        To avoid that, APE will issue a fatal error if the number of steps is greater than ODEMaxSteps.
    """

    _KEYS = [
        # SCF
        "MaximumIter", 
        "ConvAbsDens",
        "ConvRelDens",
        "ConvAbsEnergy",
        "ConvRelEnergy",
        "SmearingFunction",
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
        super(SCFControl, self).__init__(**kwargs)
                                                                   
        for k in self:
            if k not in self._KEYS:
                raise ValueError("%s is not a registered key" % k)
                                                                   
    def to_apeinput(self):
        lines = []
        for k,v in self.items():
            lines += ["%s = %s" % k, v]
        return lines

##########################################################################################

class PseudoGenerator(object):
    "A PseudoGenerator uses the data produces by the AeSolver the construct the pseudopotential"
    ape_exe = "ape"

    def __init__(self, workdir, ae_solver, format="abinit6"):
        self.workdir = os.path.abspath(workdir)
        self.ae_solver = ae_solver
        self.format = format

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

    def to_apeinput(self):
        # HACK
        with open("pp_si", "r") as fh:
            lines = fh.readlines()

        # FIXME does not work!
        #lines += ["PPTestAEDir = '%s'" % self.ae_solver.aedir]
        return lines 

    def pseudize(self):
        assert not hasattr(self, "exit_code")
        # Write the input file.
        os.makedirs(self.workdir)

        # FIXME: Have to copy data dir
        #shutil.copy(os.path.join(self.ae_solver.aedir, "data"), self.workdir)
        shutil.copytree(self.ae_solver.aedir, os.path.join(self.workdir, "ae"))

        self.ape_input = os.path.join(self.workdir, "ape.in")

        with open(self.ape_input, "w") as fh:
            fh.writelines(self.to_apeinput())

        # Launch APE in a subprocess.
        command = "%s < %s > %s 2> %s" % (self.ape_exe, self.ape_input, self.ape_output, self.ape_stderr)

        process = Popen(command, shell=True, cwd=self.workdir, stderr=PIPE)

        self.exit_code = process.wait()

        if self.exit_code != 0:
            with open(self.ape_stderr, "r") as stderr:
                raise ApeError("APE exit_code %s\n stderr:" % (self.exit_code, stderr.read()))

        with open(self.ape_stderr, "r") as stderr:
            error = stderr.read()
            if error:
                raise ApeError("APE exit_code %s\n stderr:" % (self.exit_code, error))

        # TODO Analize the output for possible warnings.

        # Check ghost-states and PP eigenvalues
                                                
        # Check logarithmic derivative.

    #def check_pseudo_quality(self):
    #    eig_err = self.check_ppeigen()
    #    has_ghosts = self.check_ghosts()
    #    lds_err = self.check_logders()

    def check_ppeigen(self):
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
                new = super(PseudoScfTest, cls).__new__(cls, **kwargs)
                return new

        out_lines = self.output

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
            #print(scf_tests[state])

    def check_ghosts(self):
        """
        Check for presence of ghost states. Reads data from ape.out.
        Example:

              Ghost state analysis:
                State: 3s
                  KB energy < 0; Eref < E0       =>  No ghost states
                  Local potential eigenvalues:   -0.1182 (E0)     0.0000 (E1)
                  Reference energy:              -0.3981 (Eref)
                State: 3p
                  KB energy < 0; Eref < E0       =>  No ghost states
                  Local potential eigenvalues:   -0.0190 (E0)     0.0000 (E1)
                  Reference energy:              -0.1533 (Eref)
                State: 3d
                  KB energy < 0; Eref = E0 = 0   =>  Unable to determine
                  Local potential eigenvalues:    0.0000 (E0)     0.0000 (E1)
                  Reference energy:               0.0000 (Eref)
                State: 4f
                  KB energy < 0; Eref = E0 = 0   =>  Unable to determine
                  Local potential eigenvalues:    0.0000 (E0)     0.0000 (E1)
                  Reference energy:               0.0000 (Eref)

                Localization radii [b]:
        """
        out_lines = self.output

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
                ans = line.split("=>")[-1].lower().strip()
                ans = ans.split()[0]
                # See states.f90
                assert ans in ["no", "ghost", "unable", "illdefined"]
                ghosts[state] = ans

        print(ghosts)

    def check_logders(self):
        "Check the quality of the log derivatives."
        for state, pp_ld in self.pp_logders.items():
            ae_ld = self.ae_logders[state]
            diff = abs(pp_ld.values - ae_ld.values)
            pprint(diff)

    def plot_waves(self, show=True, savefig=None): 
        ape_plot_waves(self.ae_waves, pp_waves=self.pp_waves, show=show, savefig=savefig)

    def plot_logders(self, show=True, savefig=None): 
        """
        Uses Matplotlib to plot the radial wavefunction (AE vs PP)
                                                                                             
        Args:
            show:
                True to show the figure
            savefig:
                'abc.png' or 'abc.eps'* to save the figure to a file.
        """
        import matplotlib.pyplot as plt
                                                                                             
        fig = plt.figure()
                                                                                             
        num_logds, spl_idx = len(self.ae_logders), 0

        for (state, pp_logd) in self.pp_logders.items():
            spl_idx += 1
            ax = fig.add_subplot(num_logds,1,spl_idx)
                                                                                             
            lines, legends = [], []
                                                                                             
            ae_logd = self.ae_logders[state]
                                                                                             
            line, = ax.plot(pp_logd.rmesh, pp_logd.values, "-->", linewidth=1.0, markersize=1)
            lines.append(line)
            legends.append("PP logder %s" % state)
                                                                                             
            line, = ax.plot(ae_logd.rmesh, ae_logd.values, "-->", linewidth=1.0, markersize=1)
            lines.append(line)
            legends.append("AE logder %s" % state)
                                                                                             
            # TODO
            #i, rsmall = ae_phi.ir_small()
            #ax.xaxis.set_view_interval(-1, 10)
            #ax.yaxis.set_view_interval(-10, emax + 0.01 * abs(emax))
                                                                                
            ax.legend(lines, legends, 'upper right', shadow=True)
                                                                                             
            # Set xticks and labels.
            #ax.grid(True)
            #ax.set_xlabel("r")
            #ax.set_ylabel("$\Delta$ Etotal [meV]")
            #ax.set_xticks(ecut_list)
                                                                                             
        if show:
            plt.show()
                                 
        if savefig is not None:
            fig.savefig(savefig)

    #def plot_potentials(self, vname):
    #def plot_density(self):
    #def plot_corecorrection(self):

##########################################################################################

def ape_plot_waves(ae_waves, pp_waves=None, show=True, savefig=None): 
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
        ax = fig.add_subplot(num_waves,1,spl_idx)

        lines, legends = [], []

        ir = len(ae_phi) + 1
        if ae_phi.isbound:
            ir = ae_phi.ir_small() + 1

        line, = ax.plot(ae_phi.rmesh[:ir], ae_phi.values[:ir], "-b", linewidth=1.0, markersize=1)

        lines.append(line)
        legends.append("AE = %s" % state)

        if pp_waves is not None:
            pp_phi = pp_waves[state]
            line, = ax.plot(pp_phi.rmesh[:ir], pp_phi.values[:ir], "*r", linewidth=2.0, markersize=2)
            lines.append(line)
            legends.append("PP = %s" % state)

        ax.legend(lines, legends, 'upper right', shadow=True)

        # Set xticks and labels.
        ax.grid(True)
        ax.set_ylabel("$\phi(r)$")
        #ax.set_xticks(ecut_list)
        if spl_idx == num_waves:
            ax.set_xlabel("r [Bohr]")

    if show:
        plt.show()
                             
    if savefig is not None:
        fig.savefig(savefig)

##########################################################################################

if __name__ == "__main__":

    ae_atom = AeAtom("Si")
    #ae_atom = AeAtom("Si+")

    aesolver = AeSolver(ae_atom)
    #print(ae.to_apeinput())
    aesolver.solve()

    pprint(aesolver.ae_waves)
    pprint(aesolver.ae_potentials)
    pprint(aesolver.ae_densities)

    aesolver.plot_waves()

    for (state, wave) in aesolver.ae_waves.items():
        #print(wave)
        print("state %s, nroots %s inodes %s" % (state, len(wave.roots),len(wave.inodes)))
        #print("state %s, integral %s" %(state, wave.integral()))
        #if state == "1s": 
            #pprint(wave.values)
            #for r, v in wave:
            #    print(v, wave.spline(r))

    ppgen = PseudoGenerator("pptest", aesolver)
    ppgen.pseudize()

    ppgen.check_ghosts()
    ppgen.check_ppeigen()
    ppgen.plot_waves()
    #ppgen.plot_logders()
