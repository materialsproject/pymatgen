# from structure_analyzer

from __future__ import division, print_function, unicode_literals

import math
from math import pi
import numpy as np
import itertools
import collections

from warnings import warn
#from pyhull.voronoi import VoronoiTess
from pymatgen import PeriodicSite
from pymatgen import Element, Specie, Composition
from pymatgen.util.num_utils import abs_cap
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder



# helper functions

def inbounds(lower, value, upper):

    """
    Returns the input value, but constraining it to lie in the interval
    between lower and upper bound.  If the original value was outside
    the interval [lower,upper] a warning is issued.

    Args:
        lower (float):
            lower bound to which to constrain value
        value (float):
            value to be kept in bounds between lower and upper
        upper (float):
            upper bound to which to constrain value
    """

    if value < lower:
        return lower
    elif value > upper:
        return upper
    else:
        return value


def gramschmidt(vin, uin):

    """
    Returns that part of the first input vector
    that is orthogonal to the second input vector.
    The output vector is not normalized.

    Args:
        vin (numpy array):
            first input vector
        uin (numpy array):
            second input vector
    """

    vin_uin = np.inner(vin, uin)
    uin_uin = np.inner(uin, uin)
    if uin_uin <= 0.0:
        raise ValueError("Zero or negative inner product!")
    vout = np.array([0, 0, 0], float)
    vout = vin - (vin_uin / uin_uin) * uin

    return vout


def normalize(vec):

    """
    Returns the input vector in a normalized form.

    Args:
        vec (numpy array):
            vector to be normalized
    """

    vout = np.copy(vec)
    norm2 = np.inner(vout, vout)
    if norm2 == 0.0:
        pass
    elif norm2 < 0.0:
        raise ValueError("Encountered a negative norm!")
    else:
        norm = math.sqrt(norm2)
        vout = vout / norm

    return vout


def vlength(v):

    """
    Returns the length of the input vector.

    Args:
        v (numpy array):
            vector for which length is to be computed
    """

    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

class OrderParameters(object):

    """
    This class permits the calculation of various types of local order
    parameters.
    """


    __supported_types = ("cn", "tet", "oct", "bcc", "q2", "q4", "q6")


    def __init__(self, types, parameters=None, cutoff=-10.0):
        """
        Create an OrderParameter analyzer instance.

        Args:
            types ([string]):
                List of strings representing the types of order parameters
                to be calculated. Note that multiple mentions of the
                same type may occur. Currently available types are
                "cn"  (simple coordination number---normalized,
                      if desired),
                "tet" [Shetty--Peters style OP recognizing tetrahedral
                      coordination (Zimmermann et al.,
                      J. Am. Chem. Soc., 137, 13352-13361, 2015)],
                "oct" [Shetty--Peters style OP recognizing octahedral
                      coordination (Zimmermann et al.,
                      J. Am. Chem. Soc., 137, 13352-13361, 2015)],
                "bcc" [Shetty--Peters style OP recognizing local
                      body-centered cubic environment (Peters,
                      J. Chem. Phys., 131, 244103, 2009)],
                "q2"  [Bond orientational order parameter (BOOP)
                      of weight l=2 (Steinhardt et al., Phys. Rev. B,
                      28, 784-805, 1983)],
                "q4"  (BOOP of weight l=4),
                "q6"  (BOOP of weight l=6).
            parameters ([[float]]):
                2D list of floating point numbers that store
                parameters associated with the different order parameters
                that are to be calculated (1st dimension = length of
                types tuple; any 2nd dimension may be zero, in which case
                default values are used). In the following, those order
                parameters q_i are listed that require further parameters
                for their computation (values in brackets denote default
                values):
                  "cn":  normalizing constant (1);
                  "tet": Gaussian width in fractions of pi (180 degrees)
                         reflecting the "speed of
                         penalizing" deviations away from the perfect
                         tetrahedral angle of any individual
                         neighbor1-center-neighbor2 configuration (0.0667);
                  "oct": threshold angle in degrees distinguishing a second
                         neighbor to be either close to the south pole or
                         close to the equator (160.0);
                         Gaussian width for penalizing deviations away
                         from south pole (0.0667);
                         Gaussian width for penalizing deviations away
                         from equator (0.0556);
                         constant for shifting q_oct toward smaller
                         values, which can be helpful when trying to fine-
                         tune the capabilities of distinguishing between
                         different environments (e.g., tet vs oct)
                         given a single mutual threshold q_thresh;
                  "bcc": south-pole threshold angle as for "oct" (160.0);
                         south-pole Gaussian width as for "oct" (0.0667).
            cutoff (float):
                Cutoff radius to determine which nearest neighbors are
                supposed to contribute to the order parameters.
                If the value is negative the neighboring sites found by
                distance and cutoff radius are further
                pruned using the get_coordinated_sites method from the
                VoronoiCoordFinder class.
        """
        if len(types) == 0:
            raise ValueError("Empty types list!")
        for t in types:
            if t not in OrderParameters.__supported_types:
                raise ValueError("Unknown order parameter type (" + \
                                 t + ")!")
        if parameters:
            if len(types) != len(parameters):
                raise ValueError("1st dimension of parameters array is not"
                                 " consistent with types list!")
            for lp in parameters:
                if len(lp) > 0:
                    for p in lp:
                        if type(p) != float and type(p) != int:
                            raise AttributeError("Expected only float and"
                                                 " integer type parameters!")
            loc_parameters = list(parameters)
        else:
            loc_parameters = [[] for t in types]
        self._types = tuple(types)
        tmpparas = []
        self._computerijs = self._geomops = self._boops = False
        self._max_trig_order = -1
        for i, t in enumerate(self._types):
            # add here any additional parameter checking and
            #     default value assignment
            tmpparas.append([])
            if t == "cn":
                if len(loc_parameters[i]) == 0:
                    tmpparas[i].append(1.0)
                else:
                    if loc_parameters[i][0] == 0.0:
                        raise ValueError("Normalizing constant for"
                                         " coordination-number based order"
                                         " parameter is zero!")
                    else:
                        tmpparas[i].append(loc_parameters[i][0])
            elif t == "tet":
                if len(loc_parameters[i]) == 0:
                    tmpparas[i].append(1.0/0.0667)
                else:
                    if loc_parameters[i][0] == 0.0:
                        raise ValueError("Gaussian width for"
                                         " tetrahedral order"
                                         " parameter is zero!")
                    else:
                        tmpparas[i].append(1.0/loc_parameters[i][0])
            elif t == "oct":
                if len(loc_parameters[i]) < 4:
                    tmpparas[i].append(8.0*pi/9.0)
                    tmpparas[i].append(1.0/0.0667)
                    tmpparas[i].append(1.0/0.0556)
                    tmpparas[i].append(0.25)
                    tmpparas[i].append(4.0/3.0)
                else:
                    if loc_parameters[i][0]<=0.0 or loc_parameters[i][0]>=180.0:
                        warn("Threshold value for south pole"
                             " configurations in octahedral order"
                             " parameter outside ]0,180[")
                    tmpparas[i].append(loc_parameters[i][0]*pi/180.0)
                    if loc_parameters[i][1] == 0.0:
                        raise ValueError("Gaussian width for south pole"
                                         " configurations in octahedral"
                                         " order parameter is zero!")
                    else:
                        tmpparas[i].append(1.0/loc_parameters[i][1])
                    if loc_parameters[i][2] == 0.0:
                        raise ValueError("Gaussian width for equatorial"
                                         " configurations in octahedral"
                                         " order parameter is zero!")
                    else:
                        tmpparas[i].append(1.0/loc_parameters[i][2])
                    if loc_parameters[i][3]-1.0 == 0.0:
                        raise ValueError("Shift constant may not be"
                                         " unity!")
                    if loc_parameters[i][3] < 0.0 or loc_parameters[i][3] > 1.0:
                        warn("Shift constant outside [0,1[.")
                    tmpparas[i].append(loc_parameters[i][3])
                    tmpparas[i].append(1.0/(1.0-loc_parameters[i][3]))
            elif t == "bcc":
                if len(loc_parameters[i]) < 2:
                    tmpparas[i].append(8.0*pi/9.0)
                    tmpparas[i].append(1.0/0.0667)
                else:
                    if loc_parameters[i][0]<=0.0 or loc_parameters[i][0]>=180.0:
                        warn("Threshold value for south pole"
                             " configurations in bcc order"
                             " parameter outside ]0,180[")
                    tmpparas[i].append(loc_parameters[i][0]*pi/180.0)
                    if loc_parameters[i][1] == 0.0:
                        raise ValueError("Gaussian width for south pole"
                                         " configurations in bcc"
                                         " order parameter is zero!")
                    else:
                        tmpparas[i].append(1.0/loc_parameters[i][1])
            # All following types should be well-defined/-implemented,
            # and they should not require parameters.
            elif t != "q2" and t != "q4" and t != "q6":
                raise ValueError("unknown order-parameter type \""+t+"\"")

            # Add here any additional flags to be used during calculation.
            if t == "tet" or t == "oct" or t == "bcc":
                self._computerijs = self._geomops = True
            if t == "q2" or t == "q4" or t == "q6":
                self._computerijs = self._boops = True
            if t == "q2" and self._max_trig_order < 2:
                self._max_trig_order = 2
            if t == "q4" and self._max_trig_order < 4:
                self._max_trig_order = 4
            if t == "q6" and self._max_trig_order < 6:
                self._max_trig_order = 6

        # Finish parameter treatment.
        self._paras = list(tmpparas)
        if cutoff < 0.0:
            self._cutoff = -cutoff
            self._voroneigh = True
        elif cutoff > 0.0:
            self._cutoff = cutoff
            self._voroneigh = False
        else:
            raise ValueError("Cutoff radius is zero!")

        # Further variable definitions.
        self._last_nneigh = -1
        self._pow_sin_t = {}
        self._pow_cos_t = {}
        self._sin_n_p = {}
        self._cos_n_p = {}


    @property
    def num_ops(self):

        """"
        Returns the number of different order parameters that are targeted
        to be calculated.
        """

        return len(self._types)


    @property
    def last_nneigh(self):

        """"
        Returns the number of neighbors encountered during the most
        recent order-parameter calculation. A value of -1 indicates that
        no such calculation has yet been performed for this instance.
        """

        return len(self._last_nneigh)


    def compute_trigonometric_terms(self, thetas=[], phis=[]):

        """"
        Computes trigonometric terms that are required to
        calculate bond orientational order parameters.

        Args:
            thetas ([float]):
                polar angles of all neighbors in radians.
            phis ([float]):
                azimuth angles of all neighbors in radians.  The list of
                azimuth angles is expected to have the same size as the list
                of polar angles; otherwise, a ValueError is raised.  Also,
                the two lists of angles have to be coherent in order. That
                is, it is expected that the order in the list of azimuth
                angles corresponds to a distinct sequence of neighbors.
                And, this sequence has to equal the sequence
                of neighbors in the list of polar angles.

        """

        if not thetas or not phis:
            raise ValueError("No entries in list/s of polar and/or azimuthal" \
                             " angles!")
        if len(thetas) != len(phis):
            raise ValueError("List of polar and azimuthal angles have to be"
                             " equal!")

        self._pow_sin_t.clear()
        self._pow_cos_t.clear()
        self._sin_n_p.clear()
        self._cos_n_p.clear()

        self._pow_sin_t[1] = [math.sin(float(t)) for t in thetas]
        self._pow_cos_t[1] = [math.cos(float(t)) for t in thetas]
        self._sin_n_p[1]   = [math.sin(float(p)) for p in phis]
        self._cos_n_p[1]   = [math.cos(float(p)) for p in phis]

        for i in range(2, self._max_trig_order+1):
            self._pow_sin_t[i] = [e[0]*e[1] for e in zip(
                self._pow_sin_t[i-1], self._pow_sin_t[1])]
            self._pow_cos_t[i] = [e[0]*e[1] for e in zip(
                self._pow_cos_t[i-1], self._pow_cos_t[1])]
            self._sin_n_p[i]  = [math.sin(float(i)*float(p)) \
                for p in phis]
            self._cos_n_p[i]  = [math.cos(float(i)*float(p)) \
                for p in phis]


    def get_q2(self, thetas=[], phis=[]):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=2.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]):
                polar angles of all neighbors in radians.
            phis ([float]):
                azimuth angles of all neighbors in radians.

        Return:
            q2 (float): bond orientational order parameter of weight l=2
                corresponding to the input angles thetas and phis.
        """

        if thetas and phis:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        sqrt_15_2pi = math.sqrt(15.0/(2.0*pi))
        sqrt_5_pi   = math.sqrt(5.0/pi)

        pre_y_2_2 = [0.25 * sqrt_15_2pi * val for val in self._pow_sin_t[2]]
        pre_y_2_1 = [0.5 * sqrt_15_2pi * val[0] * val[1]
            for val in zip(self._pow_sin_t[1],self._pow_cos_t[1])]

        acc = 0.0

        # Y_2_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        # Y_2_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_2_0
        real = imag = 0.0
        for i in nnn_range:
            real += 0.25 * sqrt_5_pi * (3.0*self._pow_cos_t[2][i] - 1.0)
        acc += (real*real)

        # Y_2_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_2_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_2_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_2_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_2_2[i] * self._cos_n_p[2][i]
            imag += pre_y_2_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        q2 = math.sqrt(4.0 * pi * acc / (5.0 * float(nnn*nnn)))
        return q2


    def get_q4(self, thetas=[], phis=[]):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=4.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]):
                polar angles of all neighbors in radians.
            phis ([float]):
                azimuth angles of all neighbors in radians.

        Return:
            q4 (float): bond orientational order parameter of weight l=4
                corresponding to the input angles thetas and phis.
        """

        if thetas and phis:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        i16_3 = 3.0 / 16.0
        i8_3 = 3.0 / 8.0

        sqrt_35_pi  = math.sqrt(35.0/pi)
        sqrt_35_2pi = math.sqrt(35.0/(2.0*pi))
        sqrt_5_pi   = math.sqrt(5.0/pi)
        sqrt_5_2pi  = math.sqrt(5.0/(2.0*pi))
        sqrt_1_pi   = math.sqrt(1.0/pi)

        pre_y_4_4 = [i16_3 * sqrt_35_2pi * val for val in self._pow_sin_t[4]]
        pre_y_4_3 = [i8_3 * sqrt_35_pi * val[0] * val[1] \
            for val in zip(self._pow_sin_t[3],self._pow_cos_t[1])]
        pre_y_4_2 = [i8_3 * sqrt_5_2pi * val[0] * (7.0*val[1] - 1.0) \
            for val in zip(self._pow_sin_t[2],self._pow_cos_t[2])]
        pre_y_4_1 = [i8_3 * sqrt_5_pi * val[0] * (7.0*val[1] - 3.0*val[2]) \
            for val in zip(self._pow_sin_t[1],self._pow_cos_t[3], \
                           self._pow_cos_t[1])]

        acc = 0.0

        # Y_4_-4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += (real*real + imag*imag)

        # Y_4_-3
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += (real*real + imag*imag)

        # Y_4_-2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        # Y_4_-1
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_4_0
        real = imag = 0.0
        for i in nnn_range:
            real += i16_3 * sqrt_1_pi * (35.0*self._pow_cos_t[4][i] - \
                30.0*self._pow_cos_t[2][i] + 3.0)
        acc += (real*real)

        # Y_4_1
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_4_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_4_2
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_2[i] * self._cos_n_p[2][i]
            imag += pre_y_4_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        # Y_4_3
        real = imag = 0.0
        for i in nnn_range:
            real -= pre_y_4_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_4_3[i] * self._sin_n_p[3][i]
        acc += (real*real + imag*imag)

        # Y_4_4
        real = imag = 0.0
        for i in nnn_range:
            real += pre_y_4_4[i] * self._cos_n_p[4][i]
            imag += pre_y_4_4[i] * self._sin_n_p[4][i]
        acc += (real*real + imag*imag)

        q4 = math.sqrt(4.0 * pi * acc / (9.0 * float(nnn*nnn)))
        return q4

    def get_q6(self, thetas=[], phis=[]):

        """
        Calculates the value of the bond orientational order parameter of
        weight l=6.  If the function is called with non-empty lists of
        polar and azimuthal angles the corresponding trigonometric terms
        are computed afresh.  Otherwise, it is expected that the
        compute_trigonometric_terms function has been just called.

        Args:
            thetas ([float]):
                polar angles of all neighbors in radians.
            phis ([float]):
                azimuth angles of all neighbors in radians.

        Return:
            q6 (float): bond orientational order parameter of weight l=6
                corresponding to the input angles thetas and phis.
        """

        if thetas and phis:
            self.compute_trigonometric_terms(thetas, phis)
        nnn = len(self._pow_sin_t[1])
        nnn_range = range(nnn)

        i64 = 1.0 / 64.0
        i32 = 1.0 / 32.0
        i32_3 = 3.0 / 32.0
        i16 = 1.0 / 16.0

        sqrt_3003_pi = math.sqrt(3003.0/pi)
        sqrt_1001_pi = math.sqrt(1001.0/pi)
        sqrt_91_2pi  = math.sqrt(91.0/(2.0*pi))
        sqrt_1365_pi = math.sqrt(1365.0/pi)
        sqrt_273_2pi = math.sqrt(273.0/(2.0*pi))
        sqrt_13_pi   = math.sqrt(13.0/pi)

        pre_y_6_6 = [i64 * sqrt_3003_pi * val for val in self._pow_sin_t[6]]
        pre_y_6_5 = [i32_3 * sqrt_1001_pi * val[0] * val[1] \
            for val in zip(self._pow_sin_t[5], self._pow_cos_t[1])]
        pre_y_6_4 = [i32_3 * sqrt_91_2pi * val[0] * (11.0*val[1] - 1.0) \
            for val in zip(self._pow_sin_t[4], self._pow_cos_t[2])]
        pre_y_6_3 = [i32 * sqrt_1365_pi * val[0] * (11.0*val[1] - 3.0*val[2]) \
            for val in zip(self._pow_sin_t[3], self._pow_cos_t[3], \
            self._pow_cos_t[1])]
        pre_y_6_2 = [i64 * sqrt_1365_pi * val[0] * (33.0*val[1] - \
            18.0*val[2] + 1.0) for val in zip(self._pow_sin_t[2], \
            self._pow_cos_t[4], self._pow_cos_t[2])]
        pre_y_6_1 = [i16 * sqrt_273_2pi * val[0] * (33.0*val[1] - \
            30.0*val[2] + 5.0*val[3]) for val in zip(self._pow_sin_t[1], \
            self._pow_cos_t[5], self._pow_cos_t[3], self._pow_cos_t[1])]

        acc = 0.0

        # Y_6_-6
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i] # cos(x) =  cos(-x)
            imag -= pre_y_6_6[i] * self._sin_n_p[6][i] # sin(x) = -sin(-x)
        acc += (real*real + imag*imag)

        # Y_6_-5
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real*real + imag*imag)

        # Y_6_-4
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag -= pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real*real + imag*imag)

        # Y_6_-3
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real*real + imag*imag)

        # Y_6_-2
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag -= pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        # Y_6_-1
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_6_0
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += i32 * sqrt_13_pi * (231.0*self._pow_cos_t[6][i] - \
                315.0*self._pow_cos_t[4][i] + 105.0*self._pow_cos_t[2][i] - 5.0)
        acc += (real*real)

        # Y_6_1
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_1[i] * self._cos_n_p[1][i]
            imag -= pre_y_6_1[i] * self._sin_n_p[1][i]
        acc += (real*real + imag*imag)

        # Y_6_2
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_2[i] * self._cos_n_p[2][i]
            imag += pre_y_6_2[i] * self._sin_n_p[2][i]
        acc += (real*real + imag*imag)

        # Y_6_3
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_3[i] * self._cos_n_p[3][i]
            imag -= pre_y_6_3[i] * self._sin_n_p[3][i]
        acc += (real*real + imag*imag)

        # Y_6_4
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_4[i] * self._cos_n_p[4][i]
            imag += pre_y_6_4[i] * self._sin_n_p[4][i]
        acc += (real*real + imag*imag)

        # Y_6_5
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real -= pre_y_6_5[i] * self._cos_n_p[5][i]
            imag -= pre_y_6_5[i] * self._sin_n_p[5][i]
        acc += (real*real + imag*imag)

        # Y_6_6
        real = imag = 0.0
        real = 0.0
        imag = 0.0
        for i in nnn_range:
            real += pre_y_6_6[i] * self._cos_n_p[6][i]
            imag += pre_y_6_6[i] * self._sin_n_p[6][i]
        acc += (real*real + imag*imag)

        q6 = math.sqrt(4.0 * pi * acc / (13.0 * float(nnn*nnn)))
        return q6


    def get_type(self, index):

        """
        Return type of order-parameter at the index provided and
        represented by a short string.

        Args:
            index (int):
                index of order-parameter for which type is to be returned
        """
        if index < 0 or index >= len(self._types):
            raise ValueError("Index for getting order-parameter type"
                             " out-of-bounds!")
        return self._types[index]


    def get_parameters(self, index):

        """
        Returns list of floats that represents
        the parameters associated with calculation of the order
        parameter that was defined at the index provided.
        Attention: the parameters do not need to equal those originally
        inputted because of processing out of efficiency reasons.

        Args:
            index (int):
                index of order-parameter for which associated parameters
                are to be returned
        """
        if index < 0 or index >= len(self._types):
            raise ValueError("Index for getting parameters associated with"
                             " order-parameter calculation out-of-bounds!")
        return self._paras[index]


    def get_order_parameters(self, structure, n, indeces_neighs=[], \
            tol=0.0, target_spec=None):

        """
        Compute all order parameters of site n.

        Args:
            structure (Structure):
                input structure.
            n (int):
                index of site in input structure, for which OPs are to be
                calculated.  Note that we do not use the sites iterator
                here, but directly access sites via struct[index].
            indeces_neighs ([int]):
                list of indeces of those neighbors in Structure object
                structure that are to be considered for OP computation.
                This optional argument overwrites the way neighbors are
                to be determined as defined in the constructor (i.e.,
                Voronoi coordination finder via negative cutoff radius
                vs constant cutoff radius if cutoff was positive).
                We do not use information about the underlying
                structure lattice if the neighbor indeces are explicitly
                provided.  This has two important consequences.  First,
                the input Structure object can, in fact, be a
                simple list of Site objects.  Second, no nearest images
                of neighbors are determined when providing an index list.
                Note furthermore that this neighbor
                determination type ignores the optional target_spec
                argument.
            tol (float):
                threshold of weight (= solid angle / maximal solid angle)
                to determine if a particular pair is
                considered neighbors; this is relevant only in the case
                when Voronoi polyhedra are used to determine coordination
            target_spec (Specie):
                target specie to be considered when calculating the order
                parameters of site n; None includes all species of input
                structure.

        Returns:
            list of floats representing order parameters.  Should it not be
            possible to compute a given OP for a conceptual reason, the
            corresponding entry is None instead of a float.  For Steinhardt
            et al.'s bond orientational OPs and the other geometric OPs
            ("tet", "oct", "bcc"), this can happen if there is a single
            neighbor around site n in the structure because that, obviously,
            does not permit calculation of angles between multiple
            neighbors.
        """

        # Do error-checking and initialization.
        if n < 0:
            raise ValueError("Site index smaller zero!")
        if n >= len(structure):
            raise ValueError("Site index beyond maximum!")
        if indeces_neighs:
            for index in indeces_neighs:
                if index >= len(structure):
                    raise ValueError("Neighbor site index beyond maximum!")
        if tol < 0.0:
            raise ValueError("Negative tolerance for weighted solid angle!")
        left_of_unity = 1.0 - 1.0e-12

        # Find central site and its neighbors.
        # Note that we adopt the same way of accessing sites here as in
        # VoronoiCoordFinder; that is, not via the sites iterator.
        centsite = structure[n]
        if indeces_neighs:
            neighsites = [structure[index] for index in indeces_neighs]
        elif self._voroneigh:
            vorocf = VoronoiCoordFinder(structure)
            neighsites = vorocf.get_coordinated_sites(n, tol, target_spec)
        else:
            # Structure.get_sites_in_sphere --> also other periodic images
            neighsitestmp = [i[0] for i in structure.get_sites_in_sphere(
                    centsite.coords, self._cutoff)]
            neighsites = []
            if centsite not in neighsitestmp:
                raise ValueError("Could not find center site!")
            else:
                neighsitestmp.remove(centsite)
            if target_spec is None:
                neighsites = list(neighsitestmp)
            else:
                neighsites[:] = [site for site in neighsitestmp \
                        if site.specie.symbol == target_spec]
        nneigh = len(neighsites)
        self._last_nneigh = nneigh

        # Prepare angle calculations, if applicable.
        rij = []
        rijnorm = []
        dist = []
        centvec = centsite.coords
        if self._computerijs:
            for j, neigh in enumerate(neighsites):
                rij.append((neigh.coords - centvec))
                dist.append(np.linalg.norm(rij[j]))
                rijnorm.append((rij[j]/dist[j]))

        # Initialize OP list and, then, calculate OPs.
        ops = [0.0 for t in self._types]

        # First, coordination number-based OPs.
        for i, t in enumerate(self._types):
            if t == "cn":
                ops[i] = nneigh / self._paras[i][0]

        # Then, bond orientational OPs based on spherical harmonics
        # according to Steinhardt et al., Phys. Rev. B, 28, 784-805, 1983.
        if self._boops:
            thetas = []
            phis = []
            for j, vec in enumerate(rijnorm):

                # z is North pole --> theta between vec and (0, 0, 1)^T.
                # Because vec is normalized, dot product is simply vec[2].
                thetas.append(math.acos(inbounds(-1.0, vec[2], 1.0)))
                tmpphi = 0.0

                # Compute phi only if it is not (almost) perfectly
                # aligned with z-axis.
                if vec[2] < left_of_unity and vec[2] > - (left_of_unity):
                    # x is prime meridian --> phi between projection of vec
                    # into x-y plane and (1, 0, 0)^T
                    tmpphi = math.acos(inbounds(
                            -1.0,
                            vec[0]/(math.sqrt(vec[0]*vec[0]+vec[1]*vec[1])),
                            1.0))
                    if vec[1] < 0.0:
                        tmpphi = -tmpphi
                phis.append(tmpphi)

            # Note that None flags that we have too few neighbors
            # for calculating BOOPS.
            for i, t in enumerate(self._types):
                if t == "q2":
                    ops[i] = self.get_q2(thetas, phis) if len(thetas) > 0 else None
                elif t == "q4":
                    ops[i] = self.get_q4(thetas, phis) if len(thetas) > 0 else None
                elif t == "q6":
                    ops[i] = self.get_q6(thetas, phis) if len(thetas) > 0 else None

        # Then, deal with the Shetty--Peters style OPs that are tailor-made
        # to recognize common structural motifs
        # (Shetty et al., J. Chem. Phys., 117, 4000-4009, 2002;
        #  Peters, J. Chem. Phys., 131, 244103, 2009;
        #  Zimmermann et al., J. Am. Chem. Soc., under revision, 2015).
        if self._geomops:
            gaussthetak = [0.0 for t in self._types] # not used by all OPs
            ipi = 1.0 / pi
            piover2 = pi / 2.0
            tetangoverpi = math.acos(-1.0 / 3.0) * ipi
            itetangminuspihalfoverpi = 1.0/(tetangoverpi-0.5)

            for j in range(nneigh): # Neighbor j is put to the North pole.
                zaxis = rijnorm[j]

                for k in range(nneigh): # From neighbor k, we construct
                    if j != k:          # the prime meridian.
                        tmp = inbounds(
                              -1.0, np.inner(zaxis, rijnorm[k]), 1.0)
                        thetak = math.acos(tmp)
                        xaxistmp = gramschmidt(rijnorm[k], zaxis)
                        xaxis = normalize(xaxistmp)
                        flag_xaxis = False
                        if vlength(xaxis) == 0.0:
                            flag_xaxis = True

                        # Contributions of j-i-k angles, where i represents the central atom
                        # and j and k two of the neighbors.
                        for i, t in enumerate(self._types):
                            if t == "tet":
                                tmp = self._paras[i][0] * (
                                      thetak * ipi - tetangoverpi)
                                gaussthetak[i] = math.exp(-0.5*tmp*tmp)
                            elif t == "oct":
                                if thetak >= self._paras[i][0]: # k is south pole to j
                                    tmp = self._paras[i][1] * (
                                          thetak * ipi - 1.0)
                                    ops[i] += 3.0 * math.exp(-0.5 * tmp*tmp)
                            elif t == "bcc" and j < k:
                                if thetak >= self._paras[i][0]: # k is south pole to j
                                    tmp     = self._paras[i][1] * (
                                              thetak * ipi - 1.0)
                                    ops[i] += 6.0 * math.exp(-0.5 * tmp*tmp)

                        for m in range(nneigh):
                            if (m != j) and (m != k):
                                tmp = inbounds(
                                    -1.0, np.inner(zaxis, rijnorm[m]), 1.0)
                                thetam = math.acos(tmp)
                                xtwoaxistmp = gramschmidt(rijnorm[m],zaxis)
                                xtwoaxis = normalize(xtwoaxistmp)
                                phi = math.acos(inbounds(
                                    -1.0, np.inner(xtwoaxis, xaxis), 1.0))
                                flag_xtwoaxis = False
                                if vlength(xtwoaxis) == 0.0:
                                    flag_xtwoaxis = True

                                # Contributions of j-i-m angle and
                                # angles between plane j-i-k and i-m vector.
                                if not flag_xaxis and not flag_xtwoaxis:
                                    for i, t in enumerate(self._types):
                                        if t == "tet":
                                            tmp     = self._paras[i][0] * (
                                                      thetam*ipi-tetangoverpi)
                                            ops[i] += gaussthetak[i] *math.exp(
                                                      -0.5*tmp*tmp) * math.cos(
                                                      3.0*phi)
                                        elif t == "oct":
                                            if thetak < self._paras[i][0] and \
                                                    thetam < self._paras[i][0]:
                                                tmp     = math.cos(2.0 * phi)
                                                tmp2    = self._paras[i][2] * (
                                                          thetam * ipi - 0.5)
                                                ops[i] += tmp*tmp * \
                                                          self._paras[i][4] * (
                                                          math.exp(
                                                          -0.5*tmp2*tmp2) - \
                                                          self._paras[i][3])
                                        elif t == "bcc" and j < k:
                                            if thetak<self._paras[i][0]:
                                                if thetak > piover2:
                                                    fac = 1.0
                                                else:
                                                    fac = -1.0
                                                tmp = (thetam - piover2) / (19.47 * pi / 180.0)
                                                ops[i] += fac * math.cos(3.0*phi)* \
                                                          1.6 * tmp * \
                                                          math.exp(-0.5*tmp*tmp)

            # Normalize Shetty--Peters style OPs.
            for i, t in enumerate(self._types):
                if t == "tet":
                    ops[i] = ops[i] / 24.0 if nneigh > 2 else None

                elif t == "oct":
                    ops[i] = ops[i] / 90.0 if nneigh > 2 else None

                elif t == "bcc":
                    # Reassured 144 by pen-and-paper for perfect bcc
                    # structure (24 from i-j-k South pole contributions
                    # and 4 * 6 * 5 = 120 from i-j-k-m non-South
                    # pole-containing quadruples.
                    # Because qbcc has a separate South pole contribution
                    # as summand, it is sufficient to have 2 neighbors to
                    # obtain a configuration that is well-defined by
                    # yielding a contribution to qbcc via thetak alone:
                    # ==> nneigh > 1.
                    ops[i] = ops[i] / 144.0 if nneigh > 1 else None

        return ops

