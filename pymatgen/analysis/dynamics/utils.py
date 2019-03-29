from __future__ import division, unicode_literals

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.signal import savgol_filter as savgol
from pymatgen.analysis.dynamics import constants

from pymatgen.io.lammps.outputs import parse_lammps_dumps

__author__ = "Nicholas Winner"
__copyright__ = "None"
__version__ = "0.2"
__maintainer__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"
__status__ = "Development"
__date__ = ""


def autocorrelation(v, mode=1, norm=False, detrend=False):
    """
    Using the function __autocorrelation__() this function calculates the autocorrelation for an input vector v. The
    only difference is that this function accounts for v being composed of vectors, e.g: v[0] = [v_x, v_y, v_z]

    Args:
        v: [array] Input vector, can be a 1D array of scalars, or an array of vectors.
        mode: (int) What mode to run the calculation (see __autocorrelation__()). 1 is for the FFT method and
                2 is for the discrete (long) calculation. Default=1
        norm: (Bool) Whether or not to normalize the result. Default=False
        detrend: (Bool) Whether or not to detrend the data, v = v - v_average. Default=False

    Returns:
        [Array] The autocorrelation data.
    """
    if isinstance(v[0], (list, tuple, np.ndarray)):
        transposed_v = list(zip(*v))
        acf = []
        for i in transposed_v:

            acf.append(__autocorrelation__(i, mode, norm, detrend))

        return np.mean(acf, axis=0)

    else:
        return __autocorrelation__(v, mode, norm, detrend)


def __autocorrelation__(v, mode=1, norm=False, detrend=False):
    """
    A simple function for performing the autocorrelation by either Wiener–Khinchin theorem using FFT, or by
    using the following expression:

            C(tau) = <v(t) . v(t+tau)> ; where <...> is ensemble averaged over time origins, t, and tau is the
            time delay.

    It is recommended to use the FFT method for speed unless you have a good reason not to do so. Generally, this
    function does not need to be implemented, because it is explicitely written for 1D arrays, and instead you
    should run autocorrelation(), which will decide if your data contains vectors, this function is the workhorse
    behind that function.

    Args:
        v: [array] the data for which to get the autocorrelation
        mode: (int) which mode to run this function. 1 is FFT based, and 2 follows the above formula.
                 default=1.
        norm: (bool) whether or not to normalize the ACF data to 1. Generally, for correlation functions in
                statistical mechanics, you do not want to normalize, because the area under the curve is
                connected to the transport coefficients (see, Green-Kubo formalism). Default=False
        detrend: (bool) Whether or not to de-trend the data by subtracting the mean from each element of v
    Returns:
         [array] the autocorrelation results
    """

    nstep = len(v)  # Number of time steps i.e. number of data points
    c = np.zeros((nstep,), dtype=float)  # empty array for the correlation data

    if detrend:
        v = np.subract(v, np.mean(v))

    # FFT Method
    if mode == 1:
        vv  = np.concatenate((v[::-1], v[1:]))  # Mirror the data (FFT assumes periodicity)
        vv0 = np.concatenate((vv, np.zeros((2*nstep,), dtype=float)))  # Pad the data with 0's (FFT assumes periodicity)
        c   = np.fft.ifft(np.abs(np.power( np.fft.fft(v), 2)))[:nstep].real  # Wiener–Khinchin, only take half mirror

    # Discrete Method
    if mode == 2:
        vv = np.concatenate((v, np.zeros((nstep,),dtype=float)))
        for t in range(nstep):
            for j in range(nstep):
                c[t] += v[j] * vv[j + t]

    c = c / nstep

    if norm:
        c = c / np.max(c)

    return c[0:int(len(v))]


def power_spectrum(v):

    """

    :param v:
    :return:
    """

    if isinstance(v[0], (list, tuple, np.ndarray)):
        transposed_v = list(zip(*v))
        acf = []
        for i in transposed_v:

            acf.append(__power_spectrum__(i))

        return np.mean(acf, axis=0)

    else:
        return __power_spectrum__(v)


def __power_spectrum__(v):

    """
    This function calculates the power spectral density for 1D array, v. It is very similar to the autocorrelation
    function, which is the more common function to use. The only difference in these is that there is no final inverse
    fast Fourier transform in the calculation, and that this function only supports the FFT method, and no "exact"
    method.

    :param v:
    :return:
    """

    nstep = len(v)  # Number of time steps i.e. number of data points
    c = np.zeros((nstep,), dtype=float)  # empty array for the correlation data

    vv = np.concatenate((v[::-1], v[1:]))  # Mirror the data (FFT assumes periodicity)
    vv0 = np.concatenate((vv, np.zeros((nstep,), dtype=float)))  # Pad the data with 0's (FFT assumes periodicity)
    c = np.abs(np.power(np.fft.fft(vv), 2))[:nstep].real  # Wiener–Khinchin, only take half mirror

    c = c / (2*nstep)

    return c[0:len(v)]


def dump_to_df(filename, write_csv=True, output="data.csv"):
    """
    A helper function that takes a lammps dump file and returns a Pandas Dataframe. It is also recommended to write
    a CSV file. This has the benefits:
        (1) CSV files take up less space than dump files.
        (2) It is far more efficient to parse massive amounts of lammps data entirely through Pandas instead of
            using a list of Pandas DataFrames, as Pymatgen currently provides, it just requires a little more knowledge
            of how to use Pandas efficiently. A single CSV can be read as a DF and then processed.
        (3) We can pre-sort the particles. LAMMPS does not retain the order of its particles in the dump file, which
            can be very annoying for post processing. When the csv is written, it sorts the Pd dataframe so that at
            each time step, the particles are listed in order of their id.

    Args:
        filename: (str) file name of the lammps dump.
        write_csv: (bool) Whether or not to write csv file
        output: (str) file name to output the csv. Include ".csv" in your filename.

    Returns:
        Pandas Dataframe of the dump
    """

    dump = parse_lammps_dumps(filename)
    dfs  = []
    for frame in dump:
        dfs.append(frame.data)
        dfs[-1]['Timestep'] = frame.timestep
    df = pd.concat(dfs).sort_values(by=['Timestep', 'id']).reset_index(drop=True)

    if write_csv:
        df.to_csv(output)

    return df


def xdatcar_to_df(filename="XDATCAR", write_csv=True, output="xdatcar.csv"):
    """
    Parses a DATCAR file into a Pandas Dataframe, and by default writes a csv file, which is quicker to read when
    doing further analysis. By default, it is equipped for reading XDATCAR, but the format is general. So, if you
    have a VDATCAR for velocities for example. You could read that as well.

    Args:
        filename: (str) Name of the file to read.
        header: [list] The values that are going to be read from the columns of the DATCAR file. Default='XDATCAR'
        write_csv: (bool) Whether or not to write the pandas dataframe to a csv file. Defaults to True.
        output: (str) Name of the csv file to write. Default='xdatcar.csv'
    Returns:
        Pandas Dataframe of the datcar
    """

    f = open(filename)
    lines = f.readlines()
    f.close()

    spec = lines[5].split()
    num  = lines[6].split()
    num = list(map(int, num))
    N    = 0

    species = []
    for i in range(len(num)):
        for j in range(num[i]):
            species.append(spec[i])

    N = sum(num) # Number of atoms/ions

    dfs = []
    j = 8
    timestep = 1

    while j+N+8 < len(lines):
        temp = [lines[i].split() for i in range(j, j+N)]
        dfs.append( pd.DataFrame().from_records(temp) )
        dfs[-1]['Timestep'] = lines[j-1].split()[-1]
        dfs[-1]['type']   = species
        j += N+8

    df = pd.concat(dfs).rename(columns={0: 'x', 1: 'y', 2: 'y'})
    df.index.names = ['Id']

    if write_csv:
        df.to_csv(output)

    return df


def ensemble_average(df, function, values, types):
    """
    For doing ensemble averaging. Using some function, it will get the ensemble average for all time, and all particle
    types that you specify

    Args:
        df: (Dataframe) Input data, needs to have the standard form with 'Id', 'type', and 'Timestep'
        function: (function) the function to whch we want to pass the ensemble data. I.e. an averaging
                    function, an autocorrelation function, a spectral density function, etc.
        values: The data that you want to evaluate with your function.
                Example: ['x', 'y', 'z']
        types: [array] Which particle types you want to evaluate.
                Example: ['H', 'O']
    :return:
    """
    ids = {}
    for i_type in types:
        ids[i_type] = df[(df.Timestep == 1) & (df.type == i_type)]['Id'].values

    results = []
    for key, value in ids.items():
        for id in value:
            data = df[df['Id'] == id][values].values
            results.append(function(data))

    return np.mean(results, axis=0)


def opt_filter(x, n):

    Nopt = 3
    N1   = 1

    while N1 != Nopt:

        N1 = int(np.floor(2*Nopt/2))
        if N1 % 2 == 0:
            N1 += 1
        print(N1)
        y  = savgol(x, N1, n)
        dy = savgol(np.diff(y,1), N1, n)
        y2 = np.diff(dy, 3)
        c1 = np.mean(np.power(y2,2))

        Nopt = np.power(((2*(n+2))*(np.power(np.math.factorial(2*n+3),2))*(np.var(x))) /
                        ((np.power(np.math.factorial(n+1),2))*(c1)), 1/(2*n+5))

    return Nopt

file = '/Users/nwinner/code/venv/dipoles.txt'

vel  = pd.read_csv('/Users/nwinner/code/venv/vdatcar.csv')

x = ensemble_average(vel, power_spectrum, ['vx','vy','vz'], ['Li','Be','F'])

x = np.power(np.abs(x), 2)/(3*98*973*constants.kB)


time = vel['Timestep'].drop_duplicates().values*.001
wavenumber = time * 100

plt.plot(wavenumber[0:10000], x[0:10000], label="")
y = savgol(x, 51, 2)
plt.plot(wavenumber[0:10000], y[0:10000], label="")

plt.legend()
plt.show()
exit()

with open(file) as f:

    lines = f.readlines()

    ionic      = []
    electronic = []


    for l in range(0, len(lines)-1, 2):

        line1 = lines[l]
        ionic_line = line1[line1.find('(')+1 : line1.find(')')]
        ionic_str = ionic_line.split()

        ionic.append([float(i) for i in ionic_str])


        line2 = lines[l+1]
        electronic_line = line2[line2.find('(')+1 : line2.find(')')]
        electronic_str = electronic_line.split()

        electronic.append([float(i) for i in electronic_str])


    total = np.add(ionic, electronic)

    v = autocorrelation(total)

    vv = np.concatenate((v[::-1], v[1:]))  # Mirror the data (FFT assumes periodicity)

    z = np.fft.fft(v)
    plt.plot(z[0:len(v)])
    plt.show()