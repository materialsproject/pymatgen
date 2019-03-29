from __future__ import division, unicode_literals

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

__author__ = "Nicholas Winner"
__copyright__ = "None"
__version__ = "0.2"
__maintainer__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"
__status__ = "Development"
__date__ = ""

from pymatgen.analysis.dynamics.utils import autocorrelation, xdatcar_to_df, ensemble_average

boltzmann = 8.6173303e-5       # eV/K
plank     = 4.135667e-15       # eV/s
plank_ps  = plank / 1e12       # eV/ps
hbar      = plank/(2*np.pi)    # eV/s
hbar_ps   = plank_ps/(2*np.pi) # eV/ps


# TODO: Lots


class Vibrations(object):

    def __init__(self, data, from_positions=False):
        if from_positions:
            self.data = np.diff(data, axis=0)
        else:
            self.data = data

        self.vdos_spectrum  = None
        self.ir_spectrum    = None
        self.raman_spectrum = None

    def calc_vdos_spectrum(self):
        v = self.data
        print(v)
        exit()
        vv = np.concatenate((v, v))  # Mirror the data (FFT assumes periodicity)

        vacf = autocorrelation(vv)

        self.vdos_spectrum = np.fft.ifft(vacf)

        time = np.arange(0, 10000)
        freq = 1 / time
        wavenumber = freq / 3e8
        plt.plot(self.vdos_spectrum[0:1000])
        plt.show()

    def calc_ir_spectrum(self):

        dacf = autocorrelation(d)

        A = np.fft.fft(dacf)

        # TODO: Create pre-factor
        prefactor = 1

        self.ir_spectrum = prefactor*A

    def calc_raman_spectrum(self):
        return

    def calc_spectra(self):

        self.calc_power_spectrum()
        self.calc_ir_spectrum()
        self.calc_raman_spectrum()

    def plot(self):
        plt.plot(self.power_spectrum)
        plt.show()
        return

    def to_csv(self):
        return


df = pd.read_csv('/Users/nwinner/code/venv/xdatcar.csv.gz', compression='infer')


def func(v):
    vv = np.diff(v, axis=0)
    mir = np.concatenate((vv[::-1], vv[1:]))  # Mirror the data (FFT assumes periodicity)
    z = np.concatenate((mir, np.zeros((2*len(mir),3), dtype=float)))  # Pad the data with 0's (FFT assumes periodicity)
    return autocorrelation(z[0:len(vv)])


positions = ensemble_average(df, func, values=['x', 'y', 'z'], types=['F'])

s = np.fft.fft(positions)

plt.plot(positions)
plt.show()

time = df['Timestep'].drop_duplicates().values[0:-1]
time = np.multiply(time, .001)


