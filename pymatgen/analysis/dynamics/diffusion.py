from __future__ import division, unicode_literals

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import linregress
from pymatgen.io.xyz import XYZ


__author__ = "Nicholas Winner"
__copyright__ = "None"
__version__ = "0.2"
__maintainer__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"
__status__ = "Development"
__date__ = ""


from pymatgen.analysis.dynamics.utils import autocorrelation, xdatcar_to_df


class Diffusion(object):

    def __init__(self):
        self.diffusion_coefficients = {}


class MSD(object):

    def __init__(self, structures=None, data=None):
        self.structures = structures
        self.sites      = []
        self.t_convert  = 1
        self.data       = data

        self.msd = {}
        self.diffusion_coef = {}

        if structures:
            for i in structures:
                self.sites.append( i.sites )
            self.N          = len(self.sites[0])

    @staticmethod
    def from_csv(filename):
        data = pd.read_csv(filename, compression='infer')
        return MSD(data=data)

    @staticmethod
    def from_xdatcar(self, filename):
        data = xdatcar_to_df(filename, write_csv=False)
        return MSD(data=data)

    @staticmethod
    def form_xyz(self, filename):

        xyz = XYZ.from_file(filename)
        mols = xyz.all_molecules

        self.total_time = len(mols)
        for i in mols:
            self.sites.append( i.sites )
        self.N          = len(self.sites[0])
        self.time = np.arange(0, self.total_time * self.t_convert, self.t_convert)
        self.msd['time'] = self.time

    @staticmethod
    def from_dump(self, filename):
        data = dump_to_df(filename, write_csv=False)
        return MSD(data=data)

    def set_time(self, t_convert):
        self.data['Time'] = self.data['Timestep'](t_convert)

    def calculate_msd(self, specie, format=['x', 'y', 'z'], time_origins=5):

        """
        Calculate the MSD for many (all) values of tau. Supports averaging over many samples, e.g. you have the coordinate
        list of 10 particles in time.

        :param  specie: (str) Species for which to calculate the MSD
                format: [list] how the data is formatted in the dataframe. Usually it is simply [x,y,z], but it could
                        be ['xu', 'yu', 'zu'] if the unwrapped coordinates were gotten from LAMMPS for example.
        :return: [array] the MSD for all tau values
        """

        df  = self.data[self.data['type'] == specie]


        msd = []
        for i in df['Id']:
            coordinates = df[format].values
            time_origins = len(coordinates)
            msds = []
            for tau in range(time_origins):
                msd_sum  = 0
                num_msds = 0
                for t in range(len(coordinates)):
                    i_tau = t + tau
                    if i_tau < len(coordinates):
                        current_msd = sum((coordinates[t] - coordinates[i_tau]) ** 2)
                        msd_sum    += current_msd
                        num_msds   += 1
                    else:
                        break
                msds.append(
                    np.divide( msd_sum, num_msds)
                )
            msd.append(msds)

        self.msd[specie] = np.mean(msd, axis=0)

    def write_csv(self, filename, filename_d):
        df = pd.DataFrame(self.msd)
        df.to_csv(filename)

        df2 = pd.DataFrame(self.diffusion_coef)
        df2.to_csv(filename_d)

    def diffusion_coefficient(self, specie):
        D = [0, 0, 0, 0]
        trimmed_time = self.time[200:]
        trimmed_msd  = self.msd[specie][200:]

        q1 = round(.125*len(trimmed_time))
        q2 = round(.25*len(trimmed_time))
        q3 = round(.375*len(trimmed_time))
        q4 = round(.5*len(trimmed_time))
        q5 = round(.625*len(trimmed_time))
        q6 = round(.75*len(trimmed_time))
        q7 = round(.875*len(trimmed_time))

        D[0] = (1/6)*(linregress(trimmed_time[q1:q2], trimmed_msd[q1:q2])[0])*(1e-16)/(1e-15)
        D[1] = (1/6)*(linregress(trimmed_time[q3:q4], trimmed_msd[q3:q4])[0])*(1e-16)/(1e-15)
        D[2] = (1/6)*(linregress(trimmed_time[q5:q6], trimmed_msd[q5:q6])[0])*(1e-16)/(1e-15)
        D[3] = (1/6)*(linregress(trimmed_time[q7:],   trimmed_msd[q7:])[0])*(1e-16)/(1e-15)

        self.diffusion_coef[specie] = [np.mean(D)*(1e5)]

    def plot(self,specie):
        plt.plot(self.time, self.msd[specie])
        plt.show()


# TODO: Reading in data from files is ambiguous. Could either be reading input data, or reading precalculated data
# TODO: Add options for both and clear up the language
class VACF(object):

    def __init__(self, vdata={}, vacf={}):
        self.vdata = vdata
        self.vacf  = vacf

    @staticmethod
    def from_csv(self, filename):
        return VACF.__init__(vdata=pd.read_csv(filename, compression='infer'))

    @staticmethod
    def from_vdatcar(self,filename="VDATCAR", header=['vx', 'vy', 'vz']):
        return VACF.__init__(xdatcar_to_df(filename, header))

    def calc_vacf(self,specie):
        v = vdata[vdata['type'] == specie]
        self.vacf['specie'] = autocorrelation(v)

    def calc_diffusion_coefficient(self, specie):
        return np.trapz(self.vacf[specie])

    def to_csv(self, output):

        if self.vacf is False:
            raise ValueError("VACF Dataframe is empty. Did you forget to initialize it before "
                             "trying to write it to a CSV? Try calculating the VACF first with "
                             "the function calc_vacf().")
        else:
            self.vacf.to_csv(output, index=False)