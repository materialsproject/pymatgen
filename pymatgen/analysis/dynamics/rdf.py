from __future__ import division, unicode_literals

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy.signal import find_peaks

__author__ = "Nicholas Winner"
__copyright__ = "None"
__version__ = "0.2"
__maintainer__ = "Nicholas Winner"
__email__ = "nwinner@berkeley.edu"
__status__ = "Development"
__date__ = ""


# TODO: Add the properties you get from RDF as properties
# TODO: Add thermodynamic summary method


class RDF(object):
    """
    Class for representing the radial distribution function (RDF) of a structure / sequence of
    structures. Can be initialized either by providing an array of structures, [structure], and
    then running "calculate_rdf," or by reading a CSV file from a previous RDF instance.
    """

    def __init__(self, structures=None, data=None):
        self.structures = structures
        if data is not None:
            self.data = data
            test = pd.DataFrame({'density': 98/1186.97},index=[0])
            self.data = self.data.append(test)
            print(self.data)
            self.specie_pairs = data.columns.tolist()[1:]
        else:
            self.data = pd.DataFrame()
            self.specie_pairs = []

    def calculate_rdf(self, R, species, sampleSize=100, sampleFreq=10, total_rdf=False):
        '''
        Main function in the RDF class. Given a sequence of structures (i.e. as obtained from
        an XDATCAR or XYZ file from an MD run) this function returns the radial distribution
        function for all the pairs of species specified.

        Args:
            R: (int) Cut-off radius. The distance up to which the radial distribution function will
                  be calculated. Units are specified by the units of the structure list you provide;
                  so, it is generally in Angstroms. The function saves the RDF to the Pandas Dataframe
                  "self.df"
            species: ([ [Species1, Species2], [Species1, Species3 ], ...]) An array where each element
                       of the array is an array/tuple of length 2 that contains the pair of elements for
                       which to calculate the RDF. The species in the array can either be species objects or
                       string symbols of the elements.
                       Example: [ ['H','O'] ] will calculate the RDF between 'H' and 'O.'
            sampleSize: (int) The number of structures to consider for the RDF calculation. Generally, you only
                            want to sample the most recent structures from an MD run, because the most recent
                            structures represent the equilibrium configuration.
            sampleFreq: (int) Of the sample size, sample every "nth" structure.
                            Example: of samplesize = 100 and sampleFreq = 10, then sample structure 0, 9, 19, ... 99
                            in a reversed list of structures of range(0,100) where 0 is the last structure from
                            an MD run
            total_rdf:
        Return:
            None
        '''

        if self.structures is None:
            raise TypeError('Structures have not been defined for this RDF Object. Re-initialize with'
                            'a proper list of structures.')
        else:
            structures = self.structures

        structures = self.structures  # Put self.structures into a local instance for convenience
        structures = structures[len(structures)-sampleSize:-1]  # Use only the last "sampleSize" structures
        structures = structures[0::sampleFreq]  # Use every "sampleFreq"th structure from the sample

        # If species was provided as a single pair [X, Y] turn it into [ [X,Y] ]
        # Count the number of RDFs needed to be calculated and store in N_rdf
        if isinstance(species[0], list):
            N_rdf = len(species)
        else:
            N_rdf = 1
            species = [species]

        # If any of the species have been provided as "specie" objects, convert them to string symbols
        for i in range(len(species)):
            if not isinstance(species[i][0], str):
                species[i][0] = species[i][0].species_string
            if not isinstance(species[i][1], str):
                species[i][1] = species[i][1].species_string

        N     = len(structures[0].sites)  # The number of species in the structure

        n = {}  # A dict to store the number of atoms of each specie
        for i in structures[0].types_of_specie:
            n[i.symbol] = 0
        for j in structures[0].species:
            n[j.symbol] += 1

        rho   = N/structures[0].volume  # Density
        V     = structures[0].volume    # Volume of the structure (simulation cell size)
        dr    = .1                      # The width of the spherical shell to consider
        r     = np.arange(0.05, R, .01) # Distance values

        bin = [ [0]*len(r) for j in range(N_rdf) ]  # Bins of the RDF. One bin for each "r" value, and an array of
                                                    # These bins for each rdf needed to be calculated

        '''
        A function for binning the RDF values. Given The list of structures, and two species for which 
        to calculate the RDF, iterate pairwise over the species in each structure and determine if 
        it should be added to the RDF bins.
        '''
        def binning_function(specie1, specie2, r_values, delta_r, rdf_bin):
            for struc in structures:
                for x in range(len(r_values)):  # Iterate over the range of distances in array "r"
                    for i in range(len(struc.sites)):  # For all sites
                        if str(struc.sites[i].specie) is specie1:  # If the site in question is the specified specie 1
                            for j in range(len(struc.sites)):  # Compare i to all other sites
                                if str(struc.sites[j].specie) is specie2:  # if j is the specified specie 2
                                    dist = struc.sites[i].distance(struc.sites[j])  # Calculate the distance between i and j
                                    if dist < r_values[x] + delta_r and dist > r_values[x] - delta_r:  # If the distance is within r +- dr
                                        rdf_bin[x] += 1  # Add it to the bin for this r value

        # Run the binning function for each RDF needed to be calculated
        for i in range(N_rdf):
            binning_function(species[i][0], species[i][1], r, dr, bin[i])

        # Create "g" which is the RDF itself. It is also a PDF like the bins, but it is normalized
        # Equation for G(r), the RDF:
        #  <N(r+dr)>     1       <--Number of particles in sampling region (within the volume of the shell)
        # ----------- * ----
        #  <V(r+dr)>    rho      <---Size of the Sampling Region (volume of shell) * the number density of the species
        g = [[] for i in range(N_rdf)]
        for i in range(N_rdf):
            for j in range(len(bin[i])):
                number_of_particles_in_shell = bin[i][j]
                volume_of_shell = 4*np.pi*r[j]*r[j]*2*dr
                normalizer = n[ species[i][1] ]**2 / V
                g[i].append(number_of_particles_in_shell/(volume_of_shell*normalizer))

        rdfs = {}
        rdfs['r'] = r
        for i in range(N_rdf):
            label = str(species[i][0]) + '-' + str(species[i][1])
            self.specie_pairs.append(label)
            rdfs[label] = g[i]

        self.data = self.data.append(pd.DataFrame(rdfs))

    def peaks(self, pair):
        num_data_points = len(self.data[pair])
        return find_peaks(self.data[pair],distance=int(.1*num_data_points))

    def coordination_number(self, pair, shell_number=1):
        """
        Calculate the coordination number
            n = 4π rho integral( g(r) * r^2 dr)
        :param pair:
        :return:
        """
        peak = self.peaks(pair)[0][shell_number-1]
        return 4*np.pi*self.density*np.trapz(
            self.data[pair][0:peak]*self.data['r'][0:peak]**2,
            self.data['r'][0:peak])

    @property
    def density(self):
        return 70/1186.97

    @staticmethod
    def from_csv(filename):
        """
        Constructor for building an rdf object from an existing CSV file. Mainly useful
        if you want to quickly plot a previously calculated RDF.

        Args:
             filename: path + name of the csv file to read
        Returns:
            RDF object with the Dataframe loaded.
        """

        return RDF(data=pd.read_csv(filename))

    def plot(self, show=True, savefig=False, filename="rdf.png", smoothing=False):
        """
        Plotter function for the RDF object using Matplotlib. Can be used to either produce a plot instantly by using
        show=True, or to save a PNG file of the plot by using savefig=True.

        Args:
            show: (bool) Whether or not to display the plot after it is calculated to the screen. Defaults to True
            savefig: (bool) Whether or not to save the plot to a file. Defaults to False.
            filename: (str) filename to which to save the plot. Requires savefig=True. Defaults to "rdf.png"
            smoothing: (bool) Whether or not to use a spline smoothing on the RDF plot. Defaults to False
        Returns:
             None
        """

        fig, ax1 = plt.subplots()

        r = self.data['r']

        for pair in self.specie_pairs:
            g = self.data[pair]

            if smoothing:

                spline = interpolate.UnivariateSpline(r, g, k=3)
                r_smooth = np.linspace(r.min(), r.max(), 300)
                g_smooth = spline(r_smooth)

                g = g_smooth

            if smoothing:
                ax1.plot(r_smooth, g, label=pair)
            else:
                ax1.plot(r, g, label=pair)

                #TODO: Remove this or make it a plotting option
                peaks = self.peaks(pair)[0]
                for p in peaks:
                    ax1.vlines(x=r[p], ymin=0, ymax=10)

            ax1.minorticks_on()
            ax1.tick_params(which='major', length=8, width=1, direction='in', top=True, right=True, labelsize=18)
            ax1.tick_params(which='minor', length=2, width=.5, direction='in', top=True, right=True)
            ax1.set_xlabel("Seperation, r (Å)", size=30)
            ax1.set_ylabel("Radial Distribution Function, g(r)", size=30)

        ax1.legend(prop={'size': 20}, frameon=False)

        fig = plt.gcf()
        if savefig:
            fig.savefig(filename, format="png")
        if show:
            plt.show()

    def to_csv(self, filename="rdf.csv"):
        """
        Saves the current RDF Pandas Dataframe to a CSV file.
        Args:
            filename: Name of the file to write
        Returns:
            None
        """

        if self.data is None:
            raise ValueError("RDF Dataframe is of type 'None.' Did you forget to initialize it before "
                             "trying to write it to a CSV? Try calculating the RDF first with "
                             "the function calculate_rdf().")
        else:
            self.data.to_csv(filename, index=False)


rdf = RDF.from_csv("/Users/nwinner/code/venv/rdf.csv")


#rdf.plot(show=True, savefig=False)

print(rdf.coordination_number('Be-F', shell_number=1))
