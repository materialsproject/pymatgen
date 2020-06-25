import os

import numpy as np
import pandas as pd
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import plotly.express as px

# Load in the form factors

with open(os.path.join(os.path.dirname(__file__),
                       "xray_factors_2.csv")) as f:
    x_ray_scatter_df = pd.read_csv(f)
    # from https://it.iucr.org/Cb/ch6o1v0001/ table 6.1.1.4
with open(os.path.join(os.path.dirname(__file__),
                       "neutron_factors.csv")) as f:
    neutron_scatter_df = pd.read_csv(f)
    # from http://www.ccp14.ac.uk/ccp/web-mirrors/neutrons/n-scatter/n-lengths/LIST~1.HTM


class FStarDiagram:
    """
    Take a list of structure objects and used them to generate an f* phase diagram.
    """

    def __init__(self, structures, scattering_type='X-ray_simple', custom_scatter=None):
        """
        Initialize the f* diagram generator with the list of structures and scattering type.

        Args:
            structures(list): List of structure objects to use in the diagram.
            scattering_type(str): Type of scattering to use in the f* calculation. Defaults to 'X-ray_simple'
                which uses the atomic number as the scattering factor. 'X-ray' and 'Neutron' are built in scattering
                types which use X-ray and neutron scattering factors, respectively. 'Custom' allows the user to
                supplement their own calculation with any set of scattering factors.
            custom_scatter(function): when using custom scattering set this equal to a global variable that is equal
                to the custom scattering function.
        """

        self._structures = structures
        self._scatter = scattering_type
        self._custscat = custom_scatter
        self._symstructs = [SpacegroupAnalyzer(structure).get_symmetrized_structure() for structure in structures]
        self._equiv_inds = [struct.equivalent_indices for struct in self._symstructs]

    def get_site_labels(self):
        """
        Generates unique site labels based on composition, order, and symetry equivalence in the structure object.
        Ex:
        Structure Summary
        Lattice
            abc : 2.851 2.851 14.275
        angles : 90.0 90.0 119.99999999999999
        volume : 100.48498759501827
            A : 2.851 0.0 1.745734012184552e-16
            B : -1.4254999999999998 2.4690384261894347 1.745734012184552e-16
            C : 0.0 0.0 14.275
        PeriodicSite: Li:0.990, Ni:0.010 (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
        PeriodicSite: Li:0.990, Ni:0.010 (1.4255, 0.8230, 4.7583) [0.6667, 0.3333, 0.3333]
        PeriodicSite: Li:0.990, Ni:0.010 (-0.0000, 1.6460, 9.5167) [0.3333, 0.6667, 0.6667]
        PeriodicSite: Li:0.010, Mn:0.333, Co:0.333, Ni:0.323 (0.0000, 0.0000, 7.1375) [0.0000, 0.0000, 0.5000]
        PeriodicSite: Li:0.010, Mn:0.333, Co:0.333, Ni:0.323 (1.4255, 0.8230, 11.8958) [0.6667, 0.3333, 0.8333]
        PeriodicSite: Li:0.010, Mn:0.333, Co:0.333, Ni:0.323 (-0.0000, 1.6460, 2.3792) [0.3333, 0.6667, 0.1667]
        PeriodicSite: O (0.0000, 0.0000, 3.5688) [0.0000, 0.0000, 0.2500]
        PeriodicSite: O (0.0000, 0.0000, 10.7063) [0.0000, 0.0000, 0.7500]
        PeriodicSite: O (1.4255, 0.8230, 8.3270) [0.6667, 0.3333, 0.5833]
        PeriodicSite: O (1.4255, 0.8230, 1.1895) [0.6667, 0.3333, 0.0833]
        PeriodicSite: O (-0.0000, 1.6460, 13.0855) [0.3333, 0.6667, 0.9167]
        PeriodicSite: O (-0.0000, 1.6460, 5.9480) [0.3333, 0.6667, 0.4167]

        results in
        labels - ['0Li','1Co','2O']
        '0Li' - PeriodicSite: Li:0.990, Ni:0.010 (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
        '1Co' - PeriodicSite: Li:0.010, Mn:0.333, Co:0.333, Ni:0.323 (0.0000, 0.0000, 7.1375) [0.0000, 0.0000, 0.5000]
        '2O' - PeriodicSite: O (0.0000, 0.0000, 3.5688) [0.0000, 0.0000, 0.2500]
        """

        site_labels = []

        for ind, site in enumerate(self._equiv_inds[0]):
            label = str(ind) + [str(sp) for sp, occ in self._structures[0][site[0]].species_and_occu.items()][0]
            if label not in site_labels:
                site_labels.append(label)
        return site_labels

    def get_fstar_coords(self):

        """
        Calculate the f* coordinates for the list of structures.
        """

        fstar_lists = []

        for ind, struct in enumerate(self._equiv_inds):
            fstar_list = []
            for ind2, site in enumerate(struct):
                occ_f_list = []
                mult = len(site)
                elements_and_occupansies = self._symstructs[ind][site[0]].species_and_occu.items()
                for sp, occ in elements_and_occupansies:
                    if self._scatter == 'X-ray_simple':
                        f_occ = sp.Z * occ
                    if self._scatter == 'X-ray':
                        for i, n in enumerate(x_ray_scatter_df['atom'].values):
                            try:
                                if n == str(sp.element):
                                    f_occ = round(
                                        sum([x_ray_scatter_df.loc[i]['a1'], x_ray_scatter_df.loc[i]['a2'],
                                             x_ray_scatter_df.loc[i]['a3'],
                                             x_ray_scatter_df.loc[i]['a4'], x_ray_scatter_df.loc[i]['c']]), 0) * occ
                                    break
                                else:
                                    continue
                            except AttributeError:
                                if n == str(sp):
                                    f_occ = round(
                                        sum([x_ray_scatter_df.loc[i]['a1'], x_ray_scatter_df.loc[i]['a2'],
                                             x_ray_scatter_df.loc[i]['a3'],
                                             x_ray_scatter_df.loc[i]['a4'], x_ray_scatter_df.loc[i]['c']]), 0) * occ
                                    break
                                else:
                                    continue

                    if self._scatter == 'Neutron':
                        for i, n in enumerate(neutron_scatter_df['Isotope'].values):
                            try:
                                if n == str(sp.element):
                                    f_occ = float(neutron_scatter_df.loc[i]['Coh b']) * occ
                                    break
                                else:
                                    continue
                            except AttributeError:
                                if n == str(sp):
                                    f_occ = float(neutron_scatter_df.loc[i]['Coh b']) * occ
                                    break
                                else:
                                    continue
                    if self._scatter == 'Custom':
                        try:
                            f_occ = self._custscat(str(sp.element), occ, ind, ind2)
                        except AttributeError:
                            f_occ = self._custscat(str(sp), occ, ind, ind2)
                    occ_f_list.append(f_occ)

                fstar = np.absolute(mult * sum(occ_f_list))
                fstar_list.append(fstar)
            tot = sum(fstar_list)
            fstar_list = [fs / tot for fs in fstar_list]
            fstar_lists.append(fstar_list)
        return fstar_lists

    def plot_fstar_diagram(self, combine_list=False, plot_list=False, **kwargs):
        """
        Plot an f* diagram using plotly express.

        Args:
            combine_list(list): This is a list of lists which indicates what unique sites need to be combined to make
                the plot ternary.
            plot_list(list): This is a list that indicates what unique sites to plot and what order to plot them in.
            kwargs: use this to add any other arguments from scatter_ternary .
        """
        site_labels = FStarDiagram(self._structures, scattering_type=self._scatter,
                                   custom_scatter=self._custscat).get_site_labels()
        coords = FStarDiagram(self._structures, scattering_type=self._scatter,
                              custom_scatter=self._custscat).get_fstar_coords()
        df = pd.DataFrame(columns=site_labels, data=coords)
        df.replace('nan', 0.0, inplace=True)
        if combine_list:
            for combo in combine_list:
                df[str(combo)] = sum([df[site] for site in combo])
                site_labels.append(str(combo))
        if plot_list:
            return px.scatter_ternary(data_frame=df, a=plot_list[0], b=plot_list[1], c=plot_list[2], **kwargs)
        return px.scatter_ternary(data_frame=df, a=site_labels[0], b=site_labels[1], c=site_labels[2], **kwargs)
