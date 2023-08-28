"""
This module implements an f* coordinate and diagram generator.
For information about f* diagrams, see the following publications:
Yin, L. et al. Extending the limits of powder diffraction analysis: diffraction parameter space, occupancy defects, and
atomic form factors. Rev. Sci. Instrum. 89, 093002 (2018). 10.1063/1.5044555
Yin, L. et al.Thermodynamics of Antisite Defects in Layered NMC Cathodes: Systematic Insights from High-Precision Powder
Diffraction Analyses Chem. Mater 2020 32 (3), 1002-1010. 10.1021/acs.chemmater.9b03646

"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import plotly.express as px

# Load in the neutron form factors
with open(os.path.join(os.path.dirname(__file__), "neutron_factors.csv")) as f:
    NEUTRON_SCATTER_DF = pd.read_csv(f)
    # from http://www.ccp14.ac.uk/ccp/web-mirrors/neutrons/n-scatter/n-lengths/LIST~1.HTM


class FStarDiagram:
    """
    Take a list of structure objects and/or cifs and use them to generate an f* phase diagram.
    """

    def __init__(self, structures, scattering_type="X-ray", custom_scatter=None):
        """
        Initialize the f* diagram generator with the list of structures and scattering type.

        Args:
            structures(list): List of structure objects to use in the diagram. These MUST be symmetrized structure
                objects.
            scattering_type(str): Type of scattering to use in the f* calculation. Defaults to 'X-ray'
                which uses the atomic number as the scattering factor. 'Neutron' is a built in scattering
                type which uses neutron scattering factors. 'Custom' allows the user to supplement their
                own calculation with any set of scattering factors.
            custom_scatter(function): when using custom scattering set this equal to a global varialble that is equal
                to the custom scattering function.
        """

        self._structures = structures
        self._scatter = scattering_type
        self._custscat = custom_scatter
        self._equiv_inds = [struct.equivalent_indices for struct in self._structures]
        self.site_labels = self.get_site_labels()
        self.coords = self.get_fstar_coords()
        self.plot = px.scatter_ternary(
            data_frame=self.coords, a=self.site_labels[0], b=self.site_labels[1], c=self.site_labels[2]
        )
        print("The labels for this structure's unique sites are")
        print(self.site_labels)

    def edit_fstar_diagram(self, combine_list=False, plot_list=False, **kwargs):
        """
        Edit the plot of the f* diagram using plotly express.

        Args:
            combine_list(list): This is a list of lists which indicates what unique sites need to be combined to make
                the plot ternary.
            plot_list(list): This is a list that indicates what unique sites or combined sites to plot and what order to
            plot them in.
            kwargs: use this to add any other arguments from scatter_ternary .
        """
        if combine_list:
            for combo in combine_list:
                self.coords[str(combo)] = sum([self.coords[site] for site in combo])
                if str(combo) not in self.site_labels:
                    self.site_labels.append(str(combo))
        if plot_list:
            self.plot = px.scatter_ternary(
                data_frame=self.coords, a=plot_list[0], b=plot_list[1], c=plot_list[2], **kwargs
            )
        else:
            self.plot = px.scatter_ternary(
                data_frame=self.coords, a=self.site_labels[0], b=self.site_labels[1], c=self.site_labels[2], **kwargs
            )

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
        labels - ['[0. 0. 0.]Li', '[0.  0.  0.5]Co', '[0.   0.   0.25]O']
        '[0. 0. 0.]Li' - PeriodicSite: Li:0.990, Ni:0.010 (0.0000, 0.0000, 0.0000) [0.0000, 0.0000, 0.0000]
        '[0.  0.  0.5]Co' - PeriodicSite: Li:0.010, Mn:0.333, Co:0.333, Ni:0.323
            (0.0000, 0.0000, 7.1375) [0.0000, 0.0000, 0.5000]
        '[0.   0.   0.25]O' - PeriodicSite: O (0.0000, 0.0000, 3.5688) [0.0000, 0.0000, 0.2500]
        """

        site_labels_fin = []
        for ind1, struct in enumerate(self._equiv_inds):
            site_labels = []
            for _ind2, site in enumerate(struct):
                label = str(self._structures[ind1][site[0]].frac_coords) + next(
                    str(sp) for sp, occ in self._structures[ind1][site[0]].species.items()
                )
                if label not in site_labels:
                    site_labels.append(label)
            if len(site_labels) > len(site_labels_fin):
                site_labels_fin = site_labels
        return site_labels_fin

    def get_fstar_coords(self):
        """
        Calculate the f* coordinates for the list of structures.
        """

        fstar_df_full = pd.DataFrame(columns=self.site_labels)

        for ind1, struct in enumerate(self._equiv_inds):
            fstar_df = pd.DataFrame(columns=self.site_labels, data=[[0.0 for i in self.site_labels]])
            for ind2, site in enumerate(struct):
                occ_f_list = []
                mult = len(site)
                site_frac_coord = str(self._structures[ind1][site[0]].frac_coords)
                column = [label for label in self.site_labels if site_frac_coord in label]
                elements_and_occupancies = self._structures[ind1][site[0]].species.items()
                for sp, occ in elements_and_occupancies:
                    if self._scatter == "X-ray":
                        f_occ = sp.Z * occ
                    if self._scatter == "Neutron":
                        for i, n in enumerate(NEUTRON_SCATTER_DF["Isotope"].values):
                            if hasattr(sp, "element"):
                                if n == str(sp.element):
                                    f_occ = float(NEUTRON_SCATTER_DF.loc[i]["Coh b"]) * occ
                                    break
                                else:
                                    continue
                            else:
                                if n == str(sp):
                                    f_occ = float(NEUTRON_SCATTER_DF.loc[i]["Coh b"]) * occ
                                    break
                                else:
                                    continue
                    if self._scatter == "Custom":
                        if hasattr(sp, "element"):
                            f_occ = self._custscat(str(sp.element), occ, ind1, ind2)
                        else:
                            f_occ = self._custscat(str(sp), occ, ind1, ind2)
                    occ_f_list.append(f_occ)

                fstar = np.absolute(mult * sum(occ_f_list))
                fstar_df.loc[0][column[0]] = round(float(fstar), 4)
            tot = sum(sum(list(fstar_df.values)))
            fstar_df = pd.DataFrame(columns=self.site_labels, data=[fs / tot for fs in list(fstar_df.values)])
            fstar_df_full = pd.concat([fstar_df_full, fstar_df], ignore_index=True)
        return fstar_df_full
