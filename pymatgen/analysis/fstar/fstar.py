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
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.express as px

if TYPE_CHECKING:
    from pymatgen.core.periodic_table import Element
    from pymatgen.core.structure import Structure

# Load in the neutron form factors
with open(f"{os.path.dirname(__file__)}/neutron_factors.csv.gz") as csv_file:
    NEUTRON_SCATTER_DF = pd.read_csv(csv_file)
    # from http://www.ccp14.ac.uk/ccp/web-mirrors/neutrons/n-scatter/n-lengths/LIST~1.HTM


class FStarDiagram:
    """
    Take a list of symmetrized structure objects and use them to generate an f* phase diagram.
    """

    def __init__(self, structures: list[Structure], scattering_type: str = "X-ray"):
        """
        Initialize the f* diagram generator with the list of structures and scattering type.

        Args:
            structures(list): List of structure objects to use in the diagram. These MUST be symmetrized structure
                objects. All structures must only vary in the occupancy of the sites.
            scattering_type(str): Type of scattering to use in the f* calculation. Defaults to 'X-ray'
                which uses the atomic number as the scattering factor. 'Neutron' is a built in scattering
                type which uses neutron scattering factors.
        """
        # check if the input structures list is valid
        for ind, struct in enumerate(structures):
            try:
                if struct.equivalent_indices == structures[ind - 1].equivalent_indices:
                    if len(struct.equivalent_indices) > 2:
                        continue
                    raise ValueError("Structures must have at least 3 unique sites")
                raise ValueError("All structures must only vary in occupancy.")
            except AttributeError:
                raise AttributeError("Must use symmeteized structure objects")
        self._structures = structures
        self._scatter = scattering_type
        self._scatter_dict = {"X-ray": self.xray_scatter, "Neutron": self.neutron_scatter}
        self.site_labels = self._get_site_labels()
        self.fstar_coords = self._get_fstar_coords()
        self.set_plot_list([self.site_labels[0], self.site_labels[1], self.site_labels[2]])
        self.make_plot()
        print("The labels for this structure's unique sites are")
        print(self.site_labels)

    def combine_sites(self, site_lists: list[list[str]]) -> None:
        """
        Many structures have more than three sites. If this is the case you may want to
            add some sites together to make a psudo-site.
        Args:
            site_lists(list): A list of lists of site label strings. This allows you to combine
                more than one set of sites at once.
        """
        for combo in site_lists:
            for site in combo:
                if site not in self.site_labels:
                    raise ValueError("All sites must be in the site_labels list")
            self.fstar_coords[str(combo)] = sum([self.fstar_coords[site] for site in combo])
            if str(combo) not in self.site_labels:
                self.site_labels.append(str(combo))

    def set_plot_list(self, site_list: list[str]) -> None:
        """
        set the list of sites to plot and the order to plot them in.
        Args:
            site_list(list): A list of site label strings. Index 0 goes on the top of the
                plot, index 1 goes on the bottom left, and index 2 goes on the bottom right.
        """
        for site in site_list:
            if site not in self.site_labels:
                raise ValueError("All sites must be in the site_labels list")
        self.plot_list = site_list

    def make_plot(self, **kwargs):
        """
        Makes a plotly express scatter_ternary plot using the fstar_coords dataframe and the
            sites in plot list.
        Args:
            **kwargs: this can be any argument that the scatter_ternary function can use.
        """
        self.plot = px.scatter_ternary(
            data_frame=self.fstar_coords, a=self.plot_list[0], b=self.plot_list[1], c=self.plot_list[2], **kwargs
        )

    def _get_site_labels(self):
        """
        Generates unique site labels based on composition, order, and symmetry equivalence in the structure object.
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
        site_labels = []
        for site in self._structures[0].equivalent_indices:
            site_labels.append(
                str(self._structures[0][site[0]].frac_coords)
                + next(str(sp) for sp, _ in self._structures[0][site[0]].species.items())
            )
        return site_labels

    def _get_fstar_coords(self):
        """
        Calculate the f* coordinates for the list of structures.
        """
        fstar_df_full = pd.DataFrame(columns=self.site_labels)
        for ind1, struct in enumerate(self._structures):
            fstar_df = pd.DataFrame(columns=self.site_labels, data=[[0.0 for i in self.site_labels]])
            for ind2, site in enumerate(struct.equivalent_indices):
                occ_f_list = []
                mult = len(site)
                site_frac_coord = str(self._structures[ind1][site[0]].frac_coords)
                column = [label for label in self.site_labels if site_frac_coord in label]
                elements_and_occupancies = self._structures[ind1][site[0]].species.items()
                for sp, occ in elements_and_occupancies:
                    # ind1 and ind2 are added in case someone wants to make a custom scatter function
                    # that uses information in the structure object
                    f_occ = self._scatter_dict[self._scatter](sp, occ, ind1, ind2)
                    occ_f_list.append(f_occ)
                fstar = np.absolute(mult * sum(occ_f_list))
                fstar_df.loc[0][column[0]] = round(float(fstar), 4)
            tot = sum(sum(list(fstar_df.values)))
            fstar_df = pd.DataFrame(columns=self.site_labels, data=[fs / tot for fs in list(fstar_df.values)])
            fstar_df_full = pd.concat([fstar_df_full, fstar_df], ignore_index=True)
        return fstar_df_full

    def xray_scatter(self, el: Element, occ: float, i1: int, i2: int) -> float:
        """
        X-ray scattering function. i2 and i2 are unused.
        """
        return el.Z * occ

    def neutron_scatter(self, el: Element, occ: float, i1: int, i2: int) -> float:
        """
        Neutron scattering function. i2 and i2 are unused.
        """
        for i, n in enumerate(NEUTRON_SCATTER_DF["Isotope"].values):
            if hasattr(el, "element"):
                if n == str(el.element):
                    f_occ = float(NEUTRON_SCATTER_DF.loc[i]["Coh b"]) * occ
                    break
            else:
                if n == str(el):
                    f_occ = float(NEUTRON_SCATTER_DF.loc[i]["Coh b"]) * occ
                    break
        return f_occ
