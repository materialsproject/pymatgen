from collections import OrderedDict
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import curve_fit
import ruamel.yaml

from monty.serialization import loadfn

from pymatgen import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.reaction_calculator import ComputedReaction

"""
This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.
"""

def func(x, *m):
    """
    Helper function for curve_fit.
    """
    return np.dot(x, m)

class CorrectionCalculator:

    """
    A CorrectionCalculator contains experimental and computed entries which it uses to compute corrections.

    It graphs residual errors after applying the computed corrections and creates the MPCompatibility.yaml
    file the Correction classes use.

    Attributes:
        species: list of species that corrections are being calculated for
        exp_compounds: list of dictionaries which each contain a compound's formula and experimental data
        calc_compounds: dictionary of ComputedEntry objects
        corrections: list of corrections in same order as species list
        corrections_std_error: list of the variances of the corrections in same order as species list
        corrections_dict: dictionary of format {'species': (value, error)} for easier correction lookup
    """

    species = [
        "oxide",
        "peroxide",
        "superoxide",
        "F",
        "Cl",
        "Br",
        "I",
        "N",
        "S",
        "Se",
        "Si",
        "Sb",
        "Te",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Mo",
    ]  # species that we're fitting corrections for

    def __init__(self, exp_gz: str, comp_gz: str) -> None:

        """
        Initializes a CorrectionCalculator.

        Args:
            exp_gz: name of gzip file that contains experimental data
            comp_gz: name of gzip file that contains computed entries
        """

        self.exp_compounds = loadfn(exp_gz)  # experimental data
        self.calc_compounds = loadfn(comp_gz)  # computed entries
        self.corrections = []
        self.corrections_std_error = []
        self.corrections_dict = {}  # {'species': (value, error)}

        # these three lists are just to help the graph_residual_error_per_species() method
        self.oxides = []
        self.peroxides = []
        self.superoxides = []

    def compute_corrections(
        self,
        allow_polyanions: bool = False,
        allow_large_errors: bool = False,
        allow_unstable: bool = False,
    ) -> dict:

        """
        Computes the corrections and fills in correction, corrections_std_error, and corrections_dict.

        Args:
            allow_polyanions: optional variable, boolean controlling whether compounds with problematic polyanions
                will be included in the fit
            allow_large_errors: optional variable, boolean controlling whether compounds with large experimental
                uncertainties will be included in the fit
            allow_unstable: optional variable, boolean controlling whether unstable compounds with large e_above_hull
                will be included in the fit

        Raises:
            ValueError: calc_compounds is missing an entry
        """

        self.names = []
        self.diffs = []
        self.coeff_mat = []
        self.exp_uncer = []

        self.mpids = []
        for cmpd_info in self.exp_compounds:
            name = cmpd_info["formula"]
            warnings = cmpd_info["warnings"]

            if allow_polyanions:
                warnings.pop("polyanion", None)
            if allow_large_errors:
                warnings.pop("large_uncertainty", None)
            if allow_unstable:
                warnings.pop("unstable", None)

            if name in self.calc_compounds and not warnings:

                comp = Composition(name)
                elems = list(comp.as_dict())

                compound = self.calc_compounds[name]

                reactants = []
                for elem in elems:
                    try:
                        reactants.append(self.calc_compounds[elem])
                    except KeyError:
                        raise ValueError("Computed entries missing " + elem)

                rxn = ComputedReaction(reactants, [compound])
                rxn.normalize_to(comp)
                energy = rxn.calculated_reaction_energy

                if compound.data["oxide_type"] == "oxide":
                    coeff = [comp["O"], 0, 0]
                    self.oxides.append(name)
                elif compound.data["oxide_type"] == "peroxide":
                    coeff = [0, comp["O"], 0]
                    self.peroxides.append(name)
                elif compound.data["oxide_type"] == "superoxide":
                    coeff = [0, 0, comp["O"]]
                    self.superoxides.append(name)
                else:
                    coeff = [0, 0, 0]
                coeff += [comp[elem] for elem in self.species[3:]]

                self.names.append(name)
                self.diffs.append((cmpd_info["exp energy"] - energy) / comp.num_atoms)
                self.coeff_mat.append([i / comp.num_atoms for i in coeff])
                self.exp_uncer.append((cmpd_info["uncertainty"]) / comp.num_atoms)

                self.mpids.append(compound.entry_id)

        # for any exp entries with no uncertainty value, assign average uncertainty value
        sigma = np.array(self.exp_uncer)
        sigma[sigma == 0] = np.nan
        mean_uncer = np.nanmean(sigma)
        sigma = np.where(np.isnan(sigma), mean_uncer, sigma)

        popt, pcov = curve_fit(
            func,
            self.coeff_mat,
            self.diffs,
            p0=np.ones(21),
            sigma=sigma,
            absolute_sigma=True,
        )
        self.corrections = popt.tolist()
        self.corrections_std_error = np.sqrt(np.diag(pcov)).tolist()
        for i in range(len(self.species)):
            self.corrections_dict[self.species[i]] = (
                round(self.corrections[i], 3),
                round(self.corrections_std_error[i], 4),
            )
        return self.corrections_dict

    def graph_residual_error(self) -> None:

        """
        Graphs the residual errors for all compounds after applying computed corrections.
        """

        if len(self.corrections) == 0:
            self.compute_corrections()

        indices = [i for i in range(len(self.diffs))]
        abs_errors = [
            abs(i) for i in (self.diffs - np.dot(self.coeff_mat, self.corrections))
        ]
        labels_graph = self.names.copy()
        abs_errors, labels_graph = (
            list(t) for t in zip(*sorted(zip(abs_errors, labels_graph)))
        )  # sort by error

        num = len(abs_errors)
        fig = go.Figure(
            data=go.Scatter(
                x=np.linspace(1, num, num),
                y=abs_errors,
                mode="markers",
                text=labels_graph,
            ),
            layout=go.Layout(
                title=go.layout.Title(text="Residual Errors"),
                yaxis=go.layout.YAxis(
                    title=go.layout.yaxis.Title(text="Residual Error (eV/atom)")
                ),
            ),
        )
        fig.show()

        print("Residual Error:")
        print("Median = " + str(np.median(np.array(abs_errors))))
        print("Mean = " + str(np.mean(np.array(abs_errors))))
        print("Std Dev = " + str(np.std(np.array(abs_errors))))
        print("Original Error:")
        print("Median = " + str(abs(np.median(np.array(self.diffs)))))
        print("Mean = " + str(abs(np.mean(np.array(self.diffs)))))
        print("Std Dev = " + str(np.std(np.array(self.diffs))))

    def graph_residual_error_per_species(self, specie: str) -> None:

        """
        Graphs the residual errors for each compound that contains specie after applying computed corrections.

        Args:
            specie: the specie/group that residual errors are being plotted for

        Raises:
            ValueError: the specie is not a valid specie that this class fits corrections for
        """

        if specie not in self.species:
            raise ValueError("not a valid specie")

        if len(self.corrections) == 0:
            self.compute_corrections()

        abs_errors = [
            abs(i) for i in (self.diffs - np.dot(self.coeff_mat, self.corrections))
        ]
        labels_species = self.names.copy()
        diffs_cpy = self.diffs.copy()
        num = len(labels_species)

        if specie == "oxide" or specie == "peroxide" or specie == "superoxide":
            if specie == "oxide":
                compounds = self.oxides
            elif specie == "peroxide":
                compounds = self.peroxides
            else:
                compounds = self.superoxides
            for i in range(num):
                if labels_species[num - i - 1] not in compounds:
                    del labels_species[num - i - 1]
                    del abs_errors[num - i - 1]
                    del diffs_cpy[num - i - 1]
        else:
            for i in range(num):
                if not Composition(labels_species[num - i - 1])[specie]:
                    del labels_species[num - i - 1]
                    del abs_errors[num - i - 1]
                    del diffs_cpy[num - i - 1]
        abs_errors, labels_species = (
            list(t) for t in zip(*sorted(zip(abs_errors, labels_species)))
        )  # sort by error

        num = len(abs_errors)
        fig = go.Figure(
            data=go.Scatter(
                x=np.linspace(1, num, num),
                y=abs_errors,
                mode="markers",
                text=labels_species,
            ),
            layout=go.Layout(
                title=go.layout.Title(text="Residual Errors for " + specie),
                yaxis=go.layout.YAxis(
                    title=go.layout.yaxis.Title(text="Residual Error (eV/atom)")
                ),
            ),
        )
        fig.show()

        print("Residual Error:")
        print("Median = " + str(np.median(np.array(abs_errors))))
        print("Mean = " + str(np.mean(np.array(abs_errors))))
        print("Std Dev = " + str(np.std(np.array(abs_errors))))
        print("Original Error:")
        print("Median = " + str(abs(np.median(np.array(diffs_cpy)))))
        print("Mean = " + str(abs(np.mean(np.array(diffs_cpy)))))
        print("Std Dev = " + str(np.std(np.array(diffs_cpy))))

    def make_yaml(self, name: str = "MP") -> None:
        """
        Creates the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityErrors.yaml
        for correction errors.

        Args:
            name: optional argument, alternate name for the outputted yaml file
        """

        if len(self.corrections) == 0:
            self.compute_corrections()

        # from old mpcompatibility
        aqueous = OrderedDict()
        aqueous["O2"] = -0.316731
        aqueous["N2"] = -0.295729
        aqueous["F2"] = -0.313025
        aqueous["Cl2"] = -0.344373
        aqueous["Br"] = -0.235039
        aqueous["Hg"] = -0.234421
        aqueous["H2"] = -3.6018845
        aqueous["H2O"] = -4.972

        compatibility = OrderedDict()
        anion_corr = OrderedDict()
        advanced = OrderedDict()
        gas_corr = OrderedDict()
        u_corr = OrderedDict()
        o = OrderedDict()
        f = OrderedDict()

        compatibility_error = OrderedDict()
        anion_corr_error = OrderedDict()
        advanced_error = OrderedDict()
        gas_corr_error = OrderedDict()
        u_corr_error = OrderedDict()
        o_error = OrderedDict()
        f_error = OrderedDict()

        anion_corr["oxide"] = self.corrections_dict["oxide"][0]
        anion_corr["peroxide"] = self.corrections_dict["peroxide"][0]
        anion_corr["superoxide"] = self.corrections_dict["superoxide"][0]
        anion_corr["ozonide"] = 0  # do i need this??

        anion_corr_error["oxide"] = self.corrections_dict["oxide"][1]
        anion_corr_error["peroxide"] = self.corrections_dict["peroxide"][1]
        anion_corr_error["superoxide"] = self.corrections_dict["superoxide"][1]
        anion_corr_error["ozonide"] = 0  # do i need this??

        anion_corr["sulfide"] = self.corrections_dict["S"][0]
        anion_corr_error["sulfide"] = self.corrections_dict["S"][1]

        for elem in ["Br", "I", "Se", "Si", "Sb", "Te"]:
            anion_corr[elem] = self.corrections_dict[elem][0]
            anion_corr_error[elem] = self.corrections_dict[elem][1]

        for elem in ["F", "Cl", "N"]:
            entry = self.calc_compounds[elem]
            key = entry.composition.reduced_formula
            val = entry.energy_per_atom - self.corrections_dict[elem][0]
            gas_corr[key] = val
            gas_corr_error[key] = self.corrections_dict[elem][1]

        # from old mpcompatibility
        gas_corr["H2"] = -3.23973666138
        gas_corr_error["H2"] = 0

        for elem in ["V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Mo"]:
            o[elem] = self.corrections_dict[elem][0]
            f[elem] = self.corrections_dict[elem][0]

            o_error[elem] = self.corrections_dict[elem][1]
            f_error[elem] = self.corrections_dict[elem][1]

        u_corr["O"] = o
        u_corr["F"] = f
        advanced["UCorrections"] = u_corr
        compatibility["Name"] = name
        compatibility["Advanced"] = advanced
        compatibility["GasCorrections"] = gas_corr
        compatibility["AnionCorrections"] = anion_corr
        compatibility["AqueuousCompoundEnergies"] = aqueous

        u_corr_error["O"] = o_error
        u_corr_error["F"] = f_error
        advanced_error["UCorrections"] = u_corr_error
        compatibility_error["Name"] = name
        compatibility_error["Advanced"] = advanced_error
        compatibility_error["GasCorrections"] = gas_corr_error
        compatibility_error["AnionCorrections"] = anion_corr_error

        fn = name + "Compatibility.yaml"

        f = open(fn, "w")
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility, f)
        f.close()

        fn = name + "CompatibilityErrors.yaml"
        f = open(fn, "w")
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility_error, f)
        f.close()
