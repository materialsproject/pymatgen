"""
This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.
"""

from collections import OrderedDict
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import curve_fit
import ruamel.yaml
import warnings as w

from typing import List, Dict, Tuple
from monty.serialization import loadfn

from pymatgen import Element, Composition
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.analysis.structure_analyzer import sulfide_type


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
        corrections_dict: dictionary of format {'species': (value, uncertainty)} for easier correction lookup
    """

    species = [
        "oxide",
        "peroxide",
        "superoxide",
        "sulfide",
        "F",
        "Cl",
        "Br",
        "I",
        "N",
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
        "W",
        "Mo",
        "H",
    ]  # species that we're fitting corrections for

    def __init__(self, exp_gz: str, comp_gz: str) -> None:

        """
        Initializes a CorrectionCalculator.

        Args:
            exp_gz: name of gzip file that contains experimental data
                    data in gzip file should be a list of dictionary objects with the following keys/values:
                    {"formula": chemical formula, "exp energy": formation energy in eV/formula unit,
                    "uncertainty": uncertainty in formation energy, "warnings": {"polyanion": str describing
                    polyanions in compound or None, "large_uncertainty": str describing if compound has a large
                    uncertainty value or None, "unstable": str describing if compound has high e above hull or None}}
            comp_gz: name of gzip file that contains computed entries
                    data in gzip file should be a dictionary of {chemical formula: ComputedEntry}
        """

        self.exp_compounds = loadfn(exp_gz)  # experimental data
        self.calc_compounds = loadfn(comp_gz)  # computed entries
        self.corrections: List[float] = []
        self.corrections_std_error: List[float] = []
        self.corrections_dict: Dict[
            str, Tuple[float, float]
        ] = {}  # {'species': (value, uncertainty)}

        # to help the graph_residual_error_per_species() method differentiate between oxygen containing compounds
        self.oxides: List[str] = []
        self.peroxides: List[str] = []
        self.superoxides: List[str] = []
        self.sulfides: List[str] = []

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

        self.names: List[str] = []
        self.diffs: List[float] = []
        self.coeff_mat: List[List[float]] = []
        self.exp_uncer: List[float] = []

        # remove any corrections in calc_compounds
        for entry in self.calc_compounds.values():
            entry.correction = 0

        for cmpd_info in self.exp_compounds:

            # to get consistent element ordering in formula
            name = Composition(cmpd_info["formula"]).reduced_formula

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
                        elem_name = Composition(elem).reduced_formula
                        reactants.append(self.calc_compounds[elem_name])
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

                if Element("S") in comp:
                    sf_type = "sulfide"
                    if compound.data.get("sulfide_type"):
                        sf_type = compound.data["sulfide_type"]
                    elif hasattr(compound, "structure"):
                        sf_type = sulfide_type(compound.structure)
                    if sf_type == "sulfide":
                        coeff += [comp["S"]]
                        self.sulfides.append(name)
                    else:
                        coeff += [0]
                else:
                    coeff += [0]

                coeff += [comp[elem] for elem in self.species[4:]]

                self.names.append(name)
                self.diffs.append((cmpd_info["exp energy"] - energy) / comp.num_atoms)
                self.coeff_mat.append([i / comp.num_atoms for i in coeff])
                self.exp_uncer.append((cmpd_info["uncertainty"]) / comp.num_atoms)

        # for any exp entries with no uncertainty value, assign average uncertainty value
        sigma = np.array(self.exp_uncer)
        sigma[sigma == 0] = np.nan

        with w.catch_warnings():
            w.simplefilter(
                "ignore", category=RuntimeWarning
            )  # numpy raises warning if the entire array is nan values
            mean_uncer = np.nanmean(sigma)

        sigma = np.where(np.isnan(sigma), mean_uncer, sigma)

        if np.isnan(mean_uncer):
            # no uncertainty values for any compounds, don't try to weight
            popt, self.pcov = curve_fit(
                func, self.coeff_mat, self.diffs, p0=np.ones(len(self.species))
            )
        else:
            popt, self.pcov = curve_fit(
                func,
                self.coeff_mat,
                self.diffs,
                p0=np.ones(len(self.species)),
                sigma=sigma,
                absolute_sigma=True,
            )
        self.corrections = popt.tolist()
        self.corrections_std_error = np.sqrt(np.diag(self.pcov)).tolist()
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

        if (
            specie == "oxide"
            or specie == "peroxide"
            or specie == "superoxide"
            or specie == "sulfide"
        ):
            if specie == "oxide":
                compounds = self.oxides
            elif specie == "peroxide":
                compounds = self.peroxides
            elif specie == "superoxides":
                compounds = self.superoxides
            else:
                compounds = self.sulfides
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

    def make_yaml(self, name: str = "MP2020") -> None:
        """
        Creates the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityUncertainties.yaml
        for correction uncertainties.

        Args:
            name: str, alternate name for the created .yaml file.
            Default: "MP2020"
        """

        if len(self.corrections) == 0:
            self.compute_corrections()

        compatibility: OrderedDict = OrderedDict()
        comp_corr: OrderedDict[str, float] = OrderedDict()
        advanced: OrderedDict[str, OrderedDict] = OrderedDict()
        u_corr: OrderedDict[str, OrderedDict] = OrderedDict()
        o: OrderedDict[str, float] = OrderedDict()
        f: OrderedDict[str, float] = OrderedDict()

        compatibility_error: OrderedDict = OrderedDict()
        comp_corr_error: OrderedDict[str, float] = OrderedDict()
        advanced_error: OrderedDict[str, OrderedDict] = OrderedDict()
        u_corr_error: OrderedDict[str, OrderedDict] = OrderedDict()
        o_error: OrderedDict[str, float] = OrderedDict()
        f_error: OrderedDict[str, float] = OrderedDict()

        comp_corr["oxide"] = self.corrections_dict["oxide"][0]
        comp_corr["peroxide"] = self.corrections_dict["peroxide"][0]
        comp_corr["superoxide"] = self.corrections_dict["superoxide"][0]
        comp_corr["ozonide"] = 0  # do i need this??

        comp_corr_error["oxide"] = self.corrections_dict["oxide"][1]
        comp_corr_error["peroxide"] = self.corrections_dict["peroxide"][1]
        comp_corr_error["superoxide"] = self.corrections_dict["superoxide"][1]
        comp_corr_error["ozonide"] = 0  # do i need this??

        comp_corr["sulfide"] = self.corrections_dict["sulfide"][0]
        comp_corr_error["sulfide"] = self.corrections_dict["sulfide"][1]

        for elem in ["Br", "I", "Se", "Si", "Sb", "Te", "F", "Cl", "N", "H"]:
            comp_corr[elem] = self.corrections_dict[elem][0]
            comp_corr_error[elem] = self.corrections_dict[elem][1]

        for elem in ["V", "Cr", "Mn", "Fe", "Co", "Ni", "W", "Mo"]:
            o[elem] = self.corrections_dict[elem][0]
            f[elem] = self.corrections_dict[elem][0]

            o_error[elem] = self.corrections_dict[elem][1]
            f_error[elem] = self.corrections_dict[elem][1]

        u_corr["O"] = o
        u_corr["F"] = f
        advanced["UCorrections"] = u_corr
        compatibility["Name"] = name
        compatibility["Advanced"] = advanced
        compatibility["CompositionCorrections"] = comp_corr

        u_corr_error["O"] = o_error
        u_corr_error["F"] = f_error
        advanced_error["UCorrections"] = u_corr_error
        compatibility_error["Name"] = name
        compatibility_error["Advanced"] = advanced_error
        compatibility_error["CompositionCorrections"] = comp_corr_error

        fn = name + "Compatibility.yaml"
        file = open(fn, "w")
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility, file)
        file.close()

        fn = name + "CompatibilityUncertainties.yaml"
        file = open(fn, "w")
        yaml = ruamel.yaml.YAML()
        yaml.Representer.add_representer(OrderedDict, yaml.Representer.represent_dict)
        yaml.default_flow_style = False
        yaml.dump(compatibility_error, file)
        file.close()
