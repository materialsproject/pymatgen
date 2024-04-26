"""This module calculates corrections for the species listed below, fitted to the experimental and computed
entries given to the CorrectionCalculator constructor.
"""

from __future__ import annotations

import os
import warnings

import numpy as np
import plotly.graph_objects as go
from monty.serialization import loadfn
from ruamel import yaml
from scipy.optimize import curve_fit

from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.analysis.structure_analyzer import sulfide_type
from pymatgen.core import Composition, Element


class CorrectionCalculator:
    """A CorrectionCalculator contains experimental and computed entries which it uses to compute corrections.

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

    def __init__(
        self,
        species: list[str] | None = None,
        max_error: float = 0.1,
        allow_unstable: float | bool = 0.1,
        exclude_polyanions: list[str] | None = None,
    ) -> None:
        """Initializes a CorrectionCalculator.

        Args:
            species: list of species to calculate corrections for
            max_error: maximum tolerable relative uncertainty in experimental energy.
                Compounds with relative uncertainty greater than this value will be excluded from the fit
            allow_unstable: whether unstable entries are to be included in the fit. If True, all compounds will
                be included regardless of their energy above hull. If False or a float, compounds with
                energy above hull greater than the given value (defaults to 0.1 eV/atom) will be
                excluded
            exclude_polyanions: a list of polyanions that contain additional sources of error that may negatively
                influence the quality of the fitted corrections. Compounds with these polyanions
                will be excluded from the fit
        """
        self.species = species or "oxide peroxide superoxide S F Cl Br I N Se Si Sb Te V Cr Mn Fe Co Ni W Mo H".split()

        self.max_error = max_error
        if not allow_unstable:
            self.allow_unstable = 0.1
        else:
            self.allow_unstable = allow_unstable
        self.exclude_polyanions = (
            exclude_polyanions
            if exclude_polyanions is not None
            else "SO4 SO3 CO3 NO3 NO2 OCl3 ClO3 ClO4 HO ClO SeO3 TiO3 TiO4 WO4 SiO3 SiO4 Si2O5 PO3 PO4 P2O7".split()
        )

        self.corrections: list[float] = []
        self.corrections_std_error: list[float] = []
        self.corrections_dict: dict[str, tuple[float, float]] = {}  # {'species': (value, uncertainty)}

        # to help the graph_residual_error_per_species() method differentiate between oxygen containing compounds
        if "oxide" in self.species:
            self.oxides: list[str] = []
        if "peroxide" in self.species:
            self.peroxides: list[str] = []
        if "superoxide" in self.species:
            self.superoxides: list[str] = []
        if "S" in self.species:
            self.sulfides: list[str] = []

    def compute_from_files(self, exp_gz: str, comp_gz: str):
        """
        Args:
            exp_gz: name of .json.gz file that contains experimental data
                    data in .json.gz file should be a list of dictionary objects with the following keys/values:
                    {"formula": chemical formula, "exp energy": formation energy in eV/formula unit,
                    "uncertainty": uncertainty in formation energy}
            comp_gz: name of .json.gz file that contains computed entries
                    data in .json.gz file should be a dictionary of {chemical formula: ComputedEntry}.
        """
        exp_entries = loadfn(exp_gz)
        calc_entries = loadfn(comp_gz)

        return self.compute_corrections(exp_entries, calc_entries)

    def compute_corrections(self, exp_entries: list, calc_entries: dict) -> dict:
        """Computes the corrections and fills in correction, corrections_std_error, and corrections_dict.

        Args:
            exp_entries: list of dictionary objects with the following keys/values:
                    {"formula": chemical formula, "exp energy": formation energy in eV/formula unit,
                    "uncertainty": uncertainty in formation energy}
            calc_entries: dictionary of computed entries, of the form {chemical formula: ComputedEntry}

        Raises:
            ValueError: calc_compounds is missing an entry
        """
        self.exp_compounds = exp_entries
        self.calc_compounds = calc_entries

        self.names: list[str] = []
        self.diffs: list[float] = []
        self.coeff_mat: list[list[float]] = []
        self.exp_uncer: list[float] = []

        # remove any corrections in calc_compounds
        for entry in self.calc_compounds.values():
            entry.correction = 0

        for cmpd_info in self.exp_compounds:
            # to get consistent element ordering in formula
            name = Composition(cmpd_info["formula"]).reduced_formula

            allow = True

            compound = self.calc_compounds.get(name)
            if not compound:
                warnings.warn(f"Compound {name} is not found in provided computed entries and is excluded from the fit")
                continue

            # filter out compounds with large uncertainties
            relative_uncertainty = abs(cmpd_info["uncertainty"] / cmpd_info["exp energy"])
            if relative_uncertainty > self.max_error:
                allow = False
                warnings.warn(
                    f"Compound {name} is excluded from the fit due to high experimental "
                    f"uncertainty ({relative_uncertainty:.1%})"
                )

            # filter out compounds containing certain polyanions
            for anion in self.exclude_polyanions:
                if anion in name or anion in cmpd_info["formula"]:
                    allow = False
                    warnings.warn(f"Compound {name} contains the poly{anion=} and is excluded from the fit")
                    break

            # filter out compounds that are unstable
            if isinstance(self.allow_unstable, float):
                try:
                    eah = compound.data["e_above_hull"]
                except KeyError:
                    raise ValueError("Missing e above hull data")
                if eah > self.allow_unstable:
                    allow = False
                    warnings.warn(f"Compound {name} is unstable and excluded from the fit (e_above_hull = {eah})")

            if allow:
                comp = Composition(name)
                elems = list(comp.as_dict())

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

                coeff = []
                for specie in self.species:
                    if specie == "oxide":
                        if compound.data["oxide_type"] == "oxide":
                            coeff.append(comp["O"])
                            self.oxides.append(name)
                        else:
                            coeff.append(0)
                    elif specie == "peroxide":
                        if compound.data["oxide_type"] == "peroxide":
                            coeff.append(comp["O"])
                            self.peroxides.append(name)
                        else:
                            coeff.append(0)
                    elif specie == "superoxide":
                        if compound.data["oxide_type"] == "superoxide":
                            coeff.append(comp["O"])
                            self.superoxides.append(name)
                        else:
                            coeff.append(0)
                    elif specie == "S":
                        if Element("S") in comp:
                            sf_type = "sulfide"
                            if compound.data.get("sulfide_type"):
                                sf_type = compound.data["sulfide_type"]
                            elif hasattr(compound, "structure"):
                                sf_type = sulfide_type(compound.structure)
                            if sf_type == "sulfide":
                                coeff.append(comp["S"])
                                self.sulfides.append(name)
                            else:
                                coeff.append(0)
                        else:
                            coeff.append(0)
                    else:
                        try:
                            coeff.append(comp[specie])
                        except ValueError:
                            raise ValueError(f"We can't detect this {specie=} in {name=}")

                self.names.append(name)
                self.diffs.append((cmpd_info["exp energy"] - energy) / comp.num_atoms)
                self.coeff_mat.append([i / comp.num_atoms for i in coeff])
                self.exp_uncer.append((cmpd_info["uncertainty"]) / comp.num_atoms)

        # for any exp entries with no uncertainty value, assign average uncertainty value
        sigma = np.array(self.exp_uncer)
        sigma[sigma == 0] = np.nan

        with warnings.catch_warnings():
            # numpy raises warning if the entire array is nan values
            warnings.simplefilter("ignore", category=RuntimeWarning)
            mean_uncert = np.nanmean(sigma)

        sigma = np.where(np.isnan(sigma), mean_uncert, sigma)

        if np.isnan(mean_uncert):
            # no uncertainty values for any compounds, don't try to weight
            p_opt, self.pcov = curve_fit(
                lambda x, *m: np.dot(x, m), self.coeff_mat, self.diffs, p0=np.ones(len(self.species))
            )
        else:
            p_opt, self.pcov = curve_fit(
                lambda x, *m: np.dot(x, m),
                self.coeff_mat,
                self.diffs,
                p0=np.ones(len(self.species)),
                sigma=sigma,
                absolute_sigma=True,
            )
        self.corrections = p_opt.tolist()
        self.corrections_std_error = np.sqrt(np.diag(self.pcov)).tolist()
        for idx, specie in enumerate(self.species):
            self.corrections_dict[specie] = (
                round(self.corrections[idx], 3),
                round(self.corrections_std_error[idx], 4),
            )

        # set ozonide correction to 0 so that this species does not receive a correction
        # while other oxide types do
        self.corrections_dict["ozonide"] = (0, 0)

        return self.corrections_dict

    def graph_residual_error(self) -> go.Figure:
        """Graphs the residual errors for all compounds after applying computed corrections."""
        if len(self.corrections) == 0:
            raise RuntimeError("Please call compute_corrections or compute_from_files to calculate corrections first")

        abs_errors = [abs(i) for i in self.diffs - np.dot(self.coeff_mat, self.corrections)]
        labels_graph = self.names.copy()
        abs_errors, labels_graph = (list(t) for t in zip(*sorted(zip(abs_errors, labels_graph))))  # sort by error

        n_err = len(abs_errors)
        fig = go.Figure(
            data=go.Scatter(
                x=np.linspace(1, n_err, n_err),
                y=abs_errors,
                mode="markers",
                text=labels_graph,
            ),
            layout=dict(title="Residual Errors", yaxis=dict(title="Residual Error (eV/atom)")),
        )

        print("Residual Error:")
        print(f"Median = {np.median(abs_errors)}")
        print(f"Mean = {np.mean(abs_errors)}")
        print(f"Std Dev = {np.std(abs_errors)}")
        print("Original Error:")
        print(f"Median = {abs(np.median(self.diffs))}")
        print(f"Mean = {abs(np.mean(self.diffs))}")
        print(f"Std Dev = {np.std(self.diffs)}")

        return fig

    def graph_residual_error_per_species(self, specie: str) -> go.Figure:
        """Graphs the residual errors for each compound that contains specie after applying computed corrections.

        Args:
            specie: the specie/group that residual errors are being plotted for

        Raises:
            ValueError: the specie is not a valid specie that this class fits corrections for
        """
        if specie not in self.species:
            raise ValueError("not a valid specie")

        if len(self.corrections) == 0:
            raise RuntimeError("Please call compute_corrections or compute_from_files to calculate corrections first")

        abs_errors = [abs(i) for i in self.diffs - np.dot(self.coeff_mat, self.corrections)]
        labels_species = self.names.copy()
        diffs_cpy = self.diffs.copy()
        n_species = len(labels_species)

        if specie in ("oxide", "peroxide", "superoxide", "S"):
            if specie == "oxide":
                compounds = self.oxides
            elif specie == "peroxide":
                compounds = self.peroxides
            elif specie == "superoxides":
                compounds = self.superoxides
            else:
                compounds = self.sulfides
            for idx in range(n_species):
                if labels_species[n_species - idx - 1] not in compounds:
                    del labels_species[n_species - idx - 1]
                    del abs_errors[n_species - idx - 1]
                    del diffs_cpy[n_species - idx - 1]
        else:
            for idx in range(n_species):
                if not Composition(labels_species[n_species - idx - 1])[specie]:
                    del labels_species[n_species - idx - 1]
                    del abs_errors[n_species - idx - 1]
                    del diffs_cpy[n_species - idx - 1]
        abs_errors, labels_species = (
            list(tup) for tup in zip(*sorted(zip(abs_errors, labels_species)))
        )  # sort by error

        n_err = len(abs_errors)
        fig = go.Figure(
            data=go.Scatter(
                x=np.linspace(1, n_err, n_err),
                y=abs_errors,
                mode="markers",
                text=labels_species,
            ),
            layout=dict(
                title=dict(text=f"Residual Errors for {specie}"),
                yaxis=dict(title="Residual Error (eV/atom)"),
            ),
        )

        print("Residual Error:")
        print(f"Median = {np.median(np.array(abs_errors))}")
        print(f"Mean = {np.mean(np.array(abs_errors))}")
        print(f"Std Dev = {np.std(np.array(abs_errors))}")
        print("Original Error:")
        print(f"Median = {abs(np.median(np.array(diffs_cpy)))}")
        print(f"Mean = {abs(np.mean(np.array(diffs_cpy)))}")
        print(f"Std Dev = {np.std(np.array(diffs_cpy))}")

        return fig

    def make_yaml(self, name: str = "MP2020", dir: str | None = None) -> None:
        """Creates the _name_Compatibility.yaml that stores corrections as well as _name_CompatibilityUncertainties.yaml
        for correction uncertainties.

        Args:
            name: str, alternate name for the created .yaml file.
                Default: "MP2020"
            dir: str, directory in which to save the file. Pass None (default) to
                save the file in the current working directory.
        """
        if len(self.corrections) == 0:
            raise RuntimeError("Please call compute_corrections or compute_from_files to calculate corrections first")

        # elements with U values
        ggau_correction_species = ["V", "Cr", "Mn", "Fe", "Co", "Ni", "W", "Mo"]

        comp_corr: dict[str, float] = {}
        o: dict[str, float] = {}
        f: dict[str, float] = {}

        comp_corr_error: dict[str, float] = {}
        o_error: dict[str, float] = {}
        f_error: dict[str, float] = {}

        for specie in [*self.species, "ozonide"]:
            if specie in ggau_correction_species:
                o[specie] = self.corrections_dict[specie][0]
                f[specie] = self.corrections_dict[specie][0]

                o_error[specie] = self.corrections_dict[specie][1]
                f_error[specie] = self.corrections_dict[specie][1]

            else:
                comp_corr[specie] = self.corrections_dict[specie][0]
                comp_corr_error[specie] = self.corrections_dict[specie][1]

        outline = """\
        Name:
        Corrections:
            GGAUMixingCorrections:
                O:
                F:
            CompositionCorrections:
        Uncertainties:
            GGAUMixingCorrections:
                O:
                F:
            CompositionCorrections:
        """
        fn = f"{name}Compatibility.yaml"
        path = os.path.join(dir, fn) if dir else fn

        yml = yaml.YAML()
        yml.default_flow_style = False
        contents = yml.load(outline)

        contents["Name"] = name

        # make CommentedMap so comments can be added
        contents["Corrections"]["GGAUMixingCorrections"]["O"] = yaml.comments.CommentedMap(o)
        contents["Corrections"]["GGAUMixingCorrections"]["F"] = yaml.comments.CommentedMap(f)
        contents["Corrections"]["CompositionCorrections"] = yaml.comments.CommentedMap(comp_corr)
        contents["Uncertainties"]["GGAUMixingCorrections"]["O"] = yaml.comments.CommentedMap(o_error)
        contents["Uncertainties"]["GGAUMixingCorrections"]["F"] = yaml.comments.CommentedMap(f_error)
        contents["Uncertainties"]["CompositionCorrections"] = yaml.comments.CommentedMap(comp_corr_error)

        contents["Corrections"].yaml_set_start_comment("Energy corrections in eV/atom", indent=2)
        contents["Corrections"]["GGAUMixingCorrections"].yaml_set_start_comment(
            "Composition-based corrections applied to transition metal oxides\nand fluorides to "
            'make GGA and GGA+U energies compatible\nwhen compat_type = "Advanced" (default)',
            indent=4,
        )
        contents["Corrections"]["CompositionCorrections"].yaml_set_start_comment(
            "Composition-based corrections applied to any compound containing\nthese species as anions",
            indent=4,
        )
        contents["Uncertainties"].yaml_set_start_comment(
            "Uncertainties corresponding to each energy correction (eV/atom)", indent=2
        )
        with open(path, mode="w") as file:
            yml.dump(contents, file)
