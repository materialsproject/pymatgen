from typing import Any, Dict, Optional

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender

from pymatgen.entries import Entry


__author__ = "Evan Spotte-Smith, Samuel Blau, Daniel Barter"
__copyright__ = "Copyright 2023, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__date__ = "02/04/2023"


class MoleculeEntry(Entry):
    """
    A class to provide easy access to Molecule properties.

    Args:
        molecule: pymatgen Molecule object
        energy: Electronic energy (units: Hartree)
        entry_id: Unique identifier for this MoleculeEntry (default: None)
        enthalpy: Enthalpy (units: kcal/mol) (default: None)
        entropy: Entropy (units: cal/mol-K) (default: None)
        mol_graph: graph representation of the Molecule (default: None)
        data: Other data associated with this entry (default: None)
    """

    def __init__(
        self,
        molecule: Molecule,
        energy: float,
        entry_id: object | None = None,
        enthalpy: Optional[float] = None,
        entropy: Optional[float] = None,
        mol_graph: Optional[MoleculeGraph] = None,
        data: Optional[Dict[str, Any]] = None
    ):
        self.molecule = molecule

        super().__init__(
            self.molecule.composition,
            self._energy
        )

        self.enthalpy = enthalpy
        self.entropy = entropy

        self.entry_id = entry_id

        if not mol_graph:
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            self.mol_graph = metal_edge_extender(mol_graph)
        else:
            self.mol_graph = mol_graph

        if data is None:
            self.data = dict()
        else:
            self.data = data

    @property
    def graph(self):
        return self.mol_graph.graph.to_undirected()

    @property
    def formula(self):
        return self.molecule.composition.alphabetical_formula

    @property
    def charge(self):
        return self.molecule.charge

    @property
    def spin_multiplicity(self):
        return self.molecule.spin_multiplicity

    @property
    def num_atoms(self):
        return len(self.molecule)

    @property
    def energy(self):
        return self._energy

    def free_energy(self, temperature: float = 298.15) -> float:
        """
        Gibbs free energy of the molecule (in eV) at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin. Default is 298.15, room temperature
        """

        if self.enthalpy is None or self.entropy is None:
            raise ValueError("Gibbs free energy requires enthalpy and entropy to be known!")

        u = self.energy * 27.2114
        h = self.enthalpy * 0.043363
        s = self.entropy * 0.000043363

        return u + h - temperature * s

    @classmethod
    def from_libe_madeira(
        cls,
        doc: Dict,
        use_thermo: str = "raw",
    ):
        """
        Initialize a MoleculeEntry from a document in the LIBE (Lithium-Ion
        Battery Electrolyte; Spotte-Smith*, Blau*, et al., Sci. Data 8(203), 2021
        ) or MADEIRA (MAgnesium Dataset of Electrolyte and
        Interphase ReAgents; manuscript in submission) datasets.

        Args:
            doc: Dictionary representing an entry from LIBE or MADEIRA
            use_thermo: One of "raw" (meaning raw, uncorrected thermo data will
                be used), "rrho_shifted" (meaning that a slightly modified
                Rigid-Rotor Harmonic Oscillator approximation will be used -
                see Ribiero et al., J. Phys. Chem. B 2011, 115, 14556-14562), or
                "qrrho" (meaning that Grimme's Quasi-Rigid Rotor Harmonic
                Oscillator - see Grimme, Chem. Eur. J. 2012, 18, 9955-9964) will
                be used.
        """

        thermo = use_thermo.lower()

        if thermo not in ["raw", "rrho_shifted", "qrrho"]:
            raise ValueError(
                "Only allowed values for use_thermo are 'raw', 'rrho_shifted', "
                "and 'qrrho'!"
            )
        try:
            if isinstance(doc["molecule"], Molecule):
                molecule = doc["molecule"]
            else:
                molecule = Molecule.from_dict(doc["molecule"])  # type: ignore

            if (
                thermo == "rrho_shifted"
                and doc["thermo"]["shifted_rrho_eV"] is not None
            ):
                energy = (
                    doc["thermo"]["shifted_rrho_eV"]["electronic_energy"] * 0.0367493
                )
                enthalpy = doc["thermo"]["shifted_rrho_eV"]["total_enthalpy"] * 23.061
                entropy = doc["thermo"]["shifted_rrho_eV"]["total_entropy"] * 23061
            elif thermo == "qrrho" and doc["thermo"]["quasi_rrho_eV"] is not None:
                energy = doc["thermo"]["quasi_rrho_eV"]["electronic_energy"] * 0.0367493
                enthalpy = doc["thermo"]["quasi_rrho_eV"]["total_enthalpy"] * 23.061
                entropy = doc["thermo"]["quasi_rrho_eV"]["total_entropy"] * 23061
            else:
                energy = doc["thermo"]["raw"]["electronic_energy_Ha"]
                enthalpy = doc["thermo"]["raw"]["total_enthalpy_kcal/mol"]
                entropy = doc["thermo"]["raw"]["total_entropy_cal/molK"]

            entry_id = doc["molecule_id"]

            if isinstance(doc["molecule_graph"], MoleculeGraph):
                mol_graph = doc["molecule_graph"]
            else:
                mol_graph = MoleculeGraph.from_dict(doc["molecule_graph"])

            data = dict()
            
            partial_charges_resp = doc.get('partial_charges', dict()).get('resp')
            if partial_charges_resp:
                data["partial_charges_resp"] = partial_charges_resp
            partial_charges_mulliken = doc.get('partial_charges', dict()).get('mulliken')
            if partial_charges_mulliken:
                data["partial_charges_mulliken"] = partial_charges_mulliken

            partial_charges_nbo = doc.get('partial_charges', dict()).get('nbo')
            if partial_charges_nbo:
                data["partial_charges_nbo"] = partial_charges_nbo
            partial_spins_nbo = doc.get('partial_spins', dict()).get('nbo')
            if partial_spins_nbo:
                data["partial_spins_nbo"] = partial_spins_nbo

            if 'redox' in doc:
                if 'electron_affinity_eV' in doc['redox']:
                    data["electron_affinity_eV"] = doc['redox']['electron_affinity_eV']

                if 'ionization_energy_eV' in doc['redox']:
                    data["ionization_energy_eV"] = doc['redox']['ionization_energy_eV']

        except KeyError as e:
            raise Exception(
                "Unable to construct molecule entry from molecule document; missing "
                f"attribute {e} in `doc`."
            )

        return cls(
            molecule=molecule,
            energy=energy,
            entry_id=entry_id,
            enthalpy=enthalpy,
            entropy=entropy,
            mol_graph=mol_graph,
            data=data
        )

    @classmethod
    def from_mpcule_summary_doc(d: Dict[str, Any],
                                solvent: str):
        """
        Initialize a MoleculeEntry from a summary document in the Materials Project
        molecules (MPcules) database.

        
        """

    @classmethod
    def from_dict(d: Dict[str, Any]):
        pass

    @classmethod
    def as_dict():
        pass

    def __repr__(self):

        output = [
            f"MoleculeEntry {self.entry_id} - {self.formula}",
            f"Charge = {self.charge}",
            f"Spin = {self.spin_multiplicity}"
        ]

        energies = [
            ("Energy", "Hartree", self.energy),
            ("Enthalpy", "kcal/mol", self.enthalpy),
            ("Entropy", "cal/mol.K", self.entropy),
            ("Free Energy (298.15 K)", "eV", self.free_energy()),
        ]
        for name, unit, value in energies:
            if value is None:
                output.append(f"{name} = {value} {unit}")
            else:
                output.append(f"{name} = {value:.4f} {unit}")

        if self.ind:
            output.append("Index: {}".format(self.ind))
        
        if len(self.data) > 0:
            output.append("Additional Data:")
        for k, v in self.data.items():
            output.append(f"\t{k} = {v}")

        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if type(self) == type(other):
            return str(self) == str(other)
        else:
            return False


class MoleculeCalculationEntry(Entry):
    pass
