from typing import Any, Dict, Optional

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender


class MoleculeEntry:
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

        self.energy = energy
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
    def free_energy():
        pass

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

    @classmethod
    def from_dataset_entry(
        cls,
        doc: Dict,
        use_thermo: str = "raw",
    ):
        """
        Initialize a MoleculeEntry from a document in the LIBE (Lithium-Ion
        Battery Electrolyte) or MADEIRA (MAgnesium Dataset of Electrolyte and
        Interphase ReAgents) datasets.

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

            partial_charges_resp = doc['partial_charges']['resp']
            partial_charges_mulliken = doc['partial_charges']['mulliken']
            spin_multiplicity = doc['spin_multiplicity']


            if doc['number_atoms'] == 1:
                partial_charges_nbo = doc['partial_charges']['mulliken']
                partial_spins_nbo = doc['partial_spins']['mulliken']
            else:
                partial_charges_nbo = doc['partial_charges']['nbo']
                partial_spins_nbo = doc['partial_spins']['nbo']

            electron_affinity_eV = None
            ionization_energy_eV = None
            if 'redox' in doc:
                if 'electron_affinity_eV' in doc['redox']:
                    electron_affinity_eV = doc['redox']['electron_affinity_eV']

                if 'ionization_energy_eV' in doc['redox']:
                    ionization_energy_eV = doc['redox']['ionization_energy_eV']

        except KeyError as e:
            raise Exception(
                "Unable to construct molecule entry from molecule document; missing "
                f"attribute {e} in `doc`."
            )



        return cls(
            molecule=molecule,
            energy=energy,
            enthalpy=enthalpy,
            entropy=entropy,
            entry_id=entry_id,
            mol_graph=mol_graph,
            partial_charges_resp=partial_charges_resp,
            partial_charges_mulliken=partial_charges_mulliken,
            partial_charges_nbo=partial_charges_nbo,
            electron_affinity=electron_affinity_eV,
            ionization_energy=ionization_energy_eV,
            spin_multiplicity=spin_multiplicity,
            partial_spins_nbo=partial_spins_nbo
        )



    def get_free_energy(self, temperature: float = ROOM_TEMP) -> Optional[float]:
        """
        Get the free energy at the give temperature.
        """
        if self.enthalpy is not None and self.entropy is not None:
            # TODO: fix these hard coded vals
            return (
                self.energy * 27.21139
                + 0.0433641 * self.enthalpy
                - temperature * self.entropy * 0.0000433641
            )
        else:
            return None

    def __repr__(self):

        output = [
            f"MoleculeEntry {self.entry_id} - {self.formula}",
            f"Total charge = {self.charge}",
        ]

        energies = [
            ("Energy", "Hartree", self.energy),
            ("Enthalpy", "kcal/mol", self.enthalpy),
            ("Entropy", "cal/mol.K", self.entropy),
            ("Free Energy (298.15 K)", "eV", self.get_free_energy()),
        ]
        for name, unit, value in energies:
            if value is None:
                output.append(f"{name} = {value} {unit}")
            else:
                output.append(f"{name} = {value:.4f} {unit}")

        if self.ind:
            output.append("index: {}".format(self.ind))

        return "\n".join(output)

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        if type(self) == type(other):
            return str(self) == str(other)
        else:
            return False
