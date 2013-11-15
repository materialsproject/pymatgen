#!/usr/bin/env python

"""
This module implements input and output processing from QChem.
"""
from pymatgen.core.structure import Molecule
from pymatgen.serializers.json_coders import MSONable

__author__ = "Xiaohui Qu"
__copyright__ = "Copyright 2013, The Electrolyte Genome Project"
__version__ = "0.1"
__maintainer__ = "Xiaohui Qu"
__email__ = "xhqu1981@gmail.com"
__date__ = "11/4/13"


class QcInput(MSONable):
    """
        An object representing a QChem input file.
    """

    optional_keywords_list = set(["basis", "ecp", "empirical_dispersion",
        "external_charges", "force_field_params", "intracule", "isotopes",
        "localized_diabatization", "multipole_field", "nbo", "occupied",
        "swap_occupied_virtual", "opt", "pcm", "pcm_solvent", "plots",
        "qm_atoms", "svp", "svpirf", "van_der_waals", "xc_functional", "cdft",
        "efp_fragments", "efp_params"])

    def __init__(self, molecule=None, charge=None, spin_multiplicity=None,
                 title=None, exchange="HF", correlation=None,
                 basis_set="6-31+G*", aux_basis_set=None, rem_params=None,
                 optional_params=None):
        """
        Args:
            molecule:
                The input molecule. If it is None of string "read", QChem will
                read geometry from checkpoint file. If it is a Molecule object,
                QcInput will convert it into Cartesian coordinates.
                Valid values: pymatgen Molecule object, "read", None
                Defaults to None.
            charge:
                Charge of the molecule. If None, charge on molecule is used.
                Defaults to None.
                Type: Integer
            spin_multiplicity:
                Spin multiplicity of molecule. Defaults to None,
                which means that the spin multiplicity is set to 1 if the
                molecule has no unpaired electrons and to 2 if there are
                unpaired electrons.
                Type: Integer
            title:
                Comments for the job. Defaults to None. Which means the
                $comment section will be discarded
                Type: str
            exchange:
                The exchange methods of the theory. Examples including:
                "B" (in pure BLYP), "PW91", "PBE", "TPSS".
                Defaults to "HF".
                This parameter can also be common names of hybrid
                functionals, such as B3LYP, TPSSh, XYGJOS. In such cases,
                the correlation parameter should be left as None.
                Type: str
            correlation:
                The correlation level of the theory. Example including:
                "MP2", "RI-MP2", "CCSD(T)", "LYP", "PBE", "TPSS"
                Defaults to None.
                Type: str
            basis_set:
                The basis set. If set to "gen", an additional parameter
                indicated by key "basis" is required in optional_params.
                Type: str
            aux_basis_set:
                Auxiliary basis set. For methods, like RI-MP2, XYG3, OXYJ-OS,
                auxiliary basis set is required.
                Type: str
            rem_params:
                The parameters supposed to write in the $rem section. Dict of
                key/value pairs.
                Type: dict
                Example: {"scf_algorithm": "diis_gdm", "scf_max_cycles": 100}
            optional_params:
                The parameter for keywords other than $rem section. Dict of
                key/value pairs.
                Type: dict
                Example: {"basis":
                          {"Li": "cc-PVTZ", "B": "aug-cc-PVTZ",
                           "F": "aug-cc-PVTZ"}
                          "ecp":
                          {"Cd": "srsc", "Br": "srlc"}
                         }

        """

        self.mol = molecule if molecule else "read"
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        if self.mol is not "read":
            if not isinstance(self.mol, Molecule):
                raise ValueError("The molecule must be a pymatgen Molecule "
                                 "object or read/None")
            self.charge = charge if charge is not None else self.mol.charge
            nelectrons = self.mol.charge + self.nelectrons - self.charge
            if spin_multiplicity is not None:
                self.spin_multiplicity = spin_multiplicity
                if (nelectrons + spin_multiplicity) % 2 != 1:
                    raise ValueError("Charge of {} and spin multiplicity of {} "
                                     "is not possible for this molecule"
                    .format(self.charge, spin_multiplicity))
            else:
                self.spin_multiplicity = 1 if nelectrons % 2 == 0 else 2
        if (self.charge is None) != (self.spin_multiplicity is None):
            raise ValueError("spin multiplicity must be set together")
        self.params = dict()
        if title is not None:
            self.params["comments"] = title
        if rem_params is not None:
            for k, v in rem_params:
                self.params["rem"][k.lower()] = v.lower()
        self.params["rem"]["basis"] = basis_set.lower()
        self.params["rem"]["exchange"] = exchange.lower()
        if correlation is not None:
            self.params["rem"]["correlation"] = correlation.lower()
        op_key = set([k.lower() for k in optional_params.keys()])
        if len(op_key - self.optional_keywords_list) > 0:
            invalid_keys = op_key - self.optional_keywords_list
            raise ValueError(','.join(['$' + k for k in invalid_keys]) +
                             'is not a valid section')
        self.params.update(optional_params)
        if aux_basis_set is None:
            if self._aux_basis_required():
                if self.params["rem"]["basis"].startswith("6-31+g"):
                    self.params["rem"]["aux_basis"] = "rimp2-aug-cc-pvdz"
                elif self.params["rem"]["basis"].startswith("6-311+g"):
                    self.params["rem"]["aux_basis"] = "rimp2-aug-cc-pvtz"
                else:
                    raise ValueError("Auxiliary basis set is missing")
        else:
            self.params["rem"]["aux_basis"] = aux_basis_set

    def _aux_basis_required(self):
        if self.params["rem"]["exchange"] in ['xygjos', 'xyg3', 'lxygjos']:
            return True
        if 'correlation' in self.params["rem"]:
            if self.params["rem"]["correlation"].startswith("ri"):
                return True