#!/usr/bin/env python

"""
This module implements input and output processing from QChem.
"""
import os
from string import Template
from pymatgen import zopen
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
                                  "external_charges", "force_field_params",
                                  "intracule", "isotopes",
                                  "localized_diabatization", "multipole_field",
                                  "nbo", "occupied",
                                  "swap_occupied_virtual", "opt", "pcm",
                                  "pcm_solvent", "plots",
                                  "qm_atoms", "svp", "svpirf", "van_der_waals",
                                  "xc_functional", "cdft",
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
        if isinstance(self.mol, str):
            self.mol = self.mol.lower()
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        if self.mol is not "read":
            if not isinstance(self.mol, Molecule):
                raise ValueError("The molecule must be a pymatgen Molecule "
                                 "object or read/None")
            self.charge = charge if charge is not None else self.mol.charge
            nelectrons = self.mol.charge + self.mol.nelectrons - self.charge
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
        if "rem" not in self.params:
            self.params["rem"] = dict()
        self.params["rem"]["basis"] = basis_set.lower()
        self.params["rem"]["exchange"] = exchange.lower()
        if correlation is not None:
            self.params["rem"]["correlation"] = correlation.lower()
        if rem_params is not None:
            for k, v in rem_params:
                self.params["rem"][k.lower()] = v.lower()
        op_key = set([k.lower() for k in optional_params.keys()])
        if len(op_key - self.optional_keywords_list) > 0:
            invalid_keys = op_key - self.optional_keywords_list
            raise ValueError(','.join(['$' + k for k in invalid_keys]) +
                             'is not a valid optional section')
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

    @property
    def molecule(self):
        return self.mol


    def __str__(self):
        sections = ["comments", "molecule", "rem"] + \
                   sorted(list(self.optional_keywords_list))
        lines = []
        for sec in sections:
            if sec in self.params:
                foramt_sec = self.__getattribute__("_format" + sec)
                lines.append("$" + sec)
                lines.extend(foramt_sec())
                lines.append("$end")
                lines.append('\n')
        return '\n'.join(lines)


    def _format_comments(self):
        lines = [' ' + self.params["comments"].strip()]
        return lines


    def _format_molecule(self):
        lines = []
        if self.charge:
            lines.append(" {charge:%d}  {multi:%d}".format(charge=self
            .charge, multi=self.spin_multiplicity))
        if self.mol == "read":
            lines.append(" read")
        else:
            for site in self.mol.sites:
                lines.append(" {element:<%4s} {x:>%12.8f} {y:>%12.8f} "
                             "{z:>%12.8f}".format(element=site.species_string,
                                                  x=site.x, y=site.y, z=site.z))
        return lines


    def _format_rem(self):
        rem_format_template = Template(" {name:>%$name_width} = "
                                       "{vaule:<%$value_with}")
        name_width = 0
        value_width = 0
        for name, value in self.params["rem"].iteritems():
            if len(name) > name_width:
                name_width = len(name)
            if len(value) > value_width:
                value_width = len(value)
        rem = rem_format_template.substitute(name_width=name_width,
                                             value_width=value_width)
        lines = []
        for name in sorted(self.params["rem"].keys()):
            value = self.params["rem"][name]
            lines.append(rem.format(name=name, value=value))
        return lines

    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "molecule": "read" if self.mol == "read" else self.mol.to_dict,
                "charge": self.charge,
                "spin_multiplicity": self.spin_multiplicity,
                "parameters": self.params}

    def from_dict(cls, d):
        mol = "read" if d["molecule"] == "read" \
                     else Molecule.from_dict(d["molecule"])
        title = d["params"].get("comments", None)
        exchange = d["params"]["rem"]["exchange"]
        correlation = d["params"]["rem"].get("correlation", None)
        basis_set = d["params"]["rem"]["basis"]
        aux_basis_set = d["params"]["rem"].get("aux_basis", None)
        optional_params = None
        op_keys = set(d["params"].keys()) - set(["comments", "rem"])
        if len(op_keys) > 0:
            optional_params = dict()
            for k in op_keys:
                optional_params[k] = d["params"][k]
        return QcInput(molecule=mol, charge=d["charge"],
                       spin_multiplicity=d["spin_multiplicity"], title=title,
                       exchange=exchange, correlation=correlation,
                       basis_set=basis_set, aux_basis_set=aux_basis_set,
                       rem_params=d["params"]["rem"], optional_params=optional_params)

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())

