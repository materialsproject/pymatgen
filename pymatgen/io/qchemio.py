#!/usr/bin/env python
# coding=utf-8

"""
This module implements input and output processing from QChem.
"""
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

    optional_keywords_list = {"basis", "ecp", "empirical_dispersion",
                              "external_charges", "force_field_params",
                              "intracule", "isotopes", "aux_basis",
                              "localized_diabatization", "multipole_field",
                              "nbo", "occupied", "swap_occupied_virtual", "opt",
                              "pcm", "pcm_solvent", "plots", "qm_atoms", "svp",
                              "svpirf", "van_der_waals", "xc_functional",
                              "cdft", "efp_fragments", "efp_params"}

    def __init__(self, molecule=None, charge=None, spin_multiplicity=None,
                 jobtype='SP', title=None, exchange="HF", correlation=None,
                 basis_set="6-31+G*", aux_basis_set=None, ecp=None,
                 rem_params=None, optional_params=None):
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
            jobtype:
                The type the QChem job. "SP" for Single Point Energy,
                "opt" for geometry optimization, "freq" for
                vibrational frequency.
                Type: str
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
                The basis set.
                If it is a dict, each element can use different basis set.
                Type: str or dict
            aux_basis_set:
                Auxiliary basis set. For methods, like RI-MP2, XYG3, OXYJ-OS,
                auxiliary basis set is required.
                If it is a dict, each element can use different auxiliary
                basis set.
                Type: str or dict
                Type: str
            ecp:
                Effective core potential (ECP) to be used.
                If it is a dict, each element can use different ECP.
                Type: str or dict
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
        self.params["rem"]["exchange"] = exchange.lower()
        available_jobtypes = {"sp", "opt", "ts", "freq", "force", "rpath",
                              "nmr", "bsse", "eda", "pes_scan", "fsm", "aimd",
                              "pimc", "makeefp"}
        if jobtype.lower() not in available_jobtypes:
            raise ValueError("Job type " + jobtype + " is not supported yet")
        self.params["rem"]["jobtype"] = jobtype.lower()
        if correlation is not None:
            self.params["rem"]["correlation"] = correlation.lower()
        if rem_params is not None:
            for k, v in rem_params.iteritems():
                if isinstance(v, str):
                    self.params["rem"][k.lower()] = v.lower()
                elif isinstance(v, int):
                    self.params["rem"][k.lower()] = v
                else:
                    raise ValueError("The value in $rem can only be Integer "
                                     "or string")
        if optional_params:
            op_key = set([k.lower() for k in optional_params.keys()])
            if len(op_key - self.optional_keywords_list) > 0:
                invalid_keys = op_key - self.optional_keywords_list
                raise ValueError(','.join(['$' + k for k in invalid_keys]) +
                                 'is not a valid optional section')
            self.params.update(optional_params)

        self.set_basis_set(basis_set)

        if aux_basis_set is None:
            if self._aux_basis_required():
                if isinstance(self.params["rem"]["basis"], str):
                    if self.params["rem"]["basis"].startswith("6-31+g"):
                        self.set_auxiliary_basis_set("rimp2-aug-cc-pvdz")
                    elif self.params["rem"]["basis"].startswith("6-311+g"):
                        self.set_auxiliary_basis_set("rimp2-aug-cc-pvtz")
                if "aux_basis" not in self.params["rem"]:
                    raise ValueError("Auxiliary basis set is missing")
        else:
            self.set_auxiliary_basis_set(aux_basis_set)

        if ecp:
            self.set_ecp(ecp)

    def _aux_basis_required(self):
        if self.params["rem"]["exchange"] in ['xygjos', 'xyg3', 'lxygjos']:
            return True
        if 'correlation' in self.params["rem"]:
            if self.params["rem"]["correlation"].startswith("ri"):
                return True

    def set_basis_set(self, basis_set):
        if isinstance(basis_set, str):
            self.params["rem"]["basis"] = basis_set.lower()
        elif isinstance(basis_set, dict):
            self.params["rem"]["basis"] = "gen"
            bs = dict()
            for element, basis in basis_set.iteritems():
                bs[element.strip().capitalize()] = basis.lower()
            self.params["basis"] = bs
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                basis_elements = set(self.params["basis"].keys())
                if len(mol_elements - basis_elements) > 0:
                    raise ValueError("The basis set for elements " +
                                     ", ".join(
                                         list(mol_elements - basis_elements)) +
                                     " is missing")
                if len(basis_elements - mol_elements) > 0:
                    raise ValueError("Basis set error: the molecule "
                                     "doesn't contain element " +
                                     ", ".join(basis_elements - mol_elements))

    def set_auxiliary_basis_set(self, aux_basis_set):
        if isinstance(aux_basis_set, str):
            self.params["rem"]["aux_basis"] = aux_basis_set.lower()
        elif isinstance(aux_basis_set, dict):
            self.params["rem"]["aux_basis"] = "gen"
            bs = dict()
            for element, basis in aux_basis_set.iteritems():
                bs[element.strip().capitalize()] = basis.lower()
            self.params["aux_basis"] = bs
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                basis_elements = set(self.params["aux_basis"].keys())
                if len(mol_elements - basis_elements) > 0:
                    raise ValueError("The auxiliary basis set for "
                                     "elements " +
                                     ", ".join(
                                         list(mol_elements - basis_elements)) +
                                     " is missing")
                if len(basis_elements - mol_elements) > 0:
                    raise ValueError("Auxiliary asis set error: the "
                                     "molecule doesn't contain element " +
                                     ", ".join(basis_elements - mol_elements))

    def set_ecp(self, ecp):
        if isinstance(ecp, str):
            self.params["rem"]["ecp"] = ecp.lower()
        elif isinstance(ecp, dict):
            self.params["rem"]["ecp"] = "gen"
            potentials = dict()
            for element, p in ecp.iteritems():
                potentials[element.strip().capitalize()] = p.lower()
            self.params["ecp"] = potentials
            if self.mol:
                mol_elements = set([site.species_string for site
                                    in self.mol.sites])
                ecp_elements = set(self.params["ecp"].keys())
                if len(ecp_elements - mol_elements) > 0:
                    raise ValueError("ECP error: the molecule "
                                     "doesn't contain element " +
                                     ", ".join(ecp_elements - mol_elements))

    @property
    def molecule(self):
        return self.mol

    def set_memory(self, total=None, static=None):
        """
        Set the maxium allowed memory.

        Args:
            total: The total memory. Integer. Unit: MBytes. If set to None,
                this parameter will be neglected.
            static: The static memory. Integer. Unit MBytes. If set to None,
                this parameterwill be neglected.
        """
        if total:
            self.params["rem"]["mem_total"] = total
        if static:
            self.params["rem"]["mem_static"] = static

    def set_max_num_of_scratch_files(self, num=16):
        """
        In QChem, the size of a single scratch is limited 2GB. By default,
        the max number of scratich is 16, which is cooresponding to 32GB
        scratch space. If you want to use more scratch disk space, you need
        to increase the number of scratch files:

        Args:
            num: The max number of the scratch files. (Integer)
        """
        self.params["rem"]["max_sub_file_num"] = num

    def set_scf_algorithm_and_iterations(self, algorithm="diis",
                                         iterations=50):
        """
        Set algorithm used for converging SCF and max number of SCF iterations.

        Args:
            algorithm: The algorithm used for converging SCF. (str)
            iterations: The max number of SCF iterations. (Integer)
        """
        available_algorithms = {"diis", "dm", "diis_dm", "diis_gdm", "gdm",
                                "rca", "rca_diis", "roothaan"}
        if algorithm.lower() not in available_algorithms:
            raise ValueError("Algorithm " + algorithm +
                             " is not available in QChem")
        self.params["rem"]["scf_algorithm"] = algorithm.lower()
        self.params["rem"]["max_scf_cycles"] = iterations

    def set_scf_convergence_threshold(self, exponent=8):
        """
        SCF is considered converged when the wavefunction error is less than
        10**(-exponent).
        In QChem, the default values are:
            5	For single point energy calculations.
            7	For geometry optimizations and vibrational analysis.
            8	For SSG calculations

        Args:
            exponent: The exponent of the threshold. (Integer)
        """
        self.params["rem"]["scf_convergence"] = exponent

    def set_integral_threshold(self, thresh=12):
        """
        Cutoff for neglect of two electron integrals. 10−THRESH (THRESH ≤ 14).
        In QChem, the default values are:
            8	For single point energies.
            10	For optimizations and frequency calculations.
            14	For coupled-cluster calculations.

        Args:
            thresh. The exponent of the threshold. (Integer)
        """
        self.params["rem"]["thresh"] = thresh

    def set_dft_grid(self, radical_points=128, angular_points=302,
                     grid_type="Lebedev"):
        """
        Set the grid for DFT numerical integrations.

        Args:
            radical_points: Radical points. (Integer)
            angular_points: Angular points. (Integer)
            grid_type: The type of of the grid. There are two standard grids:
                SG-1 and SG-0. The other two supported grids are "Lebedev" and
                "Gauss-Legendre"
        """
        available_lebedev_angular_points = {6, 18, 26, 38, 50, 74, 86, 110, 146,
                                            170, 194, 230, 266, 302, 350, 434,
                                            590, 770, 974, 1202, 1454, 1730,
                                            2030, 2354, 2702, 3074, 3470, 3890,
                                            4334, 4802, 5294}
        if grid_type.lower() == "sg-0":
            self.params["rem"]["xc_grid"] = 0
        elif grid_type.lower() == "sg-1":
            self.params["rem"]["xc_grid"] = 1
        elif grid_type.lower() == "lebedev":
            if angular_points not in available_lebedev_angular_points:
                raise ValueError(str(angular_points) + " is not a valid "
                                 "Lebedev angular points number")
            self.params["rem"]["xc_grid"] = "{rp:06d}{ap:06d}".format(
                rp=radical_points, ap=angular_points)
        elif grid_type.lower() == "gauss-legendre":
            self.params["rem"]["xc_grid"] = "-{rp:06d}{ap:06d}".format(
                rp=radical_points, ap=angular_points)
        else:
            raise ValueError("Grid type " + grid_type + " is not supported "
                                                        "currently")

    def set_scf_initial_guess(self, guess="SAD"):
        """
        Set initial guess method to be used for SCF

        Args:
            guess: The initial guess method. (str)
        """
        availabel_guesses = {"core", "sad", "gwh", "read", "fragmo"}
        if guess.lower() not in availabel_guesses:
            raise ValueError("The guess method " + guess + " is not supported "
                                                           "yet")
        self.params["rem"]["scf_guess"] = guess.lower()

    def set_geom_max_iterations(self, iterations):
        """
        Set the max iterations of geometry optimization.

        Args:
            iterations: the maximum iterations of geometry optimization.
            (Integer)
        """
        self.params["rem"]["geom_opt_max_cycles"] = iterations

    def set_geom_opt_coords_type(self, coords_type="internal_switch"):
        """
        Set the coordinates system used in geometry optimization.
        "cartesian"       --- always cartesian coordinates.
        "internal"        --- always internal coordinates.
        "internal-switch" --- try internal coordinates first, if fails, switch
            to cartesian coordinates.
        "z-matrix"        --- always z-matrix coordinates.
        "z-matrix-switch" --- try z-matrix first, if fails, switch to
            cartesian coordinates.

        Args:
            coords_type: The type of the coordinates. (str)
        """
        coords_map = {"cartesian": 0, "internal": 1, "internal-switch": -1,
                      "z-matrix": 2, "z-matrix-switch": -2}
        if coords_type.lower() not in set(coords_map.keys()):
            raise ValueError("Coodinate system " + coords_type + " is not "
                             "supported yet")
        else:
            self.params["rem"]["geom_opt_coords"] = \
                coords_map[coords_type.lower()]

    def scale_geom_opt_threshold(self, gradient=0.1, displacement=0.1,
                                 energy=0.1):
        """
        Adjust the convergence criteria of geometry optimization.

        Args:
            gradient: the scale factor for gradient criteria. If less than
                1.0, you are tightening the threshold. The base value is
                300 × 10E−6
            displacement: the scale factor for atomic displacement. If less
                then 1.0, you are tightening the threshold. The base value is
                1200 × 10E−6
            energy: the scale factor for energy change between successive
                iterations. If less than 1.0, you are tightening the
                threshold. The base value is 100 × 10E−8.
        """
        if gradient < 1.0/(300-1) or displacement < 1.0/(1200-1) or \
                energy < 1.0/(100-1):
            raise ValueError("The geometry optimization convergence criteria "
                             "is too tight")
        self.params["rem"]["geom_opt_tol_gradient"] = int(gradient * 300)
        self.params["rem"]["geom_opt_tol_displacement"] = int(displacement *
                                                              1200)
        self.params["rem"]["geom_opt_tol_energy"] = int(energy * 100)

    def set_geom_opt_use_gdiis(self, subspace_size=None):
        """
        Use GDIIS algorithm in geometry optimization.

        Args:
            subspace_size: The size of the DIIS subsapce. None for default
                value. The default value is min(NDEG, NATOMS, 4) NDEG = number
                of moleculardegrees of freedom.
        """
        subspace_size = subspace_size if subspace_size else -1
        self.params["rem"]["geom_opt_max_diis"] = subspace_size

    def disable_symmetry(self):
        """
        Turn the symmetry off.
        """
        self.params["rem"]["sym_ignore"] = True
        self.params["rem"]["symmetry"] = False

    def __str__(self):
        sections = ["comments", "molecule", "rem"] + \
            sorted(list(self.optional_keywords_list))
        lines = []
        for sec in sections:
            if sec in self.params or sec == "molecule":
                foramt_sec = self.__getattribute__("_format_" + sec)
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
        if self.charge is not None:
            lines.append(" {charge:d}  {multi:d}".format(charge=self
                         .charge, multi=self.spin_multiplicity))
        if self.mol == "read":
            lines.append(" read")
        else:
            for site in self.mol.sites:
                lines.append(" {element:<4} {x:>17.8f} {y:>17.8f} "
                             "{z:>17.8f}".format(element=site.species_string,
                                                 x=site.x, y=site.y, z=site.z))
        return lines

    def _format_rem(self):
        rem_format_template = Template("  {name:>$name_width} = "
                                       "{value}")
        name_width = 0
        for name, value in self.params["rem"].iteritems():
            if len(name) > name_width:
                name_width = len(name)
        rem = rem_format_template.substitute(name_width=name_width)
        lines = []
        all_keys = set(self.params["rem"].keys())
        priority_keys = ["jobtype", "exchange", "basis"]
        additional_keys = all_keys - set(priority_keys)
        ordered_keys = priority_keys + sorted(list(additional_keys))
        for name in ordered_keys:
            value = self.params["rem"][name]
            lines.append(rem.format(name=name, value=value))
        return lines

    def _format_basis(self):
        lines = []
        for element in sorted(self.params["basis"].keys()):
            basis = self.params["basis"][element]
            lines.append(" " + element)
            lines.append(" " + basis)
            lines.append(" ****")
        return lines

    def _format_aux_basis(self):
        lines = []
        for element in sorted(self.params["aux_basis"].keys()):
            basis = self.params["aux_basis"][element]
            lines.append(" " + element)
            lines.append(" " + basis)
            lines.append(" ****")
        return lines

    def _format_ecp(self):
        lines = []
        for element in sorted(self.params["ecp"].keys()):
            ecp = self.params["ecp"][element]
            lines.append(" " + element)
            lines.append(" " + ecp)
            lines.append(" ****")
        return lines

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "molecule": "read" if self.mol == "read" else self.mol.to_dict,
                "charge": self.charge,
                "spin_multiplicity": self.spin_multiplicity,
                "params": self.params}

    @classmethod
    def from_dict(cls, d):
        mol = "read" if d["molecule"] == "read" \
            else Molecule.from_dict(d["molecule"])
        jobtype = d["params"]["rem"]["jobtype"]
        title = d["params"].get("comments", None)
        exchange = d["params"]["rem"]["exchange"]
        correlation = d["params"]["rem"].get("correlation", None)
        basis_set = d["params"]["rem"]["basis"]
        aux_basis_set = d["params"]["rem"].get("aux_basis", None)
        ecp = d["params"]["rem"].get("ecp", None)
        optional_params = None
        op_keys = set(d["params"].keys()) - {"comments", "rem"}
        if len(op_keys) > 0:
            optional_params = dict()
            for k in op_keys:
                optional_params[k] = d["params"][k]
        return QcInput(molecule=mol, charge=d["charge"],
                       spin_multiplicity=d["spin_multiplicity"],
                       jobtype=jobtype, title=title,
                       exchange=exchange, correlation=correlation,
                       basis_set=basis_set, aux_basis_set=aux_basis_set,
                       ecp=ecp, rem_params=d["params"]["rem"],
                       optional_params=optional_params)

    def write_file(self, filename):
        with zopen(filename, "w") as f:
            f.write(self.__str__())
