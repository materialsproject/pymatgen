# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals, \
    absolute_import

from collections import OrderedDict
import itertools
import warnings
import re
from ast import literal_eval

import numpy as np
from monty.json import MSONable
from six import string_types
from ruamel.yaml import YAML

from pymatgen.util.io_utils import clean_lines
from pymatgen.core.structure import SiteCollection
from pymatgen import Molecule, Element

"""
This module implements a core class LammpsData for generating/parsing 
LAMMPS data file, and other bridging classes to build LammpsData from 
molecules. 

Only point particle styles are supported for now (atom_style in angle,
atomic, bond, charge, full and molecular only). See the pages below for 
more info.

    http://lammps.sandia.gov/doc/atom_style.html
    http://lammps.sandia.gov/doc/read_data.html

"""

__author__ = "Kiran Mathew, Zhi Deng"
__email__ = "kmathew@lbl.gov, z4deng@eng.ucsd.edu"
__credits__ = "Brandon Wood"


SECTION_KEYWORDS = {"atom": ["Atoms", "Velocities", "Masses",
                             "Ellipsoids", "Lines", "Triangles", "Bodies"],
                    "molecule": ["Bonds", "Angles", "Dihedrals", "Impropers"],
                    "ff": ["Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs",
                           "Angle Coeffs", "Dihedral Coeffs",
                           "Improper Coeffs"],
                    "class2": ["BondBond Coeffs", "BondAngle Coeffs",
                               "MiddleBondTorsion Coeffs",
                               "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
                               "AngleAngleTorsion Coeffs",
                               "BondBond13 Coeffs", "AngleAngle Coeffs"]}

ATOMS_LINE_FORMAT = {"angle": ["molecule-ID", "type", "x", "y", "z"],
                     "atomic": ["type", "x", "y", "z"],
                     "bond": ["molecule-ID", "type", "x", "y", "z"],
                     "charge": ["type", "q", "x", "y", "z"],
                     "full": ["molecule-ID", "type", "q", "x", "y", "z"],
                     "molecular": ["molecule-ID", "type", "x", "y", "z"]}

ATOMS_FLOATS = ["q", "x", "y", "z"]


class LammpsData(MSONable):
    """
    Object for representing the data in a LAMMPS data file.
    """

    def __init__(self, masses, atoms, box_bounds, box_tilt=None,
                 velocities=None, force_field=None, topology=None,
                 atom_style="full"):
        """
        This constructor is designed to work with parsed data from a
        file. Not recommended to use directly.

        Args:
            masses ([dict]): Data for Masses section.
                [{"id": 1, "mass": 1.008}, ...]
            atoms ([dict]): Data for Atoms section. Keys in dicts
                varies with atom_style.
                [{"id": 1, "type": 1,
                  "x": 0.0, "y": 0.0, "z": 0.0, ...}, ...]
            box_bounds: A (3, 2) array/list of floats setting the
                boundaries of simulation box.
            box_tilt: A (3,) array/list of floats setting the tilt of
                simulation box. Default to None, i.e., use an
                orthogonal box.
            velocities ([dict]): Data for Velocities section. Default
                to None. If not None, its length and ids should be
                consistent with atoms.
                [{"id": 1, "velocity": [0.0, 0.0, 0.0]}, ...]
            force_field (dict): Data for force field sections.
                Default to None. All six valid keys listed below are
                optional.
                {
                    "Pair Coeffs":
                        [{"id": 1, "coeffs": [coeff]}, ...],
                    "Pair IJ Coeffs":
                        [{"id1": 1, "id2": 1, "coeffs": [coeff]}, ...],
                    "Bond Coeffs":
                        [{"id": 1, "coeffs": [coeff]}, ...],
                    "Angle Coeffs":
                        [{"id": 1, "coeffs": [coeff]}, ...],
                    "Dihedral Coeffs":
                        [{"id": 1, "coeffs": [coeff]}, ...],
                    "Improper Coeffs":
                        [{"id": 1, "coeffs": [coeff]}, ...],
                }
            topology (dict): Data for topology sections. Default to
                None. All four valid keys listed below are optional.
                {
                    "Bonds":
                        [{"id": 1, "bond": [1, 2]}, ...],
                    "Angles":
                        [{"id": 1, "angle": [2, 1, 3]}, ...],
                    "Dihedrals":
                        [{"id": 1, "dihedral": [3, 1, 2, 4]}, ...],
                    "Impropes":
                        [{"id": 1, "improper": [3, 1, 2, 4]}, ...],
                }
            atom_style (str): Output atom_style. Default to "full".

        """
        bounds_arr = np.array(box_bounds)
        bounds_shape = bounds_arr.shape
        assert bounds_shape == (3, 2), \
            "Expecting a (3, 2) array for box_bounds," \
            " got {}".format(bounds_shape)
        box_bounds = bounds_arr.tolist()

        if box_tilt is not None:
            tilt_arr = np.array(box_tilt)
            tilt_shape = tilt_arr.shape
            assert tilt_shape == (3,),\
                "Expecting a (3,) array for box_tilt," \
                " got {}".format(tilt_shape)
            box_tilt = tilt_arr.tolist()

        if velocities:
            assert len(velocities) == len(atoms),\
                "Inconsistency found between atoms and velocities"

        if force_field:
            force_field = {k: force_field[k] for k in SECTION_KEYWORDS["ff"]
                           if k in force_field}
        if topology:
            topology = {k: topology[k] for k in SECTION_KEYWORDS["molecule"]
                        if k in topology}
        self.masses = masses
        self.atoms = atoms
        self.box_bounds = box_bounds
        self.box_tilt = box_tilt
        self.velocities = velocities
        self.force_field = force_field
        self.topology = topology
        self.atom_style = atom_style

    def __str__(self):
        return self.get_string()

    def get_string(self, significant_figures=6):
        """
        Returns the string representation of LammpsData, equivalent to
        the string that is to be written into a file.

        Args:
            significant_figures (int): No. of significant figures of
                (changeable) quantities to output, default ot 6.
                Quantities include box bounds and tilt, coordinates,
                charges and velocities. While other stationary
                quantities, like masses and force field coefficients,
                are displayed as-is.

        Returns:
            String representation

        """
        template = """Generated by pymatgen.io.lammps.data.LammpsData

{stats}

{box}

{masses}

{ff_sections}

{atoms}

{velocities}

{topo_sections}
"""

        contents = {}

        float_ph = "{:.%df}" % significant_figures \
            if significant_figures else "{}"

        def _pretty_section(title, str_mat):
            lens = [max(map(len, col)) for col in zip(*str_mat)]
            fmt = "  ".join("{:>%d}" % x for x in lens)
            rows = [title, ""] + [fmt.format(*row) for row in str_mat]
            sec = "\n".join(rows)
            return sec

        counts = OrderedDict([("atoms", len(self.atoms))])
        types = OrderedDict([("atom", len(self.masses))])

        # box
        box_lines = []
        for bound, d in zip(self.box_bounds, "xyz"):
            fillers = bound + [d] * 2
            bound_format = " ".join([float_ph] * 2 + ["{}lo {}hi"])
            box_lines.append(bound_format.format(*fillers))
        if self.box_tilt:
            tilt_format = " ".join([float_ph] * 3 + ["xy xz yz"])
            box_lines.append(tilt_format.format(*self.box_tilt))
        contents["box"] = "\n".join(box_lines)

        # masses
        masses_mat = [["%d" % m["id"], "{:.4f}".format(m["mass"])]
                      for m in self.masses]
        contents["masses"] = _pretty_section("Masses", masses_mat)

        # ff_sections
        contents["ff_sections"] = ""
        if self.force_field:
            ff_parts = []
            for kw in self.force_field.keys():
                if kw == "PairIJ Coeffs":
                    ff_mat = [[str(i) for i in
                               [d["id1"], d["id2"]] + d["coeffs"]]
                              for d in self.force_field[kw]]
                else:
                    ff_mat = [[str(i) for i in [d["id"]] + d["coeffs"]]
                              for d in self.force_field[kw]]
                if not kw.startswith("Pair"):
                    types[kw.lower()[:-7]] = len(ff_mat)
                ff_parts.append(_pretty_section(kw, ff_mat))
            contents["ff_sections"] = "\n\n".join(ff_parts)

        # atoms
        atom_format = ["id"] + ATOMS_LINE_FORMAT[self.atom_style]
        if "nx" in self.atoms[0].keys():
            atom_format.extend(["nx", "ny", "nz"])
        map_str = lambda t: float_ph if t in ATOMS_FLOATS else "{}"
        atoms_mat = []
        for a in self.atoms:
            atoms_mat.append([map_str(k).format(a[k]) for k in atom_format])
        contents["atoms"] = _pretty_section("Atoms", atoms_mat)

        # velocities
        contents["velocities"] = ""
        if self.velocities:
            velocities_mat = []
            for v in self.velocities:
                vs = [float_ph.format(i) for i in v["velocity"]]
                velocities_mat.append(["%d" % v["id"]] + vs)
            contents["velocities"] = _pretty_section("Velocities",
                                                     velocities_mat)

        # topo_sections
        contents["topo_sections"] = ""
        if self.topology:
            topo_parts = []
            topo_keys = [k for k in SECTION_KEYWORDS["molecule"]
                         if k in self.topology]
            for kw in topo_keys:
                skw = kw.lower()[:-1]
                topo_mat = [["%d" % v for v in [d["id"], d["type"]] + d[skw]]
                            for d in self.topology[kw]]
                counts[kw.lower()] = len(topo_mat)
                topo_parts.append(_pretty_section(kw, topo_mat))
            contents["topo_sections"] = "\n\n".join(topo_parts)

        # stats
        all_stats = list(counts.values()) + list(types.values())
        line_fmt = "{:>%d} {}" % len(str(max(all_stats)))
        count_lines = [line_fmt.format(v, k) for k, v in counts.items()]
        type_lines = [line_fmt.format(v, k + " types")
                      for k, v in types.items()]
        contents["stats"] = "\n".join(count_lines + [""] + type_lines)

        return template.format(**contents)

    def write_file(self, filename, significant_figures=6):
        """
        Writes LammpsData to file.

        Args:
            filename (str): Filename to write to.
            significant_figures (int): No. of significant figures of
                (changeable) quantities to output, default ot 6.
                Quantities include box bounds and tilt, coordinates,
                charges and velocities. While other stationary
                quantities, like masses and force field coefficients,
                are displayed as-is.

        """
        with open(filename, "w") as f:
            f.write(self.get_string(significant_figures=significant_figures))

    @classmethod
    def from_file(cls, filename, atom_style="full", sort_id=False):
        """
        Constructor from parsing a file.

        Args:
            filename (str): Filename to read.
            atom_style (str): Associated atom_style. Default to "full".
            sort_id (bool): Whether sort each section by id. Default to
                True.

        """
        with open(filename) as f:
            lines = f.readlines()
        clines = list(clean_lines(lines))
        section_marks = [i for i, l in enumerate(clines) if l
                         in itertools.chain(*SECTION_KEYWORDS.values())]
        parts = np.split(clines, section_marks)

        # First, parse header
        float_group = r'([0-9eE.+-]+)'
        header_pattern = {}
        header_pattern["counts"] = r'^\s*(\d+)\s+([a-zA-Z]+)$'
        header_pattern["types"] = r'^\s*(\d+)\s+([a-zA-Z]+)\s+types$'
        header_pattern["bounds"] = r'^\s*{}$'.format(r'\s+'.join(
            [float_group] * 2 + [r"([xyz])lo \3hi"]))
        header_pattern["tilt"] = r'^\s*{}$'.format(r'\s+'.join(
            [float_group] * 3 + ["xy xz yz"]))

        header = {"counts": {}, "types": {}}
        bounds = {}
        for l in parts[0]:
            match = None
            for k, v in header_pattern.items():
                match = re.match(v, l)
                if match:
                    break
                else:
                    continue
            if match and k in ["counts", "types"]:
                header[k][match.group(2)] = int(match.group(1))
            elif match and k == "bounds":
                g = match.groups()
                bounds[g[2]] = [float(i) for i in g[:2]]
            elif match and k == "tilt":
                header["tilt"] = [float(i) for i in match.groups()]
        header["bounds"] = [bounds.get(i, [-0.5, 0.5]) for i in "xyz"]

        # Then, parse each section
        topo_sections = SECTION_KEYWORDS["molecule"]

        def parse_section(single_section_lines):
            kw = single_section_lines[0]

            if kw in SECTION_KEYWORDS["ff"] and kw != "PairIJ Coeffs":
                parse_line = lambda l: {"coeffs": [literal_eval(x)
                                                   for x in l[1:]]}
            elif kw == "PairIJ Coeffs":
                parse_line = lambda l: {"id1": int(l[0]), "id2": int(l[1]),
                                        "coeffs": [literal_eval(x)
                                                   for x in l[2:]]}
            elif kw in topo_sections:
                n = {"Bonds": 2, "Angles": 3, "Dihedrals": 4, "Impropers": 4}
                parse_line = lambda l: {"type": int(l[1]), kw[:-1].lower():
                    [int(x) for x in l[2:n[kw] + 2]]}
            elif kw == "Atoms":
                keys = ATOMS_LINE_FORMAT[atom_style][:]
                sample_l = single_section_lines[1].split()
                if len(sample_l) == len(keys) + 1:
                    pass
                elif len(sample_l) == len(keys) + 4:
                    keys += ["nx", "ny", "nz"]
                else:
                    warnings.warn("Atoms section format might be imcompatible"
                                  " with atom_style %s." % atom_style)
                float_keys = [k for k in keys if k in ATOMS_FLOATS]
                parse_line = lambda l: {k: float(v) if k in float_keys
                else int(v) for (k, v) in zip(keys, l[1:len(keys) + 1])}
            elif kw == "Velocities":
                parse_line = lambda l: {"velocity": [float(x)
                                                     for x in l[1:4]]}
            elif kw == "Masses":
                parse_line = lambda l: {"mass": float(l[1])}
            else:
                warnings.warn("%s section parser has not been implemented. "
                              "Skipping..." % kw)
                return kw, []

            section = []
            splitted_lines = [l.split() for l in single_section_lines[1:]]
            if sort_id and kw != "PairIJ Coeffs":
                splitted_lines = sorted(splitted_lines,
                                        key=lambda l: int(l[0]))
            for l in splitted_lines:
                line_data = parse_line(l)
                if kw != "PairIJ Coeffs":
                    line_data["id"] = int(l[0])
                section.append(line_data)
            return kw, section

        err_msg = "Bad LAMMPS data format where "
        body = {}
        seen_atoms = False
        for part in parts[1:]:
            name, section = parse_section(part)
            if name == "Atoms":
                seen_atoms = True
            if name in ["Velocities"] + topo_sections and not seen_atoms:
                raise RuntimeError(err_msg + "%s section appears before"
                                             " Atoms section" % name)
            body.update({name: section})

        err_msg += "Nos. of {} do not match between header and {} section"
        assert len(body["Masses"]) == header["types"]["atom"], \
            err_msg.format("atom types", "Masses")
        atom_sections = ["Atoms", "Velocities"] \
            if body.get("Velocities") else ["Atoms"]
        for s in atom_sections:
            assert len(body[s]) == header["counts"]["atoms"], \
                err_msg.format("atoms", s)
        for s in topo_sections:
            if header["counts"].get(s.lower(), 0) > 0:
                assert len(body[s]) == header["counts"][s.lower()], \
                    err_msg.format(s.lower(), s)

        items = {k.lower(): body[k] for k in ["Masses", "Atoms"]}
        items["box_bounds"] = header["bounds"]
        items["box_tilt"] = header.get("tilt")
        items["velocities"] = body.get("Velocities")
        ff_kws = [k for k in body.keys() if k in SECTION_KEYWORDS["ff"]]
        items["force_field"] = {k: body[k] for k in ff_kws} if ff_kws \
            else None
        topo_kws = [k for k in body.keys()
                    if k in SECTION_KEYWORDS["molecule"]]
        items["topology"] = {k: body[k] for k in topo_kws} \
            if topo_kws else None
        items["atom_style"] = atom_style
        return cls(**items)

    @classmethod
    def from_ff_and_topologies(cls, ff, topologies, box_bounds, box_tilt=None,
                               atom_style="full"):
        """
        Constructor building LammpsData from a ForceField object and a
        list of Topology objects.

        Args:
            ff (ForceField): ForceField object with data for Masses and
                force field sections.
            topologies ([Topology]): List of Topology objects with data
                for Atoms, Velocities and topology sections.
            box_bounds: A (3, 2) array/list of floats setting the
                boundaries of simulation box.
            box_tilt: A (3,) array/list of floats setting the tilt of
                simulation box. Default to None, i.e., use an
                orthogonal box.
            atom_style (str): Output atom_style. Default to "full".

        """
        atom_types = set(itertools.chain(*[t.types for t in topologies]))
        assert atom_types.issubset(ff.atom_map.keys()),\
            "Unknown atom type found in topologies"

        items = {"box_bounds": box_bounds, "box_tilt": box_tilt,
                 "atom_style": atom_style}
        items["masses"] = ff.masses
        lookup = {"Atoms": ff.atom_map}

        pair_coeffs = ff.get_pair_coeffs()
        mol_coeffs = getattr(ff, "mol_coeffs")
        force_field = {} if any((pair_coeffs, mol_coeffs)) else None
        if pair_coeffs:
            force_field.update(pair_coeffs)
        if mol_coeffs:
            for kw in mol_coeffs.keys():
                coeffs, mapper = ff.get_coeffs_and_mapper(kw)
                force_field.update(coeffs)
                lookup[kw[:-7] + "s"] = mapper
        items["force_field"] = force_field

        atoms = []
        velocities = [] if topologies[0].velocities else None
        topology = {k: [] for k in SECTION_KEYWORDS["molecule"]}
        stack = {k: 0 for k in ["Atoms"] + SECTION_KEYWORDS["molecule"]}
        atom_format = ATOMS_LINE_FORMAT[atom_style]

        for mid, topo in enumerate(topologies):
            map_inds = lambda inds: tuple([topo.types[i] for i in inds])

            topo_atoms = []
            for aid, (s, t) in enumerate(zip(topo.sites, topo.types)):
                d_atom = {"id": aid + 1 + stack["Atoms"],
                          "type": lookup["Atoms"][t]}
                d_atom.update({k: getattr(s, k) for k in "xyz"})
                if "molecule-ID" in atom_format:
                    d_atom["molecule-ID"] = mid + 1
                topo_atoms.append(d_atom)
            if "q" in atom_format:
                charges = [0.0] * len(topo.sites) if not topo.charges \
                    else topo.charges
                for d_atom, q in zip(topo_atoms, charges):
                    d_atom["q"] = q
            atoms.extend(topo_atoms)

            if isinstance(velocities, list):
                velocities.extend({"id": aid + 1 + stack["Atoms"],
                                   "velocity": v}
                                  for aid, v in enumerate(topo.velocities))

            if topo.topologies:
                for kw in topo.topologies.keys():
                    topo_lookup = lookup[kw]
                    unfiltered_indices = np.array(topo.topologies[kw])
                    topo_topos = []
                    tid = stack[kw]
                    for inds in unfiltered_indices:
                        topo_type = topo_lookup.get(map_inds(inds))
                        if topo_type:
                            topo_inds = list(inds + stack["Atoms"] + 1)
                            topo_topos.append({"id": tid + 1,
                                               "type": topo_type,
                                               kw.lower()[:-1]: topo_inds})
                            tid += 1
                    topology[kw].extend(topo_topos)
                    stack[kw] = tid

            stack["Atoms"] += len(topo_atoms)

        topology = {k: v for k, v in topology.items() if len(v) > 0}
        topology = None if len(topology) == 0 else topology
        items.update({"atoms": atoms, "velocities": velocities,
                      "topology": topology})
        return cls(**items)


class Topology(MSONable):
    """
    Class carrying most data in Atoms, Velocities and molecular
    topology sections for ONE single SiteCollection or its subclasses
    (Molecule/Structure), or a plain list of Sites.

    """

    def __init__(self, sites, atom_type=None, charges=None, velocities=None,
                 topologies=None):
        """

        Args:
            sites ([Site] or SiteCollection): A group of sites in a
                list or as a Molecule/Structure.
            atom_type (str): Site property key for labeling atoms of
                different types. Default to None, i.e., use
                site.species_string.
            charges ([q, ...]): Charge of each site in a (n,)
                array/list, where n is the No. of sites. Default to
                None, i.e., search site property for charges.
            velocities ([[vx, vy, vz], ...]): Velocity of each site in
                a (n, 3) array/list, where n is the No. of sites.
                Default to None, i.e., search site property for
                velocities.
            topologies (dict): Bonds, angles, dihedrals and improper
                dihedrals defined by site indices. Default to None,
                i.e., no additional topology. All four valid keys
                listed below are optional.
                {
                    "Bonds": [[i, j], ...],
                    "Angles": [[i, j, k], ...],
                    "Dihedrals": [[i, j, k, l], ...],
                    "Impropers": [[i, j, k, l], ...]
                }

        """
        if not isinstance(sites, SiteCollection):
            sites = Molecule.from_sites(sites)

        if atom_type:
            types = sites.site_properties.get(atom_type)
        else:
            types = [site.species_string for site in sites]
        # search for site property if not override
        if charges is None:
            charges = sites.site_properties.get("charge")
        if velocities is None:
            velocities = sites.site_properties.get("velocities")
        # validate shape
        if charges is not None:
            charge_arr = np.array(charges)
            assert charge_arr.shape == (len(sites),),\
                "Wrong format for charges"
            charges = charge_arr.tolist()
        if velocities is not None:
            velocities_arr = np.array(velocities)
            assert velocities_arr.shape == (len(sites), 3), \
                "Wrong format for velocities"
            velocities = velocities_arr.tolist()

        if topologies:
            topologies = {k: topologies[k] for k in
                          SECTION_KEYWORDS["molecule"] if k in topologies}

        self.sites = sites
        self.atom_type = atom_type
        self.types = types
        self.charges = charges
        self.velocities = velocities
        self.topologies = topologies

    @classmethod
    def from_bonding(cls, molecule, bond=True, angle=True, dihedral=True,
                     atom_type=None, charges=None, velocities=None, tol=0.1):
        """
        Another constructor that creates an instance from a molecule.
        Covalent bonds and other bond-based topologies (angles and
        dihedrals) can be automatically determined. Cannot be used for
        non bond-based topologies, e.g., improper dihedrals.

        Args:
            molecule (Molecule): Input molecule.
            bond (bool): Whether find bonds. If set to False, angle and
                dihedral searching will be skipped. Default to True.
            angle (bool): Whether find angles. Default to True.
            dihedral (bool): Whether find dihedrals. Default to True.
            atom_type (str): Site property key for labeling atoms of
                different types. Default to None, i.e., use
                site.species_string.
            charges ([q, ...]): Charge of each site in a (n,)
                array/list, where n is the No. of sites. Default to
                None, i.e., search site property for charges.
            velocities ([[vx, vy, vz], ...]): Velocity of each site in
                a (n, 3) array/list, where n is the No. of sites.
                Default to None, i.e., search site property for
                velocities.
            tol (float): Bond distance tolerance. Default to 0.1.
                Not recommended to alter.

        """
        real_bonds = molecule.get_covalent_bonds(tol=tol)
        bond_list = [list(map(molecule.index, [b.site1, b.site2]))
                     for b in real_bonds]
        if not all((bond, bond_list)):
            return cls(sites=molecule, atom_type=atom_type, charges=charges,
                       velocities=velocities)
        else:
            angle_list, dihedral_list = [], []
            dests, freq = np.unique(bond_list, return_counts=True)
            hubs = dests[np.where(freq > 1)]
            bond_arr = np.array(bond_list)
            if len(hubs) > 0:
                hub_spokes = {}
                for hub in hubs:
                    ix = np.any(np.isin(bond_arr, hub), axis=1)
                    bonds = list(np.unique(bond_arr[ix]))
                    bonds.remove(hub)
                    hub_spokes[hub] = bonds
            dihedral = False if len(bond_list) < 3 or len(hubs) < 2 \
                else dihedral
            angle = False if len(bond_list) < 2 or len(hubs) < 1 else angle

            if angle:
                for k, v in hub_spokes.items():
                    angle_list.extend([[i, k, j] for i, j in
                                       itertools.combinations(v, 2)])
            if dihedral:
                hub_cons = bond_arr[np.all(np.isin(bond_arr, hubs), axis=1)]
                for i, j in hub_cons:
                    ks = [k for k in hub_spokes[i] if k != j]
                    ls = [l for l in hub_spokes[j] if l != i]
                    dihedral_list.extend([[k, i, j, l] for k,l in
                                          itertools.product(ks, ls)
                                          if k != l])

            topologies = {k: v for k, v
                          in zip(SECTION_KEYWORDS["molecule"][:3],
                                 [bond_list, angle_list, dihedral_list])
                          if len(v) > 0}
            topologies = None if len(topologies) == 0 else topologies
            return cls(sites=molecule, atom_type=atom_type, charges=charges,
                       velocities=velocities, topologies=topologies)


class ForceField(MSONable):
    """
    Class carrying most data in Masses and force field sections.
    """

    def __init__(self, mass_info, pair_coeffs=None, mol_coeffs=None):
        """

        Args:
            mass_into (list): List of atomic mass info. Each item looks
                like a key, value pair in an OrderedDict. Elements,
                strings (symbols) and floats are all acceptable for the
                values, with the first two converted to the atomic mass
                of an element.
                [("C": 12.01), ("H": Element("H")), ("O": "O"), ...]
            pair_coeffs [coeffs]: List of pair or pairij coefficients,
                of which the sequence must be sorted according to the
                species in mass_dict. Pair or PairIJ determined by the
                length of list.
            mol_coeffs (dict): Dict with force field coefficients for
                molecular topologies. Default to None, i.e., no
                additional coefficients. All four valid keys listed
                below are optional.
                {
                    "Bond Coeffs":
                        [{"coeffs": [coeffs],
                          "types": [("C", "C"), ...]}, ...],
                    "Angle Coeffs":
                        [{"coeffs": [coeffs],
                          "types": [("H", "C", "H"), ...]}, ...],
                    "Dihedral Coeffs":
                        [{"coeffs": [coeffs],
                          "types": [("H", "C", "C", "H"), ...]}, ...],
                    "Improper Coeffs":
                        [{"coeffs": [coeffs],
                          "types": [("H", "C", "C", "H"), ...]}, ...],
                }
        """
        map_mass = lambda v: v.atomic_mass.real if isinstance(v, Element) \
            else Element(v).atomic_mass.real if isinstance(v, string_types) \
            else v
        mass_info = [(k, map_mass(v)) for k, v in mass_info]

        if pair_coeffs:
            # validate No. of pair coeffs
            npc = len(pair_coeffs)
            nm = len(mass_info)
            ncomb = nm * (nm + 1) / 2
            if npc == nm:
                self.pair_type = "pair"
            elif npc == ncomb:
                self.pair_type = "pairij"
            else:
                raise ValueError("Expecting {} Pair Coeffs or "
                                 "{} PairIJ Coeffs for {} atom types,"
                                 " got {}".format(nm, ncomb, nm, npc))
        else:
            self.pair_type = None

        if mol_coeffs:
            mol_coeffs = {k: mol_coeffs[k] for k in SECTION_KEYWORDS["ff"][2:]
                          if k in mol_coeffs}
            complete_types = lambda l: l + [t[::-1] for t in l
                                            if t[::-1] not in l]
            for k, v in mol_coeffs.items():
                for d in v:
                    d["types"] = complete_types(d["types"])
                # No duplicated items under different types allowed
                distinct_types = [set(d["types"]) for d in v]
                if len(distinct_types) > 1:
                    assert set.intersection(*distinct_types) == set(),\
                        "Duplicated items found " \
                        "under different coefficients in %s" % k
                # No undefined atom types allowed
                atoms = set(np.ravel(list(itertools.chain(*distinct_types))))
                assert atoms.issubset([m[0] for m in mass_info]), \
                    "Undefined atom type found in %s" % k

        self.mass_info = mass_info
        self.pair_coeffs = pair_coeffs
        self.mol_coeffs = mol_coeffs
        masses_sec, self.atom_map = self.get_coeffs_and_mapper("Masses")
        self.masses = masses_sec["Masses"]

    def get_coeffs_and_mapper(self, section):
        """
        Returns data for Masses or a force field section for molecular
        topology, also returns a mapper dict ({type: id, ...}) for
        labeling data in Atoms or the relative topology section.

        Args:
            section (str): Section title. Choose among "Masses",
                "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", and
                "Improper Coeffs".

        Returns:
            Dict with section title as key for the usage of
            LammpsData, and a mapper dict for labeling
            {"XX Coeffs": [{"id": 1, "coeffs": coeffs}, ...]},
            {"type1": 1, ...}

        """
        data = []
        mapper = {}
        if section == "Masses":
            for i, (k, v) in enumerate(self.mass_info):
                data.append({"id": i + 1, "mass": v})
                mapper[k] = i + 1
        elif section in SECTION_KEYWORDS["ff"][2:]:
            for i, d in enumerate(self.mol_coeffs[section]):
                data.append({"id": i + 1, "coeffs": d["coeffs"]})
                mapper.update({k: i + 1 for k in d["types"]})
        else:
            raise RuntimeError("Invalid coefficient section keyword")
        return {section: data}, mapper

    def get_pair_coeffs(self):
        """
        Returns data for Pair(IJ) Coeffs section.

        Returns:
            Dict with section title as key for the usage of
            LammpsData
            {"Pair Coeffs": [{"id": 1, "coeffs": coeffs}, ...]} or
            {"PairIJ Coeffs": [{"id1": 1, "id2": 1,
                                "coeffs": coeffs}, ...]}

        """
        if self.pair_type == "pair":
            return {"Pair Coeffs": [{"id": i + 1, "coeffs": c}
                                    for i, c in enumerate(self.pair_coeffs)]}
        elif self.pair_type == "pairij":
            n = len(self.mass_info)
            ids = itertools.combinations_with_replacement(range(1, n + 1), 2)
            return {"PairIJ Coeffs": [{"id1": i[0], "id2": i[1], "coeffs": c}
                                      for i, c in zip(ids, self.pair_coeffs)]}

    def to_file(self, filename):
        """
        Saves object to a file in YAML format.

        Args:
            filename (str): File name.

        """
        d = {"mass_info": self.mass_info, "pair_coeffs": self.pair_coeffs,
             "mol_coeffs": self.mol_coeffs}
        yaml = YAML(typ="safe")
        with open(filename, "w") as f:
            yaml.dump(d, f)

    @classmethod
    def from_file(cls, filename):
        """
        Constructor that reads in a file in YAML format.

        Args:
            filename (str): File name.

        """
        yaml = YAML(typ="safe")
        with open(filename, "r") as f:
            d = yaml.load(f)
        return cls.from_dict(d)

    @classmethod
    def from_dict(cls, d):
        d["mass_info"] = [tuple(m) for m in d["mass_info"]]
        if d.get("mol_coeffs"):
            for v in d["mol_coeffs"].values():
                for c in v:
                    c["types"] = [tuple(t) for t in c["types"]]
        return cls(d["mass_info"], d["pair_coeffs"], d["mol_coeffs"])
