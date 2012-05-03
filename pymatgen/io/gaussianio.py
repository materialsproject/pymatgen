#!/usr/bin/env python

'''
This module implements input and output processing from molecules to Gaussian
Input files.
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Apr 17, 2012"

import re
import math

from pymatgen.core.operations import SymmOp
import numpy as np
from pymatgen.core.structure import Molecule

class GaussianInput(object):
    """
    An object representing a Gaussian input file.
    """
    def __init__(self, mol, charge=0, spin_multiplicity=1, title=None,
                 functional="HF", basis_set="6-31G(d)", route_parameters=None,
                 input_parameters=None):
        """
        Args:
            mol:
                Input molecule
            charge:
                Charge of the molecule. Defaults to 0.
            spin_multiplicity:
                Spin multiplicity of molecule. Defaults to 1.
            title:
                Title for run. Defaults to formula of molecule if None.
            functional:
                Functional for run.
            basis_set:
                Basis set for run.
            route_parameters:
                Additional route parameters as a dict. For example, 
                {'SP':"", "SCF":"Tight"}
            input_parameters:
                Additional input parameters for run as a dict. Used for example,
                in PCM calculations.  E.g., {"EPS":12}
        """
        self._mol = mol
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        self.functional = functional
        self.basis_set = basis_set
        self.route_parameters = route_parameters if route_parameters else {}
        self.input_parameters = input_parameters if input_parameters else {}
        self.title = title if title else self._mol.composition.formula

    @property
    def molecule(self):
        """
        Returns molecule associated with this GaussianInput.
        """
        return self._mol

    @staticmethod
    def parse_coords(coord_lines):
        """
        Helper method to parse coordinates.
        """
        paras = {}
        var_pattern = re.compile("^\s*([A-Za-z]+\S*)[\s=,]+([\d\-\.]+)\s*$")
        for l in coord_lines:
            m = var_pattern.match(l.strip())
            if m:
                paras[m.group(1)] = float(m.group(2))
        zmat_patt = re.compile("^\s*([A-Za-z]+)[\w\d\-\_]*([\s,]+(\w+)[\s,]+(\w+))*[\-\.\s,\w]*$")
        mixed_species_patt = re.compile("([A-Za-z]+)[\d\-\_]+")
        xyz_patt = re.compile("^\s*([A-Za-z]+[\w\d\-\_]*)\s+([\d\.eE\-]+)\s+([\d\.eE\-]+)\s+([\d\.eE\-]+)[\-\.\s,\w.]*$")

        parsed_species = []
        species = []
        coords = []

        for l in coord_lines:
            l = l.strip()
            if l == "":
                break
            if xyz_patt.match(l):
                m = xyz_patt.match(l)
                m2 = mixed_species_patt.match(m.group(1))
                if m2:
                    parsed_species.append(m.group(1))
                    species.append(m2.group(1))
                else:
                    species.append(m.group(1))

                toks = re.split("[,\s]+", l.strip())
                if len(toks) > 4:
                    coords.append(map(float, toks[2:5]))
                else:
                    coords.append(map(float, toks[1:4]))
            elif zmat_patt.match(l):
                toks = re.split("[,\s]+", l.strip())
                m = mixed_species_patt.match(toks[0])
                if m:
                    parsed_species.append(toks[0])
                    species.append(m.group(1))
                else:
                    species.append(toks[0])
                toks.pop(0)
                if len(toks) == 0:
                    coords.append(np.array([0, 0, 0]))
                else:
                    nn = []
                    parameters = []
                    while len(toks) > 0:
                        ind = toks.pop(0)
                        data = toks.pop(0)
                        try:
                            int(ind)
                            nn.append(int(ind))
                        except:
                            nn.append(parsed_species.index(ind))
                        parameters.append(paras[data])
                    if len(nn) == 1:
                        coords.append(np.array([0, 0, parameters[0]]))
                    elif len(nn) == 2:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        axis = [0, 1, 0]
                        coord = SymmOp.from_origin_axis_angle(coords1, axis, angle, False).operate(coords2)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)
                    elif len(nn) == 3:
                        coords1 = coords[nn[0] - 1]
                        coords2 = coords[nn[1] - 1]
                        coords3 = coords[nn[2] - 1]
                        bl = parameters[0]
                        angle = parameters[1]
                        dih = parameters[2]
                        v1 = coords3 - coords2
                        v2 = coords1 - coords2
                        axis = np.cross(v1, v2)
                        coord = SymmOp.from_origin_axis_angle(coords1, axis, angle, False).operate(coords2)
                        v1 = coord - coords1
                        v2 = coords1 - coords2
                        v3 = np.cross(v1, v2)
                        d = np.dot(v3, axis) / np.linalg.norm(v3) / np.linalg.norm(axis)
                        if d > 1:
                            d = 1
                        elif d < -1:
                            d = -1
                        adj = math.acos(d) * 180 / math.pi
                        axis = coords1 - coords2
                        coord = SymmOp.from_origin_axis_angle(coords1, axis, dih - adj, False).operate(coord)
                        vec = coord - coords1
                        coord = vec * bl / np.linalg.norm(vec) + coords1
                        coords.append(coord)

        return Molecule(species, coords)

    @staticmethod
    def from_string(contents):
        """
        Creates GaussianInput from a string.
        
        Args:
            contents:
                String representing an Gaussian input file.
        
        Returns:
            GaussianInput object
        """
        lines = contents.split("\n")
        route_patt = re.compile("^#[sSpPnN]*.*")

        route = None
        for i, l in enumerate(lines):
            if route_patt.match(l):
                route = l
                route_index = i
                break
        route_paras = {}
        if route:
            for tok in re.split("\s+", route):
                if tok.strip().startswith("#"):
                    continue
                if re.match("\w+\/.*", tok):
                    d = tok.split("/")
                    functional = d[0]
                    basis_set = d[1]
                else:
                    d = tok.split("=")
                    v = None if len(d) == 1 else d[1]
                    route_paras[d[0]] = v
        title = lines[route_index + 2]
        toks = re.split("\s+", lines[route_index + 4])
        charge = int(toks[0])
        spin_mult = int(toks[1])
        coord_lines = []
        spaces = 0
        input_paras = {}
        for i in xrange(route_index + 5, len(lines)):
            if lines[i].strip() == "":
                spaces += 1
            if spaces >= 2:
                d = lines[i].split("=")
                if len(d) == 2:
                    input_paras[d[0]] = float(d[1])
            else:
                coord_lines.append(lines[i].strip())
        mol = GaussianInput.parse_coords(coord_lines)

        return GaussianInput(mol, charge=charge, spin_multiplicity=spin_mult,
                             title=title, functional=functional,
                            basis_set=basis_set, route_parameters=route_paras,
                            input_parameters=input_paras)

    @staticmethod
    def from_file(filename):
        """
        Creates GaussianInput from a file.
        
        Args:
            filename:
                Gaussian input filename
                
        Returns:
            GaussianInput object
        """
        with open(filename, "r") as f:
            return GaussianInput.from_string(f.read())

    def _find_nn_pos_before_site(self, siteindex):
        """
        Returns index of nearest neighbor atoms.
        """
        alldist = [(self._mol.get_distance(siteindex, i), i) for i in xrange(siteindex)]
        alldist = sorted(alldist, key=lambda x:x[0])
        return [d[1] for d in alldist]

    def get_zmatrix(self):
        """
        Returns a z-matrix representation of the molecule.
        """
        output = []
        outputvar = []
        for i, site in enumerate(self._mol):
            if i == 0:
                output.append("{}".format(site.specie))
            elif i == 1:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                output.append("{} {} B{}".format(self._mol[i].specie, nn[0] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
            elif i == 2:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                angle = self._mol.get_angle(i, nn[0], nn[1])
                output.append("{} {} B{} {} A{}".format(self._mol[i].specie, nn[0] + 1, i, nn[1] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
                outputvar.append("A{}={:.6f}".format(i, angle))
            else:
                nn = self._find_nn_pos_before_site(i)
                bondlength = self._mol.get_distance(i, nn[0])
                angle = self._mol.get_angle(i, nn[0], nn[1])
                dih = self._mol.get_dihedral(i, nn[0], nn[1], nn[2])
                output.append("{} {} B{} {} A{} {} D{}".format(self._mol[i].specie, nn[0] + 1, i, nn[1] + 1, i, nn[2] + 1, i))
                outputvar.append("B{}={:.6f}".format(i, bondlength))
                outputvar.append("A{}={:.6f}".format(i, angle))
                outputvar.append("D{}={:.6f}".format(i, dih))
        return "\n".join(output) + "\n\n" + "\n".join(outputvar)

    def __str__(self):
        def para_dict_to_string(para, joiner=" "):
            para_str = []
            for k, v in para.items():
                if v:
                    para_str.append("{}={}".format(k, v))
                else:
                    para_str.append(k)
            return joiner.join(para_str)

        output = ["#P {func}/{bset} {route} Test".format(func=self.functional, bset=self.basis_set, route=para_dict_to_string(self.route_parameters))]
        output.append("")
        output.append(self.title)
        output.append("")
        output.append("{} {}".format(self.charge, self.spin_multiplicity))
        output.append(self.get_zmatrix())
        output.append("")
        output.append(para_dict_to_string(self.input_parameters, "\n"))
        output.append("")
        return "\n".join(output)

