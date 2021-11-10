# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


"""
This module provides classes to perform fitting of molecule with arbitrary
atom orders.
This module is supposed to perform exact comparisons without the atom order
correspondence prerequisite, while molecule_structure_comparator is supposed
to do rough comparisons with the atom order correspondence prerequisite.

The implementation is based on an excellent python package called `rmsd` that
you can find at https://github.com/charnley/rmsd.
"""

__author__ = "Xiaohui Qu, Adam Fekete"
__version__ = "1.0"
__email__ = "xhqu1981@gmail.com"


import abc
import copy
import itertools
import logging
import math
import re

import numpy as np
from monty.dev import requires
from monty.json import MSONable

try:
    from openbabel import openbabel as ob

    from pymatgen.io.babel import BabelMolAdaptor
except ImportError:
    ob = None

from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from pymatgen.core.structure import Molecule  # pylint: disable=ungrouped-imports

logger = logging.getLogger(__name__)


class AbstractMolAtomMapper(MSONable, metaclass=abc.ABCMeta):
    """
    Abstract molecular atom order mapping class. A mapping will be able to
    find the uniform atom order of two molecules that can pair the
    geometrically equivalent atoms.
    """

    @abc.abstractmethod
    def uniform_labels(self, mol1, mol2):
        """
        Pair the geometrically equivalent atoms of the molecules.

        Args:
            mol1: First molecule. OpenBabel OBMol or pymatgen Molecule object.
            mol2: Second molecule. OpenBabel OBMol or pymatgen Molecule object.

        Returns:
            (list1, list2) if uniform atom order is found. list1 and list2
            are for mol1 and mol2, respectively. Their length equal
            to the number of atoms. They represents the uniform atom order
            of the two molecules. The value of each element is the original
            atom index in mol1 or mol2 of the current atom in uniform atom
            order.
            (None, None) if unform atom is not available.
        """
        pass

    @abc.abstractmethod
    def get_molecule_hash(self, mol):
        """
        Defines a hash for molecules. This allows molecules to be grouped
        efficiently for comparison.

        Args:
            mol: The molecule. OpenBabel OBMol or pymatgen Molecule object

        Returns:
            A hashable object. Examples can be string formulas, etc.
        """
        pass

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (): Dict

        Returns:
            AbstractMolAtomMapper
        """
        for trans_modules in ["molecule_matcher"]:
            import sys

            if sys.version_info > (3, 0):
                level = 0  # Python 3.x
            else:
                level = -1  # Python 2.x
            mod = __import__(
                "pymatgen.analysis." + trans_modules,
                globals(),
                locals(),
                [d["@class"]],
                level,
            )
            if hasattr(mod, d["@class"]):
                class_proxy = getattr(mod, d["@class"])
                from_dict_proxy = getattr(class_proxy, "from_dict")
                return from_dict_proxy(d)
        raise ValueError("Invalid Comparator dict")


class IsomorphismMolAtomMapper(AbstractMolAtomMapper):
    """
    Pair atoms by isomorphism permutations in the OpenBabel::OBAlign class
    """

    def uniform_labels(self, mol1, mol2):
        """
        Pair the geometrically equivalent atoms of the molecules.
        Calculate RMSD on all possible isomorphism mappings and return mapping
        with the least RMSD

        Args:
            mol1: First molecule. OpenBabel OBMol or pymatgen Molecule object.
            mol2: Second molecule. OpenBabel OBMol or pymatgen Molecule object.

        Returns:
            (list1, list2) if uniform atom order is found. list1 and list2
            are for mol1 and mol2, respectively. Their length equal
            to the number of atoms. They represents the uniform atom order
            of the two molecules. The value of each element is the original
            atom index in mol1 or mol2 of the current atom in uniform atom
            order.
            (None, None) if unform atom is not available.
        """
        obmol1 = BabelMolAdaptor(mol1).openbabel_mol
        obmol2 = BabelMolAdaptor(mol2).openbabel_mol

        h1 = self.get_molecule_hash(obmol1)
        h2 = self.get_molecule_hash(obmol2)
        if h1 != h2:
            return None, None

        query = ob.CompileMoleculeQuery(obmol1)
        isomapper = ob.OBIsomorphismMapper.GetInstance(query)
        isomorph = ob.vvpairUIntUInt()
        isomapper.MapAll(obmol2, isomorph)

        sorted_isomorph = [sorted(x, key=lambda morp: morp[0]) for x in isomorph]
        label2_list = tuple(tuple(p[1] + 1 for p in x) for x in sorted_isomorph)

        vmol1 = obmol1
        aligner = ob.OBAlign(True, False)
        aligner.SetRefMol(vmol1)
        least_rmsd = float("Inf")
        best_label2 = None
        label1 = list(range(1, obmol1.NumAtoms() + 1))
        # noinspection PyProtectedMember
        elements1 = InchiMolAtomMapper._get_elements(vmol1, label1)
        for label2 in label2_list:
            # noinspection PyProtectedMember
            elements2 = InchiMolAtomMapper._get_elements(obmol2, label2)
            if elements1 != elements2:
                continue
            vmol2 = ob.OBMol()
            for i in label2:
                vmol2.AddAtom(obmol2.GetAtom(i))
            aligner.SetTargetMol(vmol2)
            aligner.Align()
            rmsd = aligner.GetRMSD()
            if rmsd < least_rmsd:
                least_rmsd = rmsd
                best_label2 = copy.copy(label2)
        return label1, best_label2

    def get_molecule_hash(self, mol):
        """
        Return inchi as molecular hash
        """
        obconv = ob.OBConversion()
        obconv.SetOutFormat(str("inchi"))
        obconv.AddOption(str("X"), ob.OBConversion.OUTOPTIONS, str("DoNotAddH"))
        inchi_text = obconv.WriteString(mol)
        match = re.search(r"InChI=(?P<inchi>.+)\n", inchi_text)
        return match.group("inchi")

    def as_dict(self):
        """
        Returns:
            Jsonable dict.
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            IsomorphismMolAtomMapper
        """
        return IsomorphismMolAtomMapper()


class InchiMolAtomMapper(AbstractMolAtomMapper):
    """
    Pair atoms by inchi labels.
    """

    def __init__(self, angle_tolerance=10.0):
        """
        Args:
            angle_tolerance (float): Angle threshold to assume linear molecule. In degrees.
        """
        self._angle_tolerance = angle_tolerance
        self._assistant_mapper = IsomorphismMolAtomMapper()

    def as_dict(self):
        """
        Returns:
            MSONAble dict.
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "angle_tolerance": self._angle_tolerance,
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict Representation

        Returns:
            InchiMolAtomMapper
        """
        return InchiMolAtomMapper(angle_tolerance=d["angle_tolerance"])

    @staticmethod
    def _inchi_labels(mol):
        """
        Get the inchi canonical labels of the heavy atoms in the molecule

        Args:
            mol: The molecule. OpenBabel OBMol object

        Returns:
            The label mappings. List of tuple of canonical label,
            original label
            List of equivalent atoms.
        """
        obconv = ob.OBConversion()
        obconv.SetOutFormat(str("inchi"))
        obconv.AddOption(str("a"), ob.OBConversion.OUTOPTIONS)
        obconv.AddOption(str("X"), ob.OBConversion.OUTOPTIONS, str("DoNotAddH"))
        inchi_text = obconv.WriteString(mol)
        match = re.search(
            r"InChI=(?P<inchi>.+)\nAuxInfo=.+" r"/N:(?P<labels>[0-9,;]+)/(E:(?P<eq_atoms>[0-9," r";\(\)]*)/)?",
            inchi_text,
        )
        inchi = match.group("inchi")
        label_text = match.group("labels")
        eq_atom_text = match.group("eq_atoms")
        heavy_atom_labels = tuple(int(i) for i in label_text.replace(";", ",").split(","))
        eq_atoms = []
        if eq_atom_text is not None:
            eq_tokens = re.findall(r"\(((?:[0-9]+,)+[0-9]+)\)", eq_atom_text.replace(";", ","))
            eq_atoms = tuple(tuple(int(i) for i in t.split(",")) for t in eq_tokens)
        return heavy_atom_labels, eq_atoms, inchi

    @staticmethod
    def _group_centroid(mol, ilabels, group_atoms):
        """
        Calculate the centroids of a group atoms indexed by the labels of inchi

        Args:
            mol: The molecule. OpenBabel OBMol object
            ilabel: inchi label map

        Returns:
            Centroid. Tuple (x, y, z)
        """
        c1x, c1y, c1z = 0.0, 0.0, 0.0
        for i in group_atoms:
            orig_idx = ilabels[i - 1]
            oa1 = mol.GetAtom(orig_idx)
            c1x += float(oa1.x())
            c1y += float(oa1.y())
            c1z += float(oa1.z())
        num_atoms = len(group_atoms)
        c1x /= num_atoms
        c1y /= num_atoms
        c1z /= num_atoms
        return c1x, c1y, c1z

    def _virtual_molecule(self, mol, ilabels, eq_atoms):
        """
        Create a virtual molecule by unique atoms, the centriods of the
        equivalent atoms

        Args:
            mol: The molecule. OpenBabel OBMol object
            ilables: inchi label map
            eq_atoms: equivalent atom labels
            farthest_group_idx: The equivalent atom group index in which
                there is the farthest atom to the centroid

        Return:
            The virtual molecule
        """
        vmol = ob.OBMol()

        non_unique_atoms = {a for g in eq_atoms for a in g}
        all_atoms = set(range(1, len(ilabels) + 1))
        unique_atom_labels = sorted(all_atoms - non_unique_atoms)

        # try to align molecules using unique atoms
        for i in unique_atom_labels:
            orig_idx = ilabels[i - 1]
            oa1 = mol.GetAtom(orig_idx)
            a1 = vmol.NewAtom()
            a1.SetAtomicNum(oa1.GetAtomicNum())
            a1.SetVector(oa1.GetVector())

        # try to align using centroids of the equivalent atoms
        if vmol.NumAtoms() < 3:
            for symm in eq_atoms:
                c1x, c1y, c1z = self._group_centroid(mol, ilabels, symm)
                min_distance = float("inf")
                for i in range(1, vmol.NumAtoms() + 1):
                    va = vmol.GetAtom(i)
                    distance = math.sqrt((c1x - va.x()) ** 2 + (c1y - va.y()) ** 2 + (c1z - va.z()) ** 2)
                    if distance < min_distance:
                        min_distance = distance
                if min_distance > 0.2:
                    a1 = vmol.NewAtom()
                    a1.SetAtomicNum(9)
                    a1.SetVector(c1x, c1y, c1z)

        return vmol

    @staticmethod
    def _align_heavy_atoms(mol1, mol2, vmol1, vmol2, ilabel1, ilabel2, eq_atoms):
        """
        Align the label of topologically identical atoms of second molecule
        towards first molecule

        Args:
            mol1: First molecule. OpenBabel OBMol object
            mol2: Second molecule. OpenBabel OBMol object
            vmol1: First virtual molecule constructed by centroids. OpenBabel
                OBMol object
            vmol2: First virtual molecule constructed by centroids. OpenBabel
                OBMol object
            ilabel1: inchi label map of the first molecule
            ilabel2: inchi label map of the second molecule
            eq_atoms: equivalent atom lables

        Return:
            corrected inchi labels of heavy atoms of the second molecule
        """

        nvirtual = vmol1.NumAtoms()
        nheavy = len(ilabel1)

        for i in ilabel2:  # add all heavy atoms
            a1 = vmol1.NewAtom()
            a1.SetAtomicNum(1)
            a1.SetVector(0.0, 0.0, 0.0)  # useless, just to pair with vmol2
            oa2 = mol2.GetAtom(i)
            a2 = vmol2.NewAtom()
            a2.SetAtomicNum(1)
            # align using the virtual atoms, these atoms are not
            # used to align, but match by positions
            a2.SetVector(oa2.GetVector())

        aligner = ob.OBAlign(False, False)
        aligner.SetRefMol(vmol1)
        aligner.SetTargetMol(vmol2)
        aligner.Align()
        aligner.UpdateCoords(vmol2)

        canon_mol1 = ob.OBMol()
        for i in ilabel1:
            oa1 = mol1.GetAtom(i)
            a1 = canon_mol1.NewAtom()
            a1.SetAtomicNum(oa1.GetAtomicNum())
            a1.SetVector(oa1.GetVector())

        aligned_mol2 = ob.OBMol()
        for i in range(nvirtual + 1, nvirtual + nheavy + 1):
            oa2 = vmol2.GetAtom(i)
            a2 = aligned_mol2.NewAtom()
            a2.SetAtomicNum(oa2.GetAtomicNum())
            a2.SetVector(oa2.GetVector())

        canon_label2 = list(range(1, nheavy + 1))
        for symm in eq_atoms:
            for i in symm:
                canon_label2[i - 1] = -1
        for symm in eq_atoms:
            candidates1 = list(symm)
            candidates2 = list(symm)
            for c2 in candidates2:
                distance = 99999.0
                canon_idx = candidates1[0]
                a2 = aligned_mol2.GetAtom(c2)
                for c1 in candidates1:
                    a1 = canon_mol1.GetAtom(c1)
                    d = a1.GetDistance(a2)
                    if d < distance:
                        distance = d
                        canon_idx = c1
                canon_label2[c2 - 1] = canon_idx
                candidates1.remove(canon_idx)

        canon_inchi_orig_map2 = list(zip(canon_label2, list(range(1, nheavy + 1)), ilabel2))
        canon_inchi_orig_map2.sort(key=lambda m: m[0])
        heavy_atom_indices2 = tuple(x[2] for x in canon_inchi_orig_map2)
        return heavy_atom_indices2

    @staticmethod
    def _align_hydrogen_atoms(mol1, mol2, heavy_indices1, heavy_indices2):
        """
        Align the label of topologically identical atoms of second molecule
        towards first molecule

        Args:
            mol1: First molecule. OpenBabel OBMol object
            mol2: Second molecule. OpenBabel OBMol object
            heavy_indices1: inchi label map of the first molecule
            heavy_indices2: label map of the second molecule

        Return:
            corrected label map of all atoms of the second molecule
        """
        num_atoms = mol2.NumAtoms()
        all_atom = set(range(1, num_atoms + 1))
        hydrogen_atoms1 = all_atom - set(heavy_indices1)
        hydrogen_atoms2 = all_atom - set(heavy_indices2)
        label1 = heavy_indices1 + tuple(hydrogen_atoms1)
        label2 = heavy_indices2 + tuple(hydrogen_atoms2)

        cmol1 = ob.OBMol()
        for i in label1:
            oa1 = mol1.GetAtom(i)
            a1 = cmol1.NewAtom()
            a1.SetAtomicNum(oa1.GetAtomicNum())
            a1.SetVector(oa1.GetVector())
        cmol2 = ob.OBMol()
        for i in label2:
            oa2 = mol2.GetAtom(i)
            a2 = cmol2.NewAtom()
            a2.SetAtomicNum(oa2.GetAtomicNum())
            a2.SetVector(oa2.GetVector())

        aligner = ob.OBAlign(False, False)
        aligner.SetRefMol(cmol1)
        aligner.SetTargetMol(cmol2)
        aligner.Align()
        aligner.UpdateCoords(cmol2)

        hydrogen_label2 = []
        hydrogen_label1 = list(range(len(heavy_indices1) + 1, num_atoms + 1))
        for h2 in range(len(heavy_indices2) + 1, num_atoms + 1):
            distance = 99999.0
            idx = hydrogen_label1[0]
            a2 = cmol2.GetAtom(h2)
            for h1 in hydrogen_label1:
                a1 = cmol1.GetAtom(h1)
                d = a1.GetDistance(a2)
                if d < distance:
                    distance = d
                    idx = h1
            hydrogen_label2.append(idx)
            hydrogen_label1.remove(idx)

        hydrogen_orig_idx2 = label2[len(heavy_indices2) :]
        hydrogen_canon_orig_map2 = list(zip(hydrogen_label2, hydrogen_orig_idx2))
        hydrogen_canon_orig_map2.sort(key=lambda m: m[0])
        hydrogen_canon_indices2 = [x[1] for x in hydrogen_canon_orig_map2]

        canon_label1 = label1
        canon_label2 = heavy_indices2 + tuple(hydrogen_canon_indices2)

        return canon_label1, canon_label2

    @staticmethod
    def _get_elements(mol, label):
        """
        The the elements of the atoms in the specified order

        Args:
            mol: The molecule. OpenBabel OBMol object.
            label: The atom indices. List of integers.

        Returns:
            Elements. List of integers.
        """
        elements = [int(mol.GetAtom(i).GetAtomicNum()) for i in label]
        return elements

    def _is_molecule_linear(self, mol):
        """
        Is the molecule a linear one

        Args:
            mol: The molecule. OpenBabel OBMol object.

        Returns:
            Boolean value.
        """
        if mol.NumAtoms() < 3:
            return True
        a1 = mol.GetAtom(1)
        a2 = mol.GetAtom(2)
        for i in range(3, mol.NumAtoms() + 1):
            angle = float(mol.GetAtom(i).GetAngle(a2, a1))
            if angle < 0.0:
                angle = -angle
            if angle > 90.0:
                angle = 180.0 - angle
            if angle > self._angle_tolerance:
                return False
        return True

    def uniform_labels(self, mol1, mol2):
        """
        Args:
            mol1 (Molecule): Molecule 1
            mol2 (Molecule): Molecule 2

        Returns:
            Labels
        """
        obmol1 = BabelMolAdaptor(mol1).openbabel_mol
        obmol2 = BabelMolAdaptor(mol2).openbabel_mol

        ilabel1, iequal_atom1, inchi1 = self._inchi_labels(obmol1)
        ilabel2, iequal_atom2, inchi2 = self._inchi_labels(obmol2)

        if inchi1 != inchi2:
            return None, None  # Topoligically different

        if iequal_atom1 != iequal_atom2:
            raise Exception("Design Error! Equavilent atoms are inconsistent")

        vmol1 = self._virtual_molecule(obmol1, ilabel1, iequal_atom1)
        vmol2 = self._virtual_molecule(obmol2, ilabel2, iequal_atom2)

        if vmol1.NumAtoms() != vmol2.NumAtoms():
            return None, None

        if vmol1.NumAtoms() < 3 or self._is_molecule_linear(vmol1) or self._is_molecule_linear(vmol2):
            # using isomorphism for difficult (actually simple) molecules
            clabel1, clabel2 = self._assistant_mapper.uniform_labels(mol1, mol2)
        else:
            heavy_atom_indices2 = self._align_heavy_atoms(obmol1, obmol2, vmol1, vmol2, ilabel1, ilabel2, iequal_atom1)
            clabel1, clabel2 = self._align_hydrogen_atoms(obmol1, obmol2, ilabel1, heavy_atom_indices2)
        if clabel1 and clabel2:
            elements1 = self._get_elements(obmol1, clabel1)
            elements2 = self._get_elements(obmol2, clabel2)

            if elements1 != elements2:
                return None, None

        return clabel1, clabel2

    def get_molecule_hash(self, mol):
        """
        Return inchi as molecular hash
        """
        obmol = BabelMolAdaptor(mol).openbabel_mol
        inchi = self._inchi_labels(obmol)[2]
        return inchi


class MoleculeMatcher(MSONable):
    """
    Class to match molecules and identify whether molecules are the same.
    """

    @requires(
        ob,
        "BabelMolAdaptor requires openbabel to be installed with "
        "Python bindings. Please get it at http://openbabel.org "
        "(version >=3.0.0).",
    )
    def __init__(self, tolerance=0.01, mapper=InchiMolAtomMapper()):
        """
        Args:
            tolerance (float): RMSD difference threshold whether two molecules are
                different
            mapper (AbstractMolAtomMapper): MolAtomMapper object that is able to map the atoms of two
                molecule to uniform order
        """
        self._tolerance = tolerance
        self._mapper = mapper

    def fit(self, mol1, mol2):
        """
        Fit two molecules.

        Args:
            mol1: First molecule. OpenBabel OBMol or pymatgen Molecule object
            mol2: Second molecule. OpenBabel OBMol or pymatgen Molecule object

        Returns:
            A boolean value indicates whether two molecules are the same.
        """
        return self.get_rmsd(mol1, mol2) < self._tolerance

    def get_rmsd(self, mol1, mol2):
        """
        Get RMSD between two molecule with arbitrary atom order.

        Returns:
            RMSD if topology of the two molecules are the same
            Infinite if  the topology is different
        """
        label1, label2 = self._mapper.uniform_labels(mol1, mol2)
        if label1 is None or label2 is None:
            return float("Inf")
        return self._calc_rms(mol1, mol2, label1, label2)

    @staticmethod
    def _calc_rms(mol1, mol2, clabel1, clabel2):
        """
        Calculate the RMSD.

        Args:
            mol1: The first molecule. OpenBabel OBMol or pymatgen Molecule
                object
            mol2: The second molecule. OpenBabel OBMol or pymatgen Molecule
                object
            clabel1: The atom indices that can reorder the first molecule to
                uniform atom order
            clabel1: The atom indices that can reorder the second molecule to
                uniform atom order

        Returns:
            The RMSD.
        """
        obmol1 = BabelMolAdaptor(mol1).openbabel_mol
        obmol2 = BabelMolAdaptor(mol2).openbabel_mol

        cmol1 = ob.OBMol()
        for i in clabel1:
            oa1 = obmol1.GetAtom(i)
            a1 = cmol1.NewAtom()
            a1.SetAtomicNum(oa1.GetAtomicNum())
            a1.SetVector(oa1.GetVector())
        cmol2 = ob.OBMol()
        for i in clabel2:
            oa2 = obmol2.GetAtom(i)
            a2 = cmol2.NewAtom()
            a2.SetAtomicNum(oa2.GetAtomicNum())
            a2.SetVector(oa2.GetVector())

        aligner = ob.OBAlign(True, False)
        aligner.SetRefMol(cmol1)
        aligner.SetTargetMol(cmol2)
        aligner.Align()
        return aligner.GetRMSD()

    def group_molecules(self, mol_list):
        """
        Group molecules by structural equality.

        Args:
            mol_list: List of OpenBabel OBMol or pymatgen objects

        Returns:
            A list of lists of matched molecules
            Assumption: if s1=s2 and s2=s3, then s1=s3
            This may not be true for small tolerances.
        """
        mol_hash = [(i, self._mapper.get_molecule_hash(m)) for i, m in enumerate(mol_list)]
        mol_hash.sort(key=lambda x: x[1])

        # Use molecular hash to pre-group molecules.
        raw_groups = tuple(tuple(m[0] for m in g) for k, g in itertools.groupby(mol_hash, key=lambda x: x[1]))

        group_indices = []
        for rg in raw_groups:
            mol_eq_test = [
                (p[0], p[1], self.fit(mol_list[p[0]], mol_list[p[1]])) for p in itertools.combinations(sorted(rg), 2)
            ]
            mol_eq = {(p[0], p[1]) for p in mol_eq_test if p[2]}
            not_alone_mols = set(itertools.chain.from_iterable(mol_eq))
            alone_mols = set(rg) - not_alone_mols
            group_indices.extend([[m] for m in alone_mols])
            while len(not_alone_mols) > 0:
                current_group = {not_alone_mols.pop()}
                while len(not_alone_mols) > 0:
                    candidate_pairs = {tuple(sorted(p)) for p in itertools.product(current_group, not_alone_mols)}
                    mutual_pairs = candidate_pairs & mol_eq
                    if len(mutual_pairs) == 0:
                        break
                    mutual_mols = set(itertools.chain.from_iterable(mutual_pairs))
                    current_group |= mutual_mols
                    not_alone_mols -= mutual_mols
                group_indices.append(sorted(current_group))

        group_indices.sort(key=lambda x: (len(x), -x[0]), reverse=True)
        all_groups = [[mol_list[i] for i in g] for g in group_indices]
        return all_groups

    def as_dict(self):
        """
        Returns:
            MSONAble dict.
        """
        return {
            "version": __version__,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "tolerance": self._tolerance,
            "mapper": self._mapper.as_dict(),
        }

    @classmethod
    def from_dict(cls, d):
        """
        Args:
            d (dict): Dict representation

        Returns:
            MoleculeMatcher
        """
        return MoleculeMatcher(
            tolerance=d["tolerance"],
            mapper=AbstractMolAtomMapper.from_dict(d["mapper"]),
        )


class KabschMatcher(MSONable):
    """Molecule matcher using Kabsch algorithm

    The Kabsch algorithm capable aligning two molecules by finding the parameters
    (translation, rotation) which minimize the root-mean-square-deviation (RMSD) of
    two molecules which are topologically (atom types, geometry) similar two each other.

    Notes:
        When aligning molecules, the atoms of the two molecules **must** be in the same
        order for the results to be sensible.
    """

    def __init__(self, target: Molecule):
        """Constructor of the matcher object.

        Args:
            target: a `Molecule` object used as a target during the alignment
        """
        self.target = target

    def match(self, p: Molecule):
        """Using the Kabsch algorithm the alignment of two molecules (P, Q)
        happens in three steps:
        - translate the P and Q into their centroid
        - compute of the optimal rotation matrix (U) using Kabsch algorithm
        - compute the translation (V) and rmsd

        The function returns the rotation matrix (U), translation vector (V),
        and RMSD between Q and P', where P' is:

            P' = P * U + V

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            U: Rotation matrix (D,D)
            V: Translation vector (D)
            RMSD : Root mean squared deviation between P and Q
        """
        if self.target.atomic_numbers != p.atomic_numbers:
            raise ValueError("The order of the species aren't matching! Please try using `PermInvMatcher`.")

        p_coord, q_coord = p.cart_coords, self.target.cart_coords

        # Both sets of coordinates must be translated first, so that their
        # centroid coincides with the origin of the coordinate system.
        p_trans, q_trans = p_coord.mean(axis=0), q_coord.mean(axis=0)
        p_centroid, q_centroid = p_coord - p_trans, q_coord - q_trans

        # The optimal rotation matrix U using Kabsch algorithm
        U = self.kabsch(p_centroid, q_centroid)

        p_prime_centroid = np.dot(p_centroid, U)
        rmsd = np.sqrt(np.mean(np.square(p_prime_centroid - q_centroid)))

        V = q_trans - np.dot(p_trans, U)

        return U, V, rmsd

    def fit(self, p: Molecule):
        """Rotate and transform `p` molecule according to the best match.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            p_prime: Rotated and translated of the `p` `Molecule` object
            rmsd: Root-mean-square-deviation between `p_prime` and the `target`
        """
        U, V, rmsd = self.match(p)

        # Rotate and translate matrix `p` onto the target molecule.
        # P' = P * U + V
        p_prime = p.copy()
        for site in p_prime:
            site.coords = np.dot(site.coords, U) + V

        return p_prime, rmsd

    @staticmethod
    def kabsch(P: np.ndarray, Q: np.ndarray):
        """The Kabsch algorithm is a method for calculating the optimal rotation matrix
        that minimizes the root mean squared deviation (RMSD) between two paired sets of points
        P and Q, centered around the their centroid.

        For more info see:
        - http://en.wikipedia.org/wiki/Kabsch_algorithm and
        - https://cnx.org/contents/HV-RsdwL@23/Molecular-Distance-Measures

        Args:
            P: Nx3 matrix, where N is the number of points.
            Q: Nx3 matrix, where N is the number of points.

        Returns:
            U: 3x3 rotation matrix
        """

        # Computation of the cross-covariance matrix
        C = np.dot(P.T, Q)

        # Computation of the optimal rotation matrix
        # using singular value decomposition (SVD).
        V, S, WT = np.linalg.svd(C)

        # Getting the sign of the det(V*Wt) to decide whether
        d = np.linalg.det(np.dot(V, WT))

        # And finally calculating the optimal rotation matrix R
        # we need to correct our rotation matrix to ensure a right-handed coordinate system.
        U = np.dot(np.dot(V, np.diag([1, 1, d])), WT)

        return U


class BruteForceOrderMatcher(KabschMatcher):
    """Finding the best match between molecules by selecting molecule order
    with the smallest RMSD from all the possible order combinations.

    Notes:
        When aligning molecules, the atoms of the two molecules **must** have same number
        of atoms from the same species.
    """

    def match(self, p: Molecule, ignore_warning=False):
        """Similar as `KabschMatcher.match` but this method also finds the order of
        atoms which belongs to the best match.

        A `ValueError` will be raised when the total number of possible combinations
        become unfeasible (more than a million combination).

        Args:
            p: a `Molecule` object what will be matched with the target one.
            ignore_warning: ignoring error when the number of combination is too large

        Returns:
            inds: The indices of atoms
            U: 3x3 rotation matrix
            V: Translation vector
            rmsd: Root mean squared deviation between P and Q
        """

        q = self.target

        if sorted(p.atomic_numbers) != sorted(q.atomic_numbers):
            raise ValueError("The number of the same species aren't matching!")

        _, count = np.unique(p.atomic_numbers, return_counts=True)
        total_permutations = 1
        for c in count:
            total_permutations *= np.math.factorial(c)  # type: ignore

        if not ignore_warning and total_permutations > 1_000_000:
            raise ValueError(
                "The number of all possible permutations "
                "({}) is not feasible to run this method!".format(total_permutations)
            )

        p_coord, q_coord = p.cart_coords, q.cart_coords
        p_atoms, q_atoms = np.array(p.atomic_numbers), np.array(q.atomic_numbers)

        # Both sets of coordinates must be translated first, so that
        # their centroid coincides with the origin of the coordinate system.
        p_trans, q_trans = p_coord.mean(axis=0), q_coord.mean(axis=0)
        p_centroid, q_centroid = p_coord - p_trans, q_coord - q_trans

        # Sort the order of the target molecule by the elements
        q_inds = np.argsort(q_atoms)
        q_centroid = q_centroid[q_inds]

        # Initializing return values
        rmsd = np.inf

        # Generate all permutation grouped/sorted by the elements
        for p_inds_test in self.permutations(p_atoms):

            p_centroid_test = p_centroid[p_inds_test]
            U_test = self.kabsch(p_centroid_test, q_centroid)

            p_centroid_prime_test = np.dot(p_centroid_test, U_test)
            rmsd_test = np.sqrt(np.mean(np.square(p_centroid_prime_test - q_centroid)))

            if rmsd_test < rmsd:
                p_inds, U, rmsd = p_inds_test, U_test, rmsd_test

        # Rotate and translate matrix P unto matrix Q using Kabsch algorithm.
        # P' = P * U + V
        V = q_trans - np.dot(p_trans, U)

        # Using the original order of the indices
        inds = p_inds[np.argsort(q_inds)]

        return inds, U, V, rmsd

    def fit(self, p: Molecule, ignore_warning=False):
        """Order, rotate and transform `p` molecule according to the best match.

        A `ValueError` will be raised when the total number of possible combinations
        become unfeasible (more than a million combinations).

        Args:
            p: a `Molecule` object what will be matched with the target one.
            ignore_warning: ignoring error when the number of combination is too large

        Returns:
            p_prime: Rotated and translated of the `p` `Molecule` object
            rmsd: Root-mean-square-deviation between `p_prime` and the `target`
        """

        inds, U, V, rmsd = self.match(p, ignore_warning=ignore_warning)

        p_prime = Molecule.from_sites([p[i] for i in inds])
        for site in p_prime:
            site.coords = np.dot(site.coords, U) + V

        return p_prime, rmsd

    @staticmethod
    def permutations(atoms):
        """Generates all the possible permutations of atom order. To achieve better
        performance all tha cases where the atoms are different has been ignored.
        """
        element_iterators = [itertools.permutations(np.where(atoms == element)[0]) for element in np.unique(atoms)]

        for inds in itertools.product(*element_iterators):
            yield np.array(list(itertools.chain(*inds)))


class HungarianOrderMatcher(KabschMatcher):
    """This method pre-aligns the molecules based on their principal inertia
    axis and then re-orders the input atom list using the Hungarian method.

    Notes:
        This method cannot guarantee the best match but is very fast.

        When aligning molecules, the atoms of the two molecules **must** have same number
        of atoms from the same species.
    """

    def match(self, p: Molecule):
        """Similar as `KabschMatcher.match` but this method also finds the order of
        atoms which belongs to the best match.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            inds: The indices of atoms
            U: 3x3 rotation matrix
            V: Translation vector
            rmsd: Root mean squared deviation between P and Q
        """

        if sorted(p.atomic_numbers) != sorted(self.target.atomic_numbers):
            raise ValueError("The number of the same species aren't matching!")

        p_coord, q_coord = p.cart_coords, self.target.cart_coords
        p_atoms, q_atoms = (
            np.array(p.atomic_numbers),
            np.array(self.target.atomic_numbers),
        )

        p_weights = np.array([site.species.weight for site in p])
        q_weights = np.array([site.species.weight for site in self.target])

        # Both sets of coordinates must be translated first, so that
        # their center of mass with the origin of the coordinate system.
        p_trans, q_trans = p.center_of_mass, self.target.center_of_mass
        p_centroid, q_centroid = p_coord - p_trans, q_coord - q_trans

        # Initializing return values
        rmsd = np.inf

        # Generate all permutation grouped/sorted by the elements
        for p_inds_test in self.permutations(p_atoms, p_centroid, p_weights, q_atoms, q_centroid, q_weights):

            p_centroid_test = p_centroid[p_inds_test]
            U_test = self.kabsch(p_centroid_test, q_centroid)

            p_centroid_prime_test = np.dot(p_centroid_test, U_test)
            rmsd_test = np.sqrt(np.mean(np.square(p_centroid_prime_test - q_centroid)))

            if rmsd_test < rmsd:
                inds, U, rmsd = p_inds_test, U_test, rmsd_test

        # Rotate and translate matrix P unto matrix Q using Kabsch algorithm.
        # P' = P * U + V
        V = q_trans - np.dot(p_trans, U)

        return inds, U, V, rmsd

    def fit(self, p: Molecule):
        """Order, rotate and transform `p` molecule according to the best match.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            p_prime: Rotated and translated of the `p` `Molecule` object
            rmsd: Root-mean-square-deviation between `p_prime` and the `target`
        """

        inds, U, V, rmsd = self.match(p)

        # Translate and rotate `mol1` unto `mol2` using Kabsch algorithm.
        p_prime = Molecule.from_sites([p[i] for i in inds])
        for site in p_prime:
            site.coords = np.dot(site.coords, U) + V

        return p_prime, rmsd

    @staticmethod
    def permutations(p_atoms, p_centroid, p_weights, q_atoms, q_centroid, q_weights):
        """Generates two possible permutations of atom order. This method uses the principle component
        of the inertia tensor to prealign the molecules and hungarian method to determine the order.
        There are always two possible permutation depending on the way to pre-aligning the molecules.

        Args:
            p_atoms: atom numbers
            p_centroid: array of atom positions
            p_weights: array of atom weights
            q_atoms: atom numbers
            q_centroid: array of atom positions
            q_weights: array of atom weights

        Yield:
            perm_inds: array of atoms' order
        """
        # get the principal axis of P and Q
        p_axis = HungarianOrderMatcher.get_principal_axis(p_centroid, p_weights)
        q_axis = HungarianOrderMatcher.get_principal_axis(q_centroid, q_weights)

        # rotate Q onto P considering that the axis are parallel and antiparallel
        U = HungarianOrderMatcher.rotation_matrix_vectors(q_axis, p_axis)
        p_centroid_test = np.dot(p_centroid, U)

        # generate full view from q shape to fill in atom view on the fly
        perm_inds = np.zeros(len(p_atoms), dtype=int)

        # Find unique atoms
        species = np.unique(p_atoms)

        for specie in species:
            p_atom_inds = np.where(p_atoms == specie)[0]
            q_atom_inds = np.where(q_atoms == specie)[0]
            A = q_centroid[q_atom_inds]
            B = p_centroid_test[p_atom_inds]

            # Perform Hungarian analysis on distance matrix between atoms of 1st
            # structure and trial structure
            distances = cdist(A, B, "euclidean")
            a_inds, b_inds = linear_sum_assignment(distances)

            perm_inds[q_atom_inds] = p_atom_inds[b_inds]

        yield perm_inds

        # rotate Q onto P considering that the axis are parallel and antiparallel
        U = HungarianOrderMatcher.rotation_matrix_vectors(q_axis, -p_axis)
        p_centroid_test = np.dot(p_centroid, U)

        # generate full view from q shape to fill in atom view on the fly
        perm_inds = np.zeros(len(p_atoms), dtype=int)

        # Find unique atoms
        species = np.unique(p_atoms)

        for specie in species:
            p_atom_inds = np.where(p_atoms == specie)[0]
            q_atom_inds = np.where(q_atoms == specie)[0]
            A = q_centroid[q_atom_inds]
            B = p_centroid_test[p_atom_inds]

            # Perform Hungarian analysis on distance matrix between atoms of 1st
            # structure and trial structure
            distances = cdist(A, B, "euclidean")
            a_inds, b_inds = linear_sum_assignment(distances)

            perm_inds[q_atom_inds] = p_atom_inds[b_inds]

        yield perm_inds

    @staticmethod
    def get_principal_axis(coords, weights):
        """Get the molecule's principal axis.

        Args:
            coords: coordinates of atoms
            weights: the weight use for calculating the inertia tensor

        Returns:
            Array of dim 3 containing the principal axis
        """

        Ixx = Iyy = Izz = Ixy = Ixz = Iyz = 0.0

        for (x, y, z), wt in zip(coords, weights):

            Ixx += wt * (y * y + z * z)
            Iyy += wt * (x * x + z * z)
            Izz += wt * (x * x + y * y)

            Ixy += -wt * x * y
            Ixz += -wt * x * z
            Iyz += -wt * y * z

        inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

        eigvals, eigvecs = np.linalg.eigh(inertia_tensor)

        principal_axis = eigvecs[:, 0]
        return principal_axis

    @staticmethod
    def rotation_matrix_vectors(v1, v2):
        """Returns the rotation matrix that rotates v1 onto v2 using
        Rodrigues' rotation formula.

        See more: https://math.stackexchange.com/a/476311

        Args:
            v1: initial vector
            v2: target vector

        Returns:
            3x3 rotation matrix
        """

        if np.allclose(v1, v2):
            # same direction
            return np.eye(3)

        if np.allclose(v1, -v2):
            # opposite direction: return a rotation of pi around the y-axis
            return np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])

        v = np.cross(v1, v2)
        s = np.linalg.norm(v)
        c = np.vdot(v1, v2)

        vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])

        return np.eye(3) + vx + np.dot(vx, vx) * ((1.0 - c) / (s * s))


class GeneticOrderMatcher(KabschMatcher):
    """This method was inspired by genetic algorithms and tries to match molecules
    based on their already matched fragments.

    It uses the fact that when two molecule is matching their sub-structures have to match as well.
    The main idea here is that in each iteration (generation) we can check the match of all possible
    fragments and ignore those which are not feasible.

    Although in the worst case this method has N! complexity (same as the brute force one),
    in practice it performs much faster because many of the combination can be eliminated
    during the fragment matching.

    Notes:
        This method very robust and returns with all the possible orders.

        There is a well known weakness/corner case: The case when there is
        a outlier with large deviation with a small index might be ignored.
        This happens due to the nature of the average function
        used to calculate the RMSD for the fragments.

        When aligning molecules, the atoms of the two molecules **must** have the
        same number of atoms from the same species.
    """

    def __init__(self, target: Molecule, threshold: float):
        """Constructor of the matcher object.

        Args:
            target: a `Molecule` object used as a target during the alignment
            threshold: value used to match fragments and prune configuration
        """
        super().__init__(target)
        self.threshold = threshold
        self.N = len(target)

    def match(self, p: Molecule):
        """Similar as `KabschMatcher.match` but this method also finds all of the
        possible atomic orders according to the `threshold`.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            Array of the possible matches where the elements are:
                inds: The indices of atoms
                U: 3x3 rotation matrix
                V: Translation vector
                rmsd: Root mean squared deviation between P and Q
        """
        out = []
        for inds in self.permutations(p):
            p_prime = p.copy()
            p_prime._sites = [p_prime[i] for i in inds]

            U, V, rmsd = super().match(p_prime)

            out.append((inds, U, V, rmsd))

        return out

    def fit(self, p: Molecule):
        """Order, rotate and transform all of the matched `p` molecule
        according to the given `threshold`.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            Array of the possible matches where the elements are:
                p_prime: Rotated and translated of the `p` `Molecule` object
                rmsd: Root-mean-square-deviation between `p_prime` and the `target`
        """
        out = []
        for inds in self.permutations(p):
            p_prime = p.copy()
            p_prime._sites = [p_prime[i] for i in inds]

            U, V, rmsd = super().match(p_prime)

            # Rotate and translate matrix `p` onto the target molecule.
            # P' = P * U + V
            for site in p_prime:
                site.coords = np.dot(site.coords, U) + V

            out.append((p_prime, rmsd))

        return out

    def permutations(self, p: Molecule):
        """Generates all of possible permutations of atom order according the threshold.

        Args:
            p: a `Molecule` object what will be matched with the target one.

        Returns:
            Array of index arrays
        """

        # caching atomic numbers and coordinates
        p_atoms, q_atoms = p.atomic_numbers, self.target.atomic_numbers
        p_coords, q_coords = p.cart_coords, self.target.cart_coords

        if sorted(p_atoms) != sorted(q_atoms):
            raise ValueError("The number of the same species aren't matching!")

        # starting maches (only based on element)
        partial_matches = [[j] for j in range(self.N) if p_atoms[j] == q_atoms[0]]

        for i in range(1, self.N):
            # extending the target fragment with then next atom
            f_coords = q_coords[: i + 1]
            f_atom = q_atoms[i]

            f_trans = f_coords.mean(axis=0)
            f_centroid = f_coords - f_trans

            matches = []
            for indices in partial_matches:

                for j in range(self.N):

                    # skipping if the this index is already matched
                    if j in indices:
                        continue

                    # skipping if they are different species
                    if p_atoms[j] != f_atom:
                        continue

                    inds = indices + [j]
                    P = p_coords[inds]

                    # Both sets of coordinates must be translated first, so that
                    # their centroid coincides with the origin of the coordinate system.
                    p_trans = P.mean(axis=0)
                    p_centroid = P - p_trans

                    # The optimal rotation matrix U using Kabsch algorithm
                    U = self.kabsch(p_centroid, f_centroid)

                    p_prime_centroid = np.dot(p_centroid, U)
                    rmsd = np.sqrt(np.mean(np.square(p_prime_centroid - f_centroid)))

                    # rejecting if the deviation is too large
                    if rmsd > self.threshold:
                        continue

                    logger.debug("match - rmsd: {}, inds: {}".format(rmsd, inds))
                    matches.append(inds)

            partial_matches = matches

            logger.info(
                "number of atom in the fragment: {}, number of possible matches: {}".format(i + 1, len(matches))
            )

        return matches
