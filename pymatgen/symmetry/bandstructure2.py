from __future__ import division, unicode_literals

import numpy as np
import networkx as nx
import warnings
import operator
from math import ceil
from math import floor
from math import cos
from math import sin
from math import tan
from math import pi
from math import e
from scipy.optimize import linprog
from scipy.linalg import sqrtm
from warnings import warn
from pymatgen.core.operations import SymmOp, MagSymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bz_labels import LABELS, LABEL_POINTS

"""
Created on July 3, 2018

@author: K Latimer
"""

class HighSymmKpath2(object):
    """
    This class looks for a path along high symmetry lines in the 
    Brillouin zone.
    It is based on the derived symmetry of the energy spectrum
    for a crystalline solid given by A P Cracknell in J. Phys. C:
    Solid State Phys. Vol. 6 (1973), pp. 826-840 ('Van Hove 
    singularities and zero-slope points in crystals') and pp. 841-
    854 ('Van Hove singularities and zero-slope points in magnetic 
    crystals').    
    The user should ensure that the lattice of the input structure
    is as reduced as possible, ie that there is no linear 
    combination of lattice vectors which can produce a vector of
    lesser magnitude than the given set (this is required to
    obtain the correct Brillouin zone within the current 
    implementaiton). This is checked during initialization and a
    warning is issued if the condition is not fulfilled.
    In the case of magnetic structures, care must also be taken to
    provide the magnetic primitive cell (i.e. that which reproduces
    the entire crystal, including the correct magnetic ordering, 
    upon application of lattice translations). There is no way to 
    programatically check for this, so if the input structure is 
    incorrect, the class will output the incorrect kpath without
    any warning being issued.

    Args:
        structure (Structure): Structure object
        symprec (float): Tolerance for symmetry finding
        angle_tolerance (float): Angle tolerance for symmetry finding.
        atol (float): Absolute tolerance used to determine symmetric
        equivalence of points and lines on the BZ.
        has_magmoms (boolean): Whether the input structure contains
            magnetic moments as site properties with the key 'magmom.'
            Values may be in the form of 3-component vectors given in
            the basis of the input lattice vectors, or as scalars, in
            which case the spin axis will default to a_3, the third
            real-space lattice vector (this triggers a warning).
    """
    
    def __init__(self, structure, has_magmoms=False, symprec=0.01, angle_tolerance=5, atol=1e-3):
        self._struct = structure
        self._latt = self._struct.lattice
        
        # Check to see if input lattice is reducible. Ref: B Gruber in Acta. Cryst. Vol. A29,
        # pp. 433-440 ('The Relationship between Reduced Cells in a General Bravais lattice'). 
        # The correct BZ will still be obtained if the lattice vectors are reducible by any
        # linear combination of themselves with coefficients of absolute value less than 2, 
        # hence a missing factor of 2 as compared to the reference.       
        reducible = []        
        for i in range(3):
            for j in range(3):
                if i != j:
                    if np.absolute(np.dot(self._latt.matrix[i], self._latt.matrix[j])) \
                    > np.dot(self._latt.matrix[i], self._latt.matrix[i]) and \
                    np.absolute(np.dot(self._latt.matrix[i], self._latt.matrix[j]) - 
                        np.dot(self._latt.matrix[i], self._latt.matrix[i])) > atol:
                        reducible.append(True)
                    else:
                        reducible.append(False)
        if np.any(reducible):
            print('reducible')
            warnings.warn("The lattice of the input structure is not fully reduced!"
                          "The path can be incorrect. Use at your own risk.")

        self._rec_lattice = self._struct.lattice.reciprocal_lattice
        self._kpath = self._get_kpath(has_magmoms, symprec, angle_tolerance, atol)
        
    @property
    def structure(self):
        """
        Returns:
            The input structure
        """
        return self._struct
    
    @property
    def lattice(self):
        """
        Returns:
            The real space lattice
        """
        return self._latt

    @property
    def rec_lattice(self):
        """
        Returns:
            The reciprocal space lattice
        """
        return self._rec_lattice

    @property
    def kpath(self):
        """
        Returns:
            The symmetry line path in reciprocal space
        """
        return self._kpath

    def get_kpoints(self, line_density=20, coords_are_cartesian=True):
        """
        Returns:
            the kpoints along the paths in cartesian coordinates
            together with the labels for symmetry points -Wei
        """
        list_k_points = []
        sym_point_labels = []
        for b in self.kpath['path']:
            for i in range(1, len(b)):
                start = np.array(self.kpath['kpoints'][b[i - 1]])
                end = np.array(self.kpath['kpoints'][b[i]])
                distance = np.linalg.norm(
                    self._rec_lattice.get_cartesian_coords(start) -
                    self._rec_lattice.get_cartesian_coords(end))
                nb = int(ceil(distance * line_density))
                sym_point_labels.extend([b[i - 1]] + [''] * (nb - 1) + [b[i]])
                list_k_points.extend(
                    [self._rec_lattice.get_cartesian_coords(start)
                     + float(i) / float(nb) *
                     (self._rec_lattice.get_cartesian_coords(end)
                      - self._rec_lattice.get_cartesian_coords(start))
                     for i in range(0, nb + 1)])
        if coords_are_cartesian:
            return list_k_points, sym_point_labels
        else:
            frac_k_points = [self._rec_lattice.get_fractional_coords(k)
                             for k in list_k_points]
            return frac_k_points, sym_point_labels
        
    def _get_kpath(self, has_magmoms, symprec, angle_tolerance, atol):
        decimals = ceil(-1*np.log10(atol)) - 1
        
        ID = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        PAR = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]) #parity, aka the inversion operation (not calling it  
                                                             #INV to avoid confusion with np.linalg.inv() function)        

        ### 1: Get lattices of real and reciprocal structures, and reciprocal point group, and Brillouin zone (BZ) ###
        
        V = self._latt.matrix.T #fractional real space to cartesian real space
        W = self._rec_lattice.matrix.T #fractional reciprocal space to cartesian reciprocal space
        A = np.dot(np.linalg.inv(W), V) #fractional real space to fractional reciprocal space
        
        if has_magmoms:
            grey_struct = self._struct.copy()
            grey_struct.remove_site_property('magmom')
            sga = SpacegroupAnalyzer(grey_struct, symprec=symprec, angle_tolerance=angle_tolerance) 
            grey_ops = sga.get_symmetry_operations()
            struct = self._convert_all_magmoms_to_vectors(self._struct)
            mag_ops = self._get_magnetic_symmetry_operations(self, struct, grey_ops)

            D = [SymmOp.from_rotation_and_translation(rotation_matrix=op.rotation_matrix, 
            translation_vec=op.translation_vector) for op in mag_ops if op.time_reversal == 1]
            
            fD = [SymmOp.from_rotation_and_translation(rotation_matrix=op.rotation_matrix, 
            translation_vec=op.translation_vector) for op in mag_ops if op.time_reversal == -1]

            if np.array([m == np.array([0, 0, 0]) for m in struct.site_properties['magmom']]).all():
                fD = D
                D = []

            if len(fD) == 0: # no operations contain time reversal; type 1
                isomorphic_point_group = [d.rotation_matrix for d in D]
                recip_point_group = self._get_reciprocal_point_group(isomorphic_point_group, ID, A)
            elif len(D) == 0: # all operations contain time reversal / all magmoms zero; type 2
                isomorphic_point_group = [d.rotation_matrix for d in fD]
                recip_point_group =  self._get_reciprocal_point_group(isomorphic_point_group, PAR, A)
            else: # half and half; type 3 or 4
                f = get_coset_factor(D + fD, D)
                isomorphic_point_group = [d.rotation_matrix for d in D]
                recip_point_group =  self._get_reciprocal_point_group(isomorphic_point_group, 
                                                                      np.dot(PAR, f.rotation_matrix), A)

        else:
            sga = SpacegroupAnalyzer(self._struct)
            ops = sga.get_symmetry_operations()
            isomorphic_point_group = [op.rotation_matrix for op in ops]
            recip_point_group =  self._get_reciprocal_point_group(isomorphic_point_group, PAR, A) 
        
        bz = self._rec_lattice.get_wigner_seitz_cell()
        self._recip_point_group = recip_point_group
        
        ### 2: Get all vertices, edge- and face- center points of BZ ("key points") ###
        key_points = []
        face_center_inds = []
        bz_as_key_point_inds = []

        # pymatgen gives BZ in cartesian coordinates; convert to fractional in the primitive basis for reciprocal space
        for (i, facet) in enumerate(bz):
            for (j, vert) in enumerate(facet):
                vert = self._rec_lattice.get_fractional_coords(vert)
                bz[i][j] = vert
        pop = []        
        for i, facet in enumerate(bz):
            rounded_facet = np.around(facet, decimals=decimals)
            u, indices = np.unique(rounded_facet, axis=0, return_index=True)
            if len(u) in [1, 2]:
                pop.append(i)
            else:
                bz[i] = [facet[j] for j in np.sort(indices)]
        bz = [bz[i] for i in range(len(bz)) if i not in pop]
        self._bz = bz
        
        #use vertex points to calculate edge- and face- centers        
        for (i, facet) in enumerate(bz):
            bz_as_key_point_inds.append([])
            for (j, vert) in enumerate(facet):
                edge_center = (vert + facet[j + 1])/2.0 if j != len(facet) - 1 else (vert + facet[0])/2.0
                duplicatevert = False
                duplicateedge = False
                for (k, point) in enumerate(key_points):
                    if np.allclose(vert, point, atol=atol):
                        bz_as_key_point_inds[i].append(k)
                        duplicatevert = True
                        break
                for (k, point) in enumerate(key_points):
                    if np.allclose(edge_center, point, atol=atol):
                        bz_as_key_point_inds[i].append(k)
                        duplicateedge = True
                        break
                if not duplicatevert:
                    key_points.append(vert)
                    bz_as_key_point_inds[i].append(len(key_points) - 1)
                if not duplicateedge:
                    key_points.append(edge_center)
                    bz_as_key_point_inds[i].append(len(key_points) - 1)
            if len(facet) == 4: #parallelogram facet
                face_center = (facet[0] + facet[1] + facet[2] + facet[3])/4.0
                key_points.append(face_center)
                face_center_inds.append(len(key_points) - 1)
                bz_as_key_point_inds[i].append(len(key_points) - 1)
            else: #hexagonal facet
                face_center = (facet[0] + facet[1] + facet[2] + facet[3] + facet[4] + facet[5])/6.0
                key_points.append(face_center)
                face_center_inds.append(len(key_points) - 1)
                bz_as_key_point_inds[i].append(len(key_points) - 1)

        #add gamma point
        key_points.append(np.array([0, 0, 0])) 
        self._key_points = key_points
        ### 3: Find symmetry-equivalent points, which can be mapped to each other by a combination of point group ###
        ### operations and integer translations by lattice vectors. The integers will only be -1, 0, or 1, since  ###
        ### we are restricting to the BZ.                                                                         ###

        key_points_copy = dict(zip(range(len(key_points) - 1), key_points[0:len(key_points) - 1])) 
        # gamma not equivalent to any on BZ and is last point added to key_points
        key_points_inds_orbits = []

        i = 0
        while len(key_points_copy) > 0: 
            key_points_inds_orbits.append([])
            k0ind = list(key_points_copy.keys())[0]
            k0 = key_points_copy[k0ind]
            key_points_inds_orbits[i].append(k0ind)
            key_points_copy.pop(k0ind)

            for op in recip_point_group:
                to_pop = []
                k1 = np.dot(op, k0)
                for ind_key in key_points_copy:
                    diff = k1 - key_points_copy[ind_key]
                    if self._all_ints(diff, atol=atol):
                        key_points_inds_orbits[i].append(ind_key)
                        to_pop.append(ind_key)

                for key in to_pop:
                    key_points_copy.pop(key)
            i += 1

        key_points_inds_orbits.append([len(key_points) - 1])
        self._key_points_inds_orbits = key_points_inds_orbits
        ### 4: Get all lines on BZ between adjacent key points and between gamma and key points ("key lines") ###

        key_lines = []
        gamma_ind = len(key_points) - 1 

        for (i, facet_as_key_point_inds) in enumerate(bz_as_key_point_inds):
            facet_as_key_point_inds_bndy = facet_as_key_point_inds[: len(facet_as_key_point_inds) - 1] 
            # not the face center point (don't need to check it since it's not shared with other facets)
            face_center_ind = facet_as_key_point_inds[-1]
            for (j, ind) in enumerate(facet_as_key_point_inds_bndy):
                if (min(ind, facet_as_key_point_inds_bndy[j - 1]), \
                max(ind, facet_as_key_point_inds_bndy[j - 1])) not in key_lines:
                    key_lines.append((min(ind, facet_as_key_point_inds_bndy[j - 1]), 
                                      max(ind, facet_as_key_point_inds_bndy[j - 1])))
                k = j + 1 if j != len(facet_as_key_point_inds_bndy) - 1 else 0
                if (min(ind, facet_as_key_point_inds_bndy[k]), max(ind, facet_as_key_point_inds_bndy[k])) not in key_lines:
                    key_lines.append((min(ind, facet_as_key_point_inds_bndy[k]), max(ind, facet_as_key_point_inds_bndy[k])))
                if (ind, gamma_ind) not in key_lines:
                    key_lines.append((ind, gamma_ind))
                key_lines.append((min(ind, face_center_ind), max(ind, face_center_ind)))
            key_lines.append((face_center_ind, gamma_ind)) 

        ### 5: Find symmetry-equivalent key lines, defined as endpoints of first line being equivalent ###
        ### to end points of second line, and a random point in between being equivalent to the mapped ###
        ### random point.                                                                              ###

        key_lines_copy = dict(zip(range(len(key_lines)), key_lines))
        key_lines_inds_orbits = []

        i = 0
        while len(key_lines_copy) > 0:
            key_lines_inds_orbits.append([])
            l0ind = list(key_lines_copy.keys())[0]
            l0 = key_lines_copy[l0ind]
            key_lines_inds_orbits[i].append(l0)
            key_lines_copy.pop(l0ind)
            to_pop = []
            p00 = key_points[l0[0]]
            p01 = key_points[l0[1]]
            pmid0 = p00 + e/pi*(p01 - p00)
            for ind_key in key_lines_copy:

                l1 = key_lines_copy[ind_key]
                p10 = key_points[l1[0]]
                p11 = key_points[l1[1]]                       
                equivptspar = False
                equivptsperp = False
                equivline = False

                if np.array([l0[0] in orbit and l1[0] in orbit for orbit in key_points_inds_orbits]).any() \
                and np.array([l0[1] in orbit and l1[1] in orbit for orbit in key_points_inds_orbits]).any():
                    equivptspar = True
                elif np.array([l0[1] in orbit and l1[0] in orbit for orbit in key_points_inds_orbits]).any() \
                and np.array([l0[0] in orbit and l1[1] in orbit for orbit in key_points_inds_orbits]).any():
                    equivptsperp = True  

                if equivptspar:
                    pmid1 = p10 + e/pi*(p11 - p10)
                    for op in recip_point_group:
                        if not equivline:
                            p00pr = np.dot(op, p00)
                            diff0 = p10 - p00pr
                            if self._all_ints(diff0, atol=atol):                                     
                                pmid0pr = np.dot(op, pmid0) + diff0
                                p01pr = np.dot(op, p01) + diff0
                                if np.allclose(p11, p01pr, atol=atol) and np.allclose(pmid1, pmid0pr, atol=atol):
                                    equivline = True

                elif equivptsperp:
                    pmid1 = p11 + e/pi*(p10 - p11)
                    for op in recip_point_group:
                        if not equivline:
                            p00pr = np.dot(op, p00)
                            diff0 = p11 - p00pr
                            if self._all_ints(diff0, atol=atol):                               
                                pmid0pr = np.dot(op, pmid0) + diff0
                                p01pr = np.dot(op, p01) + diff0
                                if np.allclose(p10, p01pr, atol=atol) and np.allclose(pmid1, pmid0pr, atol=atol):
                                    equivline = True

                if equivline:
                    key_lines_inds_orbits[i].append(l1)
                    to_pop.append(ind_key)

            for key in to_pop:
                key_lines_copy.pop(key)
            i += 1
        self._key_lines_inds_orbits = key_lines_inds_orbits
        ### 6: Get little groups for key points (group of symmetry elements present at that point) ###

        little_groups_points = [] # elements are lists of indicies of recip_point_group. the 
                                  # list little_groups_points[i] is the little group for the 
                                  # orbit key_points_inds_orbits[i]
        v1 = W
        for (i, orbit) in enumerate(key_points_inds_orbits):
            k0 = key_points[orbit[0]]
            little_groups_points.append([])
            for (j, op) in enumerate(recip_point_group):
                gamma_to = np.dot(op, -1*k0) + k0     
                check_gamma = True
                if not self._all_ints(gamma_to, atol=atol):
                    check_gamma = False
                if check_gamma:
                    little_groups_points[i].append(j)

        ### 7: Get little groups for key lines (group of symmetry elements present at every point ###
        ### along the line). This is implemented by testing the symmetry at a point e/pi of the   ###
        ### way between the two endpoints.                                                        ###

        little_groups_lines = [] #elements are lists of indicies of recip_point_group. the list little_groups_lines[i] is 
                                 #the little group for the orbit key_points_inds_lines[i]

        for (i, orbit) in enumerate(key_lines_inds_orbits):
            l0 = orbit[0]
            v = key_points[l0[1]] - key_points[l0[0]]
            k0 = key_points[l0[0]] + np.e/pi*v
            little_groups_lines.append([])
            for (j, op) in enumerate(recip_point_group):
                gamma_to = np.dot(op, -1*k0) + k0
                check_gamma = True
                if not self._all_ints(gamma_to, atol=atol):
                    check_gamma = False
                if check_gamma:
                    little_groups_lines[i].append(j)

        ### 8: Choose key lines for k-path. Current implementation: include all lines which are two-fold symmetry ###
        ### axes (includes 4-fold and 6-fold), or three-fold symmetry axes and the line of intersection of three  ###
        ### mirror planes (C3v/C6v pt gp) ###

        line_orbits_in_path = []
        point_orbits_in_path = []
        for (i, little_group) in enumerate(little_groups_lines):
            add_rep = False
            nC2 = 0
            nC3 = 0
            nsig = 0
            for j, opind in enumerate(little_group):
                op = recip_point_group[opind]
                if not (op == ID).all():
                    if (np.dot(op, op) == ID).all():
                        if np.linalg.det(op) == 1:
                            nC2 += 1
                        elif not (op == PAR).all():
                            nsig += 1
                    elif (np.dot(op, np.dot(op, op)) == ID).all() and np.linalg.det(op) == 1:
                        nC3 += 1
            if nC2 > 0 or (nC3 > 0 and nsig > 2):
                add_rep = True

            if add_rep:
                line_orbits_in_path.append(i)
                l = key_lines_inds_orbits[i][0]
                ind0 = l[0]
                ind1 = l[1]
                found0 = False
                found1 = False
                for (j, orbit) in enumerate(key_points_inds_orbits):
                    if ind0 in orbit:
                        point_orbits_in_path.append(j)
                        found0 = True
                    if ind1 in orbit:
                        point_orbits_in_path.append(j)
                        found1 = True
                    if found0 and found1:
                        break
        point_orbits_in_path = list(set(point_orbits_in_path))

        ### 9: Choose remaining unconnected key points for k-path. Under SC criteria, these would be points which ###
        ### either contain inversion symmetry, or which have a rotation/rotoinversion axis passing through the     ###
        ### point but not on the BZ boundary or going through gamma (not sure if this occurs so I'm adding a flag  ###
        ### for it). Connect them to gamma if selected (this is arbitrary but there needs to be some connection    ###
        ### for a path to exist).                                                                                  ###

        unconnected = []

        for i in range(len(key_points_inds_orbits)):
            if i not in point_orbits_in_path:
                unconnected.append(i)

        for ind in unconnected:
            connect = False
            for op_ind in little_groups_points[ind]:
                op = recip_point_group[op_ind]
                if (op == ID).all():
                    pass
                elif (op == PAR).all():
                    connect = True
                    break
                elif np.linalg.det(op) == 1:
                    if (np.dot(op, np.dot(op, op)) == ID).all():
                        pass
                    else:
                        connect = True
                        break
                else:
                    pass
            if connect:
                l = (key_points_inds_orbits[ind][0], gamma_ind)
                for (j, orbit) in enumerate(key_lines_inds_orbits):
                    if l in orbit:
                        line_orbits_in_path.append(j)
                        break
                if gamma_ind not in point_orbits_in_path:
                    point_orbits_in_path.append(gamma_ind)
                point_orbits_in_path.append(ind)
                
        ### 10: Consolidate selected segments into a single irreducible section of the Brilouin zone (as determined ###
        ### by the reciprocal point and lattice symmetries). This is accomplished by identifying the boundary       ###
        ### planes of the IRBZ.`                                                                                    ###
        # return None
        IRBZ_points_inds = self._get_IRBZ(recip_point_group, W, key_points, face_center_inds, atol)
        lines_in_path_inds = []
        for ind in line_orbits_in_path:
            for tup in key_lines_inds_orbits[ind]:
                if tup[0] in IRBZ_points_inds and tup[1] in IRBZ_points_inds:
                    lines_in_path_inds.append(tup)
                    break
        G = nx.Graph(lines_in_path_inds)
        lines_in_path_inds = list(nx.edge_dfs(G))
        points_in_path_inds = [ind for tup in lines_in_path_inds for ind in tup]
        points_in_path_inds_unique = list(set(points_in_path_inds))

        orbit_labels = self._get_orbit_labels(key_points_inds_orbits, key_points, atol)
        key_points_labels = ['' for i in range(len(key_points))]
        for i, orbit in enumerate(key_points_inds_orbits):
            for point_ind in orbit:
                key_points_labels[point_ind] = LABELS[int(orbit_labels[i])]

        kpoints = {}
        reverse_kpoints = {}
        for point_ind in points_in_path_inds_unique:
            point_label = key_points_labels[point_ind]
            if point_label not in kpoints.keys():
                kpoints[point_label] = key_points[point_ind]
                reverse_kpoints[point_ind] = point_label
            else:
                existing_labels = [key for key in kpoints.keys() if point_label in key]
                if '\'' not in point_label:
                    existing_labels[:] = [label for label in existing_labels if '\'' not in label]
                if len(existing_labels) == 1:
                    max_occurence = 0
                else:
                    max_occurence = max([int(label[3:-1]) for label in existing_labels[1:]])
                kpoints[point_label + '_{' + str(max_occurence + 1) + '}'] = key_points[point_ind]
                reverse_kpoints[point_ind] = point_label + '_{' + str(max_occurence + 1) + '}'
        print('hello world')
        path = []
        i = 0
        start_of_subpath = True
        while i < len(points_in_path_inds):
            if start_of_subpath:
                path.append([reverse_kpoints[points_in_path_inds[i]]])
                i += 1
                start_of_subpath = False
            elif points_in_path_inds[i] == points_in_path_inds[i+1]:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 2
            else:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 1
                start_of_subpath = True
            if i == len(points_in_path_inds) - 1:
                path[-1].append(reverse_kpoints[points_in_path_inds[i]])
                i += 1

        return {'kpoints': kpoints, 'path': path}

    def _convert_all_magmoms_to_vectors(self, struct):
        if not 'magmom' in struct.site_properties:
            warnings.warn('The \'magmom\' property is not set in the structure\'s site properties.' 
                  'All magnetic moments are being set to zero.')
            struct.add_site_property('magmom', [np.array([0, 0, 0]) for i in range(len(struct.sites))])
            return struct

        old_magmoms = struct.site_properties['magmom']
        new_magmoms = []
        found_scalar = False

        for magmom in old_magmoms:
            if type(magmom) == np.ndarray:
                new_magmoms.append(magmom)
            elif type(magmom) == list:
                new_magmoms.append(np.array(magmom))
            else:
                found_scalar = True
                new_magmoms.append(np.array([0, 0, magmom]))

        if found_scalar:
            warnings.warn('At least one magmom had a scalar value. Defaulted to z+ spinor.')

        struct.remove_site_property('magmom')
        struct.add_site_property('magmom', new_magmoms)
        return struct
        
    def _get_magnetic_symmetry_operations(self, struct, grey_ops):
        mag_ops = []
        magmoms = struct.site_properties['magmom']
        nonzero_magmom_inds = [i for i in range(len(struct.sites)) if not (magmoms[i] == np.array([0, 0, 0])).all()]
        init_magmoms = [site.properties['magmom'] for (i, site) in enumerate(struct.sites) if i in nonzero_magmom_inds] 
        sites = [site for (i, site) in enumerate(struct.sites) if i in nonzero_magmom_inds]
        init_site_coords = [site.frac_coords for site in sites]
        for op in grey_ops:
            r = op.rotation_matrix
            t = op.translation_vector
            xformed_magmoms = [self._apply_op_to_magmom(self, r, magmom) for magmom in init_magmoms] 
            xformed_site_coords = [np.dot(r, site.frac_coords) + t for site in sites]
            permutation = ['a' for i in range(len(sites))]
            not_found = list(range(len(sites)))
            for i in range(len(sites)):
                xformed = xformed_site_coords[i]
                for k, j in enumerate(not_found):
                    init = init_site_coords[j]
                    diff = xformed - init
                    if self._all_ints(diff, atol=atol):
                        permutation[i] = j
                        not_found.pop(k)
                        break

            same = np.zeros(len(sites))
            flipped = np.zeros(len(sites))
            for i, magmom in enumerate(xformed_magmoms):
                if (magmom == init_magmoms[permutation[i]]).all():
                    same[i] = 1
                elif (magmom == -1*init_magmoms[permutation[i]]).all():
                    flipped[i] = 1

            if same.all(): #add symm op without tr
                mag_ops.append(MagSymmOp.from_rotation_and_translation_and_time_reversal(rotation_matrix=op.rotation_matrix,
                                                                                        translation_vec=op.translation_vector,
                                                                                        time_reversal=1))
            if flipped.all(): #add symm op with tr
                mag_ops.append(MagSymmOp.from_rotation_and_translation_and_time_reversal(rotation_matrix=op.rotation_matrix,
                                                                                        translation_vec=op.translation_vector,
                                                                                        time_reversal=-1))

        return mag_ops
        
    def _get_reciprocal_point_group(self, ops, R, A):
        Ainv = np.linalg.inv(A)
        recip_point_group = [R]
        for op in ops:
            op = np.around(np.dot(A, np.dot(op, Ainv)), decimals=2) #convert to reciprocal primitive basis
            new = True
            new_coset = True
            for thing in recip_point_group:
                if (thing == op).all():
                    new = False
                if (thing == np.dot(R, op)).all():
                    new_coset = False

            if new:
                recip_point_group.append(op)
            if new_coset:
                recip_point_group.append(np.dot(R, op))

        return recip_point_group
        
    def _get_coset_factor(self, G, H):
         # finds g for left coset decomposition G = H + gH (H must be subgroup of G with index two.) 
         # in this implementation, G and H are lists of objects of type SymmOp
        gH = []
        for i, op1 in enumerate(G):
            in_H = False
            for op2 in H:
                if op1.equiv(op2): # TODO: equiv fn in symmop to account for integer lattice translations
                    in_H = True
                    break
            if not in_H:
                gH.append(op1)

        for op in gH:
            opH = [op.__mul__(h) for h in H]
            is_coset_factor = True
            for op1 in opH:
                for op2 in H:
                    if op1.equiv(op2):
                        is_coset_factor = False
                        break
                if not is_coset_factor:
                    break
            if is_coset_factor:
                return op

        return "No coset factor found."
        
    def _apply_op_to_magmom(self, r, magmom):
        if np.linalg.det(r) == 1:
            return np.dot(r, magmom)
        else:
            return -1*np.dot(r, magmom)
        
    def _all_ints(self, arr, atol):
        rounded_arr = np.around(arr, decimals=0)
        return np.allclose(rounded_arr, arr, atol=atol)
    
    def _get_IRBZ(self, recip_point_group, W, key_points, face_center_inds, atol):
        PAR = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        rpgdict = self._get_reciprocal_point_group_dict(recip_point_group, atol)
        self._rpgdict = rpgdict

        g = np.dot(W.T, W) # just using change of basis matrix rather than 
                           # Lattice.get_cartesian_coordinates for conciseness
        ginv = np.linalg.inv(g)
        D = np.linalg.det(W)

        primary_orientation = None
        secondary_orientation = None
        tertiary_orientation = None

        planar_boundaries = []
        IRBZ_points = [(i, point) for i, point in enumerate(key_points)]

        for sigma in rpgdict['reflections']:
            norm = sigma['normal']
            if primary_orientation is None:
                primary_orientation = norm
                planar_boundaries.append(norm)
            elif np.isclose(np.dot(primary_orientation, np.dot(g, norm)), 0, atol=atol):
                if secondary_orientation is None:
                    secondary_orientation = norm
                    planar_boundaries.append(norm)
                elif np.isclose(np.dot(secondary_orientation, np.dot(g, norm)), 0, atol=atol):
                    if tertiary_orientation is None:
                        tertiary_orientation = norm
                        planar_boundaries.append(norm)
                    elif np.allclose(norm, -1*tertiary_orientation, atol=atol):
                        pass
                elif np.dot(secondary_orientation, np.dot(g, norm)) < 0:
                    planar_boundaries.append(-1*norm)
                else:
                    planar_boundaries.append(norm)
            elif np.dot(primary_orientation, np.dot(g, norm)) < 0:
                planar_boundaries.append(-1*norm)
            else:
                planar_boundaries.append(norm)
        
        IRBZ_points = self._reduce_IRBZ(IRBZ_points, planar_boundaries, g, atol)

        used_axes = []

        for rotn in rpgdict['rotations']['six-fold']: #six-fold rotoinversion always comes with horizontal mirror so don't need to check
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict['rotations']['rotoinv-four-fold']:
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict['rotations']['four-fold']:
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict['rotations']['rotoinv-three-fold']:
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(sqrtm(-1*op), cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict['rotations']['three-fold']:
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        for rotn in rpgdict['rotations']['two-fold']:
            ax = rotn['axis']
            op = rotn['op']
            if not np.any([np.allclose(ax, usedax, atol) for usedax in used_axes]):
                if self._op_maps_IRBZ_to_self(op, IRBZ_points, atol):
                    face_center_found = False
                    for point in IRBZ_points:
                        if point[0] in face_center_inds:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                face_center_found = True
                                used_axes.append(ax)
                                break
                    if not face_center_found:
                        print('face center not found')
                        for point in IRBZ_points:
                            cross = D*np.dot(ginv, np.cross(ax, point[1]))
                            if not np.allclose(cross, 0, atol=atol):
                                rot_boundaries = [cross, -1*np.dot(op, cross)]
                                used_axes.append(ax)
                                break
                    IRBZ_points = self._reduce_IRBZ(IRBZ_points, rot_boundaries, g, atol)

        return [point[0] for point in IRBZ_points]
    
    def _get_reciprocal_point_group_dict(self, recip_point_group, atol):        
        ID = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        PAR = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        
        d = {'reflections': [],
             'rotations': {'two-fold': [], 'three-fold': [], 'four-fold': [], 'six-fold': [], 
             'rotoinv-three-fold': [], 'rotoinv-four-fold': [], 'rotoinv-six-fold': []},
             'inversion': []}

        for i, op in enumerate(recip_point_group):
            evals, evects = np.linalg.eig(op)
            if np.all(op == ID):
                continue
            elif np.all(op == PAR):
                d['inversion'].append({'ind': i, 'op': PAR})
            elif np.all(np.dot(op, op) == ID):
                if np.linalg.det(op) == 1: #two-fold rotation
                    for j in range(3):
                        if np.isclose(evals[j], 1, atol=atol):
                            ax = evects[:, j]
                    d['rotations']['two-fold'].append({'ind': i, 'axis': ax, 'op': op})
                else: #reflections
                    for j in range(3):
                        if np.isclose(evals[j], -1, atol=atol):
                            norm = evects[:, j]
                    d['reflections'].append({'ind': i, 'normal': norm, 'op': op})
            elif np.all(np.dot(op, np.dot(op, op)) == ID): #three-fold rotation
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['three-fold'].append({'ind': i, 'axis': ax, 'op': op})
            elif  np.all(np.dot(op, np.dot(op, op)) == PAR): #three-fold rotoinversion
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['rotoinv-three-fold'].append({'ind': i, 'axis': ax, 'op': op})
            elif np.all(np.dot(op, np.dot(op, np.dot(op, op))) == ID) and np.linalg.det(op) == 1: #four-fold rotation
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['four-fold'].append({'ind': i, 'axis': ax, 'op': op})
            elif np.all(np.dot(op, np.dot(op, np.dot(op, op))) == ID) and np.linalg.det(op) == -1: #four-fold rotoinversion
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['rotoinv-four-fold'].append({'ind': i, 'axis': ax, 'op': op})
            elif np.all(np.dot(op, np.dot(op, np.dot(op, np.dot(op, np.dot(op, op))))) == ID) and np.linalg.det(op) == 1: #six-fold rotation
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['six-fold'].append({'ind': i, 'axis': ax, 'op': op})
            elif np.all(np.dot(op, np.dot(op, np.dot(op, np.dot(op, np.dot(op, op))))) == ID) and np.linalg.det(op) == -1: #six-fold rotoinversion
                for j in range(3):
                    if np.isreal(evals[j]) and np.isclose(np.absolute(evals[j]), 1, atol=atol):
                        ax = evects[:, j]
                d['rotations']['rotoinv-six-fold'].append({'ind': i, 'axis': ax, 'op': op})

        return d

    def _op_maps_IRBZ_to_self(self, op, IRBZ_points, atol):
        point_coords = [point[1] for point in IRBZ_points]
        for point in point_coords:
            point_prime = np.dot(op, point)
            mapped_back = False
            for checkpoint in point_coords:
                if np.allclose(point_prime, checkpoint, atol):
                    mapped_back = True
                    break
            if not mapped_back:
                return False

        return True


    def _reduce_IRBZ(self, IRBZ_points, boundaries, g, atol):
        in_reduced_section = []
        for point in IRBZ_points:
            in_reduced_section.append(np.all([(np.dot(point[1], np.dot(g, boundary)) >= 0 or 
                np.isclose(np.dot(point[1], np.dot(g, boundary)), 0, atol=atol)) for boundary in boundaries]))

        return [IRBZ_points[i] for i in range(len(IRBZ_points)) if in_reduced_section[i]]

    def _get_orbit_labels(self, key_points_inds_orbits, key_points, atol):
        orbit_cosines = []
        for i, orbit in enumerate(key_points_inds_orbits):
            orbit_cosines.append(sorted(sorted([(j, np.round(np.dot(key_points[k], LABEL_POINTS[j])/ \
            (np.linalg.norm(key_points[k])*np.linalg.norm(LABEL_POINTS[j])), decimals=3)) \
            for k in orbit for j in range(26)], key=operator.itemgetter(0)), key=operator.itemgetter(1), reverse=True))

        orbit_labels_unsorted = [(len(key_points_inds_orbits) - 1, 26)]
        pop_orbits = []
        pop_labels = []

        for i, orb_cos in enumerate(orbit_cosines):
            if np.isclose(orb_cos[0][1], 1.0, atol=atol):
                orbit_labels_unsorted.append((i, orb_cos[0][0])) #(point orbit index, label index)
                pop_orbits.append(i)
                pop_labels.append(orb_cos[0][0])

        orbit_cosines = self._reduce_cosines_array(orbit_cosines, pop_orbits, pop_labels)

        while len(orbit_labels_unsorted) < len(key_points_inds_orbits):
            pop_orbits = []
            pop_labels = []
            max_cosine_value = max([orb_cos[0][1] for orb_cos in orbit_cosines])
            max_cosine_value_inds = [j for j in range(len(orbit_cosines)) if orbit_cosines[j][0][1] == max_cosine_value]
            max_cosine_label_inds = self._get_max_cosine_labels([orbit_cosines[j] for j in max_cosine_value_inds])

            for j, label_ind in enumerate(max_cosine_label_inds):
                orbit_labels_unsorted.append((max_cosine_value_inds[j], label_ind))
                pop_orbits.append(max_cosine_value_inds[j])
                pop_labels.append(label_ind)
            orbit_cosines = self._reduce_cosines_array(orbit_cosines, pop_orbits, pop_labels)

        orbit_labels = np.zeros(len(key_points_inds_orbits))
        for tup in orbit_labels_unsorted:
            orbit_labels[tup[0]] = tup[1]

        return orbit_labels

    def _reduce_cosines_array(self, orbit_cosines, pop_orbits, pop_labels):

        return [[orb_cos[i] for i in range(len(orb_cos)) if orb_cos[i][0] not in pop_labels] \
        for j, orb_cos in enumerate(orbit_cosines) if j not in pop_labels]

    def _get_max_cosine_labels(self, max_cosine_orbits):
        max_cosine_label_inds = np.zeros(len(max_cosine_orbits))
        initial_max_cosine_label_inds = [max_cos_orb[0][0] for max_cos_orb in max_cosine_orbits]
        u, inds, counts = np.unique(initial_max_cosine_label_inds, return_index=True, return_counts=True)

        for i, ind in enumerate(inds):
            if counts[i] == 1:
                max_cosine_label_inds[ind] = initial_max_cosine_label_inds[ind]
            else:
                recurs_max_cosine_label_inds = self._get_max_cosine_labels([max_cos_orb[1:] for \
                max_cos_orb in max_cosine_orbits if max_cos_orb[0][0] == initial_max_cosine_label_inds[ind]])
                k = 0
                for j, max_cos_orb in enumerate(max_cosine_orbits):
                    if max_cos_orb[0][0] == initial_max_cosine_label_inds[ind]:
                        max_cosine_label_inds[j] = recurs_max_cosine_label_inds[k]
                        k += 1

        return max_cosine_label_inds