__author__ = 'HongDing'

from pymatgen.core.surface import get_symmetrically_distinct_miller_indices as get_indices
from pymatgen.core.surface import SlabGenerator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.matproj.rest import MPRester
from monty.serialization import dumpfn,loadfn
from Elastic_Def_Energy_Class import DefGradientMatrix,DefEnergyDensity_from_file
from fractions import gcd
from MP_API import my_MP_API_Key
from math import pi
import numpy as np
import time
import json


def structure_interface_matching(str_film,str_sub,criteria):
    """
    Calculate the reduced matching lattice vectors for heterostructure interfaces
    The method is proposed in by Zur and McGill:
    Journal of Applied Physics 55 (1984), 378 ; doi: 10.1063/1.333084

    :param str_film: film structure object
    :param str_sub: substrate structure object
    :param criteria: A dict of mathcing creteria
    :return all_matching_results:
                a list of matching result tuple and each tuple contains:
                    [0] Film Miller index
                    [1] Substrate Miller index
                    [2] Matching vector sets (tuple)
                        (0) film vector 0
                        (1) index of vector 0
                        (2) mismatch strain of vecotr 0
                        (3) film vector 1
                        (4) index of vector 1
                        (5) mismatch strain fo vecotr 1
                    [3] Matching area general info (tuple)
                        (0) Match area
                        (1) lengh of super lattice vector 0 for film
                        (2) lengh of super lattice vector 1 for film
                        (3) angle between super lattice vecotrs 0 and 1 for film
                        (4) lengh of super lattice vector 0 for substrate
                        (5) lengh of super lattice vector 1 for substrate
                        (6) angle between super lattice vecotrs 0 and 1 for substrate



    """

    def vec_Zur_def(a,b):
        """
        :param a: one two-dimensional vector
        :param b: another two-dimensional vector
        :return: a order of vectors of (a, b) that satisfy the Zur vector definition
                 for the surface vector pairs
        """
        if np.dot(a, b) < 0:
            b = -b
            return vec_Zur_def(a, b)

        if np.linalg.norm(a) > np.linalg.norm(b):
            a, b = b, a
            return vec_Zur_def(a, b)

        if np.linalg.norm(b)> np.linalg.norm(np.add(b, a)):
            b = np.add(b, a)
            return vec_Zur_def(a, b)

        if np.linalg.norm(b)> np.linalg.norm(np.subtract(b, a)):
            b = np.subtract(b, a)
            return vec_Zur_def(a, b)

        return [a,b]

    def surf_area(vec):
        """
        :param vec: a list consist of two vectors
        :return: The area of the parallelogram that the two vectors make
        """
        return np.linalg.norm(np.cross(vec[0],vec[1]))

    def factor(n):
        """
        :param n: integer number
        :return:  a list of factors of integer n
        """
        return [x for x in range(1, n+1) if n % x == 0]

    def gen_all_matrices(area_num):
        """
        :param area_num: Area input
        :return: all possible matrices of
                | i ,  j      |
                | 0 ,  area/i |
                equal to Area
        """
        matrices_list = list()
        for i in factor(area_num):
            for j in range(area_num/i):
                matrix = np.matrix(((i, j),(0, area_num/i)))
                matrices_list.append(matrix)

        return matrices_list
    def py_ang(v1, v2):
        """
        :param v1: one vector
        :param v2: another vector
        :return: the angle between the two factors
        """
        cosang = np.dot(v1, v2)
        sinang = np.linalg.norm(np.cross(v1, v2))
        return np.arctan2(sinang, cosang)

    def rel_len_strain(vec_comp1, vec_comp2):
        """
        :param vec_comp1: One lattice vector
        :param vec_comp2: Another lattice vector
        :return: relative length difference
        """
        return np.linalg.norm(vec_comp2)/np.linalg.norm(vec_comp1)-1

    def ref_angle(vec_set1, vec_set2):
        """
        :param vec_comp1: One lattice vector set (consist a list of two vectors)
        :param vec_comp2: Another lattice vector set (consist a list of two vectors)
        :return: relative angle difference between the two set
        """
        return py_ang(vec_set2[0],vec_set2[1])/py_ang(vec_set1[0],vec_set1[1])-1

    def Tol(vec_set1, vec_set2,tol_value_length,tol_value_angle):
        """
        :param vec_comp1: One lattice vector set (consist a list of two vectors)
        :param vec_comp2: Another lattice vector set (consist a list of two vectors)
        :param tol_value_length: lattice vector length difference tolerance
        :param tol_value_angle: angle difference tolerance
        :return: Whether the two vector sets satisify length and angle tolerance
        """
        return np.absolute(rel_len_strain(vec_set1[0],vec_set2[0]))< tol_value_length\
            and np.absolute(rel_len_strain(vec_set1[1],vec_set2[1]))< tol_value_length\
            and np.absolute(ref_angle(vec_set1,vec_set2)) < tol_value_angle


    # Specify film miller index under consideration
    if not criteria["film_miller"] == None:
        film_miller_index_list = criteria["film_miller"]
    else:
        film_miller_index_list = sorted(get_indices(str_film,criteria["film_max_miller_index"]))

    # Specify substrate miller index under consideration
    if not criteria["sub_miller"] == None:
        sub_miller_index_list = criteria["sub_miller"]
    else:
        sub_miller_index_list = sorted(get_indices(str_sub,criteria["sub_max_miller_index"]))

    # Final output data
    all_matching_results = list()


    # Loop through miller indices between substrate and film
    for sub_miller in sub_miller_index_list:
        for film_miller in film_miller_index_list:

            # Generate film and substrate slabs
            film_slab = SlabGenerator(str_film, film_miller, 20, 15,primitive=False).get_slab()
            sub_slab  = SlabGenerator(str_sub, sub_miller, 20, 15,primitive=False).get_slab()

            # Obtain the slab vector info.
            film_vectors, sub_vectors = film_slab.lattice_vectors(), sub_slab.lattice_vectors()

            # Transform the slab vectors to satifiy Zur vector definition
            film_surf_vectors = vec_Zur_def(film_vectors[0],film_vectors[1])
            sub_surf_vectors = vec_Zur_def(sub_vectors[0], sub_vectors[1])

            # Areas of surface for film and substrate
            film_area, sub_area = surf_area(film_surf_vectors), surf_area(sub_surf_vectors)

            sub_miller_result = list()
            reduced_vector_result = list()

            # Loop through possible area for matching
            for i in range(1, int(criteria["max_area"]/sub_area)):
                for j in range(1, int(criteria["max_area"]/film_area)): # j*film_area = i*sub_area
                    # i and j should be coprime integers   & Area difference are within max_area_ratio_tol
                    if gcd(i,j) == 1 and np.absolute(film_area/sub_area  - float(i)/j) < criteria["max_area_ratio_tol"]:

                        # Consider all combination matrices
                        film_matrices_list, sub_matrices_list = gen_all_matrices(j), gen_all_matrices(i)

                        # Loop through combination matrices
                        for sub_matrix in sub_matrices_list:

                            # Searching the smallest mathcing area
                            min_area_recorder = criteria["max_area"]

                            for film_matrix in film_matrices_list:

                                # Simply the super lattice vectors
                                film_SL_vec_unZur   = np.squeeze(np.asarray(film_matrix*((film_surf_vectors[0]),\
                                                                                         (film_surf_vectors[1]))))
                                sub_SL_vec_unZur    = np.squeeze(np.asarray(sub_matrix*((sub_surf_vectors[0]),\
                                                                                        (sub_surf_vectors[1]))))

                                # Make the super lattice vector satisfy Zur vector definition
                                film_SL_vec = vec_Zur_def(film_SL_vec_unZur[0],film_SL_vec_unZur[1])
                                sub_SL_vec  = vec_Zur_def(sub_SL_vec_unZur[0],sub_SL_vec_unZur[1])


                                # Find smaller super lattice mathing and satify deformation tolerance
                                if Tol(film_SL_vec,sub_SL_vec,criteria["max_length_tol"],criteria["max_angle_tol"]) \
                                        and j*film_area<min_area_recorder:

                                    # Pack matching area general info
                                    new_entry = np.array([  j*film_area,\
                                                            np.linalg.norm(film_SL_vec[0]), \
                                                            np.linalg.norm(film_SL_vec[1]), \
                                                            py_ang(film_SL_vec[0],film_SL_vec[1])/pi*180,\
                                                            np.linalg.norm(sub_SL_vec[0]), \
                                                            np.linalg.norm(sub_SL_vec[1]),\
                                                            py_ang(sub_SL_vec[0],sub_SL_vec[1])/pi*180 ]
                                                        )

                                    similar_find = False

                                    # Check whether it is a duplication in the result
                                    if not similar_find:
                                        min_area_recorder = j*film_area
                                        film_vec_reduced_cell, sub_vec_reduced_cell = {}, {}

                                        # Solve the vector index
                                        for k in range(2):
                                            film_vec_reduced_cell[k] = np.linalg.solve(str_film.lattice.matrix, np.array(film_SL_vec[k]))
                                            sub_vec_reduced_cell[k] = np.linalg.solve(str_sub.lattice.matrix, np.array(sub_SL_vec[k]))


                                        sub_miller_result.append(new_entry)

                                        # Pack matching vector info.
                                        reduced_vector_result.append(
                                            tuple([film_SL_vec[0],\
                                                 ([int(round(k)) for k in sub_vec_reduced_cell[0]]),\
                                                 (rel_len_strain(film_SL_vec[0],sub_SL_vec[0])),\
                                                  film_SL_vec[1],\
                                                 ([int(round(k)) for k in sub_vec_reduced_cell[1]]),\
                                                 (rel_len_strain(film_SL_vec[1],sub_SL_vec[1]))\
                                               ]\
                                            )
                                        )

            # Reorganize the final output
            for element, vectors in zip(sub_miller_result,reduced_vector_result):
                all_matching_results.append((film_miller,sub_miller,vectors,element))

    return all_matching_results

def str_prepartion(Source='MP',mpid_list=[],final=True,str_file_list_by_user= []):
    """
    :param Source: if "MP", we will obtain the structure object from Materials Project
                    otherwise, we will read from user_input_str_file_list
    :param mpid_list: mp_id list to query
    :param final: (bool) whether choose the relaxed  (True) or initial (False) structure
    :param user_input_str_file_list: user selected structure file
    :return: a list of structure objects   and    a list of structure id
    """
    output_str_list, output_id_list= [], []
    if Source == 'MP':
        for mpid in mpid_list:
            try:
                structure_import = mprest.get_structure_by_material_id(mpid,final=final)
            except:
                print "Fail to query structure of material_id %s" % mpid
                raise

            structure_conv = SpacegroupAnalyzer(structure_import,symprec=0.1).get_conventional_standard_structure()
            output_str_list.append(structure_conv)
        return output_str_list, mpid_list

    elif str_file_list_by_user != []:
        count_str = 0
        for str_file in str_file_list_by_user:
            output_id_list.append(count_str)

            try:
                structure_import = Structure.from_file(str_file)
            except:
                print "Fail to obtain structure from file %s" % str_file
                raise
            structure_conv = SpacegroupAnalyzer(structure_import,symprec=0.1).get_conventional_standard_structure()
            output_str_list.append(structure_conv)
            count_str += 1

        return output_str_list, output_id_list
    else:
        print "Please specify the method to input structure"
        exit()

def DefEnergyDensity(matching_element,film_str,ElaC_file):

    c_axis = np.cross(matching_element[2][0],matching_element[2][3])
    Matrix_Undef, Matrix_Def= np.zeros((3,3)),np.zeros((3,3))

    for i in range(3):
        Matrix_Undef[i][0] = matching_element[2][0][i]
        Matrix_Def[i][0] = matching_element[2][0][i]*(matching_element[2][2]+1)
        Matrix_Undef[i][1] = matching_element[2][3][i]
        Matrix_Def[i][1] = matching_element[2][3][i]*(matching_element[2][5]+1)
        Matrix_Undef[i][2] = c_axis[i]
        Matrix_Def[i][2] = c_axis[i]

    F = np.transpose(np.linalg.solve(np.transpose(Matrix_Undef),np.transpose(Matrix_Def)))

    Def_Gradient = DefGradientMatrix(F)

    return DefEnergyDensity_from_file(Def_Gradient,ElaC_file)*film_str.volume / len(film_str.sites)



def Zur_Matching(film_id_list, film_str_list, sub_id_list, sub_str_list, ElaC_film_file_list):

    # loop over all substructure and film combinations
    for sub_id, sub_str  in zip(sub_id_list,sub_str_list):

        area_result, energy_result  = dict(), dict()

        for film_id, film_str in zip(film_id_list,film_str_list):

            for matching_element in structure_interface_matching(film_str,sub_str,matching_creteria):
                try:
                    energy_density = DefEnergyDensity(matching_element,film_str, ElaC_file= ElaC_film_file_list[film_id])
                except:
                    energy_density = 999

                # Sort the data by magnitude of superlattice area
                if matching_element[1] in area_result.keys():
                    dict_list = area_result[matching_element[1]]
                    energy_dict_list = energy_result[matching_element[1]]
                    if film_id in dict_list:
                        if matching_element[3][0] < dict_list[film_id]:
                            dict_list[film_id] = matching_element[3][0]
                            energy_dict_list[film_id] = energy_density
                    else:
                        dict_list[film_id] = matching_element[3][0]
                        energy_dict_list[film_id] = energy_density
                else:
                    area_result[matching_element[1]] = {film_id:matching_element[3][0]}
                    energy_result[matching_element[1]] = {film_id:energy_density}

            # organizing the oupput
            new_area_result = dict()
            for key, value in area_result.iteritems():
                new_area_result[str(key)] = value

            new_energy_result = dict()
            for key, value in energy_result.iteritems():
                new_energy_result[str(key)] = value

        # Output the data
        dumpfn(new_area_result,sub_id+'-matrix',indent=4)
        dumpfn(new_energy_result,sub_id+'-energy',indent=4)
        print sub_id, "done"


if __name__ == "__main__":

    print time.time()

    #mapi_key  = my_MP_API_Key
    mprest = MPRester()

    # Matching creteria
    matching_creteria = loadfn('matching_creteria.json')

    film_str_file_list_by_user = ['POSCAR_0']
    ElaC_film_file_list = []

    for file in film_str_file_list_by_user:
        try:
            ElaC_film_file_list.append(file+'_Cij')
        except:
            ElaC_film_file_list.append(None)

    film_str_list, film_id_list = str_prepartion(Source='User',str_file_list_by_user=film_str_file_list_by_user)
    #film_str_list = str_prepartion(Source='MP',mpid_list=['mp-149','mp-66'],final=False)

    sub_str_list, sub_id_list= str_prepartion(Source='MP',mpid_list=['mp-5229'])#'mp-1143',,'mp-1265','mp-2920','mp-5229'])


    Zur_Matching(film_id_list, film_str_list, sub_id_list, sub_str_list,ElaC_film_file_list)

    print  time.time()
