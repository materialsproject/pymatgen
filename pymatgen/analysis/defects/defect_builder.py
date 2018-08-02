"""
This class takes a list of database task objects
(for a fixed bulk composition) and then returns
Defect Entries which are ready for thermo analysis.
Includes use of Ccmpatibility class if desired.
"""

import os
import numpy as np
from monty.json import MontyDecoder, MontyEncoder

from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen import MPRester

from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.analysis.defects.core import DefectEntry
# from pymatgen.analysis.defects.corrections import find_optimal_gamma, generate_g_sum


def standardize_sc_matrix(sc_mat):
    #takes a scaling matrix = [3,3,3] or 3 or [[3,0,0],[0,3,0],[0,0,3]]
    #and standardizes it to the longer format
    if (type(sc_mat) == float) or (type(sc_mat) == int):
        return [[float(sc_mat), 0., 0.], [0., float(sc_mat), 0.], [0., 0., float(sc_mat)] ]
    elif (type( sc_mat[0]) == float) or (type( sc_mat[0]) == int):
        return [[float(sc_mat[0]),0., 0.], [0., float(sc_mat[1]), 0.], [0.,0., float(sc_mat[2])]]
    else:
        return sc_mat


def add_vr_eigenvalue_dict( task):
    """
    This takes a task object from drone output and
    does an additional vr_eigenvalue_dict key
    addition to the output key of the most
    recent calcs_reversed

    This is likely a temporary hack until
    I add this functionality to the VaspDrone?

    """
    task_dir = task['calcs_reversed'][0]['dir_name']
    vrpath = os.path.join( task_dir, 'vasprun.xml')
    if not os.path.exists(vrpath):
        vrpath += '.gz'
    if not os.path.exists(vrpath):
        vrpath = os.path.join(task_dir, 'vasprun.xml.relax2')
    if not os.path.exists(vrpath):
        vrpath += '.gz'
    if not os.path.exists(vrpath):
        print("ERROR: Could not find Vasprun file path at {}. This is needed for "
              "doing band filling correction.".format(task_dir))
        return task

    vr = Vasprun( vrpath)
    eigenvalues = {cnt:v for cnt, v in enumerate(vr.eigenvalues.values())}
    kpoint_weights = vr.actual_kpoints_weights

    vr_eigenvalue_dict = {'eigenvalues': eigenvalues, 'kpoint_weights': kpoint_weights}

    task['calcs_reversed'][0]['output']['vr_eigenvalue_dict'] = vr_eigenvalue_dict
    return task


class TaskDefectBuilder(object):
    """
    This does not have same format as a standard Builder, but it does everything we would want
    a defect builder to do.

    Input is a list of tasks from database...
    """
    def __init__(self, list_of_tasks, compatibility=DefectCompatibility()):
        self.list_of_tasks = list_of_tasks
        self.compatibility = compatibility

    def process_entries(self, run_compatibility=False):
        #first load up bulk and dielectric information
        bulk_data = {} #has keys based on scaling_matrix, values contain essential information needed for each defect entry
        diel_data = {}
        hybrid_level_data = {}
        defect_task_list = []
        mpid_checked = False

        print('Running builder for {} tasks'.format(len(self.list_of_tasks)))
        for task in self.list_of_tasks:
            #call it a dielectric calculation if LEPSILON and LPEAD are in incar
            if ('LEPSILON' in task['input']['incar'].keys()) and ('LPEAD' in task['input']['incar'].keys()):
                #process dielectric data
                eps_ionic = task['output']['epsilon_ionic']
                eps_static = task['output']['epsilon_static']
                eps_total = []
                for i in range(3):
                    eps_total.append([e[0]+e[1] for e in zip(eps_ionic[i], eps_static[i])])

                #this is a check up to see if dielectric task already loaded? for danny...
                if 'epsilon_ionic' in diel_data.keys():
                    raise ValueError("Multiple dielectric tasks exists? How to handle??")

                diel_data.update( {'epsilon_ionic': eps_ionic, 'epsilon_static':  eps_static,
                                   'dielectric': eps_total})

            #call it a bulk calculation if transformation class is just a SupercellTranformation
            elif 'SupercellTransformation' == task['transformations']['history'][0]['@class']:
                #process bulk task information
                bstruct = Structure.from_dict(task['transformations']['history'][0]['input_structure'])

                #check if mpid exists for this structure (only need to do once)
                if not mpid_checked:
                    mpid_checked = True
                    with MPRester() as mp:
                        mplist = mp.find_structure(bstruct)

                    mpid = None
                    for trial_mpid in mplist:
                        with MPRester() as mp:
                            mpstruct = mp.get_structure_by_material_id(trial_mpid)
                        if StructureMatcher(scale=False).fit(bstruct, mpstruct):
                            if mpid:
                                ValueError("Trouble Matching mp-ids. Multiple possible mpids exist?")
                            else:
                                mpid = trial_mpid

                    if mpid:
                        with MPRester() as mp:
                            bs = mp.get_bandstructure_by_material_id(mpid)
                        cbm = bs.get_cbm()['energy']
                        vbm = bs.get_vbm()['energy']
                        gga_gap = bs.get_band_gap()['energy']
                        band_stats_from_MP = True
                    else:
                        band_stats_from_MP = False

                #fetch band data if no information from MP database...
                #TODO: allow for this information to be replaced by a user's band structure calculation?
                if not band_stats_from_MP:
                    cbm = task['output']['cbm']
                    vbm = task['output']['vbm']
                    gga_gap = task['output']['bandgap']

                scaling_matrix = standardize_sc_matrix( task['transformations']['history'][0]['scaling_matrix'])
                bulk_energy = task['output']['energy']

                #check to see if bulk task of this size already loaded
                if repr(scaling_matrix) in list(bulk_data.keys()):
                    raise ValueError("Multiple bulk tasks with scaling matrix {} exist? "
                                     "How to handle??".format( scaling_matrix))

                bulk_data.update( {repr(scaling_matrix): {'bulk_energy':bulk_energy,
                                                    'bulk_structure': bstruct.copy(), 'mpid': mpid,
                                                    "cbm": cbm, "vbm": vbm, "gga_gap": gga_gap,
                                                    "band_stats_from_MP": band_stats_from_MP} } )

                if 'locpot' in task['calcs_reversed'][0]['output'].keys():
                    bulklpt = task['calcs_reversed'][0]['output']['locpot']
                    bulk_planar_averages = [bulklpt[ax] for ax in range(3)]
                    bulk_data[repr(scaling_matrix)]['bulk_planar_averages'] = bulk_planar_averages
                else:
                    print('bulk supercell {} does not have locpot values for parsing'.format(scaling_matrix))


                if 'outcar' in task['calcs_reversed'][0]['output'].keys():
                    bulkoutcar = task['calcs_reversed'][0]['output']['outcar']
                    bulk_atomic_site_averages = bulkoutcar['electrostatic_potential']
                    bulk_data[repr(scaling_matrix)]['bulk_atomic_site_averages'] = bulk_atomic_site_averages
                else:
                    print('bulk supercell {} does not have outcar values for parsing'.format(scaling_matrix))



            elif 'defect' in task['transformations']['history'][0].keys():
                defect_task_list.append(task)

            #assume that if not already identified as a defect , then a hybrid level calculation will yield additional
            # band structure information is given by task which has 'HFSCREEN' in incar...
            #TODO: figure out a better way to see if a hybrid BS caclulation is being suppled for band edge corrections?
            elif 'HFSCREEN' in task['input']['incar'].keys():
                hybrid_cbm = task['output']['cbm']
                hybrid_vbm = task['output']['vbm']
                hybrid_gap = task['output']['bandgap']

                #this is a check up to see if hybrid task already loaded? for danny...
                if hybrid_level_data:
                    raise ValueError("Hybrid level data already parsed? How to deal with this?")

                hybrid_level_data.update( {'hybrid_cbm': hybrid_cbm, 'hybrid_vbm': hybrid_vbm,
                                           'hybrid_gap': hybrid_gap})

        # """TO BE REMOVED WHEN FULL DATABASE APPROACH USED"""
        # if os.path.exists('corr_history_temp.json'):
        #     corr_hist = loadfn('corr_history_temp.json', cls=MontyDecoder)
        # else:
        #     corr_hist = {}
        #
        # if os.path.exists('kumagai_g_sums.json'):
        #     kumagai_helper = loadfn('kumagai_g_sums.json', cls=MontyDecoder)
        # else:
        #     kumagai_helper = {}
        # """END BEING REMOVED WHEN FULL DATABASE APPROACH USED"""

        #now set up and load defect entries, running compatibility as desired.
        defect_entry_list = []
        for defect_task in defect_task_list:
            #figure out size of bulk calculation and initialize it for parameters, along with dielectric data
            scaling_matrix = standardize_sc_matrix( defect_task['transformations']['history'][0]['scaling_matrix'])

            parameters = bulk_data[repr(scaling_matrix)].copy()
            parameters.update( diel_data)
            parameters.update( hybrid_level_data)
            parameters.update( {'scaling_matrix': scaling_matrix})

            #get defect object, energy and related structures
            defect = MontyDecoder().process_decoded( defect_task['transformations']['history'][0]['defect'])
            defect_energy = defect_task['output']['energy']
            bulk_struct_sc = parameters['bulk_structure'].copy()
            bulk_struct_sc.make_supercell( scaling_matrix)
            initial_defect_structure = defect.generate_defect_structure( scaling_matrix)
            final_defect_structure = Structure.from_dict( defect_task['output']['structure'])

            parameters.update( {'defect_energy': defect_energy,
                                'final_defect_structure': final_defect_structure,
                                'initial_defect_structure': initial_defect_structure} )


            # check to make sure that the bulk system which the defect is modeled after is the same as the bulk
            # calculation... this sometimes is inconsistent because a second calculation was set up with a structure
            # which was symmetrized slightly different (for example some atoms are shifted to the origin)
            bulk_calc_fin = parameters['bulk_structure'].copy() #non-supercell structure
            bulk_def_compare = defect.bulk_structure.copy() #non-supercell structure
            #compare site coords and compare
            bad_bulk_flag = False
            bad_site_compare_list = []
            for bc_site, bd_site in zip(bulk_calc_fin.sites, bulk_def_compare.sites):
                if np.linalg.norm( np.subtract(bc_site.coords, bd_site.coords)) > 0.001: #then something is wrong...
                    bad_bulk_flag = True
                    bad_site_compare_list.append( [bc_site, bd_site])
            sm = StructureMatcher( primitive_cell=False, scale=False, attempt_supercell=False, allow_subset=False)
            if bad_bulk_flag and sm.fit( bulk_calc_fin, bulk_def_compare):
                print('\tWARNING bulk structure in defect object is correct fit but appears to have shifted basis'
                      'for {}_chg{}...will not append to list. '
                      'See bad bulk sites:\n{}'.format( defect.name, defect.charge, bad_site_compare_list))
                continue
            elif not sm.fit( bulk_calc_fin, bulk_def_compare):
                print('\tFATAL ERROR bulk structure is not equivalent type between defect and bulk calculation '
                      'for {}_chg{}...will not append to list'.format( defect.name, defect.charge))
                continue

            #Load information for Freysoldt related parsing
            if 'locpot' in defect_task['calcs_reversed'][0]['output'].keys():
                deflpt = defect_task['calcs_reversed'][0]['output']['locpot']
                defect_planar_averages = [deflpt[ax] for ax in range(3)]
                abc = bulk_struct_sc.lattice.abc
                axis_grid = []
                for ax in range(3):
                    num_pts = len(defect_planar_averages[ax])
                    axis_grid.append( [i / num_pts * abc[ax] for i in range(num_pts)] )

                parameters.update( {'axis_grid': axis_grid,
                                    'defect_planar_averages': defect_planar_averages} )
            else:
                print('ERR: defect  {}_{} does not have locpot values for parsing Freysoldt'.format(defect.name, defect.charge))


            #Load information for Kumagai related parsing
            if 'outcar' in defect_task['calcs_reversed'][0]['output'].keys():
                defoutcar = defect_task['calcs_reversed'][0]['output']['outcar']
                defect_atomic_site_averages = defoutcar['electrostatic_potential']
                #NOTE: since I am not adding dim then Kumagai correction is not being performed in compatibility.
                #TODO: once Kumagai correction is fixed can modify this to allow for Kumagai correction be done.
                dim = defoutcar['ngf']

                #create list that maps site indices from bulk structure to defect structure (needed by Kumagai correction)
                site_matching_indices = []
                for dindex, dsite in enumerate(initial_defect_structure.sites):
                    if np.linalg.norm( np.subtract(dsite.coords, defect.site.coords)) > 0.001: #exclude the defect site..
                        poss_deflist = sorted(bulk_struct_sc.get_sites_in_sphere(dsite.coords, 1, include_index=True), key=lambda x: x[1])
                        bulkindex = poss_deflist[0][2]
                        site_matching_indices.append( [bulkindex, dindex])

                # assuming Wigner-Seitz radius for sampling radius
                wz = initial_defect_structure.lattice.get_wigner_seitz_cell()
                dist = []
                for facet in wz:
                    midpt = np.mean(np.array(facet), axis=0)
                    dist.append(np.linalg.norm(midpt))
                sampling_radius = min(dist)

                parameters.update( {#'dim': dim,
                                    'defect_atomic_site_averages': defect_atomic_site_averages,
                                    'site_matching_indices': site_matching_indices,
                                    'sampling_radius': sampling_radius} )
            else:
                print('ERR: defect {}_{} does not have outcar values for parsing Kumagai'.format(defect.name, defect.charge))


            #Load information for Bandfilling related parsing
            if 'vr_eigenvalue_dict' not in defect_task['calcs_reversed'][0]['output'].keys():
                defect_task = add_vr_eigenvalue_dict( defect_task)

            if 'vr_eigenvalue_dict' in defect_task['calcs_reversed'][0]['output'].keys():
                eigenvalues = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues'][:]
                kpoint_weights = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['kpoint_weights']
                parameters.update( {'eigenvalues': eigenvalues,
                                    'kpoint_weights': kpoint_weights} )
            else:
                print('ERR: defect {}_{} does not have eigenvalue data for parsing bandfilling.'.format(defect.name, defect.charge))



            defect_entry = DefectEntry( defect, parameters['defect_energy'] - parameters['bulk_energy'],
                                        corrections = {}, parameters = parameters)

            if run_compatibility:
                #first load any useful information that might exists in local json file...
                #TODO: replace this approach with loading from database (if information exists)
                """TO BE REPLACED BY DATABASE APPROACH"""
                # if task['task_id'] in corr_hist.keys():
                #     defect_entry.parameters.update( corr_hist( task['task_id']))
                # redform = bulk_struct_sc.composition.reduced_formula
                # gamma = None
                # g_sum = None
                # if redform in kumagai_helper.keys():
                #     if repr(scaling_matrix) in kumagai_helper[redform].keys():
                #         gamma = kumagai_helper[redform][repr(scaling_matrix)]['gamma']
                #         g_sum = kumagai_helper[redform][repr(scaling_matrix)]['g_sum']
                # if not gamma:
                #     print("RUNNING kumagai setup...")
                #     gamma = find_optimal_gamma(defect_struct_sc.lattice, defect_entry.parameters["dielectric"])
                #     g_sum = generate_g_sum(defect_struct_sc.lattice, defect_entry.parameters["dielectric"],
                #                            defect_entry.parameters['dim'], gamma)
                #     kumagai_helper.update( {redform: {repr(scaling_matrix): {'gamma': gamma, 'g_sum': g_sum}}})
                #
                # defect_entry.parameters.update( {'gamma': gamma, 'g_sum': g_sum})
                """END BEING REPLACED BY DATABASE APPROACH"""

                defect_entry = self.compatibility.process_entry( defect_entry)

                """TO BE REPLACED BY DATABASE APPROACH"""
                # corr_hist.update( {task['task_id']: defect_entry.parameters} )
                """END BEING REPLACED BY DATABASE APPROACH"""


            defect_entry_list.append( defect_entry)

        """TO BE REMOVED WHEN FULL DATABASE APPROACH USED"""
        # dumpfn(corr_hist, 'corr_history_temp.json', cls=MontyEncoder, indent=2)
        # dumpfn(kumagai_helper, 'kumagai_g_sums.json', cls=MontyEncoder, indent=2)
        """END BEING REMOVED WHEN FULL DATABASE APPROACH USED"""

        return defect_entry_list





