"""
This class takes a list of database task objects
(for a fixed bulk composition) and then returns
Defect Entries which are ready for thermo analysis.
Includes use of Compatibility class if desired.
"""

import os
import numpy as np
from monty.json import MontyDecoder, MontyEncoder

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen import MPRester

from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.analysis.defects.core import DefectEntry, Interstitial


class TaskDefectBuilder(object):
    """
    This does not have same format as a standard Builder, but it does everything we would want
    a defect builder to do.

    Input is a list of tasks from database...
    """
    def __init__(self, list_of_tasks,  gga_run=True, compatibility=DefectCompatibility()):
        self.list_of_tasks = list_of_tasks
        self.compatibility = compatibility
        self.gga_run = gga_run

    def process_entries(self, run_compatibility=False):
        """
        :param run_compatibility (bool): Whether to run defect_compatibility or not
        :param gga_run:  needed to determine whether MP entry is relevant for chemical potentiials or not
        :return:

        """
        #first load up bulk and dielectric information
        bulk_data = {} #has keys based on number of atoms, values contain essential information needed for each defect entry
        diel_data = {}
        hybrid_level_data = {}
        defect_task_list = []
        mpid = None
        #TODO: add confirmation that all tasks are derived from same approach (gga vs. hse) / same structure

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
                continue

            elif ('history' in task['transformations'].keys()):
                if ('DefectTransformation' in task['transformations']['history'][0]['@class']):
                    defect_task_list.append(task)

                elif 'SupercellTransformation' == task['transformations']['history'][0]['@class']:
                    #process bulk task information
                    bulk_sc_structure = task['input']['structure']
                    if type(bulk_sc_structure) == dict:
                        bulk_sc_structure = Structure.from_dict( bulk_sc_structure)

                    #check if mpid exists for this structure (only need to do once)
                    if mpid is None:
                        mpid = self.check_mpid( bulk_sc_structure)
                        if (mpid == 'DNE') or (not self.gga_run):
                            # if no Mp-id exists (or if calculation is not a GGA type calcualtion)
                            # then use cbm and vbm from bulk calc
                            #TODO: allow for this information to be replaced by a user's band structure calculation?
                            cbm = task['output']['cbm']
                            vbm = task['output']['vbm']
                            gga_gap = task['output']['bandgap']
                        elif mpid:
                            with MPRester() as mp:
                                bs = mp.get_bandstructure_by_material_id(mpid)
                            cbm = bs.get_cbm()['energy']
                            vbm = bs.get_vbm()['energy']
                            gga_gap = bs.get_band_gap()['energy']

                    bulk_energy = task['output']['energy']

                    #check to see if bulk task of this size already loaded
                    if len(bulk_sc_structure) in list(bulk_data.keys()):
                        raise ValueError("Multiple bulk tasks with size {} atoms exist? "
                                         "How to handle??".format( len(bulk_sc_structure)))

                    bulk_data.update( {len(bulk_sc_structure): {'bulk_energy':bulk_energy, 'mpid': mpid,
                                                                "cbm": cbm, "vbm": vbm, "gga_gap": gga_gap,
                                                                'bulk_sc_structure': bulk_sc_structure} } )

                    if 'locpot' in task['calcs_reversed'][0]['output'].keys():
                        bulklpt = task['calcs_reversed'][0]['output']['locpot']
                        axes = list(bulklpt.keys())
                        axes.sort()
                        bulk_data[len(bulk_sc_structure)]['bulk_planar_averages'] = [bulklpt[ax] for ax in axes]
                    else:
                        print('bulk supercell with {}atoms does not have locpot values for parsing'.format(len(bulk_sc_structure)))


                    if 'outcar' in task['calcs_reversed'][0]['output'].keys():
                        bulkoutcar = task['calcs_reversed'][0]['output']['outcar']
                        bulk_atomic_site_averages = bulkoutcar['electrostatic_potential']
                        bulk_data[len(bulk_sc_structure)]['bulk_atomic_site_averages'] = bulk_atomic_site_averages
                    else:
                        print('bulk supercell {} does not have outcar values for parsing'.format(len(bulk_sc_structure)))

            # assume that if not already identified as a defect or bulk supercel, then a new
            # hybrid level calculation is for adjusting band edge correction
            #TODO: figure out a better way to see if a hybrid BS caclulation
            #    is being suppled for band edge corrections?
            elif 'HFSCREEN' in task['input']['incar'].keys():
                hybrid_cbm = task['output']['cbm']
                hybrid_vbm = task['output']['vbm']
                hybrid_gap = task['output']['bandgap']

                #this is a check up to see if hybrid task already loaded? for danny...
                if hybrid_level_data:
                    raise ValueError("Hybrid level data already parsed? How to deal with this?")

                hybrid_level_data.update( {'hybrid_cbm': hybrid_cbm, 'hybrid_vbm': hybrid_vbm,
                                           'hybrid_gap': hybrid_gap})


        #now set up and load defect entries, running compatibility as desired.
        defect_entry_list = []
        for defect_task in defect_task_list:
            #figure out size of bulk calculation and initialize it for parameters, along with dielectric data
            final_defect_structure = Structure.from_dict( defect_task['output']['structure'])

            map_bulk_struct = {abs(size - len(final_defect_structure)): size for size in bulk_data.keys()}
            bulk_size = map_bulk_struct[min(map_bulk_struct.keys())]

            if (min(map_bulk_struct.keys()) > 1):
                raise ValueError("ERROR no bulk tasks exist for defect sc "
                                 "size {}".format(len(final_defect_structure)))

            parameters = bulk_data[bulk_size].copy()
            parameters.update( diel_data)
            parameters.update( hybrid_level_data)

            #get defect object, energy and related structures
            defect = MontyDecoder().process_decoded( defect_task['transformations']['history'][0]['defect'])
            scaling_matrix = MontyDecoder().process_decoded( defect_task['transformations']['history'][0]['scaling_matrix'])
            initial_defect_structure = defect.generate_defect_structure( scaling_matrix)
            initial_defect_structure = self.reorder_structure( initial_defect_structure, final_defect_structure)
            defect_energy = defect_task['output']['energy']

            parameters.update( {'defect_energy': defect_energy,
                                'final_defect_structure': final_defect_structure,
                                'initial_defect_structure': initial_defect_structure} )

            #Load information for Freysoldt related parsing
            if 'locpot' in defect_task['calcs_reversed'][0]['output'].keys():
                deflpt = defect_task['calcs_reversed'][0]['output']['locpot']
                axes = list(deflpt.keys())
                axes.sort()
                defect_planar_averages = [deflpt[ax] for ax in axes]

                abc = initial_defect_structure.lattice.abc
                axis_grid = []
                for ax in range(3):
                    num_pts = len(defect_planar_averages[ax])
                    axis_grid.append( [i / num_pts * abc[ax] for i in range(num_pts)] )

                parameters.update( {'axis_grid': axis_grid,
                                    'defect_planar_averages': defect_planar_averages} )
            else:
                print('ERR: defect  {}_{} does not have locpot values for parsing Freysoldt correction'.format(defect.name, defect.charge))


            #Load information for Kumagai related parsing
            if 'outcar' in defect_task['calcs_reversed'][0]['output'].keys():
                defoutcar = defect_task['calcs_reversed'][0]['output']['outcar']
                defect_atomic_site_averages = defoutcar['electrostatic_potential']
                bulk_sc_structure = parameters['bulk_sc_structure']

                #find fractional coordinates of defect within the final_defect_structure
                if type(defect) != Interstitial:
                    poss_deflist = sorted(
                        bulk_sc_structure.get_sites_in_sphere(defect.site.coords, 2, include_index=True), key=lambda x: x[1])
                    defect_frac_sc_coords = bulk_sc_structure[poss_deflist[0][2]].frac_coords
                else:
                    poss_deflist = sorted(
                        initial_defect_structure.get_sites_in_sphere(defect.site.coords, 2, include_index=True), key=lambda x: x[1])
                    defect_frac_sc_coords = initial_defect_structure[poss_deflist[0][2]].frac_coords

                #create list that maps site indices from bulk structure to defect structure (needed by Kumagai correction)
                site_matching_indices = []
                for dindex, dsite in enumerate(initial_defect_structure.sites):
                    if dsite.distance_and_image_from_frac_coords( defect_frac_sc_coords)[0] > 0.001:  #exclude the defect site from site_matching...
                        poss_deflist = sorted(bulk_sc_structure.get_sites_in_sphere(dsite.coords, 1, include_index=True), key=lambda x: x[1])
                        bulkindex = poss_deflist[0][2]
                        site_matching_indices.append( [bulkindex, dindex])

                # assuming Wigner-Seitz radius for sampling radius
                wz = initial_defect_structure.lattice.get_wigner_seitz_cell()
                dist = []
                for facet in wz:
                    midpt = np.mean(np.array(facet), axis=0)
                    dist.append(np.linalg.norm(midpt))
                sampling_radius = min(dist)

                parameters.update( {'defect_atomic_site_averages': defect_atomic_site_averages,
                                    'site_matching_indices': site_matching_indices,
                                    'sampling_radius': sampling_radius,
                                    'defect_frac_sc_coords': defect_frac_sc_coords} )
            else:
                print('ERR: defect {}_{} does not have outcar values for parsing Kumagai'.format(defect.name, defect.charge))


            if 'vr_eigenvalue_dict' in defect_task['calcs_reversed'][0]['output'].keys():
                eigenvalues = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues']
                kpoint_weights = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['kpoint_weights']
                parameters.update( {'eigenvalues': eigenvalues,
                                    'kpoint_weights': kpoint_weights} )
            else:
                print('ERR: defect {}_{} does not have eigenvalue data for parsing '
                      'bandfilling.'.format(defect.name, defect.charge))

            if 'defect' in defect_task['calcs_reversed'][0]['output'].keys():
                parameters.update( {'defect_ks_delocal_data': defect_task['calcs_reversed'][0]['output']['defect']})
            else:
                print('ERR: defect {}_{} does not have defect data for parsing '
                      'delocalization.'.format(defect.name, defect.charge))


            defect_entry = DefectEntry( defect, parameters['defect_energy'] - parameters['bulk_energy'],
                                        corrections = {}, parameters = parameters)

            if run_compatibility:
                defect_entry = self.compatibility.process_entry( defect_entry)

            dir_name = defect_task['calcs_reversed'][0]['dir_name']
            defect_entry.parameters.update( {'dir_name': dir_name, 'task_db_task_id': defect_task['_id']})

            defect_entry_list.append( defect_entry)


        return defect_entry_list

    def check_mpid(self, bulk_sc_structure):
        with MPRester() as mp:
            mplist = mp.find_structure(bulk_sc_structure)

        mpid_fit_list = []
        for trial_mpid in mplist:
            with MPRester() as mp:
                mpstruct = mp.get_structure_by_material_id(trial_mpid)
            if StructureMatcher(scale=False).fit(bulk_sc_structure, mpstruct):
                mpid_fit_list.append( trial_mpid)

        if len(mpid_fit_list) == 1:
            mpid = mpid_fit_list[0]
        elif len(mpid_fit_list) > 1:
            num_mpid_list = [int(mp.split('-')[1]) for mp in mpid_fit_list]
            num_mpid_list.sort()
            mpid  = 'mp-'+str(num_mpid_list[0])
            print("Multiple mp-ids found for bulk structure:{}\nWill use lowest number mpid "
                  "for bulk band structure = {}.".format(str(mpid_fit_list), mpid))
        else:
            print("Could not find bulk structure in MP database after tying the following list:\n{}"
                  "Continuing with band_status_from_MP=False.".format( mplist))
            mpid = 'DNE'

        return mpid

    def reorder_structure( self, s1, s2):
        """
        This takes s1 and attempts to reorder the structure in same site-indexed fashion as s2

        If substantial relaxation occurs this might cause errors...

        returns: reordered s1 structure
        """
        if s1.composition != s2.composition:
            raise ValueError('Compositions are not equivalent: {} vs. {}'.format( s1.composition, s2.composition))

        list_matcher = []
        for s2ind, s2site in enumerate(s2):
            match_this_site = [[s2site.distance_and_image_from_frac_coords( s1site.frac_coords)[0],
                               s1ind] for s1ind, s1site in enumerate( s1)]
            match_this_site.sort() #put closest site first
            for match in match_this_site:
                #take first site matching that matches specie
                if s1[match[1]].specie == s2site.specie:
                    list_matcher.append( [ match[1], s2ind])
                    break

        s1_matcher_set = list(np.array( list_matcher)[:, 0])
        num_unique_s1_sites_matched = len(set(s1_matcher_set))
        s2_matcher_set = list(np.array( list_matcher)[:, 1])
        num_unique_s2_sites_matched = len(set(s2_matcher_set))

        if (len(s1) != num_unique_s1_sites_matched) or \
            (len(s2) != num_unique_s2_sites_matched):
            raise ValueError("Number of unique sites differs in defect site matching routine." \
                  "\n\ts1 {} vs. {}\n\ts2 {} vs. {}".format(len(s1), num_unique_s1_sites_matched,
                                                            len(s2), num_unique_s2_sites_matched))

        new_species_list = [s1[ind].specie for ind in s1_matcher_set]
        new_coords_list = [s1[ind].frac_coords for ind in s1_matcher_set]
        new_struct = Structure( s1.lattice, new_species_list, new_coords_list,
                               validate_proximity=True, to_unit_cell=True,
                               coords_are_cartesian=False)

        return new_struct

    def process_and_push_to_db(self, db, run_compatibility=False, gga_run=True):
        """
        For running process_entries and then immediately pushing to database
        :param db:
        :param tags:
        :param run_compatibility: (bool) whether to run_compatibility or not...not really needed?
        :param gga_run: (whether this run is a gga run or not - impacts HSE bs parsing...)
        :return:
        """
        entries = self.process_entries(run_compatibility=run_compatibility, gga_run=gga_run)

        for entry in entries:
            #TODO: probably dont want to do this with storage of the massive delocalization data...
            t_id = db.insert_task( entry) #note that this updates previous entries too...



