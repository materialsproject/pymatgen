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
            # mplist = mp.find_structure(bulk_sc_structure) #had to hack this because this wasnt working??
            tmp_mplist = mp.get_entries_in_chemsys(list(bulk_sc_structure.symbol_set))

        mplist = [ment.entry_id for ment in tmp_mplist if ment.composition.reduced_composition == \
                  bulk_sc_structure.composition.reduced_composition]
        #TODO: this is a hack because find_structure was data intensive. simplify the hack to do less queries...

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



"""Begin real builders (not hacks)"""

from datetime import datetime
from itertools import chain
import numpy as np

from monty.json import MontyDecoder, jsanitize

from pymatgen import Structure
from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.defects.core import Interstitial, DefectEntry
from pymatgen import MPRester
from pymatgen.electronic_structure.bandstructure import BandStructure

from maggma.builder import Builder



__author__ = "Danny Broberg <dpbroberg@lbl.gov>"


class DefectBuilder(Builder):
    def __init__(self,
                 tasks,
                 defects,
                 query=None,
                 compatibility=DefectCompatibility(),
                 ltol=0.2,
                 stol=0.3,
                 angle_tol=5,
                 max_items_size=0,
                 **kwargs):
        """
        Creates DefectEntry from vasp task docs

        Args:
            tasks (Store): Store of tasks documents
            defects (Store): Store of defect entries with all metadata required for followup decisions on defect thermo
            query (dict): dictionary to limit materials to be analyzed
            compatibility (PymatgenCompatability): Compatability module to ensure defect calculations are compatible
            ltol (float): StructureMatcher tuning parameter for matching tasks to materials
            stol (float): StructureMatcher tuning parameter for matching tasks to materials
            angle_tol (float): StructureMatcher tuning parameter for matching tasks to materials
            max_items_size (int): limits number of items approached from tasks (zero places no limit on number of items)
        """

        self.tasks = tasks
        self.defects = defects
        self.query = query if query else {}
        self.compatibility = compatibility
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        self.max_items_size = max_items_size
        super().__init__(sources=[tasks], targets=[defects], **kwargs)

    def get_items(self):
        """
        Gets sets of entries from chemical systems that need to be processed

        Returns:
            generator of relevant entries from one chemical system

        """
        self.logger.info("Defect Builder Started")

        self.logger.info("Setting indexes")
        self.ensure_indicies() #TODO: is this neccessary? Is there a better way to use it?

        # Save timestamp for update operation
        self.time_stamp = datetime.utcnow()

        # Get all successful defect tasks that have been updated since
        # defect_store was last updated
        q = dict(self.query)
        q["state"] = "successful"
        # q.update(self.defects.lu_filter(self.defects)) #This wasnt working because DefectEntries dont have obvious way of doing this? (tried parameters.last_updated but this broke because parameters is a property and last_updated is a key)
        #TODO: does self.tasks.lu_filter(self.defects) work better?
        q.update({'task_id': {'$nin': self.defects.distinct('entry_id')}}) #dont redo previous tasks...
        # q.update({'transformations.history.0.@module':
        q.update({'transformations.history.@module':
                      {'$in': ['pymatgen.transformations.defect_transformations']}})
        defect_tasks = list(self.tasks.query(criteria=q,
                                             properties=['task_id', 'transformations', 'input',
                                                         'task_label', 'last_updated',
                                                         'output', 'calcs_reversed', 'chemsys']))
        if self.max_items_size and len(defect_tasks) > self.max_items_size:
            defect_tasks = [dtask for dind, dtask in enumerate(defect_tasks) if dind < self.max_items_size]
        task_ids = [dtask['task_id'] for dtask in defect_tasks]
        self.logger.info("Found {} new defect tasks to consider:\n{}".format( len(defect_tasks), task_ids))
        log_defect_bulk_types = [frozenset(Structure.from_dict(dt['transformations']['history'][0]['defect']['structure']).symbol_set)
                                 for dt in defect_tasks]
        log_defect_bulk_types = list(set( log_defect_bulk_types))

        #get a few other tasks which are needed for defect entries (regardless of when they were last updated):
        #   bulk_supercell, dielectric calc, BS calc, HSE-BS calc
        bulksc = {"state": "successful", 'transformations.history.0.@module':
            {'$in': ['pymatgen.transformations.standard_transformations']}} #ALSO -> confirm identical INCAR settings are set?
        dielq = {"state": "successful", "input.incar.LEPSILON": True, "input.incar.LPEAD": True}
        HSE_BSq = {"state": "successful", 'calcs_reversed.0.input.incar.LHFCALC': True,
                   'transformations.history.0.@module':
                        {'$nin': ['pymatgen.transformations.defect_transformations',
                                  'pymatgen.transformations.standard_transformations']}}
        # TODO: add smarter capability for getting HSE bandstructure from tasks
        # TODO: add capability for getting GGA bandstructure from tasks? (Currently uses bulk supercell cbm/vbm when no MPID exists...)

        all_bulk_tasks = list(self.tasks.query(criteria=bulksc, properties=['task_id', 'chemsys', 'last_updated',
                                                                            'input', 'task_label']))
        self.logger.info('Queried {} bulk calculations'.format( len(all_bulk_tasks)))
        log_additional_tasks = dict() #organize based on symbols to save a bit of parsing time
        for blktask in all_bulk_tasks:
            sym_set = frozenset(blktask['chemsys'].split('-'))
            if sym_set not in log_defect_bulk_types:
                continue

            if sym_set not in log_additional_tasks.keys():
                log_additional_tasks[sym_set] = {'bulksc': [blktask.copy()]}
            else:
                log_additional_tasks[sym_set]['bulksc'].append( blktask.copy())

            #grab diel
            if 'diel' not in log_additional_tasks[sym_set].keys():
                q = dielq.copy()
                q.update( {'elements': {'$all': list( sym_set)}})
                diel_tasks = list(self.tasks.query(criteria=q,
                                                   properties=['task_id', 'task_label', 'last_updated',
                                                               'input', 'output']))
                log_additional_tasks[sym_set]['diel'] = diel_tasks[:]

            #grab hse bs
            if 'hsebs' not in log_additional_tasks[sym_set].keys():
                q = HSE_BSq.copy()
                q.update( {'elements': {'$all': list( sym_set)}})
                hybrid_tasks = list(self.tasks.query(criteria=q,
                                                     properties=['task_id', 'input', 'output', 'task_label']))
                log_additional_tasks[sym_set]['hsebs'] = hybrid_tasks[:]

        self.logger.info('Populated bulk, diel, and hse bs lists')
        #now load up all defect tasks with relevant information required for analysis
        temp_log_bs_bulk = dict() #to minimize number of band structure queries to MP, log by element sets
        for d_task in defect_tasks:
            sym_set = frozenset(Structure.from_dict(d_task['transformations']['history'][0]['defect']['structure']).symbol_set)
            defect_task = self.load_defect_task( d_task, log_additional_tasks[sym_set])

            #import additional bulk task information and get mp-id as appropriate
            bulk_task_id = defect_task['bulk_task']['task_id']
            if bulk_task_id not in temp_log_bs_bulk.keys():
                full_bulk_task = list(self.tasks.query(criteria={'task_id': bulk_task_id},
                                                  properties=['task_id', 'chemsys', 'task_label',
                                                              'transformations',
                                                              'input', 'output', 'calcs_reversed']))
                if len(full_bulk_task) != 1:
                    raise ValueError("got {} bulk task for id {}?".format( len(full_bulk_task), bulk_task_id))
                else:
                    full_bulk_task = full_bulk_task[0]

                bulk_structure = Structure.from_dict(full_bulk_task['input']['structure'])
                mpid = self.get_bulk_mpid( bulk_structure)

                if mpid:
                    with MPRester() as mp:
                        bs = mp.get_bandstructure_by_material_id(mpid)
                    full_bulk_task['MP-gga-BScalc'] = bs.as_dict().copy()
                    full_bulk_task['mpid'] = mpid
                else:
                    full_bulk_task['MP-gga-BScalc'] = None
                    full_bulk_task['mpid'] = None
                temp_log_bs_bulk[bulk_task_id] = full_bulk_task.copy()
            else:
                full_bulk_task = temp_log_bs_bulk[bulk_task_id].copy()

            defect_task['bulk_task'] = full_bulk_task.copy()

            if defect_task:
                yield defect_task


    def process_item(self, item):
        """
        Process a defect item (containing defect, bulk and dielectric information as processed in get_items)

        Args:
            item (defect_task): a defect_task to process into a DefectEntry object

        Returns:
            dict: a DefectEntry dictionary to update defect database with
        """

        if 'bulk_task' not in item:
            raise ValueError("bulk_task is not in item! Cannot parse this.")
        elif 'diel_task_meta' not in item:
            raise ValueError("diel_task_meta is not in item! Cannot parse this.")

        #initialize parameters with dielectric data
        eps_ionic = item['diel_task_meta']['epsilon_ionic']
        eps_static = item['diel_task_meta']['epsilon_static']
        eps_total = []
        for i in range(3):
            eps_total.append([e[0]+e[1] for e in zip(eps_ionic[i], eps_static[i])])
        parameters = {'epsilon_ionic': eps_ionic, 'epsilon_static':  eps_static,
                      'dielectric': eps_total,
                      'task_level_metadata':
                          {'diel_taskdb_task_id': item['diel_task_meta']['diel_taskdb_task_id']}}

        #initialize bulk data in parameters
        bulk_task = item['bulk_task']
        bulk_energy = bulk_task['output']['energy']
        mpid = bulk_task['mpid']
        bulk_sc_structure = bulk_task['input']['structure']
        if type(bulk_sc_structure) == dict:
            bulk_sc_structure = Structure.from_dict(bulk_sc_structure)

        parameters.update( {'bulk_energy': bulk_energy, 'mpid': mpid,
                            'bulk_sc_structure': bulk_sc_structure} )

        if 'locpot' in bulk_task['calcs_reversed'][0]['output'].keys():
            bulklpt = bulk_task['calcs_reversed'][0]['output']['locpot']
            axes = list(bulklpt.keys())
            axes.sort()
            parameters.update( {'bulk_planar_averages': [bulklpt[ax] for ax in axes]} )
        else:
            self.logger.error('BULKTYPEcalc: {} (task-id {}) does not '
                              'have locpot values for parsing'.format( bulk_task['task_label'],
                                                                       bulk_task['task_id']))

        if 'outcar' in bulk_task['calcs_reversed'][0]['output'].keys():
            bulkoutcar = bulk_task['calcs_reversed'][0]['output']['outcar']
            bulk_atomic_site_averages = bulkoutcar['electrostatic_potential']
            parameters.update( {'bulk_atomic_site_averages': bulk_atomic_site_averages})
        else:
            self.logger.error('BULKTYPEcalc: {} (task-id {}) does not '
                              'have outcar values for parsing'.format( bulk_task['task_label'],
                                                                       bulk_task['task_id']))

        #load INCAR, KPOINTS, POTCAR and task_id (for both bulk and defect)
        potcar_summary = {'pot_spec': list([potelt["titel"] for potelt in item['input']['potcar_spec']]),
                            'pot_labels': list(item['input']['pseudo_potential']['labels'][:]),
                            'pot_type': item['input']['pseudo_potential']['pot_type'],
                            'functional': item['input']['pseudo_potential']['functional']} #note bulk has these potcar values also, other wise it would not get to process_items
        dincar = item["input"]["incar"].copy()
        dincar_reduced = {k: dincar.get(k, None) for k in ["LHFCALC", "HFSCREEN", "IVDW", "LUSE_VDW",
                                                           "LDAU", "METAGGA"]} #same as bulk
        bincar = item["input"]["incar"].copy()
        d_kpoints = item['calcs_reversed'][0]['input']['kpoints']
        if type(d_kpoints) != dict:
            d_kpoints = d_kpoints.as_dict()
        b_kpoints = bulk_task['calcs_reversed'][0]['input']['kpoints']
        if type(b_kpoints) != dict:
            b_kpoints = b_kpoints.as_dict()

        dir_name = item['calcs_reversed'][0]['dir_name']
        bulk_dir_name = bulk_task['calcs_reversed'][0]['dir_name']

        parameters['task_level_metadata'].update( {'defect_dir_name': dir_name, 'bulk_dir_name': bulk_dir_name,
                                                   'bulk_taskdb_task_id': bulk_task['task_id'],
                 'potcar_summary': potcar_summary.copy(), 'incar_calctype_summary': dincar_reduced.copy(),
                 'defect_incar': dincar.copy(), 'bulk_incar': bincar.copy(),
                 'defect_kpoints': d_kpoints.copy(), 'bulk_kpoints': b_kpoints.copy(),
                                                   'defect_task_last_updated': item['last_updated']})

        #load band edge characteristics
        if mpid:
            #TODO: NEED to be smarter about use of +U etc in the MP gga band structure calculations...
            bs = bulk_task['MP-gga-BScalc']
            if type(bs) == dict:
                bs = BandStructure.from_dict( bs)

            parameters['task_level_metadata'].update( {'MP_gga_BScalc_data':
                                                           bs.get_band_gap().copy()} ) #contains gap kpt transition
            cbm = bs.get_cbm()['energy']
            vbm = bs.get_vbm()['energy']
            gap = bs.get_band_gap()['energy']
        else:
            cbm = bulk_task['output']['cbm']
            vbm = bulk_task['output']['vbm']
            gap = bulk_task['output']['bandgap']

        parameters.update( {"cbm": cbm, "vbm": vbm, "gap": gap})

        if (dincar_reduced['HFSCREEN'] in [None, False, 'False']) and ("hybrid_bs_meta" in item):
            parameters.update({k: item["hybrid_bs_meta"][k] for k in ['hybrid_cbm', 'hybrid_vbm']})
            parameters.update({'hybrid_gap': parameters['hybrid_cbm'] - parameters['hybrid_vbm']})
            parameters['task_level_metadata'].update( {k: item["hybrid_bs_meta"][k] for k in ['hybrid_CBM_task_id',
                                                                                              'hybrid_VBM_task_id']})

        # get defect object, energy and related structures
        final_defect_structure = item['output']['structure']
        if type(final_defect_structure) != Structure:
            final_defect_structure = Structure.from_dict( final_defect_structure)
        if type( item['transformations']) != dict:
            item['transformations'] = item['transformations'].as_dict()
        defect = item['transformations']['history'][0]['defect']
        try:
            defect = MontyDecoder().process_decoded( defect)
        except:
            self.logger.info("Error with defect loading for task-id {}. "
                             "Trying again after removing unnecessary keys.".format( item['task_id']))
            needed_keys = ['@module', '@class', 'structure', 'defect_site', 'charge']
            defect = MontyDecoder().process_decoded( {k:v for k,v in defect.items() if k in needed_keys})

        scaling_matrix = MontyDecoder().process_decoded( item['transformations']['history'][0]['scaling_matrix'])
        initial_defect_structure = defect.generate_defect_structure(scaling_matrix)
        defect_energy = item['output']['energy']
        try:
            initial_defect_structure = self.reorder_structure(initial_defect_structure, final_defect_structure)
            parameters.update({'stdrd_init_to_final_structure_match': True})
        except:
            #above can fail in cases of large relaxation. If this is case, then can sometimes rely on input.structure
            #for initial_defect_structure...this can cause minor problems if the defect required re-running...
            self.logger.debug("WARNING: had issue with reordering_structure from the defect structure "
                              "description. Switching to use of input.structure and confirming that "
                              "species are identical and proceeding...Note that this can possibly cause "
                              "problems with potential alignment or structure analysis later on")

            initial_defect_structure = item['input']['structure']
            if type( initial_defect_structure) != Structure:
                initial_defect_structure = Structure.from_dict(initial_defect_structure)
            for index, u, v in zip(range(len(initial_defect_structure)),
                                   initial_defect_structure, final_defect_structure):
                if u.specie != v.specie:
                    raise ValueError("Could not match index {}. {} != {}".format(index, u, v))

            parameters.update({'stdrd_init_to_final_structure_match': False})

        parameters.update({'defect_energy': defect_energy,
                           'final_defect_structure': final_defect_structure,
                           'initial_defect_structure': initial_defect_structure})

        if 'locpot' in item['calcs_reversed'][0]['output'].keys():
            #Load information for Freysoldt related parsing
            deflpt = item['calcs_reversed'][0]['output']['locpot']
            axes = list(deflpt.keys())
            axes.sort()
            defect_planar_averages = [deflpt[ax] for ax in axes]
            abc = initial_defect_structure.lattice.abc
            axis_grid = []
            for ax in range(3):
                num_pts = len(defect_planar_averages[ax])
                axis_grid.append([i / num_pts * abc[ax] for i in range(num_pts)])
            parameters.update({'axis_grid': axis_grid,
                               'defect_planar_averages': defect_planar_averages})
        else:
            self.logger.error('DEFECTTYPEcalc: {} (task-id {}) does not have locpot values for '
                              'parsing Freysoldt correction'.format(item['task_label'], item['task_id']))


        if 'outcar' in item['calcs_reversed'][0]['output'].keys():
            #Load information for Kumagai related parsing
            defoutcar = item['calcs_reversed'][0]['output']['outcar']
            defect_atomic_site_averages = defoutcar['electrostatic_potential']
            bulk_sc_structure = parameters['bulk_sc_structure']

            if type(defect) != Interstitial:
                poss_deflist = sorted(
                    bulk_sc_structure.get_sites_in_sphere(defect.site.coords, 2,
                                                          include_index=True), key=lambda x: x[1])
                defect_frac_sc_coords = bulk_sc_structure[poss_deflist[0][2]].frac_coords
            else:
                poss_deflist = sorted(
                    initial_defect_structure.get_sites_in_sphere(defect.site.coords, 2,
                                                                 include_index=True), key=lambda x: x[1])
                defect_frac_sc_coords = initial_defect_structure[poss_deflist[0][2]].frac_coords

            #create list that maps site indices from bulk structure to defect structure
            site_matching_indices = []
            for dindex, dsite in enumerate(initial_defect_structure.sites):
                if dsite.distance_and_image_from_frac_coords( defect_frac_sc_coords)[0] > 0.001:  #exclude the defect site from site_matching...
                    poss_deflist = sorted(bulk_sc_structure.get_sites_in_sphere(dsite.coords, 1,
                                                                                include_index=True), key=lambda x: x[1])
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
            self.logger.error('DEFECTTYPEcalc: {} (task-id {}) does not have outcar values for '
                              'parsing Kumagai'.format(item['task_label'], item['task_id']))

        if 'vr_eigenvalue_dict' in item['calcs_reversed'][0]['output'].keys():
            eigenvalues = item['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues']
            kpoint_weights = item['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['kpoint_weights']
            parameters.update( {'eigenvalues': eigenvalues,
                                'kpoint_weights': kpoint_weights} )
        else:
            self.logger.error('DEFECTTYPEcalc: {} (task-id {}) does not have eigenvalue data for parsing '
                  'bandfilling.'.format(item['task_label'], item['task_id']))

        if 'defect' in item['calcs_reversed'][0]['output'].keys():
            parameters.update( {'defect_ks_delocal_data': item['calcs_reversed'][0]['output']['defect']})
        else:
            self.logger.error('DEFECTTYPEcalc: {} (task-id {}) does not have defect data for parsing '
                  'delocalization.'.format(item['task_label'], item['task_id']))


        defect_entry = DefectEntry( defect, parameters['defect_energy'] - parameters['bulk_energy'],
                                    corrections = {}, parameters = parameters, entry_id= item['task_id'])

        defect_entry = self.compatibility.process_entry( defect_entry)
        defect_entry.parameters['last_updated'] = datetime.utcnow()

        #add additional tags as desired...
        dentry_as_dict = defect_entry.as_dict()
        dentry_as_dict['task_id'] = item['task_id'] #this will need to be deleted when loading DefectEntry.as_dict()

        return dentry_as_dict

    def update_targets(self, items):
        """
        Inserts the defect docs into the defect collection

        Args:
            items ([dict]): a list of defect entries as dictionaries
        """

        self.logger.info("Updating {} defect documents".format(len(items)))

        self.defects.update(items, update_lu=True, key='entry_id')

    def ensure_indicies(self):
        """
        Ensures indicies on the tasks and defects collections
        :return:
        """
        # Search indicies for tasks
        self.tasks.ensure_index(self.tasks.key, unique=True)
        # self.tasks.ensure_index(self.tasks.lu_field)
        self.tasks.ensure_index("chemsys")

        # Search indicies for defects
        self.defects.ensure_index(self.defects.key, unique=True)
        # self.defects.ensure_index(self.defects.lu_field)
        self.defects.ensure_index("chemsys")

    def load_defect_task(self, defect_task, additional_tasks):
        """
        This takes defect_task and finds bulk task, diel task, and other things as appropriate (example hse_data...)

        needs to make sure INCAR settings and structure matching works for all

        :param defect_task:
        :param additional_tasks:
        :return:
        """
        out_defect_task = defect_task.copy()

        # get suitable BULK calc by matching structures (right lattice
        # constant + same supercell size...) and matching essential INCAR + POTCAR settings...
        bulk_tasks = additional_tasks['bulksc']
        bulk_sm = StructureMatcher( ltol=self.ltol, stol=self.stol, angle_tol=self.angle_tol,
            primitive_cell=False, scale=False, attempt_supercell=False, allow_subset=False)
        bulk_matched = []
        dstruct_withoutdefect = Structure.from_dict(out_defect_task['transformations']['history'][0]['defect']['structure'])
        scaling_matrix = out_defect_task['transformations']['history'][0]['scaling_matrix']
        dstruct_withoutdefect.make_supercell( scaling_matrix)
        dincar = out_defect_task["input"]["incar"] #identify essential INCAR properties which differentiate different calcs
        dincar_reduced = {k: dincar.get(k, None) for k in ["LHFCALC", "HFSCREEN", "IVDW", "LUSE_VDW",
                                                           "LDAU", "METAGGA"]}
        d_potcar_base = {'pot_spec': [potelt["titel"] for potelt in out_defect_task['input']['potcar_spec']],
                         'pot_labels': out_defect_task['input']['pseudo_potential']['labels'][:],
                         'pot_type': out_defect_task['input']['pseudo_potential']['pot_type'],
                         'functional': out_defect_task['input']['pseudo_potential']['functional']}

        for b_task in bulk_tasks:
            bstruct = Structure.from_dict(b_task['input']['structure'])

            if bulk_sm.fit( bstruct, dstruct_withoutdefect):
                #also match essential INCAR and POTCAR settings
                bincar = b_task["input"]["incar"]
                bincar_reduced = {k: bincar.get(k, None) for k in dincar_reduced.keys()}

                b_potcar = {'pot_spec': set([potelt["titel"] for potelt in b_task['input']['potcar_spec']]),
                            'pot_labels': set(b_task['input']['pseudo_potential']['labels'][:]),
                            'pot_type': b_task['input']['pseudo_potential']['pot_type'],
                            'functional': b_task['input']['pseudo_potential']['functional']}
                d_potcar = d_potcar_base.copy() #need to reduce in cases of extrinsic species or reordering
                d_potcar['pot_spec'] = set([d for d in d_potcar['pot_spec'] if d in b_potcar['pot_spec']])
                d_potcar['pot_labels'] = set([d for d in d_potcar['pot_labels'] if d in b_potcar['pot_labels']])

                #track to make sure that cartesian coords are same (important for several levels of analysis in defect builder)
                same_cart_positions = True
                for bsite_coords in bstruct.cart_coords:
                    if not len( dstruct_withoutdefect.get_sites_in_sphere(bsite_coords, 1)):
                        same_cart_positions = False

                if bincar_reduced == dincar_reduced and b_potcar == d_potcar and same_cart_positions:
                    bulk_matched.append( b_task.copy())
                else:
                    self.logger.debug("Bulk structure match was found for {} with {}, "
                                      "but:".format( b_task['task_label'], out_defect_task['task_label']))
                    if bincar_reduced != dincar_reduced:
                        out_inc = {k:[v, bincar_reduced[k]] for k,v in dincar_reduced.items() if v != bincar_reduced[k]}
                        self.logger.debug("\tIncars were different: {} ".format( out_inc))
                    if b_potcar != d_potcar:
                        out_pot = {k:[v, b_potcar[k]] for k,v in d_potcar.items() if v != b_potcar[k]}
                        self.logger.debug("\tPotcar specs were different: {} ".format( out_pot))
                    if not same_cart_positions:
                        self.logger.debug("\tBulk site coords were different")

        #if bulk_task found then take most recently updated bulk_task for defect
        if len( bulk_matched):
            bulk_matched = sorted( bulk_matched, key=lambda x: x["last_updated"], reverse=True)
            out_defect_task["bulk_task"] = bulk_matched[0].copy()
            self.logger.debug("Found {} possible bulk supercell structures. Taking most recent entry updated "
                  "on: {}".format(len(bulk_matched), bulk_matched[0]['last_updated']))
        else:
            self.logger.error("Bulk task doesnt exist for: {}! Cant create defect "
                             "object...\nMetadata: {}\n{}".format( out_defect_task['task_label'],
                                                                   d_potcar, dincar_reduced))
            return None


        # get suitable dielectric calc by matching structures (only lattice
        # constant fitting needed - not supercell) and POTCAR settings...
        diel_task_list = additional_tasks['diel']
        diel_sm = StructureMatcher( ltol=self.ltol, stol=self.stol, angle_tol=self.angle_tol,
            primitive_cell=True, scale=False, attempt_supercell=True, allow_subset=False)
        diel_matched = []
        for diel_task in diel_task_list:
            diel_struct = Structure.from_dict( diel_task['input']['structure'])
            if diel_sm.fit( diel_struct, dstruct_withoutdefect):
                #also match essential POTCAR settings and confirm LVTOT = True and LVHAR = True

                diel_potcar = {'pot_spec': set([potelt["titel"] for potelt in diel_task['input']['potcar_spec']]),
                            'pot_labels': set(diel_task['input']['pseudo_potential']['labels'][:]),
                            'pot_type': diel_task['input']['pseudo_potential']['pot_type'],
                            'functional': diel_task['input']['pseudo_potential']['functional']}
                d_potcar = d_potcar_base.copy() #need to reduce in cases of extrinsic species or reordering
                d_potcar['pot_spec'] = set([d for d in d_potcar['pot_spec'] if d in diel_potcar['pot_spec']])
                d_potcar['pot_labels'] = set([d for d in d_potcar['pot_labels'] if d in diel_potcar['pot_labels']])

                if diel_potcar == d_potcar:
                    diel_matched.append( diel_task.copy())
                else:
                    self.logger.debug("Dielectric structure match was found for {} with {}, "
                                      "but:".format( diel_task['task_label'], out_defect_task['task_label']))
                    out_pot = {k:[v, diel_potcar[k]] for k,v in d_potcar.items() if v != diel_potcar[k]}
                    self.logger.debug("\tPotcar specs were different: {} ".format( out_pot))
            # else:
            #     self.logger.debug("{} ({}) had a structure which did not match {} for use "
            #                       "as a dielectric calculation".format( diel_task['task_label'],
            #                                                          diel_task['task_id'],
            #                                                          out_defect_task['task_label']))

        #if diel_tasks found then take most recently updated bulk_task for defect
        if len( diel_matched):
            diel_matched = sorted( diel_matched, key=lambda x: x["last_updated"], reverse=True)
            diel_dict = {'diel_taskdb_task_id': diel_matched[0]['task_id'],
                         'epsilon_static': diel_matched[0]['output']['epsilon_static'],
                         'epsilon_ionic': diel_matched[0]['output']['epsilon_ionic']}
            out_defect_task["diel_task_meta"] = diel_dict.copy()
            self.logger.debug("Found {} possible dieletric calcs. Taking most recent entry updated "
                  "on: {}".format(len(diel_matched), diel_matched[0]['last_updated']))
        else:
            self.logger.error("Dielectric task doesnt exist for: {}! Cant create defect "
                             "object...\nMetadata for defect: {}\n{}".format( out_defect_task['task_label'],
                                                                   d_potcar, dincar_reduced))
            return None

        # FINALLY consider grabbing extra hybrid BS information...
        # first confirm from the INCAR setting that this defect is NOT an HSE calculation itself...
        if dincar_reduced['LHFCALC'] in [None, False, 'False'] and len(additional_tasks['hsebs']):
            hse_bs_matched = []
            for hse_bs_task in additional_tasks['hsebs']:
                hse_bs_struct = Structure.from_dict(hse_bs_task['input']['structure'])
                if diel_sm.fit(hse_bs_struct, dstruct_withoutdefect): #can use same matching scheme as the dielectric structure matcher
                    hse_bs_matched.append( hse_bs_task.copy())
                else:
                    self.logger.debug("{} ({}) had a structure which did not match {} for use "
                                      "as an HSE BS calculation".format( hse_bs_task['task_label'],
                                                                         hse_bs_task['task_id'],
                                                                         out_defect_task['task_label']))

            if len(hse_bs_matched): #match the lowest CBM  and highest VBM values, keeping track of their task_ids
                hybrid_cbm_data = min([[htask['output']['cbm'], htask['task_id']] for htask in hse_bs_matched])
                hybrid_vbm_data = max([[htask['output']['vbm'], htask['task_id']] for htask in hse_bs_matched])
                hybrid_meta = {'hybrid_cbm': hybrid_cbm_data[0], 'hybrid_CBM_task_id': hybrid_cbm_data[1],
                               'hybrid_vbm': hybrid_vbm_data[0], 'hybrid_VBM_task_id': hybrid_vbm_data[1]}
                out_defect_task["hybrid_bs_meta"] = hybrid_meta.copy()
                self.logger.debug("Found hybrid band structure properties for {}:\n\t{}".format( out_defect_task['task_label'],
                                                                                               hybrid_meta))
            else:
                self.logger.debug("Could NOT find hybrid band structure properties for {} despite "
                                  "there being {} eligible hse calculations".format( out_defect_task['task_label'],
                                                                                     len(additional_tasks['hsebs'])))

        return out_defect_task

    def get_bulk_mpid(self, bulk_structure):
        try:
            with MPRester() as mp:
                # mplist = mp.find_structure(bulk_structure) #had to hack this because this wasnt working??
                tmp_mplist = mp.get_entries_in_chemsys(list(bulk_structure.symbol_set))
            mplist = [ment.entry_id for ment in tmp_mplist if ment.composition.reduced_composition == \
                      bulk_structure.composition.reduced_composition]
            #TODO: this is a hack because find_structure was data intensive. simplify the hack to do less queries...
        except:
            raise ValueError("Error with querying MPRester for {}".format( bulk_structure.composition.reduced_formula))

        mpid_fit_list = []
        for trial_mpid in mplist:
            with MPRester() as mp:
                mpstruct = mp.get_structure_by_material_id(trial_mpid)
            if StructureMatcher(ltol=self.ltol, stol=self.stol, angle_tol=self.angle_tol,
                                primitive_cell=True, scale=False, attempt_supercell=True,
                                allow_subset=False).fit(bulk_structure, mpstruct):
                mpid_fit_list.append( trial_mpid)

        if len(mpid_fit_list) == 1:
            mpid = mpid_fit_list[0]
            self.logger.debug("Single mp-id found for bulk structure:{}.".format( mpid))
        elif len(mpid_fit_list) > 1:
            num_mpid_list = [int(mp.split('-')[1]) for mp in mpid_fit_list]
            num_mpid_list.sort()
            mpid  = 'mp-'+str(num_mpid_list[0])
            self.logger.debug("Multiple mp-ids found for bulk structure:{}\nWill use lowest number mpid "
                  "for bulk band structure = {}.".format(str(mpid_fit_list), mpid))
        else:
            self.logger.debug("Could not find bulk structure in MP database after tying the "
                              "following list:\n{}".format( mplist))
            mpid = None

        return mpid

    def reorder_structure( self, s1, s2):
        """
        This takes s1 and reorders the structure in same site-indexed fashion as s2

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


class DefectThermoBuilder(Builder):
    def __init__(self,
                 defects,
                 defectthermo,
                 query=None,
                 compatibility=DefectCompatibility(),
                 ltol=0.2,
                 stol=0.3,
                 angle_tol=5,
                 **kwargs):
        """
        Creates DefectEntry from vasp task docs

        Args:
            defects (Store): Store of defect entries
            defectthermo (Store): Store of DefectPhaseDiagram documents
            query (dict): dictionary to limit materials to be analyzed
            compatibility (PymatgenCompatability): Compatability module to ensure defect calculations are compatible
            ltol (float): StructureMatcher tuning parameter for matching tasks to materials
            stol (float): StructureMatcher tuning parameter for matching tasks to materials
            angle_tol (float): StructureMatcher tuning parameter for matching tasks to materials
        """

        self.defects = defects
        self.defectthermo = defectthermo
        self.query = query if query else {}
        self.compatibility = compatibility #TODO: how it this going to be used? For cut off parameters?
        self.ltol = ltol
        self.stol = stol
        self.angle_tol = angle_tol
        super().__init__(sources=[defects], targets=[defectthermo], **kwargs)

    def get_items(self):
        self.logger.info("DefectThermo Builder Started")

        # Save timestamp for update operation
        self.time_stamp = datetime.utcnow()

        #get all new Defect Entries since last time DefectThermo was updated...
        q = dict(self.query)
        q.update(self.defects.lu_filter(self.defectthermo))  #TODO: does this work??

        defect_entries = list(self.defects.query(criteria=q))
        self.logger.info("Found {} new defect entries to consider".format( len(defect_entries)))

        #group them based on bulk composition element types and task level metadata info
        grpd_entry_list = {}
        for entry_dict in defect_entries:
            #get bulk symbol set
            # #TODO = create the bulk_sym_set...
            bulk_struct = Structure.from_dict( entry_dict['defect']['bulk_structure'])
            bulk_sym_set = frozenset( bulk_struct.symbol_set)

            if bulk_sym_set not in grpd_entry_list.keys():
                grpd_entry_list[bulk_sym_set] = {}

            #get run metadata info
            run_metadata = entry_dict['parameters']['potcar_summary'].copy()
            run_metadata.update( entry_dict['parameters']['incar_calctype_summary'].copy())
            run_metadata['pot_spec'] = frozenset(run_metadata['pot_spec'])
            run_metadata['pot_labels'] = frozenset(run_metadata['pot_labels'])

            if run_metadata not in grpd_entry_list[bulk_sym_set].keys():
                grpd_entry_list[bulk_sym_set][run_metadata] = []

            for prev_run_metadata in grpd_entry_list[bulk_sym_set].keys():
                if prev_run_metadata == run_metadata:
                    grpd_entry_list[bulk_sym_set][prev_run_metadata].append( entry_dict)

        for bulk_set, metadatadict in grpd_entry_list.items():
            self.logger.info("Symbol set {} categorized".format(bulk_set))
            for metakey, metalist in metadatadict.items():
                self.logger.info("\t {} DefectEntries found for md:\n{}".format(len(metalist), metakey))

        flatten_entry_list = [grpd_entry_set for grpd_entry_by_elts in grpd_entry_list.values()
                              for grpd_entry_set in grpd_entry_by_elts.values()]

        self.logger.info('Found {} new thermo entries to update.'.format(len(flatten_entry_list)))
        for entry_list in flatten_entry_list:
            #TODO: also output the defectthermo object if it already exists
            yield entry_list

    def process_items(self, item):
        #group defect entries into same bulk structure type set

        #see if thermo object already exists (add to it if it already does...)

        #store run_meta_data type (hse / scan /gga) and other relevant metadata worth keeping
        # 'task_level_metadata' key in parameters has keys:
        #     defect_dir_name, bulk_dir_name, bulk_taskdb_task_id, potcar_summary, incar_calctype_summary,
        #     defect_incar, bulk_incar, defect_kpoints, bulk_kpoints, defect_task_last_updated

        pass

    def update_targets(self, items):

        self.logger.info("Updating {} DefectThermo documents".format(len(items)))

        self.defectthermo.update(items, update_lu=True, key='entry_id')

