"""
This class takes a list of database task objects
(for a fixed bulk composition) and then returns
Defect Entries which are ready for thermo analysis.
Includes use of Ccmpatibility class if desired.
"""

import os
import numpy as np
from itertools import product
from monty.json import MontyDecoder, MontyEncoder
from monty.shutil import decompress_file

from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen import MPRester
from pymatgen.io.vasp.outputs import Wavecar, Procar, Vasprun

from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.analysis.defects.core import DefectEntry, Interstitial


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
    vrpath = os.path.join( task_dir, 'vasprun.xml.relax2')
    if not os.path.exists(vrpath):
        vrpath += '.gz'
    if not os.path.exists(vrpath):
        vrpath = os.path.join(task_dir, 'vasprun.xml.relax1')
    if not os.path.exists(vrpath):
        vrpath += '.gz'
    if not os.path.exists(vrpath):
        vrpath = os.path.join(task_dir, 'vasprun.xml')
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


def get_wf_data( task_doc):
    """
    :param wf_file_path:
    :param kpt_wgts: listof kpt weights (as far as I can tell, this list is always in same order as parsed kpoints)
    :return:
    """
    #TODO: add Kyle's approach?
    #see if file is zipped and if it is then decompress
    wf_file_path = os.path.join( task_doc['calcs_reversed'][0]['dir_name'],
                               task_doc['calcs_reversed'][0]['output_file_paths']['wavecar'])
    if '.gz' in task_doc['calcs_reversed'][0]['output_file_paths']['wavecar']:
        print('Decompressing wavecar...')
        decompress_file( wf_file_path)
        wf_file_path = wf_file_path[:-3] #this path relabeling only works for .gz labeled files...TODO: fix this

    print('\t--> loading wavefunction data (a little time intensive)')
    #TODO: this does excessive loading..could fix by loading just the bands I want to look at?
    w = Wavecar(wf_file_path, verbose=False) #note this file path must be decompressed to load

    print('\twavecar file loaded')
    #bands to consider: w.nb, w.efermi, w.band_energy
    fermi_band_set = []
    if w.spin == 1:
        for kpt_ind in range( w.nk):
            for band_ind in range( w.nb):
                if not w.band_energy[kpt_ind][band_ind][2]:
                    fermi_band_set.append(band_ind)
                    break
    elif w.spin == 2:
        for spin_ind in range(len(w.band_energy)):
            for kpt_ind in range( w.nk):
                for band_ind in range( w.nb):
                    if not w.band_energy[spin_ind][kpt_ind][band_ind][2]:
                        fermi_band_set.append(band_ind)
                        break
    fermi_band_set = list(set(fermi_band_set))
    sample_bands = np.arange(min(fermi_band_set) - 10, max(fermi_band_set) + 11)

    spin_str = 'False' if w.spin == 1 else 'True'
    print('\tsampling {} kpts and {} bands for spin={}'.format( w.nk, len(sample_bands), spin_str))
    output = [] #first key = spin index, secondkey = kpt index, third key = band index, values = [x,y] for plotting
    #determine conversion factor for nx,ny.nz
    conv = np.array([np.linalg.norm(w.a[ind])/w.ng[ind] for ind in range(3)])
    for spin_ind in range( w.spin):
        print('\t spin_index {}'.format(spin_ind))
        output.append( [])
        for kpt_ind in range( w.nk):
            print('\t\tkpt {}'.format(kpt_ind))
            output[spin_ind].append( {})
            for band_ind in sample_bands:
                """THIS approach can be used for histogram of spatial extent in units of bins"""
                # # coefficients after FFT to real space
                # c = np.fft.fftn(w.fft_mesh(kpt_ind, band_ind, spin=spin_ind, shift=False))
                # # magnitude squared of wavefunction
                # n = np.abs(np.conj(c)*c)
                # #find charge center
                # avg = np.zeros(3)
                # for p in product(range(n.shape[0]), range(n.shape[1]), range(n.shape[2])):
                #     avg += np.array(p)*n[p]
                # chg_center = avg/np.sum(n)
                # #collect data and bin for reduced storage
                # x, y = ([], [])
                # for pind, p in enumerate(product(range(n.shape[0]), range(n.shape[1]), range(n.shape[2]))):
                #     x.append(np.linalg.norm(chg_center - np.array(p)))
                #     y.append(n[p])
                # H, xedges, yedges = np.histogram2d(x,y, bins=(50, 30))
                # output[spin_ind][kpt_ind][band_ind] = [H, xedges, yedges, chg_center]
                """THIS approach can be used for simple radial spatial extent (max y val) with x in Angstrom"""
                # coefficients after FFT to real space
                c = np.fft.fftn(w.fft_mesh(kpt_ind, band_ind, spin=spin_ind, shift=False))
                # magnitude squared of wavefunction
                n = np.abs(np.conj(c)*c)
                n /= np.sum(n) #normalize the single particle wavefunction
                #find charge center
                avg = np.zeros(3)
                for p in product(range(n.shape[0]), range(n.shape[1]), range(n.shape[2])):
                    avg += np.array(p) * n[p]
                chg_center = avg #units of bins
                #collect data, bin and reduce storage to a simple radial plot of wavefunction extent
                x, y = ([], [])
                for p in product(range(n.shape[0]), range(n.shape[1]), range(n.shape[2])):
                    x.append(np.linalg.norm((chg_center - np.array(p)) * conv)) #converted to angstrom distance
                    y.append(n[p])
                H, xedges, yedges = np.histogram2d(x,y, bins=(50, 20))
                reduced_x, reduced_y = ([], []) #only store the max of each x-value
                for xind, xset in enumerate(H):
                    maxy = 0.
                    for yind, yval in enumerate(xset):
                        if yval:
                            maxy = yedges[yind]
                    reduced_x.append( xedges[xind])
                    reduced_y.append( maxy)
                output[spin_ind][kpt_ind][band_ind] = [reduced_x, reduced_y, chg_center]

    task_doc['calcs_reversed'][0]['output']['wf_data']  = output

    return task_doc


def get_procar_data( task_doc, drone):

    procar_files = drone.filter_files( task_doc['calcs_reversed'][0]['dir_name'],
                                       file_pattern="PROCAR")
    path_to_procar = os.path.join( task_doc['calcs_reversed'][0]['dir_name'],
                                   list(procar_files.values())[-1]) #take most recent procar file
    p = Procar(path_to_procar)
    p_data = {'orbital_keys': p.orbitals}

    #find bands for efermi
    efermi = task_doc['calcs_reversed'][0]['output']['efermi']
    fermi_band_set = []
    for spinset in task_doc['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues'].values():
        for bandset_at_kpt in spinset:
            for band_ind, band in enumerate(bandset_at_kpt):
                if abs(band[0] - efermi ) < 0.3:
                    fermi_band_set.append(band_ind)
    fermi_band_set = list(set(fermi_band_set))
    sample_bands = np.arange(min(fermi_band_set) - 10, max(fermi_band_set) + 11)
    p_data['bands_sampling'] = sample_bands

    for spin, spinval in p.data.items(): #NOTE: this sometimes has spin zero contribution?
        spin_dat = []
        for kpt_set in spinval:
            kpt_dat = []
            for band_ind, band in enumerate(kpt_set):
                if band_ind in sample_bands:
                    reduced_band_data = []
                    for atomdata in band:
                        reduce_band_atom_data = {orb: val for orb,val in zip(p.orbitals, atomdata) if val} #removes all zero orbital contributions
                        reduced_band_data.append( reduce_band_atom_data)
                    kpt_dat.append(reduced_band_data)
            spin_dat.append( kpt_dat) #this contains list with atoms, and then atomic projection
        p_data[spin.value] = spin_dat

    task_doc['calcs_reversed'][0]['output']['procar_data'] = p_data
    return task_doc



class TaskDefectBuilder(object):
    """
    This does not have same format as a standard Builder, but it does everything we would want
    a defect builder to do.

    Input is a list of tasks from database...
    """
    def __init__(self, list_of_tasks, compatibility=DefectCompatibility()):
        self.list_of_tasks = list_of_tasks
        self.compatibility = compatibility

    def process_entries(self, run_compatibility=False, gga_run=True):
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
        #TODO: add confirmation that all tasks are derived from same approach / same structure

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

            elif 'DefectTransformation' in task['transformations']['history'][0]['@class']:
                defect_task_list.append(task)

            elif 'SupercellTransformation' == task['transformations']['history'][0]['@class']:
                #process bulk task information
                bulk_sc_structure = task['input']['structure']

                #check if mpid exists for this structure (only need to do once)
                if (mpid is None):
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

                if (mpid == 'DNE') or (not gga_run):
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

            # assume that if not already identified as a defect , then a new
            # hybrid level calculation is for adjusting band edge correction
            #TODO: figure out a better way to see if a hybrid BS caclulation
            #    is being suppled for band edge corrections?
            elif 'history' not in task['transformations'].keys() and 'HFSCREEN' in task['input']['incar'].keys():
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
                        bulk_sc_structure.get_sites_in_sphere(defect.coords, 2, include_index=True), key=lambda x: x[1])
                    defect_frac_sc_coords = bulk_sc_structure[poss_deflist[0][2]].frac_coords
                else:
                    poss_deflist = sorted(
                        initial_defect_structure.get_sites_in_sphere(defect.coords, 2, include_index=True), key=lambda x: x[1])
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


            #Load information for Bandfilling related parsing
            if 'vr_eigenvalue_dict' not in defect_task['calcs_reversed'][0]['output'].keys():
                defect_task = add_vr_eigenvalue_dict( defect_task)

            if 'vr_eigenvalue_dict' in defect_task['calcs_reversed'][0]['output'].keys():
                eigenvalues = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['eigenvalues']
                kpoint_weights = defect_task['calcs_reversed'][0]['output']['vr_eigenvalue_dict']['kpoint_weights']
                parameters.update( {'eigenvalues': eigenvalues,
                                    'kpoint_weights': kpoint_weights} )
            else:
                print('ERR: defect {}_{} does not have eigenvalue data for parsing '
                      'bandfilling.'.format(defect.name, defect.charge))

            #try wavefunction loading for delocalization analysis
            if 'wf_data' not in ['calcs_reversed'][0]['output'].keys():
                try:
                    task_doc = get_wf_data( task_doc)
                except:
                    print(" ERROR IN WAVEFUNCTION DATA LOADING... continuing without this data")

            #try procar loading for delocalization analysis
            if 'procar_data' not in ['calcs_reversed'][0]['output'].keys():
                try:
                    task_doc = get_procar_data( task_doc, self.drone)
                except:
                    print(" ERROR IN PROCAR DATA LOADING... continuing without this data")


            defect_entry = DefectEntry( defect, parameters['defect_energy'] - parameters['bulk_energy'],
                                        corrections = {}, parameters = parameters)

            if run_compatibility:
                defect_entry = self.compatibility.process_entry( defect_entry)


            defect_entry_list.append( defect_entry)


        return defect_entry_list


    def analyze_KS_level(self, task_doc):
        """
        do single particle (KS) level analysis:
        - get band width of bands near the bandgap
        - divy up procar contributions into indivudal bands
        - find spatial extent of KS levels (from wavecar)

        """



