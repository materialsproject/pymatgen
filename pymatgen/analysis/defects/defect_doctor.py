"""
A set of classes and functions that are useful for defect workflow management
"""

from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.analysis.defects.defect_compatibility import DefectCompatibility
from pymatgen.analysis.defects.defect_builder import TaskDefectBuilder
from pymatgen.analysis.defects.thermodynamics import DefectPhaseDiagram
from pymatgen.analysis.structure_matcher import StructureMatcher

from atomate.vasp.fireworks.core import TransmuterFW

from fireworks import Workflow





def get_fw_from_defect( defect, supercell_size,
                        defect_input_set=None,
                        job_type='normal', db_file='>>db_file<<', vasp_cmd='>>vasp_cmd<<'):

    chgdef_trans = ["DefectTransformation"]
    chgdef_trans_params = [{"scaling_matrix": supercell_size,
                            "defect": defect.copy()}]

    test_bulk_struct = defect.bulk_structure.copy()
    test_bulk_struct.make_supercell( supercell_size)
    num_atoms = len(test_bulk_struct)

    def_tag = "{}:{}_{}_{}_{}atoms".format(test_bulk_struct.composition.reduced_formula,
                                           job_type, defect.name, defect.charge, num_atoms)


    if (job_type == 'normal') and (defect_input_set is None):
        defect_sc = defect.generate_defect_structure(supercell=supercell_size)

        reciprocal_density = 50 if job_type == 'hse' else 100
        kpoints_settings = {"reciprocal_density": reciprocal_density}
        stdrd_defect_incar_settings = {"EDIFF": 0.0001, "EDIFFG": 0.001, "IBRION": 2, "ISMEAR": 0, "SIGMA": 0.05,
                                       "ISPIN": 2, "ISYM": 2, "LVHAR": True, "LVTOT": True, "NSW": 100,
                                       "NELM": 60, "ISIF": 2, "LAECHG": False, "LWAVE": True}
        defect_input_set = MPRelaxSet(defect_sc,
                                      user_incar_settings=stdrd_defect_incar_settings.copy(),
                                      user_kpoints_settings=kpoints_settings,
                                      use_structure_charge=True)

    elif defect_input_set is None:
        raise ValueError("job_type = {} and no defect input set is specified...need to specify input set".format( job_type))


    fw = TransmuterFW(name=def_tag, structure=defect.bulk_structure.copy(),
                      transformations=chgdef_trans,
                      transformation_params=chgdef_trans_params,
                      vasp_input_set=defect_input_set,
                      vasp_cmd=vasp_cmd,
                      copy_vasp_outputs=False,
                      db_file=db_file,
                      job_type=job_type,
                      bandstructure_mode="auto")

    return fw




class DefectResubber(object):
    """
    Takes a defect builder result and creates follow up
    workflow outputs as needed (from atomate)
    for completing the defect thermodynamics desired.

    Does this based on supercell size scaling, delocalization metrics...
    """

    def __init__(self, tasks, defects, lpad):
        self.tasks = tasks
        self.defects = defects
        self.lpad = lpad

    def resubmit(self, base_structure, name_wf):
        fws = []

        # get relevant defect tasks from database
        from pymatgen.core import Structure
        defect_task_list = []
        #get defects and bulk transformations
        poss_defect_tasks = self.tasks.find({'transformations.history.0.@module':
                                                  {'$in': ['pymatgen.transformations.standard_transformations',
                                                           'pymatgen.transformations.defect_transformations']}})
        sm = StructureMatcher(primitive_cell=True, scale=False,
                              attempt_supercell=False, allow_subset=False)
        for def_task in poss_defect_tasks:
            init_struct = Structure.from_dict( def_task['transformations']['history'][0]['input_structure'])
            if sm.fit( base_structure, init_struct):
                defect_task_list.append( def_task)

        #get dielectric
        eltset = [elt.symbol for elt in base_structure.composition.elements]
        eltset.sort()
        poss_diel_tasks = self.tasks.find( {'elements': eltset, 'input.incar.LEPSILON': True})
        for diel_task in poss_diel_tasks:
            init_struct = Structure.from_dict( diel_task['input']['structure'])
            if sm.fit( base_structure, init_struct):
                defect_task_list.append( diel_task)



        # make tasks into defect entry objects and run them through compatibility
        dcompat = DefectCompatibility( plnr_avg_var_tol=0.1, plnr_avg_minmax_tol=0.1,
                             atomic_site_var_tol=0.1, atomic_site_minmax_tol=0.1,
                             tot_relax_tol=1.0, perc_relax_tol=20.,
                             defect_tot_relax_tol=0.1, preferred_cc='freysoldt',
                             free_chg_cutoff=4.,  use_bandfilling=True, use_bandedgeshift=False)
        tdb = TaskDefectBuilder( defect_task_list, gga_run=True, compatibility=dcompat)
        defect_entry_list = tdb.process_entries(self, run_compatibility=True)



        # load to thermodynamics and figure out which additional charge states need to be calculated
        vbm = defect_entry_list[0].parameters['vbm']
        band_gap = defect_entry_list[0].parameters['gga_gap']
        defect_thermo = DefectPhaseDiagram( defect_entry_list, vbm, band_gap, filter_compatible=False)
        print('--> Defects loaded to DefectPhaseDiagram.\n\nFinished charges:')
        for k, v in defect_thermo.finished_charges.items():
            print('\t{}: {}'.format( k, v))
        print('\nSTABLE charges:')
        for k, v in defect_thermo.stable_charges.items():
            print('\t{}: {}\t(t.l.: {} ) '.format( k, v, defect_thermo.transition_levels[k]))

        print('\nNow consider charges for follow up..')
        rec_dict = defect_thermo.suggest_charges()
        for defname, charge_list in rec_dict.items():
            defect_template = defect_thermo.stable_entries[defname][0].copy()
            for charge in charge_list:
                defect = defect_template.copy()
                defect.set_charge( charge)

                for task in defect_task_list:
                    if defect_template.parameters['dir_name'] in task['dir_name']:
                        supercell_size = task['transformations']['history'][0]['scaling_matrix'][:]
                        break

                fw = get_fw_from_defect( defect, supercell_size)
                fws.append( fw)
                print('\trerunning ', defect.name, defect.charge)



        # load Workflow to lpad
        # wf = Workflow( fws, name=name_wf)
        # self.lpad.add_wf( wf)


        return



if __name__ == "__main__":
    #get database
    from atomate.vasp.database import VaspCalcDb
    db = VaspCalcDb.from_db_file("db.json", admin=True)
    tasks = db.get_collections('tasks')
    defects?? #TODO make a defect collection for storing defect entries...

    drs = DefectResubber(tasks, defects, lpad)
    drs.resubmit( bulk_struct, name_wf)



