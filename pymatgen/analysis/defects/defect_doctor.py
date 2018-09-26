"""
A set of classes and functions that are useful for defect workflow management
"""

import os

from atomate.vasp.drones import VaspDrone


class DefectFirstAid(object):
    """
    NOTE: this may not be universally useful for MP defects.
    More for merging new data into existing databases...

    Takes a file path and confirms that it is parsable by atomate (for purposes of defect calculations)

    confirms that:
        - calculation can be loaded with atomate drone
        - Procar exists
        - Wavecar exists

    If these tests fail, then attempts to return a set of files which can
    """
    def __init__(self, file_path, check_for_procar = True,
                 check_for_wavecar = True, runs=None):
        self.file_path = file_path
        self.check_for_procar = check_for_procar
        self.check_for_wavecar = check_for_wavecar

        self.task_set = None
        self.db_ready = True

        self.drone = VaspDrone(runs=runs, parse_dos=True, parse_locpot=True)

        #run drone assimilation
        try:
            self.task_doc = self.drone.assimilate(self.file_path)
        except:
            print('assimilation not possible for ', self.file_path)
            self.db_ready = False

        #check if procar and wavecar exist
        if self.check_for_procar and self.db_ready:
            if 'procar' not in self.task_doc['calcs_reversed'][0]['output_file_paths'].keys():
                print( 'assimilation worked, but no procar exists...')
                self.db_ready = False

        if self.check_for_wavecar and self.db_ready:
            if 'wavecar' not in self.task_doc['calcs_reversed'][0]['output_file_paths'].keys():
                print( 'assimilation worked, but no wavecar exists...')
                self.db_ready = False


class DefectResubber(object):
    """
    Takes a defect builder result and creates follow up
    workflow outputs as needed (from atomate)
    for completing the defect thermodynamics desired.

    Does this based on supercell size scaling, delocalization metrics...
    """

    def __init__(self, db):
        self.db = db

    def test(self):
        return

