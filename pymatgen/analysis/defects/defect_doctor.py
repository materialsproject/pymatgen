"""
A set of classes and functions that are useful for defect workflow management
"""



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
    def __init__(self, file_path, check_for_procar = True, check_for_wavecar = True):
        self.file_path = file_path
        self.check_for_procar = check_for_procar
        self.check_for_wavecar = check_for_wavecar

        self.task_set = None
        self.db_ready = True

        self.initial_atomate_check()

        if self.check_for_procar and self.db_ready:
            # do procar checking routine


        if self.check_for_wavecar and self.db_ready:
            # do wavecar checking routine

    def initial_atomate_check(self):
        """
        run initial atomate drone to see if ready for database

        If this fails, then can run a simple checkup - approach
        which tries to get as much information as possible...
        """


        return

    def get_rerun_info(self):
        """
        If calculation is not ready then this returns information
        for setting up and rerunning atomate workflow...
        """

        if self.db_ready:
            print("This calculation is finished! No rerun is neccessary")
            return
        else:
            # TODO: create this functionality.
            # Note that a similar funcitonality is neccessary for DefectResubber, but should not be identical
            # Want to make this easy for doing a single workflow resubmission of data in a way that is easy to track

            return



class DefectResubber(object):
    """
    Takes a defect builder result and creates follow up
    workflow outputs as needed (from atomate)
    for completing the defect thermodynamics desired
    """

    def __init__(self, var):
        self.var = var

    def test(self):
        return

