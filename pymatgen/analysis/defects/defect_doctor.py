"""
A set of classes and functions that are useful for defect workflow management
"""

import os


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

