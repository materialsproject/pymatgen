# coding: utf-8
"""Mock objects for unit tests."""
from __future__ import division, print_function, unicode_literals

from .tasks import AbinitTask

def mock_task_start(task):
    """
    """
    task.__class__ = ErroredAbinitTask
    return task


class ErroredAbinitTask(AbinitTask): 
    def start(self, **kwargs):
        self.set_status(self.S_ERROR, info_msg="Mocking Error")
        return 1


