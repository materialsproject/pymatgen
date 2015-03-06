# coding: utf-8
"""Mock objects for unit tests."""
from __future__ import division, print_function, unicode_literals

from .nodes import Status
from .tasks import AbinitTask


def mock_task_start(task, mocked_status="Error"):
    """
    """
    task.__class__ = AbinitTaskMockedStart
    task.mocked_status = Status.as_status(mocked_status)
    return task


class AbinitTaskMockedStart(AbinitTask): 
    def start(self, **kwargs):
        self.set_status(self.mocked_status, msg="Mocking status with %s" % self.mocked_status)
        return 1
