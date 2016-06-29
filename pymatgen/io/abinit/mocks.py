# coding: utf-8
"""Mock objects for unit tests."""
from __future__ import division, print_function, unicode_literals

from .nodes import Status
from .tasks import AbinitTask
from .flows import Flow


def change_task_start(task, mocked_status="Error"):
    """Return a AbinitTaskMockedStart object."""
    task.__class__ = AbinitTaskMockedStart
    task.mocked_status = Status.as_status(mocked_status)
    return task


class AbinitTaskMockedStart(AbinitTask):
    """A Task whose status is always self.mocked_status."""
    def start(self, **kwargs):
        self.set_status(self.mocked_status, msg="Mocking status with %s" % self.mocked_status)
        return 1


def infinite_flow(flow):
    """Return an InfiniteFlow."""
    flow.__class__ = InfiniteFlow
    return flow


class InfiniteFlow(Flow):
    """A Flow that will never reach `all_ok`"""
    def check_status(self, **kwargs):
        super(InfiniteFlow, self).check_status(**kwargs)

        for task in self.iflat_tasks(status=self.S_OK):
            task.set_status(task.S_INIT)
            task.reset_from_scratch()
