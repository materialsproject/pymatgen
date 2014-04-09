"""
Error handlers for errors originating from the Submission systems.
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Mar 24, 2014"


import re
from abc import ABCMeta


class AbstractCorrection():
    """
    Correction base class
    """
    __metaclass__ = ABCMeta

    def __init__(self, errmsg, meta_data):
        self.errmsg = errmsg
        self.meta_data = meta_data if meta_data is not None else {}

    def __str__(self):
        _message = '%s  %s\n' \
                        '  error message : %s \n' \
                        '  meta data     : %s' % (self.name, self.__doc__, self.errmsg, str(self.meta_data))
        return _message

    @property
    def name(self):
        return self.__class__.__name__


class SubmitCorrection(AbstractCorrection):
    """
    Errors occurring at submission. The limits on the cluster may have changed.
    """


class FullQueueCorrection(AbstractCorrection):
    """
    Errors occurring at submission. To many jobs in the queue / total cpus / .. .
    """


class DiskCorrection(AbstractCorrection):
    """
    Errors involving problems writing to disk.
    """


class TimeCancelCorrection(AbstractCorrection):
    """
    Error due to exceeding the time limit for the job.
      .limit will return the limit that was broken, None if it could not be determined.
    """

    @property
    def limit(self):
        return self.meta_data.get('broken_limit')


class MemoryCancelCorrection(AbstractCorrection):
    """
    Error due to exceeding the memory limit for the job.
      .limit will return the limit that was broken, None if it could not be determined.
    """

    @property
    def limit(self):
        return self.meta_data.get('broken_limit')


class NodeFailureCorrection(AbstractCorrection):
    """
    Error due the hardware failure of a specific node.
     .node will return the problematic node, None if it could not be determined.
    """

    @property
    def node(self):
        return self.meta_data.get('node')


class AbstractCorrecter():
    """
    Abstract class for parsing errors originating from the scheduler system and error that are not reported by the
    program itself, i.e. segmentation faults.

    A concrete implementation of this class for a specific sceduler needs a class attribute ERRORS for containing a
    dictionary specifying error:

    ERRORS = {ErrorClass: {'filespecifier' : {'string': "the string to be lookedfor",
                                              'metafilter': "string specifing the regular exprecion to obtian the meta data"}}

    """
    __metaclass__ = ABCMeta

    def __init__(self, err_file, out_file=None, run_err_file=None, batch_err_file=None):
        self.files = {'err': err_file, 'out': out_file, 'run_err': run_err_file, 'batch_err': batch_err_file}
        self.errors = []
        return

    @staticmethod
    def extract_metadata(lines, metafilter):
        meta_dict = {}
        for key in metafilter.keys():
            values = []
            for line in lines:
                match = re.match(metafilter[key][0], line)
                if match is not None:
                    values.append(re.match(metafilter[key][0], lines).group(metafilter[key][1]))
            values = set(values)
            meta_dict.update({key: values})
        return meta_dict

    def parse_single(self, errmsg):
        """
        Parse the provided files for the corresponding strings.
        """
        found = False
        message = None
        metadata = None
        for k in errmsg.keys():
            print 'parsing ', self.files[k], ' for ', errmsg[k]['string']
            try:
                with open(self.files[k], mode='r') as f:
                    lines = f.read().split('\n')
                for line in lines:
                    if errmsg[k]['string'] in line:
                        found = True
                if found:
                    metadata = self.extract_metadata(lines, errmsg[k]['metafilter'])
            except (IOError, OSError):
                print self.files[k], 'not found'

        return found, message, metadata

    def parse(self):
        """
        Parse for the occurens of all errors defined in ERRORS
        """
        for my_error in self.ERRORS:
            result = self.parse_single(self.ERRORS[my_error])
            if result[0]:
                self.errors.append(my_error(result[1], result[2]))


class SlurmCorrecter(AbstractCorrecter):
    """
    Implementation for the Slurm scheduler
    """


class PBSCorrecter(AbstractCorrecter):
    """
    Implementation for the PBS scheduler
    """

ALL_CORRECTERS = {'slurm': SlurmCorrecter, 'pbs': PBSCorrecter}


def get_correcter(scheduler, code):
    """
    Factory function to provide the correcter for the specified scheduler. If the scheduler is not implemented None is
    returned. The code is passed to the correcter.
    """
    cls = ALL_CORRECTERS.get(scheduler)

    return cls(code)


if __name__ == "__main__":
    my_correcter = get_correcter('slurm', 'ABINIT')
    my_correcter.parse()
    print 'parser.errors', my_correcter.errors
    for error in my_correcter.errors:
        print error