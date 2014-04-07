"""
Error handelers of errors originating from the Submission systems.
"""

from __future__ import division

__author__ = "Michiel van Setten"
__copyright__ = " "
__version__ = "0.9"
__maintainer__ = "Michiel van Setten"
__email__ = "mjvansetten@gmail.com"
__date__ = "Mar 24, 2014"

import os
import os.path
import abc
import re

from abc import ABCMeta

from custodian.custodian import Custodian, ErrorHandler


class SchedulerErrorHandler(ErrorHandler):
    """
    Custodian error handler for scheduler related errors
      scheduler takes the schduler, currently: slurm is supported
    """
    def __init__(self, scheduler, err_file='queue.err', out_file='queue.out', run_err_file='run.err', batch_err_file='batch.err'):
        self.scheduler = scheduler
        self.err_file = err_file
        self.out_file = out_file
        self.run_err_file = run_err_file
        self.batch_err_file = batch_err_file
        self.errors = []

    def check(self):
        parser = AbstractErrorParser.factory('slurm', err_file=self.err_file, out_file=self.out_file, run_err_file=self.run_err_file, batch_err_file=self.batch_err_file)
        parser.parse()
        self.errors = parser.errors
        if len(self.errors) == 0:
            return False
        else:
            return True

    def correct(self):
        pass


class AbstractError():
    """
    Error base class
    """
    __metaclass__ = ABCMeta

    def __init__(self, errmsg, meta_data):
        self.errmsg = errmsg
        self.meta_data = meta_data
        self._message = str

    def __str__(self):
        self._message = '%s  %s\n' \
                        '  error message : %s \n' \
                        '  meta data     : %s' % (self.name, self.__doc__, self.errmsg, str(self.meta_data))
        return self._message

    @property
    def name(self):
        return self.__class__.__name__


class SubmitError(AbstractError):
    """
    Errors occurring at submission. The limits on the cluster may have changed.
    """
    def __init__(self, errmsg, meta_data):
        super(SubmitError, self).__init__(errmsg, meta_data)


class FullQueueError(AbstractError):
    """
    Errors occurring at submission. To many jobs in the queue / total cpus / .. .
    """
    def __init__(self, errmsg, meta_data):
        super(FullQueueError, self).__init__(errmsg, meta_data)


class DiskError(AbstractError):
    """
    Errors involving problems writing to disk.
    """
    def __init__(self, errmsg, meta_data):
        super(DiskError, self).__init__(errmsg, meta_data)


class TimeCancelError(AbstractError):
    """
    Error due to exceeding the time limit for the job.
      .limit will return the limit that was broken, None if it could not be determined.
    """
    def __init__(self, errmsg, meta_data):
        super(TimeCancelError, self).__init__(errmsg, meta_data)
        self.limit = None
        if meta_data is not None and 'broken_limit' in meta_data.keys():
                self.limit = meta_data['broken_limit']


class MemoryCancelError(AbstractError):
    """
    Error due to exceeding the memory limit for the job.
      .limit will return the limit that was broken, None if it could not be determined.
    """
    def __init__(self, errmsg, meta_data):
        super(MemoryCancelError, self).__init__(errmsg, meta_data)
        self.limit = None
        if meta_data is not None and 'broken_limit' in meta_data.keys():
                self.limit = meta_data['broken_limit']


class NodeFailureError(AbstractError):
    """
    Error due the hardware failure of a specific node.
     .node will return the problematic node, None if it could not be determined.
    """
    def __init__(self, errmsg, meta_data):
        super(NodeFailureError, self).__init__(errmsg, meta_data)
        self.node = None
        if len(meta_data) > None and 'node' in meta_data.keys():
                self.node = meta_data['node']


class AbstractErrorParser():
    """
    Abstract class for parsing errors originating from the scheduler system and error that are not reported by the
    program itself, i.e. segmentation faults.

    A concrete implementation of this class for a specific sceduler needs a class attribute ERRORS for containing a
    dictionary specifying error:

    ERRORS = {ErrorClass: {'filespecifier' : {'string': "the string to be lookedfor",
                                              'metafilter': "string specifing the regular exprecion to obtian the meta data"}}

    """
    __metaclass__ = ABCMeta

    def __init__(self, err_file, out_file, run_err_file, batch_err_file):
        self.files = {'err': err_file, 'out': out_file, 'run_err': run_err_file, 'batch_err': batch_err_file}
        self.errors = []
        return

    def factory(scheduler, *kwargs):
        if scheduler == "slurm":
            return SlurmErrorParser(*kwargs)
        elif scheduler == "pbs":
            return PBSErrorParse(*kwargs)
        assert 0, "undefined scheduler " + scheduler
    factory = staticmethod(factory)

    @staticmethod
    def extract_metadata(message, metafilter):
        meta_dict = {}
        for key in metafilter.keys():
            print key
            print metafilter[key]
            print metafilter[key][0]
            print " !!!!!!!  next line"
            print message
            match = re.match(metafilter[key][0], message)
            if match is not None:
                print ' !!! match found'
                print match
                meta_dict.update({key: re.match(metafilter[key][0], message).group(metafilter[key][1])})
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
                        message = line.strip()
                        metadata = self.extract_metadata(message, errmsg[k]['metafilter'])
            except (IOError, OSError):
                print self.files[k], not found

        return found, message, metadata

    def parse(self):
        """
        Parse for the occurce of all errors defined in ERRORS
        """
        for my_error in self.ERRORS:
            result = self.parse_single(self.ERRORS[my_error])
            if result[0]:
                self.errors.append(my_error(result[1], result[2]))


class SlurmErrorParser(AbstractErrorParser):
    """
    Implementation for the Slurm scheduler
    """
    def __init__(self, err_file, out_file, run_err_file, batch_err_file):
        super(SlurmErrorParser, self).__init__(err_file, out_file, run_err_file, batch_err_file)

    ERRORS = {SubmitError: {'batch_err': {'string': "sbatch: error: Batch job submission failed:",
                                          'metafilter': ""}},
              FullQueueError: {'batch_err': {'string': "sbatch: error: Batch job submission failed: Job violates accounting/QOS policy",
                                             'metafilter': ""}},
              MemoryCancelError: {'err': {'string': "Exceeded job memory limit",
                                          'metafilter': ""}},
              NodeFailureError: {'run_err': {'string': "can't open /dev/ipath, network down",
                                      'metafilter': {'node': [r"node(\d+)\.(\d+)can't open (\S*), network down \(err=26\)", 1]}}},
              AbstractError: {'out': {'string': "a string to be found",
                                      'metafilter': ""}},
              }


class PBSErrorParse(AbstractErrorParser):
    """
    Implementation for the PBS scheduler
    """


if __name__ == "__main__":
    my_parser = SlurmErrorParser(err_file='queue.err', out_file='queue.out', run_err_file='run.err', batch_err_file='batch.err')
    my_parser.parse()
    print 'parser.errors', my_parser.errors
    for error in my_parser.errors:
        print error
        print error.node