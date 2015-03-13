# coding: utf-8
from __future__ import unicode_literals, division, print_function

import os
import datetime

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio import events 

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', "abinitio")

def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


class EventsParserTest(PymatgenTest):
    def test_mgb2_outputs(self):
        """Testing MgB2 output files."""
        # Analyze scf log
        parser = events.EventsParser()
        report = parser.parse(ref_file("mgb2_scf.log"), verbose=1)
        self.assertPMGSONable(report)

        print(report)
        assert (report.num_errors, report.num_warnings, report.num_comments) == (0, 0, 0)
        assert report.run_completed
        fmt = "%a %b %d %H:%M:%S %Y"
        assert report.start_datetime ==  datetime.datetime.strptime("Fri Mar 13 20:08:51 2015", fmt)
        assert report.end_datetime ==  datetime.datetime.strptime("Fri Mar 13 20:08:57 2015", fmt)

        # Analyze nscf log
        report = events.EventsParser().parse(ref_file("mgb2_nscf.log"), verbose=0)
        assert (report.num_errors, report.num_warnings, report.num_comments) == (0, 2, 0)
        self.assertPMGSONable(report)

        #d = report.as_dict()
        #print(d)
        #assert 0

        for warning in report.warnings:
            print(warning)
            # Msonable is conflict with YAMLObject
            #self.assertPMGSONable(warning, check_inst=False)


class EventHandlersTest(PymatgenTest):
    def test_events(self):
        # Autodoc
        events.autodoc_event_handlers()

        for cls in events.get_event_handler_classes():
            # Test pickle
            handler = cls()
            self.serialize_with_pickle(handler, test_eq=False)

        assert events.as_event_class(events.AbinitWarning) == events.AbinitWarning
        assert events.as_event_class('!WARNING') == events.AbinitWarning


if __name__ == '__main__':
    import unittest
    unittest.main()
