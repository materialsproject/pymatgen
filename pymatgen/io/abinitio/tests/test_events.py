# coding: utf-8
from __future__ import unicode_literals, division, print_function

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio import events


#class EventsParserTest(PymatgenTest):
#    report = EventsParser().parse(verbose=1)
#
#    print(report)
#    assert report.run_completed
#    assert report.num_errors, report.num_warnings, report.num_comments == (0, 0, 0)


class EventHandlersTest(PymatgenTest):
    def test_events(self):
        # Autodoc
        events.autodoc_event_handlers()

        for cls in events.get_event_handler_classes():
            # Test pickle
            handler = cls()
            self.serialize_with_pickle(handler, test_eq=False)


if __name__ == '__main__':
    import unittest
    unittest.main()
