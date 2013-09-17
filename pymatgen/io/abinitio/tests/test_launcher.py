#!/usr/bin/env python
from __future__ import division, print_function

import unittest

from pymatgen.io.abinitio.launcher import ScriptEditor

class ScriptEditorTest(unittest.TestCase):

    def test_base(self):
        "base test"
        se = ScriptEditor()
        se.shebang()
        se.declare_var("FOO", "BAR")
        se.add_emptyline()
        se.add_comment("This is a comment")
        se.declare_vars({"FOO1": "BAR1"})
        se.load_modules(["module1", "module2"])
        print(se.get_script_str())



if __name__ == '__main__':
    unittest.main()
