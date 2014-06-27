
#!/usr/bin/env python

"""Compute the deltafactor for a given pseudopotential."""

from __future__ import division, print_function

__author__ = 'setten'


import os
import sys
# import abipy.data as data
import abipy.abilab as abilab

from abipy.data.runs import AbipyTest, MixinTest
from pseudo_dojo.dojo.deltaworks import DeltaFactory


class DeltaFactorFlowTest(AbipyTest, MixinTest):
    """
    Unit test for the flow defined in this module.
    Users who just want to learn how to use this flow can ignore this section.
    """
    def setUp(self):
        super(DeltaFactorFlowTest, self).setUp()
        self.init_dirs()
        self.flow = build_flow(workdir=self.workdir)


def build_flow(options):
    # Path of the pseudopotential to test.
    #pseudo = data.pseudo("14si.pspnc")
    #pseudo = data.pseudo("Si.GGA_PBE-JTH-paw.xml")
    here = os.path.curdir
    pseudo = os.path.join(here, "pseudo_to_test")
    print(pseudo)


    # Working directory (default is the name of the script with '.py' removed and "run_" replaced by "flow_")
    # workdir = options.workdir
    # if not options.workdir:
    #    #workdir = os.path.basename(__file__).replace(".py", "").replace("run_", "flow_")

    workdir = 'df_run'

    # Instantiate the TaskManager.
    manager = abilab.TaskManager.from_user_config()  # if not options.manager else options.manager

    # Build the workflow for the computation of the deltafactor.
    # The calculation is done with the parameters and the cif files
    # used in the original paper. We only have to specify
    # the cutoff energy ecut (Ha) for the pseudopotential.
    # The workflow will produce a pdf file with the equation of state
    # and a file deltafactor.txt with the final results in the
    # outdir directory DELTAFACTOR/work_0/outdir.
    factory = DeltaFactory()

    #extra = {}

    if options['test']:
        kppa = 50    # this value is for testing purpose.
        ecut = 8
        pawecutdg = ecut * 2

        # Initialize the flow.
        # flow = abilab.AbinitFlow(workdir=workdir, manager=manager, pickle_protocol=0)

        work = factory.work_for_pseudo(pseudo, accuracy="normal", kppa=kppa,
                                       ecut=ecut, pawecutdg=pawecutdg,
                                       toldfe=1.e-8, smearing="fermi_dirac:0.0005")

        work.start()
        work.wait()
        return
    else:
        flow = abilab.AbinitFlow(workdir=workdir, manager=manager, pickle_protocol=0)
        kppa = 6750  # Use this to have the official k-point sampling
        for ecut in [8, 12, 16, 20, 24, 28, 32, 36]:
            pawecutdg = ecut * 2
            work = factory.work_for_pseudo(pseudo, accuracy="normal", kppa=kppa,
                                           ecut=ecut, pawecutdg=pawecutdg,
                                           toldfe=1.e-8, smearing="fermi_dirac:0.0005")
            # Register the workflow.
            flow.register_work(work, workdir='W'+str(ecut))
        workflow = flow.allocate()
        workflow.build_and_pickle_dump()
        return


#abilab.flow_main
def main(options):
    print(options)
    build_flow(options)


if __name__ == "__main__":
    my_options = {'test': False}

    for arg in sys.argv:
        my_options.update({arg: True})

    print(my_options)
    main(options=my_options)
