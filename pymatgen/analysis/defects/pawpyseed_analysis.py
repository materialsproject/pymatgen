import numpy as np
from pawpyseed.core.projector import Projector, Wavefunction, make_c_ops, copy_wf, cfunc_call, PAWC


def get_level_projection_amounts(defect_path, bulk_path, spinpol = False, bands_to_project=None,
                                  return_bulk_object=False, unsym = True):
    """
    A simple container for using pawpyseed to project bands onto host bands

    :param defect_path:
    :param bulk_path:
    :param bands_to_project:
    :param return_bulk_object:
    :return:
    """

    mb = True
    if type(bulk_path) == str:
        #MAKE BULK basis thing
        bulk_basis = Wavefunction.from_directory(bulk_path, False)
        mb = False
    else:
        bulk_basis = bulk_path
    defect = Wavefunction.from_directory(defect_path, False)

    projector = Projector(wf = defect, basis = bulk_basis, projector_list = None,
                         unsym_basis = unsym and mb, unsym_wf = unsym, pseudo = False)
    if bands_to_project == None:
        bands_to_project = np.arange(defect.nband)
    #GET projection amounts
    proj_amounts = {}
    for band in bands_to_project:
        proj_amounts[band] = projector.proportion_conduction(band, spinpol)

    projector.wf.free_all()
    projector.free_all()

    if return_bulk_object:
        return proj_amounts, projector.basis
    else:
        projector.basis.free_all()
        return proj_amounts
