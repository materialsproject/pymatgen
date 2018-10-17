


def get_level_projection_amounts( defect_wavecar_path, bulk_wavecar_path,
                                  bands_to_project=None, return_bulk_object=False):
    """
    A simple container for using pawpyseed to project bands onto host bands

    :param defect_wavecar_path:
    :param bulk_wavecar_path:
    :param bands_to_project:
    :param return_bulk_object:
    :return:
    """

    if type(bulk_wavecar_path ) == str:
        #MAKE BULK basis thing
        bulk_basis = ???
    else:
        bulk_basis = bulk_wavecar_path


    #GET projection amounts
    proj_amounts = ??

    if return_bulk_object:
        return proj_amounts, bulk_basis
    else:
        return proj_amounts