


def get_projection_amounts( defect_wavecar_path, bulk_wavecar_path):
    """
    A simple container for using pawpyseed to project bands onto host bands

    :param defect_wavecar_path:
    :param bulk_wavcar_path:
    :return:
    """

    if type(bulk_wavcar_path) == str:
        #MAKE BULK basis thing
        bulk_basis = ???
    else:
        bulk_basis = bulk_wavcar_path


    #GET projection amounts
    proj_amounts = ??


    return proj_amounts, bulk_wavecar_path