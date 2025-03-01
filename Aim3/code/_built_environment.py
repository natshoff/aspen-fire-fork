"""
Prepare the WUI grid, COMBUST, population, and structure counts
Author: maxwell.cook@colorado.edu
"""

import os, sys

# Environment variables
maindir = ""

# Functions


def reclassify_bin(img, save_file=True, folder=None):
    """
    :param img:
    :param folder:
    :param save_file: whether the raster should be saved
    :return:
    """
    in_img = rxr.open_rasterio(img,masked=True,lock=False,chunks='auto').squeeze()
    out_img = xr.where(in_img > 0,1,0)

    if save_file is True:
        out_file = os.path.basename(str(img))[:-4]+'_bin.tif'
        out_img.rio.to_raster(os.path.join(folder,out_file))

    return out_img