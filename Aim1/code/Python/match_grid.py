"""
Match a raster to another raster

maxwell.cook@colorado.edu
"""

# Packages
import os,sys,time
import threading

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

begin = time.time()

# Globals

proj = 'EPSG:5070'

projdir = os.path.join(os.getcwd(),'data/spatial/mod/')

# Define the two regions
rois = ['srme','wrnf']

# Reference grid (binary, matched)
ref_path = os.path.join(projdir,'reference/lc16_evt_200_bin_srme_10m.tif')
# Testing grid (to be matched)
test_path = os.path.join(projdir,'reference/spatial_block_grid_50km2.tif')

################################
# Match the grids if necessary #
################################

test = rxr.open_rasterio(test_path, masked=True, cache=False).squeeze()
ref = rxr.open_rasterio(ref_path, masked=True, cache=False).squeeze()

if test.rio.resolution() == ref.rio.resolution() and \
        test.rio.bounds() == ref.rio.bounds() and \
        test.shape == ref.shape:

    print("Ref and Test match ...")

else:
    print("Mismatch between ref and test ...")

    print(f"Shape of test: {test.shape}\nBounds of ref: {ref.shape}")
    print(f"Resolution of test: {test.rio.resolution()}\nResolution of ref: {ref.rio.resolution()}")
    print(f"Bounds of test: {test.rio.bounds()}\nBounds of ref: {ref.rio.bounds()}")

    print(f"Matching reference image to test image for {os.path.basename(ref_path)}")

    test_ = test.fillna(0).astype(np.uint16)
    img_match = test_.rio.reproject_match(test)

    print(img_match.shape)

    out_path = os.path.basename(test_path)[:-4]+"_match.tif"
    print(out_path)

    img_match.rio.to_raster(
        out_path, compress='zstd', zstd_level=9,
        dtype='uint16', driver='GTiff'
    )

print(time.time() - begin)