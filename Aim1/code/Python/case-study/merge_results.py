
globals().clear()

# Packages

import os,sys
from rioxarray.merge import merge_arrays
import rasterio.crs

import time
begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/case-study/'))
from functions import *

# Globals

maindir = os.path.join(os.getcwd(),'data/spatial/mod/case-study/')
dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:32613'  # project CRS (UTM 13N)


##########################
# Prep the probability map
##########################

# Load, mosaic, and clip the SRME aspen probability

probs = list_files(os.path.join(datamod,'results/probability/export/prop/'),'*.tif')
out_img = os.path.join(datamod, 'results/probability/srme_skcv_probability_p.tif')

tiles = []
for pr in probs:
    print(os.path.basename(pr))
    tile = rxr.open_rasterio(
        pr,masked=True,cache=False,chunks=True,lock=False
    ).squeeze().astype(rasterio.uint64)
    tiles.append(tile)
    del tile

# Merge the rasters
print("Merging arrays ...")

probs_merge = merge_arrays(
    dataarrays=tiles,
    res=(10, 10),
    crs=rasterio.crs.CRS.from_epsg(32613),
    nodata=0
)

time.sleep(1)
print(time.time() - begin)

# Write to disk
print("Writing merged raster image ...")

probs_merge.rio.to_raster(
    out_img, tiled=True, lock=threading.Lock(), windowed=True,
    compress='zstd', zstd_level=1, num_threads='all_cpus',
    dtype='uint64', driver='GTiff')
tiles = []  # clear the list

time.sleep(1)
print(time.time() - begin)

del probs_merge  # clear memory
