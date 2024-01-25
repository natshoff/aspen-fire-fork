"""
This script mosaics the results of the image classification into a single composite

The results are the probability of aspen occurrence

maxwell.cook@colorado.edu
"""

# Packages

import os,sys
from rioxarray.merge import merge_arrays
import rasterio.crs

import time
begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

# Globals

dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:5070'  # NAD83 CONUS Albers

# Southern Rocky Mountains Ecoregion (SRME)
srme = gpd.read_file(os.path.join(dataraw,'boundaries/us_eco_l3_srme.gpkg')).to_crs(proj)

############################
# Prep the probability map #
############################

# Load, mosaic, and clip the SRME aspen probability

probs = list_files(os.path.join(datamod,'results/probability/export/prop/'),'*.tif')
out_img = os.path.join(datamod, 'results/probability/s2aspen_prob_10m.tif')

tiles = []
for pr in probs:
    print(os.path.basename(pr))
    tile = rxr.open_rasterio(
        pr,masked=True,cache=False,chunks=True
    ).squeeze().astype(rasterio.uint16)
    tiles.append(tile)
    del tile, pr

# Merge the rasters
print("Merging arrays.")

probs_merge = merge_arrays(
    dataarrays=tiles,
    res=(10, 10),
    crs=rasterio.crs.CRS.from_string(proj),
    nodata=0
)

del tiles

print(f"Merging arrays took {(time.time() - begin)/60} minutes")
begin2 = time.time()

# Write to disk
print("Writing merged raster image ...")

probs_merge.rio.to_raster(
    out_img, tiled=True, windowed=True,
    compress='zstd', zstd_level=9,
    dtype='uint16', driver='GTiff')

del probs_merge

print(f"Writing to disk took {(time.time() - begin2)/60} minutes.")

print(f"Total elapsed time: {(time.time() - begin)/60} minutes.")
