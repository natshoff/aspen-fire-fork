
"""
This script reclassifies the Sentinel-based aspen probability map to a binary distribution map
Reclassification is based on the optimum threshold on probability based on results from the accuracy assessment

maxwell.cook@colorado.edu
"""

# Packages

import os,sys
import time

begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

# Globals

dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:5070'  # NAD83 CONUS Albers

# Load the mosaic probability surface for the case study ROI
prob_path = os.path.join(datamod, 'results/probability/s2aspen_prob_10m.tif')

# Create the 10m binary raster based on the optimum threshold
# Starting with the 10m probability map
thresh = 424  # 'optimum' threshold based on model averages and F1 score
# Read in the probability surface
prob = rxr.open_rasterio(prob_path,masked=True,cache=False).squeeze()

print("Calculating binary resample ...")
# Reclassify to binary based on threshold
bin_out = os.path.join(datamod,'results/classification/s2aspen_prob_10m_binOpt.tif')
bin10 = xr.where(prob >= thresh, 1, 0).astype(rasterio.uint8).rio.reproject(proj)

del prob

print("Exporting raster grid")
bin10.rio.to_raster(
    bin_out,
    compress='zstd', zstd_level=9,
    dtype='uint8', driver='GTiff'
)

print(f"Time elapsed: {(time.time() - begin)/60} minutes.")
