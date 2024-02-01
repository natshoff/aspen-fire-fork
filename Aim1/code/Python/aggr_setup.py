
"""
This script preps the reference and testing images for the agreement assessment

    - Resamples existing products to common 10m resolution matching the Sentinel-based aspen map
    - Creates a rasterized grid from the spatial blocks, matched to test image

LANDFIRE Existing Vegetation Type (EVT, c. 2016 Remap)
USFS National Individual Tree Species Atlas (c. 2014)
USFS TreeMap (c. 2016)

maxwell.cook@colorado.edu
"""

# Packages

import os,sys
from rasterio.enums import Resampling

import time
begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

# Globals

dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:5070'  # NAD83 CONUS Albers

# Load the study area boundaries
# Southern Rocky Mountains Ecoregion (SRME)
srme = gpd.read_file(os.path.join(dataraw,'boundaries/us_eco_l3_srme.gpkg')).to_crs(proj)
# White River National Forest
wrnf = gpd.read_file(os.path.join(dataraw,'boundaries/wrnf_boundary_srme.gpkg')).to_crs(proj)
# Combine in a list
gdfs = [srme,wrnf]
names = ["srme","wrnf"]

#############################
# Prep the reference images #
#############################

print("Starting prep of reference images ...")

# Target grid (Sentinel-based aspen distribution map)
tests = [
    os.path.join(datamod,'results/classification/s2aspen_prob_10m_binOpt_srme.tif'),
    os.path.join(datamod,'results/classification/s2aspen_prob_10m_binOpt_wrnf.tif')
]

# Load the "reference" images of aspen distribution
# USFS TreeMap, ITSP (basal area), LANDFIRE EVT
refs = [
    os.path.join(datamod,'LANDFIRE/lc16_evt_200_bin.tif'),
    os.path.join(datamod,'USFS/usfs_itsp_aspen_ba_gt10.tif'),
    os.path.join(datamod,'USFS/usfs_treemap16_balive_int_bin.tif')
]

# Loop the reference images, exporting a matched 10-meter grid

for r in refs:
    name = os.path.basename(r)[:-4]
    print(f"Processing: {name}")

    # Open the raster file
    ref = rxr.open_rasterio(r,masked=True,cache=False).squeeze()

    for i in range(0,len(names)):
        print(f"Clipping to the '{names[i]}'")

        gdf = gdfs[i]

        # Open the test image to match
        test = rxr.open_rasterio(tests[i],masked=True,cache=False)

        # Clip to study region
        cl = ref.rio.clip(gdf.geometry)

        # Upscale to 10m and match to the test image
        print("Upsampling to 10m")

        out = os.path.join(datamod, f'reference/{name}_10m.tif')

        ups = resample_match_grid(
            in_img=cl, to_img=test,
            scale_factor=3, crs=proj, method=Resampling.nearest,
            dtype="uint8", out_path=out
        )

        del test

    del ref

print(f"Time elapsed: {(time.time() - begin) / 60} minutes.")
