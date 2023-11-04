
globals().clear()

# Packages

import os,sys
import time
from rasterio.enums import Resampling

begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/case-study/'))
from functions import *

# Globals

maindir = os.path.join(os.getcwd(),'data/spatial/mod/case-study/')
dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:32613'  # project CRS (UTM 13N)

# Load the mosaic probability surface for the case study ROI
prob_path = os.path.join(datamod, 'results/probability/srme_skcv_probability_mosaic_prop.tif')
# Bring in a reference image for 30-meter matching (LANDFIRE)
ref = os.path.join(datamod,'landfire/lc16_evt_srme_aspen_r01_utm.tif')

# Generate a 30m probability map matched to LANDFIRE
# (mean & median (probability), max & sum (binary))
print("Calculating mean and median resample ...")
# Mean
mean30 = resample_match_grid(
    in_img=prob_path, to_img=ref,
    scale_factor=1/3, crs=proj, method=Resampling.average, dtype="float32",
    out_path=os.path.join(datamod,'results/probability/aspen_prob_30m_mn.tif')
)
# Median
med30 = resample_match_grid(
    in_img=prob_path, to_img=ref,
    scale_factor=1/3, crs=proj, method=Resampling.med, dtype="float32",
    out_path=os.path.join(datamod,'results/probability/aspen_prob_30m_md.tif')
)

del mean30
del med30

# Create the 10m binary raster based on the optimum threshold
# Starting with the 10m probability map
thresh = 424  # 'optimum' threshold based on model averages and F1 score
# Read in the probability surface
prob = rxr.open_rasterio(prob_path,masked=True,cache=False,chunks=True,lock=False).squeeze()

print("Calculating binary resample ...")
# Reclassify to binary based on threshold
bin_out = os.path.join(datamod,'results/classification/aspen_prob_10m_binOpt.tif')
bin10 = xr.where(prob >= thresh, 1, 0).astype(rasterio.uint8).rio.reproject(proj)
bin10.rio.to_raster(
    bin_out, tiled=True, lock=threading.Lock(), windowed=True,
    compress='zstd', zstd_level=1, num_threads='all_cpus',
    dtype='uint8', driver='GTiff'
)

# Create the binary surfaces at 30m

# Maximum within 30m
max30 = resample_match_grid(
    in_img=bin_out, to_img=ref, scale_factor=1/3, crs=proj,
    method=Resampling.max, dtype="uint8",
    out_path=os.path.join(datamod,'results/classification/aspen_prob_10m_binOpt_max30m.tif')
)
# Sum total (how many s2 pixels in one 30m)
sum30 = resample_match_grid(
    in_img=bin_out, to_img=ref, scale_factor=1/3, crs=proj,
    method=Resampling.sum, dtype="uint8",
    out_path=os.path.join(datamod,'results/classification/aspen_prob_10m_binOpt_sum30m.tif')
)

time.sleep(1)
print(time.time() - begin)
