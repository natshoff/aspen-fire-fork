
globals().clear()

# Packages

import os,sys
from rasterio.enums import Resampling

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

# Bring in the regional case studies,
# White River National Forest
wrnf = os.path.join(maindir,'wrnf_boundary.gpkg')
# San Juan National Forest
sjnf = os.path.join(maindir,'sjnf_boundary.gpkg')
# Rocky Mountain National Park
rmnp = os.path.join(maindir,'rmnp_boundary.gpkg')

# Combine in a list
gdfs = [wrnf,sjnf,rmnp]


###########################
# Prep the reference images
###########################

print("Starting prep of ref images ...")

# Load the "reference" images of aspen distribution
# USFS TreeMap, ITSP (basal area), LANDFIRE EVT
refs = [
    os.path.join(datamod,'landfire/lc16_evt_srme_aspen_r01_utm.tif'),
    os.path.join(datamod,'USFS/TreeMap2016_clip_spcd746_bin.tif'),
    os.path.join(datamod,'USFS/ITSP_Aspen_30m_srme.tif')
]

match = rxr.open_rasterio(refs[0],masked=True,cache=False).squeeze().rio.reproject(proj)  # for reproject match

# Prep each reference image
i = 0
for r in refs:
    print("Starting: ", os.path.basename(r))
    for gdf in gdfs:
        name = os.path.basename(gdf)[:4]
        outs = [f'{name}_lf_evt_30m.tif', f'{name}_treemap_30m.tif', f'{name}_itsp_30m.tif']

        # Open the reference image
        img = rxr.open_rasterio(r,cache=False).squeeze().rio.reproject(proj)
        rm = img.rio.reproject_match(match)  # Reproject match

        # Clip to the ROI
        print("Clipping ...")
        cl = rm.rio.clip(gdf.geometry)

        # Save out
        print("Saving out ...")

        out_img = os.path.join(maindir, 'ref/', outs[i])

        # Reclassify the ITSP first
        if "ITSP" in os.path.basename(r):
            print("Reclassifying the ITSP ...")
            clr = xr.where(cl > 0, 1, 0).astype(rasterio.uint8)
            print("Reproject and clip ...")
            clr = clr.rio.reproject_match(match)  # need to run this again after reclassify
            clrc = clr.rio.clip(gdf.geometry)
            print("Saving ...")
            clrc.rio.to_raster(
                out_img, tiled=True, lock=threading.Lock(), windowed=True,
                compress='zstd', zstd_level=1, num_threads='all_cpus',
                dtype='uint8', driver='GTiff'
            )
            del clr
            del clrc
        else:
            cl.rio.to_raster(
                out_img, tiled=True, lock=threading.Lock(), windowed=True,
                compress='zstd', zstd_level=1, num_threads='all_cpus',
                dtype='uint8', driver='GTiff'
            )

        del img
        del rm
        del cl

        i = i + 1

    del match

time.sleep(1)
print(time.time() - begin)


###############################################
# Resample / reclassify the probability surface
###############################################

print("Starting prep of the aspen probability surface ...")

# Load the mosaic probability surface for the case study ROI
prob_path = os.path.join(maindir, 'wrnf_aspen_prob_10m.tif')

# Bring in a reference image for 30-meter matching (LANDFIRE)
ref = os.path.join(maindir,'ref/wrnf_lf_evt_30m.tif')

# Generate a 30m probability map matched to LANDFIRE
# (mean & median (probability), max & sum (binary))

# Probability surfaces

print("Calculating mean and median resample ...")

# Average
mean30 = resample_match_grid(
    in_img=prob_path, to_img=ref,
    scale_factor=1/3, crs=proj, method=Resampling.average, dtype="float32",
    out_path=os.path.join(maindir,'wrnf_aspen_prob_30m_avg.tif')
)

# Median
med30 = resample_match_grid(
    in_img=prob_path, to_img=ref,
    scale_factor=1/3, crs=proj, method=Resampling.med, dtype="float32",
    out_path=os.path.join(maindir,'wrnf_aspen_prob_30m_med.tif')
)

del mean30
del med30


# Create the 10m binary raster based on the optimum threshold
# Starting with the 10m probability map
thresh = 434  # 'optimum' threshold based on model averages and F1 score

# Read in the probability surface
prob = rxr.open_rasterio(prob_path,masked=True,cache=False,chunks=True,lock=False).squeeze()

print("Calculating binary resample ...")

# Reclassify to binary based on threshold
bin_out = os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt.tif')
bin10 = xr.where(prob >= thresh, 1, 0).astype(rasterio.uint8).rio.reproject(proj)
bin10 = bin10.rio.clip(gdf.geometry)
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
    out_path=os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt_max30m.tif')
)

# Sum total (how many s2 pixels in one 30m)
sum30 = resample_match_grid(
    in_img=bin_out, to_img=ref, scale_factor=1/3, crs=proj,
    method=Resampling.sum, dtype="uint8",
    out_path=os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt_sum30m.tif')
)
