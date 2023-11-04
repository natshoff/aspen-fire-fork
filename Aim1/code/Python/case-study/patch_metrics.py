
globals().clear()

# Imports

import os, sys
import pylandstats as pls
from xrspatial.classify import reclassify
from rasterio.enums import Resampling

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/case-study/'))
from functions import *

maindir = os.path.join(os.getcwd(),'data/spatial/mod/case-study/')

proj = 'EPSG:32613'  # project CRS (UTM 13N)

# Data

tiffs = [
    os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt.tif'),  # 10m binary grid
    os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt_max30m.tif'),  # 30m max
    os.path.join(maindir,'ref/wrnf_lf_evt_30m.tif'),  # landfire
    os.path.join(maindir,'ref/wrnf_treemap_30m.tif'),  # treemap
    os.path.join(maindir,'ref/wrnf_itsp_30m.tif')  # itsp
]

# Create 10m grids from the reference images
ref_tiffs = tiffs[2:5]
print(ref_tiffs)

ref_tiffs10m = []
for rt in ref_tiffs:

    name = os.path.basename(rt)[:-8]+"_10m.tif"
    print(name)

    # Resample to 10m and match to Sentinel-based map
    out_tif = resample_match_grid(
        in_img=rt, to_img=tiffs[0], scale_factor=3,
        crs=proj, dtype='uint8', method=Resampling.max,
        out_path=os.path.join(maindir,'ref/{}'.format(name)))

    print_raster(os.path.join(maindir,'ref/{}'.format(name)),open_file=True)
    ref_tiffs10m.append(out_tif)

# Specify a list of metrics to compute

metrics_pch = [
    'area', 'perimeter', 'perimeter_area_ratio',
    'shape_index', 'fractal_dimension', 'euclidean_nearest_neighbor'
]
metrics_cls = [
    'total_area', 'proportion_of_landscape', 'number_of_patches',
    'patch_density', 'largest_patch_index', 'total_edge',
    'edge_density', 'landscape_shape_index', 'effective_mesh_size'
]

# Loop through tiffs and create the dataframes
for tif in tiffs:
    name = os.path.basename(tif)[:-4]
    print(name)

    if "30m" in tif:
        res = (30,30)
    else:
        res = (10,10)
    print(res)

    # Open the image
    img = rxr.open_rasterio(tif,masked=True,cache=False).squeeze()
    # Reclassify 0s
    img_r = reclassify(img, bins=[0, 1], new_values=[1, 2]).astype(int)

    # Initiate the landscape
    print("Converting to landscape object ...")
    ls = pls.Landscape(img_r.values, res=res)

    # Calculate patch metrics
    print("Calculating patch metrics ...")
    patch = ls.compute_patch_metrics_df(metrics=metrics_pch, class_val=2)

    # Calculate class metrics
    print("Calculating class metrics ...")
    clss = ls.compute_class_metrics_df(metrics=metrics_cls)

    # Save out
    patch.to_csv(os.path.join(maindir,'lsmets/wrnf_aspen_{}_patch.csv').format(name))
    clss.to_csv(os.path.join(maindir, 'lsmets/wrnf_aspen_{}_class.csv').format(name))

print("Success ! ...")