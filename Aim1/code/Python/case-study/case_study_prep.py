
# Packages
import os,sys
import time
begin = time.time()

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/case-study/'))
from functions import *

# Globals

maindir = os.path.join(os.getcwd(),'data/case-study/')
dataraw = os.path.join(os.getcwd(),'data/spatial/raw/')
datamod = os.path.join(os.getcwd(),'data/spatial/mod/')

proj = 'EPSG:32613'  # project CRS (UTM 13N)

# Read in the probability raster
prob_path = os.path.join(datamod, 'results/probability/srme_skcv_probability_mosaic_prop.tif')
prob_img = rxr.open_rasterio(
    prob_path, masked=True, cache=False, chunks=True, lock=False
).squeeze().astype(rasterio.uint16)

# Load the SRME boundary
srme = os.path.join(datamod,'boundaries/us_eco_l3_srme.gpkg')
# White River National Forest
wrnf = os.path.join(maindir,'wrnf_boundary.gpkg')

# Combine in a list
gdfs = [srme,wrnf]
# Clip to case study ROI
print("Clipping for case study ROIs ...")
for g in gdfs:
    name = os.path.basename(g)[:4]
    print(f"Starting {name} ...")
    gdf = gpd.read_file(g).to_crs(proj).geometry
    # Clip the probability raster
    clipped = prob_img.rio.clip(gdf)
    # Save out the image
    print("Saving case study image ...")
    out_img = os.path.join(maindir,f'test/{name}_aspen_prob_10m.tif')
    clipped.rio.to_raster(
        out_img, tiled=True, lock=threading.Lock(), windowed=True,
        compress='zstd', zstd_level=1, num_threads='all_cpus',
        dtype='uint32', driver='GTiff')
    # Clean up
    del clipped
    del gdf

time.sleep(1)
print(time.time() - begin)