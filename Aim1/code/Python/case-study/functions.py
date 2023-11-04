# Packages

import os, glob
import geopandas as gpd
import rioxarray as rxr
import rasterio
import xarray as xr
import threading


# Functions

def list_files(path, ext):
    return glob.glob(os.path.join(path, '**', '*{}'.format(ext)), recursive=True)


def print_raster(raster,open_file):
    if open_file is True:
        img = rxr.open_rasterio(raster,masked=True).squeeze()
    else:
        img = raster
    print(
        f"shape: {img.rio.shape}\n"
        f"resolution: {img.rio.resolution()}\n"
        f"bounds: {img.rio.bounds()}\n"
        f"sum: {img.sum().item()}\n"
        f"CRS: {img.rio.crs}\n"
        f"NoData: {img.rio.nodata}"
        f"Array: {img}"
    )
    del img


# Clip an image to input shape
def clip_to_shape(img, shp, proj4, buffer=0):
    image = rxr.open_rasterio(
        img, masked=True, cache=False, chunks=True, lock=False
    ).squeeze()
    # Reproject
    image_repr = image.rio.reproject(proj4)
    # read in the roi, match projection, grab geometry
    geom = gpd.read_file(shp).to_crs(crs=proj4).buffer(buffer).geometry
    # clip the image
    return image_repr.rio.clip(geom).fillna(0)


# Function to resample and snap grids
def resample_match_grid(
        in_img, to_img, scale_factor,
        crs, dtype, method='', out_path=os.getcwd()):

    toimg = rxr.open_rasterio(to_img,masked=True,cache=False).squeeze()
    inimg = rxr.open_rasterio(in_img,masked=True,cache=False).squeeze()

    # Define the resample dimensions
    new_height = int(inimg.rio.height * scale_factor)
    new_width = int(inimg.rio.width * scale_factor)

    # Reproject w/ resampling
    resamp = inimg.rio.reproject(
        crs,
        shape=(new_height, new_width),
        resampling=method
    )

    # Reproject match, save out
    out_grid = resamp.rio.reproject_match(toimg)
    out_grid.rio.to_raster(
        out_path, tiled=True, lock=threading.Lock(), windowed=True,
        compress='zstd', zstd_level=1, num_threads='all_cpus',
        dtype=dtype, driver='GTiff'
    )


# Function to reclassify input image to binary
def reclassify_bin(img,to_img,geom,folder,proj4):

    toimg = rxr.open_rasterio(to_img,masked=True).squeeze()
    roi = gpd.read_file(geom).to_crs(crs=proj4).geometry

    in_img = rxr.open_rasterio(img, masked=True).squeeze()
    out_img = xr.where(in_img > 0, 1, 0).astype(rasterio.uint8)
    out_img_match = out_img.rio.reproject_match(toimg)
    out_img_clip = out_img_match.rio.clip(roi)
    out_file = os.path.basename(str(img))[:-4] + '_bin.tif'

    out_img_clip.rio.to_raster(os.path.join(folder, out_file))