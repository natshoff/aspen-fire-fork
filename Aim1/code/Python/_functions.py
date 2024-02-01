# Packages

import os, glob, sys
import geopandas as gpd
import rioxarray as rxr
import rasterio
import xarray as xr
import numpy as np
import pandas as pd
from osgeo import gdal

# Bring in Johannes' accuracy scripts

sys.path.append(
    '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/earth-lab/opp-urban-fuels/local_accuracy-main')
import accmeas


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

    # Define the resample dimensions
    new_height = int(in_img.rio.height * scale_factor)
    new_width = int(in_img.rio.width * scale_factor)

    # Reproject w/ resampling
    resamp = in_img.rio.reproject(
        crs,
        shape=(new_height, new_width),
        resampling=method
    )

    del in_img

    # Reproject match, clip, save out
    out_grid = resamp.rio.reproject_match(to_img)

    del to_img, resamp

    out_grid.rio.to_raster(
        out_path, compress='zstd', zstd_level=9,
        dtype=dtype, driver='GTiff'
    )

    del out_grid


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


def blockmax(inarr, _blocksize_):
    n = _blocksize_  # Height of window
    m = _blocksize_  # Width of window
    modulo = inarr.shape[0] % _blocksize_
    if modulo > 0:
        padby = _blocksize_ - modulo
        inarr_pad = np.pad(inarr, ((0, padby), (0, 0)), mode='constant', constant_values=0)
    else:
        inarr_pad = inarr
    modulo = inarr.shape[1] % _blocksize_
    if modulo > 0:
        padby = _blocksize_ - modulo
        inarr_pad = np.pad(inarr_pad, ((0, 0), (0, padby)), mode='constant', constant_values=0)
    k = int(inarr_pad.shape[0] / n)  # Must divide evenly
    l = int(inarr_pad.shape[1] / m)  # Must divide evenly
    inarr_pad_blockmax = inarr_pad.reshape(k, n, l, m).max(axis=(-1, -3))  # Numpy >= 1.7.1
    return inarr_pad_blockmax


def get_agreement(test_bin_path,ref_bin_path,
                  block_df=None,zones_path=None):

    # Open the test and ref binary arrays
    test_bin_arr = gdal.Open(test_bin_path).ReadAsArray().flatten()
    ref_bin_arr = gdal.Open(ref_bin_path).ReadAsArray().flatten()

    print("Setting blank arrays ...")
    tps = np.zeros(ref_bin_arr.shape)
    fps = np.zeros(ref_bin_arr.shape)
    tns = np.zeros(ref_bin_arr.shape)
    fns = np.zeros(ref_bin_arr.shape)

    print("Assigning confusion matrix ...")
    tps[np.logical_and(ref_bin_arr == 1, test_bin_arr == 1)] = 1
    fps[np.logical_and(ref_bin_arr == 0, test_bin_arr == 1)] = 1
    tns[np.logical_and(ref_bin_arr == 0, test_bin_arr == 0)] = 1
    fns[np.logical_and(ref_bin_arr == 1, test_bin_arr == 0)] = 1

    print("Creating DataFrame ...")
    acc_df = pd.DataFrame()
    acc_df['tp'] = tps.astype(np.int8)
    acc_df['fp'] = fps.astype(np.int8)
    acc_df['tn'] = tns.astype(np.int8)
    acc_df['fn'] = fns.astype(np.int8)

    if zones_path is not None:
        # Open the zones array
        zone_arr = gdal.Open(zones_path).ReadAsArray().flatten()
        # Gather some information
        acc_df['block_id'] = zone_arr
        acc_df = acc_df.dropna()

        acc_df['block_id'] = acc_df['block_id'].apply(str)

        print("Summarizing by spatial block ...")
        aggr_cat_sum_df = acc_df.groupby('block_id')[['tp', 'fp', 'tn', 'fn']].sum().reset_index()
        aggr_cat_sum_df = aggr_cat_sum_df[aggr_cat_sum_df['block_id'] != '-9999']

        acc_df_out = calc_accmeas(aggr_cat_sum_df)

        # Join to the RUCC table
        block_df = block_df[['id','Block_ID']]
        acc_df_out = pd.merge(acc_df_out, block_df, how='inner', on='id')
        print(acc_df_out.head())

        return acc_df_out

    else:
        acc_df = acc_df.dropna()
        acc_df_out = calc_accmeas(acc_df)

        return acc_df_out


def calc_accmeas(df):

    df['recall'] = df.apply(lambda row: accmeas.recall(row.tp, row.tn, row.fp, row.fn), axis=1)
    df['precision'] = df.apply(lambda row: accmeas.precision(row.tp, row.tn, row.fp, row.fn), axis=1)
    df['f1'] = df.apply(lambda row: accmeas.f1(row.tp, row.tn, row.fp, row.fn), axis=1)

    return df