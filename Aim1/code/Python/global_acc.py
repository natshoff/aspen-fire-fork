
"""
This script calculates the agreement between our Sentinel-based map and three reference datasets
Returns precision, recall and F1 for 1x1 and 3x3 windows

maxwell.cook@colorado.edu
"""

# Packages
import os,sys,time
import numpy as np
import pandas as pd

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

begin = time.time()

# Globals

projdir = os.path.join(os.getcwd(),'data/spatial/mod/')

rois = ['srme','wrnf']

# Target grid (resampled 10-meter map using maximum resampling)
tests = [
    os.path.join(projdir,'results/classification/s2aspen_prob_10m_binOpt_srme.tif'),
    os.path.join(projdir,'results/classification/s2aspen_prob_10m_binOpt_wrnf.tif')
]

# Reference grids (binary, matched)
refs = [
    os.path.join(projdir,'reference/lc16_evt_200_bin_srme_10m.tif'),
    os.path.join(projdir,'reference/lc16_evt_200_bin_wrnf_10m.tif'),
    os.path.join(projdir,'reference/usfs_treemap16_bin_srme_10m.tif'),
    os.path.join(projdir,'reference/usfs_treemap16_bin_wrnf_10m.tif'),
    os.path.join(projdir,'reference/usfs_itsp_aspen_ba_gt10_srme_10m.tif'),
    os.path.join(projdir,'reference/usfs_itsp_aspen_ba_gt10_wrnf_10m.tif')
]

# Loop through ROIs
for i in range(len(rois)):

    roi = rois[i]

    test_file_paths = [test for test in tests if str(roi)+".tif" in test]
    print(test_file_paths[0])
    test = rxr.open_rasterio(test_file_paths[0], cache=False).squeeze()

    ref_file_paths = [ref for ref in refs if str(roi)+"_10m.tif" in ref]
    print(ref_file_paths)

    # Check that they match with the aspen surfaces

    for ref in ref_file_paths:
        print(os.path.basename(ref))

        ref_ = rxr.open_rasterio(ref_file_paths[0], cache=False).squeeze()

        if test.rio.resolution() == ref_.rio.resolution() and \
                test.rio.bounds() == ref_.rio.bounds() and \
                test.shape == ref_.shape:

            print("Ref and Test match ...")

            del ref_

        else:
            print("Mismatch between ref and test ...")

            print(f"Shape of test: {test.shape}\nBounds of ref: {ref_.shape}")
            print(f"Resolution of test: {test.rio.resolution()}\nResolution of ref: {ref_.rio.resolution()}")
            print(f"Bounds of test: {test.rio.bounds()}\nBounds of ref: {ref_.rio.bounds()}")

            del ref_

            print(f"Matching reference image to test image for {os.path.basename(ref)}")
            img = rxr.open_rasterio(ref,masked=True,cache=False).squeeze()
            img = img.fillna(0).astype(np.uint16)
            img_match = img.rio.reproject_match(test)
            out_path = ref[:-4]+".tif"
            print(out_path)
            img_match.rio.to_raster(
                out_path, tiled=True, lock=threading.Lock(), windowed=True,
                compress='zstd', zstd_level=9, num_threads='all_cpus',
                dtype='uint16', driver='GTiff'
            )

            del img, img_match

    del test, ref


###########
# Workflow
###########

def blockmax(inarr, blocksize):
    n = blocksize  # Height of window
    m = blocksize  # Width of window
    modulo = inarr.shape[0] % blocksize
    if modulo > 0:
        padby = blocksize - modulo
        inarr_pad = np.pad(inarr, ((0, padby), (0, 0)), mode='constant', constant_values=0)
    else:
        inarr_pad = inarr
    modulo = inarr.shape[1] % blocksize
    if modulo > 0:
        padby = blocksize - modulo
        inarr_pad = np.pad(inarr_pad, ((0, 0), (0, padby)), mode='constant', constant_values=0)
    k = int(inarr_pad.shape[0] / n)  # Must divide evenly
    l = int(inarr_pad.shape[1] / m)  # Must divide evenly
    inarr_pad_blockmax = inarr_pad.reshape(k, n, l, m).max(axis=(-1, -3))  # Numpy >= 1.7.1
    return inarr_pad_blockmax


for roi in rois:

    print(f"Starting for {roi}")

    # Load the test image
    test_path = [test for test in tests if str(roi)+".tif" in test]
    print(test_path[0])

    file_paths = [ref for ref in refs if str(roi)+"_10m.tif" in ref]
    print(file_paths)

    # Loop through reference images
    out_refs = []
    for ref_tif in file_paths:

        test_arr = rxr.open_rasterio(test_path[0], masked=True, cache=False).squeeze().values.astype(np.uint16)

        name = os.path.basename(ref_tif)[:-4]
        print(name)

        blocksizes = [1, 3]  # block sizes (in pixel) used as analytical units.

        outdata = []
        for blocksize in blocksizes:

            ref_arr = rxr.open_rasterio(ref_tif, masked=True, cache=False).squeeze().values.astype(np.uint16)

            if blocksize > 1:
                arr_ref_res = blockmax(ref_arr, blocksize)
                arr_test_res = blockmax(test_arr, blocksize)
            else:
                arr_ref_res = ref_arr
                arr_test_res = test_arr

            # Free up some space
            del ref_arr

            # Print the shapes for debugging
            print(
                f"Blocksize {blocksize}: Reference - {arr_ref_res.shape}, Test - {arr_test_res.shape}")

            print("Creating data frame")

            # Check if the reshaped arrays have the same shape
            if arr_ref_res.shape != arr_test_res.shape:
                raise ValueError(
                    f"Reference and test arrays have different shapes: {arr_ref_res.shape} vs {arr_test_res.shape}")

            currdf = pd.DataFrame({
                'ref': arr_ref_res.flatten(),
                'test': arr_test_res.flatten()
            })

            # Free up some more space
            del arr_ref_res, arr_test_res

            currdf = currdf[-np.logical_and(currdf.ref == 0, currdf.test == 0)]
            tp = len(currdf[np.logical_and(currdf.ref == 1, currdf.test == 1)])
            fp = len(currdf[np.logical_and(currdf.ref == 0, currdf.test == 1)])
            fn = len(currdf[np.logical_and(currdf.ref == 1, currdf.test == 0)])
            print(blocksize, tp, fp, fn)
            outdata.append([blocksize, tp, fp, fn])

        del test_arr

        outdatadf = pd.DataFrame(outdata, columns=['blocksize', 'tp', 'fp', 'fn'])
        outdatadf['prec'] = outdatadf.tp / (outdatadf.tp + outdatadf.fp).astype(np.float64)
        outdatadf['rec'] = outdatadf.tp / (outdatadf.tp + outdatadf.fn).astype(np.float64)
        outdatadf['source'] = name
        outdatadf.to_csv(
            os.path.join(
                os.getcwd(),f'data/tabular/mod/results/accmeas/global_accmeas_multi_blocks_{name}.csv'),index=False)
        out_refs.append(outdatadf)
        outdata = []

    # Bind the results together for plotting
    outdfs = pd.concat(out_refs).reset_index(drop=True)
    outdfs.to_csv(os.path.join(os.getcwd(),f'data/tabular/mod/results/global_accmeas_multi_blocks_full_{roi}.csv'),
                  index=False)

print("Complete!")

print(time.time() - begin)