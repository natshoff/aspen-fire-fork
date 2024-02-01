
"""
This script calculates the agreement between our Sentinel-based map and three reference datasets
Returns precision, recall and F1 for 1x1 and 3x3 windows

maxwell.cook@colorado.edu
"""

# Packages
import os,sys,time

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/'))
from _functions import *

begin = time.time()

# Globals

proj = 'EPSG:5070'

projdir = os.path.join(os.getcwd(),'data/spatial/mod/')

# Define the two regions
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
    os.path.join(projdir,'reference/usfs_treemap16_balive_int_bin_srme_10m.tif'),
    os.path.join(projdir,'reference/usfs_treemap16_balive_int_bin_wrnf_10m.tif'),
    os.path.join(projdir,'reference/usfs_itsp_aspen_ba_gt10_srme_10m.tif'),
    os.path.join(projdir,'reference/usfs_itsp_aspen_ba_gt10_wrnf_10m.tif')
]

# Load the rasterized spatial grid
blocks_img_path = os.path.join(projdir,'reference/spatial_block_grid_50km2_match.tif')
blocks_img = rxr.open_rasterio(blocks_img_path,masked=True,cache=False).squeeze()
blocks_arr = blocks_img.values.astype(np.uint16)
del blocks_img


#####################################
# Confirm raster grids are matching #
#####################################

# Loop through ROIs
for i in range(len(rois)):
    roi = rois[i]
    print(roi)

    # Sentinel-based map
    test_file_paths = [test for test in tests if str(roi) + ".tif" in test]
    print(os.path.basename(test_file_paths[0]))

    # Reference images
    test = rxr.open_rasterio(test_file_paths[0], cache=False).squeeze().fillna(0)
    print(test.shape)

    ref_file_paths = [ref for ref in refs if str(roi) + "_10m.tif" in ref]
    ref_file_paths.append(blocks_img_path)  # also test the spatial blocks
    print([os.path.basename(ref) for ref in ref_file_paths])

    # Check that they match with the aspen surfaces
    for ref in ref_file_paths:
        print(os.path.basename(ref))

        ref_ = rxr.open_rasterio(ref, cache=False).squeeze()

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
                out_path, compress='zstd', zstd_level=9,
                dtype='uint16', driver='GTiff'
            )

            del img, img_match

        del ref

    del test


################################################
# Workflow to calculate the confusion matrices #
################################################

# Loop regions, perform the analysis
for roi in rois:

    print(f"Starting for {roi}")

    # Load the test image paths for the region
    test_path = [test for test in tests if str(roi) + ".tif" in test][0]
    print(os.path.basename(test_path[0]))
    # Load the reference image paths
    ref_paths = [ref for ref in refs if str(roi) + "_10m.tif" in ref]
    print([os.path.basename(ref) for ref in ref_paths])

    # Open the Sentinel-based map
    test_img = rxr.open_rasterio(test_path,masked=True,cache=False).squeeze().fillna(0)
    test_arr = test_img.values.astype(np.uint16)

    del test_img  # clear up space

    # Loop through reference images
    out_refs = []
    for ref_tif in ref_paths:

        # Open the reference image
        ref_img = rxr.open_rasterio(ref_tif, masked=True, cache=False).squeeze()
        name = os.path.basename(ref_tif)[:-4]  # name of the reference dataset
        print(name)

        blocksizes = [1, 3]  # block sizes (in pixel) used as analytical units.

        # Convert to arrays
        ref_arr = ref_img.values.astype(np.uint16)

        del ref_img  # clean up

        # Check what region (if Southern Rockies, include Block_ID)
        outdata = []
        for blocksize in blocksizes:
            if blocksize > 1:
                arr_ref_res = blockmax(ref_arr, blocksize)
                arr_test_res = blockmax(test_arr, blocksize)
                arr_block_res = blockmax(blocks_arr, blocksize)
            else:
                arr_ref_res = ref_arr
                arr_test_res = test_arr
                arr_block_res = blocks_arr

            # Print the shapes for debugging
            print(
                f"Blocksize {blocksize}: "
                f"Reference - {arr_ref_res.shape}, "
                f"Test - {arr_test_res.shape}, "
                f"Blocks - {arr_block_res.shape}")

            # Check if the reshaped arrays have the same shape
            if arr_ref_res.shape != arr_test_res.shape:
                raise ValueError(
                    f"Reference and test arrays have different shapes: "
                    f"{arr_ref_res.shape} vs {arr_test_res.shape}")
            elif arr_block_res != arr_ref_res:
                raise ValueError(
                    f"Reference and block arrays have different shapes: "
                    f"{arr_block_res.shape} vs {arr_ref_res.shape}")

            print("Creating data frame ...")

            if roi == 'srme':
                currdf = pd.DataFrame({
                    'block': blocks_arr.flatten(),  # if SRME, add the blocks
                    'ref': arr_ref_res.flatten(),
                    'test': arr_test_res.flatten()
                }).query('ref != 0 or test != 0')

                # Calculate confusion matrix by Block_ID
                currdf = currdf.groupby('block').apply(lambda x: pd.Series({
                    'tp': ((x['ref'] == 1) & (x['test'] == 1)).sum(),
                    'fp': ((x['ref'] == 0) & (x['test'] == 1)).sum(),
                    'fn': ((x['ref'] == 1) & (x['test'] == 0)).sum()
                })).reset_index()

                currdf['blocksize'] = blocksize
                outdata.append(currdf)

                del currdf, arr_ref_res, arr_test_res

            else:
                currdf = pd.DataFrame({
                    'ref': arr_ref_res.flatten(),
                    'test': arr_test_res.flatten()
                }).query('ref != 0 or test != 0')

                # Calculate tp, fp, fn directly
                tp = (currdf['ref'] == 1) & (currdf['test'] == 1).sum()
                fp = (currdf['ref'] == 0) & (currdf['test'] == 1).sum()
                fn = (currdf['ref'] == 1) & (currdf['test'] == 0).sum()

                # Create a new DataFrame to store these values
                currdf = pd.DataFrame([{
                    'tp': tp,
                    'fp': fp,
                    'fn': fn,
                    'blocksize': blocksize
                }])

                outdata.append(currdf)

                del currdf, tp, fp, fn, arr_ref_res, arr_test_res

            # Free up some more space
            del arr_ref_res, arr_test_res

        del test_arr, ref_arr, ref_img

        outdatadf = pd.concat(outdata)
        outdatadf = calc_accmeas(outdatadf)
        outdatadf['source'] = name

        # Save the results to a CSV
        outdatadf.to_csv(
            os.path.join(
                os.getcwd(), f'data/tabular/mod/results/accmeas/global_accmeas_multi_blocks_{name}.csv'),
            index=False)
        out_refs.append(outdatadf)

        del outdata, outdatadf

    # Bind the results together for plotting
    outdfs = pd.concat(out_refs).reset_index(drop=True)
    outdfs.to_csv(
        os.path.join(os.getcwd(),f'data/tabular/mod/results/global_accmeas_multi_blocks_full_{roi}.csv'),
        index=False)

    del test_img, outdfs

print("Complete!")

print(time.time() - begin)