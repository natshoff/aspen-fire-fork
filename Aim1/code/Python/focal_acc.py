
globals().clear()

# Packages

import os,sys,time
import numpy as np
import pandas as pd

# Custom functions
sys.path.append(os.path.join(os.getcwd(),'code/Python/case-study/'))
from _functions import *

begin = time.time()

# Globals

projdir = os.path.join(os.getcwd(),'data/spatial/mod/')

rois = ['srme','wrnf']

# Target grid (resampled 10-meter map using maximum resampling)
tests = [
    os.path.join(projdir,'results/classification/aspen_prob_10m_binOpt_max30m_srme.tif'),  # 30m max
    os.path.join(projdir,'results/classification/aspen_prob_10m_binOpt_max30m_wrnf.tif')
]

# Reference grids (binary, matched)
refs = [
    os.path.join(projdir,'landfire/lc16_evt_srme_aspen_r01_utm_match_srme.tif'),
    os.path.join(projdir,'landfire/lc16_evt_srme_aspen_r01_utm_match_wrnf.tif'),
    os.path.join(projdir,'USFS/treemap16_spcd_746_int_match_srme.tif'),
    os.path.join(projdir,'USFS/treemap16_spcd_746_int_match_wrnf.tif'),
    os.path.join(projdir,'USFS/itsp_aspen_srme_r1__match_srme.tif'),
    os.path.join(projdir, 'USFS/itsp_aspen_srme_r1__match_wrnf.tif')
]

# Check that they match with the aspen surfaces
img1 = rxr.open_rasterio(tests[0],cache=False).squeeze()
img2 = rxr.open_rasterio(refs[0],cache=False).squeeze()
if img1.rio.resolution() == img2.rio.resolution() and \
        img1.rio.bounds() == img2.rio.bounds():
    print("Ref and Test match ...")
else:
    print("Mismatch between ref and test ...")
# Clear the memory
del img1, img2


###########
# Workflow
###########


def blockmax(inarr, blocksize):
    n = blocksize  # Height of window
    m = blocksize  # Width of window
    print(inarr.shape[0])
    modulo = inarr.shape[0] % blocksize
    print(modulo)
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
    print(k,n,l,m)
    # print(inarr_pad)
    inarr_pad_blockmax = inarr_pad.reshape(k, n, l, m).max(axis=(-1, -3))  # Numpy >= 1.7.1
    # print(inarr_pad_blockmax)
    return inarr_pad_blockmax


geog_scales = [400,600,800]  # spatial support in m
blocksizes = [1]  # block size in px. /
# For a block size of 3 or larger, /
# positive instances in any cell within the block are considered true positives.
orig_res = 30  # in m or linear unit of the raster data CRS
stride = 1  # as measured in 0.5 times the geog_scale.
# e.g. stride=2 means we shift a window of geog_scale=1000 by 500m in x,y.
extract = True  # perform extraction
vis = False  # focal precision recall scatterplots.

# Loop through ROIs
for roi in rois:

    print(f"Starting for {roi}")

    # Load the test image
    test_path = [test for test in tests if str(roi) + ".tif" in test]
    print(test_path[0])

    ref_file_paths = [ref for ref in refs if str(roi) + ".tif" in ref]
    print(ref_file_paths)

    for ref in ref_file_paths:

        print(os.path.basename(ref))
        name = os.path.basename(ref)[:-4]

        arr_test = rxr.open_rasterio(test_path[0], masked=True, cache=False).squeeze().values.astype(np.uint16)
        arr_ref = rxr.open_rasterio(ref, masked=True, cache=False).squeeze().values.astype(np.uint16)

        allfocaldata = []
        for geog_scale in geog_scales:
            for blocksize in blocksizes:

                if blocksize > 1:
                    arr_ref_res = blockmax(arr_ref, blocksize)
                    arr_test_res = blockmax(arr_test, blocksize)
                else:
                    arr_ref_res = arr_ref
                    arr_test_res = arr_test

                curr_blocksize_m = blocksize * orig_res
                curr_bocksize_resampled_px = int(geog_scale / float(curr_blocksize_m))

                shift_m = geog_scale / float(stride)
                shift_orig_cells = shift_m / float(orig_res)
                shift_target_cells = int(shift_orig_cells / blocksize)
                shifts = np.arange(0, curr_bocksize_resampled_px, shift_target_cells)
                reftiles = []
                for i in shifts:
                    for j in shifts:
                        tiles_y = np.array_split(arr_ref_res[i:, j:],
                                                 int(arr_ref_res[i:, j:].shape[0] / float(curr_bocksize_resampled_px)),
                                                 axis=0)
                        for tile_y in tiles_y:
                            tiles_x = np.array_split(tile_y, int(tile_y.shape[1] / float(curr_bocksize_resampled_px)),
                                                     axis=1)
                            for tile_x in tiles_x:
                                reftiles.append(tile_x.copy())
                                if np.sum(tile_x) > 0:
                                    print(i, 'ref', max(shifts), geog_scale, blocksize, np.sum(tile_y), np.sum(tile_x))

                testtiles = []
                for i in shifts:
                    for j in shifts:
                        tiles_y = np.array_split(arr_test_res[i:, j:],
                                                 int(arr_test_res[i:, j:].shape[0] / float(curr_bocksize_resampled_px)),
                                                 axis=0)
                        for tile_y in tiles_y:
                            tiles_x = np.array_split(tile_y, int(tile_y.shape[1] / float(curr_bocksize_resampled_px)),
                                                     axis=1)
                            for tile_x in tiles_x:
                                testtiles.append(tile_x.copy())
                                if np.sum(tile_x) > 0:
                                    print(i, 'test', max(shifts), geog_scale, blocksize, np.sum(tile_y), np.sum(tile_x))

                alltiles = list(zip(reftiles, testtiles))

                del reftiles, testtiles

                numtiles = len(alltiles)
                tilecount = 0
                for tilepair in alltiles:
                    if np.nansum(tilepair[0]) == 0 and np.nansum(tilepair[1]) == 0:
                        continue
                    refvec = tilepair[0].flatten()
                    testvec = tilepair[1].flatten()

                    currdf = pd.DataFrame()
                    currdf['ref'] = refvec
                    currdf['test'] = testvec

                    tp = len(currdf[np.logical_and(currdf.ref == 1, currdf.test == 1)])
                    fp = len(currdf[np.logical_and(currdf.ref == 0, currdf.test == 1)])
                    fn = len(currdf[np.logical_and(currdf.ref == 1, currdf.test == 0)])

                    if float(tp+fp) == 0:
                        prec = 0
                    else:
                        prec = np.around(np.divide(tp, float(tp + fp)), decimals=5)

                    if float(tp+fn) == 0:
                        rec = 0
                    else:
                        rec = np.around(np.divide(tp, float(tp + fn)), decimals=5)

                    allfocaldata.append([geog_scale, blocksize, tp, fp, fn, prec, rec])
                    print(tilecount, numtiles, geog_scale, blocksize, tp, fp, fn, prec, rec)
                    tilecount += 1

                    del tilepair, refvec, testvec, currdf

                del alltiles

        del arr_ref, arr_test

        allfocaldatadf = pd.DataFrame(allfocaldata)
        allfocaldatadf.columns = ['geog_scale', 'blocksize', 'tp', 'fp', 'fn', 'prec', 'rec']
        allfocaldatadf.to_csv(
            os.path.join(os.getcwd(),f'data/tabular/mod/results/focal_accmeas_multi_blocks_{name}_bs1.csv'),index=False)

        del allfocaldatadf

        # if vis:
        #     allfocaldatadf = allfocaldatadf
        #     fig,axs = plt.subplots(len(geog_scales),len(blocksizes),sharex='all',sharey='all',figsize=(5,5))
        #     scalecount = 0
        #     for geog_scale in geog_scales:
        #         blocksizecount = 0
        #         for blocksize in blocksizes:
        #             plotdf = allfocaldatadf[allfocaldatadf.geog_scale == geog_scale]
        #             plotdf = plotdf[plotdf.blocksize == blocksize]
        #             plotdf['refbudens'] = np.log(1+(plotdf.tp+plotdf.fn))
        #             ax=axs[scalecount,blocksizecount]
        #             ax.scatter(plotdf.prec.values,plotdf.rec.values,s=2,alpha=0.9,c=plotdf.refbudens.values,cmap='viridis')
        #             ax.set_xlim([0,1])
        #             ax.set_ylim([0,1])
        #             ax.set_yticks([0,0.5,1])
        #             ax.set_yticklabels([0.0,0.5,1.0])
        #             if scalecount == 2:
        #                 ax.set_xlabel('Block size = %s' % blocksize)
        #             if blocksizecount == 0:
        #                 ax.set_ylabel('Support = %s' % geog_scale)
        #             blocksizecount += 1
        #         scalecount += 1
        #     plt.suptitle('Spatially explicit accuracy: Precision (x) vs. Recall (y)\n'
        #                  'for multiple analytical units and spatial support levels')
        #     plt.show()
        #     fig.savefig(os.path.join(maindir,'acc/prec_rec_scat_sensitivity_{}.png'.format(name)), dpi=150)

print("Complete!")
time.sleep(1)
print(time.time() - begin)