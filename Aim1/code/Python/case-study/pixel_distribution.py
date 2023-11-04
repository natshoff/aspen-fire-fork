
globals().clear()

# Packages

import os,sys
import pandas as pd
import numpy as np

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

# Read in the downsampled grid (sum resampling)

sumgrid = os.path.join(maindir,'wrnf_aspen_prob_10m_binOpt_sum30m.tif')

refgrid1 = os.path.join(maindir,'ref/wrnf_lf_evt_30m.tif')
refgrid2 = os.path.join(maindir,'ref/wrnf_itsp_30m.tif')
refgrid3 = os.path.join(maindir,'ref/wrnf_treemap_30m.tif')

# Create a stack of arrays

a = np.asarray(rxr.open_rasterio(sumgrid,masked=True).squeeze())
b = np.asarray(rxr.open_rasterio(refgrid1,masked=True).squeeze())
c = np.asarray(rxr.open_rasterio(refgrid2,masked=True).squeeze())
d = np.asarray(rxr.open_rasterio(refgrid3,masked=True).squeeze())

arrays = [a,b,c,d]

stack = pd.DataFrame(
    # concatenate column vectors
    np.hstack([
        # first flatten, then convert row vectors to columns
        ar.ravel().reshape(-1, 1)
        # for each array in your list
        for ar in arrays
    ])
)

stack = stack.dropna()
stack = stack.rename(columns={0: 'sum', 1: 'lfevt', 2: 'itsp', 3: 'treemap'})
print(stack.head(20))

stack.to_csv(os.path.join(maindir,'acc/sum_within_out_ref_all.csv'))

