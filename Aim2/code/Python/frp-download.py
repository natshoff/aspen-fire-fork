
# Imports

import os, sys
import geopandas as gpd

sys.path.insert(0, '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/code/Python/')
from laads_data_download import *

########
# Data #
########

# Fire perimeters w/ aspen cover >= 5%
mtbs = gpd.read_file('data/spatial/mod/mtbs/mtbs_perims_west_w_aspen.gpkg')

