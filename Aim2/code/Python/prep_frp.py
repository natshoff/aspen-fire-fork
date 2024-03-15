
"""
Read in the archived FRP data from the SUOMI VIIRS C2 (375m Obs.) as shapefile for the western U.S. (2017-2022)
Downloaded from the LAADS archive: https://firms.modaps.eosdis.nasa.gov/download/

Extract the total and percent area of aspen forest cover within observations
"""

import os
import pandas as pd
import geopandas as gpd

# Load the VIIRS archive shapefile
frp_path = 'data/spatial/raw/VIIRS/DL_FIRE_SC-C2_361441/fire_archive_SV-C2_361441.shp'
frp = gpd.read_file(frp_path)
frp.head()

frp = pd.read_csv("data/tabular/mod/TreeMap/treemap16_west_foresttype_vnp_plot.csv")
print(frp.head())

# Calculate the percent cover for each class per fire
frp['FTypeM2'] = frp['Count'] * 900  # convert to m2
frp['FTypeAcres'] = frp['FTypeM2'] * 0.000247105
frp['FTypePct'] = frp['FTypeM2'] / 140625 * 100
frp = \
    frp[
        ['FID','ForTypName','FTypeAcres','FTypePct']
    ]
print(frp.columns.values)
print(frp.head())

# Save a copy of the table with percentages attached
frp.to_csv("data/tabular/mod/TreeMap/frp_plot_foresttype.csv")

# Attach some data from the TreeMap 2016 (forest type, live basal area, and canopy percent)

...

# Tidy and export for upload to GEE