
"""
Preparation of wildfire footprint data in the western U.S.

"""

# Imports
import pandas as pd
import geopandas as gpd

########
# Data #
########

# MTBS footprints (western US 2017-2021)
mtbs_gpkg = gpd.read_file('data/spatial/raw/MTBS/mtbs_perims_west_2017to2021.shp')
mtbs = pd.DataFrame(mtbs_gpkg.drop(columns='geometry'))

# Load the TreeMap 2016 summary tables (for MTBS perimeters)

# Forest Type
ftype = pd.read_csv('data/tabular/mod/TreeMap/treemap16_MTBS_west_foresttype.csv')
print(ftype.columns.values)
# Tidy the data frame
ftype = ftype[['Event_ID','ForTypName','Count']]  # retain needed columns
ftype['ForTypName'] = ftype['ForTypName'].str.lower()
print(ftype.head())
# Find aspen classes
types = list(ftype['ForTypName'].unique())
aspen_classes = list(filter(lambda x: 'aspen' in x, types))
print(aspen_classes)

# Live Basal Area
balive = pd.read_csv('data/tabular/mod/TreeMap/treemap16_MTBS_west_balive.csv')
print(balive.columns.values)
# Tidy the data frame
balive = balive[['Event_ID','SUM','MEAN','STD']]

# Join to the MTBS footprints
mtbs_treemap = \
    mtbs[['Event_ID','Incid_Name','BurnBndAc','Ig_Date','Ig_Year']]\
    .merge(ftype, on='Event_ID')\
    .dropna(subset='ForTypName')

# Calculate the percent cover for each class per fire
mtbs_treemap['FTypeM2'] = mtbs_treemap['Count'] * 900  # convert to m2
mtbs_treemap['FTypeAcres'] = mtbs_treemap['FTypeM2'] * 0.000247105
mtbs_treemap['FTypePct'] = mtbs_treemap['FTypeAcres'] / mtbs_treemap['BurnBndAc'] * 100
mtbs_treemap = \
    mtbs_treemap[
        ['Event_ID','Incid_Name','Ig_Date','Ig_Year','BurnBndAc',
         'ForTypName','FTypeAcres','FTypePct']]
print(mtbs_treemap.columns.values)

# Subset to fires which burned in aspen forests
mtbs_aspen_fires = mtbs_treemap[mtbs_treemap['ForTypName'] == 'aspen']
print(len(mtbs_aspen_fires[mtbs_aspen_fires['FTypePct'] >= 5.0]))

# Join back to the spatial data, export
mtbs_aspen_fires = mtbs_gpkg[['Event_ID','geometry']].merge(mtbs_aspen_fires, on="Event_ID")
mtbs_aspen_fires_o5pct = mtbs_aspen_fires[mtbs_aspen_fires['FTypePct'] >= 5.0]
mtbs_aspen_fires_o5pct.to_file("data/spatial/mod/MTBS/mtbs_perims_west_w_aspen.gpkg", driver="GPKG")
