
"""
Calculate the per-block vegetation types and area form the LANDFIRE Existing Vegetation Type (EVT) c. 2016
Southern Rockies ecoregion

"""

import os
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import numpy as np
import random
import rasterio as rio
from shapely.geometry import Point

import time
begin = time.time()

# Load the LF-EVT raster, clipped to spatial blocks
lfevt = rxr.open_rasterio('data/spatial/mod/LANDFIRE/lc16_evt_200_blocks.tif')
# Read in the LANDFIRE lookup table
lookup = pd.read_csv('data/tabular/raw/USFS/LF16_EVT_200.csv')
# Load the spatial block grid
blocks = gpd.read_file('data/spatial/mod/boundaries/spatial_block_grid_50km2.gpkg')
blocks['Block_M2'] = blocks['block_area'] * 1e6  # convert area back to m2
block_ids = list(blocks['grid_id'].unique())  # list of block ids

###########################################################################
# Mask the LF-EVT where aspen live basal area from USFS TreeMap c.2016 > 0

# Generate the LF-EVT summaries by Block

out_csv = os.path.join(os.getcwd(),'data/tabular/mod/training/background/lc16_evt_summary_by_block.csv')

if not os.path.exists(out_csv):
    print("Creating the LF-EVT Data Frame")

    results = []
    for i in range(len(block_ids)):
        block_id = block_ids[i]
        print(f'Processing block {i}')

        # Clip the EVT raster to the current block
        block = blocks[blocks['grid_id'] == block_id]
        block_geom = block.geometry.iloc[0]  # grab the geometry
        clipped = lfevt.rio.clip([block_geom])

        # Count the frequency of each pixel value
        unique, counts = np.unique(clipped.data, return_counts=True)
        freq = dict(zip(unique, counts))
        # Handle NoData value
        freq.pop(-32768, None)

        # Calculate the area for each land cover type
        block_area = block['block_area'].iloc[0]
        pixel_area = 900  # 30-meter spatial resolution

        for value, count in freq.items():
            area = count * pixel_area
            # Find the corresponding land cover type in the lookup table
            lc_row = lookup[lookup['VALUE'] == value].iloc[0]
            # Create the dataframe with EVT attributes
            results.append({
                'Block_ID': block_id,
                'Block_M2': block_area,
                'EVT_NAME': lc_row['EVT_NAME'],
                'EVT_GP_N': lc_row['EVT_GP_N'],
                'EVT_CLASS': lc_row['EVT_CLASS'],
                'EVT_SBCLS': lc_row['EVT_SBCLS'],
                'EVT_AREA_M2': area
            })

        del clipped, block, block_geom, block_area, unique, counts, \
            freq, pixel_area, value, block_id, area, count, lc_row

    del lfevt, lookup

    print(f"Total elapsed time: {(time.time() - begin)/60} minutes.")

    # Create the final data frame
    background = pd.DataFrame(results)
    background.head()

    del results

    # Save the file
    background.to_csv(out_csv)

else:
    # Load the results
    background = pd.read_csv(out_csv).drop('Unnamed: 0', axis=1)

    del lfevt

print(list(background.columns))

# Remove the aspen classes from the dataframe
background = background[~background['EVT_NAME'].str.contains("aspen", case=False, na=False)]


####################################################
# Get the count of presence points within each block

# Load the aspen presence data
presence = gpd.read_file('data/spatial/mod/training/points/gee/pi_points_srme_m500_.shp')
presence = presence.to_crs(blocks.crs)  # Ensure the same CRS

# Join to the block grid
presence = gpd.sjoin(presence,blocks,how="left",predicate="within")
presence = presence.drop(['index_right'], axis=1)
print(presence.columns)

# Get a count of presence points per block
count = presence.groupby('grid_id').size()
count.name = 'n_presence'
count = count.reset_index()
count = count.rename(columns={'grid_id': 'Block_ID'})
print(count.head())

# Join the presence count to the background data frame, handle NA
background = pd.merge(background, count, on='Block_ID', how='left')
background['n_presence'] = background['n_presence'].fillna(0).astype(int)

# Join to the blocks as well
blocks = blocks.rename(columns={'grid_id': 'Block_ID'})
blocks = pd.merge(blocks,count,on='Block_ID',how='left')

del presence

print(background.columns)
print(blocks.columns)

blocks.to_file('data/spatial/mod/boundaries/spatial_block_grid_50km2_count.gpkg')

###########################################################################
# Calculate the number of background samples per block for the three levels
# Proportionally stratified (area) and related to the number of presence points in the block

blocks_f = blocks[blocks['n_presence'] >= 100].copy()  # at least 100 presence samples in the block
blocks_f['n_background'] = blocks_f['n_presence'] * 10  # ten times the number of presence samples (10:1 ratio)

# Remove blocks with < 100 presence samples
background_f = background[background['n_presence'] >= 100].copy()
# Calculate the number of background samples, keep it to a 10:1 ratio for moderate class imbalance
background_f['n_background'] = background_f['n_presence'] * 10  # ten times the number of presence samples (10:1 ratio)
background_f['n_background'] = background_f['n_background'].astype(int)
background_f['n_background'].describe()

# Create a dataframe which is just Block_ID, n_presence, n_background
print(blocks.columns)

cols = ['EVT_NAME','EVT_GP_N','EVT_SBCLS']  # levels

results = []
dfs_totals = []

MIN_SAMPLES_PER_CLASS = 10

for c in cols:
    print(f'Processing background sample distribution for {c}')

    area_nm = c + '_M2'
    prop_nm = c + '_PROP'
    n_samples = 'n_samples_' + c

    # Calculate the total area for each class/level by block
    df = background_f.groupby(['Block_ID', c])['EVT_AREA_M2'].sum().reset_index(name=area_nm)
    print(len(df))
    print(df.head())

    # Re-attach the attributes
    df_m = pd.merge(df, blocks_f[['Block_ID','Block_M2','n_background']],on='Block_ID',how='left')
    print(len(df_m))
    print(df_m.head())

    # Calculate proportions and number of samples
    df_m[prop_nm] = df_m[area_nm] / df_m['Block_M2']
    df_m[n_samples] = df_m[prop_nm] * df_m['n_background']

    # Ensure a minimum sample
    df_m[n_samples] = df_m[n_samples].apply(lambda x: max(x, MIN_SAMPLES_PER_CLASS))

    results.append(pd.DataFrame(df_m))

    pd.DataFrame(df_m).to_csv(f'data/tabular/mod/training/lc16_evt_{c}_by_block.csv')

    df_total = df_m.groupby(c)[n_samples].sum().reset_index()
    dfs_totals.append(pd.DataFrame(df_total))

    pd.DataFrame(df_total).to_csv(f'data/tabular/mod/training/lc16_evt_{c}_totals.csv')

print(results[0])
print(results[1])
print(results[2])

for df in results:
    print(len(df))


#############################################################################
# Generate spatially balanced points within the LANDFIRE EVT Subclass pixels


def is_far_enough(point, points, min_distance):
    return all(point.distance(other) > min_distance for other in points)


# Load the EVT Subclass raster, with aspen areas from TreeMap masked out (values of 0)
sbcls = rxr.open_rasterio('data/spatial/mod/LANDFIRE/lc16_evt_200_blocks_sbcls_mask_.tif')
# Grab the EVT_SBCLS number of samples by block
sbcls_samples = results[2]  # number of samples per subclass per block
# Initiate a geodataframe to store the output
sample_gdfs = gpd.GeoDataFrame(columns=['geometry'], crs=sbcls.rio.crs)

# Loop blocks and generate the random points
block_ids = list(blocks_f['Block_ID'].unique())  # list of block ids
for block_id in block_ids:
    print(f'Processing background sampling for block {block_id}')

    # Clip the EVT raster to the current block
    block = blocks_f[blocks_f['Block_ID'] == block_id]
    block_geom = block.geometry.iloc[0]  # grab the geometry
    clipped = sbcls.rio.clip([block_geom])  # clip the raster
    transform = clipped.rio.transform()

    # Create a grid of pixel centers
    height, width = clipped.squeeze().shape
    x, y = np.meshgrid(np.arange(width), np.arange(height))
    x, y = rio.transform.xy(transform, y, x, offset='center')
    x = np.array(x)
    y = np.array(y)

    # Get the sample data for this block
    block_samples = sbcls_samples[sbcls_samples['Block_ID'] == str(block_id)]

    for _, lc_row in block_samples.iterrows():
        subclass_id = lc_row['EVT_SBCLS']

        if subclass_id == "Deciduous shrubland":
            continue  # this class has too few samples

        N = int(lc_row['n_samples_EVT_SBCLS'])
        N_ = int(N + 100)  # add 100 because we will take a random sample with minimum distance constraint

        # Create a grid mask
        class_mask = np.isin(clipped.data[0], list(lookup[lookup['EVT_SBCLS'] == subclass_id]['VALUE'].unique()))
        # Grab the pixel centers corresponding to the class
        # Use np.where to find indices where class_mask is True
        y_indices, x_indices = np.where(class_mask)
        # Use these indices to get the corresponding x and y values
        x_selected, y_selected = x[y_indices, x_indices], y[y_indices, x_indices]

        min_distance = 100
        max_attempts = 1000
        i = 0

        sample_points = []
        while len(sample_points) < N and i < max_attempts:
            i += 1
            idx = random.randint(0, len(x_selected) - 1)
            candidate = Point(x_selected[idx], y_selected[idx])

            if not sample_points or is_far_enough(candidate, sample_points, min_distance):
                sample_points.append(candidate)

        if i == max_attempts:
            print(f"Reached a maximum for subclass: {subclass_id}")

        gdf = gpd.GeoDataFrame(geometry=sample_points, crs=clipped.rio.crs)
        gdf['EVT_SBCLS'] = subclass_id
        gdf['Block_ID'] = block_id
        gdf['label'] = 0

        sample_gdfs = pd.concat([sample_gdfs, gdf], ignore_index=True)

print(f"Total elapsed time: {(time.time() - begin)/60} minutes.")

print(sample_gdfs)

# Export to geopackage

# Get the unique classes and set a numeric class code for modelling
classes = {code: idx+1 for idx, code in enumerate(sample_gdfs['EVT_SBCLS'].unique())}  # create numeric label
sample_gdfs['SBCLS_CODE'] = sample_gdfs['EVT_SBCLS'].apply(lambda x: classes.get(x, x))

out_file = os.path.join(os.getcwd(),'data/spatial/mod/training/points/background_samples_evt_sbcls.gpkg')
sample_gdfs.to_file(out_file, driver="GPKG")

final_counts = sample_gdfs['EVT_SBCLS'].value_counts()  # check the counts

# Export a summary table
final_counts.to_csv(os.path.join(os.getcwd(),'data/tabular/mod/training/background_samples_evt_sbcls_counts.csv'))
