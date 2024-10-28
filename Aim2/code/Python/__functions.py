"""
Helper functions for aspen intensity/severity work
maxwell.cook@colorado.edu
"""

import gc, time
import pandas as pd
import numpy as np
import pytz

from datetime import datetime
from shapely.geometry import box
from shapely.geometry import Polygon, MultiPolygon
from rasterstats import zonal_stats

import warnings
warnings.filterwarnings("ignore")  # suppresses annoying geopandas warning

 
def compute_band_stats(geoms, image_da, id_col, attr=None, stats=None, ztype='categorical'):
    """
    Function to compute band statistics for geometries and a single raster band.
    Args:
        geoms: the geometries for which to calculate zonal statistics
        image_da: categorical raster image array
        id_col: the unique identifier for geometries
        attr: the attribute to calculate (example, 'CBH_') for naming
        stats: statistics to calculate (if continuous input data) as list of strings
        ztype: whether to treat raster data as categorical or continuous
    """
    affine = image_da.rio.transform()
    nodataval = image_da.rio.nodata
    arr = image_da.values

    if ztype == 'categorical':

        if attr is None:
            attr = 'evt'

        zs = zonal_stats(
            vectors=geoms[[id_col, 'geometry']],
            raster=arr,
            affine=affine,
            nodata=nodataval,
            categorical=True,
            all_touched=True,
            geojson_out=True
        )

        # Extract the results (properties)
        stats_df = pd.DataFrame(zs)
        stats_df[id_col] = stats_df['properties'].apply(lambda x: x.get(id_col))
        stats_df['properties'] = stats_df['properties'].apply(
            lambda x: {key: val for key, val in x.items() if key != id_col})
        stats_df['props_list'] = stats_df['properties'].apply(lambda x: list(x.items()))

        # Explode the properties to column
        props = stats_df.explode('props_list').reset_index(drop=True)
        props[[attr,'count']] = pd.DataFrame(props['props_list'].tolist(), index=props.index)

        # Handle NaN values
        props.dropna(subset=[attr], inplace=True)  # handle cases where EVT is NaN

        # Tidy the columns.
        props[attr] = props[attr].astype(int)
        props = props[[id_col,attr,'count']].reset_index(drop=True)

        # Calculate the total pixels and percent cover
        total_pixels = props.groupby(props[id_col])['count'].transform('sum')
        props['total_pixels'] = total_pixels
        props['pct_cover'] = (props['count'] / props['total_pixels']) * 100

        del arr, stats, stats_df  # clean up
        gc.collect()

        return props

    elif ztype == 'continuous':
        # Make sure 'stats' is defined
        if stats is None:
            print("! Please provide list of statistics to calculate !")
            return None
        else:
            zs = zonal_stats(
                vectors=geoms[[id_col, 'geometry']],
                raster=arr,
                affine=affine,
                nodata=nodataval,
                stats=stats,
                categorical=False,
                all_touched=True,
                geojson_out=True
            )

            # Extract the dataframe
            stats_df = pd.DataFrame(zs)
            stats_df[id_col] = stats_df['properties'].apply(lambda x: x.get(id_col))
            for stat in stats:
                stats_df[stat] = stats_df['properties'].apply(lambda x: x.get(stat))
                stats_df.rename(columns={stat: f'{attr}_{stat}'}, inplace=True)

            # Tidy the columns
            cols_to_keep = [id_col] + [f'{attr}_{stat}' for stat in stats]
            stats_df = stats_df[cols_to_keep]

            return stats_df


def create_bounds(gdf, buffer=None, method='bounds'):
    """
    Calculate a bounding rectangle for a given geometry and buffer
    Args:
        gdf: perimeter geometry
        buffer: buffer distance to be applied
        method: one of ['bounds','convex_hull','exact']
    """
    if method == 'bounds':
        geom = gdf.geometry.apply(lambda geo: box(*geo.bounds))
        if buffer is not None:
            geom = geom.buffer(buffer)
        # Apply the new geometry
        gdf_ = gdf.copy()
        gdf_.geometry = geom.geometry.apply(
            lambda geo: Polygon(geo) if geo.geom_type == 'Polygon' else MultiPolygon([geo]))
        return gdf_
    elif method == 'convex_hull':
        gdf_ = gdf.copy()
        gdf_['geometry'] = gdf_.geometry.convex_hull
        if buffer is not None:
            gdf_['geometry'] = gdf_.geometry.buffer(buffer)
        return gdf_
    elif method == 'exact':
        if buffer is not None:
            gdf_ = gdf.geometry.buffer(buffer)
            return gdf_
        else:
            gdf_ = gdf.copy()
            gdf_.geometry = gdf_.geometry.buffer(buffer)
            return gdf_


def convert_datetime(acq_date, acq_time):
    """ Function to convert ACQ_DATE and ACQ_TIME to a datetime object in UTC """
    # Ensure ACQ_TIME is in HHMM format
    acq_time = str(acq_time) # force to string
    if len(acq_time) == 3:
        acq_time = '0' + acq_time
    elif len(acq_time) == 2:
        acq_time = '00' + acq_time
    elif len(acq_time) == 1:
        acq_time = '000' + acq_time

    acq_date_str = acq_date.strftime('%Y-%m-%d')
    dt = datetime.strptime(acq_date_str + acq_time, '%Y-%m-%d%H%M')
    dt_utc = pytz.utc.localize(dt)  # Localize the datetime object to UTC
    return dt_utc


def weighted_variance(values, weights):
    """ Calculate weighted variance. """
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)
    return variance


def find_best_match(perim, neighbors, max_size_diff):
    """
    Identifies the 'best match' based on spatial distance and size difference
    Args:
        perim: the polygon geometries
        neighbors: nearest neighbors identified (see 'find_nearest' function)
        max_size_diff: the maximum allowed size difference (in % difference)
    """
    # Get the fire size from the perimeter data
    perim_size = perim['GIS_ACRES']

    # Initialize best score and match
    best_score = float('inf')
    best_match = None

    for _, point in neighbors.iterrows():
        # Calculate the size difference
        ics_size = point.get('FINAL_ACRES', np.nan)

        if perim_size != 0:
            size_diff = abs((ics_size - perim_size) / perim_size) * 100
        else:
            size_diff = float('inf')  # If perim_size is 0, treat it as infinite difference

        if max_size_diff is not None:
            if size_diff > max_size_diff:
                continue  # Skip if size difference exceeds max_size_diff

        # Calculate the spatial distance (assuming it's precomputed)
        spatial_dist = point.get('spatial_dist', 0)

        # Composite score (you can adjust the weights as needed)
        score = spatial_dist + size_diff

        # Check if this is the best match
        if score < best_score:
            best_score = score
            best_match = point

    return best_match


def find_nearest(perims, points, nn, max_dist=50000, max_size_diff=150):
    """
    Finds the nearest points based on spatial proximity, size, and temporal alignment.
    """

    out_nns = []  # storing the resulting nearest neighbors for each perimeter
    no_matches = []  # to store fires with no matches

    for _, perim in perims.iterrows():
        fire_id = perim['OBJECTID']
        perim_geom = perim.geometry
        fire_year = perim['FIRE_YEAR']

        # Filter incident points to the fire year (filtered once per perimeter)
        inci_points = points[points['START_YEAR'] == fire_year]

        # Early check if no points match the fire year
        if inci_points.empty:
            print(f"No matching points for fire year and fire id: {fire_year} / {fire_id}")
            no_matches.append(perim.to_frame().T)  # Append as DataFrame
            continue

        # Calculate distances from the fire perimeter to the incident points
        distances = inci_points.geometry.apply(lambda x: perim_geom.distance(x))

        # Filter by the maximum distance if provided
        if max_dist is not None:
            inci_points = inci_points[distances <= max_dist]
            distances = distances[distances <= max_dist]

        # Check if there are still points left after filtering
        if inci_points.empty:
            no_matches.append(perim.to_frame().T)  # Convert row to DataFrame and append
            continue

        # Sort by distance and retain the nearest NN points
        nearest_points = inci_points.iloc[distances.argsort()[:nn]].copy()

        # Calculate the best match based on size and distance
        best_ = find_best_match(perim, nearest_points, max_size_diff=max_size_diff)

        if best_ is not None:
            best_['NIFC_ID'] = perim['OBJECTID']
            best_['NIFC_NAME'] = perim['INCIDENT']
            best_['NIFC_ACRES'] = perim['GIS_ACRES']
            out_nns.append(best_.to_frame().T)  # Convert best match to DataFrame before appending

    # Concatenate the no_matches
    print(f"There were [{len(no_matches)}/{len(perims)}] fires with no matches.")
    if len(no_matches) > 0:
        no_matches = pd.concat(no_matches, ignore_index=True)
    else:
        no_matches = pd.DataFrame()

    # Concatenate the matches
    if len(out_nns) > 0:
        out_nns = pd.concat(out_nns, ignore_index=True)
    else:
        out_nns = pd.DataFrame()

    return out_nns, no_matches


def monitor_export(task, timeout=30):
    """ Monitors EE export task """
    while task.active():
        print('Waiting for export to finish..\n\tPatience young padawan.')
        time.sleep(timeout)  # Check every 30 seconds

    # Get the status of the task
    status = task.status()

    # Check if the task failed or succeeded
    if status['state'] == 'COMPLETED':
        print("Export completed successfully !!!!")
    elif status['state'] == 'FAILED':
        print(f"Export failed! Bummer. Reason: {status.get('error_message', 'Unknown error')}")
    else:
        print(f"Export ended with state: {status['state']}")