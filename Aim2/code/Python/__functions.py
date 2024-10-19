"""
Helper functions for aspen intensity/severity work
maxwell.cook@colorado.edu
"""

import pandas as pd
import numpy as np
import gc
import pytz

from datetime import datetime
from shapely.geometry import box
from shapely.geometry import Polygon, MultiPolygon
from rasterstats import zonal_stats

import warnings
warnings.filterwarnings("ignore")  # suppresses annoying geopandas warning


def compute_band_stats(geoms, image_da, id_col):
    """
    Function to compute band statistics for geometries and a single raster band.
    Args:
        geoms: the geometries for which to calculate zonal statistics
        image_da: categorical raster image array
        id_col: the unique identifier for geometries
    """
    affine = image_da.rio.transform()
    nodataval = image_da.rio.nodata
    arr = image_da.values

    stats = zonal_stats(
        vectors=geoms[[id_col, 'geometry']],
        raster=arr,
        affine=affine,
        nodata=nodataval,
        categorical=True,
        all_touched=True,
        geojson_out=True
    )

    # Extract the results (properties)
    stats_df = pd.DataFrame(stats)
    stats_df[id_col] = stats_df['properties'].apply(lambda x: x.get(id_col))
    stats_df['properties'] = stats_df['properties'].apply(
        lambda x: {key: val for key, val in x.items() if key != id_col})
    stats_df['props_list'] = stats_df['properties'].apply(lambda x: list(x.items()))

    # Explode the properties to column
    props = stats_df.explode('props_list').reset_index(drop=True)
    props[['evt', 'count']] = pd.DataFrame(props['props_list'].tolist(), index=props.index)

    # Handle NaN values
    props.dropna(subset=['evt'], inplace=True)  # handle cases where EVT is NaN
    props['evt'] = props['evt'].astype(int)
    props = props[[id_col, 'evt', 'count']].reset_index(drop=True)

    # Calculate the total pixels and percent cover
    total_pixels = props.groupby(props[id_col])['count'].transform('sum')
    props['total_pixels'] = total_pixels
    props['pct_cover'] = (props['count'] / props['total_pixels']) * 100

    del arr, stats, stats_df  # clean up
    gc.collect()

    return props


def create_bounds(gdf, buffer=None, by_bounds=True):
    """ Calculate a bounding rectangle for a given geometry and buffer """
    if by_bounds is True:
        geom = gdf.geometry.apply(lambda geo: box(*geo.bounds))
        if buffer is not None:
            geom = geom.buffer(buffer)
        # Apply the new geometry
        gdf_ = gdf.copy()
        gdf_.geometry = geom.geometry.apply(
            lambda geo: Polygon(geo) if geo.geom_type == 'Polygon' else MultiPolygon([geo]))
        return gdf_
    elif by_bounds is False:
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


def best_match(fired_perim, neighbors, max_date_diff=14, max_size_diff=50):
    fired_perim['ig_date'] = pd.to_datetime(fired_perim['ig_date'])
    ig_date = fired_perim['ig_date']
    perim_size = fired_perim['tot_ar_km2'] * 247.105  # Convert from km^2 to acres

    # Initialize best score and match
    best_score = float('inf')
    best_match = None

    for _, point in neighbors.iterrows():
        # Calculate the date difference
        point['DISCOVERY_DATE'] = pd.to_datetime(point['DISCOVERY_DATE'])
        ics_start_date = point['DISCOVERY_DATE']
        date_diff = abs((ig_date - ics_start_date).days)

        if date_diff > max_date_diff:
            continue  # Skip if date difference exceeds max_date_diff

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
        score = spatial_dist + date_diff + size_diff

        # Check if this is the best match
        if score < best_score:
            best_score = score
            best_match = point

    return best_match


def find_nearest(perims, points, nn=10, max_dist=None, date_range=14, max_size_diff=None):
    """
    Finds the nearest points based on spatial proximity and temporal alignment.
    Args:
        - perims: GeoDataFrame of wildfire perimeters
        - points: GeoDataFrame of incident points (ICS-209-PLUS)
        - NN: the number of neighbors to return
        - max_dist: the maximum distance to search for nearest points
        - date_range: range of days to include in temporal filtering
    Returns: NN nearest points
    """

    # Ensure the same projection
    points = points.to_crs(perims.crs)

    # Convert date columns to datetime
    perims['ig_date'] = pd.to_datetime(perims['ig_date'])
    points['DISCOVERY_DATE'] = pd.to_datetime(points['DISCOVERY_DATE'])

    out_nns = {}  # storing the resulting nearest neighbors for each perimeter
    for _, perim in perims.iterrows():
        fired_id = perim['fired_id']
        centroid = perim.geometry.centroid  # centroid of the fire perimeter

        # Use the entire geometry instead of just the centroid for spatial matching
        perim_geom = perim.geometry

        # Grab the ignition date information
        ig_date = perim['ig_date']

        # Create a date range filter
        date_filter = (points['DISCOVERY_DATE'] >= ig_date - pd.Timedelta(days=date_range)) & \
                      (points['DISCOVERY_DATE'] <= ig_date + pd.Timedelta(days=date_range))
        inci_points = points[date_filter]

        # Check if there are any matches first
        if inci_points.empty:
            print(f"No matching points based on ignition date found for fire {fired_id}")
            out_nns[fired_id] = None
            continue

        # Calculate distances from the fire perimeter to the incident points
        distances = inci_points.geometry.apply(lambda x: perim_geom.distance(x))

        # Filter by the maximum distance if provided
        if max_dist is not None:
            inci_points = inci_points[distances <= max_dist]
            distances = distances[distances <= max_dist]

        # Sort by distance and retain the nearest NN points
        nearest_points = inci_points.iloc[distances.argsort()[:nn]].copy() if not inci_points.empty else None

        if nearest_points is not None:
            best_ = best_match(perim, nearest_points, max_size_diff=max_size_diff)
            out_nns[fired_id] = best_
        else:
            print(f"No matching ICS-209-PLUS points for {fired_id}, skipping ...")
            continue

    return out_nns