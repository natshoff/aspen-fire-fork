"""
Helper functions for aspen intensity/severity work
maxwell.cook@colorado.edu
"""

import pandas as pd
import rasterio as rio
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import gc

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


def create_bounds(gdf, buffer=None):
    """
    Calculate a bounding rectangle for a given geometry and buffer

    Args:
        gdf: the geometry for which to create bounds
        buffer (optional): buffer distance to apply to bounds
    """
    bounds = gdf.geometry.apply(lambda geom: box(*geom.bounds))
    if buffer is not None:
        bounds = bounds.buffer(buffer)
    # Assign the geometry to the geodataframe
    gdf_ = gdf.copy()
    gdf_.geometry = bounds.geometry.apply(
        lambda geom: Polygon(geom) if geom.geom_type == 'Polygon' else MultiPolygon([geom])
    )
    return gdf_