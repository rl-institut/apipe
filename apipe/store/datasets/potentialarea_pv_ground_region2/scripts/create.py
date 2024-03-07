import os

import geopandas as gpd
import numpy as np
import rasterio
from rasterstats import zonal_stats

from apipe.scripts.geo import (
    convert_to_multipolygon,
    overlay,
    raster_to_vector,
    write_geofile,
)


def process_dataset(
    clipped_raster_path,
    region_muns_path,
    output_gpkg_path,
    area_threshold,
    raster_value_threshold,
):
    """Vectorize and calc raster mean value per polygon

    Parameters
    ----------
    clipped_raster_path : pathlib.Path
        Path to raster file with data
    region_muns_path : pathlib.Path
        Path to regions' municipalities
    output_gpkg_path : pathlib.Path
        Path to output vector file
    area_threshold : float
        Threshold for min. area to be vectorized (in ha/cells)
    raster_value_threshold : float
        Threshold for raster value to be vectorized (0..1)

    Returns
    -------
    None
    """

    # Step 1: Set raster values > 0 to 1 and save
    with rasterio.open(clipped_raster_path) as src:
        raster = src.read(1)
        binary_raster = np.where(raster > raster_value_threshold, 1, 0)
        meta = src.meta
        meta["dtype"] = "int32"
        binary_raster_path = clipped_raster_path.replace(".tif", "_binary.tif")
        with rasterio.open(binary_raster_path, "w", **meta) as dst:
            dst.write(binary_raster, 1)

    # Step 2: Overlay with municipal boundaries and rename
    municipalities = gpd.read_file(region_muns_path)
    gdf = raster_to_vector(binary_raster_path)
    gdf_overlayed = overlay(gdf, municipalities, {"id": "municipality_id"})

    # Calculate the area for each polygon and add as a new column 'area'
    gdf_overlayed["area"] = gdf_overlayed.geometry.area
    # Step 3: Merge adjacent polygons if area meets the threshold
    gdf_overlayed = convert_to_multipolygon(gdf_overlayed)
    gdf_overlayed = gdf_overlayed[
        gdf_overlayed["area"] >= (area_threshold * 10000)
    ]

    # Step 4: Calculate the mean of the original raster values per polygon
    stats = zonal_stats(
        gdf_overlayed,
        clipped_raster_path,
        stats="mean",
        nodata=0,
        all_touched=True,
        geojson_out=False,
    )
    # Add the calculated means as a new column to the GeoDataFrame
    mean_values = [stat["mean"] for stat in stats]
    gdf_overlayed["mean_raster_value"] = mean_values

    # Step 5: Save the result
    write_geofile(
        gdf_overlayed[
            ["municipality_id", "mean_raster_value", "geometry", "area"]
        ],
        output_gpkg_path,
        driver="GPKG",
    )

    # Remove binary raster temp-files
    os.remove(binary_raster_path)


process_dataset(
    clipped_raster_path=snakemake.input.clipped_raster,
    region_muns_path=snakemake.input.region_muns,
    output_gpkg_path=snakemake.output.vector_overlay_gpkg,
    area_threshold=float(snakemake.params.area_threshold),
    raster_value_threshold=float(snakemake.params.raster_value_threshold),
)
