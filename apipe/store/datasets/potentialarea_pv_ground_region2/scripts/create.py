import os

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import shapes
from rasterstats import zonal_stats
from shapely.geometry import shape

from apipe.scripts.geo import convert_to_multipolygon, overlay, write_geofile


def process_dataset(
    clipped_raster_path, region_muns_path, output_gpkg_path, area_threshold
):
    # Step 1: Set raster values > 0 to 1 and save
    with rasterio.open(clipped_raster_path) as src:
        raster = src.read(1)
        binary_raster = np.where(raster > 0, 1, 0)
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


def raster_to_vector(raster_path):
    with rasterio.open(raster_path) as src:
        image = src.read(1)
        results = [
            {"properties": {"value": v}, "geometry": shape(s)}
            for s, v in shapes(image, transform=src.transform)
            if v > 0
        ]
        gdf = gpd.GeoDataFrame.from_features(results, crs=src.crs)
    return gdf


clipped_raster_path = snakemake.input.clipped_raster
region_muns_path = snakemake.input.region_muns
output_gpkg_path = snakemake.output.vector_overlay_gpkg
area_threshold = float(snakemake.params.area_threshold)

process_dataset(
    clipped_raster_path, region_muns_path, output_gpkg_path, area_threshold
)
