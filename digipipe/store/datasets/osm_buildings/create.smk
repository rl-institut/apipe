"""
Snakefile for this dataset

Note: To include the file in the main workflow, it must be added to the respective module.smk .
"""

import json
import geopandas as gpd

from digipipe.store.utils import (
    get_abs_dataset_path,
    create_tag_string_ogr
)

DATASET_PATH = get_abs_dataset_path(
    "datasets", "osm_buildings", data_dir=True)

rule extract_buildings:
    """
    Create layers from converted OSM file using requested tags per layer. Only
    requested attributes are retained.
    """
    input:
        get_abs_dataset_path(
            "preprocessed", "osm_filtered", data_dir=True
        ) / "germany-230101_filtered.osm.gpkg"
    output: DATASET_PATH / "osm_buildings.gpkg"
    params:
        tags=create_tag_string_ogr(config["tags"]),
        layer_name=config["layer_name"]
    run:
        conditions = (
            params.tags['conditions']
            if config["use_conditions"]
            else ""
        )
        shell(
            f"ogr2ogr -f GPKG -select {params.tags['tags']} "
            f"{conditions} {output} {input} -nln {params.layer_name} "
            f"multipolygons"
        )

rule create_centroids:
    """
    Create centroids for buildings and attach ground area (in sqm)
    """
    input: DATASET_PATH / "osm_buildings.gpkg"
    output: DATASET_PATH / "osm_buildings_centroids.gpkg"
    params:
        layer_name=config["layer_name"]
    run:
        shell(
            f"ogr2ogr -f GPKG -sql "
            f"'SELECT ST_Centroid(geom) as geom, ST_Area(geom) AS area_sqm "
            f"FROM {params.layer_name}' -dialect sqlite "
            f"{output} {input} -nln '{params.layer_name}'"
        )

rule merge_with_region_file:
    """
    Merge centroids file with region file
    """
    input:
        centroids=DATASET_PATH / "osm_buildings_centroids.gpkg",
        region=rules.datasets_bkg_vg250_region_create.output
    output: DATASET_PATH / "osm_buildings_centroids_regionmerge.gpkg"
    params:
        layer_name=config["layer_name"]
    run:
        shell(
            f"ogrmerge.py -f GPKG -o {output} {input.centroids} {input.region} "
            "-nln '{{DS_BASENAME}}'"
        )

rule intersect_with_region:
    """
    Intersect centroids with region
    """
    input: DATASET_PATH / "osm_buildings_centroids_regionmerge.gpkg"
    output: DATASET_PATH / "osm_buildings_centroids_region.gpkg"
    params:
        layer_name=config["layer_name"]
    run:
        shell(
            f"ogr2ogr -sql 'SELECT a.area_sqm, a.geom "
            f"FROM osm_buildings_centroids AS a, bkg_vg250_region AS b "
            f"WHERE ST_INTERSECTS(a.geom, b.geom)' -dialect sqlite "
            f"{output} {input} -nln '{params.layer_name}'"
        )

rule calc_area_totals:
    """
    Calculate total ground area of buildings for region and country
    """
    input:
        region=DATASET_PATH / "osm_buildings_centroids_region.gpkg",
        country=DATASET_PATH / "osm_buildings_centroids.gpkg"
    output:
        region=DATASET_PATH / "osm_buildings_ground_area_region.gpkg",
        country=DATASET_PATH/ "osm_buildings_ground_area_country.gpkg"
    params:
        layer_name=config["layer_name"]
    run:
        shell(
            f"ogr2ogr -sql 'SELECT sum(area_sqm) AS area_sum_sqm "
            f"FROM {params.layer_name}' -dialect sqlite "
            f"{output.region} {input.region} -nln '{params.layer_name}'"
        )
        shell(
            f"ogr2ogr -sql 'SELECT sum(area_sqm) AS area_sum_sqm "
            f"FROM {params.layer_name}' -dialect sqlite "
            f"{output.country} {input.country} -nln '{params.layer_name}'"
        )

rule calc_building_ground_area_share:
    """
    Calculate share of region's total building ground area in country's total
    building ground area.
    """
    input:
        region=DATASET_PATH / "osm_buildings_ground_area_region.gpkg",
        country=DATASET_PATH/ "osm_buildings_ground_area_country.gpkg"
    output:
        DATASET_PATH / "osm_buildings_ground_area.json"
    run:
        area_region = round(
            float(gpd.read_file(input.region).area_sum_sqm.sum())
        )
        area_country = round(
            float(gpd.read_file(input.country).area_sum_sqm.sum())
        )
        area_share = round(area_region / area_country, 4)
        print(
            f"Share of region's total building ground area in country's total "
            f"building ground area: {area_share}"
        )
        with open(output[0], "w", encoding="utf8") as f:
            json.dump(
                {
                    "building_ground_area_country": area_country,
                    "building_ground_area_region": area_region,
                    "building_ground_area_share_region": area_share
                },
                f,
                indent=4
            )
